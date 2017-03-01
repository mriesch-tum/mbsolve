#include <solver_2lvl_pc_red.hpp>
#include <iostream>
#include <boost/foreach.hpp>

namespace mbsolve{

static SolverFactory<SolverOMP_2lvl_pc_red> factory("openmp-2lvl-pc-red");

struct sim_constants gsc[MaxRegions];
unsigned int OL;

SolverOMP_2lvl_pc_red::SolverOMP_2lvl_pc_red(const Device& device,
					     const Scenario& scenario) :
    ISolver(device, scenario)
{
    /* total device length */
    Quantity length = device.XDim();

    /* minimum relative permittivity */
    Quantity minRelPermittivity = device.MinRelPermittivity();

    /* determine grid point and time step size */
    real C = 0.5; /* courant number */
    real velocity = sqrt(MU0() * EPS0() * minRelPermittivity());
    m_scenario.GridPointSize = length()/(m_scenario.NumGridPoints - 1);
    real timestep  = C * m_scenario.GridPointSize * velocity;
    m_scenario.NumTimeSteps = ceil(m_scenario.SimEndTime/timestep) + 1;
    m_scenario.TimeStepSize = m_scenario.SimEndTime /
	(m_scenario.NumTimeSteps - 1);

    /* determine border indices and initialize region settings */
    if (device.Regions.size() > MaxRegions) {
	throw std::invalid_argument("Too many regions requested");
    }

    unsigned int i = 0;
    BOOST_FOREACH(Region reg, device.Regions) {
	if (i > 0) {
	    gsc[i - 1].idx_end = round(reg.X0()/m_scenario.GridPointSize) - 1;
	}
	gsc[i].idx_start = round(reg.X0()/m_scenario.GridPointSize);
	gsc[i].M_CE = m_scenario.TimeStepSize/(EPS0() * reg.RelPermittivity());
	gsc[i].M_CH = m_scenario.TimeStepSize/(MU0() *
					      m_scenario.GridPointSize);
	/* TODO: overlap factor? */
	gsc[i].M_CP = -2.0 * reg.DopingDensity * HBAR();
	gsc[i].sigma = 2.0 * sqrt(EPS0 * reg.RelPermittivity/MU0) * reg.Losses;

	gsc[i].w12 = (reg.TransitionFrequencies.size() < 1) ? 0.0 :
	    reg.TransitionFrequencies[0]();
	/* TODO rename to rabi freqs or something */
	gsc[i].d12 = (reg.DipoleMoments.size() < 1) ? 0.0 :
	    reg.DipoleMoments[0]() * E0()/HBAR();
	gsc[i].tau1 = (reg.ScatteringRates.size() < 1) ? 0.0 :
	    reg.ScatteringRates[0]();
	gsc[i].gamma12 = (reg.DephasingRates.size() < 1) ? 0.0 :
	    reg.DephasingRates[0]();

	gsc[i].d_x_inv = 1.0/m_scenario.GridPointSize;
	gsc[i].d_t = m_scenario.TimeStepSize;

	if (reg.DopingDensity() < 1.0) {
	    gsc[i].dm11_init = 0.0;
	    gsc[i].dm22_init = 0.0;
	} else {
	    gsc[i].dm11_init = 0.0;
	    gsc[i].dm22_init = 1.0;
	}

	i++;
    }
    if (i > 0) {
	gsc[i - 1].idx_end = m_scenario.NumGridPoints - 1;
    }

    /* redundant calculation overlap */
    OL = 1;

    unsigned int P = omp_get_max_threads();
    unsigned int chunk = m_scenario.NumGridPoints/P;

    std::cout << "Number of threads: " << P << std::endl;
    m_dm11 = new real*[P];
    m_dm12r = new real*[P];
    m_dm12i = new real*[P];
    m_dm22 = new real*[P];
    m_e = new real*[P];
    m_h = new real*[P];
    region_indices = new unsigned int*[P];

    for (int tid = 0; tid < P; tid++) {
	unsigned int size = chunk + 2 * OL;

	if (tid == P - 1) {
	    size += m_scenario.NumGridPoints % P;
	}

	/* allocation */
	m_dm11[tid] = new real[size];
	m_dm12r[tid] = new real[size];
	m_dm12i[tid] = new real[size];
	m_dm22[tid] = new real[size];

	m_h[tid] = new real[size];
	m_e[tid] = new real[size];

	region_indices[tid] = new unsigned int[size];
    }

    std::cout << "m_dm11 " << m_dm11 << std::endl;
    std::cout << "m_dm22 " << m_dm22 << std::endl;
    std::cout << "m_e " << m_e << std::endl;

    /* set up results transfer data structures */
    BOOST_FOREACH(Record rec, m_scenario.Records) {
	unsigned int row_ct = m_scenario.SimEndTime/rec.Interval;
	unsigned int interval = ceil(1.0 * m_scenario.NumTimeSteps/row_ct);
	unsigned int position_idx;
	unsigned int col_ct;

	if (rec.Position() < 0.0) {
	    /* copy complete grid */
	    position_idx = 0;
	    col_ct = m_scenario.NumGridPoints;
	} else {
	    position_idx = round(rec.Position()/m_scenario.GridPointSize);
	    col_ct = 1;
	}

	/* allocate result memory */
	Result *res = new Result(rec.Name, col_ct, row_ct);
	m_results.push_back(res);

	/* create copy list entry */
	CopyListEntry *entry;
	if (rec.Type == EField) {
	    entry = new CLEField_red(m_e, res, col_ct, position_idx,
				 interval);
	    m_copyList.push_back(entry);
	} else if (rec.Type == HField) {
	    /* TODO: numGridPoints + 1 */
	    entry = new CLEField_red(m_h, res, col_ct, position_idx,
				 interval);
	    m_copyList.push_back(entry);
	} else if (rec.Type == Density) {
	    if ((rec.I - 1 < 2) && (rec.J - 1 < 2)) {
		if (rec.I == rec.J) {
		    /* main diagonal entry */
		    if (rec.I == 1) {
			entry = new CLEField_red(m_dm11, res, col_ct, position_idx,
					     interval);
		    } else {
			entry = new CLEField_red(m_dm22, res, col_ct, position_idx,
					     interval);
		    }
		    m_copyList.push_back(entry);
		} else {
		    /* off-diagonal entry */
		    /* TODO */
		    /* if complex */
		    /* create two list entries */
		    /* create two Results, or one complex Result */

		    /* real part: GetSrcDensity(&dm, rec.I, rec.J); */
		    /* imag part: GetSrcDensity(&dm, rec.J, rec.I); */
		}
	    } else {
	    // throw exc
	    }
	} else {
	    // throw exc
	}
    }
}

SolverOMP_2lvl_pc_red::~SolverOMP_2lvl_pc_red()
{
    /* delete copy lists */
    BOOST_FOREACH(CopyListEntry *entry, m_copyList) {
	delete entry;
    }

    for (int tid = 0; tid < omp_get_max_threads(); tid++) {
	delete[] m_h[tid];
	delete[] m_e[tid];

	delete[] m_dm11[tid];
	delete[] m_dm12r[tid];
	delete[] m_dm12i[tid];
	delete[] m_dm22[tid];

	delete[] region_indices[tid];
    }

    delete[] m_dm11;
    delete[] m_dm12r;
    delete[] m_dm12i;
    delete[] m_dm22;

    delete[] m_e;
    delete[] m_h;
    delete[] region_indices;

}

std::string
SolverOMP_2lvl_pc_red::getName() const
{
    return factory.getName();
}

void
SolverOMP_2lvl_pc_red::run() const
{
#pragma omp parallel
    {
	struct sim_constants sc[MaxRegions];

	for (unsigned int j = 0; j < MaxRegions; j++) {
	    sc[j] = gsc[j];
	}

	unsigned int P = omp_get_max_threads();
	unsigned int tid = omp_get_thread_num();
	unsigned int chunk_base = m_scenario.NumGridPoints/P;
	unsigned int chunk = chunk_base;

	if (tid == P - 1) {
	    chunk += m_scenario.NumGridPoints % P;
	}

	/* initialization */
	for (int i = 0; i < chunk + 2 * OL; i++) {
	    unsigned int region = 0;
	    int global_idx = tid * chunk_base + (i - OL);

	    if ((global_idx >= 0) && (global_idx < m_scenario.NumGridPoints)) {
		for (unsigned int j = 0; j < MaxRegions; j++) {
		    if (global_idx <= sc[j].idx_end) {
			region = j;
			break;
		    }
		}
		m_dm11[tid][i] = sc[region].dm11_init;
		m_dm22[tid][i] = sc[region].dm22_init;
		region_indices[tid][i] = region;
	    } else {
		region_indices[tid][i] = 0;
	    }

	    m_dm12r[tid][i] = 0.0;
	    m_dm12i[tid][i] = 0.0;
	    m_e[tid][i] = 0.0;
	    m_h[tid][i] = 0.0;
	}
#pragma omp barrier

	/* main loop */
	for (unsigned int n = 0; n < m_scenario.NumTimeSteps/OL; n++) {

	    /* exchange data */
	    if (tid > 0) {
		for (int i = 0; i < OL; i++) {
		    m_dm11[tid][i] = m_dm11[tid - 1][chunk + i];
		    m_dm12r[tid][i] = m_dm12r[tid - 1][chunk + i];
		    m_dm12i[tid][i] = m_dm12i[tid - 1][chunk + i];
		    m_dm22[tid][i] = m_dm22[tid - 1][chunk + i];

		    m_e[tid][i] = m_e[tid - 1][chunk + i];
		    m_h[tid][i] = m_h[tid - 1][chunk + i];
		}
	    }

	    if (tid < P - 1) {
		for (int i = 0; i < OL; i++) {
		    m_dm11[tid][OL + chunk + i] = m_dm11[tid + 1][OL + i];
		    m_dm12r[tid][OL + chunk + i] = m_dm12r[tid + 1][OL + i];
		    m_dm12i[tid][OL + chunk + i] = m_dm12i[tid + 1][OL + i];
		    m_dm22[tid][OL + chunk + i] = m_dm22[tid + 1][OL + i];

		    m_e[tid][OL + chunk + i] = m_e[tid + 1][OL + i];
		    m_h[tid][OL + chunk + i] = m_h[tid + 1][OL + i];
		}
	    }

	    /* sync after communication */
#pragma omp barrier

	    /* sub-loop */
	    for (unsigned int m = 0; m < OL; m++) {

		/* update dm and e */
		for (int i = m; i < chunk + 2 * OL - m; i++) {
		    int region = region_indices[tid][i];

		    real rho11_e = m_dm11[tid][i];
		    real rho12r_e = m_dm12r[tid][i];
		    real rho12i_e = m_dm12i[tid][i];
		    real rho22_e = m_dm22[tid][i];
		    real e_e = m_e[tid][i];

		    for (int pc_step = 0; pc_step < 4; pc_step++) {
			/* execute prediction - correction steps */

			real rho11  = 0.5 * (m_dm11[tid][i] + rho11_e);
			real rho22  = 0.5 * (m_dm22[tid][i] + rho22_e);
			real rho12r = 0.5 * (m_dm12r[tid][i] + rho12r_e);
			real rho12i = 0.5 * (m_dm12i[tid][i] + rho12i_e);
			real OmRabi = 0.5 * sc[region].d12 *
			    (m_e[tid][i] + e_e);

			rho11_e = m_dm11[tid][i] + sc[region].d_t *
			    (- 2.0 * OmRabi * rho12i -
			     sc[region].tau1 * rho11);

			rho12i_e = m_dm12i[tid][i] + sc[region].d_t *
			    (- sc[region].w12 * rho12r
			     + OmRabi * (rho11 - rho22)
			     - sc[region].gamma12 * rho12i);

			rho12r_e = m_dm12r[tid][i] + sc[region].d_t *
			    (+ sc[region].w12 * rho12i
			     - sc[region].gamma12 * rho12r);

			rho22_e = m_dm22[tid][i] + sc[region].d_t *
			    (+ 2.0 * OmRabi * rho12i
			     + sc[region].tau1 * rho11);

			real j = 0;
			real p_t = sc[region].M_CP * sc[region].d12 *
			    (sc[region].w12 * rho12i -
			     sc[region].gamma12 * rho12r);

			e_e = m_e[tid][i] + sc[region].M_CE *
			    (-j - p_t + (m_h[tid][i + 1] - m_h[tid][i]) *
			     sc[region].d_x_inv);

			/* TODO fix source condition */
			if (tid * chunk_base + (i - OL) == 0) {
			    /* calculate source value */
			    real f_0 = 2e14;
			    real t = (n * OL + m) * m_scenario.TimeStepSize;
			    real T_p = 20/f_0;
			    real gamma = 2 * t/T_p - 1;
			    real E_0 = 4.2186e9;
			    real src = E_0 * 1/std::cosh(10 * gamma) *
				sin(2 * M_PI * f_0 * t);
			    src = src/2; // pi pulse

			    e_e = src; /* hard source */
			}
		    }

		    /* final update step */
		    m_dm11[tid][i] = rho11_e;
		    m_dm12i[tid][i] = rho12i_e;
		    m_dm12r[tid][i] = rho12r_e;
		    m_dm22[tid][i] = rho22_e;

		    m_e[tid][i] = e_e;
		}

		/* update h */
		for (int i = m; i < chunk + 2 * OL - m; i++) {
		    int region = region_indices[tid][i];

		    m_h[tid][i] += sc[region].M_CH *
			(m_e[tid][i] - m_e[tid][i - 1]);

		    if (tid * chunk_base + (i - OL) == 0) {
			m_h[tid][i] = 0;
		    }
		}

		/* copy result data */
		BOOST_FOREACH(CopyListEntry *entry, m_copyList) {
		    if (entry->record(n * OL + m)) {
			std::copy(entry->getSrc(tid) + OL,
				  entry->getSrc(tid) + OL + chunk,
				  entry->getDst(n * OL + m)
				  + tid * chunk_base);
		    }
		}
	    }

	    /* sync after computation */
#pragma omp barrier
	}
    }
}

}
