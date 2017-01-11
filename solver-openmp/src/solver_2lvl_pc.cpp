#include <solver_2lvl_pc.hpp>

#include <boost/foreach.hpp>

namespace mbsolve{

static SolverFactory<SolverOMP_2lvl_pc> factory("openmp-2lvl-pc");

static struct sim_constants gsc[MaxRegions];

unsigned int get_region(unsigned int idx)
{
    for (unsigned int i = 0; i < MaxRegions; i++) {
	if (idx <= gsc[i].idx_end) {
	    return i;
	}
    }
    return 0;
}

SolverOMP_2lvl_pc::SolverOMP_2lvl_pc(const Device& device,
				     const Scenario& scenario) :
    ISolver(device, scenario)
{
    /* total device length */
    Quantity length = device.XDim();

    /* minimum relative permittivity */
    Quantity minRelPermittivity = device.MinRelPermittivity();

    /* TODO: sanity check scenario? */
    if (m_scenario.NumGridPoints % 32 != 0) {
	throw std::invalid_argument("Number of grid points must be multiple"
				    " of 32");
    }

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
    struct sim_constants sc[MaxRegions];

    unsigned int i = 0;
    BOOST_FOREACH(Region reg, device.Regions) {
	if (i > 0) {
	    sc[i - 1].idx_end = round(reg.X0()/m_scenario.GridPointSize) - 1;
	}
	sc[i].idx_start = round(reg.X0()/m_scenario.GridPointSize);
	sc[i].M_CE = m_scenario.TimeStepSize/(EPS0() * reg.RelPermittivity());
	sc[i].M_CH = m_scenario.TimeStepSize/(MU0() *
					      m_scenario.GridPointSize);
	/* TODO: overlap factor? */
	sc[i].M_CP = -2.0 * reg.DopingDensity * HBAR();
	sc[i].sigma = 2.0 * sqrt(EPS0 * reg.RelPermittivity/MU0) * reg.Losses;

	sc[i].w12 = (reg.TransitionFrequencies.size() < 1) ? 0.0 :
	    reg.TransitionFrequencies[0]();
	/* TODO rename to rabi freqs or something */
	sc[i].d12 = (reg.DipoleMoments.size() < 1) ? 0.0 :
	    reg.DipoleMoments[0]() * E0()/HBAR();
	sc[i].tau1 = (reg.ScatteringRates.size() < 1) ? 0.0 :
	    reg.ScatteringRates[0]();
	sc[i].gamma12 = (reg.DephasingRates.size() < 1) ? 0.0 :
	    reg.DephasingRates[0]();

	sc[i].d_x = m_scenario.GridPointSize;
	sc[i].d_t = m_scenario.TimeStepSize;

	if (reg.DopingDensity() < 1.0) {
	    sc[i].dm11_init = 0.0;
	    sc[i].dm22_init = 0.0;
	} else {
	    sc[i].dm11_init = 0.0;
	    sc[i].dm22_init = 1.0;
	}

	i++;
    }
    if (i > 0) {
	sc[i - 1].idx_end = m_scenario.NumGridPoints - 1;
    }

m_dm = new DensityMatrixData(m_scenario.NumGridPoints, 2, 5);

}

SolverOMP_2lvl_pc::~SolverOMP_2lvl_pc()
{
delete m_dm;
}

std::string
SolverOMP_2lvl_pc::getName() const
{
    return factory.getName();
}

void
SolverOMP_2lvl_pc::run() const
{
    /* main loop */
    for (unsigned int i = 0; i < m_scenario.NumTimeSteps; i++) {


	/* calculate source value */
	/* TODO */
	real f_0 = 2e14;
	real t = i * m_scenario.TimeStepSize;
	real T_p = 20/f_0;
	real gamma = 2 * t/T_p - 1;
	real E_0 = 4.2186e9;
	real src = E_0 * 1/std::cosh(10 * gamma) * sin(2 * M_PI * f_0 * t);
	src = src/2; // pi pulse


	/* execute prediction - correction steps */
	for (int pc_step = 0; pc_step < 4; pc_step++) {
	    #pragma omp parallel sections
	    {
	    #pragma omp section
	    {
	    // TODO: estimate dm in parallel

	}
	    #pragma omp section
	    {
	    // TODO: estimate e in parallel
	}
	}

}

	// TODO save result

	// TODO OMP parallel section
	// TODO update dm in parallel
	// TODO update e in parallel

	// TODO update h in parallel

	// TODO save result
    }
}

}
