/*
 * mbsolve: Framework for solving the Maxwell-Bloch/-Lioville equations
 *
 * Copyright (c) 2016, Computational Photonics Group, Technical University of
 * Munich.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

/**
 * \defgroup MBSOLVE_TOOL mbsolve-tool
 * Runs different simulation setups.
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <mbsolve.hpp>

namespace po = boost::program_options;
namespace ti = boost::timer;

static std::string device_file;
static std::string output_file;
static std::string scenario_file;
static std::string solver_method;
static std::string writer_method;
static mbsolve::real sim_endtime;
static unsigned int num_gridpoints;

static void parse_args(int argc, char **argv)
{
    po::options_description desc("Allowed options");
    desc.add_options()
	("help,h", "Print usage")
	("device,d", po::value<std::string>(&device_file)->required(),
	 "Set device settings file")
	("output,o", po::value<std::string>(&output_file), "Set output file")
	("scenario,s", po::value<std::string>(&scenario_file),
	 "Set scenario settings file")
	("method,m", po::value<std::string>(&solver_method)->required(),
	 "Set solver method")
	("writer,w", po::value<std::string>(&writer_method)->required(),
	 "Set writer")
        ("endtime,e", po::value<mbsolve::real>(&sim_endtime),
         "Set simulation end time")
        ("gridpoints,g", po::value<unsigned int>(&num_gridpoints),
         "Set number of spatial grid points");

    po::variables_map vm;
    try {
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
	    std::cout << desc;
	    exit(0);
	}
	po::notify(vm);
    } catch (po::error& e) {
	std::cerr << "Error: " << e.what() << std::endl;
	std::cerr << desc;
	exit(1);
    }

    if (vm.count("device")) {
	device_file = vm["device"].as<std::string>();
	std::cout << "Using device file " << device_file << std::endl;
    }
    if (vm.count("scenario")) {
	scenario_file = vm["scenario"].as<std::string>();
	std::cout << "Using scenario file " << scenario_file << std::endl;
    }
}

/**
 * mbsolve-tool main function.
 * \ingroup MBSOLVE_TOOL
 *
 * Specify the simulation setup using the -d parameter (The available setups
 * are described below.) and the solver method using the -m parameter
 * (Currently, the only option available is openmp-Xlvl-os-red. Replace X
 * with the number of levels.) For the complete list of parameters, run
 * mbsolve-tool -h.
 */
int main(int argc, char **argv)
{
    /* parse command line arguments */
    parse_args(argc, argv);

    try {
        ti::cpu_timer timer;
        double total_time = 0;

        std::shared_ptr<mbsolve::device> dev;
        std::shared_ptr<mbsolve::scenario> scen;

        if (device_file == "song2005") {
            /**
             * The song2005 setup features a three-level active region which
             * is excited by a sech pulse. For details see literature:
             * Song X. et al., Propagation of a Few-Cycle Laser Pulse
             * in a V-Type Three-Level System, Optics and Spectroscopy,
             * Vol. 99, No. 4, 2005, pp. 517–521
             * https://doi.org/10.1134/1.2113361
             */

            /* set up quantum mechanical description */
            std::vector<mbsolve::real> energies = {
                0,
                2.3717 * mbsolve::HBAR * 1e15,
                2.4165 * mbsolve::HBAR * 1e15
            };

            mbsolve::qm_operator H(energies);

            std::vector<mbsolve::complex> dipoles = {
                mbsolve::E0 * 9.2374e-11,
                mbsolve::E0 * 9.2374e-11 * sqrt(2),
                0
            };

            mbsolve::qm_operator u({0, 0, 0}, dipoles);

            mbsolve::real rate = 1e10;
            std::vector<std::vector<mbsolve::real> > scattering_rates = {
                { 0, rate, rate },
                { rate, 0, rate },
                { rate, rate, 0 } };

            auto relax_sop = std::make_shared<mbsolve::qm_lindblad_relaxation>
                (scattering_rates);

            mbsolve::qm_operator rho_init({1, 0, 0});

            auto qm = std::make_shared<mbsolve::qm_description>
                (6e24, H, u, relax_sop, rho_init);

            auto mat_vac = std::make_shared<mbsolve::material>("Vacuum");
            mbsolve::material::add_to_library(mat_vac);
            auto mat_ar = std::make_shared<mbsolve::material>("AR_Song", qm);
            mbsolve::material::add_to_library(mat_ar);

            dev = std::make_shared<mbsolve::device>("Song");
            dev->add_region(std::make_shared<mbsolve::region>
                            ("Active region", mat_ar, 0, 150e-6));

            /* default settings */
            if (num_gridpoints == 0) {
                num_gridpoints = 32768;
            }
            if (sim_endtime < 1e-21) {
                sim_endtime = 80e-15;
            }

            /* Song basic scenario */
            scen = std::make_shared<mbsolve::scenario>
                ("Basic", num_gridpoints, sim_endtime);

            auto sech_pulse = std::make_shared<mbsolve::sech_pulse>
                ("sech", 0.0, mbsolve::source::hard_source, 3.5471e9,
                 3.8118e14, 17.248, 1.76/5e-15, -M_PI/2);
            scen->add_source(sech_pulse);

            scen->add_record(std::make_shared<mbsolve::record>
                             ("d11", mbsolve::record::type::density, 1, 1,
                              0, 0));
            scen->add_record(std::make_shared<mbsolve::record>
                             ("d22", mbsolve::record::type::density, 2, 2,
                              0, 0));
            scen->add_record(std::make_shared<mbsolve::record>
                             ("d33", mbsolve::record::type::density, 3, 3,
                              0, 0));
            scen->add_record(std::make_shared<mbsolve::record>("e", 0, 0.0));


        } else if (device_file == "ziolkowski1995") {
            /**
             * The ziolkowski1995 setup is a self induced transparency (SIT)
             * setup that consists of a two-level active region embedded in
             * two vacuum section. For details see literature:
             * https://doi.org/10.1103/PhysRevA.52.3082
             */

            /* set up quantum mechanical description */
            auto qm = std::make_shared<mbsolve::qm_desc_2lvl>
                (1e24, 2 * M_PI * 2e14, 6.24e-11, 1.0e10, 1.0e10);

            /* materials */
            auto mat_vac = std::make_shared<mbsolve::material>("Vacuum");
            auto mat_ar = std::make_shared<mbsolve::material>
                ("AR_Ziolkowski", qm);
            mbsolve::material::add_to_library(mat_vac);
            mbsolve::material::add_to_library(mat_ar);

            /* set up device */
            dev = std::make_shared<mbsolve::device>("Ziolkowski");
            dev->add_region(std::make_shared<mbsolve::region>
                            ("Vacuum left", mat_vac, 0, 7.5e-6));
            dev->add_region(std::make_shared<mbsolve::region>
                            ("Active region", mat_ar, 7.5e-6, 142.5e-6));
            dev->add_region(std::make_shared<mbsolve::region>
                            ("Vacuum right", mat_vac, 142.5e-6, 150e-6));

            /* default settings */
            if (num_gridpoints == 0) {
                num_gridpoints = 32768;
            }
            if (sim_endtime < 1e-21) {
                sim_endtime = 200e-15;
            }

            /* Ziolkowski basic scenario */
            scen = std::make_shared<mbsolve::scenario>
                ("Basic", num_gridpoints, sim_endtime);

            auto sech_pulse = std::make_shared<mbsolve::sech_pulse>
                //("sech", 0.0, mbsolve::source::hard_source, 4.2186e9/2, 2e14,
                ("sech", 0.0, mbsolve::source::hard_source, 4.2186e9, 2e14,
                 10, 2e14);
            scen->add_source(sech_pulse);

            scen->add_record(std::make_shared<mbsolve::record>
                             ("inv12", 2.5e-15));
            scen->add_record(std::make_shared<mbsolve::record>("e", 2.5e-15));

        } else if (device_file == "tzenov2018-cpml") {
            /**
             * The tzenov2018-cpml setup consists of an absorber region
             * embedded in two gain regions. Each region is modeled as a
             * two-level system. The results show that short pulses are
             * generated due to colliding pulse mode-locking (CPML).
             * For details see literature:
             * https://doi.org/10.1088/1367-2630/aac12a
             */

            /* set up quantum mechanical descriptions */
            auto qm_gain = std::make_shared<mbsolve::qm_desc_2lvl>
                (5e21, 2 * M_PI * 3.4e12, 2e-9, 1.0/10e-12, 1.0/200e-15, 1.0);

            auto qm_absorber = std::make_shared<mbsolve::qm_desc_2lvl>
                (1e21, 2 * M_PI * 3.4e12, 6e-9, 1.0/3e-12, 1.0/160e-15);


#if 0
            /* initial value density matrix */
            d_init << 0.5, 0.001,
                0.001, 0.5;
#endif
            /* materials */
            auto mat_absorber = std::make_shared<mbsolve::material>
                ("Absorber", qm_absorber, 12.96, 1.0, 500);
            auto mat_gain = std::make_shared<mbsolve::material>
                ("Gain", qm_gain, 12.96, 1.0, 500);
            mbsolve::material::add_to_library(mat_absorber);
            mbsolve::material::add_to_library(mat_gain);

            /* set up device */
            dev = std::make_shared<mbsolve::device>("tzenov-cpml");
            dev->add_region(std::make_shared<mbsolve::region>
                            ("Gain R", mat_gain, 0, 0.5e-3));
            dev->add_region(std::make_shared<mbsolve::region>
                            ("Absorber", mat_absorber, 0.5e-3, 0.625e-3));
            dev->add_region(std::make_shared<mbsolve::region>
                            ("Gain L", mat_gain, 0.625e-3, 1.125e-3));

            /* default settings */
            if (num_gridpoints == 0) {
                num_gridpoints = 8192;
            }
            if (sim_endtime < 1e-21) {
                sim_endtime = 2e-9;
            }

            /* basic scenario */
            scen = std::make_shared<mbsolve::scenario>
                ("basic", num_gridpoints, sim_endtime);

            scen->add_record(std::make_shared<mbsolve::record>
                             ("inv12", 1e-12));
            scen->add_record(std::make_shared<mbsolve::record>("e", 1e-12));
            scen->add_record(std::make_shared<mbsolve::record>
                             ("e0", mbsolve::record::electric, 1, 1, 0.0,
                              0.0));
            scen->add_record(std::make_shared<mbsolve::record>
                             ("h0", mbsolve::record::magnetic, 1, 1, 0.0,
                              1.373e-7));
        } else {
            throw std::invalid_argument("Specified device not found!");
        }
        /* tic */
        timer.start();

        mbsolve::writer writer(writer_method);
        mbsolve::solver solver(solver_method, dev, scen);

        /* toc */
        timer.stop();
        ti::cpu_times times = timer.elapsed();
        std::cout << "Time required (setup): " << 1e-9 * times.wall
                  << std::endl;
        total_time +=1e-9 * times.wall;

        std::cout << solver.get_name() << std::endl;

        /* tic */
        timer.start();

        /* execute solver */
        solver.run();

        /* toc */
        timer.stop();
        times = timer.elapsed();
        std::cout << "Time required (run): " << 1e-9 * times.wall << std::endl;
        total_time +=1e-9 * times.wall;

        /* grid point updates per second */
        double gpups = 1e-6 * 1e9/times.wall * scen->get_num_gridpoints() *
            scen->get_endtime()/scen->get_timestep_size();
        std::cout << "Performance: " << gpups << " MGPU/s" << std::endl;

        /* tic */
        timer.start();

        /* write results */
        writer.write(output_file, solver.get_results(), dev, scen);

        /* toc */
        timer.stop();
        times = timer.elapsed();
        std::cout << "Time required (write): " << 1e-9 * times.wall
                  << std::endl;
        total_time +=1e-9 * times.wall;

        std::cout << "Time required (total): " << total_time << std::endl;

    } catch (std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        exit(1);
    }

    exit(0);
}
