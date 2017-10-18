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

static void parse_args(int argc, char **argv)
{
    po::options_description desc("Allowed options");
    desc.add_options()
	("help,h", "Print usage")
	("device,d", po::value<std::string>(&device_file),
	 "Set device settings file")
	("output,o", po::value<std::string>(&output_file), "Set output file")
	("scenario,s", po::value<std::string>(&scenario_file),
	 "Set scenario settings file")
	("method,m", po::value<std::string>(&solver_method)->required(),
	 "Set solver method")
	("writer,w", po::value<std::string>(&writer_method)->required(),
	 "Set writer");

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

Eigen::Matrix<mbsolve::complex, 3, 3>
relax_sop_song3(const Eigen::Matrix<mbsolve::complex, 3, 3>& arg)
{
    Eigen::Matrix<mbsolve::complex, 3, 3> ret =
        Eigen::Matrix<mbsolve::complex, 3, 3>::Zero();

    mbsolve::real d7_eq = 0.0;
    mbsolve::real d8_eq = 0.0;
    mbsolve::real T_1 = 1e-10;

    ret(0, 0) = 1.0/3 * (1/T_1 * (arg(2, 2) + arg(1, 1) - 2 * arg(0, 0)
                                  - d7_eq - d8_eq));
    ret(1, 1) = 1.0/3 * (1/T_1 * (arg(2, 2) + arg(1, 1) - 2 * arg(0, 0)
                                  - d7_eq - d8_eq))
        - 1/T_1 * (arg(1, 1) - arg(0, 0) - d7_eq);
    ret(2, 2) = 1.0/3 * (1/T_1 * (arg(2, 2) + arg(1, 1) - 2 * arg(0, 0)
                                  - d7_eq - d8_eq))
        - 1/T_1 * (arg(2, 2) - arg(0, 0) - d8_eq);
    ret(1, 0) = -1/T_1 * arg(1, 0);
    ret(0, 1) = -1/T_1 * arg(0, 1);
    ret(2, 0) = -1/T_1 * arg(2, 0);
    ret(0, 2) = -1/T_1 * arg(0, 2);
    ret(2, 1) = -1/T_1 * arg(2, 1);
    ret(1, 2) = -1/T_1 * arg(1, 2);

    return ret;
}

Eigen::Matrix<mbsolve::complex, 3, 3>
relax_sop_ziolk3(const Eigen::Matrix<mbsolve::complex, 3, 3>& arg)
{
    Eigen::Matrix<mbsolve::complex, 3, 3> ret =
        Eigen::Matrix<mbsolve::complex, 3, 3>::Zero();

    ret(0, 0) = +1e10 * arg(1, 1);
    ret(1, 1) = -1e10 * arg(1, 1);
    ret(0, 1) = -1e10 * arg(0, 1);
    ret(1, 0) = -1e10 * arg(1, 0);

    return ret;
}

Eigen::Matrix<mbsolve::complex, 2, 2>
relax_sop_ziolk2(const Eigen::Matrix<mbsolve::complex, 2, 2>& arg)
{
    Eigen::Matrix<mbsolve::complex, 2, 2> ret =
        Eigen::Matrix<mbsolve::complex, 2, 2>::Zero();

    ret(0, 0) = +1e10 * arg(1, 1);
    ret(1, 1) = -1e10 * arg(1, 1);
    ret(0, 1) = -1e10 * arg(0, 1);
    ret(1, 0) = -1e10 * arg(1, 0);

    return ret;
}

int main(int argc, char **argv)
{
    /* parse command line arguments */
    parse_args(argc, argv);

    try {
	ti::cpu_timer timer;
	double total_time = 0;

#define DEVICE 4
#define SCENARIO 2

#if DEVICE==1
        /* Song setup */

        Eigen::Matrix<mbsolve::complex, 3, 3> H, u, d_init;

        H <<0, 0, 0,
            0, 2.3717, 0,
            0, 0, 2.4165;
        H = H * mbsolve::HBAR * 1e15;

        // mbsolve::real g = 1.0;
        mbsolve::real g = sqrt(2);

        u <<0, 1.0, g,
            1.0, 0, 0,
            g, 0, 0;
        u = u * mbsolve::E0 * 9.2374e-11;

        d_init << 1, 0, 0,
            0, 0, 0,
            0, 0, 0;

        auto qm = std::make_shared<mbsolve::qm_desc_3lvl>
            (6e24, H, u, &relax_sop_song3, d_init);

        auto mat_vac = std::make_shared<mbsolve::material>("Vacuum");
        mbsolve::material::add_to_library(mat_vac);
        auto mat_ar = std::make_shared<mbsolve::material>("AR_Song", qm);
        mbsolve::material::add_to_library(mat_ar);

        auto dev = std::make_shared<mbsolve::device>("Song");
        //        dev->add_region(std::make_shared<mbsolve::region>
        //                        ("Vacuum left", mat_vac, 0, 7.5e-6));
        dev->add_region(std::make_shared<mbsolve::region>
                        //("Active region", mat_ar, 7.5e-6, 150e-6));
                        ("Active region", mat_ar, 0, 150e-6));

#elif DEVICE==2
        /* Ziolkowski setup in 3-lvl desc */

        Eigen::Matrix<mbsolve::complex, 3, 3> H, u, d_init;

        H <<-0.5, 0, 0,
            0, 0.5, 0,
            0, 0, 0;
        H = H * mbsolve::HBAR * 2 * M_PI * 2e14;
        u <<0, 1.0, 0,
            1.0, 0, 0,
            0, 0, 0;
        u = u * mbsolve::E0 * 6.24e-11 * (-1);
        d_init << 1, 0, 0,
            0, 0, 0,
            0, 0, 0;

        auto qm = std::make_shared<mbsolve::qm_desc_3lvl>
            (1e24, H, u, &relax_sop_ziolk3, d_init);

        auto mat_vac = std::make_shared<mbsolve::material>("Vacuum");
        auto mat_ar = std::make_shared<mbsolve::material>("AR_Ziolkowski", qm);
        mbsolve::material::add_to_library(mat_vac);
        mbsolve::material::add_to_library(mat_ar);

        /* set up device */
        auto dev = std::make_shared<mbsolve::device>("Ziolkowski");
        dev->add_region(std::make_shared<mbsolve::region>
                        ("Vacuum left", mat_vac, 0, 7.5e-6));
        dev->add_region(std::make_shared<mbsolve::region>
                        ("Active region", mat_ar, 7.5e-6, 142.5e-6));
        dev->add_region(std::make_shared<mbsolve::region>
                        ("Vacuum right", mat_vac, 142.5e-6, 150e-6));
        /* end */
#elif DEVICE==3
        /* Ziolkowski setup in old 2-lvl desc */
        auto qm = std::make_shared<mbsolve::qm_desc_2lvl>
            (1e24, 2 * M_PI * 2e14, 6.24e-11, 0.5e10, 1.0e10);

        auto mat_vac = std::make_shared<mbsolve::material>("Vacuum");
        auto mat_ar = std::make_shared<mbsolve::material>("AR_Ziolkowski", qm);
        mbsolve::material::add_to_library(mat_vac);
        mbsolve::material::add_to_library(mat_ar);

        /* set up device */
        auto dev = std::make_shared<mbsolve::device>("Ziolkowski");
        dev->add_region(std::make_shared<mbsolve::region>
                        ("Vacuum left", mat_vac, 0, 7.5e-6));
        dev->add_region(std::make_shared<mbsolve::region>
                        ("Active region", mat_ar, 7.5e-6, 142.5e-6));
        dev->add_region(std::make_shared<mbsolve::region>
                        ("Vacuum right", mat_vac, 142.5e-6, 150e-6));
        /* end */
#elif DEVICE==4
        /* Ziolkowski setup in new 2-lvl desc */

        Eigen::Matrix<mbsolve::complex, 2, 2> H, u, d_init;

        H <<-0.5, 0,
            0, 0.5;
        H = H * mbsolve::HBAR * 2 * M_PI * 2e14;
        u <<0, 1.0,
            1.0, 0;
        u = u * mbsolve::E0 * 6.24e-11 * (-1.0);
        d_init << 1, 0,
            0, 0;
        auto qm = std::make_shared<mbsolve::qm_desc_clvl<2> >
            (1e24, H, u, &relax_sop_ziolk2, d_init);

        auto mat_vac = std::make_shared<mbsolve::material>("Vacuum");
        auto mat_ar = std::make_shared<mbsolve::material>("AR_Ziolkowski", qm);
        mbsolve::material::add_to_library(mat_vac);
        mbsolve::material::add_to_library(mat_ar);

        /* set up device */
        auto dev = std::make_shared<mbsolve::device>("Ziolkowski");
        dev->add_region(std::make_shared<mbsolve::region>
                        ("Vacuum left", mat_vac, 0, 7.5e-6));
        dev->add_region(std::make_shared<mbsolve::region>
                        ("Active region", mat_ar, 7.5e-6, 142.5e-6));
        dev->add_region(std::make_shared<mbsolve::region>
                        ("Vacuum right", mat_vac, 142.5e-6, 150e-6));
        /* end */
#else

#endif

#if SCENARIO==1
        /* Song basic scenario */
        auto scen = std::make_shared<mbsolve::scenario>
            //("Basic", 32768, 10e-15);
            ("Basic", 32768, 200e-15);
            //("Basic", 65536, 200e-15);
            //("Basic", 131072, 200e-15);
            //("Basic", 262144, 200e-15);

        auto sech_pulse = std::make_shared<mbsolve::sech_pulse>
            ("sech", 0.0, mbsolve::source::hard_source, 3.5471e9, 3.8118e14,
             17.248, 1.76/5e-15, -M_PI/2);
        scen->add_source(sech_pulse);

        scen->add_record(std::make_shared<mbsolve::record>("inv12", 0, 0.0));
        scen->add_record(std::make_shared<mbsolve::record>("e", 0, 0.0));
        //scen->add_record(std::make_shared<mbsolve::record>("e", 2.5e-15));
        //scen->add_record(std::make_shared<mbsolve::record>("inv12", 2.5e-15));

#elif SCENARIO==2
        /* Ziolkowski basic scenario */

        auto scen = std::make_shared<mbsolve::scenario>
            //("Basic", 32768, 10e-15);
            ("Basic", 32768, 200e-15);
            //("Basic", 65536, 200e-15);
            //("Basic", 131072, 200e-15);
            //("Basic", 262144, 200e-15);

        auto sech_pulse = std::make_shared<mbsolve::sech_pulse>
            //("sech", 0.0, mbsolve::source::hard_source, 4.2186e9/2, 2e14,
            ("sech", 0.0, mbsolve::source::hard_source, 4.2186e9, 2e14,
             10, 2e14);
        scen->add_source(sech_pulse);

        scen->add_record(std::make_shared<mbsolve::record>("inv12", 2.5e-15));
        scen->add_record(std::make_shared<mbsolve::record>("e", 2.5e-15));
#else

#endif


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
