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

int main(int argc, char **argv)
{
    /* parse command line arguments */
    parse_args(argc, argv);

    try {
	ti::cpu_timer timer;
	double total_time = 0;

        /* Ziolkowski setup */
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

        /* Ziolkowski basic scenario */
        auto scen = std::make_shared<mbsolve::scenario>
            ("Basic", 32768, 200e-15);
        //32768;
        //65536;
        //131072;
        //262144;

        scen->add_record(std::make_shared<mbsolve::record>("d11", 2e-15));
        scen->add_record(std::make_shared<mbsolve::record>("d22", 2e-15));
        scen->add_record(std::make_shared<mbsolve::record>("e", 2e-15));

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
