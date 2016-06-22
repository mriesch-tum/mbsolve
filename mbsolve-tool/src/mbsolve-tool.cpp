#include <cstdlib>
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <Device.hpp>
#include <Scenario.hpp>
#include <Solver.hpp>
#include <SolverCUDA.hpp>

namespace po = boost::program_options;

static void parse_args(int argc, char **argv)
{
    po::options_description desc("Allowed options");
    desc.add_options()
	("help,h", "Print usage")
	("device,d", po::value<std::string>(), "Set device settings file")
	("scenario,s", po::value<std::string>(), "Set scenario settings file")
	("method,m", po::value<std::string>(), "Set solver method");

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
	std::cout << "Using device file " << vm["device"].as<std::string>()
		  << std::endl;
    }
    if (vm.count("scenario")) {
	std::cout << "Using scenario file " << vm["scenario"].as<std::string>()
		  << std::endl;
    }
    if (vm.count("method")) {
	std::cout << "Using solver method " << vm["method"].as<std::string>()
		  << std::endl;
    }
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    /* read inputs */

    /* setup solver */
    mbsolve::SolverCUDA solver;

    std::cout << solver.name() << std::endl;

    /* tic */

    /* execute solver */

    /* toc */

    /* write results */

    /* cleanup solver */

    /* cleanup */

    exit(0);
}

