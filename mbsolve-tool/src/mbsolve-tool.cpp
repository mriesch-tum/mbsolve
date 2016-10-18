#include <cstdlib>
#include <iostream>
#include <string>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <Device.hpp>
#include <Scenario.hpp>
#include <Solver.hpp>
#include <Writer.hpp>

namespace po = boost::program_options;
namespace ti = boost::timer;

static std::string device_file;
static std::string scenario_file;
static std::string solver_method;

/* TODO: output file parameter ? */
/* TODO: writer parameter ? */

static void parse_args(int argc, char **argv)
{
    po::options_description desc("Allowed options");
    desc.add_options()
	("help,h", "Print usage")
	("device,d", po::value<std::string>(), "Set device settings file")
	("scenario,s", po::value<std::string>(), "Set scenario settings file")
	("method,m", po::value<std::string>(&solver_method)->required(),
	 "Set solver method");

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

mbsolve::Device parse_device(const std::string& file)
{
    mbsolve::Device dev;

    /* TODO: read xml file */

    /* Ziolkowski settings */
    mbsolve::Region vacuum;
    vacuum.XDim = 7.5e-6;
    vacuum.Overlap = 1;
    vacuum.Losses = 0;
    vacuum.DopingDensity = 0;

    vacuum.Name = "Vacuum left";
    vacuum.X0 = 0;
    dev.Regions.push_back(vacuum);

    mbsolve::Region active;
    active.Name = "Active Region";
    active.X0 = 7.5e-6;
    active.XDim = 135e-6;
    active.Overlap = 1;
    active.Losses = 0;
    active.DopingDensity = 1e24;
    active.TransitionFrequencies.push_back(mbsolve::Quantity(M_PI * 4e14));
    active.DipoleMoments.push_back(mbsolve::Quantity(6.24e-9));
    //active.ScatteringRates.push_back(mbsolve::Quantity(-1.0e10));
    active.ScatteringRates.push_back(mbsolve::Quantity(+1.0e10));
    active.DephasingRates.push_back(mbsolve::Quantity(1.0e10));
    dev.Regions.push_back(active);

    vacuum.Name = "Vacuum right";
    vacuum.X0 = 142.5e-6;
    dev.Regions.push_back(vacuum);

    return dev;
}

mbsolve::Scenario parse_scenario(const std::string& file)
{
    mbsolve::Scenario scen;

    /* TODO: read xml file */

    /* Ziolkowski settings */
    scen.Name = "Ziolkowski";
    scen.SimEndTime = 5e-9;
    scen.NumGridPoints = 5760;

    return scen;
}

// TODO: setup - cleanup sequence? make sure everything is cleanup in case of an error
/* TODO: add Solver mysolver(args) that does the setup/RAII scheme*/
/* TODO: analogue for writer */
/* TODO: abstract IWriter/ISolver etc */
int main(int argc, char **argv)
{
    mbsolve::Device device;
    mbsolve::Scenario scenario;
    mbsolve::Solver *solver;
    mbsolve::Writer *writer;
    ti::cpu_timer timer;
    std::vector<mbsolve::Result *> results;

    /* parse command line arguments */
    parse_args(argc, argv);

    /* parse device file */
    try {
	device = parse_device(device_file);
    } catch (std::exception& e) {
	std::cout << "Error: Could not parse device file " << device_file
		  << std::endl << e.what() << std::endl;
	exit(1);
    }

    /* parse scenario file */
    try {
	scenario = parse_scenario(scenario_file);
    } catch (std::exception& e) {
	std::cout << "Error: Could not parse scenario file " << scenario_file
		  << std::endl << e.what() << std::endl;
	exit(1);
    }

    /* select writer */
    try {
	writer = mbsolve::Writer::create("MATLAB");
    } catch (std::exception& e) {
	std::cout << "Error: " << e.what() << std::endl;
	exit(1);
    }

    /* select solver */
    try {
	solver = mbsolve::Solver::create(solver_method);
    } catch (std::exception& e) {
	std::cout << "Error: " << e.what() << std::endl;
	exit(1);
    }

    /* setup solver */
    solver->setup(device, scenario);

    std::cout << solver->name() << std::endl;

    /* tic */
    timer.start();

    /* execute solver */
    solver->run(results);

    /* toc */
    timer.stop();
    ti::cpu_times times = timer.elapsed();
    std::cout << "Time required: " << 1e-9 * times.wall << std::endl;

    /* write results */
    try {
	writer->write("test.mat", results, device, scenario);
    } catch (std::exception& e) {
	std::cout << "Error: Could not write results."
		  << std::endl << e.what() << std::endl;
    }

     /* cleanup */
    BOOST_FOREACH(mbsolve::Result *result, results) {
	delete result;
    }
    delete solver;
    delete writer;

    exit(0);
}
