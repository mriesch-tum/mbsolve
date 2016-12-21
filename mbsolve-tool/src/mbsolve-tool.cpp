#include <cstdlib>
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <Device.hpp>
#include <Record.hpp>
#include <Scenario.hpp>
#include <Solver.hpp>
#include <Writer.hpp>

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

mbsolve::Device parse_device(const std::string& file)
{
    mbsolve::Device dev;
    dev.Name = "Ziolkowski";

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
    //    active.X0 = 0;
    //    active.XDim = 150e-6;
    active.X0 = 7.5e-6;
    active.XDim = 135e-6;
    active.Overlap = 1;
    active.Losses = 0;
    active.DopingDensity = 1e24;
    //active.DopingDensity = 1e16;
    active.TransitionFrequencies.push_back(mbsolve::Quantity(M_PI * 4e14));
    active.DipoleMoments.push_back(mbsolve::Quantity(6.24e-9));
    active.ScatteringRates.push_back(mbsolve::Quantity(0.5e10));
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
    scen.Name = "Basic";
    scen.SimEndTime = 100e-15;
    //scen.SimEndTime = 500e-15;
    //scen.SimEndTime = 500e-12;
    //scen.NumGridPoints = 23040;
    scen.NumGridPoints = 32768;
    //scen.NumGridPoints = 65536;
    //scen.NumGridPoints = 131072;
    //scen.NumGridPoints = 262144;

    scen.Records.push_back(mbsolve::Record("dm11", mbsolve::Density, 1, 1,
					   1e-15));

    scen.Records.push_back(mbsolve::Record("dm22", mbsolve::Density, 2, 2,
					   1e-15));
    scen.Records.push_back(mbsolve::Record("e", mbsolve::EField, 1, 1,
					   1e-15));
    return scen;
}

int main(int argc, char **argv)
{
    mbsolve::Device device;
    mbsolve::Scenario scenario;

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

    try {
	mbsolve::Writer writer(writer_method);
	mbsolve::Solver solver(solver_method, device, scenario);

	std::cout << solver.getName() << std::endl;

	/* tic */
	ti::cpu_timer timer;
	timer.start();

	/* execute solver */
	solver.run();

	/* toc */
	timer.stop();
	ti::cpu_times times = timer.elapsed();
	std::cout << "Time required: " << 1e-9 * times.wall << std::endl;

	/* write results */
	writer.write(output_file, solver.getResults(), device,
		     solver.getScenario());

    } catch (std::exception& e) {
	std::cout << "Error: " << e.what() << std::endl;
	exit(1);
    }

    exit(0);
}
