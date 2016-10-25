#include <Writer.hpp>

namespace mbsolve {

std::map<std::string, IWriterFactory *> Writer::m_factories;

void
Writer::registerFactory(const std::string& name, IWriterFactory *factory)
{
    if (m_factories[name]) {
	throw std::invalid_argument("Writer already registered.");
    }
    m_factories[name] = factory;
}

Writer::Writer(const std::string& name)
{
    std::map<std::string, IWriterFactory *>::iterator it;
    it = m_factories.find(name);
    if (it == m_factories.end()) {
	throw std::invalid_argument("Unknown writer " + name);
    }
    m_writer = it->second->createInstance();
}

Writer::~Writer()
{
    delete m_writer;
}

void
Writer::write(const std::string& file, const std::vector<Result *>& results,
	      const Device& device, const Scenario& scenario) const
{
    std::string def = device.Name + "-" + scenario.Name + "." +
	m_writer->getExtension();

    m_writer->write(file.empty() ? def : file, results, device, scenario);
}

}
