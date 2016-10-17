#include <Writer.hpp>

namespace mbsolve {

std::map<std::string, IWriterFactory *> Writer::m_writers;

void
Writer::register_new(const std::string& name, IWriterFactory *factory)
{
    // TODO: check doubles
    m_writers[name] = factory;
}

Writer *
Writer::create(const std::string& name)
{
    std::map<std::string, IWriterFactory *>::iterator it;
    it = m_writers.find(name);
    if (it == m_writers.end()) {
	throw std::invalid_argument("Unknown writer " + name);
    }
    return it->second->create();
}

}
