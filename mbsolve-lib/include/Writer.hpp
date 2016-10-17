#ifndef WRITER_H
#define WRITER_H

#include <map>
#include <string>
#include <vector>
#include <Device.hpp>
#include <Scenario.hpp>
#include <Types.hpp>

namespace mbsolve {

class IWriterFactory;

class Writer
{
private:
    static std::map<std::string, IWriterFactory *> m_writers;

public:
    virtual void write(const std::string& file,
		       const std::vector<Result *>& results,
		       const Device& device,
		       const Scenario& scenario)
    {
    }

    static void register_new(const std::string& name, IWriterFactory *factory);

    static Writer *create(const std::string& name);
};

class IWriterFactory
{
public:
    virtual Writer* create() = 0;
};

template<typename T>
class WriterFactory : IWriterFactory
{
public:
    explicit WriterFactory(const std::string& name) {
	Writer::register_new(name, this);
    }

    Writer* create() { return new T; }
};

}

#endif
