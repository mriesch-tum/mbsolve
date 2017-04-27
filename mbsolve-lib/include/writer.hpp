#ifndef WRITER_H
#define WRITER_H

#include <map>
#include <string>
#include <vector>
#include <device.hpp>
#include <scenario.hpp>
#include <types.hpp>

namespace mbsolve {

class IWriterFactory;

class IWriter
{
public:
    IWriter() { }
    virtual ~IWriter() { }
    virtual std::string getExtension() const = 0;
    virtual void write(const std::string& file,
                       const std::vector<Result *>& results,
                       const Device& device,
                       const Scenario& scenario) const = 0;
};


class Writer
{
private:
    static std::map<std::string, IWriterFactory *> m_factories;
    IWriter *m_writer;

public:
    Writer(const std::string& name);

    ~Writer();

    void write(const std::string& file,
               const std::vector<Result *>& results,
               const Device& device,
               const Scenario& scenario) const;

    static void registerFactory(const std::string& name,
				IWriterFactory *factory);
};

class IWriterFactory
{
public:
    IWriterFactory() { }
    virtual ~IWriterFactory() { }
    virtual IWriter* createInstance() const = 0;
};

template<typename T>
class WriterFactory : IWriterFactory
{
public:
    explicit WriterFactory(const std::string& name) {
        Writer::registerFactory(name, this);
    }

    IWriter* createInstance() const { return new T; }

};

}

#endif
