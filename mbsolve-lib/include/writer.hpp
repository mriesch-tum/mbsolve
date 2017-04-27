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
