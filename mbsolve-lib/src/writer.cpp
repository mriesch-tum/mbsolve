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

#include <writer.hpp>

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
    std::string def = device.get_name() + "-" + scenario.Name + "." +
        m_writer->getExtension();

    m_writer->write(file.empty() ? def : file, results, device, scenario);
}

}
