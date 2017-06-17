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

std::map<std::string, writer_factory_int *>
writer_int::m_factories;

writer_int::writer_int()
{
}

writer_int::~writer_int()
{
}

void
writer_int::register_factory(const std::string& name,
                             writer_factory_int *factory)
{
    if (m_factories[name]) {
	throw std::invalid_argument("Writer already registered.");
    }
    m_factories[name] = factory;
}

writer_factory_int *
writer_int::find_factory(const std::string& name)
{
    auto it = m_factories.find(name);
    if (it == m_factories.end()) {
        throw std::invalid_argument("Unknown writer " + name);
    }
    return it->second;
}

writer::writer(const std::string& name)
{
    /* create writer */
    writer_factory_int *factory = writer_int::find_factory(name);
    m_writer = factory->create_instance();
}

writer::~writer()
{
}

void
writer::write(const std::string& file,
              const std::vector<std::shared_ptr<result> >& results,
              std::shared_ptr<const device> dev,
              std::shared_ptr<const scenario> scen) const
{
    std::string def = dev->get_name() + "_" + scen->get_name() + "." +
        m_writer->get_extension();

    m_writer->write(file.empty() ? def : file, results, dev, scen);
}

const std::string&
writer::get_extension() const
{
    return m_writer->get_extension();
}

}
