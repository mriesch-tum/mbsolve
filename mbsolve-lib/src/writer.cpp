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

std::map<std::string, writer::bootstrap_t> writer::m_bootstraps;

std::shared_ptr<writer>
writer::create_instance(const std::string& name)
{
    auto it = m_bootstraps.find(name);
    if (it == m_bootstraps.end()) {
        throw std::invalid_argument("Unknown writer " + name);
    }
    return it->second();
}

void
writer::register_bootstrap(const std::string& name, bootstrap_t b)
{
    if (m_bootstraps[name]) {
        throw std::invalid_argument("Writer already registered.");
    }
    m_bootstraps[name] = b;
}

writer::writer(const std::string& name, const std::string& extension)
  : m_name(name), m_ext(extension)
{
}

writer::~writer()
{
}

const std::string&
writer::get_name() const
{
    return m_name;
}

const std::string&
writer::get_extension() const
{
    return m_ext;
}

std::vector<std::string>
writer::get_avail_writers()
{
    std::vector<std::string> writers;

    for (const auto& s : m_bootstraps) {
        writers.push_back(s.first);
    }

    return writers;
}
}
