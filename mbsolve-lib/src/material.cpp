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

#include <material.hpp>

namespace mbsolve {

std::map<std::string, std::shared_ptr<material> >
material::m_materials;

void
material::add_to_library(std::shared_ptr<material> mat)
{
    auto it = m_materials.find(mat->get_id());
    if (it != m_materials.end()) {
        throw std::invalid_argument("Material already available in library.");
    }
    m_materials[mat->get_id()] = mat;
}

void
material::add_to_library(const material& mat)
{
    auto ptr = std::make_shared<material>(mat);
    material::add_to_library(ptr);
}

std::shared_ptr<material>
material::get_from_library(const std::string& id)
{
    auto it = m_materials.find(id);
    if (it == m_materials.end()) {
        throw std::invalid_argument("Material not available in library.");
    }
    return it->second;
}

/*
void
material::add_to_library() const
{

}*/

}
