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

#include <solver_generic.hpp>

namespace mbsolve {

static solver_factory<solver_generic> factory("generic");

solver_generic::solver_generic(std::shared_ptr<const device> dev,
                               std::shared_ptr<scenario> scen) :
    solver_int(dev, scen)
{

    for (auto rec : scen->get_records()) {

        unsigned int cols = 20;
        unsigned int rows = 10;

        auto res = std::make_shared<result>(rec->get_name(), cols, rows);

        m_results.push_back(res);
    }
}

solver_generic::~solver_generic()
{
}

const std::string&
solver_generic::get_name() const
{
    return factory.get_name();
}

void
solver_generic::run() const
{
    /* fill test results */

    std::vector<real> rpart(20);
    std::vector<real> ipart(20);

    for (auto res : m_results) {
        for (unsigned int row_ct = 0; row_ct < 10; row_ct++) {
            for (auto& r : rpart) {
                r = 10.0 * row_ct;
            }

            for (auto& i : ipart) {
                i = 42.0;
            }

            std::copy(rpart.cbegin(), rpart.cend(),
                      res->get_data_real(row_ct, 0));
            std::copy(ipart.cbegin(), ipart.cend(),
                      res->get_data_imag(row_ct, 0));
        }
    }
}

}
