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

#include <H5Cpp.h>
#include <writer_hdf5.hpp>

namespace mbsolve {

static writer_factory<writer_hdf5> factory("hdf5");

writer_hdf5::writer_hdf5() : m_ext("hdf")
{
}

writer_hdf5::~writer_hdf5()
{
}

const std::string&
writer_hdf5::get_name() const
{
    return factory.get_name();
}

void
writer_hdf5::write(const std::string& filename,
                   const std::vector<std::shared_ptr<result> >& results,
                   std::shared_ptr<const device> dev,
                   std::shared_ptr<const scenario> scen) const
{
    auto h5_dbl = H5::PredType::NATIVE_DOUBLE;

    /* create file */
    H5::H5File file(filename, H5F_ACC_TRUNC);

    /* attributes */
    H5::DataSpace space_scalar;
    H5::Attribute a_timestep = file.createAttribute("timestep_size", h5_dbl,
                                                    space_scalar);
    double timestep = scen->get_timestep_size();
    a_timestep.write(h5_dbl, &timestep);

    H5::Attribute a_gridpoint = file.createAttribute("gridpoint_size", h5_dbl,
                                                     space_scalar);
    double gridpoint = scen->get_gridpoint_size();
    a_gridpoint.write(h5_dbl, &gridpoint);

    /*
     * ...
     */

    for (auto result : results) {
        /* dataset dimensions */
        const int rank = 2;
        hsize_t dims[rank];
        dims[0] = result->get_rows();
        dims[1] = result->get_cols();

        H5::DataSpace dataspace(rank, dims);

        H5::DataSet dataset = file.createDataSet(result->get_name(), h5_dbl,
                                                 dataspace);
        dataset.write(result->get_data_real_raw(), h5_dbl);

        /* TODO: create group, attribute is_complex,
         *  two datasets real/complex
         */
    }
}

const std::string&
writer_hdf5::get_extension() const
{
    return m_ext;
}

}
