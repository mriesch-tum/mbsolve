/*
 * mbsolve: An open-source solver tool for the Maxwell-Bloch equations.
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

#ifndef MBSOLVE_WRITER_HDF5_H
#define MBSOLVE_WRITER_HDF5_H

#include <mbsolve/lib/writer.hpp>

namespace mbsolve {

/**
 * \defgroup MBSOLVE_WRITER_HDF5 writer-hdf5
 * Writer for Hierarchic Data Format (HDF) 5 format.
 */

/**
 * Writer class for Hierarchic Data Format (HDF) 5 format.
 * \ingroup MBSOLVE_WRITER_HDF5
 */
class writer_hdf5 : public writer
{
public:
    /**
     * Constructs writer for Hierarchic Data Format (HDF) 5 format.
     */
    writer_hdf5();

    /**
     * Destructs writer for Hierarchic Data Format (HDF) 5 format.
     */
    ~writer_hdf5();

    /**
     * Writes results to HDF5 file.
     */
    void write(
        const std::string& file,
        const std::vector<std::shared_ptr<result> >& results,
        std::shared_ptr<const device> dev,
        std::shared_ptr<const scenario> scen) const;
};

/**
 * Helper class that enables usage of the HDF5 writer.
 */
class writer_hdf5_loader
{
public:
    writer_hdf5_loader();
};
}

#endif
