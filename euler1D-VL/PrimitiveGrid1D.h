/*!
 * \file PrimitiveGrid1D.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains a simple struct for a 1D grid of primitive variables along
 * with the constructor and destructor
 * \version 0.1
 * \date 2020-09-29
 *
 * \copyright Copyright (c) 2020
 *
 */
#pragma once

#include <vector>

/*!
 * \brief A simple struct to for a 1D grid of primitive (density, velocity,
 * pressure) variables
 *
 */
struct PrimitiveGrid1D
{
public:
    /// The density in a cell. Measure in kilograms/cubic meter
    std::vector<double> density;
    /// The velocity in a specific cell. Measured in meters per second
    std::vector<double> velocity;
    /// The pressure in a cell. Measured in Pascals
    std::vector<double> pressure;

    /*!
     * \brief Construct a new Primitive Grid 1 D object
     *
     * \param[in] reals Number of real cells
     * \param[in] ghosts Number of ghost cells
     */
    PrimitiveGrid1D(size_t const &reals,
                    size_t const &ghosts)
    {
        size_t numTotCells = reals + ghosts;
        density.resize(numTotCells);
        velocity.resize(numTotCells);
        pressure.resize(numTotCells);
    };

    /*!
     * \brief Destroy the Primitive Grid 1 D object with default destructor
     */
    ~PrimitiveGrid1D() = default;
};