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
    /// The total number of cells including ghost cells
    size_t numTotCells;

    /// The density in a cell. Measure in kilograms/cubic meter
    stdVector1D density;
    /// The velocity in a specific cell. Measured in meters per second
    stdVector2D velocity;
    /// The magnetic on the i-1/2 face of a specific cell.
    stdVector2D magnetic;
    /// The pressure in a cell. Measured in Pascals
    stdVector1D pressure;

    /*!
     * \brief Construct a new Primitive Grid 1 D object
     *
     * \param[in] reals Number of real cells
     * \param[in] ghosts Number of ghost cells
     */
    PrimitiveGrid1D(size_t const &reals,
                    size_t const &ghosts)
        :
        numTotCells(reals + 2 * ghosts),
        density(numTotCells),
        velocity(numTotCells, std::vector<double> (3, 0)),
        magnetic(numTotCells, std::vector<double> (3, 0)),
        pressure(numTotCells)
    {};

    /*!
     * \brief Destroy the Primitive Grid 1 D object with default destructor
     */
    ~PrimitiveGrid1D() = default;
};