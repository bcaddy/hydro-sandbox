/*!
 * \file HydroHelper.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains some basic utility functions commonly used in hydrodynamics.
 * All functions are inlined.
 * \version 0.1
 * \date 2020-08-11
 * 
 * \copyright Copyright (c) 2020
 * 
 */
#pragma once

#include <stdexcept>
#include <cmath>

namespace HydroHelper
{
// =============================================================================
    /*!
     * \brief Compute the local sound speed
     * 
     * \param pressure The pressure
     * \param density The density
     * \param gamma The ratio of specific heats
     * \return double The local speed of sound, c
     */
    inline double soundSpeed(double const & pressure,
                      double const & density,
                      double const & gamma)
    {
        // Compute the sound speeds
        double c = std::sqrt(gamma * pressure / density);

        // Check for Nan values in the speeds
        if (std::isnan(c))
            {throw std::runtime_error("Complex valued sound speed detected. Exiting.");}
        else
            {return c;}
    }
// =============================================================================

// =============================================================================
    /*!
     * \brief Compute the velocity in a cell
     *
     * \param momentum The momentum
     * \param density The density
     * \return double The velocity
     */
    inline double computeVelocity(double const &momentum,
                           double const &density)
        {return momentum / density;}
// =============================================================================

// =============================================================================
    /*!
     * \brief Compute the momentum in a cell
     *
     * \param velocity The velocity
     * \param density The density
     * \return double The momentum
     */
    inline double computeMomentum(double const &velocity,
                                  double const &density)
        {return velocity * density;}
// =============================================================================

// =============================================================================
    /*!
     * \brief Compute the pressure in a cell
     *
     * \param energy The energy
     * \param density The density
     * \param velocity The velocity
     * \param gamma The ratio of specific heats
     * \return double The pressure
     */
    inline double computePressure(double const &energy,
                                  double const &density,
                                  double const &velocity,
                                  double const &gamma)
        {return (gamma - 1) * ( energy - 0.5 * density * std::pow(velocity, 2) );}
// =============================================================================

// =============================================================================
    /*!
     * \brief Compute the energy in a cell
     *
     * \param pressure The pressure
     * \param density The density
     * \param velocity The velocity
     * \param gamma The ratio of specific heats
     * \return double The energy
     */
    inline double computeEnergy(double const &pressure,
                                double const &density,
                                double const &velocity,
                                double const &gamma)
        {return (pressure/(gamma - 1)) + 0.5 * density * std::pow(velocity,2);}
// =============================================================================
}
