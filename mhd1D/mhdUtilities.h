/*!
 * \file mhdUtilities.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains some basic utility functions commonly used in
 * magnetohydrodynamics. All functions are inlined.
 * \version 0.1
 * \date 2020-08-11
 *
 * \copyright Copyright (c) 2020
 *
 */
#pragma once

#include <stdexcept>
#include <cmath>
#include <numeric>
#include <vector>

/*!
 * \brief A namespace for common functions used in magnetohydrodynamics
 *
 */
namespace mhdUtilities
{
// =============================================================================
    /*!
     * \brief Compute the speed of the fast magnetosonic wave in the cell.
     *
     * \param pressure The pressure
     * \param density The density
     * \param magnetic The magnetic field
     * \param gamma The ratio of specific heats
     * \return double The local speed of sound, c
     */
    inline double magnetosonicSpeed(double const & pressure,
                                    double const & density,
                                    std::vector<double> magnetic,
                                    double const & gamma)
    {
        // Compute the sound speed
        double bSquared, term1, term2, cF;
        bSquared = std::inner_product(magnetic.begin(), magnetic.end(), magnetic.begin(), 0);
        term1 = gamma * pressure + bSquared;
        // Make term2 negative to compute the slow magnetic sonic speed.
        // Uncomment line below and add suitable arguments to the function to
        // enable computing either the slow or fast MS waves
        term2 = std::sqrt( std::pow(gamma * pressure + bSquared,2) - 4. * gamma * pressure * std::pow(magnetic[0], 2));
        // term2 = (waveType == "fast")? term2: -term2;

        cF = std::sqrt( (term1 + term2) / (0.5 * density) );

        // Check for Nan values in the speeds
        if (std::isnan(cF))
            throw std::runtime_error("Complex valued magnetosonic speed detected. Exiting.");

        return cF;
    }
// =============================================================================

// =============================================================================
    /*!
     * \brief Compute the speed of the Alfven wave in a cell
     *
     * \param[in] magneticX The magnetic field in the x direction, ie the direction
     * along with the Riemann solver is acting
     * \param[in] density The density in the cell
     * \return double The Alfven wave speed
     */
    inline double alfvenSpeed(double const &magneticX,
                              double const &density)
    {
        // Compute the Alfven wave speed
        double cA = std::abs(magneticX) / std::sqrt(density);

        // Check for Nan values in the speeds
        if (std::isnan(cA))
            throw std::runtime_error("Complex valued Alfven velocity detected. Exiting.");

        return cA;
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
     * \brief Compute the pressure in a cell, this is NOT the total pressure
     *
     * \param[in] energy The energy
     * \param[in] density The density
     * \param[in] velocity The velocity
     * \param[in] magnetic The magnetic field
     * \param[in] gamma The ratio of specific heats
     * \return double The pressure
     */
    inline double computePressure(double const &energy,
                                  double const &density,
                                  std::vector<double> const &velocity,
                                  std::vector<double> const &magnetic,
                                  double const &gamma)
        {return (gamma - 1)
              * ( energy
              - 0.5 * density * std::inner_product(velocity.begin(), velocity.end(), velocity.begin(), 0)
              - 0.5 * std::inner_product(magnetic.begin(), magnetic.end(), magnetic.begin(), 0));}
// =============================================================================

// =============================================================================
    /*!
     * \brief Compute the total pressure in a cell
     *
     * \param[in] energy The energy
     * \param[in] density The density
     * \param[in] velocity The velocity
     * \param[in] magnetic The magnetic field
     * \param[in] gamma The ratio of specific heats
     * \return double The total pressure
     */
    inline double computeTotalPressure(double const &energy,
                                       double const &density,
                                       std::vector<double> const &velocity,
                                       std::vector<double> const &magnetic,
                                       double const &gamma)
        {return computePressure(energy, density, velocity, magnetic, gamma)
                + 0.5 * std::inner_product(magnetic.begin(), magnetic.end(), magnetic.begin(), 0);}
// =============================================================================

// =============================================================================
    /*!
     * \brief Compute the energy in a cell
     *
     * \param[in] pressure The pressure
     * \param[in] density The density
     * \param[in] velocity The velocity
     * \param[in] magnetic The magnetic field
     * \param[in] gamma The ratio of specific heats
     * \return double The energy
     */
    inline double computeEnergy(double const &pressure,
                                double const &density,
                                std::vector<double> const &velocity,
                                std::vector<double> const &magnetic,
                                double const &gamma)
        {return (pressure/(gamma - 1))
         + 0.5 *
         ( density
         * std::inner_product(velocity.begin(), velocity.end(), velocity.begin(), 0)
         + std::inner_product(magnetic.begin(), magnetic.end(), magnetic.begin(), 0));}
// =============================================================================
}