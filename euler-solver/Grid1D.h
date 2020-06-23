/*!
 * \file Grid1D.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief The Grid1D class stores the grid and provides member functions to
 * manipulate the grid
 *
 * \version 0.1
 * \date 2020-06-23
 *
 * \copyright Copyright (c) 2020
 *
 * The Grid1D class stores the grid as a set of four arrays, one each for velocity,
 * density, pressure, and specific internal energy
 */
#pragma once

#include <vector>
#include <string>
#include <fstream>

struct Grid1D
{
private:
    // Output files
    std::ofstream VelocitySaveFile;
    std::ofstream DensitySaveFile;
    std::ofstream PressureSaveFile;
    std::ofstream siEnergySaveFile;

public:
    // Properties of the grid
    size_t numGhostCells;
    size_t numRealCells;
    size_t numTotCells;

    // Grid values
    std::vector<double> velocity;
    std::vector<double> density;
    std::vector<double> pressure;
    std::vector<double> siEnergy; //Specific Internal Energy

    // Get the conserved variables
    double ComputeMomentumElement(size_t const &i);
    double ComputeTotalSpecEnergyElement(size_t const &i);
    std::vector<double> ComputeMomentumVec();
    std::vector<double> ComputeTotalSpecEnergyVec();

    // Update the boundary conditions
    void UpdateBoundaries();

    // Save the array
    void SaveState();

    // Constructors and Destructor =============================================
    Grid1D(size_t const &reals,
           size_t const &ghosts,
           std::string const &saveloc);
    Grid1D(size_t const &reals,  // This constructor is for instances that will
           size_t const &ghosts);  // never be saved
    ~Grid1D();
};


