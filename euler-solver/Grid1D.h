/*!*****************************************************************************
 * \file Grid1D.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief The Grid1D struct stores the grid and provides member functions to
 * manipulate the grid
 *
 * \date 2020-06-23
 *
 * \copyright Copyright (c) 2020
 *
 * \details The Grid1D struct stores the grid as a set of four arrays, one each
 * for velocity, density, pressure, and specific internal energy. There are
 * member functions for finding the momentum, total specific energy, saving the 
 * grid, and updating the ghost cells.
 ******************************************************************************/
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
    /// The number of ghost cells on either side of the grid
    size_t numGhostCells;
    /// The number of real cells in the grid
    size_t numRealCells;
    /// The total number of cells in the grid. Equal to 2 * numGhostCells + numRealCells
    size_t numTotCells;

    /// The velcity in a specific cell. Measured in meters per second
    std::vector<double> velocity;
    /// The density in a cell. Measure in kilograms/cubic meter
    std::vector<double> density;
    /// The pressure in a cell. Measured in Pascals (Newtons per square meter)
    std::vector<double> pressure;
    /// The Specific internal energy in a cell. Measured in Joules per kilogram
    std::vector<double> siEnergy; //Specific Internal Energy

    /*!
     * \brief Compute the momentum in a cell
     *
     * \param i The index of the cell to find the momentum of
     * \return double The value of the momentum in the ith cell.
     */
    double ComputeMomentumElement(size_t const &i);

    /*!
     * \brief Compute the specific total energy in a cell.
     *
     * \param i The index of the cell in which to find the specific total energy.
     * \return double The value of the specific total energy in the ith cell.
     */
    double ComputeTotalSpecEnergyElement(size_t const &i);

    /*!
     * \brief Compute the momentum for every cell in the grid. Uses
     * Grid1D::ComputeMomentumElement to compute the momentum in each cell
     *
     * \return std::vector<double> An array of the momentum in each cell
     */
    std::vector<double> ComputeMomentumVec();

    /*!
     * \brief Compute the total specific energy for every cell in the grid. Uses
     * Grid1D::ComputeTotalSpecEnergyElement to compute the total specific energy in each cell
     *
     * \return std::vector<double> An array of the momentum in each cell
     */
    std::vector<double> ComputeTotalSpecEnergyVec();

    /*!
     * \brief Applies the boundary conditions by updating the ghost cells
     *
     * \todo Currently this only uses periodic boundary conditions and I would
     * like to add outflow, reflective, inflow, and hydrostatic boundary conditions.
     */
    void UpdateBoundaries();

    /*!
     * \brief Saves all the grid variables to their own csv files
     * \details Calling this functions saves Grid1D::velocity, Grid1D::density,
     * Grid1D::pressure, and Grid1D::siEnergy each to their own CSV files. The
     * files are stored in the directory given to the constructor
     * Grid1D::Grid1D(size_t const &reals,size_t const &ghosts,std::string const &saveloc)
     * which opens a file for each vector and saves them in that directory.
     *
     */
    void SaveState();

    // Constructors and Destructor =============================================
    /*!
     * \brief Construct a new Grid1D object that can save the primitive variables
     *
     * \details This constructor does the following:
     * - initializes Grid1D::numGhostCells, Grid1D::numRealCells, and 
     *   Grid1D::numTotCells.
     * - Reserves the memory for Grid1D::velocity, Grid1D::density,
     *   Grid1D::pressure, Grid1D::siEnergy
     * - Opens save files for Grid1D::velocity, Grid1D::density,
     *   Grid1D::pressure, Grid1D::siEnergy
     *
     * \param reals The number of real cells in the grid. Assigned to Grid1D numRealCells
     * \param ghosts The number of ghost cells. Assigned to Grid1D::numGhostCells
     * \param saveloc The directory to save files to
     *
     * \see Grid1D::Grid1D(size_t const &reals, size_t const &ghosts)
     *      This version instantiates variables but does not provide the ability
     *      to save the grid
     */
    Grid1D(size_t const &reals,
           size_t const &ghosts,
           std::string const &saveloc);

    /*!
     * \brief Construct a new Grid1D object
     *
     * \details This constructor does the following:
     * - initializes Grid1D::numGhostCells, Grid1D::numRealCells, and 
     *   Grid1D::numTotCells.
     * - Reserves the memory for Grid1D::velocity, Grid1D::density,
     *   Grid1D::pressure, Grid1D::siEnergy
     *
     * \param reals The number of real cells in the grid. Assigned to Grid1D numRealCells
     * \param ghosts The number of ghost cells. Assigned to Grid1D::numGhostCells
     *
     * \see Grid1D::Grid1D(size_t const &reals, size_t const &ghosts, std::string const &saveloc)
     *      This version can also save the grid.
     */
    Grid1D(size_t const &reals,  // This constructor is for instances that will
           size_t const &ghosts);  // never be saved

    /*!
     * \brief Destroy the Grid1D object
     * 
     * \details Closes any open files and then allows all other objects to use
     * their destructors.
     */
    ~Grid1D();
};