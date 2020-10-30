/*!*****************************************************************************
 * \file Grid1D.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains the declaration of the Grid1D struct
 *
 * \date 2020-06-23
 *
 * \copyright Copyright (c) 2020
 *
 ******************************************************************************/
#pragma once

#include <vector>
#include <string>
#include <fstream>

/*!
 * \brief The Grid1D struct stores the grid and provides member functions to
 * manipulate the grid
 *
 * \details The Grid1D struct stores the grid as a set of three arrays, one each
 * for density, momentum, and energy. There are member functions for saving the
 * grid and updating the ghost cells/boundary conditions.
 */
struct Grid1D
{
private:
    // Output files

    /// Save file object for the density
    std::ofstream _densitySaveFile;

    /// Save file object for the momentum in the X direction
    std::ofstream _momentumXSaveFile;
    /// Save file object for the momentum in the Y direction
    std::ofstream _momentumYSaveFile;
    /// Save file object for the momentum in the Z direction
    std::ofstream _momentumZSaveFile;

    /// Save file object for the magnetic field in the X direction
    std::ofstream _magneticXSaveFile;
    /// Save file object for the magnetic field in the Y direction
    std::ofstream _magneticYSaveFile;
    /// Save file object for the magnetic field in the Z direction
    std::ofstream _magneticZSaveFile;

    /// Save file object for the energy
    std::ofstream _energySaveFile;

public:
    /// The number of ghost cells on either side of the grid
    size_t numGhostCells;
    /// The number of real cells in the grid
    size_t numRealCells;
    /// The total number of cells in the grid. Equal to 2 * numGhostCells + numRealCells
    size_t numTotCells;

    /// The density in a cell. Measure in kilograms/cubic meter
    std::vector<double> density;
    /// The momentum in a specific cell. Measured in kilogram meters per second
    std::vector<std::vector<double>> momentum;
    /// The magnetic field on the i-1/2 face of a specific cell.
    std::vector<std::vector<double>> magnetic;
    /// The energy in a cell. Measured in Joules
    std::vector<double> energy;

    /// The type of boundary conditions to use
    std::string boundaryConditionKind;

    /*!
     * \brief Applies the boundary conditions by updating the ghost cells
     *
     * \param[in] gamma The ratio of specific heats
     */
    void updateBoundaries(double const &gamma);

    /*!
     * \brief Saves all the grid variables to their own csv files
     * \details Calling this functions saves the entire grid, each conserved
     * variable to their own CSV files. The files are stored in the directory
     * given to the constructor Grid1D::Grid1D(size_t const &reals,size_t const
     * &ghosts,std::string const &saveDir) which opens a file for each vector
     * and saves them in that directory.
     *
     */
    void saveState();

    // Constructors and Destructor =============================================
    /*!
     * \brief Construct a new Grid1D object
     *
     * \details Constructs the Grid1D object. It does the following:
     * - initializes Grid1D::numGhostCells, Grid1D::numRealCells, and
     *   Grid1D::numTotCells.
     * - Reserves the memory for  Grid1D::density, Grid1D::momentum, and
     *   Grid1D::energy
     * - Optionally: Opens save files for Grid1D::density, Grid1D::momentum, and
     *   Grid1D::energy
     *
     * \param[in] reals The number of real cells in the grid. Assigned to Grid1D::numRealCells
     * \param[in] ghosts The number of ghost cells. Assigned to Grid1D::numGhostCells
     * \param[in] boundaryConditions What kind of boundary conditions to use. Defaults to "pass"
     * \param[in] saveDir Optional: The directory to save files to. If no
     *            directory is provided then saving is disabled
     */
    Grid1D(size_t const &reals,
           size_t const &ghosts,
           std::string const &boundaryConditions = "pass",
           std::string const &saveDir = "0");

    /*!
     * \brief Destroy the Grid1D object
     *
     * \details Closes any open files and then allows all other objects to use
     * their destructors.
     */
    ~Grid1D();
};