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
    /// Save file object for the momentum
    std::ofstream _momentumSaveFile;
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
    std::vector<double> momentum;
    /// The energy in a cell. Measured in Joules
    std::vector<double> energy;

    /// The type of boundary conditions to use
    std::string boundaryConditionKind;

    /*!
     * \brief Applies the boundary conditions by updating the ghost cells
     *
     * \todo Currently this only uses periodic boundary conditions and I would
     * like to add outflow, reflective, inflow, and hydrostatic boundary conditions.
     *
     * \param[in] gamma The ratio of specific heats
     */
    void updateBoundaries(double const &gamma);

    /*!
     * \brief Saves all the grid variables to their own csv files
     * \details Calling this functions saves Grid1D::density, Grid1D::momentum,
     * and Grid1D::energy each to their own CSV files. The
     * files are stored in the directory given to the constructor
     * Grid1D::Grid1D(size_t const &reals,size_t const &ghosts,std::string const &saveDir)
     * which opens a file for each vector and saves them in that directory.
     *
     */
    void saveState();

    /*!
     * \brief Construct a new Grid1D object
     *
     * \details This method provides the actual construction of the object so
     * that the object can either be constructed immediately upon declaration or
     * at a later time. It does the following:
     * - initializes Grid1D::numGhostCells, Grid1D::numRealCells, and
     *   Grid1D::numTotCells.
     * - Reserves the memory for  Grid1D::density, Grid1D::momentum, and
     *   Grid1D::energy
     * - Opens save files for Grid1D::density, Grid1D::momentum, and
     *   Grid1D::energy
     *
     * \param[in] reals The number of real cells in the grid. Assigned to Grid1D numRealCells
     * \param[in] ghosts The number of ghost cells. Assigned to Grid1D::numGhostCells
     * \param[in] saveDir The directory to save files to. If set to "no saving" then the
     * initialized grid will not be able to save itself to a file.
     * \param[in] boundaryConditions What kind of boundary conditions to use
     *
     * \see Grid1D::Grid1D(size_t const &reals, size_t const &ghosts)
     *      Constructor that instantiates variables but does not provide the
     *      ability to save the grid
     * \see Grid1D::Grid1D(size_t const &reals, size_t const &ghosts, std::string const &saveDir)
     *      Constructor that instantiates variables and does provide the ability
     *      to save the grid
     */
    void init(size_t const &reals,
              size_t const &ghosts,
              std::string const &saveDir,
              std::string const &boundaryConditions);

    // Constructors and Destructor =============================================
    /*!
     * \brief Construct a new uninitialized Grid1D object for initialization
     * later using the Grid1D::Init() method
     */
    Grid1D() = default;

    /*!
     * \brief Construct a new Grid1D object that can save the conserved variables
     *
     * \details This constructor utilizes the Grid1D::Init() method to initialize
     *          a Grid1D object with the ability to save the gride to a file.
     *
     * \param[in] reals The number of real cells in the grid. Assigned to Grid1D::numRealCells
     * \param[in] ghosts The number of ghost cells. Assigned to Grid1D::numGhostCells
     * \param[in] saveDir The directory to save files to
     * \param[in] boundaryConditions What kind of boundary conditions to use
     *
     * \see Grid1D::Init() The initializing function used by this constructor
     * \see Grid1D::Grid1D(size_t const &reals, size_t const &ghosts)
     *      Constructor that instantiates variables but does not provide the
     *      ability to save the grid
     *
     */
    Grid1D(size_t const &reals,
           size_t const &ghosts,
           std::string const &saveDir,
           std::string const &boundaryConditions);

    /*!
     * \brief Construct a new Grid1D object that cannot save itself to a file
     *
     * \details This constructor utilizes the Grid1D::Init() method to initialize
     *          a Grid1D object with the ability to save the gride to a file.
     *
     * \param[in] reals The number of real cells in the grid. Assigned to Grid1D::numRealCells
     * \param[in] ghosts The number of ghost cells. Assigned to Grid1D::numGhostCells
     *
     * \see Grid1D::Init() The initializing function used by this constructor
     * \see Grid1D::Grid1D(size_t const &reals, size_t const &ghosts, std::string const &saveDir)
     *      Constructor that instantiates variables and does provide the ability
     *      to save the grid
     *
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