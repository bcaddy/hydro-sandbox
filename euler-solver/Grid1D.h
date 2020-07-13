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
 * for velocity, density, and pressure. There are member functions for finding
 * the momentum, saving the grid, and updating the ghost
 * cells.
 ******************************************************************************/
#pragma once

#include <vector>
#include <string>
#include <fstream>

struct Grid1D
{
private:
    // Output files
    std::ofstream _velocitySaveFile;
    std::ofstream _densitySaveFile;
    std::ofstream _pressureSaveFile;

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

    /*!
     * \brief Compute the momentum in a cell
     *
     * \param[in] i The index of the cell to find the momentum of
     * \return double The value of the momentum in the ith cell.
     */
    double computeMomentumElement(size_t const &i);

    /*!
     * \brief Compute the momentum for every cell in the grid. Uses
     * Grid1D::ComputeMomentumElement to compute the momentum in each cell
     *
     * \return std::vector<double> An array of the momentum in each cell
     */
    std::vector<double> computeMomentumVec();

    /*!
     * \brief Applies the boundary conditions by updating the ghost cells
     *
     * \todo Currently this only uses periodic boundary conditions and I would
     * like to add outflow, reflective, inflow, and hydrostatic boundary conditions.
     */
    void updateBoundaries();

    /*!
     * \brief Saves all the grid variables to their own csv files
     * \details Calling this functions saves Grid1D::velocity, Grid1D::density,
     * and Grid1D::pressure each to their own CSV files. The
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
     * - Reserves the memory for Grid1D::velocity, Grid1D::density, and
     *   Grid1D::pressure
     * - Opens save files for Grid1D::velocity, Grid1D::density, and
     *   Grid1D::pressure
     *
     * \param[in] reals The number of real cells in the grid. Assigned to Grid1D numRealCells
     * \param[in] ghosts The number of ghost cells. Assigned to Grid1D::numGhostCells
     * \param[in] saveDir The directory to save files to. If set to "no saving" then the
     * initialized grid will not be able to save itself to a file.
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
              std::string const &saveDir);

    // Constructors and Destructor =============================================
    /*!
     * \brief Construct a new uninitialized Grid1D object for initialization
     * later using the Grid1D::Init() method
     */
    Grid1D() = default;

    /*!
     * \brief Construct a new Grid1D object that can save the primitive variables
     *
     * \details This constructor utilizes the Grid1D::Init() method to initialize
     *          a Grid1D object with the ability to save the gride to a file.
     *
     * \param[in] reals The number of real cells in the grid. Assigned to Grid1D::numRealCells
     * \param[in] ghosts The number of ghost cells. Assigned to Grid1D::numGhostCells
     * \param[in] saveDir The directory to save files to
     *
     * \see Grid1D::Init() The initializing function used by this constructor
     * \see Grid1D::Grid1D(size_t const &reals, size_t const &ghosts)
     *      Constructor that instantiates variables but does not provide the
     *      ability to save the grid
     *
     */
    Grid1D(size_t const &reals,
           size_t const &ghosts,
           std::string const &saveDir);

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