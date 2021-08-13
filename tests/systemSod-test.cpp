/*!
 * \file systemSod-test.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Perform a system test of the Euler1D-VL program with the sod shock tube
 * \version 0.1
 * \date 2021-07-29
 *
 * \copyright Copyright (c) 2021
 *
 */


// STL includes
#include <stdlib.h>
#include <future>
#include <string>
#include <algorithm>

// Include GoogleTest and related libraries/headers
#include <gtest/gtest.h>

// Local includes
#include "testingUtilities.h"
using namespace testingUtilities;

// Lets start testing
TEST(SystemTest_VanLeer, sod_HLLC_1k_cells)
{
    // First let's start our cholla analogue running asynchronously. Note that
    // this dumps all console I/O to /dev/nul
    // auto sodProcess = std::async(std::launch::async, // Launch operation asynchronously
    //                              system,             // Choose the function to launch
    //                              "../euler1D-VL/euler-solver.exe >/dev/null 2>&1"); // Args to send to "system" call
    system("../euler1D-VL/euler-solver.exe >/dev/null 2>&1");
    // While that's running we're going to load the fiducial data and find the
    // number of fiducial time steps
    std::string fidDensity  = file2String("System-Test-Data/sod-VL-HLLC/Density.csv");
    std::string fidEnergy   = file2String("System-Test-Data/sod-VL-HLLC/Energy.csv");
    std::string fidMomentum = file2String("System-Test-Data/sod-VL-HLLC/Momentum.csv");
    size_t numFidTimesteps  = std::count(fidDensity.begin(), fidDensity.end(), '\n') - 1;
    std::cout << "Finished fiducial import" << std::endl;

    // Wait for sodProcess to finish
    // sodProcess.get();
    system("pwd");
    system("ls ../data/");
    system("ls ../");
    // Load the data and compute the time steps from the test data
    // std::string testDensity  = file2String("../data/Density.csv");
    // std::string testEnergy   = file2String("../data/Energy.csv");
    // std::string testMomentum = file2String("../data/Momentum.csv");
    // size_t numTestTimesteps  = std::count(testDensity.begin(), testDensity.end(), '\n') - 1;
    // std::cout << "Finished test data import" << std::endl;

    // Now let's do the actual testings
    // This could be broken up into two tests. One for time steps and one for
    // content equality but I think it would be better to keep them as is and
    // stop the string equality tests from running if the time step equality fails
    // =========================================================================
    // First assert that the number of time steps is equal. If not there's no
    // point in running the string comparison
    // ASSERT_EQ(numTestTimesteps, numFidTimesteps) << "Time step equality failed";

    // EXPECT_EQ(testDensity,  fidDensity);
    // EXPECT_EQ(testEnergy,   fidEnergy);
    // EXPECT_EQ(testMomentum, fidMomentum);
}