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
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

// Include GoogleTest and related libraries/headers
#include <gtest/gtest.h>


std::string file2String(const char *filename)
{
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in)
  {
    std::ostringstream contents;
    contents << in.rdbuf();
    in.close();
    return(contents.str());
  }
  throw std::invalid_argument("File not found");
}

// Lets start testing
TEST(SystemTest_VanLeer, sod_HLLC_1k_cells)
{
    // First let's start our cholla analogue running asynchronously. Note that
    // this dumps all console I/O to /dev/nul
    auto sodProcess = std::async(std::launch::async, // Launch operation asynchronously
                                 system,             // Choose the function to launch
                                 "../euler1D-VL/euler-solver.exe >/dev/null 2>&1"); // Args to send to "system" call

    // While that's running we're going to load the fiducial data and find the
    // number of fiducial time steps
    std::string fidDensity  = file2String("System-Test-Data/sod-VL-HLLC/Density.csv");
    std::string fidEnergy   = file2String("System-Test-Data/sod-VL-HLLC/Energy.csv");
    std::string fidMomentum = file2String("System-Test-Data/sod-VL-HLLC/Momentum.csv");

    sodProcess.get();  // Wait for sodProcess to finish
    EXPECT_EQ(1,1);
}