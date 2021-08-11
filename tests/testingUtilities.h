/*!
 * \file testUtilities.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Declares some basic utility functions commonly used in
 * testing.
 * \version 0.1
 * \date 2021-06-29
 *
 * \copyright Copyright (c) 2020
 *
 */
#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

/*!
 * \brief A namespace for common functions used in testing
 *
 */
namespace testingUtilities
{
// =============================================================================
/*!
 * \brief Reads an entire plain text file into a single std::string and returns
 * it
 *
 * \param filename[in] The path and name of the file to be read
 * \return std::string The contents of the file as a string
 */
std::string file2String(const char *filename);
// =============================================================================

// =============================================================================
// Functions for testing ifdef statements in Cholla
#ifdef  FUNC_1
int ifdefTester();
#endif  //FUNC_1

#ifdef  FUNC_2
int ifdefTester();
#endif  //FUNC_2
// =============================================================================
}