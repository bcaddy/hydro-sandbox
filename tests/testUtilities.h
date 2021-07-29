/*!
 * \file testUtilities.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains some basic utility functions commonly used in
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
namespace testUtilities
{
// =============================================================================
/*!
 * \brief Reads an entire plain text file into a single std::string and returns
 * it
 *
 * \param filename[in] The path and name of the file to be read
 * \return std::string The contents of the file as a string
 */
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
// =============================================================================
}