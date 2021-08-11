/*!
 * \file testUtilities.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Implements some basic utility functions commonly used in
 * testing.
 * \version 0.1
 * \date 2021-06-29
 *
 * \copyright Copyright (c) 2020
 *
 */

#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace testingUtilities
{
// =============================================================================
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

// =============================================================================
// Functions for testing ifdef statements in Cholla
#ifdef  FUNC_1
int ifdefTester()
{
  return 1;
}
#endif  //FUNC_1

#ifdef  FUNC_2
int ifdefTester()
{
  return 2;
}
#endif  //FUNC_2
// =============================================================================
}