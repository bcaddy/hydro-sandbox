/*!
 * \file ifdef-test.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Prototyping testing with if def statements
 * \version 0.1
 * \date 2021-07-29
 *
 * \copyright Copyright (c) 2021
 *
 */

// Include GoogleTest and related libraries/headers
#include <gtest/gtest.h>

// Local includes
#include "testingUtilities.h"
using namespace testingUtilities;

// Lets start testing.

// The tests that test a function that might not exist have to either be hidden
// with ifdef or the function they're testing has to have a null return option.
// Tests that might return different results depending on the combo if ifdefs
// also either have to be conditionally compiled OR somehow gather the ifdef
// values and then choose the correct value to compare against
#ifdef  FUNC_1
TEST(ifdef, testsFuncOne)
{

    // EXPECT_EQ does a bitwise comparison. For FP numbers this should be
    // equivalent to EXPECT_NEAR with the margin set to zero
    EXPECT_EQ(ifdefTester(), 1);
}
#endif  //FUNC_1

#ifdef  FUNC_2
TEST(ifdef, testsFuncTwo)
{

    // EXPECT_EQ does a bitwise comparison. For FP numbers this should be
    // equivalent to EXPECT_NEAR with the margin set to zero
    EXPECT_EQ(ifdefTester(), 2);
}
#endif  //FUNC_2