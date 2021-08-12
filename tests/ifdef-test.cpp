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

// Lets start testing
TEST(ifdef, testsFuncOne)
{

    // EXPECT_EQ does a bitwise comparison. For FP numbers this should be
    // equivalent to EXPECT_NEAR with the margin set to zero
    EXPECT_EQ(ifdefTester(), 1);
}

TEST(ifdef, testsFuncTwo)
{

    // EXPECT_EQ does a bitwise comparison. For FP numbers this should be
    // equivalent to EXPECT_NEAR with the margin set to zero
    EXPECT_EQ(ifdefTester(), 2);
}