/*!
* \file mhdUtilities-tests.cpp
* \author Robert 'Bob' Caddy (rvc@pitt.edu)
* \brief Test the functions within mhd1D/mhdoUtilities.h
* \version 0.1
* \date 2021-07-27
*
* \copyright Copyright (c) 2021
*
*/

// Include GoogleTest and related libraries/headers
#include <gtest/gtest.h>

// Include code to test
#include "../mhd1D/mhdUtilities.h"
using namespace mhdUtilities;

// Lets start testing
TEST(AlfvenSpeedTest, HandlesProperInput)
{
    // EXPECT_DOUBLE_EQ and EXPECT_FLOAT_EQ compare the floats and return a pass
    // if they're within 4 ULPs. A ULP is found by converting the floating point
    // numbers to integers then subtracting. The resulting number is the ULP
    // (Units in the Last Place)
    // https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
    // Correct Value is 93435318430231351
    EXPECT_DOUBLE_EQ(alfvenSpeed(1.1, 2.1), 0.75907211527658969);

    // EXPECT_NEAR(val1, val2, absolute error) returns a pass if the two values
    // are within the absolute error margin. That margin can be zero
    EXPECT_NEAR(alfvenSpeed(1.1, 2.1), 0.75907211527658969,
                                       0.0);
}

TEST(AlfvenSpeedTest, HandlesNegativeValues)
{
    // Check that the code fails properly when given negative values
    // This test just tests that the proper type of exception is thrown, not
    // that the message is also correct.
    EXPECT_THROW(alfvenSpeed(1,-1), std::runtime_error);

    // If you want to check the message as well you can do this
    // this tests _that_ the expected exception is thrown
    EXPECT_THROW({
        try
        {
            alfvenSpeed(1,-1);
        }
        catch( const std::runtime_error& errMessage )
        {
            // and this tests that it has the correct message
            EXPECT_STREQ("Complex valued Alfven velocity detected. Exiting.", errMessage.what() );
            throw;
        }
    }, std::runtime_error );
}
