/*!
 * \file gpu_reductions_tests.cu
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Tests the contents of gpu_reduction.h and gpu_reductions.cpp
 *
 */

// STL Includes
#include <vector>
#include <string>
#include <iostream>

// External Includes
#include <gtest/gtest.h>    // Include GoogleTest and related libraries/headers

// Local Includes
#include "../global/global.h"
#include "../global/global_cuda.h"
#include "../utils/testing_utilities.h"
#include "../utils/gpu_reductions.h"

// =============================================================================
TEST(tHYDROReduceMaxGPU,
     PowerOf2RealInputExpectCorrectOutput)
{
    // Test Parameters
    size_t size = TPB * 4;
    Real const maxVal = 123.4;
    std::vector<Real> testArray(size);

    // Populate vector
    for (size_t i = 0; i < size; i++)
    {
        testArray.at(i) = -1.0 * i;
    }

    // Assign maximum to random spot
    srand(time(NULL));
    size_t location = rand() % size;
    testArray.at(location) = maxVal;

    // Copy Data
    Real* device_testArray;
    CudaSafeCall(cudaMalloc(&device_testArray, sizeof(Real)*size));
    CudaSafeCall(cudaMemcpy(device_testArray, testArray.data(), sizeof(Real)*size, cudaMemcpyHostToDevice));

    // Run kernel
    Real reducedValue = gpuReductions::reduceMaxGPU(device_testArray, size);

    // Perform check
    testingUtilities::checkResults(maxVal,
        reducedValue,
        "maximum value. Stored at location = " + std::to_string(location));
}
// =============================================================================

// =============================================================================
TEST(tHYDROReduceMaxGPU,
     NonPowerOf2RealInputExpectCorrectOutput)
{
   // Test Parameters
   size_t size = (TPB * 4) + 7;
   Real const maxVal = 123.4;
   std::vector<Real> testArray(size);

   // Populate vector
   for (size_t i = 0; i < size; i++)
   {
       testArray.at(i) = -1.0 * i;
   }

   // Assign maximum to the last spot
   size_t location = size-1;
   testArray.at(location) = maxVal;

   // Copy Data
   Real* device_testArray;
   CudaSafeCall(cudaMalloc(&device_testArray, sizeof(Real)*size));
   CudaSafeCall(cudaMemcpy(device_testArray, testArray.data(), sizeof(Real)*size, cudaMemcpyHostToDevice));

   // Run kernel
   Real reducedValue = gpuReductions::reduceMaxGPU(device_testArray, size);

   // Perform check
   testingUtilities::checkResults(maxVal,
       reducedValue,
       "maximum value. Stored at location = " + std::to_string(location));
}
// =============================================================================

// =============================================================================
TEST(tHYDROReduceMinGPU,
     PowerOf2RealInputExpectCorrectOutput)
{
   // Test Parameters
   size_t size = TPB * 4;
   Real const minVal = -123.4;
   std::vector<Real> testArray(size);

   // Populate vector
   for (size_t i = 0; i < size; i++)
   {
       testArray.at(i) = 1.0 * i;
   }

   // Assign maximum to random spot
   srand(time(NULL));
   size_t location = rand() % size;
   testArray.at(location) = minVal;

   // Copy Data
   Real* device_testArray;
   CudaSafeCall(cudaMalloc(&device_testArray, sizeof(Real)*size));
   CudaSafeCall(cudaMemcpy(device_testArray, testArray.data(), sizeof(Real)*size, cudaMemcpyHostToDevice));

   // Run kernel
   Real reducedValue = gpuReductions::reduceMinGPU(device_testArray, size);

   // Perform check
   testingUtilities::checkResults(minVal,
       reducedValue,
       "minimum value. Stored at location = " + std::to_string(location));
}
// =============================================================================

// =============================================================================
TEST(tHYDROReduceMinGPU,
     NonPowerOf2RealInputExpectCorrectOutput)
{
  // Test Parameters
  size_t size = (TPB * 4) + 7;
  Real const minVal = -123.4;
  std::vector<Real> testArray(size);

  // Populate vector
  for (size_t i = 0; i < size; i++)
  {
      testArray.at(i) = 1.0 * i;
  }

  // Assign maximum to the last spot
  size_t location = size-1;
  testArray.at(location) = minVal;

  // Copy Data
  Real* device_testArray;
  CudaSafeCall(cudaMalloc(&device_testArray, sizeof(Real)*size));
  CudaSafeCall(cudaMemcpy(device_testArray, testArray.data(), sizeof(Real)*size, cudaMemcpyHostToDevice));

  // Run kernel
  Real reducedValue = gpuReductions::reduceMinGPU(device_testArray, size);

  // Perform check
  testingUtilities::checkResults(minVal,
      reducedValue,
      "minimum value. Stored at location = " + std::to_string(location));
}
// =============================================================================
