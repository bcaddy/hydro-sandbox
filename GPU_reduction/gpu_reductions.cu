/*!
 * \file gpu_reductions.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief
 * \version 0.1
 * \date 2022-03-04
 *
 * \copyright Copyright (c) 2022
 *
 */

// Globals
#define WARPSIZE 32
static constexpr int maxWarpsPerBlock = 1024/WARPSIZE; // outside kernel
typedef double Real;

 // External Includes

// STL Includes
#include <stdio.h>
#include <cstdint>
#include <iostream>
#include <vector>
#include <string>

// Local Includes
#include "sumReduction.h"
#include "maxReduction.h"

__global__ void checkDims()
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx == 0)
    {
        printf("\nthreadIdx.x = %i", threadIdx.x);
        printf("\nblockIdx.x  = %i", blockIdx.x);
        printf("\nblockDim.x  = %i", blockDim.x);
        printf("\ngridDim.x   = %i", gridDim.x);
    }
}

int	main()
{
    // cudaDeviceProp prop;
    // cudaGetDeviceProperties(&prop, 0);

    // std::cout << "prop.maxThreadsPerMultiProcessor = " << prop.maxThreadsPerMultiProcessor << std::endl;
    // std::cout << "prop.multiProcessorCount         = " << prop.multiProcessorCount << std::endl;
    // std::cout << "prop.maxThreadsPerBlock          = " << prop.maxThreadsPerBlock << std::endl;

    // int numThreads = prop.maxThreadsPerMultiProcessor * prop.multiProcessorCount;
    // int numBlocks  = numThreads / prop.maxThreadsPerBlock;

    // std::cout << std::endl;
    // std::cout << "numThreads = " << numThreads << std::endl;
    // std::cout << "numBlocks  = " << numBlocks << std::endl;


    // int sumReduced = gpuSumReduction();

    Real maxReduced = gpuAtomicMaxReduction();
    return 0;
}
