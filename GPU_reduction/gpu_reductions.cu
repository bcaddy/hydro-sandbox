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

#ifdef  CUDA_BUILD
    #define WARPSIZE 32
#endif  //CUDA_BUILD
#ifdef  HIP_BUILD
    #define WARPSIZE 64
#endif  //HIP_BUILD

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
// #include "sumReduction.h"
#include "maxReduction.h"
#include "timeStepTestOriginal.h"
#include "timeStepTestNew.h"

__global__ void checkDims()
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx == 0)
    {
        int const tid  = threadIdx.x;
        int const bid  = blockIdx.x;
        int const bdim = blockDim.x;
        int const gdim = gridDim.x;

        printf("\nthreadIdx.x = %i", tid);
        printf("\nblockIdx.x  = %i", bid);
        printf("\nblockDim.x  = %i", bdim);
        printf("\ngridDim.x   = %i", gdim);
    }
}

int	main()
{
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);

    std::cout << "prop.maxThreadsPerMultiProcessor = " << prop.maxThreadsPerMultiProcessor << std::endl;
    std::cout << "prop.multiProcessorCount         = " << prop.multiProcessorCount << std::endl;
    std::cout << "prop.maxThreadsPerBlock          = " << prop.maxThreadsPerBlock << std::endl;

    int numThreads = prop.maxThreadsPerMultiProcessor * prop.multiProcessorCount;
    int numBlocks  = numThreads / prop.maxThreadsPerBlock;

    std::cout << std::endl;
    std::cout << "numThreads = " << numThreads << std::endl;
    std::cout << "numBlocks  = " << numBlocks << std::endl;


    int trials = 1000;

    std::cout << std::endl;
    Real maxReducedAtomic = gpuAtomicMaxReduction(trials);

    std::cout << std::endl;
    Real maxReduced = gpuMaxReduction(trials);

    std::cout << std::endl;
    calcDtiOriginal(trials);

    std::cout << std::endl;
    calcDtiNEW(trials);

    // Mark variables as used to avoid warnings
    (void)maxReducedAtomic;
    (void)maxReduced;

    return 0;
}
