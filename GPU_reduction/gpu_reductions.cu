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

 // External Includes

// STL Includes
#include <stdio.h>
#include <cstdint>
#include <iostream>
#include <vector>
#include <string>

// Local Includes
#include "sumReduction.h"

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
    int sumReduced = gpuSumReduction();
    return 0;
}
