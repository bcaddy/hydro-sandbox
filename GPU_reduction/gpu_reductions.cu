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

// External Includes

// STL Includes
#include <stdio.h>
#include <cstdint>

// Local Includes

// Done editing for CUDA/HIP
__inline__ __device__ int warpReduceSum(int val)
{
    #ifdef O_HIP
        uint64_t const mask = 0xFFFFFFFFFFFFFFFF;
    #else
        uint32_t const mask = 0xFFFFFFFF;
    #endif  //HIP

    for (int offset = warpSize/2; offset > 0; offset /= 2)
    {
        val += __shfl_down_sync(mask, val, offset);
    }

    return val;
}


__inline__ __device__ int blockReduceSum(int val)
{

    int const numWarps = blockDim.x/warpSize;

    __shared__ int shared[numWarps]; // Shared memory for storing the results of each warp-wise partial reduction

    // TODO STOPPED HERE LAST NIGHT
    int lane = threadIdx.x % warpSize;
    int warpId = threadIdx.x / warpSize;

    val = warpReduceSum(val);     // Each warp performs partial reduction

    if (lane==0) shared[warpId]=val; // Write reduced value to shared memory

    __syncthreads();              // Wait for all partial reductions

    //read from shared memory only if that warp existed
    val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;

    if (warpId==0) val = warpReduceSum(val); //Final reduce within first warp

    return val;
}


__global__ void deviceReduceBlockAtomicKernel(int *in, int* out, int N)
{
    // Initialize sum
    int sum = int(0);

    // Grid stride loop to perform as much of the reduction as possible
    for(int i = blockIdx.x * blockDim.x + threadIdx.x; i < N; i += blockDim.x * gridDim.x)
    {
        // A transformation could go here

        // Grid stride reduction
        sum += in[i];
    }

    // Reduce the entire block in parallel
    sum = blockReduceSum(sum);

    // Write block level reduced value to the output scalar atomically
    if (threadIdx.x == 0) atomicAdd(out, sum);
}


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
    int const numBlocks       = 2;
    int const threadsPerBlock = 1024;
    checkDims<<<numBlocks, threadsPerBlock>>>();
    cudaDeviceSynchronize();

    return 0;
}
