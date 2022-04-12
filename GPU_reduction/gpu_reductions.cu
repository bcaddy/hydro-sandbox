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
#include <iostream>
#include <vector>

// Local Includes

#define WARPSIZE 32
static constexpr int maxWarpsPerBlock = 1024/WARPSIZE; // outside kernel

// Done editing for CUDA/HIP
__inline__ __device__ int warpReduceSum(int val)
{
    for (int offset = warpSize/2; offset > 0; offset /= 2)
    {
        val += __shfl_down_sync(0xFFFFFFFF,val, offset);
    }
    return val;
}

__inline__ __device__ int blockReduceSum(int val)
{
    __shared__ int shared[::maxWarpsPerBlock]; // Shared memory for storing the results of each warp-wise partial reduction

    int lane   = threadIdx.x % warpSize;  // thread ID within the warp,
    int warpId = threadIdx.x / warpSize;  // ID of the warp itself

    val = warpReduceSum(val);     // Each warp performs partial reduction

    if (lane==0) shared[warpId]=val; // Write reduced value to shared memory

    __syncthreads();              // Wait for all partial reductions

    //read from shared memory only if that warp existed
    val = (threadIdx.x < ::maxWarpsPerBlock) ? shared[lane] : 0;

    if (warpId==0) val = warpReduceSum(val); //Final reduce within first warp

    return val;
}

/*!
 * \brief Find the sum of the array
 *
 * \param[in] in The input array
 * \param[out] out The output scalar variable
 * \param[in] N the number of elements in the `in` array
 * \return __global__
 */
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
    // checkDims<<<numBlocks, threadsPerBlock>>>();

    size_t const size = 256*256*256;
    int host_sum;
    int *dev_vec, *dev_sum;
    std::vector<int> host_vec(size, 1);

    cudaMalloc(&dev_vec, host_vec.size() * sizeof(int));
    cudaMalloc(&dev_sum, sizeof(int));
    cudaMemcpy(dev_vec, host_vec.data(), host_vec.size() * sizeof(int), cudaMemcpyHostToDevice);

    deviceReduceBlockAtomicKernel<<<numBlocks, threadsPerBlock>>>(dev_vec, dev_sum, host_vec.size());

    cudaMemcpy(&host_sum, dev_sum, sizeof(int), cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();


    std::cout << "The final result should be " << host_vec.size() << " and is " << host_sum << std::endl;

    return 0;
}
