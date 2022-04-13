#pragma once
#include "gpu.hpp"

// =============================================================================
// NOTE: This sum reduction is not modified to work properly with HIP and
// Cholla. It should work but you would be better off basing a new sum reduction
// off of the max reduction
// =============================================================================

// Done editing for CUDA/HIP
__inline__ __device__ int warpReduceSum(int val)
{
    for (int offset = warpSize/2; offset > 0; offset /= 2)
    {
        val += __shfl_down(val, offset);
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
__global__ void deviceReduceBlockAtomicKernelSum(int *in, int* out, int N)
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

int gpuSumReduction()
{
    int const numBlocks       = 2;
    int const threadsPerBlock = 1024;

    size_t const size = 256*256*256;
    int host_sum;
    int *dev_vec, *dev_sum;
    std::vector<int> host_vec(size, 1);

    cudaMalloc(&dev_vec, host_vec.size() * sizeof(int));
    cudaMalloc(&dev_sum, sizeof(int));
    cudaMemcpy(dev_vec, host_vec.data(), host_vec.size() * sizeof(int), cudaMemcpyHostToDevice);

    hipLaunchKernelGGL(deviceReduceBlockAtomicKernelSum, numBlocks, threadsPerBlock, 0, 0, dev_vec, dev_sum, host_vec.size());

    cudaMemcpy(&host_sum, dev_sum, sizeof(int), cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();

    std::string result = (host_sum == host_vec.size())? "correct": "incorrect";
    std::cout << std::endl << "The final result should be " << host_vec.size() << " and is " << host_sum << " which is " << result << std::endl;

    return host_sum;
}