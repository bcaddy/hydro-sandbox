#include <limits>

#pragma once

// Done editing for CUDA/HIP
__inline__ __device__ int warpReduceMax(int val)
{
    for (int offset = warpSize/2; offset > 0; offset /= 2)
    {
        val = max(val, __shfl_down_sync(0xFFFFFFFF,val, offset));
    }
    return val;
}

__inline__ __device__ int blockReduceMax(int val)
{
    __shared__ int shared[::maxWarpsPerBlock]; // Shared memory for storing the results of each warp-wise partial reduction

    int lane   = threadIdx.x % warpSize;  // thread ID within the warp,
    int warpId = threadIdx.x / warpSize;  // ID of the warp itself

    val = warpReduceMax(val);     // Each warp performs partial reduction

    if (lane==0) shared[warpId]=val; // Write reduced value to shared memory

    __syncthreads();              // Wait for all partial reductions

    //read from shared memory only if that warp existed
    val = (threadIdx.x < ::maxWarpsPerBlock) ? shared[lane] : 0;

    if (warpId==0) val = warpReduceMax(val); //Final reduce within first warp

    return val;
}

/*!
 * \brief Find the max of the array
 *
 * \param[in] in The input array
 * \param[out] out The output scalar variable
 * \param[in] N the number of elements in the `in` array
 * \return __global__
 */
__global__ void deviceReduceBlockAtomicKernelMax(int *in, int* out, int N)
{
    // Initialize variable to store the max value
    int maxVal = INT_MIN;
    // -DBL_MAX from float.h is what we want eventually

    // Grid stride loop to perform as much of the reduction as possible
    for(int i = blockIdx.x * blockDim.x + threadIdx.x; i < N; i += blockDim.x * gridDim.x)
    {
        // A transformation could go here

        // Grid stride reduction
        maxVal = max(maxVal,in[i]);
    }

    // Reduce the entire block in parallel
    maxVal = blockReduceMax(maxVal);

    // Write block level reduced value to the output scalar atomically
    if (threadIdx.x == 0) atomicMax(out, maxVal);
}

int gpuMaxReduction()
{
    int const numBlocks       = 2;
    int const threadsPerBlock = 1024;

    size_t const size = 256*256*256;
    int const maxValue = 3;

    int host_max;
    int *dev_vec, *dev_max;
    std::vector<int> host_vec(size, 1);
    host_vec.at(137) = maxValue;

    cudaMalloc(&dev_vec, host_vec.size() * sizeof(int));
    cudaMalloc(&dev_max, sizeof(int));
    cudaMemcpy(dev_vec, host_vec.data(), host_vec.size() * sizeof(int), cudaMemcpyHostToDevice);

    deviceReduceBlockAtomicKernelMax<<<numBlocks, threadsPerBlock>>>(dev_vec, dev_max, host_vec.size());

    cudaMemcpy(&host_max, dev_max, sizeof(int), cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();

    std::string result = (host_max == maxValue)? "correct": "incorrect";
    std::cout << std::endl << "The final result should be " << maxValue << " and is " << host_max << " which is " << result << std::endl;

    return host_max;
}
