#include <limits>
#include <float.h>
#include <random>
#include "PerfTimer.h"
#include "gpu.hpp"

#pragma once


__inline__ __device__ Real warpReduceMax(Real val)
{
    for (int offset = warpSize/2; offset > 0; offset /= 2)
    {
        val = max(val, __shfl_down(val, offset));
    }
    return val;
}

__inline__ __device__ Real blockReduceMax(Real val)
{
    __shared__ Real shared[::maxWarpsPerBlock]; // Shared memory for storing the results of each warp-wise partial reduction

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

__device__ double atomicMax_double(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*) address;
    unsigned long long int old = *address_as_ull, assumed;
    // Explanation of loop here:
    // https://stackoverflow.com/questions/16077464/atomicadd-for-double-on-gpu
    // The loop is to make sure the value at address doesn't change between the
    // load at the atomic since the entire operation isn't atomic

    // While it appears that this could result in many times more atomic
    // operations than required, in practice it's only a handful of extra
    // operation even in the worst case. Running with 16,000 blocks gives ~8-37
    // atomics after brief testing
    do {
        assumed = old;
        old = atomicCAS(address_as_ull,
                        assumed,
                        __double_as_longlong(fmax(val, __longlong_as_double(assumed))));
    } while (assumed != old);
    return __longlong_as_double(old);
}

/*!
 * \brief Find the max of the array
 *
 * \param[in] in The input array
 * \param[out] out The output scalar variable
 * \param[in] N the number of elements in the `in` array
 * \return __global__
 */
__global__ void deviceReduceAtomicMax(Real *in, Real* out, int N)
{
    // Initialize variable to store the max value
    Real maxVal = -DBL_MAX;

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
    if (threadIdx.x == 0) atomicMax_double(out, maxVal);
}

/*!
 * \brief Find the max of the array
 *
 * \param[in] in The input array
 * \param[out] out The output scalar variable
 * \param[in] N the number of elements in the `in` array
 * \return __global__
 */
 __global__ void deviceReduceMax(Real *in, Real* out, int N)
 {
     // Initialize variable to store the max value
     Real maxVal = -DBL_MAX;

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
     if (threadIdx.x == 0)
     {
        out[blockIdx.x] = maxVal;
     }
 }

Real gpuAtomicMaxReduction(int numTrials = 100)
{
    // Launch parameters
    // =================
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);

    // Divide the total number of allowed threads by the number of threads per block
    int const numBlocks  = (prop.maxThreadsPerMultiProcessor * prop.multiProcessorCount) / prop.maxThreadsPerBlock;
    int const threadsPerBlock = prop.maxThreadsPerBlock;

    // Grid Parameters & testing parameters
    // ====================================
    size_t const size     = std::pow(512, 3);;
    Real   const maxValue = 4;
    std::vector<Real> host_grid(size);
    Real host_max;

    int const warmUps = 5;
    PerfTimer atomicTimer("AtomicMax Reduction Timer");

    // Fill grid with random values and randomly assign maximum value
    std::random_device rd;
    std::mt19937 prng(rd());
    std::uniform_real_distribution<double> doubleRand(-std::abs(maxValue)-1, std::abs(maxValue) - 1);
    std::uniform_int_distribution<int> intRand(0, host_grid.size()-1);
    for (size_t i = 0; i < host_grid.size(); i++)
    {
        host_grid.at(i) = doubleRand(prng);
    }
    host_grid.at(intRand(prng)) = maxValue;


    // Allocating and copying to device
    // ================================
    Real *dev_grid, *dev_max;
    cudaMalloc(&dev_grid, host_grid.size() * sizeof(Real));
    cudaMalloc(&dev_max, sizeof(Real));
    cudaMemcpy(dev_grid, host_grid.data(), host_grid.size() * sizeof(Real), cudaMemcpyHostToDevice);

    for (size_t trial = 0; trial < numTrials + warmUps; trial++)
    {
        if (trial >= warmUps)
        {
            atomicTimer.startTimer();
        }
        // Do the reduction
        // ================
        hipLaunchKernelGGL(deviceReduceAtomicMax, numBlocks, threadsPerBlock, 0, 0, dev_grid, dev_max, host_grid.size());

        // Copy back and sync
        // ==================
        cudaMemcpy(&host_max, dev_max, sizeof(Real), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();

        if (trial >= warmUps)
        {
            atomicTimer.stopTimer();
        }

        // Check Results
        // =============
        if (host_max != maxValue)
        {
            std::cout << "The final result should be " << maxValue
                      << " but is " << host_max << " which is incorrect." << std::endl;
        }
    }

    // Report Performance Results
    // ==========================
    atomicTimer.reportStats();

    return host_max;
}

Real gpuMaxReduction(int numTrials = 100)
{
    // Launch parameters
    // =================
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);

    // Divide the total number of allowed threads by the number of threads per block
    int const numBlocks  = (prop.maxThreadsPerMultiProcessor * prop.multiProcessorCount) / prop.maxThreadsPerBlock;
    int const threadsPerBlock = prop.maxThreadsPerBlock;

    // Grid Parameters & testing parameters
    // ====================================
    size_t const size     = std::pow(512, 3);;
    Real   const maxValue = 4;
    std::vector<Real> host_grid(size);
    Real host_max;

    int const warmUps = 5;
    PerfTimer atomicTimer("Max Reduction Timer");

    // Fill grid with random values and randomly assign maximum value
    std::random_device rd;
    std::mt19937 prng(rd());
    std::uniform_real_distribution<double> doubleRand(-std::abs(maxValue)-1, std::abs(maxValue) - 1);
    std::uniform_int_distribution<int> intRand(0, host_grid.size()-1);
    for (size_t i = 0; i < host_grid.size(); i++)
    {
        host_grid.at(i) = doubleRand(prng);
    }
    host_grid.at(intRand(prng)) = maxValue;


    // Allocating and copying to device
    // ================================
    Real *dev_grid, *dev_max;
    cudaMalloc(&dev_grid, host_grid.size() * sizeof(Real));
    cudaMalloc(&dev_max, numBlocks*sizeof(Real));
    cudaMemcpy(dev_grid, host_grid.data(), host_grid.size() * sizeof(Real), cudaMemcpyHostToDevice);

    for (size_t trial = 0; trial < numTrials + warmUps; trial++)
    {
        if (trial >= warmUps)
        {
            atomicTimer.startTimer();
        }
        // Do the reduction
        // ================
        hipLaunchKernelGGL(deviceReduceMax, numBlocks, threadsPerBlock, 0, 0, dev_grid, dev_max, host_grid.size());
        hipLaunchKernelGGL(deviceReduceMax, 1,         threadsPerBlock, 0, 0, dev_max,  dev_max, threadsPerBlock);

        // Copy back and sync
        // ==================
        cudaMemcpy(&host_max, dev_max, sizeof(Real), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();

        if (trial >= warmUps)
        {
            atomicTimer.stopTimer();
        }

        // Check Results
        // =============
        if (host_max != maxValue)
        {
            std::cout << "The final result should be " << maxValue
                      << " but is " << host_max << " which is incorrect." << std::endl;
        }
    }

    // Report Performance Results
    // ==========================
    atomicTimer.reportStats();

    return host_max;
}
