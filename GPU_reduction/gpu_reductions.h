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
#pragma once

// External Includes
#include <cmath>

// Local Includes
#include "../global/global.h"
#include "../global/global_cuda.h"
#include "../utils/gpu.hpp"

/*!
 * \brief This namespace contains functions that perform various GPU resident
 * reductions. The algorithm is taken and modified from
 * [Mark Harris](https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf)
 * and
 * [Sean Bone](http://seanbone.ch/cuda-efficient-parallel-reduction/)
 *
 */
namespace gpuReductions
{
    // PRE:
    // dA is an array allocated on the GPU
    // N <= len(dA) is a power of two (N >= BLOCKSIZE)
    // POST: the sum of the first N elements of dA is returned


    // =========================================================================
    // GPU resident reduction to find the maximum of an array
    // =========================================================================
    /*!
     * \brief Performs a reduction within a single warp. Called by the
       gpuReductions::reduceMaxKernel
     *
     * \tparam T Any numerical type
     * \param[in,out] sdata The __shared__ array to find the maximum of
     * \param[in] tid The thread ID
     */
    template <typename T>
    __device__ void warpReduceMax(volatile T *sdata, size_t tid)
    {
        if (TPB >= 64) sdata[tid] = max(sdata[tid], sdata[tid + 32]);
        if (TPB >= 32) sdata[tid] = max(sdata[tid], sdata[tid + 16]);
        if (TPB >= 16) sdata[tid] = max(sdata[tid], sdata[tid +  8]);
        if (TPB >=  8) sdata[tid] = max(sdata[tid], sdata[tid +  4]);
        if (TPB >=  4) sdata[tid] = max(sdata[tid], sdata[tid +  2]);
        if (TPB >=  2) sdata[tid] = max(sdata[tid], sdata[tid +  1]);
    }

    /*!
     * \brief Kernel to compute the maximum of an array. Note that this probably
       shouldn't be called on its own. Instead call gpuReductions::reduceMaxGPU
     *
     * \tparam T Any numerical type
     * \param[in] inputArray Input array to find maximum of
     * \param[out] outputArray Output array to store the maximums that each
       thread finds
     * \param[in] n How many elements of the input array to search in
     */
    template <typename T>
    __global__ void reduceMaxKernel(T* inputArray, T* outputArray, size_t n)
    {
        __shared__ T sdata[TPB];

        size_t tid = threadIdx.x;
        size_t i = blockIdx.x*(TPB) + tid;
        size_t gridSize = TPB*gridDim.x;
        sdata[tid] = -1E30;

        while (i < n) { sdata[tid] = max(sdata[tid],inputArray[i]); i += gridSize; }

        __syncthreads();

        if (TPB >= 1024) { if (tid < 512) { sdata[tid] = max(sdata[tid],sdata[tid + 512]); } __syncthreads(); }
        if (TPB >=  512) { if (tid < 256) { sdata[tid] = max(sdata[tid],sdata[tid + 256]); } __syncthreads(); }
        if (TPB >=  256) { if (tid < 128) { sdata[tid] = max(sdata[tid],sdata[tid + 128]); } __syncthreads(); }
        if (TPB >=  128) { if (tid <  64) { sdata[tid] = max(sdata[tid],sdata[tid +  64]); } __syncthreads(); }

        if (tid < 32) warpReduceMax(sdata, tid);
        if (tid == 0) outputArray[blockIdx.x] = sdata[0];
    }

    /*!
     * \brief Determine the maximum of an array using an entirely GPU resident
       reduction algorithm. Does not require that the array be a power of 2 in
       length (I think)
     *
     * \tparam T Any numerical type
     * \param inputArray The pointer to the device array to be reduced
     * \param N The size of inputArray
     * \return T The maximum value in that array
     */
    template<typename T>
    T reduceMaxGPU(T* inputArray, size_t N)
    {
        // Set initial values
        T max = 0.;
        size_t n = N;
        size_t blocksPerGrid = std::ceil((1.*n) / TPB);

        // Allocate temporary array used for sorting
        T* tmp;
        CudaSafeCall(cudaMalloc(&tmp, sizeof(T) * blocksPerGrid));

        // Create a new pointer that we can repoint to avoid overwriting
        // inputArray
        T* from = inputArray;

        // Launch reductions in a loop. Each iteration reduces the number of
        // elements by a factor of 1/TPB.
        do
        {

            blocksPerGrid   = std::ceil((1.*n) / TPB);
            hipLaunchKernelGGL(reduceMaxKernel, blocksPerGrid, TPB, 0, 0, from, tmp, n);
            from = tmp;
            n = blocksPerGrid;
        } while (n > TPB);

        if (n > 1)
            hipLaunchKernelGGL(reduceMaxKernel, 1, TPB, 0, 0, tmp, tmp, n);

        // Sync now that reduction is done
        cudaDeviceSynchronize();

        // Copy the reduced value from the host
        CudaSafeCall(cudaMemcpy(&max, tmp, sizeof(T), cudaMemcpyDeviceToHost));

        // Free device memory
        CudaSafeCall(cudaFree(tmp));

        return max;
    }
    // =========================================================================
    // End of GPU resident reduction to find the maximum of an array
    // =========================================================================

    // =========================================================================
    // GPU resident reduction to find the minimum of an array
    // =========================================================================
    /*!
     * \brief Performs a reduction within a single warp. Called by the
       gpuReductions::reduceMinKernel
     *
     * \tparam T Any numerical type
     * \param[in,out] sdata The __shared__ array to find the minimum of
     * \param[in] tid The thread ID
     */
    template <typename T>
    __device__ void warpReduceMin(volatile T *sdata, size_t tid)
    {
        if (TPB >= 64) sdata[tid] = min(sdata[tid], sdata[tid + 32]);
        if (TPB >= 32) sdata[tid] = min(sdata[tid], sdata[tid + 16]);
        if (TPB >= 16) sdata[tid] = min(sdata[tid], sdata[tid +  8]);
        if (TPB >=  8) sdata[tid] = min(sdata[tid], sdata[tid +  4]);
        if (TPB >=  4) sdata[tid] = min(sdata[tid], sdata[tid +  2]);
        if (TPB >=  2) sdata[tid] = min(sdata[tid], sdata[tid +  1]);
    }

    /*!
     * \brief Kernel to compute the minimum of an array. Note that this probably
       shouldn't be called on its own. Instead call gpuReductions::reduceMinGPU
     *
     * \tparam T Any numerical type
     * \param[in] inputArray Input array to find minimum of
     * \param[out] outputArray Output array to store the minimums that each
       thread finds
     * \param[in] n How many elements of the input array to search in
     */
    template <typename T>
    __global__ void reduceMinKernel(T* inputArray, T* outputArray, size_t n)
    {
        __shared__ T sdata[TPB];

        size_t tid = threadIdx.x;
        size_t i = blockIdx.x*(TPB) + tid;
        size_t gridSize = TPB*gridDim.x;
        sdata[tid] = 1E30;

        while (i < n) { sdata[tid] = min(sdata[tid],inputArray[i]); i += gridSize; }

        __syncthreads();

        if (TPB >= 1024) { if (tid < 512) { sdata[tid] = min(sdata[tid],sdata[tid + 512]); } __syncthreads(); }
        if (TPB >=  512) { if (tid < 256) { sdata[tid] = min(sdata[tid],sdata[tid + 256]); } __syncthreads(); }
        if (TPB >=  256) { if (tid < 128) { sdata[tid] = min(sdata[tid],sdata[tid + 128]); } __syncthreads(); }
        if (TPB >=  128) { if (tid <  64) { sdata[tid] = min(sdata[tid],sdata[tid +  64]); } __syncthreads(); }

        if (tid < 32) warpReduceMin(sdata, tid);
        if (tid == 0) outputArray[blockIdx.x] = sdata[0];
    }

    /*!
     * \brief Determine the minimum of an array using an entirely GPU resident
       reduction algorithm. Does not require that the array be a power of 2 in
       length (I think)
     *
     * \tparam T Any numerical type
     * \param inputArray The pointer to the device array to be reduced
     * \param N The size of inputArray
     * \return T The minimum value in that array
     */
    template<typename T>
    T reduceMinGPU(T* inputArray, size_t N)
    {
        T min = 0.;
        size_t n = N;
        size_t blocksPerGrid = std::ceil((1.*n) / TPB);

        T* tmp;
        CudaSafeCall(cudaMalloc(&tmp, sizeof(T) * blocksPerGrid));

        T* from = inputArray;

        do
        {
            blocksPerGrid   = std::ceil((1.*n) / TPB);
            hipLaunchKernelGGL(reduceMinKernel, blocksPerGrid, TPB, 0, 0, from, tmp, n);
            from = tmp;
            n = blocksPerGrid;
        } while (n > TPB);

        if (n > 1)
            hipLaunchKernelGGL(reduceMinKernel, 1, TPB, 0, 0, tmp, tmp, n);

        cudaDeviceSynchronize();

        CudaSafeCall(cudaMemcpy(&min, tmp, sizeof(T), cudaMemcpyDeviceToHost));
        CudaSafeCall(cudaFree(tmp));
        return min;
    }
    // =========================================================================
    // End of GPU resident reduction to find the minimum of an array
    // =========================================================================
}