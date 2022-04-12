#include <limits>
#include <float.h>

#pragma once

// Done editing for CUDA/HIP
__inline__ __device__ Real warpReduceMax(Real val)
{
    for (int offset = warpSize/2; offset > 0; offset /= 2)
    {
        val = max(val, __shfl_down_sync(0xFFFFFFFF,val, offset));
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
__global__ void deviceReduceBlockAtomicKernelMax(Real *in, Real* out, int N)
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

Real gpuMaxReduction()
{
    int const numBlocks       = 100;
    int const threadsPerBlock = 1024;

    size_t const size = 256*256*256;
    Real const maxValue = 3;

    Real host_max;
    Real *dev_vec, *dev_max;
    std::vector<Real> host_vec(size, 1);
    host_vec.at(256*123*185) = maxValue;

    cudaMalloc(&dev_vec, host_vec.size() * sizeof(Real));
    cudaMalloc(&dev_max, sizeof(Real));
    cudaMemcpy(dev_vec, host_vec.data(), host_vec.size() * sizeof(Real), cudaMemcpyHostToDevice);

    deviceReduceBlockAtomicKernelMax<<<numBlocks, threadsPerBlock>>>(dev_vec, dev_max, host_vec.size());

    cudaMemcpy(&host_max, dev_max, sizeof(Real), cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();

    std::string result = (host_max == maxValue)? "correct": "incorrect";
    std::cout << std::endl << "The final result should be " << maxValue << " and is " << host_max << " which is " << result << std::endl;

    return host_max;
}
