#pragma once

#ifdef HIP_BUILD
    #include <hip/hip_runtime.h>

    #define cudaDeviceSynchronize hipDeviceSynchronize
    #define cudaError hipError_t
    #define cudaError_t hipError_t
    #define cudaErrorInsufficientDriver hipErrorInsufficientDriver
    #define cudaErrorNoDevice hipErrorNoDevice
    #define cudaEvent_t hipEvent_t
    #define cudaEventCreate hipEventCreate
    #define cudaEventElapsedTime hipEventElapsedTime
    #define cudaEventRecord hipEventRecord
    #define cudaEventSynchronize hipEventSynchronize
    #define cudaFree hipFree
    #define cudaFreeHost hipHostFree
    #define cudaGetDevice hipGetDevice
    #define cudaGetDeviceCount hipGetDeviceCount
    #define cudaGetErrorString hipGetErrorString
    #define cudaGetLastError hipGetLastError
    #define cudaHostAlloc hipHostMalloc
    #define cudaHostAllocDefault hipHostMallocDefault
    #define cudaMalloc hipMalloc
    #define cudaMemcpy hipMemcpy
    #define cudaMemcpyAsync hipMemcpyAsync
    #define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
    #define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
    #define cudaMemcpyHostToDevice hipMemcpyHostToDevice
    #define cudaMemGetInfo hipMemGetInfo
    #define cudaMemset hipMemset
    #define cudaReadModeElementType hipReadModeElementType
    #define cudaSetDevice hipSetDevice
    #define cudaSuccess hipSuccess

    // New stuff
    #define cudaDeviceProp hipDeviceProp_t
    #define cudaGetDeviceProperties hipGetDeviceProperties

#endif  // HIP_BUILD


#ifdef CUDA_BUILD
    // #define hipLaunchKernelGGL(F,G,B,M,S,...) F<<<G,B,M,S>>>(__VA_ARGS__)  // FIXME not new
    #define __shfl_down(...) __shfl_down_sync(0xFFFFFFFF, __VA_ARGS__)
#endif  // CUDA_BUILD