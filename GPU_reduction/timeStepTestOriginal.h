#pragma once

#include "maxReduction.h"

#define TPB 256

__device__ __host__ Real hydroInverseCrossingTime(Real const &E, Real const &d, Real const &d_inv, Real const &vx, Real const &vy, Real const &vz, Real const &dx, Real const &dy, Real const &dz, Real const &gamma)
{
    // Compute pressure and sound speed
    Real P  = (E - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    Real cs = sqrt(d_inv * gamma * P);

    // Find maximum inverse crossing time in the cell (i.e. minimum crossing time)
    Real cellMaxInverseDt = fmax((fabs(vx)+cs)/dx, (fabs(vy)+cs)/dy);
    cellMaxInverseDt      = fmax(cellMaxInverseDt, (fabs(vz)+cs)/dz);
    cellMaxInverseDt      = fmax(cellMaxInverseDt, 0.0);

    return cellMaxInverseDt;
}


__global__ void Calc_dt_3D(Real *dev_conserved, int nx, int ny, int nz, int n_ghost, int n_fields, Real dx, Real dy, Real dz, Real *dti_array, Real gamma)
{
    __shared__ Real max_dti[TPB];

    Real d, d_inv, vx, vy, vz, E;
    int id, xid, yid, zid, n_cells;
    int tid;

    n_cells = nx*ny*nz;

    // get a global thread ID
    id = threadIdx.x + blockIdx.x * blockDim.x;
    zid = id / (nx*ny);
    yid = (id - zid*nx*ny) / nx;
    xid = id - zid*nx*ny - yid*nx;
    // and a thread id within the block
    tid = threadIdx.x;

    // set shared memory to 0
    max_dti[tid] = 0;
    __syncthreads();

    // threads corresponding to real cells do the calculation
    if (xid > n_ghost-1 && xid < nx-n_ghost && yid > n_ghost-1 && yid < ny-n_ghost && zid > n_ghost-1 && zid < nz-n_ghost)
    {
        // every thread collects the conserved variables it needs from global memory
        d  =  dev_conserved[            id];
        d_inv = 1.0 / d;
        vx =  dev_conserved[1*n_cells + id] * d_inv;
        vy =  dev_conserved[2*n_cells + id] * d_inv;
        vz =  dev_conserved[3*n_cells + id] * d_inv;
        E  = dev_conserved[4*n_cells + id];

        // Compute the maximum inverse crossing time in the cell
        max_dti[tid] = hydroInverseCrossingTime(E, d, d_inv, vx, vy, vz, dx, dy, dz, gamma);
    }
    __syncthreads();

    // do the reduction in shared memory (find the max inverse timestep in the block)
    for (unsigned int s=1; s<blockDim.x; s*=2) {
        if (tid % (2*s) == 0) {
            max_dti[tid] = fmax(max_dti[tid], max_dti[tid + s]);
        }
        __syncthreads();
    }

    // write the result for this block to global memory
    if (tid == 0) dti_array[blockIdx.x] = max_dti[0];
}


Real Calc_dt_GPU_ORIGINAL(Real *dev_conserved, Real *dev_dti_array, Real *host_dti_array, int nx, int ny, int nz, int n_ghost, int n_fields, Real dx, Real dy, Real dz, Real gamma, int n_cells)
{
    // set values for GPU kernels
    int const ngrid = (n_cells + TPB - 1) / TPB;
    // number of blocks per 1D grid
    dim3 dim1dGrid(ngrid, 1, 1);
    //  number of threads per 1D block
    dim3 dim1dBlock(TPB, 1, 1);


    // compute dt and store in dev_dti_array
    hipLaunchKernelGGL(Calc_dt_3D, dim1dGrid, dim1dBlock, 0, 0, dev_conserved, nx, ny, nz, n_ghost, n_fields, dx, dy, dz, dev_dti_array, gamma);

    // copy dev_dti_array to host_dti_array
    cudaMemcpy(host_dti_array, dev_dti_array, ngrid*sizeof(Real), cudaMemcpyDeviceToHost);

    Real max_dti = 0.0;
    for (int i=0; i<ngrid; i++) {
        max_dti = fmax(max_dti, host_dti_array[i]);
    }

    return max_dti;
}


void calcDtiOriginal(int numTrials=100)
{
    // Grid Parameters & testing parameters
    // ====================================
    size_t const nx = 512, ny = nx, nz = nx;
    size_t const n_ghost = 4;
    size_t const n_cells  = (nx+n_ghost)*(ny+n_ghost)*(nz+n_ghost);
    size_t const n_fields = 5;
    Real dx = 3, dy = dx, dz = dx;
    Real gamma = 5./3.;
    Real dti;
    int const ngrid = (n_cells + TPB - 1) / TPB;

    std::vector<Real> host_grid(n_cells*n_fields);

    int const warmUps = 5;
    PerfTimer timer("DTI Original Timer");

    // Fill grid with random values and randomly assign maximum value
    std::random_device rd;
    std::mt19937 prng(rd());
    std::uniform_real_distribution<double> doubleRand(1, 5);
    std::uniform_int_distribution<int> intRand(0, host_grid.size()-1);
    for (size_t i = 0; i < host_grid.size(); i++)
    {
        host_grid.at(i) = doubleRand(prng);
    }

    // Allocating and copying to device
    // ================================

    Real *dev_grid, *dev_max, *dev_dti_array, *host_dti_array;
    cudaMalloc(&dev_grid, host_grid.size() * sizeof(Real));
    cudaMalloc(&dev_max, sizeof(Real));
    cudaMalloc(&dev_dti_array, ngrid*sizeof(Real));
    cudaHostAlloc(&host_dti_array, ngrid*sizeof(Real));
    cudaMemcpy(dev_grid, host_grid.data(), host_grid.size() * sizeof(Real), cudaMemcpyHostToDevice);

    for (size_t trial = 0; trial < numTrials + warmUps; trial++)
    {
        if (trial >= warmUps)
        {
            timer.startTimer();
        }
        // Do the reduction
        // ================
        dti = Calc_dt_GPU_ORIGINAL(dev_grid, dev_dti_array, host_dti_array, nx, ny, nz, n_ghost, n_fields, dx, dy, dz, gamma, n_cells);

        cudaDeviceSynchronize();

        if (trial >= warmUps)
        {
            timer.stopTimer();
        }
    }

    // Report Performance Results
    // ==========================
    timer.reportStats();
}