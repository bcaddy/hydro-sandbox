#pragma once

#include "maxReduction.h"

#define TPB 256

__device__ __host__ Real hydroInverseCrossingTimeNew(Real const &E, Real const &d, Real const &d_inv, Real const &vx, Real const &vy, Real const &vz, Real const &dx, Real const &dy, Real const &dz, Real const &gamma)
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


__global__ void Calc_dt_3D_new(Real *dev_conserved, int nx, int ny, int nz, int n_ghost, int n_fields, Real dx, Real dy, Real dz, Real *dev_dti, Real gamma)
{
    Real max_dti;

    Real d, d_inv, vx, vy, vz, E;
    int xid, yid, zid, n_cells;

    n_cells = nx*ny*nz;

    // Grid stride loop to perform as much of the reduction as possible
    for(int id = threadIdx.x + blockIdx.x * blockDim.x; id < n_cells; id += blockDim.x * gridDim.x)
    {
        // Compute the real indices
        zid = id / (nx*ny);
        yid = (id - zid*nx*ny) / nx;
        xid = id - zid*nx*ny - yid*nx;

        // threads corresponding to real cells do the calculation
        if (    xid > n_ghost-1
            and xid < nx-n_ghost
            and yid > n_ghost-1
            and yid < ny-n_ghost
            and zid > n_ghost-1
            and zid < nz-n_ghost)
        {
            // every thread collects the conserved variables it needs from global memory
            d  =  dev_conserved[            id];
            d_inv = 1.0 / d;
            vx =  dev_conserved[1*n_cells + id] * d_inv;
            vy =  dev_conserved[2*n_cells + id] * d_inv;
            vz =  dev_conserved[3*n_cells + id] * d_inv;
            E  = dev_conserved[4*n_cells + id];

            // Compute the maximum inverse crossing time in the cell
            max_dti = fmax(max_dti, hydroInverseCrossingTimeNew(E, d, d_inv, vx, vy, vz, dx, dy, dz, gamma));
        }
    }

    // Perform reduction across the entire grid
    gridReduceMax(max_dti, dev_dti);
}


Real Calc_dt_GPU_NEW(Real *dev_conserved, Real *dev_dti, int nx, int ny, int nz, int n_ghost, int n_fields, Real dx, Real dy, Real dz, Real gamma, int n_cells)
{
    // Set launch parameters
    int threadsPerBlock, numBlocks;
    reductionLaunchParams(numBlocks, threadsPerBlock);

    // compute dt and store in dev_dti_array
    hipLaunchKernelGGL(Calc_dt_3D, numBlocks, threadsPerBlock, 0, 0, dev_conserved, nx, ny, nz, n_ghost, n_fields, dx, dy, dz, dev_dti, gamma);

    // copy dev_dti_array to host_dti_array
    Real host_dti;
    cudaMemcpy(&host_dti, dev_dti, sizeof(Real), cudaMemcpyDeviceToHost);

    return host_dti;
}


void calcDtiNEW(int numTrials=100)
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

    std::vector<Real> host_grid(n_cells*n_fields);

    int const warmUps = 5;
    PerfTimer timer("DTI New Timer");

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

    Real *dev_grid, *dev_dti;
    cudaMalloc(&dev_grid, host_grid.size() * sizeof(Real));
    cudaMalloc(&dev_dti, sizeof(Real));
    cudaMemcpy(dev_grid, host_grid.data(), host_grid.size() * sizeof(Real), cudaMemcpyHostToDevice);

    for (size_t trial = 0; trial < numTrials + warmUps; trial++)
    {
        if (trial >= warmUps)
        {
            timer.startTimer();
        }
        // Do the reduction
        // ================
        dti = Calc_dt_GPU_NEW(dev_grid, dev_dti, nx, ny, nz, n_ghost, n_fields, dx, dy, dz, gamma, n_cells);

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
