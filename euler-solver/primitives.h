#pragma once

#include <vector>
#include <cmath>

struct primitives
{
public:
    int size;  // Length of each array

    std::vector<double> velocity;
    std::vector<double> density;
    std::vector<double> pressure;
    std::vector<double> siEnergy; //Specific Internal Energy

    // // Get the primitive variables
    // std::vector<double> GetVelocityVec() { return velocity; };
    // std::vector<double> GetDensityVec () { return density;  };
    // std::vector<double> GetPressureVec() { return pressure; };
    // std::vector<double> GetsiEnergyVec() { return siEnergy; };
    // double GetVelocityElement(int const &i) { return &velocity[i]; };
    // double GetDensityElement (int const &i) { return density[i];  };
    // double GetPressureElement(int const &i) { return pressure[i]; };
    // double GetsiEnergyElement(int const &i) { return siEnergy[i]; };

    // // Get the conserved variables
    // double GetMomentumElement(int const &i) { return velocity[i] * density[i]; };
    // double GetTotalEnergyElement(int const &i) { return siEnergy[i] + 0.5 * density[i] * std::pow(velocity[i], 2); };

    // std::vector<double> GetMomentumVec(int const &size);
    // std::vector<double> GetTotalEnergyVec(int const &size);

    primitives(int const &length);
    ~primitives();
};

primitives::primitives(int const &length)
{
    size = length;
    velocity.reserve(size);
    density.reserve(size);
    pressure.reserve(size);
    siEnergy.reserve(size);
}

primitives::~primitives()
{
}
