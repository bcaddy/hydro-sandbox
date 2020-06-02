#pragma once

#include <vector>
#include <fstream>
#include <string>

void saveArray(const std::vector<double> &arr,
               std::ofstream &fileName,
               const int &numGhosts);

int setInitialConditions(std::vector<double> &arr, 
                          const int &size,
                          const std::string kind);

double minModLimiter(const double &a0,
                     const double &a1,
                     const double &a2,
                     const double &deltax);

double VelInterface(const double &a,
                    const double &b,
                    const double &c,
                    const double &d,
                    const double &deltax,
                    const double &deltat);