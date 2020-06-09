#pragma once
#include <string>

double VelInterface(const double &a,
                    const double &b,
                    const double &c,
                    const double &d,
                    const double &deltax,
                    const double &deltat,
                    const std::string &LimType);

double AdvectionInterface(const double &a,
                          const double &b,
                          const double &c,
                          const double &d,
                          const double &vel,
                          const double &deltax,
                          const double &deltat,
                          const std::string &LimType);