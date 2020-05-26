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