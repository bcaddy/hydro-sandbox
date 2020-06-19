#pragma once

#include <vector>
#include <fstream>

void saveArray(const std::vector<double> &arr,
               std::ofstream &fileName,
               const int &numGhosts);
