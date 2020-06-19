#include <iostream>

#include "helper.h"

using std::cout;
using std::endl;

// =============================================================================
void saveArray(const std::vector<double> &arr,
               std::ofstream &fileName,
               const int &numGhosts)
{
    int size = arr.size();

    if (fileName.is_open())
    {
        fileName << arr[numGhosts];
        for (int i = numGhosts+1; i < (size - numGhosts); i++)
        {
            fileName << "," << arr[i];
        }
        fileName << endl;
    }
    else
    {
        cout << "File not open";
    }
}
// =============================================================================