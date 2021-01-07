# Description

This solver solves the ideal MHD equations in their Eulerian form in 1 dimension
using the the VL+CT (Van-Leer + Constrained Transport) algorithm from Stone &
Gardiner 2009 and the HLLD Riemann solver. It includes support for multiple
types of initial conditions and boundary conditions. See the documentation for
MhdSimulation1D::_setInitialConditions and Grid1D::updateBoundaries respectively
for details.