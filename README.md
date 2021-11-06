# Usage
```console
mkdir build && cd build
cmake ..
cmake --build .
./Ising2d_relaxation.out
./Ising2d_equilibrium.out
```
# Implementation
Modern Fortran implementation of Metropolis method for 2-dimensional ferromagnetic Ising model with skew boundary condition.
## 2D-Ising
You can initilize `Ising2d` by kbt (temperature), x and y (number of sites on the x,y-axis), update it by Metropolis method and calculate its magnetism and energy.
