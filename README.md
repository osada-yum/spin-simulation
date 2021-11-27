# Usage
```console
mkdir _build && cd _build
cmake ..
cmake --build .
./Ising2d_relaxation.out
./Ising2d_equilibrium.out
```
or
```console
fpm build --flag '-cpp'
fpm run --flag '-cpp' Ising2d_equilibirum
fpm run --flag '-cpp' Ising2d_relaxation
```
# Implementation
Modern Fortran implementation of Metropolis method for 2-dimensional ferromagnetic Ising model with skew boundary condition.
## 2D-Ising
You can initilize `Ising2d` by kbt (temperature), x and y (number of sites on the x,y-axis), update it by Metropolis method and calculate its magnetism and energy.
