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
fpm run --flag '-cpp' Ising2d_equilibrium
fpm run --flag '-cpp' Ising2d_relaxation
```

# TODO
- [ ] compare speed of `f77 style`, `f90 with module style` and  `f90 with OOP`.
  - [x] f77 style
  - [x] f77 style with function
  - [ ] f90 with module style
  - [x] f90 with OOP
- [ ] implement several spin models and algorithms.
  - [ ] 2-dimensional Ising model
	- [x] with Metropolis method

# Implementation
Modern Fortran implementation of Metropolis method for 2-dimensional ferromagnetic Ising model with skew boundary condition.

## 2D-Ising
You can initialize `Ising2d` by kbt (temperature), x and y (number of sites on the x,y-axis), update it by Metropolis method and calculate its magnetism and energy.
