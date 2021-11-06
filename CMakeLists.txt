cmake_minimum_required(VERSION 3.0)

enable_language(Fortran)
enable_testing()

add_compile_options(-O3 -cpp)
add_compile_options(-Wall -cpp -O -fbounds-check -fbacktrace -g -pg -ffpe-trap=invalid,zero,overflow -std=f2008)

project(monte-carlo-simulation Fortran)

add_subdirectory(src)
