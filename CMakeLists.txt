cmake_minimum_required(VERSION 3.0)

enable_language(Fortran)
enable_testing()

set(CMAKE_Fortran_FLAGS           "-cpp")
set(CMAKE_Fortran_FLAGS_DEBUG     "-Og -DDEBUG -Wall -Wextra -fbounds-check -fbacktrace -g -pg -ffpe-trap=invalid,zero,overflow -std=f2008")
set(CMAKE_Fortran_FLAGS_RELEASE   "-O3 -DNDEBUG")

project(monte-carlo-simulation Fortran)

set(LIB_SPINS spins) # utility functions.

add_subdirectory(src)
add_subdirectory(app)
