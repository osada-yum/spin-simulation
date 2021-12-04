set(LIB_UTIL utils) # utility functions.

add_library(${LIB_UTIL}
  STATIC
  utility_m.f90
  )

set(LIB_RANDOM rand_gen) # random number genrator base class.
add_subdirectory(random)
set(LIB_MODEL model)     # library for spin models.
add_subdirectory(model)

add_library(${LIB_SPINS}
  INTERFACE
  )
target_link_libraries(${LIB_SPINS}
  INTERFACE
  ${LIB_UTIL}
  ${LIB_RANDOM}
  ${LIB_MODEL}
  )
add_library(${LIB_BENCH}
  STATIC
  benchmark_m.f90
  )
