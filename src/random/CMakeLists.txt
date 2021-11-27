add_library(${LIB_RANDOM}
  STATIC
  random_base_m.f90
  builtin_rand_m.f90
  xorshift_m.f90
  )
add_subdirectory(test)
