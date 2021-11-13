set(EXE_2D-Ising_EMC Ising2d_equilibrium.out) # Equilibrium Monte-Carlo
set(EXE_2D-Ising_NER Ising2d_relaxation.out)  # Non-Equilibrium relaxation
set(LIB_2D-Ising     Ising2d)

add_library(${LIB_2D-Ising}
  STATIC
  utility_m.f90
  Ising2d_m.f90
  )

add_executable(${EXE_2D-Ising_EMC}
  Ising2d_equilibirum.f90
  )
target_link_libraries(${EXE_2D-Ising_EMC}
  PRIVATE
  ${LIB_2D-Ising}
  )
set_target_properties(${EXE_2D-Ising_EMC}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )

add_executable(${EXE_2D-Ising_NER}
  Ising2d_relaxation.f90
  )
target_link_libraries(${EXE_2D-Ising_NER}
  PRIVATE
  ${LIB_2D-Ising}
  )

set_target_properties(${EXE_2D-Ising_NER}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )

add_subdirectory(test)
