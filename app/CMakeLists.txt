# Equilibrium Monte-Carlo
set(EXE_2D-Ising_EMC_BUILTIN_F77 Ising2d_equilibrium_builtin_rand_f77.out)
set(EXE_2D-Ising_EMC_BUILTIN Ising2d_equilibrium_builtin_rand.out)
set(EXE_2D-Ising_EMC_XOR96   Ising2d_equilibrium_xor96.out)
# Non-Equilibrium relaxation
set(EXE_2D-Ising_NER_BUILTIN Ising2d_relaxation_builtin_rand.out)
set(EXE_2D-Ising_NER_XOR96   Ising2d_relaxation_xor96.out)

add_executable(${EXE_2D-Ising_EMC_BUILTIN_F77}
  Ising2d_equilibrium_f77.f
  )
set_target_properties(${EXE_2D-Ising_EMC_BUILTIN_F77}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )
target_compile_options(${EXE_2D-Ising_EMC_BUILTIN_F77}
  PRIVATE
  )
target_link_libraries(${EXE_2D-Ising_EMC_BUILTIN_F77}
  PRIVATE
  ${LIB_BENCH}
  )
target_include_directories(${EXE_2D-Ising_EMC_BUILTIN_F77}
  PRIVATE
  ${CMAKE_BINARY_DIR}/src
  )

foreach(emc IN ITEMS ${EXE_2D-Ising_EMC_BUILTIN} ${EXE_2D-Ising_EMC_XOR96})
  add_executable(${emc}
    Ising2d_equilibrium.f90
    )
  target_link_libraries(${emc}
    PRIVATE
    ${LIB_SPINS}
    ${LIB_BENCH}
    )
  target_include_directories(${emc}
    PRIVATE
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_BINARY_DIR}/src/model
    ${CMAKE_BINARY_DIR}/src/random
    )
  set_target_properties(${emc}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endforeach()

foreach(ner IN ITEMS ${EXE_2D-Ising_NER_BUILTIN} ${EXE_2D-Ising_NER_XOR96})
  add_executable(${ner}
    Ising2d_relaxation.f90
    )
  target_link_libraries(${ner}
    PRIVATE
    ${LIB_SPINS}
    ${LIB_BENCH}
    )
  target_include_directories(${ner}
    PRIVATE
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_BINARY_DIR}/src/model
    ${CMAKE_BINARY_DIR}/src/random
    )
  set_target_properties(${ner}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endforeach()

foreach(builtin IN ITEMS ${EXE_2D-Ising_EMC_BUILTIN} ${EXE_2D-Ising_NER_BUILTIN})
  target_compile_options(${builtin}
    PRIVATE
    -DRAND_GEN_TYPE=builtin_rand_wrapper -cpp
    )
endforeach()
foreach(xor96 IN ITEMS ${EXE_2D-Ising_EMC_XOR96} ${EXE_2D-Ising_NER_XOR96})
  target_compile_options(${xor96}
    PRIVATE
    -DRAND_GEN_TYPE=xor96 -cpp
    )
endforeach()
