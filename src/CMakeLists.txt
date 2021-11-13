# Equilibrium Monte-Carlo
set(EXE_2D-Ising_EMC_BUILTIN Ising2d_equilibrium_builtin_rand.out)
set(EXE_2D-Ising_EMC_XOR96   Ising2d_equilibrium_xor96.out)
# Non-Equilibrium relaxation
set(EXE_2D-Ising_NER_BUILTIN Ising2d_relaxation_builtin_rand.out)
set(EXE_2D-Ising_NER_XOR96   Ising2d_relaxation_xor96.out)
set(LIB_2D-Ising     Ising2d)

add_library(${LIB_2D-Ising}
  STATIC
  utility_m.f90
  Ising2d_m.f90
  )
target_link_libraries(${LIB_2D-Ising}
  PUBLIC
  ${LIB_RANDOM}
  )
target_include_directories(${LIB_2D-Ising}
  PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}/../random
  )
target_compile_options(${LIB_2D-Ising}
  PRIVATE
  -cpp
  )

foreach(emc IN ITEMS ${EXE_2D-Ising_EMC_BUILTIN} ${EXE_2D-Ising_EMC_XOR96})
  add_executable(${emc}
    Ising2d_equilibirum.f90
    )
  target_link_libraries(${emc}
    PRIVATE
    ${LIB_2D-Ising}
    )
  target_include_directories(${emc}
    PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/../random
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
    ${LIB_2D-Ising}
    )
  target_include_directories(${ner}
    PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/../random
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

add_subdirectory(test)
