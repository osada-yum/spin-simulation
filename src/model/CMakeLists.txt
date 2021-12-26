set(LIB_2D-Ising Ising2d)

foreach(model IN ITEMS ${LIB_2D-Ising})
  add_library(${model}
    STATIC
    Ising2d_m.f90
    )
  target_link_libraries(${model}
    PRIVATE
    ${LIB_UTIL}
    ${LIB_RANDOM}
    )
  target_include_directories(${model}
    PRIVATE
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_BINARY_DIR}/src/random
    )
  target_compile_options(${model}
    PRIVATE
    -cpp
    )
endforeach()

add_library(${LIB_MODEL}
  INTERFACE
  )
target_link_libraries(${LIB_MODEL}
  INTERFACE
  ${LIB_2D-Ising}
  )
