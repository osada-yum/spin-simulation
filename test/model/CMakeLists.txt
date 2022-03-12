set(TEST_2D-Ising_element test_Ising2d_element.out)

add_executable(${TEST_2D-Ising_element} test_Ising2d_element.f90)

foreach(TEST IN ITEMS ${TEST_2D-Ising_element})
target_include_directories(${TEST}
  PRIVATE
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_BINARY_DIR}/src/random
  ${CMAKE_BINARY_DIR}/src/model
  )
target_link_libraries(${TEST}
  PRIVATE
  ${LIB_SPINS}
  ${LIB_BENCH}
  )
add_test(NAME ${TEST}
  COMMAND $<TARGET_FILE:${TEST}>
  )
endforeach()
