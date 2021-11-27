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
# set_property(TARGET ${LIB_SPINS} PROPERTY
#   IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/src/lib${LIB_UTIL}.a
#   IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/src/random/lib${LIB_RANDOM}.a
#   IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/src/model/lib${LIB_MODEL}.a
#   )
# target_link_libraries(${LIB_SPINS}
#   PRIVATE
#   ${LIB_UTIL}
#   ${LIB_RANDOM}
#   ${LIB_MODEL}
#   )
# target_include_directories(${LIB_SPINS}
#   PRIVATE
#   ${CMAKE_CURRENT_BINARY_DIR}/src
#   ${CMAKE_CURRENT_BINARY_DIR}/src/random
#   ${CMAKE_CURRENT_BINARY_DIR}/src/model
#   )
