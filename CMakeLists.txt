cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(AllHexLayerTopoSolver VERSION 0.1 LANGUAGES C CXX)

# Shared libraries
set(BUILD_SHARED_LIBS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
    message(STATUS "CMAKE_BUILD_TYPE set to ${CMAKE_BUILD_TYPE} (default)")
endif()

# Module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

set (CMAKE_CXX_STANDARD 17)

# Path to executable directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set (ENABLE_SANITY_CHECKS OFF)
  if(ENABLE_SANITY_CHECKS)
    message(STATUS "using sanitizer options")
    add_compile_options("-Wall"
      "-Wextra"
      "-fno-omit-frame-pointer"
      "-Wno-missing-field-initializers"
      "-Wpointer-arith"
      "-Wcast-align"
      "-Wcast-qual"
      "-Wstrict-aliasing"
      "-Wpointer-arith"
      "-Winit-self"
      "-Wredundant-decls"
      "-fsanitize=address"
      "-fsanitize=leak"
      "-fsanitize=undefined")

    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} \
    -fsanitize=address -fsanitize=leak -fsanitize=undefined")
  endif()
endif()

# Print useful information
MESSAGE(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
MESSAGE(STATUS "Library Type: " ${LIB_TYPE})
MESSAGE(STATUS "Compiler flags:" ${CMAKE_CXX_COMPILE_FLAGS})
MESSAGE(STATUS "Compiler cxx debug flags:" ${CMAKE_CXX_FLAGS_DEBUG})
MESSAGE(STATUS "Compiler cxx release flags:" ${CMAKE_CXX_FLAGS_RELEASE})
MESSAGE(STATUS "Compiler cxx min size flags:" ${CMAKE_CXX_FLAGS_MINSIZEREL})
MESSAGE(STATUS "Compiler cxx flags:" ${CMAKE_CXX_FLAGS})

# gecode dependency
find_package(Gecode 6.2 REQUIRED COMPONENTS Driver Int Kernel Search)
find_package(OpenMP)

set(AHLTS_SRC_FILES
  src/solver.cpp
  src/disk_triangulations.cpp
  )
add_library(ahlts_lib ${AHLTS_SRC_FILES})
target_include_directories(ahlts_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(ahlts_lib PUBLIC Gecode::Driver Gecode::Int Gecode::Kernel Gecode::Search
  Gecode::Minimodel Gecode::Support Gecode::Float)
# add_custom_command(
#   TARGET ahlts_lib POST_BUILD
#   COMMAND ${CMAKE_COMMAND} -E copy
#   ${CMAKE_SOURCE_DIR}/assets/disk_trgls.data
#   ${CMAKE_CURRENT_BINARY_DIR}/assets/disk_trgls.data)
# add_custom_command(
#   TARGET ahlts_lib POST_BUILD
#   COMMAND ${CMAKE_COMMAND} -E copy
#   ${CMAKE_SOURCE_DIR}/assets/disk_trgls_2.data
#   ${CMAKE_CURRENT_BINARY_DIR}/assets/disk_trgls_2.data)
add_custom_command(
  TARGET ahlts_lib POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_SOURCE_DIR}/assets/disk_trgls_3.data
  ${CMAKE_CURRENT_BINARY_DIR}/assets/disk_trgls_3.data)

add_executable(basic_tests tests/basic_tests.cpp)
target_link_libraries(basic_tests ahlts_lib)

add_executable(experimental_proof tests/eproof.cpp)
target_link_libraries(experimental_proof ahlts_lib)
