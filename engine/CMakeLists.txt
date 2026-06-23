# CMake build for the alloy chemical-ordering Monte Carlo engine.
#
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build -j
#
# Produces build/alloy_mc.  Requires an MPI C++ compiler.
cmake_minimum_required(VERSION 3.12)
project(alloy_mc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

find_package(MPI REQUIRED COMPONENTS CXX)

add_executable(alloy_mc
  src/config.cpp
  src/lattice.cpp
  src/alloy_model.cpp
  src/parallel_tempering.cpp
  src/main.cpp)

target_include_directories(alloy_mc PRIVATE include)
target_link_libraries(alloy_mc PRIVATE MPI::MPI_CXX)
target_compile_options(alloy_mc PRIVATE -Wall -Wextra)
