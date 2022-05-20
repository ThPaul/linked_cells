cmake_minimum_required(VERSION 3.10)
project(particle_timing)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++17")
add_executable(linked_cells linked_cells.cpp)
target_include_directories (linked_cells PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/utils)
