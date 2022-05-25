cmake_minimum_required(VERSION 3.10)
project(particle_timing)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++2a")
add_executable(atomic_test atomic_test.cpp)
#target_include_directories (atomic_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/utils)
