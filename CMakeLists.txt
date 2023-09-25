cmake_minimum_required(VERSION 3.25)
enable_testing()

project(FEM)

add_subdirectory(tests)
add_subdirectory(src)

set(CMAKE_CXX_STANDARD 20)



