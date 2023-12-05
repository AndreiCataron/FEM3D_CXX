cmake_minimum_required(VERSION 3.25)
enable_testing()

project(FEM)

set(CMAKE_CTEST_ARGUMENTS "--verbose")

#set(GCC_COVERAGE_COMPILE_FLAGS "-march=native")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )

#if(NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE Release)
#endif()
#
#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_subdirectory(tests)
add_subdirectory(src)

set(CMAKE_CXX_STANDARD 20)



