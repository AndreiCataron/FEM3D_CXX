project(TESTS)

add_executable(TESTS test_setBoundaryConditions.cpp)
add_test(NAME testSetBoundsryConditions COMMAND TESTS)