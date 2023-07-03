set(CMAKE_CXX_STANDARD 20)

find_library(FOO_LIB gmsh)

find_path(FOO_PATH gmsh.h)

add_executable(sbc test_setBoundaryConditions.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp)
add_executable(cnsbc test_checkNodeSatisfiesBoundaryCondition.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp)
add_executable(cpi test_checkParamsInitialization.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp)

target_link_libraries(sbc LINK_PUBLIC "${FOO_LIB}")

target_include_directories(sbc PUBLIC "${FOO_PATH}")

target_link_libraries(cnsbc LINK_PUBLIC "${FOO_LIB}")

target_include_directories(cnsbc PUBLIC "${FOO_PATH}")

target_link_libraries(cpi LINK_PUBLIC "${FOO_LIB}")

target_include_directories(cpi PUBLIC "${FOO_PATH}")

add_test(NAME testSetBoundaryConditions COMMAND sbc)
add_test(NAME testCheckNodeSatisfiesBoundaryCondition COMMAND cnsbc)
add_test(NAME testCheckParamsInitialization COMMAND cpi)