set(CMAKE_CXX_STANDARD 20)

find_library(FOO_LIB gmsh)

find_path(FOO_PATH gmsh.h)

find_package(OpenMP REQUIRED)

add_executable(sbc test_setBoundaryConditions.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp ../include/FEM3Dvector.hpp ../src/FEM3Dvector.cpp ../include/LinearElasticity3D.hpp ../src/LinearElasticity3D.cpp ../include/utils.hpp ../src/utils.cpp)
add_executable(cnsbc test_checkNodeSatisfiesBoundaryCondition.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp ../include/FEM3Dvector.hpp ../src/FEM3Dvector.cpp ../include/LinearElasticity3D.hpp ../src/LinearElasticity3D.cpp ../include/utils.hpp ../src/utils.cpp)
add_executable(cpi test_checkParamsInitialization.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp ../include/FEM3Dvector.hpp ../src/FEM3Dvector.cpp ../include/LinearElasticity3D.hpp ../src/LinearElasticity3D.cpp ../include/utils.hpp ../src/utils.cpp)
add_executable(lj learn_jabobians.cpp ../include/utils.hpp ../src/utils.cpp)
add_executable(lnbag learn_isNodalBasisAutomaticallyGenerated.cpp)
add_executable(rtst randomCXXTests.cpp)
add_executable(ret randomEigenTests.cpp)
add_executable(stiff test_stiffnessAndLoad.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp ../include/FEM3Dvector.hpp ../src/FEM3Dvector.cpp ../include/LinearElasticity3D.hpp ../src/LinearElasticity3D.cpp ../include/utils.hpp ../src/utils.cpp)
add_executable(exam testExampleSolutions.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp ../include/FEM3Dvector.hpp ../src/FEM3Dvector.cpp ../include/LinearElasticity3D.hpp ../src/LinearElasticity3D.cpp ../include/utils.hpp ../src/utils.cpp)

target_link_libraries(sbc LINK_PUBLIC "${FOO_LIB}")

target_link_libraries(stiff PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(rtst PUBLIC OpenMP::OpenMP_CXX)

target_include_directories(sbc PUBLIC "${FOO_PATH}")

target_link_libraries(cnsbc LINK_PUBLIC "${FOO_LIB}")

target_include_directories(cnsbc PUBLIC "${FOO_PATH}")

target_link_libraries(cpi LINK_PUBLIC "${FOO_LIB}")

target_include_directories(cpi PUBLIC "${FOO_PATH}")

target_link_libraries(lj LINK_PUBLIC "${FOO_LIB}")

target_include_directories(lj PUBLIC "${FOO_PATH}")

target_link_libraries(lnbag LINK_PUBLIC "${FOO_LIB}")

target_include_directories(lnbag PUBLIC "${FOO_PATH}")

target_link_libraries(stiff LINK_PUBLIC "${FOO_LIB}")

target_include_directories(stiff PUBLIC "${FOO_PATH}")

target_link_libraries(exam LINK_PUBLIC "${FOO_LIB}")

target_include_directories(exam PUBLIC "${FOO_PATH}")

add_test(NAME testSetBoundaryConditions COMMAND sbc)
add_test(NAME testCheckNodeSatisfiesBoundaryCondition COMMAND cnsbc)
add_test(NAME testCheckParamsInitialization COMMAND cpi)
add_test(NAME learnJacobians COMMAND lj)
add_test(NAME learnNodalBasisAutomaticallyGenerated COMMAND lnbag)
add_test(NAME randomTst COMMAND rtst)
add_test(NAME randomEigenTst COMMAND ret)
add_test(NAME testStiffnessAndLoad COMMAND stiff)
add_test(NAME testExamples COMMAND exam)
