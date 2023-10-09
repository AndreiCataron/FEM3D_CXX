#include <iostream>
#include "../include/LinearElasticity3D.hpp"
#include "../include/utils.hpp"
#include <gmsh.h>
#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Dense>
#include "../include/params.hpp"
#include <iostream>
#include "../include/Mesh.hpp"
#include <memory>

int main(int argc, char **argv) {
    auto par = std::make_shared<ParamsLE>(ParamsLE{
            0.3, // h,
            10,
            1,
            "0 == 0", // dirichlet BC
            2, // quadrature precision
            2, // order of lagrange polynomials
            { "0", "0", "0" }, // f
            { "x", "y", "z" }, // g
            -1, // lambda
            -1, // mu
            0.34, // nu
            12.864 // E
    });

//    ParamsLE par2 = {
//            0.1, // h
//            10,
//            1,
//            "0 == 0",
//            1,
//            1,
//            {"-1.9831", "-1.9831", "50.016"},
//            {"x * z + 1 / 56 * x ^ 2", "y * z + 1 / 56 * y ^ 2", "- z ^ 2 + 1 / 56 * z ^ 2"},
//            56,
//            36,
//            -1,
//            -1
//    };

    utils::checkParamsLE(*par);

    std::cout << "Salut" << '\n';


    Mesh msh(argc, argv, par);

    msh.cubeMesh();
    msh.initMesh();

    LinearElasticity3D fem(*par, msh);

    fem.setBoundaryConditions();
    auto start = std::chrono::steady_clock::now();

    fem.computeStiffnessMatrixAndLoadVector();
    fem.solveDisplacements();

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;

    std::cout << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

    Eigen::SparseMatrix<double> sm = fem.getStiffnessMatrix();
    Eigen::VectorXd lv = fem.getLoadVector();

    std::cout << "Rows: " << sm.rows() << "; Cols: " << sm.cols() << "; Non-zero elements: " << sm.nonZeros() << '\n';
    std::cout << "Load Vector Size: " << lv.size() << '\n';


    fem.outputData("/Users/andrei/CLionProjects/FEM/outputs/out.txt");
//
//    std::set<std::string> args(argv, argv + argc);
//    if(!args.count("-nopopup")) gmsh::fltk::run();
//

}