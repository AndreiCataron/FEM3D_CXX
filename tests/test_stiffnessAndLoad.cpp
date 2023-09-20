#include <iostream>
#include "../include/LinearElasticity3D.hpp"
#include <gmsh.h>
#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Dense>
#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Sparse>
#include <cassert>

int main(int argc, char **argv) {
    LinearElasticity3D::ParamsLE par = {
            0, // h
            "0 == 0", // dirichlet BC
            3, // quadrature precision
            1, // order of lagrage polynomials
            {"2 * x", "3 * y + 1", "x + z"}, // g
            {"0", "0", "0"}, // f
            1, // lambda
            2, // mu
    };

    LinearElasticity3D fem(par);

    gmsh::initialize(argc, argv);

    gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1, 1000);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.5);
    gmsh::option::setNumber("Mesh.MaxNumThreads3D", 10);
    gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::model::mesh::generate(3);

    fem.setupMesh();
    fem.setBoundaryConditions();
    fem.computeStiffnessMatrixAndLoadVector();
    Eigen::SparseMatrix<double> sm = fem.getStiffnessMatrix();
    Eigen::VectorXd lv = fem.getLoadVector();

    std::cout << "Rows: " << sm.rows() << "; Cols: " << sm.cols() << '\n';
    std::cout << "Load Vector Size: " << lv.size() << '\n';

    // check number of nodes in mesh

    std::vector<std::size_t> nodeTags;
    std::vector<double> coord, paramCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, paramCoord, -1, -1, true, false);

    std::cout << "No nodes: " << nodeTags.size();

    //Eigen::VectorXd prod = sm * lv;

    gmsh::finalize();
}