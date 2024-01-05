#include <iostream>
#include "../include/LinearElasticity3D.hpp"
#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Dense>
#include <memory>
#include <vector>

int main(int argc, char **argv) {
    auto exact = [] (double x, double y, double z) {return std::vector<double>{x, y, z};};
    auto grad = [] (double x, double y, double z) {
        Eigen::Matrix3d grd;
        grd << 1, 0, 0,
                0, 1, 0,
                0, 0, 1;
        return grd;
    };
    auto f = [] (double x, double y, double z) {return std::vector<double>{0, 0, 0};};
    auto g = [] (double x, double y, double z) {return std::vector<double>{x, y, z};};

    auto par = std::make_shared<ParamsLE>(ParamsLE{
            0.1 , // h,
            10,
            1,
            "z > 0", // dirichlet BC
            "z == 0", // neumann BC
            3, // quadrature precision
            2, // triangle quadrature precision
            1, // order of lagrange polynomials
            exact, // exact
            grad, // solution gradient
            f, // f
            g, // g
            -1, // lambda
            -1, // mu
            0.33, // nu
            70.0 // E
    });

    utils::checkParamsLE(*par);

    Eigen::initParallel();

    std::cout << "Salut " << par -> mu << '\n';


    auto msh = std::make_shared<Mesh>(Mesh(argc, argv, par));

    Mesh::cubeMesh();
    msh -> initMesh();

    LinearElasticity3D fem(par, msh);

    auto rez = fem.checkNodeSatisfiesBoundaryEquation(0.2, 1, 0);
    std::cout << rez;

    return 0;
}