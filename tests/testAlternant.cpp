#include <iostream>
#include "../include/LinearElasticity3D.hpp"
#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Dense>
#include <memory>
#include <vector>
#include "../include/Algoritm.hpp"

int main(int argc, char **argv) {
    auto exact = [] (double x, double y, double z) {
        double nu = 0.33, E = 70.0, sig = 1.5;
        return std::vector<double>{sig / (2 * E) * (x * x + nu * y * y + nu * z * z),
                                   - nu * sig / E * x * y,
                                   - nu * sig / E * x * z};
    };

    auto grad = [] (double x, double y, double z) {
        double nu = 0.33, E = 70.0, sig = 1.5;
        Eigen::Matrix3d grd;
        grd << x, nu * y, nu * z,
                -nu * y, -nu * x, 0,
                -nu * z, 0, -nu * x;
        grd = sig / E * grd;
        return grd;
    };

    auto f = [] (double x, double y, double z) {return std::vector<double>{-1.5, 0, 0};};
    auto g = [] (double x, double y, double z) {
        double nu = 0.33, E = 70.0, sig = 1.5;
        return std::vector<double>{sig / (2 * E) * (x * x + nu * y * y + nu * z * z),
                                   - nu * sig / E * x * y,
                                   - nu * sig / E * x * z};
    };

    auto initialg = [] (double x, double y, double z) {
        return std::vector<double>{0, 0, 0};
    };

    auto par = std::make_shared<ParamsLE>(ParamsLE{
            0.2 , // h,
            10,
            1,
            "x != 0.5", // dirichlet BC
            "x == 0.5", // neumann BC
            3, // quadrature precision
            3, // triangle quadrature precision
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

    auto msh = std::make_shared<Mesh>(argc, argv, par);

    Mesh::cubeMesh(0, 0, 0, 0.5);
    msh -> initMesh();

    auto fem = std::make_shared<LinearElasticity3D>(par, msh);

    Algoritm alg(fem);

    alg.iterations(7, 0, true);

}