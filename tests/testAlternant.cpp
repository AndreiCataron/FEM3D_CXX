#include <iostream>
#include "../include/LinearElasticity3D.hpp"
#include <memory>
#include <vector>
#include "../include/AlternantKozlov.hpp"
#include <limits>

int main(int argc, char **argv) {
    std::cout << std::numeric_limits<float>::denorm_min() << '\n';

    auto exact = [] (double x, double y, double z) {
        double nu = 0.33, E = 70.0, sig = 1.5;
        return std::vector<double>{sig / (2 * E) * (x * x + nu * y * y + nu * z * z),
                                   - nu * sig / E * x * y,
                                   - nu * sig / E * x * z};
    };

    auto grad = [] (double x, double y, double z) -> Eigen::Matrix3d{
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

    auto par = std::make_shared<ParamsLE>(ParamsLE{
            0.05 , // h,
            1,
            1,
//            "x > 0 and (y == 0 or y == 0.5 or z == 0 or z == 0.5 or x == 0.5)", // dirichlet BC
//            "x == 0", // neumann BC
//            "x == 0.5",
//            "x < 0.5 and (y == 0 or y == 0.5 or z == 0 or z == 0.5 or x == 0)",
            "(x > 0.25 and (y == 0 or y == 0.5 or z == 0 or z == 0.5 or x == 0.5))",
            "(x < 0.1 and (y == 0 or y == 0.5 or z == 0 or z == 0.5 or x == 0))",
            // x >=0.2 and x <= 0.3 and (y == 0 or y == 0.5 or z == 0 or z == 0.5 or x == 0.5)
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

    //Eigen::initParallel();

    auto msh = std::make_shared<Mesh>(argc, argv, par);

    Mesh::cubeMesh(0, 0, 0, 0.5);
    msh -> initMesh();

    auto fem = std::make_shared<LinearElasticity3D>(par, msh);

    std::cout << omp_get_num_threads() << " Threads \n";

    std::string buffer = "(x >= 0.1 and x <= 0.25 and (y == 0 or y == 0.5 or z == 0 or z == 0.5))";
    AlternantKozlov alg(fem, buffer);

    //std::cout << fem -> checkNodeSatisfiesCustomEquation(buffer, 0.4, 0.3, 1);
    alg.iterations(2, 0, true);

}