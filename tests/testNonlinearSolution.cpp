#include <iostream>
#include "../include/LinearElasticity3D.hpp"
#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Dense>
#include <memory>
#include <vector>

int main(int argc, char **argv) {
//    auto exact = [] (double x, double y, double z) {
//        double nu = 0.33, E = 70.0, sig = 1.5;
//        return std::vector<double>{sig / (2 * E) * (x * x + nu * y * y + nu * z * z),
//                                   - nu * sig / E * x * y,
//                                   - nu * sig / E * x * z};
//    };
//
//    auto grad = [] (double x, double y, double z) {
//        double nu = 0.33, E = 70.0, sig = 1.5;
//        Eigen::Matrix3d grd;
//        grd << x, nu * y, nu * z,
//               -nu * y, -nu * x, 0,
//               -nu * z, 0, -nu * x;
//        grd = sig / E * grd;
//        return grd;
//    };
//
//    auto f = [] (double x, double y, double z) {return std::vector<double>{-1.5, 0, 0};};
//    auto g = [] (double x, double y, double z) {
//        double nu = 0.33, E = 70.0, sig = 1.5;
//        return std::vector<double>{sig / (2 * E) * (x * x + nu * y * y + nu * z * z),
//                                   - nu * sig / E * x * y,
//                                   - nu * sig / E * x * z};
//    };
//
//    auto par = std::make_shared<ParamsLE>(ParamsLE{
//            0.15 , // h,
//            10,
//            1,
//            "x != 1", // dirichlet BC
//            "x == 1", // neumann BC
//            3, // quadrature precision
//            3, // triangle quadrature precision
//            1, // order of lagrange polynomials
//            exact, // exact
//            grad, // solution gradient
//            f, // f
//            g, // g
//            -1, // lambda
//            -1, // mu
//            0.33, // nu
//            70.0 // E
//    });

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
//    auto g = [] (double x, double y, double z) {
//        double nu = 0.33, E = 70.0, sig = 1.5;
//        return std::vector<double>{sig / (2 * E) * (x * x + nu * y * y + nu * z * z),
//                                   - nu * sig / E * x * y,
//                                   - nu * sig / E * x * z};
//    };

    auto g = [] (double x, double y, double z) {
        return std::vector<double>{0, 0, 0};
    };

    auto par = std::make_shared<ParamsLE>(ParamsLE{
            0.05 , // h,
            10,
            1,
            "x > 0", // dirichlet BC
            "x == 0", // neumann BC
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

    std::cout << "Salut neliniar " << '\n';


    auto msh = std::make_shared<Mesh>(argc, argv, par);

    Mesh::cubeMesh(0, 0, 0, 0.5);
    msh -> initMesh();

    LinearElasticity3D fem(par, msh);

    fem.setBoundaryConditions();
    fem.computeIntegrationPointsStresses();
    auto start = std::chrono::steady_clock::now();

    fem.computeStiffnessMatrixAndLoadVector();

    auto end = std::chrono::steady_clock::now();

    fem.solveDirectProblem();

    auto diff = end - start;

    std::cout << "Stiff:" << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

    Eigen::SparseMatrix<double> sm = fem.getStiffnessMatrix();
    Eigen::VectorXd lv = fem.getLoadVector();

    std::cout << "Rows: " << sm.rows() << "; Cols: " << sm.cols() << "; Non-zero elements: " << sm.nonZeros() << '\n';
    std::cout << "Load Vector Size: " << lv.size() << '\n';

    start = std::chrono::steady_clock::now();
    fem.computeL2Error();
    std::cout << "L2 error: " << fem.getL2Error();
    end = std::chrono::steady_clock::now();
    diff = end - start;

    std::cout << "\nL2: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

    start = std::chrono::steady_clock::now();
    fem.computeH1Error();
    std::cout << "H1 error: " << fem.getH1Error();
    end = std::chrono::steady_clock::now();
    diff = end - start;

    std::cout << "\nH1: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

    //std::vector<double> plane = {0, 1, 0, 0};
    std::vector<double> plane = {1, 0, 0, 0};
    fem.outputData("/Users/andrei/CLionProjects/FEM/outputs/out.txt", true, plane);

}