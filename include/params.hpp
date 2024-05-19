#ifndef FEM_PARAMS_HPP
#define FEM_PARAMS_HPP

#include <string>
#include <functional>
#include "eigen_common.hpp"

struct Params{
    // mesh size
    double h;
    // maximum number of threads for meshing
    unsigned int max_num_threads = 1;
    // gmsh verbosity
    unsigned int verbosity = 1;
    // equation of boundary surfaces where Dirichlet BC are imposed
    std::string dirichlet_bc;
    // equation of boundary surfaces where Neumann BC are imposed
    std::string neumann_bc;
    // Quadrature precision
    unsigned int quadrature_precision;
    // Quadrature precision on triangles
    unsigned int triangle_quadrature_precision;
    // Order of Lagrange polynomials
    int element_order;
public:
    void swapBC() {
        std::swap(dirichlet_bc, neumann_bc);
    };
};

struct ParamsVector : Params{
    std::function<std::vector<double>(double, double, double)> exact_solution;
    std::function<Eigen::Matrix3d(double, double, double)> solution_gradient;
    std::function<std::vector<double>(double, double, double)> f;
    std::function<std::vector<double>(double, double, double)> g;
};

struct ParamsLE : ParamsVector {
    // Lame constants
    double lambda = -1;
    double mu = -1;
    // Poisson ratio
    double nu = -1;
    // Young Modulus
    double E = -1;
};

#endif
