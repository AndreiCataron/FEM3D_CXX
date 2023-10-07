#ifndef FEM_PARAMS_HPP
#define FEM_PARAMS_HPP

#include <string>

struct Params{
    // mesh size
    double h;
    // maximum number of threads for meshing
    unsigned int max_num_threads = 10;
    // gmsh verbosity
    unsigned int verbosity = 1;
    // equation of boundary surfaces where Dirichlet BC are imposed
    std::string dirichlet_bc;
    // Quadrature precision
    unsigned int quadrature_precision;
    // Order of Lagrange polynomials
    int element_order;
};

struct ParamsVector : Params{
    std::vector<std::string> f;
    std::vector<std::string> g;
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
