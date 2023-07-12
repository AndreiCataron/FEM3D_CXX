#ifndef FEM_FEM3D_HPP
#define FEM_FEM3D_HPP

#endif //FEM_FEM3D_HPP

#include <unordered_map>
#include <string>
#include "../include/exprtk/exprtk.hpp"
#include <eigen3/Eigen/Dense>

class FEM3D{
public:
    struct Params{
        // mesh size
        double h;
        // equation of boundary surfaces where Dirichlet BC are imposed
        std::string dirichlet_bc;
        // Quadrature precision
        unsigned int quadrature_precision;
        // Order of Lagrange polynomials
        unsigned int element_order;
    };

private:
    // struct of parameters
    const Params params_;

protected:
    // stiffness matrix
    Eigen::MatrixXd stiffness_matrix;

public:
    // constructor
    FEM3D(const Params&);

    // methods
    double parseExpression(std::string, double, double, double);

    virtual void setBoundaryConditions() = 0;
    virtual void computeStiffnessMatrix() = 0;

    void setupMesh();

    // check if a node is on a part of the boundary where boundary conditions are imposed
    // return 0 if no bc are imposed, 1 for dirichlet, 2 for neumann
    int checkNodeSatisfiesBoundaryEquation(const std::size_t, double, double, double);

    //getters
    FEM3D::Params getParams();
};
