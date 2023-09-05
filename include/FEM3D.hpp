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
    // dictionary of type tag : index
    std::unordered_map<std::size_t, int> nodeIndexes;
    // indexes of constrained nodes
    std::vector<int> constrainedNodes;
    // indexes of free nodes
    std::vector<int> freeNodes;
    // stiffness matrix
    Eigen::MatrixXd stiffness_matrix;

public:
    // constructor
    explicit FEM3D(const Params&);

    // methods
    double parseExpression(const std::string&, double, double, double);

    virtual void setBoundaryConditions() = 0;
    virtual void computeStiffnessMatrix() = 0;

    void setupMesh();
    virtual void indexConstrainedNodes() = 0;
    void indexFreeNodes();

    // check if a node is on a part of the boundary where boundary conditions are imposed
    // return 0 if no bc are imposed, 1 for dirichlet, 2 for neumann
    int checkNodeSatisfiesBoundaryEquation(double, double, double);

    //getters
    FEM3D::Params getParams();
    std::unordered_map<std::size_t, int> getNodeIndexes();
    std::vector<int> getConstrainedNodes();
    std::vector<int> getFreeNodes();
};
