#ifndef FEM_FEM3D_HPP
#define FEM_FEM3D_HPP

#include <unordered_map>
#include <string>
#include "../include/exprtk/exprtk.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

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
    // for parseExpression
    double _x{}, _y{}, _z{};
    symbol_table_t symbol_table;
    expression_t expression;
    parser_t parser;

protected:
    // dictionary of type tag : index
    std::unordered_map<std::size_t, int> nodeIndexes;
    // indexes of constrained nodes
    std::vector<int> constrainedNodes;
    // indexes of free nodes
    std::vector<int> freeNodes;
    // stiffness matrix
    Eigen::SparseMatrix<double> stiffness_matrix;
    // load vector
    Eigen::VectorXd load_vector;
    // solution in displacements
    Eigen::VectorXd displacements;

public:
    // constructor
    explicit FEM3D(Params );

    // methods
    double parseExpression(const std::string&, double, double, double);

    void setupMesh();
    virtual void indexConstrainedNodes() = 0;
    void indexFreeNodes();
    virtual void getNodesCoordinates() = 0;

    virtual void setBoundaryConditions() = 0;
    // check if a node is on a part of the boundary where boundary conditions are imposed
    // return 0 if no bc are imposed, 1 for dirichlet, 2 for neumann
    int checkNodeSatisfiesBoundaryEquation(double, double, double);
    virtual void computeStiffnessMatrixAndLoadVector() = 0;
    virtual void solveDisplacements() = 0;

    virtual void outputData(std::string) = 0;

    //getters
    FEM3D::Params getParams();
    std::unordered_map<std::size_t, int> getNodeIndexes();
    std::vector<int> getConstrainedNodes();
    std::vector<int> getFreeNodes();
    Eigen::SparseMatrix<double> getStiffnessMatrix();
    Eigen::VectorXd getLoadVector();
    Eigen::VectorXd getDisplacements();
};

#endif