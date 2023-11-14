#ifndef FEM_FEM3D_HPP
#define FEM_FEM3D_HPP

#include <unordered_map>
#include <string>
#include "../include/exprtk/exprtk.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "Mesh.hpp"
#include "params.hpp"
#include <memory>

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

// todo-idea Integration on threads
//   check below
//   https://gitlab.onelab.info/gmsh/fem/-/blob/master/src/post/Integrate.cpp

// todo-idea Move boundary node information to Mesh class
//   am 3 locuri in care fac acelasi for (auto b : boundary)
//   super rau !!!

class FEM3D{
private:
    // struct of parameters
    // const Params params_;
    std::shared_ptr<Params> params_;
    // for parseExpression
    double _x{}, _y{}, _z{};
    symbol_table_t symbol_table;
    expression_t expression;
    parser_t parser;

protected:
    // the mesh
    Mesh mesh;
    // dictionary of type tag : index
    std::unordered_map<std::size_t, int> nodeIndexes;
    // indexes of constrained nodes
    std::vector<int> constrainedNodes;
    // indexes of free nodes
    std::vector<int> freeNodes;
    // faces where neumann boundary conditions are imposed
    std::vector<std::size_t> neumannFaceNodes = {};
    std::vector<std::size_t> neumannFaceTags = {};
    // stiffness matrix
    Eigen::SparseMatrix<double> stiffness_matrix;
    // load vector
    Eigen::VectorXd load_vector;
    // solution in displacements
    Eigen::VectorXd displacements;

public:
    // constructors
    explicit FEM3D(std::shared_ptr<Params> const&);
    FEM3D(std::shared_ptr<Params> const&, Mesh&);

    // methods
    double parseExpression(const std::string&, double, double, double);

    virtual void indexConstrainedNodes() = 0;
    void indexFreeNodes();

    virtual void setBoundaryConditions() = 0;
    void setNeumannBoundaryConditions();
    // check if a node is on a part of the boundary where boundary conditions are imposed
    // return 0 if no bc are imposed, 1 for dirichlet, 2 for neumann
    int checkNodeSatisfiesBoundaryEquation(double, double, double);
    virtual void computeStiffnessMatrixAndLoadVector() = 0;
    virtual void solveDisplacements() = 0;

    virtual void outputData(std::string) = 0;

    virtual double computeL2Error() = 0;

    //getters
    std::unordered_map<std::size_t, int> getNodeIndexes();
    std::vector<int> getConstrainedNodes();
    std::vector<int> getFreeNodes();
    Eigen::SparseMatrix<double> getStiffnessMatrix();
    Eigen::VectorXd getLoadVector();
    Eigen::VectorXd getDisplacements();
};

#endif