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
#include "/opt/homebrew/Cellar/libomp/17.0.4/include/omp.h"

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

// todo-idea Integration on threads
//   check below
//   https://gitlab.onelab.info/gmsh/fem/-/blob/master/src/post/Integrate.cpp

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
    std::shared_ptr<Mesh> mesh;
    // dictionary of type tag : index
    std::unordered_map<std::size_t, int> nodeIndexes;
    // indexes of constrained nodes
    std::vector<int> constrainedNodes;
    // indexes of free nodes
    std::vector<int> freeNodes;
    // faces where neumann boundary conditions are imposed
    // store for each neumann triangle on the boundary the tag of the respective boundary surface
    std::unordered_map<std::size_t, std::size_t> neumannBoundaryTriangles = {};
    // stiffness matrix
    Eigen::SparseMatrix<double> stiffness_matrix;
    // load vector
    Eigen::VectorXd load_vector;
    // solution in displacements
    Eigen::VectorXd displacements;
    // errors
    double l2_error = -1, h1_error = -1;

public:
    // constructors
    explicit FEM3D(std::shared_ptr<Params> const&);
    FEM3D(std::shared_ptr<Params> const&, std::shared_ptr<Mesh> const&);

    // methods
    void setBoundaryConditions(auto&&...) noexcept;

    virtual void resetBoundaryConditions() noexcept = 0;

    virtual void computeStiffnessMatrixAndLoadVector() = 0;
    virtual void solveDirectProblem() = 0;

    virtual void outputData(const std::string& file, bool boundaryError, const std::vector<double>& plane) = 0;

    virtual void computeL2Error() = 0;
    virtual void computeH1Error() = 0;

protected:
    double parseExpression(const std::string&, double, double, double);

    // check if a node is on a part of the boundary where boundary conditions are imposed
    // return 0 if no bc are imposed, 1 for dirichlet, 2 for neumann
    int checkNodeSatisfiesBoundaryEquation(double, double, double);

    void setNeumannBoundaryConditions() noexcept;
    virtual void setDirichletBoundaryConditions() noexcept = 0;

    virtual void indexConstrainedNodes() noexcept = 0;
    void indexFreeNodes() noexcept;

public:
    // getters
    [[nodiscard]] const Mesh& getMesh() const;
    std::shared_ptr<Params> getParamsPointer();
    std::unordered_map<std::size_t, int> getNodeIndexes();
    std::vector<int> getConstrainedNodes();
    std::vector<int> getFreeNodes();
    Eigen::SparseMatrix<double> getStiffnessMatrix();
    Eigen::VectorXd getLoadVector();
    Eigen::VectorXd getDisplacements();
    [[nodiscard]] double getL2Error() const;
    [[nodiscard]] double getH1Error() const;
};

#endif