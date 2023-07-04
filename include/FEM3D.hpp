#ifndef FEM_FEM3D_HPP
#define FEM_FEM3D_HPP

#endif //FEM_FEM3D_HPP

#include <unordered_map>
#include <string>
#include "../include/exprtk/exprtk.hpp"

class FEM3D{
public:
    struct Params{
        double h;
        std::string dirichlet_bc;
    };

private:
    // struct of parameters
    const Params params_;

public:
    // constructor
    FEM3D(const Params&);

    // methods
    double parseExpression(std::string, double, double, double);

    virtual void setBoundaryConditions() = 0;

    // check if a node is on a part of the boundary where boundary conditions are imposed
    // return 0 if no bc are imposed, 1 for dirichlet, 2 for neumann
    int checkNodeSatisfiesBoundaryEquation(const std::size_t, double, double, double);



    //getters
    FEM3D::Params getParams();
};
