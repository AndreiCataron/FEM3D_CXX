#ifndef FEM_FEM3D_HPP
#define FEM_FEM3D_HPP

#endif //FEM_FEM3D_HPP

#include <unordered_map>
#include <string>

class FEM3D{
public:
    struct Params{
        double h;
        std::string dirichlet_bc;
    };

private:
    // struct of parameters
    const Params params_;

    // store for each node on the boundary [tag : dirichlet bc]
    std::unordered_map<std::size_t, double> dirichlet_bc;
    // store for each node on the boundary [tag : satisfies neumann boundary equation]
    // temporary
    std::unordered_map<std::size_t, bool> neumann_bc;

public:
    // constructor
    FEM3D(const Params&);

    // methods
    virtual void setBoundaryConditions() = 0;
    bool checkNodeSatisfiesBoundaryEquation(const std::size_t, double, double, double);

    //getters
    std::unordered_map<std::size_t, double> getDirichletBC();
    FEM3D::Params getParams();
};
