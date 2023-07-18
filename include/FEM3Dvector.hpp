#ifndef FEM_FEM3DVECTOR_HPP
#define FEM_FEM3DVECTOR_HPP

#endif //FEM_FEM3DVECTOR_HPP

#include "FEM3D.hpp"
#include <string>

class FEM3DVector : public FEM3D {
public:
    struct ParamsVector : FEM3D::Params{
        std::vector<std::string> g;
    };
private:
    const ParamsVector params3d_;

    // store for each node on the boundary [tag : dirichlet bc]
    std::unordered_map<std::size_t, std::vector<double> > dirichlet_bc;
    // store for each node on the boundary [tag : satisfies neumann boundary equation]
    // temporary
    std::unordered_map<std::size_t, bool> neumann_bc;

public:
    explicit FEM3DVector(const ParamsVector&);
    void setBoundaryConditions() override;

    void indexConstrainedNodes() override;

    //getters
    std::unordered_map<std::size_t, std::vector<double> > getDirichletBC();
};