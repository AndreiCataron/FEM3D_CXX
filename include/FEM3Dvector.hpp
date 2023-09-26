#ifndef FEM_FEM3DVECTOR_HPP
#define FEM_FEM3DVECTOR_HPP

#include "FEM3D.hpp"
#include <string>

class FEM3DVector : public FEM3D {
public:
    struct ParamsVector : FEM3D::Params{
        std::vector<std::string> f;
        std::vector<std::string> g;
    };

private:
    const ParamsVector params3d_;

protected:
    // store for each node on the boundary [tag : dirichlet bc]
    std::unordered_map<std::size_t, std::vector<double> > dirichlet_bc;
    // store for each node on the boundary [tag : satisfies neumann boundary equation]
    // temporary
    std::unordered_map<std::size_t, bool> neumann_bc;
    // store for each node [tag : coordinates]
    std::unordered_map<std::size_t, std::tuple<double, double, double> > node_coordinates;

public:
    explicit FEM3DVector(const ParamsVector&);

    void setBoundaryConditions() override;
    void indexConstrainedNodes() override;

    void getNodesCoordinates() override;

    void outputData(std::string) override;

    //getters
    std::unordered_map<std::size_t, std::vector<double> > getDirichletBC();
    std::unordered_map<std::size_t, std::tuple<double, double, double> > getNodeCoordinates();
};

#endif