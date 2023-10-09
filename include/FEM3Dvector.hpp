#ifndef FEM_FEM3DVECTOR_HPP
#define FEM_FEM3DVECTOR_HPP

#include "FEM3D.hpp"
#include <string>
#include "params.hpp"

class FEM3DVector : public FEM3D {
private:
    const ParamsVector params3d_;

protected:
    // store for each node on the boundary [tag : dirichlet bc]
    std::unordered_map<std::size_t, std::vector<double> > dirichlet_bc;
    // store for each node on the boundary [tag : satisfies neumann boundary equation]
    // temporary
    std::unordered_map<std::size_t, bool> neumann_bc;

public:
    explicit FEM3DVector(const ParamsVector&);
    FEM3DVector(const ParamsVector&, Mesh&);

    void setBoundaryConditions() override;
    void indexConstrainedNodes() override;

    void outputData(std::string) override;

    double computeL2Error() override;

    //getters
    std::unordered_map<std::size_t, std::vector<double> > getDirichletBC();
};

#endif