#ifndef FEM_FEM3DVECTOR_HPP
#define FEM_FEM3DVECTOR_HPP

#include "FEM3D.hpp"
#include <string>
#include "params.hpp"

class FEM3DVector : public FEM3D {
private:
    //const ParamsVector params3d_;
    std::shared_ptr<ParamsVector> params3d_;

protected:
    // store for each node on the boundary [tag : dirichlet bc]
    std::unordered_map<std::size_t, std::vector<double> > dirichlet_bc;

public:
    explicit FEM3DVector(std::shared_ptr<ParamsVector> const&);
    FEM3DVector(std::shared_ptr<ParamsVector> const&, Mesh&);

    void setDirichletBoundaryConditions() noexcept override;
    void indexConstrainedNodes() noexcept override;

    void outputData(std::string, bool, std::vector<double>) override;

    void computeL2Error() override;
    void computeH1Error() override;

    //getters
    std::unordered_map<std::size_t, std::vector<double> > getDirichletBC();
};

#endif