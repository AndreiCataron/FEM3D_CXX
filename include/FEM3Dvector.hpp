#ifndef FEM_FEM3DVECTOR_HPP
#define FEM_FEM3DVECTOR_HPP

#include "FEM3D.hpp"
#include <string>
#include "params.hpp"

using DirichletMap = std::unordered_map<std::size_t, std::vector<double> >;

class FEM3DVector : public FEM3D {
private:
    //const ParamsVector params3d_;
    std::shared_ptr<ParamsVector> params3d_;

protected:
    // store for each node on the boundary [tag : dirichlet bc]
    DirichletMap dirichlet_bc;

public:
    explicit FEM3DVector(std::shared_ptr<ParamsVector> const&);
    FEM3DVector(std::shared_ptr<ParamsVector> const&, std::shared_ptr<Mesh> const&);

    void outputData(std::string, bool, std::vector<double>) override;

    void computeL2Error() override;
    void computeH1Error() override;

private:
    void setDirichletBoundaryConditions() noexcept override;
    void setDirichletBoundaryConditions(const DirichletMap&) noexcept;
    void indexConstrainedNodes() noexcept override;

public:
    //getters
    DirichletMap getDirichletBC();
};

#endif