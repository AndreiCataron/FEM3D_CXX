#ifndef FEM_LINEARELASTICITY3D_HPP
#define FEM_LINEARELASTICITY3D_HPP

#include "FEM3Dvector.hpp"
#include "params.hpp"

class LinearElasticity3D : public FEM3DVector {
private:
    //const ParamsLE paramsLE_;
    std::shared_ptr<ParamsLE> paramsLE_;
public:
    explicit LinearElasticity3D(std::shared_ptr<ParamsLE> const&);
    LinearElasticity3D(std::shared_ptr<ParamsLE> const&, Mesh&);

    void computeStiffnessMatrixAndLoadVector() override;
    void solveDisplacements() override;
};

#endif
