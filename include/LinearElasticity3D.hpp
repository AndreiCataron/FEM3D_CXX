#ifndef FEM_LINEARELASTICITY3D_HPP
#define FEM_LINEARELASTICITY3D_HPP

#include "FEM3Dvector.hpp"
#include "params.hpp"

class LinearElasticity3D : public FEM3DVector {
private:
    const ParamsLE paramsLE_;
public:
    explicit LinearElasticity3D(const ParamsLE&);
    LinearElasticity3D(const ParamsLE&, Mesh&);

    void computeStiffnessMatrixAndLoadVector() override;
    void solveDisplacements() override;
};

#endif
