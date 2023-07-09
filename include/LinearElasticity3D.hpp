#ifndef FEM_LINEARELASTICITY3D_HPP
#define FEM_LINEARELASTICITY3D_HPP

#endif //FEM_LINEARELASTICITY3D_HPP

#include "FEM3Dvector.hpp"

class LinearElasticity3D : public FEM3DVector {
public:
    LinearElasticity3D(const FEM3DVector::ParamsVector&);

    void computeStiffnessMatrix();
};