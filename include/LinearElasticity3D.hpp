#ifndef FEM_LINEARELASTICITY3D_HPP
#define FEM_LINEARELASTICITY3D_HPP

#include "FEM3Dvector.hpp"

class LinearElasticity3D : public FEM3DVector {
public:
    struct ParamsLE : FEM3DVector::ParamsVector {
        // Lame constants
        double lambda;
        double mu;
        // Poisson ratio
        double nu;
        // Young Modulus
        double E;
    };
private:
    const ParamsLE paramsLE_;
public:
    explicit LinearElasticity3D(const ParamsLE&, const Mesh&);

    void computeStiffnessMatrixAndLoadVector() override;
    void solveDisplacements() override;
};

#endif
