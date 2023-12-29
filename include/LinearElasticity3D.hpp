#ifndef FEM_LINEARELASTICITY3D_HPP
#define FEM_LINEARELASTICITY3D_HPP

#include "FEM3Dvector.hpp"
#include "params.hpp"

class LinearElasticity3D : public FEM3DVector {
private:
    // ParamsLE paramsLE_;
    std::shared_ptr<ParamsLE> paramsLE_;
    // neumann boundary condition
    //std::function<Eigen::Vector3d(double, double, double)> h;
public:
    explicit LinearElasticity3D(std::shared_ptr<ParamsLE> const&);
    LinearElasticity3D(std::shared_ptr<ParamsLE> const&, Mesh&);

    Eigen::Vector3d h(std::vector<double>&, int);
    void computeStiffnessMatrixAndLoadVector() override;
    void solveDirectProblem() override;
};

#endif
