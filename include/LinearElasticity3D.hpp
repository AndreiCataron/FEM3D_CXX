#ifndef FEM_LINEARELASTICITY3D_HPP
#define FEM_LINEARELASTICITY3D_HPP

#include "FEM3Dvector.hpp"
#include "params.hpp"
#include <mutex>

using stressVector = std::vector<Eigen::Matrix3d>;

class LinearElasticity3D : public FEM3DVector {
private:
    // ParamsLE paramsLE_;
    std::shared_ptr<ParamsLE> paramsLE_;
    // matrices used for stiffness computation;
    Eigen::MatrixXd D, Mk;
    std::vector<Eigen::Triplet<double> > tripletList;
    // stresses at face integration points
    stressVector integrationPointsStresses = {};
    // std::vector<Eigen::Matrix3d> approximatedStresses = {};
    // for parallelism
    std::mutex mtx_stiff, mtx_neu;

public:
    explicit LinearElasticity3D(std::shared_ptr<ParamsLE> const&);
    LinearElasticity3D(std::shared_ptr<ParamsLE> const&, std::shared_ptr<Mesh> const&);

    void resetBoundaryConditions() noexcept override;

    void computeIntegrationPointsStresses() noexcept;
    void computeApproximatedStresses() noexcept;

    void computeStiffnessMatrixAndLoadVector() override;

    void solveDirectProblem() override;

private:
    Eigen::Vector3d h(int);

    void stiffnessIterations(int, int);
    void neumannIterations(int, int);

public:
    void setStresses(stressVector);
    stressVector getStresses();
};

#endif
