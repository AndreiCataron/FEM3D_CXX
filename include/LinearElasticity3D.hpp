#ifndef FEM_LINEARELASTICITY3D_HPP
#define FEM_LINEARELASTICITY3D_HPP

#include "FEM3Dvector.hpp"
#include "params.hpp"
#include <mutex>

class LinearElasticity3D : public FEM3DVector {
private:
    // ParamsLE paramsLE_;
    std::shared_ptr<ParamsLE> paramsLE_;
    // matrices used for stiffness computation;
    Eigen::MatrixXd D, Mk;
    std::vector<Eigen::Triplet<double> > tripletList;
    // for parallelism
    std::mutex mtx_triplet, mtx_load, mtx_neu, mtx_h;

public:
    explicit LinearElasticity3D(std::shared_ptr<ParamsLE> const&);
    LinearElasticity3D(std::shared_ptr<ParamsLE> const&, std::shared_ptr<Mesh> const&);

    Eigen::Vector3d h(std::vector<double>&, int);
    void computeStiffnessMatrixAndLoadVector() override;

    void stiffnessIterations(int, int);
    void neumannIterations(int, int);

    void solveDirectProblem() override;
};

#endif
