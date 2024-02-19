#ifndef FEM_ALGORITM_HPP
#define FEM_ALGORITM_HPP

#include "FEM3D.hpp"
#include "FEM3Dvector.hpp"
#include "LinearElasticity3D.hpp"

class Algoritm {
private:
    std::shared_ptr<LinearElasticity3D> LE;
    std::vector<double> projectionPlane = {1, 0, 0, 0};
    DirichletMap initialDirichletCondition = {};
    std::vector<Eigen::Matrix3d> initialStresses = {};
    std::vector<std::size_t> accesibleBdryTags, unaccessibleBdryTags;

    DirichletMap prevDirichletCondition = {};

    // relaxation factor
    double theta = 1.5;
public:
    explicit Algoritm(std::shared_ptr<LinearElasticity3D> const&);
    Algoritm(std::shared_ptr<LinearElasticity3D> const&, std::vector<double>&);

    void iteration(int);
    void iterations(int, double, bool);

private:
    void unaccessibleDirichletPrep();
    void unaccessibleNeumannPrep();
    void computeSolution();
};


#endif
