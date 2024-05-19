#ifndef FEM_ALTERNANTKOZLOV_HPP
#define FEM_ALTERNANTKOZLOV_HPP

#include "FEM3D.hpp"
#include "FEM3Dvector.hpp"
#include "LinearElasticity3D.hpp"

class AlternantKozlov {
private:
    std::shared_ptr<LinearElasticity3D> LE;

    // corespunzatoare portiunii de frontiera Gamma_1 din Nachaoui
    std::string dirichletBuffer, initialDirichletBdry, initialNeumannBdry;

    DirichletMap buffer = {};

    std::vector<double> projectionPlane = {1, 0, 0, 0};
    DirichletMap initialDirichletCondition = {};
    std::vector<Eigen::Matrix3d> initialStresses = {};
    std::vector<std::size_t> accesibleBdryTags, unaccessibleBdryTags;

    DirichletMap prevDirichletCondition = {};

    // relaxation factor
    double theta = 1.5;
    stressVector previousStresses = {};
public:
    AlternantKozlov(std::shared_ptr<LinearElasticity3D> const&, std::string);
    AlternantKozlov(std::shared_ptr<LinearElasticity3D> const&, std::vector<double>&, std::string);

    void iteration(int);
    void iterations(int, double, bool);

private:
    void computeRelaxedStresses(stressVector&);
    void unaccessibleDirichletPrep();
    void unaccessibleNeumannPrep();
    void computeSolution();
};


#endif