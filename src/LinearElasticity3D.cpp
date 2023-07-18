#include "../include/LinearElasticity3D.hpp"
#include "../include/utils.hpp"
#include <gmsh.h>
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;

LinearElasticity3D::LinearElasticity3D(const LinearElasticity3D::ParamsLE &params) : paramsLE_(params), FEM3DVector(params){}

void LinearElasticity3D::computeStiffnessMatrix() {
    // initialize D matrix
    MatrixXd D(6, 6);
    double lam = paramsLE_.lambda;
    double mu = paramsLE_.mu;
    D << lam + 2 * mu, lam, lam, 0, 0, 0,
         lam, lam + 2 * mu, lam, 0, 0, 0,
         lam, lam, lam + 2 * mu, 0, 0, 0,
         0, 0, 0, mu, 0, 0,
         0, 0, 0, 0, mu, 0,
         0, 0, 0, 0, 0, mu;

    // get elements
    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t> > elemTags, nTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, nTags, 3);

    // get element type
    int elementType = elemTypes[0];
    // element tags
    std::vector<std::size_t> elementTags = elemTags[0];
    // node tags
    std::vector<std::size_t> nodeTags = nTags[0];
    // no of nodes per element
    int noNodes = binom(int(paramsLE_.element_order) + 3, 3);

    // get integration points
    std::vector<double> localCoord, weights;
    std::string intRule = "Gauss" + std::to_string(paramsLE_.quadrature_precision);

    gmsh::model::mesh::getIntegrationPoints(elementType, intRule, localCoord, weights);

    int noInterpolationPoints = int(localCoord.size()) / 3;

    //get gradients at interpolation points
    std::string functionSpaceType = "GradLagrange" + std::to_string(paramsLE_.element_order);
    std::vector<double> basisFunctions;
    int numOrient, numComp;

    gmsh::model::mesh::getBasisFunctions(elementType, localCoord, functionSpaceType, numComp, basisFunctions, numOrient);

    int noBasisFunctions = noNodes;
    int bNoCols = 3 * noBasisFunctions;

    MatrixXd K(bNoCols, bNoCols);

    for(int i = 0; i < noInterpolationPoints; i++) {
        MatrixXd B = MatrixXd::Zero(6, bNoCols);
        for(int k = 0; k < noBasisFunctions; k++) {
            double a, b, c;
            a = basisFunctions[i * noInterpolationPoints + 3 * k];
            b = basisFunctions[i * noInterpolationPoints + 3 * k + 1];
            c = basisFunctions[i * noInterpolationPoints + 3 * k + 2];
            B(0, 3 * k) = a;
            B(1, 3 * k + 1) = b;
            B(2, 3 * k + 2) = c;
            B(3, 3 * k) = b;
            B(3, 3 * k + 1) = a;
            B(4, 3 * k + 1) = c;
            B(4, 3 * k + 2) = b;
            B(5, 3 * k) = c;
            B(5, 3 * k + 2) = a;
        }
        MatrixXd Ktemp = B.transpose() * D * B;
        K = K + weights[i] * Ktemp;
    }

    // get jacobians
    std::vector<double> jacobians, determinants, coord;
    std::vector<double> jacobianCoords = {0.25, 0.25, 0.25};

    gmsh::model::mesh::preallocateJacobians(elementType, 1, false, true, false, jacobians, determinants, coord);
    gmsh::model::mesh::getJacobians(elementType, jacobianCoords, jacobians, determinants, coord);

    // loop through all elements
    for(int i = 0; i < nodeTags.size(); i++) {
        std::size_t elemTag = elementTags[i];
        std::vector<std::size_t> elementNodeTags = std::vector<std::size_t>(nodeTags.begin() + i * noNodes, nodeTags.begin() + (i + 1) * noNodes);
        double det = determinants[i];
    }

}