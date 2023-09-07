#include "../include/LinearElasticity3D.hpp"
#include "../include/utils.hpp"
#include <gmsh.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

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
    int noNodes = int(nodeTags.size());
    // no of nodes per element
    // nodeTags.size() = elementTags.size() * noNodesPerElement
    int noNodesPerElement = binom(int(paramsLE_.element_order) + 3, 3);

    // get integration points
    std::vector<double> localCoord, weights;
    std::string intRule = "Gauss" + std::to_string(paramsLE_.quadrature_precision);

    gmsh::model::mesh::getIntegrationPoints(elementType, intRule, localCoord, weights);

    int noInterpolationPoints = int(localCoord.size()) / 3;

    // get gradients at interpolation points
    std::string functionSpaceType = "GradLagrange" + std::to_string(paramsLE_.element_order);
    std::vector<double> basisFunctionsGradients;
    int numOrient, numComp;

    gmsh::model::mesh::getBasisFunctions(elementType, localCoord, functionSpaceType, numComp, basisFunctionsGradients, numOrient);

    int noBasisFunctions = noNodesPerElement;
    // number of columns of the strain matrix B introduced in Larson page 267
    int bNoCols = 3 * noBasisFunctions;

    // get jacobians
    std::vector<double> jacobians, determinants, coord;
    std::vector<double> jacobianCoords = {0.25, 0.25, 0.25};

    gmsh::model::mesh::preallocateJacobians(elementType, 1, false, true, false, jacobians, determinants, coord);
    gmsh::model::mesh::getJacobians(elementType, jacobianCoords, jacobians, determinants, coord);

    // initialize stiffness matrix as Eigen Sparse Matrix

    std::vector<Eigen::Triplet<double> > tripletList;
    //tripletList.reserve();
    Eigen::SparseMatrix<double> stiffnessMatrix (3 * noNodes, 3 * noNodes);

    // assemble tripletList by looping through all elements and computing the element stiffness matrix

    MatrixXd K(bNoCols, bNoCols);

    for (int i = 0; i < noInterpolationPoints; i++) {
        MatrixXd B = MatrixXd::Zero(6, bNoCols);
        for (int k = 0; k < noBasisFunctions; k++) {
            double a, b, c;
            a = basisFunctionsGradients[i * noInterpolationPoints + 3 * k];
            b = basisFunctionsGradients[i * noInterpolationPoints + 3 * k + 1];
            c = basisFunctionsGradients[i * noInterpolationPoints + 3 * k + 2];
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

    // loop through all elements
    for (int i = 0; i < elementTags.size(); i++) {
        std::size_t elemTag = elementTags[i];
        std::vector<std::size_t> elementNodeTags = std::vector<std::size_t>(nodeTags.begin() + i * noNodesPerElement, nodeTags.begin() + (i + 1) * noNodesPerElement);
        double det = determinants[i];

        MatrixXd elementStiffness = det / 6 * K;
        std::vector<int> elementNodeIndexes;
        elementNodeIndexes.reserve(3 * noNodesPerElement);

        for (const auto& tag : elementNodeTags) {
            int index = nodeIndexes[tag];

            elementNodeIndexes.emplace_back(3 * index);
            elementNodeIndexes.emplace_back(3 * index + 1);
            elementNodeIndexes.emplace_back(3 * index + 2);
        }

        for (int k = 0; k < bNoCols; k++) {
            for (int j = 0; j < bNoCols; j++) {
                tripletList.emplace_back(elementNodeIndexes[k], elementNodeIndexes[j], elementStiffness(k, j));
            }
        }
    }

    stiffnessMatrix.setFromTriplets(tripletList.begin(), tripletList.end());

}

void LinearElasticity3D::computeLoadVector() {
    /////////////////////
    // bucata asta e la fel ca la stifness
    // ar trebui cumva rezolvata problema ca sa nu refolosesc cod
    ////////////////////
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
    int noNodes = int(nodeTags.size());
    // no of nodes per element
    // nodeTags.size() = elementTags.size() * noNodesPerElement
    int noNodesPerElement = binom(int(paramsLE_.element_order) + 3, 3);

    // get integration points
    std::vector<double> localCoord, weights;
    std::string intRule = "Gauss" + std::to_string(paramsLE_.quadrature_precision);

    gmsh::model::mesh::getIntegrationPoints(elementType, intRule, localCoord, weights);

    int noInterpolationPoints = int(localCoord.size()) / 3;

    ////////////////////
    // final bucata dublata
    ////////////////////

    // get basis function values at interpolation points
    std::string functionSpaceType = "Lagrange" + std::to_string(paramsLE_.element_order);
    std::vector<double> basisFunctionsValues;
    int numOrient, numComp;

    gmsh::model::mesh::getBasisFunctions(elementType, localCoord, functionSpaceType, numComp, basisFunctionsValues, numOrient);

    int noBasisFunctions = noNodesPerElement;

    // computation of the matrix M^k from Larson page 270

    MatrixXd Mk(3 * noNodesPerElement, 3 * noNodesPerElement);

    for (int i = 0; i < noInterpolationPoints; i++) {
        MatrixXd Mktemp = MatrixXd::Zero(3 * noNodesPerElement, 3);
        for (int j = 0; j < noBasisFunctions; j++) {
            MatrixXd basisFunctionMatrix(3, 3);
            basisFunctionMatrix(0, 0) = basisFunctionsValues[i * noInterpolationPoints + j];
            basisFunctionMatrix(1, 1) = basisFunctionsValues[i * noInterpolationPoints + j];
            basisFunctionMatrix(2, 2) = basisFunctionsValues[i * noInterpolationPoints + j];
            Mktemp.block<3, 3>(0, 3 * j) = basisFunctionMatrix;
        }
        Mk = Mk + weights[i] * (Mktemp * Mktemp.transpose());
    }

}