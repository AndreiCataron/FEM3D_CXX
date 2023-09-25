#include "../include/LinearElasticity3D.hpp"
#include "../include/utils.hpp"
#include <iostream>
#include <gmsh.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/IterativeLinearSolvers>

LinearElasticity3D::LinearElasticity3D(const LinearElasticity3D::ParamsLE &params) : paramsLE_(params), FEM3DVector(params){}

void LinearElasticity3D::computeStiffnessMatrixAndLoadVector() {
    // initialize D matrix
    Eigen::MatrixXd D(6, 6);
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
    // number of nodes in the mesh
    std::vector<std::size_t> tags;
    std::vector<double> co, pc;
    gmsh::model::mesh::getNodes(tags, co, pc, -1, -1, true, false);
    utils::deleteDuplicatesFromVector(tags);
    int noNodes = int(tags.size());
    // no of nodes per element
    // nodeTags.size() = elementTags.size() * noNodesPerElement
    int noNodesPerElement = utils::binom(int(paramsLE_.element_order) + 3, 3);

    // get integration points
    std::vector<double> localCoord, weights;
    std::string intRule = "Gauss" + std::to_string(paramsLE_.quadrature_precision);

    gmsh::model::mesh::getIntegrationPoints(elementType, intRule, localCoord, weights);

    int noInterpolationPoints = int(localCoord.size()) / 3;

    // get basis function values at interpolation points
    std::string functionSpaceType = "Lagrange" + std::to_string(paramsLE_.element_order);
    std::vector<double> basisFunctionsValues;
    int numOrient, numComp;

    gmsh::model::mesh::getBasisFunctions(elementType, localCoord, functionSpaceType, numComp, basisFunctionsValues, numOrient);

    // get gradients at interpolation points
    functionSpaceType = "GradLagrange" + std::to_string(paramsLE_.element_order);
    std::vector<double> basisFunctionsGradients;

    gmsh::model::mesh::getBasisFunctions(elementType, localCoord, functionSpaceType, numComp, basisFunctionsGradients, numOrient);

    int noBasisFunctions = noNodesPerElement;
    // number of columns of the strain matrix B introduced in Larson page 267
    int bNoCols = 3 * noBasisFunctions;

    // get jacobians
    std::vector<double> jacobians, determinants, coord;
    std::vector<double> jacobianCoords = {0.25, 0.25, 0.25};

    gmsh::model::mesh::preallocateJacobians(elementType, 1, false, true, false, jacobians, determinants, coord);
    gmsh::model::mesh::getJacobians(elementType, jacobianCoords, jacobians, determinants, coord);

    // computation of the matrix M^k from Larson page 270

    Eigen::MatrixXd Mk(3 * noNodesPerElement, 3 * noNodesPerElement);

    for (int i = 0; i < noInterpolationPoints; i++) {
        Eigen::MatrixXd Mktemp = Eigen::MatrixXd::Zero(3 * noNodesPerElement, 3);
        for (int j = 0; j < noBasisFunctions; j++) {
            Eigen::MatrixXd basisFunctionMatrix(3, 3);
            basisFunctionMatrix(0, 0) = basisFunctionsValues[i * noInterpolationPoints + j];
            basisFunctionMatrix(1, 1) = basisFunctionsValues[i * noInterpolationPoints + j];
            basisFunctionMatrix(2, 2) = basisFunctionsValues[i * noInterpolationPoints + j];
            Mktemp.block<3, 3>(3 * j, 0) = basisFunctionMatrix;
        }
        Mk = Mk + weights[i] * (Mktemp * Mktemp.transpose());
    }

    // initialize load vector
    load_vector = Eigen::VectorXd::Zero(3 * noNodes);

    // initialize stiffness matrix as Eigen Sparse Matrix
    std::vector<Eigen::Triplet<double> > tripletList;
    tripletList.reserve(elementTags.size() * bNoCols * bNoCols);
    stiffness_matrix.resize(3 * noNodes, 3 * noNodes);

    // matrix used to compute local stiffness matrix for each element
    Eigen::MatrixXd K(bNoCols, bNoCols);

    for (int i = 0; i < noInterpolationPoints; i++) {
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, bNoCols);
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
        Eigen::MatrixXd Ktemp = B.transpose() * D * B;
        K = K + weights[i] * Ktemp;
    }

    // assemble tripletList and load vector by looping through all elements and computing the element stiffness matrix and load vector
    for (int i = 0; i < elementTags.size(); i++) {
        // get tags of nodes in current element
        std::vector<std::size_t> elementNodeTags = std::vector<std::size_t>(nodeTags.begin() + i * noNodesPerElement, nodeTags.begin() + (i + 1) * noNodesPerElement);
        double det = determinants[i];

        Eigen::MatrixXd elementStiffness = det / 6 * K;
        Eigen::MatrixXd elementMass = det / 6 * Mk;

        // find stiffness and load indexes for the displacements in the element
        std::vector<int> elementNodeIndexes;
        elementNodeIndexes.reserve(3 * noNodesPerElement);

        for (const auto& tag : elementNodeTags) {
            int index = nodeIndexes[tag];

            elementNodeIndexes.emplace_back(3 * index);
            elementNodeIndexes.emplace_back(3 * index + 1);
            elementNodeIndexes.emplace_back(3 * index + 2);
        }

        // add triplets to tripletList
        // repeated pairs of indexes are summed up when initializing the sparse stiffness matrix
        for (int k = 0; k < bNoCols; k++) {
            for (int j = 0; j < bNoCols; j++) {
                tripletList.emplace_back(elementNodeIndexes[k], elementNodeIndexes[j], elementStiffness(k, j));
            }
        }

        // the vector f^k from Larson
        std::vector<double> fk;
        fk.reserve(3 * noNodesPerElement);
        for (const auto& tag : elementNodeTags) {
            std::tuple<double, double, double> nodeCoord = node_coordinates[tag];

            for (const auto& component : paramsLE_.f) {
                fk.emplace_back(parseExpression(component, std::get<0>(nodeCoord), std::get<1>(nodeCoord), std::get<2>(nodeCoord)));
            }
        }

        double *ptr = &fk[0];
        Eigen::Map<Eigen::VectorXd> fkConverted(ptr, int(fk.size()));

        Eigen::VectorXd localLoad = elementMass * fkConverted;

        load_vector(elementNodeIndexes) += localLoad;
    }

    stiffness_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

}

void LinearElasticity3D::solveDisplacements() {
    displacements = Eigen::VectorXd(nodeIndexes.size());

    // assemble vector of constrained values
    Eigen::VectorXd constrainedValues(3 * constrainedNodes.size());

    std::unordered_map<int, std::size_t> reverseNodeIndexes = utils::inverseMap(nodeIndexes);

    std::vector<double> tempDirBC;

    for (int i = 0; i < constrainedNodes.size(); i++) {
        std::size_t idx = reverseNodeIndexes[constrainedNodes[i]];

        tempDirBC = dirichlet_bc[idx];
        constrainedValues(3 * i) = tempDirBC[0];
        constrainedValues(3 * i + 1) = tempDirBC[1];
        constrainedValues(3 * i + 2) = tempDirBC[2];
    }

    std::vector<int> freeIndexes;
    freeIndexes.reserve(3 * freeNodes.size());
    for (auto idx : freeNodes) {
        freeIndexes.emplace_back(3 * idx);
        freeIndexes.emplace_back(3 * idx + 1);
        freeIndexes.emplace_back(3 * idx + 2);
    }

    load_vector = load_vector(freeIndexes) -
            stiffness_matrix.bottomRows(3 * freeNodes.size()).leftCols(3 * constrainedNodes.size()) * constrainedValues;

    stiffness_matrix = stiffness_matrix.bottomRightCorner(3 * freeNodes.size(), 3 * freeNodes.size());

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;

    displacements = solver.compute(stiffness_matrix).solve(load_vector);

    try {
        if(solver.info() != Eigen::Success) {
                // solving failed
                throw 1;
            }
        }
    catch (int err) {
        std::cout << "solving the sparse system failed";
        return;
    }
}