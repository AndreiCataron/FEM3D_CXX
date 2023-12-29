#include "../include/LinearElasticity3D.hpp"
#include "../include/utils.hpp"
#include <iostream>
#include <set>
#include <gmsh.h>
#include <eigen3/Eigen/IterativeLinearSolvers>

LinearElasticity3D::LinearElasticity3D(std::shared_ptr<ParamsLE> const& params) : FEM3DVector(params), paramsLE_(params) {}

LinearElasticity3D::LinearElasticity3D(std::shared_ptr<ParamsLE> const& params, Mesh &msh) : FEM3DVector(params, msh), paramsLE_(params) {}

Eigen::Vector3d LinearElasticity3D::h(std::vector<double>& coord, const int tag) {
    // compute strain tensor
    double x = coord[0], y = coord[1], z = coord[2];
    Eigen::Matrix3d strain = 0.5 * (paramsLE_ -> solution_gradient(x, y, z) + (paramsLE_ -> solution_gradient(x, y, z)).transpose());

    // compute stress tensor
    Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d stress = (paramsLE_ -> E) / (1 + paramsLE_ -> nu) * (strain + (paramsLE_ -> nu) / (1 - 2 * paramsLE_ -> nu) * strain.trace() * identity);

    // get normal
    std::vector<double> parametricCoord, normal;

    gmsh::model::getParametrization(2, tag, coord, parametricCoord);

    gmsh::model::getNormal(tag, parametricCoord, normal);

    // convert normal to Eigen::Vector3d
    double *ptr = &normal[0];
    Eigen::Map<Eigen::Vector3d> normalConverted(ptr, 3);

    Eigen::Vector3d rez = stress * normalConverted;

    return rez;
}

void LinearElasticity3D::computeStiffnessMatrixAndLoadVector() {
    // initialize D matrix
    Eigen::MatrixXd D(6, 6);
    double lam = paramsLE_ -> lambda;
    double mu = paramsLE_ -> mu;
    D << lam + 2 * mu, lam, lam, 0, 0, 0,
         lam, lam + 2 * mu, lam, 0, 0, 0,
         lam, lam, lam + 2 * mu, 0, 0, 0,
         0, 0, 0, mu, 0, 0,
         0, 0, 0, 0, mu, 0,
         0, 0, 0, 0, 0, mu;

    // number of nodes in the mesh
    // first, delete duplicate nodes (which appear because they are shared by multiple surfaces)
    std::vector<std::size_t> tags;
    std::vector<double> co, pc;
    gmsh::model::mesh::getNodes(tags, co, pc, -1, -1, true, false);
    utils::deleteDuplicatesFromVector(tags);
    int noNodes = int(tags.size());

    // number of columns of the strain matrix B introduced in Larson page 267
    int bNoCols = 3 * mesh.elems.noNodesPerElement;

    // computation of the matrix M^k from Larson page 270

    Eigen::MatrixXd Mk(3 * mesh.elems.noNodesPerElement, 3 * mesh.elems.noNodesPerElement);

    for (int i = 0; i < mesh.elems.noIntegrationPoints; i++) {
        Eigen::MatrixXd Mktemp = Eigen::MatrixXd::Zero(3 * mesh.elems.noNodesPerElement, 3);
        for (int j = 0; j < mesh.elems.noNodesPerElement; j++) {
            Eigen::Matrix3d basisFunctionMatrix = Eigen::Matrix3d::Zero();
            basisFunctionMatrix(0, 0) = mesh.elems.basisFunctionsValues[i * mesh.elems.noNodesPerElement + j];
            basisFunctionMatrix(1, 1) = mesh.elems.basisFunctionsValues[i * mesh.elems.noNodesPerElement + j];
            basisFunctionMatrix(2, 2) = mesh.elems.basisFunctionsValues[i * mesh.elems.noNodesPerElement + j];
            Mktemp.block<3, 3>(3 * j, 0) = basisFunctionMatrix;
        }
        Mk = Mk + mesh.elems.weights[i] * (Mktemp * Mktemp.transpose());
    }

    // initialize load vector
    load_vector = Eigen::VectorXd::Zero(3 * noNodes);

    // initialize stiffness matrix as Eigen Sparse Matrix
    std::vector<Eigen::Triplet<double> > tripletList;
    tripletList.reserve(mesh.elems.elementTags.size() * bNoCols * bNoCols);
    stiffness_matrix.resize(3 * noNodes, 3 * noNodes);

    auto start = std::chrono::steady_clock::now();

    // assemble tripletList and load vector by looping through all elements and computing the element stiffness matrix and load vector
    for (int i = 0; i < mesh.elems.elementTags.size(); i++) {
        // get tags of nodes in current element
        std::vector<std::size_t> elementNodeTags = std::vector<std::size_t>(
                mesh.elems.nodeTags.begin() + i * mesh.elems.noNodesPerElement,
                mesh.elems.nodeTags.begin() + (i + 1) * mesh.elems.noNodesPerElement);

        // compute element stiffness matrix
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero(bNoCols, bNoCols);

        for (int j = 0; j < mesh.elems.noIntegrationPoints; j++) {
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, bNoCols);

            for (int k = 0; k < mesh.elems.noNodesPerElement; k++) {
                Eigen::Vector3d referenceGrads, elementGrads;
                referenceGrads << mesh.elems.basisFunctionsGradients[j * 3 * mesh.elems.noNodesPerElement + 3 * k],
                                  mesh.elems.basisFunctionsGradients[j * 3 * mesh.elems.noNodesPerElement + 3 * k + 1],
                                  mesh.elems.basisFunctionsGradients[j * 3 * mesh.elems.noNodesPerElement + 3 * k + 2];

                elementGrads = mesh.elems.inverse_jacobians[i] * referenceGrads;

                B(0, 3 * k) = elementGrads(0);
                B(1, 3 * k + 1) = elementGrads(1);
                B(2, 3 * k + 2) = elementGrads(2);
                B(3, 3 * k) = elementGrads(1);
                B(3, 3 * k + 1) = elementGrads(0);
                B(4, 3 * k + 1) = elementGrads(2);
                B(4, 3 * k + 2) = elementGrads(1);
                B(5, 3 * k) = elementGrads(2);
                B(5, 3 * k + 2) = elementGrads(0);
            }

            K = K + mesh.elems.weights[j] * B.transpose() * D * B;
        }

        double det = mesh.elems.determinants[i];

        // Eigen::MatrixXd elementStiffness = det * K;
        // Eigen::MatrixXd elementMass = det * Mk;

        // find stiffness and load indexes for the displacements in the element
        std::vector<int> elementNodeIndexes;
        elementNodeIndexes.reserve(3 * mesh.elems.noNodesPerElement);

        for (const auto &tag: elementNodeTags) {
            int index = nodeIndexes[tag];

            elementNodeIndexes.emplace_back(3 * index);
            elementNodeIndexes.emplace_back(3 * index + 1);
            elementNodeIndexes.emplace_back(3 * index + 2);
        }

        // add triplets to tripletList
        // repeated pairs of indexes are summed up when initializing the sparse stiffness matrix
        for (int k = 0; k < bNoCols; k++) {
            for (int j = 0; j < bNoCols; j++) {
                tripletList.emplace_back(elementNodeIndexes[k], elementNodeIndexes[j], det * K(k, j));
            }
        }

        // the vector f^k from Larson
        std::vector<double> fk;
        fk.reserve(3 * mesh.elems.noNodesPerElement);
        for (const auto &tag: elementNodeTags) {
            // std::tuple<double, double, double> nodeCoord = mesh.elems.node_coordinates[tag];

            for (auto component: paramsLE_ -> f(std::get<0>(mesh.elems.node_coordinates[tag]),
                                                std::get<1>(mesh.elems.node_coordinates[tag]),
                                                std::get<2>(mesh.elems.node_coordinates[tag]))) {
                fk.emplace_back(component);
            }
        }

        double *ptr = &fk[0];
        Eigen::Map<Eigen::VectorXd> fkConverted(ptr, int(fk.size()));

        // Eigen::VectorXd localLoad = elementMass * fkConverted;


        load_vector(elementNodeIndexes) = load_vector(elementNodeIndexes) + det * Mk * fkConverted;
    }

    auto end = std::chrono::steady_clock::now();

    auto diff = end - start;

    std::cout << "STIFF 1: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

    // neumann contribution
    for (int i = 0; i < mesh.elems.boundaryFacesTags.size(); i++) {
        std::size_t faceTag = mesh.elems.boundaryFacesTags[i];
        if (neumannBoundaryTriangles.find(faceTag) != neumannBoundaryTriangles.end()) {
            std::vector<std::size_t> faceNodeTags = std::vector<std::size_t>(
                    mesh.elems.boundaryFacesNodes.begin() + i * mesh.elems.noNodesPerTriangle,
                    mesh.elems.boundaryFacesNodes.begin() + (i + 1) * mesh.elems.noNodesPerTriangle);

            // compute integral of h * phi for all basis functions
            for (int k = 0; k < mesh.elems.noNodesPerTriangle; k++) {
                Eigen::Vector3d integral = Eigen::Vector3d::Zero();

                for (int j = 0; j < mesh.elems.triangleNoIntegrationPoints; j++) {
                    //std::cout << mesh.elems.trianglesGlobalCoord.size() << '\n';
                    std::vector<double> coord = std::vector<double>(
                            mesh.elems.trianglesGlobalCoord.begin() +
                            i * 3 * mesh.elems.triangleNoIntegrationPoints + 3 * j,
                            mesh.elems.trianglesGlobalCoord.begin() +
                            i * 3 * mesh.elems.triangleNoIntegrationPoints + 3 * j + 3);

                    integral += mesh.elems.triangleWeights[j] *
                                mesh.elems.triangleBasisFunctionsValues[j * mesh.elems.noNodesPerTriangle + k] *
                                h(coord, int(neumannBoundaryTriangles[faceTag]));
                }

                integral *= mesh.elems.trianglesDeterminants[i];

                int index = nodeIndexes[faceNodeTags[k]];
                std::vector<int> nodeIndexes = {3 * index, 3 * index + 1, 3 * index + 2};

//#pragma omp critical
//                {
                    load_vector(nodeIndexes) = load_vector(nodeIndexes) + integral;
//                }
            }
        }
    }
//}

    stiffness_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

}

void LinearElasticity3D::solveDirectProblem() {
    displacements = Eigen::VectorXd(3 * nodeIndexes.size());

    // assemble vector of constrained values
    Eigen::VectorXd constrainedValues(3 * constrainedNodes.size());

    std::unordered_map<int, std::size_t> reverseNodeIndexes = utils::inverseMap(nodeIndexes);

    std::vector<double> tempDirBC;

    for (int idx : constrainedNodes) {
        std::size_t tag = reverseNodeIndexes[idx];

        tempDirBC = dirichlet_bc[tag];

        constrainedValues(3 * idx) = tempDirBC[0];
        constrainedValues(3 * idx + 1) = tempDirBC[1];
        constrainedValues(3 * idx + 2) = tempDirBC[2];
    }

    displacements.head(3 * constrainedNodes.size()) = constrainedValues;

    load_vector = (load_vector.tail(3 * freeNodes.size()) -
            stiffness_matrix.bottomRows(3 * freeNodes.size()).leftCols(3 * constrainedNodes.size()) * constrainedValues).eval();

    stiffness_matrix = stiffness_matrix.bottomRightCorner(3 * freeNodes.size(), 3 * freeNodes.size()).eval();

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;

    auto start = std::chrono::steady_clock::now();

    displacements.tail(3 * freeNodes.size()) = solver.compute(stiffness_matrix).solve(load_vector);

    auto end = std::chrono::steady_clock::now();

    auto diff = end - start;

    std::cout << "SOLVING TIME:" << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;


    try {
        if (solver.info() != Eigen::Success) {
                // solving failed
                throw 1;
        }
    }
    catch (int err) {
        std::cout << "solving the sparse system failed";
        return;
    }
}