#include "../include/Mesh.hpp"
#include <gmsh.h>
#include "../include/utils.hpp"
#include <iostream>

Mesh::Mesh() {
    if (gmsh::isInitialized() == 0) {
        gmsh::initialize();
    }
    params = nullptr;
}

Mesh::Mesh(int argc, char **argv) {
    if (gmsh::isInitialized() == 0) {
        gmsh::initialize(argc, argv);
    }
    params = nullptr;
}

Mesh::Mesh(int argc, char **argv, const std::string& meshType) : Mesh(argc, argv) {
    if (meshType.compare("unit_cube") == 1) {
        cubeMesh();
    }
}

Mesh::Mesh(int argc, char **argv, std::shared_ptr<Params> const &par) : Mesh(argc, argv) {
    params = par;
}

Mesh::~Mesh() {
    if (gmsh::isInitialized() == 1) {
        gmsh::finalize();
    }
}

void Mesh::getNodesCoordinates() {
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord, parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, -1, -1, true, false);

    for (int i = 0; i < nodeTags.size(); i++) {
        global.node_coordinates[nodeTags[i]] = {coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]};
    }
}

void Mesh::cubeMesh() {
    gmsh::model::occ::addBox(-1, -1, -1, 2, 2, 2, 1000);
    gmsh::model::occ::synchronize();
}

void Mesh::computeInverseJacobians() {
    global.inverse_jacobians.reserve(int(global.jacobians.size() / 9));

    for (int i = 0; i < int(global.jacobians.size() / 9); i++) {
        // jacobian is global.jacobians[9 * i, 9 * i + 3, 9 * i + 6
        //                             9 * i + 1, 9 * i + 4, 9 * i + 7
        //                             9 * i + 2, 9 * i + 5, 9 * i + 8]
        Eigen::MatrixXd jacobian(3, 3);
        jacobian << global.jacobians[9 * i], global.jacobians[9 * i + 3], global.jacobians[9 * i + 6],
                    global.jacobians[9 * i + 1], global.jacobians[9 * i + 4], global.jacobians[9 * i + 7],
                    global.jacobians[9 * i + 2], global.jacobians[9 * i + 5], global.jacobians[9 * i + 8];

        Eigen::MatrixXd inverseJacobian(3, 3);
        inverseJacobian << - jacobian(1, 2) * jacobian(2, 1) + jacobian(1, 1) * jacobian(2, 2),
                           jacobian(0, 2) * jacobian(2, 1) - jacobian(0, 1) * jacobian(2, 2),
                           - jacobian(0, 2) * jacobian(1, 1) + jacobian(0, 1) * jacobian(1, 2),
                           jacobian(1, 2) * jacobian(2, 0) - jacobian(1, 0) * jacobian(2, 2),
                           - jacobian(0, 2) * jacobian(2, 0) + jacobian(0, 0) * jacobian(2, 2),
                           jacobian(0, 2) * jacobian(1, 0) - jacobian(0, 0) * jacobian(1, 2),
                           - jacobian(1, 1) * jacobian(2, 0) + jacobian(1, 0) * jacobian(2, 1),
                           jacobian(0, 1) * jacobian(2, 0) - jacobian(0, 0) * jacobian(2, 1),
                           - jacobian(0, 1) * jacobian(1, 0) + jacobian(0, 0) * jacobian(1, 1);

        inverseJacobian = 1 / global.determinants[i] * inverseJacobian;
        inverseJacobian.transposeInPlace();

        global.inverse_jacobians.emplace_back(inverseJacobian);
    }
}

void Mesh::computeNormalsAtFaceIntegrationPoints() {
    auto it = global.trianglesGlobalCoord.cbegin();
    for (const auto& b : bdry.boundary) {
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t> > elementTags, nodeTags;
        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, b.first, b.second);

        auto coords = std::vector<double>(it, it + 3 * bdry.triangleNoIntegrationPoints * int(elementTags[0].size()));
        it += bdry.triangleNoIntegrationPoints * int(elementTags[0].size());

        std::vector<double> parametricCoords, normals;

        gmsh::model::getParametrization(2, b.second, coords, parametricCoords);

        gmsh::model::getNormal(b.second, parametricCoords, normals);

        for (int i = 0; i < normals.size() / 3; i++) {
            double *ptr = &normals[3 * i];
            Eigen::Map<Eigen::Vector3d> normalConverted(ptr, 3);
            global.normals.emplace_back(normalConverted);
        }
    }
}

void Mesh::initMesh() {
    if (params == nullptr) {
        gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.1);
        gmsh::option::setNumber("Mesh.MaxNumThreads3D", 10);
        gmsh::option::setNumber("General.Verbosity", 1);

        gmsh::model::mesh::generate(3);

        gmsh::model::mesh::setOrder(1);
    }
    else {
        gmsh::option::setNumber("Mesh.CharacteristicLengthMax", params -> h);
        gmsh::option::setNumber("Mesh.MaxNumThreads3D", params -> max_num_threads);
        gmsh::option::setNumber("General.Verbosity", params -> verbosity);

        gmsh::model::mesh::generate(3);

        gmsh::model::mesh::setOrder(params -> element_order);
    }

    getNodesCoordinates();

    // get elements
    gmsh::model::mesh::getElements(elems.elementTypes, elems.elemTags, elems.nTags, 3);
    // get element type
    elems.elementType = elems.elementTypes[0];
    // element tags
    elems.elementTags = elems.elemTags[0];
    // node tags
    elems.nodeTags = elems.nTags[0];

    // create faces
    gmsh::model::mesh::createFaces();

    // get boundary
    std::vector<std::pair<int, int> > domain_entity;
    gmsh::model::getEntities(domain_entity, 3);

    gmsh::model::getBoundary(domain_entity, bdry.boundary, true, false, false);

    // get face element to identify triangle type
    auto b = bdry.boundary[0];

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t> > elementTags, nodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, b.first, b.second);

    bdry.triangleElementType = elementTypes[0];

    // get integration points for triangles
    std::string triangleIntRule = "Gauss" + std::to_string(params -> triangle_quadrature_precision);

    gmsh::model::mesh::getIntegrationPoints(bdry.triangleElementType, triangleIntRule, bdry.triangleLocalCoord, bdry.triangleWeights);
    bdry.triangleNoIntegrationPoints = int(bdry.triangleLocalCoord.size() / 3);

    // get integration points for tetrahedrons
    std::string intRule = "Gauss" + std::to_string(params -> quadrature_precision);

    gmsh::model::mesh::getIntegrationPoints(elems.elementType, intRule, elems.localCoord, elems.weights);
    elems.noIntegrationPoints = int(elems.localCoord.size()) / 3;

    // get basis function values at integration points
    std::string functionSpaceType = "Lagrange" + std::to_string(params -> element_order);
    int numOrient, numComp;

    gmsh::model::mesh::getBasisFunctions(elems.elementType, elems.localCoord, functionSpaceType, numComp, elems.basisFunctionsValues, numOrient);

    // get basis function values at triangle integration points
    gmsh::model::mesh::getBasisFunctions(bdry.triangleElementType, bdry.triangleLocalCoord, functionSpaceType, numComp, bdry.triangleBasisFunctionsValues, numOrient);

    elems.noNodesPerElement = utils::binomialCoefficient(int(params -> element_order) + 3, 3);
    bdry.noNodesPerTriangle = utils::binomialCoefficient(int(params -> element_order) + 2, 2);

    // get gradients at integration points
    functionSpaceType = "GradLagrange" + std::to_string(params -> element_order);

    gmsh::model::mesh::getBasisFunctions(elems.elementType, elems.localCoord, functionSpaceType, numComp, elems.basisFunctionsGradients, numOrient);

    // get jacobians of triangles
    std::vector<double> triangleJacobianCoords = {0.25, 0.25, 0}, triangleDeterminants, triangleJacobians, triangleGlobalCoords;
    gmsh::model::mesh::preallocateJacobians(bdry.triangleElementType, 1, false, true, false, triangleJacobians, global.trianglesDeterminants, triangleGlobalCoords);
    gmsh::model::mesh::getJacobians(bdry.triangleElementType, triangleJacobianCoords, triangleJacobians, global.trianglesDeterminants, triangleGlobalCoords);
    gmsh::model::mesh::preallocateJacobians(bdry.triangleElementType, int(bdry.triangleLocalCoord.size()), false, true, true, triangleJacobians, triangleDeterminants, global.trianglesGlobalCoord);
    gmsh::model::mesh::getJacobians(bdry.triangleElementType, bdry.triangleLocalCoord, triangleJacobians, triangleDeterminants, global.trianglesGlobalCoord);

    // get the determinant of the jacobian of each element
    // todo-idea Modification
    //   check line 197 below for parallelization
    //   https://gitlab.onelab.info/gmsh/fem/-/blob/master/src/term/Assembler.cpp?ref_type=heads
    std::vector<double> jacobianCoords = {0.25, 0.25, 0.25};

    std::vector<double> dummyCoords;
    gmsh::model::mesh::preallocateJacobians(elems.elementType, 1, true, true, false, global.jacobians, global.determinants, dummyCoords);
    gmsh::model::mesh::getJacobians(elems.elementType, jacobianCoords, global.jacobians, global.determinants, dummyCoords);

    // get the global coordinates of the integration points in each element

    std::vector<double> dummyDeterminants;
    std::vector<double> dummyJacobians;
    gmsh::model::mesh::preallocateJacobians(elems.elementType, int(elems.localCoord.size()), false, true, true, dummyJacobians, dummyDeterminants, global.globalCoord);
    gmsh::model::mesh::getJacobians(elems.elementType, elems.localCoord, dummyJacobians, dummyDeterminants, global.globalCoord);

    computeInverseJacobians();
    computeNormalsAtFaceIntegrationPoints();
}

void Mesh::showMesh(int argc, char **argv) {
    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();
}