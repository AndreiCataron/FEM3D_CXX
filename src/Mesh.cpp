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
        elems.node_coordinates[nodeTags[i]] = {coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]};
    }
}

void Mesh::cubeMesh() {
    gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1, 1000);
    gmsh::model::occ::synchronize();
}

void Mesh::computeInverseJacobians() {
    elems.inverse_jacobians.reserve(int(elems.jacobians.size() / 9));

    for (int i = 0; i < int(elems.jacobians.size() / 9); i++) {
        // jacobian is elems.jacobians[9 * i, 9 * i + 3, 9 * i + 6
        //                             9 * i + 1, 9 * i + 4, 9 * i + 7
        //                             9 * i + 2, 9 * i + 5, 9 * i + 8]
        Eigen::MatrixXd jacobian(3, 3);
        jacobian << elems.jacobians[9 * i], elems.jacobians[9 * i + 3], elems.jacobians[9 * i + 6],
                    elems.jacobians[9 * i + 1], elems.jacobians[9 * i + 4], elems.jacobians[9 * i + 7],
                    elems.jacobians[9 * i + 2], elems.jacobians[9 * i + 5], elems.jacobians[9 * i + 8];

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

        elems.inverse_jacobians.emplace_back(1 / elems.determinants[i] * inverseJacobian);
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

    // get vertices of tetrahedron on faces
    gmsh::model::mesh::getElementFaceNodes(elems.elementType, 3, elems.facesVerticesTags, -1, true);

    // get boundary
    std::vector<std::pair<int, int> > domain_entity;
    gmsh::model::getEntities(domain_entity, 3);

    gmsh::model::getBoundary(domain_entity, elems.boundary, true, false, false);

    for (auto b : elems.boundary) {
        std::vector<std::size_t> tags;
        std::vector<double> coord, param_coords;
        gmsh::model::mesh::getNodes(tags, coord, param_coords, b.first, b.second, true, false);

        for (auto tag : tags) {
            elems.boundaryTags[tag] = b.second;
        }
    }

    // get integration points
    std::string intRule = "Gauss" + std::to_string(params -> quadrature_precision);

    gmsh::model::mesh::getIntegrationPoints(elems.elementType, intRule, elems.localCoord, elems.weights);
    elems.noIntegrationPoints = int(elems.localCoord.size()) / 3;

    // get basis function values at integration points
    std::string functionSpaceType = "Lagrange" + std::to_string(params -> element_order);
    int numOrient, numComp;

    gmsh::model::mesh::getBasisFunctions(elems.elementType, elems.localCoord, functionSpaceType, numComp, elems.basisFunctionsValues, numOrient);

    elems.noNodesPerElement = utils::binom(int(params -> element_order) + 3, 3);

    // get gradients at integration points
    functionSpaceType = "GradLagrange" + std::to_string(params -> element_order);

    gmsh::model::mesh::getBasisFunctions(elems.elementType, elems.localCoord, functionSpaceType, numComp, elems.basisFunctionsGradients, numOrient);


    // get the determinant of the jacobian of each element
    // todo-idea Modification
    //   check line 197 below for parallelization
    //   https://gitlab.onelab.info/gmsh/fem/-/blob/master/src/term/Assembler.cpp?ref_type=heads
    std::vector<double> jacobianCoords = {0.25, 0.25, 0.25};

    std::vector<double> dummyCoords;
    gmsh::model::mesh::preallocateJacobians(elems.elementType, 1, true, true, false, elems.jacobians, elems.determinants, dummyCoords);
    gmsh::model::mesh::getJacobians(elems.elementType, jacobianCoords, elems.jacobians, elems.determinants, dummyCoords);

    // get the global coordinates of the integration points in each element

    std::vector<double> dummyDeterminants;
    std::vector<double> dummyJacobians;
    gmsh::model::mesh::preallocateJacobians(elems.elementType, int(elems.localCoord.size()), false, true, true, dummyJacobians, dummyDeterminants, elems.globalCoord);
    gmsh::model::mesh::getJacobians(elems.elementType, elems.localCoord, dummyJacobians, dummyDeterminants, elems.globalCoord);

    computeInverseJacobians();
}