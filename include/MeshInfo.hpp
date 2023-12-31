#ifndef FEM_MESHINFO_HPP
#define FEM_MESHINFO_HPP

#include <unordered_map>
#include <vector>
#include <eigen3/Eigen/Dense>

using CoordTuple = std::tuple<double, double, double>;

struct GlobalData {
    // store for each node [tag : coordinates]
    std::unordered_map<std::size_t, CoordTuple> node_coordinates = {};

    // jacobians
    std::vector<double> jacobians = {}, determinants = {}, globalCoord = {};
    // inverse jacobians
    std::vector<Eigen::MatrixXd> inverse_jacobians;

    // boundary triangles jacobians
    std::vector<double> trianglesDeterminants = {}, trianglesGlobalCoord = {};
};

struct BoundaryData {
    // boundary
    std::vector<std::pair<int, int> > boundary = {};

    std::vector<std::size_t> boundaryFacesTags = {}, boundaryFacesNodes = {};

    int triangleElementType = -1;

    // faces integration points
    std::vector<double> triangleLocalCoord = {}, triangleWeights = {};
    int triangleNoIntegrationPoints = -1;
    // basis functions values at triangleLocalCoord
    std::vector<double> triangleBasisFunctionsValues = {};
    int noNodesPerTriangle = -1;
};

struct ElementData {
    // tetrahedron elements
    std::vector<int> elementTypes = {};
    std::vector<std::vector<std::size_t> > elemTags = {}, nTags = {};
    int elementType = -1;
    std::vector<std::size_t> elementTags = {};
    std::vector<std::size_t> nodeTags = {};

    // integration points
    std::vector<double> localCoord = {}, weights = {};
    int noIntegrationPoints = -1;

    // basis function values at integration  points
    std::vector<double> basisFunctionsValues = {};
    std::vector<double> basisFunctionsGradients = {};
    int noNodesPerElement = -1;
};

#endif //FEM_MESHINFO_HPP
