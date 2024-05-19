#ifndef FEM_MESH_HPP
#define FEM_MESH_HPP

#include <string>
#include <memory>
#include <set>
#include "params.hpp"
#include "utils.hpp"
#include "MeshInfo.hpp"

using CoordTuple = std::tuple<double, double, double>;

class Mesh {
protected:
    std::shared_ptr<Params> params;

    GlobalData global;
    BoundaryData bdry;
    ElementData elems;

public:
    // constructors
    Mesh();
    Mesh(int, char**);
    Mesh(int, char**, const std::string&);
    Mesh(int, char**, std::shared_ptr<Params> const&);
    ~Mesh();

    friend class FEM3D;
    friend class FEM3DVector;
    friend class LinearElasticity3D;
    friend class AlternantKozlov;

    void getNodesCoordinates();
    static void cubeMesh(double, double, double, double);
    void initMesh();
    static void showMesh(int, char**);

private:
    void computeInverseJacobians();
    void computeNormalsAtFaceIntegrationPoints();
    void findBoundaryAdjacentElements();

public:
    std::unordered_map<std::size_t, CoordTuple> nodeCoords();
};


#endif
