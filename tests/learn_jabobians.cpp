#include <gmsh.h>
#include <iostream>
#include <set>

int main(int argc, char **argv) {
    gmsh::initialize(argc, argv);

    gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1, 1000);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.3);
    gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::model::mesh::generate(3);

    gmsh::model::mesh::setOrder(4);

//    std::set<std::string> args(argv, argv + argc);
//    if(!args.count("-nopopup")) gmsh::fltk::run();

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t> > elementTags;
    std::vector<std::vector<std::size_t> > nodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 3);

    for(auto i : elementTypes){
        std::cout << i << " ";
    }
    std::cout << std::endl;
    for(auto j : elementTags){
        for(auto k : j){
            std::cout << k << " ";
        }
    }
    std::cout << "\n";

    for(auto elemType : elementTypes) {
        std::string name;
        int d, order, numv, numpv;
        std::vector<double> param;
        gmsh::model::mesh::getElementProperties(elemType, name, d, order, numv,
                                                param, numpv);
        std::cout << " - Element type: " << name << ", order " << order << "\n";
        std::cout << "   with " << numv << " nodes in param coord: (";
        for(auto p : param) std::cout << p << " ";
        std::cout << ")\n";

        // Integration Points

        std::string intRule = "Gauss1";
        std::vector<double> localCoord, weights;

        std::cout << "\n";
        gmsh::model::mesh::getIntegrationPoints(elemType, intRule, localCoord, weights);
        std::cout << localCoord.size() << "\n";
        for(auto c : localCoord) std::cout << c << " ";
        std::cout << "\n";
        for(auto w : weights) std::cout << w << " ";

        std::string functionSpaceType = "Lagrange2";
        std::vector<double> basisFunctions;
        int numOrient, numComp = 3;

        // Basis functions

        gmsh::model::mesh::getBasisFunctions(elemType, localCoord, functionSpaceType, numComp, basisFunctions, numOrient);

        std::cout << "\n" << basisFunctions.size() << "\n";
        for(auto b : basisFunctions){
            std::cout << b << " ";
        }
        std::cout << "\n" << numOrient;

        // Jacobians

        std::vector<double> jacobians, determinants, coord;

        gmsh::model::mesh::getJacobians(elemType, localCoord, jacobians, determinants, coord);

        std::cout << "\n";
        for(auto i : determinants){
            std::cout << i << " ";
        }

    }

    gmsh::finalize();

    return 0;
}