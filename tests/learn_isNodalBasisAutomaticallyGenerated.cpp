#include <gmsh.h>
#include <set>
#include <iostream>

int main(int argc, char **argv) {
    gmsh::initialize(argc, argv);

    gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1, 1000);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.5);
    gmsh::option::setNumber("General.Verbosity", 1);

    gmsh::model::mesh::generate(3);
    gmsh::model::mesh::setOrder(3);

    std::vector<std::size_t> nodeTags;
    std::vector<double> coord, parametricCoord;

    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, -1, -1, false, false);

    std::cout << nodeTags.size();

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();
}