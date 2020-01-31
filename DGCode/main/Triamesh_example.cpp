/*!
    @file   Triamesh_example.cpp
    @author Nicola Melas
    @brief  Generator of Triangular mesh
*/

#include <iostream>
#include "Mesh.hpp"
#include "Points_Utilities.hpp"
#include "muParserXInterface.hpp"
#include <vector>
#include <string>
#include "TriaGenerator.hpp"
#include <fstream>
#include "GetPot.hpp"


using namespace MuParserInterface;
using namespace PolygonDG;
/*!
  @brief It test the generator of triangular mesh. The mesh can be exported
  to be visualized in paraview. However it is possible to print the points
  of all triangles.
  */

int main(int argc, char** argv) {

GetPot file(argc, argv);

std::string filename = file("filename", "Triamesh.txt");

std::cout<<"Reading parameters from "<<filename<<std::endl;
//reading parameters
GetPot g2(filename.c_str());
std::size_t ne = g2("ne", 3);


std::vector<double> domain = {0,1,0,1};

//default empty mesh
Mesh Th;
//Triangle mesh generator
TriaGenerator generator(domain, ne, ne);
generator.generate(Th);

//output
Th.print_mesh(std::cout, 0);
//paraview compatible file
Th.printvtk("Triamesh");

return 0;



}
