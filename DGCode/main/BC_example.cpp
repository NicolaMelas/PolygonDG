/*!
    @file   BC_example.cpp
    @author Nicola Melas
    @brief  Check validity of BoundaryCondition
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
#include "ElastodynamicsProblem.hpp"
#include <unordered_map>

using namespace MuParserInterface;
using namespace PolygonDG;

/*!
  @brief The example consists of passing the boundary condition to the problem and check that
  the problem is correctly modified in a way such that the different boundaries
  where boundary conditions are imposed can be correctly recognized
  */

int main(int argc, char** argv) {

GetPot file(argc, argv);

std::string filename = file("filename", "BC_example.txt");

std::cout<<"Reading parameters from "<<filename<<std::endl;

GetPot g2(filename.c_str());

std::size_t ne = g2("ne", 8);
double rho = g2("rho", 1);
double mu = g2("mu", 1);
double lam = g2("lam", 1);
double dt = g2("dt", 0.0001);
std::size_t fem_degree = g2("fem_degree", 2);
double penalty_coeff = g2("penalty_coeff", 10);
std::string source_1 = g2("source_1", "0.0");
std::string source_2 = g2("source_2", "0.0");

muParserXInterface<3> fx;
fx.set_expression(source_1);
fx.DefineConst("mu", mu);
fx.DefineConst("lam", lam);
fx.DefineConst("rho", rho);
muParserXInterface<3> fy;
fy.set_expression(source_2);
fy.DefineConst("mu", mu);
fy.DefineConst("lam", lam);
fy.DefineConst("rho", rho);

std::string exact_sol_1 = g2("exact_sol_1", "0.0");
std::string exact_sol_2 = g2("exact_sol_2", "0.0");

muParserXInterface<3> exact_sol_x;
exact_sol_x.set_expression(exact_sol_1);
exact_sol_x.DefineConst("mu", mu);
exact_sol_x.DefineConst("lam", lam);
exact_sol_x.DefineConst("rho", rho);
muParserXInterface<3> exact_sol_y;
exact_sol_y.set_expression(exact_sol_2);
exact_sol_y.DefineConst("mu", mu);
exact_sol_y.DefineConst("lam", lam);
exact_sol_y.DefineConst("rho", rho);

std::string neum_1 = g2("neumann_1", "0.0");
std::string neum_2 = g2("neumann_2", "0.0");


muParserXInterface<3> neum_x;
neum_x.set_expression(neum_1);
neum_x.DefineConst("mu", mu);
neum_x.DefineConst("lam", lam);
neum_x.DefineConst("rho", rho);
muParserXInterface<3> neum_y;
neum_y.set_expression(neum_2);
neum_y.DefineConst("mu", mu);
neum_y.DefineConst("lam", lam);
neum_y.DefineConst("rho", rho);

std::unordered_map<std::string, double> coeffs = {{"rho",rho},{"mu",mu},{"lam",lam},{"dt",dt},
                                                  {"penalty",penalty_coeff}};
Coefficients coefficients(coeffs, fx, fy);

Point2D p1,p2;

p1 << 0.0,0.0;
p2 << 1.0,1.0;

Mesh Th(p1,p2,ne,1,"Polygon");

Neighbour neighbour(Th);


FeSpace femTh(Th, neighbour, fem_degree);


DirichletBC diri1([](Point2D x)-> bool {return x(0) < 0.0000000001;},exact_sol_x, exact_sol_y, -1, coefficients);
DirichletBC diri2([](Point2D x)-> bool {return x(1) < 0.0000000001;},exact_sol_x, exact_sol_y, -2,coefficients);
DirichletBC diri3([](Point2D x)-> bool {return x(0) > 1-0.0000000001;},exact_sol_x, exact_sol_y, -3, coefficients);
NeumannBC neum1([](Point2D x)-> bool {return x(1) > 1-0.0000000001;},neum_x, neum_y, -4);

std::vector<std::shared_ptr<BoundaryCondition>> bcs;
bcs.push_back(std::make_shared<DirichletBC>(diri1));
bcs.push_back(std::make_shared<DirichletBC>(diri2));
bcs.push_back(std::make_shared<DirichletBC>(diri3));
bcs.push_back(std::make_shared<NeumannBC>(neum1));
ElastodynamicsProblem problem(femTh, coefficients, bcs);

const Neighbour & neighbourp = problem.get_Neighbour();
auto polygons = Th.get_coords_element();

//check if tags are set correctly
std::cout<< std::endl;
std::cout << "                            Tag List" << std::endl;
std::cout<< std::endl;
const VecVecXi & neigh = neighbourp.get_neigh();
std::size_t n_elements = Th.get_ne();

for (std::size_t ii = 0; ii < n_elements; ii++) {
  const auto & loc_coords = polygons[ii]->theVertexes();
  std::size_t int_size = neigh[ii].size();
  for (std::size_t jj = 0; jj < int_size; jj++) {
    if (neigh[ii](jj) < 0) {
      std::cout << "Tag: " << neigh[ii](jj) << " for the edge of vertexes: " << loc_coords[jj](0) <<" "<< loc_coords[jj](1)
        << " and " <<  loc_coords[(jj+1)%int_size](0) << " "<<loc_coords[(jj+1)%int_size](1) << std::endl;
    }
  }
}

return 0;


}
