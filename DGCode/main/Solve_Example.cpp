/*!
    @file   Solve_Example.cpp
    @author Nicola Melas
    @brief  Solve the ElastodynamicsProblem with mixed BC
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
#include "EvalFunction.hpp"
#include "ExportFunction.hpp"
#include <unordered_map>
#include <chrono>

using namespace MuParserInterface;
using namespace PolygonDG;


/*!
  @brief example that shows how to solve the problem. And print also
  the L2 error.
  */
int main(int argc, char** argv) {

GetPot file(argc, argv);

std::string filename = file("filename", "Solve_example.txt");

std::cout<<"Reading parameters from "<<filename<<std::endl;

GetPot g2(filename.c_str());

std::size_t ne = g2("ne", 8);
double rho = g2("rho", 1);
double mu = g2("mu", 1);
double lam = g2("lam", 1);
double dt = g2("dt", 0.0001);
double T = g2("T", 0.2);
double t0 = g2("t0", 0.0);
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

std::unordered_map<std::string, double> coeffs = {{"rho",rho},{"mu",mu},{"lam",lam},{"dt",dt},{"penalty",penalty_coeff}};
Coefficients coefficients(coeffs, fx, fy);

Point2D p1,p2;

p1 << 0.0,0.0;
p2 << 1.0,1.0;

Mesh Th(p1,p2,ne,1,"Polygon");

Neighbour neighbour(Th);


FeSpace femregion(Th, neighbour, fem_degree);


DirichletBC diri1([](Point2D x)-> bool {return x(0) < 0.0000000001;},exact_sol_x, exact_sol_y, -1, coefficients);
DirichletBC diri2([](Point2D x)-> bool {return x(1) < 0.0000000001;},exact_sol_x, exact_sol_y, -2, coefficients);
DirichletBC diri3([](Point2D x)-> bool {return x(0) > 1-0.0000000001;},exact_sol_x, exact_sol_y, -3, coefficients);
NeumannBC neum1([](Point2D x)-> bool {return x(1) > 1-0.0000000001;},neum_x, neum_y, -4);

std::vector<std::shared_ptr<BoundaryCondition>> bcs;
bcs.push_back(std::make_shared<DirichletBC>(diri1));
bcs.push_back(std::make_shared<DirichletBC>(diri2));
bcs.push_back(std::make_shared<DirichletBC>(diri3));
bcs.push_back(std::make_shared<NeumannBC>(neum1));
ElastodynamicsProblem problem(femregion, coefficients, bcs);
//problem.set_method(1);
problem.assemble();

double time = t0-dt;
VectorD uu2 = evaluate_function(femregion, exact_sol_x, exact_sol_y, time);
time  = t0;
VectorD uu1 = evaluate_function(femregion, exact_sol_x, exact_sol_y, time );

std::size_t nstep = T/dt;
VectorD solution;
problem.update_rhs(time);

auto MM = problem.get_M();
Eigen::SparseLU<SparseMatrixXd> solver;
solver.compute(MM);
VectorD u1old =   solver.solve(uu1);
VectorD u2old =   solver.solve(uu2);

for (std::size_t mnb = 1; mnb <= nstep;mnb++){
  time = mnb * dt ;
  std::cout << "Computing time: " << time << std::endl;
  auto start = std::chrono::system_clock::now();
  problem.solve(solution, u1old, u2old);
  auto middle = std::chrono::system_clock::now();
  problem.update_rhs(time);
  auto end = std::chrono::system_clock::now();
  u2old = u1old;
  u1old = solution;
  std::chrono::duration<double> elapsed_seconds1 = middle-start;
  std::cout << "Solving the system needs " << elapsed_seconds1.count() << std::endl;
  std::chrono::duration<double> elapsed_seconds2 = end-middle;
  std::cout << "Updating the system needs " << elapsed_seconds2.count() << std::endl;
}

exportfunctionvtk(femregion,solution,"Solvexample");

double L2_error = problem.L2_error(solution, exact_sol_x, exact_sol_y, time);
//! Print error and diameter
std::cout << "L2_error: " << L2_error << std::endl;
std::cout << "Diameter of the mesh: " << Th.get_diameter() << std::endl;
return 0;

}
