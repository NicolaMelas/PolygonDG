/*!
    @file   Polygonal_convrate_example.cpp
    @author Nicola Melas
    @brief  Compute the L2-convergence rate
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

using namespace MuParserInterface;
using namespace PolygonDG;

/*!
  @brief Compute the L2 convergence rate in a Polygonal mesh.
  The expect theorical convergence rate in L2 norm is the number of
  polynomial degree + 1.
  */
int main(int argc, char** argv) {

GetPot file(argc, argv);

std::string filename = file("filename", "PolyL2convergence.txt");

std::cout<<"Reading parameters from "<<filename<<std::endl;

GetPot g2(filename.c_str());


std::size_t n_repetitions = g2("reps",3);

std::vector<std::size_t> ne_s(n_repetitions);

for (std::size_t c = 0; c < n_repetitions; c++){
  std::string toberead = "ne"+std::to_string(c);
  std::size_t ne = g2(toberead.c_str(), 5);
  ne_s[c] = ne;
}

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

// Source
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

//Exact solution/Dirichlet
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

// Neumann
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

//Gradient exact sol for the error
std::string exact_grad_1x = g2("exact_grad_1x", "0.0");
std::string exact_grad_1y = g2("exact_grad_1y", "0.0");
std::string exact_grad_2x = g2("exact_grad_2x", "0.0");
std::string exact_grad_2y = g2("exact_grad_2y", "0.0");

muParserXInterface<3> grad_1_x;
grad_1_x.set_expression(exact_grad_1x);
grad_1_x.DefineConst("mu", mu);
grad_1_x.DefineConst("lam", lam);
grad_1_x.DefineConst("rho", rho);
muParserXInterface<3> grad_1_y;
grad_1_y.set_expression(exact_grad_1y);
grad_1_y.DefineConst("mu", mu);
grad_1_y.DefineConst("lam", lam);
grad_1_y.DefineConst("rho", rho);
muParserXInterface<3> grad_2_x;
grad_2_x.set_expression(exact_grad_2x);
grad_2_x.DefineConst("mu", mu);
grad_2_x.DefineConst("lam", lam);
grad_2_x.DefineConst("rho", rho);
muParserXInterface<3> grad_2_y;
grad_2_y.set_expression(exact_grad_2y);
grad_2_y.DefineConst("mu", mu);
grad_2_y.DefineConst("lam", lam);
grad_2_y.DefineConst("rho", rho);

//Solving the problem
std::unordered_map<std::string, double> coeffs = {{"rho",rho},{"mu",mu},{"lam",lam},{"dt",dt},
                                                  {"penalty",penalty_coeff}};
Coefficients coefficients(coeffs, fx, fy);
Point2D p1,p2;

p1 << 0.0,0.0;
p2 << 1.0,1.0;

std::vector<double> L2errors;
std::vector<double> H1errors;

std::vector<double> diameters;

for (std::size_t j = 0; j < n_repetitions; j++) {

  Mesh Th(p1,p2,ne_s[j],1,"Polygon");

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
    problem.solve(solution, u1old, u2old);
    problem.update_rhs(time);
    u2old = u1old;
    u1old = solution;
  }

  double L2_error = problem.L2_error(solution, exact_sol_x, exact_sol_y, time);
  double H1_error = problem.H1_error(solution, grad_1_x, grad_1_y, grad_2_x, grad_2_y, time);

  L2errors.push_back(L2_error);
  H1errors.push_back(H1_error);

  diameters.push_back(Th.get_diameter());
  exportfunctionvtk(femregion,solution,"Displacement" + std::to_string(j));

}

std::vector<double> p;
std::vector<double> pp;


std::cout << "Ne   Diameter   L2error    H1error" <<std::endl;
for (std::size_t k = 0; k < n_repetitions;k++)
  std::cout << ne_s[k] << "  " << diameters[k] << "  " << L2errors[k]  << "  " << H1errors[k]<< std::endl;

for (std::size_t k = 0; k < n_repetitions-1;k++) {
  p.push_back(std::log(L2errors[k]/L2errors[k+1])/std::log(diameters[k]/diameters[k+1]));
  pp.push_back(std::log(H1errors[k]/H1errors[k+1])/std::log(diameters[k]/diameters[k+1]));
}

std::cout << std::endl;
std::cout << "L2 convergence rate" << std::endl;
for (std::size_t k = 0; k < n_repetitions-1;k++)
  std::cout << p[k] << std::endl;

std::cout << std::endl;
std::cout << "H1 convergence rate" << std::endl;
for (std::size_t k = 0; k < n_repetitions-1;k++)
  std::cout << pp[k] << std::endl;

return 0;

}
