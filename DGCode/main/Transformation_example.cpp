/*!
    @file   Transformation_example.cpp
    @author Nicola Melas
    @brief  Check the map from the reference element.
*/


#include <iostream>
#include "Mesh.hpp"
#include "Points_Utilities.hpp"
#include "muParserXInterface.hpp"
#include <vector>
#include <string>
#include "Transformation.hpp"
#include "Triangle.hpp"
#include "Polygon.hpp"
#include <fstream>
#include "GetPot.hpp"


using namespace PolygonDG;

/*!
  @brief It checks the transformation from reference triangle/square to the
  physical points. Done both for triangles and polygons.
  */

int main() {

  Point2D p1,p2,p3;

  p1 << 10.0, 3.0;
  p2 << 31.54, 52.74;
  p3 << 20,70.56;

  auto tria = std::make_shared<Triangle>(p1,p2,p3);
  auto BJ_trans = get_jacobian_and_translation(tria);

  MatrixXd BJ = BJ_trans.first;
  VectorD translation = BJ_trans.second;

  Point2D result, partial;
  tria->showMe();

  std::cout << std::endl;
  partial << 0,0;
  result = BJ *  partial + translation;
  std::cout << BJ(0,0) << "  " << BJ(0,1) << "   *   " << partial(0) << "   +   " << translation(0) << "   =   " << result(0) << std::endl;
  std::cout << BJ(1,0) << "  " << BJ(1,1) << "       " << partial(1) << "   +   " << translation(1) << "   =   " << result(1) << std::endl;

  std::cout << std::endl;
  partial << 1,0;
  result = BJ *  partial + translation;
  std::cout << BJ(0,0) << "  " << BJ(0,1) << "   *   " << partial(0) << "   +   " << translation(0) << "   =   " << result(0) << std::endl;
  std::cout << BJ(1,0) << "  " << BJ(1,1) << "       " << partial(1) << "   +   " << translation(1) << "   =   " << result(1) << std::endl;

  std::cout << std::endl;
  partial << 0,1;
  result = BJ *  partial + translation;
  std::cout << BJ(0,0) << "  " << BJ(0,1) << "   *   " << partial(0) << "   +   " << translation(0) << "   =   " << result(0) << std::endl;
  std::cout << BJ(1,0) << "  " << BJ(1,1) << "       " << partial(1) << "   +   " << translation(1) << "   =   " << result(1) << std::endl;

  Point2D p4,p5; ;

  p4 << 15.0, 50.0;
  p5 << -5, 20;

  Points vertexes = {p1,p2,p3,p4,p5};
  auto poly = std::make_shared<Polygon>(vertexes);
  BJ_trans = get_jacobian_and_translation(poly);

  std::cout << std::endl;
  poly->showMe();

  auto BBox = poly->get_BBox();
  std::cout << std::endl;
  std::cout << "Polygon BBox" << std::endl;
  std::cout << BBox[0] << "  " << BBox[2] << std::endl;
  std::cout << BBox[1] << "  " << BBox[2] << std::endl;
  std::cout << BBox[0] << "  " << BBox[3] << std::endl;
  std::cout << BBox[1] << "  " << BBox[3] << std::endl;

  BJ = BJ_trans.first;
  translation = BJ_trans.second;

  std::cout << std::endl;
  partial << -1,-1;
  result = BJ *  partial + translation;
  std::cout << BJ(0,0) << "  " << BJ(0,1) << "   *   " << partial(0) << "   +   " << translation(0) << "   =   " << result(0) << std::endl;
  std::cout << BJ(1,0) << "  " << BJ(1,1) << "       " << partial(1) << "   +   " << translation(1) << "   =   " << result(1) << std::endl;

  std::cout << std::endl;
  partial << 1,-1;
  result = BJ *  partial + translation;
  std::cout << BJ(0,0) << "  " << BJ(0,1) << "   *   " << partial(0) << "   +   " << translation(0) << "   =   " << result(0) << std::endl;
  std::cout << BJ(1,0) << "  " << BJ(1,1) << "       " << partial(1) << "   +   " << translation(1) << "   =   " << result(1) << std::endl;

  std::cout << std::endl;
  partial << -1, 1;
  result = BJ *  partial + translation;
  std::cout << BJ(0,0) << "  " << BJ(0,1) << "   *   " << partial(0) << "   +   " << translation(0) << "   =   " << result(0) << std::endl;
  std::cout << BJ(1,0) << "  " << BJ(1,1) << "       " << partial(1) << "   +   " << translation(1) << "   =   " << result(1) << std::endl;

  std::cout << std::endl;
  partial << 1,1;
  result = BJ *  partial + translation;
  std::cout << BJ(0,0) << "  " << BJ(0,1) << "   *   " << partial(0) << "   +   " << translation(0) << "   =   " << result(0) << std::endl;
  std::cout << BJ(1,0) << "  " << BJ(1,1) << "       " << partial(1) << "   +   " << translation(1) << "   =   " << result(1) << std::endl;


  return 0;
}
