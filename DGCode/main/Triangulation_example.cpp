/*!
    @file   Triangulation_example.cpp
    @author Nicola Melas
    @brief  Check the CGAL function to triangulate a polygon.
*/

#include "Polygon.hpp"
#include <vector>
#include "Triangulation.hpp"

using namespace PolygonDG;

/*!
  @brief A Polygon can be subdivided in triangles. To do it the CGAL library is
  used. This example verify the output of calling this type of function.
  */

int main(){


  // Define the points
  Point2D p1,p2,p3,p4,p5;
  p1 << 0,0;
  p2 << 1,0;
  p3 << 2,1;
  p4 << 2,2;
  p5 << 1.5,3;
  Points ppp = {p1,p2,p3,p4,p5};

  //build the polygon
  auto point = std::make_shared<Polygon>(ppp);

  //output
  point->showMe();
  std::cout << std::endl;

  //Build the triangulation
  auto triangul = Triangulation(point);


  //Output
  std::cout << "Number of subtriangles: " << triangul.size()<< "." << std::endl;

  std::cout << "The elements: " << std::endl;
  std::cout << std::endl;

  for (auto it : triangul) {
    it->showMe();
    std::cout << std::endl;
  }


return 0;



}
