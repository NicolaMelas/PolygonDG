/*!
    @file   TriaGenerator.cpp
    @author Nicola Melas
    @brief  Implementation of class TriaGenerator
*/


#include "TriaGenerator.hpp"

namespace PolygonDG {

  void TriaGenerator::generate(Mesh & mesh) {

    mesh.ne = ne;

    double dx = (domain[1]-domain[0])/ne_x;
    double dy = (domain[1]-domain[0])/ne_y;

    for(std::size_t jj = 0; jj < ne_y + 1; jj++){
      for(std::size_t ii = 0; ii < ne_x + 1; ii++){
        //building points
        Point2D tmp;
        tmp << ii*dx, jj*dy;
        mesh.coord.push_back(tmp);

        if (ii != ne_x && jj != ne_y) {
          //building connectivity
          VectorI conn1(3), conn2(3);
          conn1 << jj*(ne_x + 1) + ii, (jj+1)*(ne_x + 1) + ii + 1, (jj+1)*(ne_x + 1) + ii;
          conn2 << jj*(ne_x + 1) + ii, jj*(ne_x + 1) + ii + 1, (jj+1)*(ne_x + 1) + ii + 1;
          mesh.connectivity.push_back(conn1);
          mesh.connectivity.push_back(conn2);
        }
      }
    }


    for (std::size_t kk = 0; kk < mesh.connectivity.size(); kk++){
      //building triangles
      Point2D v1,v2,v3;
      v1 = mesh.coord[mesh.connectivity[kk](0)];
      v2 = mesh.coord[mesh.connectivity[kk](1)];
      v3 = mesh.coord[mesh.connectivity[kk](2)];

      Triangle tria(v1,v2,v3);
      mesh.polygons.push_back(std::make_shared<Triangle>(tria));

    }
  }

}
