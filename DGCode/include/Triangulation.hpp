/*!
    @file   Triangulation.hpp
    @author Nicola Melas
    @brief  Function that perform a triangulation of a polygon
*/


#ifndef TRIANGULATION
#define TRIANGULATION

#include "AbstractPolygon.hpp"
#include "Triangle.hpp"

namespace PolygonDG {

  /*! Returns the triangulation of the polygon. This function
    relies on the CGAL library. The goal is to subdivide
    a polygon in Triangles in a safe manner. The polygon is passed through
    pointer to handle both cases, triangular and polygonal*/
  std::vector<std::shared_ptr<Triangle>> Triangulation(const std::shared_ptr<AbstractPolygon>);

}

#endif
