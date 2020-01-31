/*!
    @file   Triangulation.cpp
    @author Nicola Melas
    @brief  Implementation of triangulation routine
*/

#include "Triangulation.hpp"
#include <iterator>

namespace PolygonDG {



  std::vector<std::shared_ptr<Triangle>> Triangulation(const std::shared_ptr<AbstractPolygon> polygon) {

    //converting for boost library
    Points loc_coord = polygon->theVertexes();
    std::vector<Point_3> polyline;

    for (std::size_t ii = 0; ii < loc_coord.size(); ii++)
      polyline.push_back(Point_3(loc_coord[ii](0),loc_coord[ii](1),0.));

    std::vector<Triangle_int> patch;
    std::size_t n_tria = polyline.size() -2;
    patch.reserve(n_tria); // there will be exactly n-2 triangles in the patch
    //trinagulation boost routine
    CGAL::Polygon_mesh_processing::triangulate_hole_polyline(polyline,std::back_inserter(patch));

    std::vector<std::shared_ptr<Triangle>> result;
    result.reserve(n_tria);
    //return Triangles
    for (std::size_t ii = 0; ii < n_tria; ii++)
      result.push_back(std::make_shared<Triangle>(loc_coord[patch[ii].second],
                                                  loc_coord[patch[ii].third],
                                                  loc_coord[patch[ii].first]));


  return result;

  }


}
