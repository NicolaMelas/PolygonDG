/*!
    @file   AbstractPolygon.cpp
    @author Nicola Melas
    @brief  Implementation of class AbstractPolygons
*/

#include "Polygon.hpp"
#include <limits>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/foreach.hpp>
#include "Points_Utilities.hpp"

namespace PolygonDG{

  AbstractPolygon::AbstractPolygon(Points const & v, bool check):vertexes(v)
  {
    if (check) this->checkConvexity();
  }

  void AbstractPolygon::showMe(std::ostream & out)const
  {
    if (this->size()==0){
      out<< "empty polygon" << std::endl;
    }
    else {
      out<<"Points:"<<std::endl;
      out<<"    X     " << "   Y    "<<std::endl;
      for (auto const & i : this->vertexes)
        out<<std::setprecision(4) << i(0)<<" " << i(1)<<std::endl;
    }
  }

  void AbstractPolygon::checkConvexity()
  {
    Points const & myV(this->vertexes);
    auto mysize=this->size();
    // We consider segments and triangles as convex
    if (mysize <= 3) {
      this->isconvex=true;
      return;
    }
    //! Since we are dealing with floating points it is better to have
    //  a small number so that |a| < smallNumber means for us a==0
    double const smallNumber(1000*std::numeric_limits<double>::min());
    Point2D p;
    Point2D v;
    Point2D u;
    double res(0.0);
    double newres(0.0);
    //! C++11 sintax. decltype(expr) returns the type of the expression
    for ( decltype(mysize) i=0; i < mysize; ++i)  {
      p = myV[i];
      // ! next point
      Point2D tmp = myV[(i+1) % myV.size()];
      v = tmp - p;
      //! next next point
      u = myV[(i+2) % myV.size()];
      if (i == 0) // in first loop direction is unknown, so save it in res
        res = u(0) * v(1) - u(1) * v(0) + v(0) * p(1) - v(1) * p(0);
	    else{
        newres = u(0) * v(1) - u(1) * v(0) + v(0) * p(1) - v(1) * p(0);
        if (std::abs(res)<smallNumber){
        // The two edges are aligned, skip test and update res
          res=newres;
        }
        else if ((newres > 0 && res < 0) || (newres < 0 && res > 0) ) {
          this->isconvex=false;
          return;
        }
      }
    }// end for
    this->isconvex=true;
    return;
  }

  Point2D AbstractPolygon::compute_centroid() const {

    // Converting Polygon class in the boost polygon class
    namespace bg = boost::geometry;
    typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
    typedef bg::model::polygon<point_t> polygon_t;

    polygon_t poly;

    // Fill the Outer Ring
    for (std::size_t kk = 0; kk < vertexes.size();kk++)
      bg::append(poly.outer(),point_t(vertexes[kk](0),vertexes[kk](1)));
    bg::append(poly.outer(),point_t(vertexes[0](0),vertexes[0](1)));

    point_t pp(0.0,0.0);
    // Computing the centroid of poly inside pp
    bg::centroid(poly,pp);

    // Converting Boost sintax to Eigen sintax
    Point2D centroid;
    centroid << pp.get<0>(), pp.get<1>();
    return centroid;

  }

  double AbstractPolygon::x_min() const {

    double result = vertexes[0](0);
    std::size_t psize = vertexes.size();

    for (std::size_t jj = 1; jj < psize; jj++) {
      if (vertexes[jj](0) < result)
        result = vertexes[jj](0);
    }
    return result;
  }

  double AbstractPolygon::y_min() const {

    double result = vertexes[0](1);
    std::size_t psize = vertexes.size();

    for (std::size_t jj = 1; jj < psize; jj++) {
      if (vertexes[jj](1) < result)
        result = vertexes[jj](1);
    }
    return result;
  }

  double AbstractPolygon::x_max() const {

    double result = vertexes[0](0);
    std::size_t psize = vertexes.size();

    for (std::size_t jj = 1; jj < psize; jj++) {
      if (vertexes[jj](0) > result)
        result = vertexes[jj](0);
    }
    return result;
  }

  double AbstractPolygon::y_max() const {
    double result = vertexes[0](1);
    std::size_t psize = vertexes.size();

    for (std::size_t jj = 1; jj < psize; jj++) {
      if (vertexes[jj](1) > result)
        result = vertexes[jj](1);
    }
    return result;
  }

  std::vector<double> AbstractPolygon::get_x() const {

    return PolygonDG::get_x(vertexes);
  }

  std::vector<double> AbstractPolygon::get_y() const {

    return PolygonDG::get_y(vertexes);

  }

  std::vector<Edge> AbstractPolygon::get_edges() const {

    std::vector<Edge> edges;
    std::size_t nv = vertexes.size();
    edges.reserve(nv);

    for (std::size_t i = 0; i < nv; i++ )
      edges.push_back(Edge(vertexes[i],vertexes[(i+1)%nv]));

    return edges;

  }

  Points AbstractPolygon::get_normals() const {

    std::vector<Edge> loc_edges = get_edges();
    Points normal;
    normal.reserve(loc_edges.size());

    for (const Edge & ed : loc_edges)
      normal.push_back(ed.get_normal());

    return normal;
  }

  std::vector<double> AbstractPolygon::get_edges_length() const{

    std::vector<Edge> loc_edges = get_edges();
    std::vector<double> meshsize;
    meshsize.reserve(loc_edges.size());

    for (const Edge & ed : loc_edges)
      meshsize.push_back(ed.get_size());

    return meshsize;

  }

  std::vector<std::shared_ptr<AbstractPolygon>> AbstractPolygon::polyintersection(const std::shared_ptr<AbstractPolygon> polygon) const {

    //Converting the polygons to boost polygon sintax
    namespace bg = boost::geometry;
    typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
    typedef bg::model::polygon<point_t> polygon_t;

    auto verticesp1 = this->theVertexes();
    auto verticesp2 = polygon->theVertexes();
    std::size_t p1size = verticesp1.size();
    std::size_t p2size = verticesp2.size();

    polygon_t poly1, poly2;

    // Fill outer rings
    for (std::size_t kk = 0; kk < p1size;kk++)
      bg::append(poly1.outer(),point_t(verticesp1[kk](0),verticesp1[kk](1)));
    bg::append(poly1.outer(),point_t(verticesp1[0](0),verticesp1[0](1)));

    for (std::size_t kk = 0; kk < p2size;kk++)
      bg::append(poly2.outer(),point_t(verticesp2[kk](0),verticesp2[kk](1)));
    bg::append(poly2.outer(),point_t(verticesp2[0](0),verticesp2[0](1)));

    // initialization of the output
    std::vector<polygon_t> output;
    bg::intersection(poly1, poly2, output);

    // Converting to Mypolygon class
    std::vector<std::shared_ptr<AbstractPolygon>> out;
    out.reserve(output.size());

    for (std::size_t jj = 0; jj < output.size(); jj++) {
      auto tmppoly = std::make_shared<Polygon>();
      for(auto it = boost::begin(boost::geometry::exterior_ring(output[jj])); it != boost::end(boost::geometry::exterior_ring(output[jj])); ++it)
      {
        double x = bg::get<0>(*it);
        double y = bg::get<1>(*it);
        Point2D tmp;
        tmp << x,y;
        (*tmppoly).addvertex(tmp);
      }
      out.push_back(tmppoly);
    }
  return out;
  }

  std::vector<double> AbstractPolygon::get_BBox() const {
    std::vector<double> BBox(4);
    BBox = {x_min(),x_max(),y_min(),y_max()};
    return BBox;
  }

}
