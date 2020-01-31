/*!
    @file   AbstractPolygon.hpp
    @author Nicola Melas
    @brief  Class that defines the abstract class of polygons
*/

#ifndef HH_POLYGON_HH
#define HH_POLYGON_HH
#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <memory>
#include <Eigen/Dense>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/foreach.hpp>
#include "Edge.hpp"
#include "types.hpp"



namespace PolygonDG
{
  /*!
      @brief Class that defines the abstract class of polygons
      This pure virtual base class defines the properties shared by all polygons.@n
      It is represented by the protected member vertexes that contains the points of the polygon.
      There is also a function that assembles the points in edges and return them.@n
      For example the first edge has the vertex 0 and the vertex 1, the last edge has
      the last vertex and the vertex 0. It is pure virtual due to Specialised versions
      of some functions, in particular how area is compute and the addvertex function.

  */

  class AbstractPolygon
  {
  public:

    //! Constructor taking vertices as Points
    /*!
      It checks convexity if check=true
     */
    AbstractPolygon(Points const & v, bool check=true);

    //! Default constructor is defaulted

    AbstractPolygon()=default;

    //! Default Assignment operator
    AbstractPolygon & operator=(AbstractPolygon const&)=default;

    //! Default Copy constructor
    AbstractPolygon(AbstractPolygon const &)=default;

    //! Default Move constructor
    AbstractPolygon(AbstractPolygon&&)=default;

    //! Default Move assignement
    AbstractPolygon & operator=(AbstractPolygon&&)=default;

    //! Defautl virtual destructor
    virtual ~AbstractPolygon()=default;

    //! Is the polygon convex
    bool isConvex() const {return isconvex;}

    //! Returns the number of vertices.
    std::size_t size() const {return vertexes.size();}

    //! Intersection between polygons
    /*!
      It relies on the boost library. @n
      The @param polygon is passed through pointer to handle correctly the children class
    */
    std::vector<std::shared_ptr<AbstractPolygon>> polyintersection(const std::shared_ptr<AbstractPolygon> polygon) const ;

    //! Returns the vertices (read only)
    Points const & theVertexes()const {return vertexes;}

    //! Outputs some info on the polygon
    virtual void showMe(std::ostream & out=std::cout) const;

    //! Compute Centroid of polygon
    /*! It relies on the boost library */
    Point2D compute_centroid() const;

    //! Returns the minimum of the x-component of the vertices
    double x_min() const;

    //! Returns the maximum of the x-component of the vertices
    double x_max() const;

    //! Returns the minimum of the y-component of the vertices
    double y_min() const;

    //! Returns the maximum of the y-component of the vertices
    double y_max() const;

    //! Returns the vector of all the x-components
    std::vector<double> get_x() const;

    //! Returns the vector of all the y-components
    std::vector<double> get_y() const;

    /*!
      @brief Get a Vertex
      This functions returns the i-th Vertex.
      @param i The index of the Vertex required.
  */
    Point2D get_vertex(std::size_t i) const {return vertexes[i];};

    //! Return vector of Edge
    /*! The first edge has vertex 0 and vertex 1. @n
    The last one has vertex size()-1 and vertex 0. */
    std::vector<Edge> get_edges() const;

    //! Return the normals(external) to the Edges
    Points get_normals() const;

    //! Return the size of each edge of the polygon
    std::vector<double> get_edges_length() const;

    //! Return the Bounding Box
    /*! It's the smallest rectangle containing the polygon.*/
    std::vector<double> get_BBox() const ;

    //! Return max_kb
    /*! This variable contains, for each edge (following the vertexes order)
    the maximum area of a triangle containing the edge and one of the other
    vertexes of the polygon. @n
    It's pure virtual because there could be an easy version
    in some special cases.*/
    virtual std::vector<double> get_max_kb() const = 0;

    //! The area of the polygon (with sign!).
    virtual double area() const=0;

    //! Add a vertex to the polygon
    /*! It's pure virtual, the children classes has to eventually check if it's
    ok or not to add a vertex.*/
    virtual void addvertex(const Point2D & vertex) = 0;

  protected:
    //! Vertexes of the polygon
    Points vertexes;

    //! True if polygon is convex
    bool isconvex;

    //! Test convexity of the polygon
    void checkConvexity();

  };


}

#endif
