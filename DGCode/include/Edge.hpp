/*!
    @file   Edge.hpp
    @author Nicola Melas
    @brief  Class that defines the Edge as segment of two point
*/

#ifndef EDGE_H
#define EDGE_H

#include <array>
#include <vector>
#include "Coefficients.hpp"
#include "types.hpp"

namespace PolygonDG {
  /*!
      @brief Class that defines the Edge as segment of two point. @n
      It's a simple class. The edge is composed by two point. The order of
      the point is important. For example the function that return the normal
      to the edge change orientation if the vertexes are stored in the
      reverse order.
  */
  class Edge {

  private:
    std::array<Point2D,2> edge;

  public:
    //!Constructor from two points, p1 and p2
    Edge(const Point2D & p1, const Point2D & p2) {edge = {{p1,p2}};};

    //!Constructor from array
    Edge(const std::array<Point2D,2> & edge_) : edge(edge_){};

    //! Default Assignment operator
    Edge & operator=(Edge const&)=default;

    //! Default Copy constructor
    Edge(Edge const &)=default;

    //! Default Move constructor
    Edge(Edge&&)=default;

    //! Default Move assignement
    Edge & operator=(Edge&&)=default;

    //! Return the edge points
    inline const std::array<Point2D,2> & get_points() const {return edge;};

    //! return the first point
    inline const Point2D & first() const {return edge[0];};

    //! return the second point
    inline const Point2D & second() const {return edge[1];};

    //! Return the normal vector to the edge (w.r.t. the order the points are stored)
    Point2D get_normal() const;

    //! Distance between the two vectors
    double get_size() const;

  };

  //! Returns the distance between two points
  double distance(const Point2D & p1, const Point2D & p2);
}

#endif
