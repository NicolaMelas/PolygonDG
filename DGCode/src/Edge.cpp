/*!
    @file   Edge.cpp
    @author Nicola Melas
    @brief  Implementation of class Edge
*/

#include "Edge.hpp"

namespace PolygonDG {

  double distance(Point2D const & a, Point2D const & b){

    Point2D tmp(a-b);
    return std::sqrt(tmp(0)*tmp(0)+
		     tmp(1)*tmp(1)
		     );
  }

  Point2D Edge::get_normal() const {

    Point2D nn = VectorD::Zero(2);

    nn << -edge[0](1) +  edge[1](1), edge[0](0) -  edge[1](0);  // corresponding normal

    double nn_norm = nn.norm();
    //normalizing
    nn/=nn_norm;

    return nn;

  }

  double Edge::get_size() const {

    const Point2D & p1 = first();
    const Point2D & p2 = second();

    return distance(p1,p2);

  }
}
