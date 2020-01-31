/*!
    @file   Triangle.cpp
    @author Nicola Melas
    @brief  Implementation of class Triangle
*/


#include "Triangle.hpp"


namespace PolygonDG {

  Triangle::Triangle(Point2D const & v1,Point2D const & v2,Point2D const & v3) {

    this->isconvex=true;
    this->addvertex(v1);
    this->addvertex(v2);
    this->addvertex(v3);

  }
  Triangle::Triangle(Points const & v):AbstractPolygon(v,false){
    this->isconvex=true;
    // Check if we give 3 vertices
    if(this->size() != 3){
      throw std::runtime_error(" A triangle must be created giving three vertices");
    }
  }


  double Triangle::area() const
  {
    if(this->size()==0) return 0.0;
    // I use the cross product since this is a
    // signed area!
    Point2D v(this->vertexes[1]-this->vertexes[0]);
    Point2D w(this->vertexes[2]-this->vertexes[0]);
    // area = 0.5* v \times w. Positive if triangle counterclockwise oriented
    return 0.5*(v(0)*w(1)-v(1)*w(0));
    ;}

  void Triangle::showMe(std::ostream & out) const{
    out<<"A Triangle"<<std::endl;
    AbstractPolygon::showMe(out);
  }

  int Triangle::orientation() const
  {
    const auto & p1 = vertexes[0];
    const auto & p2 = vertexes[1];
    const auto & p3 = vertexes[2];

    double val = (p2(1) - p1(1)) * (p3(0) - p2(0)) -  (p2(0) - p1(0)) * (p3(1) - p2(1));

    if (val == 0) return 0;  // colinear

    return (val > 0)? 1: 2; // clock or counterclock wise
  }

  bool Triangle::clockwise() {
    int orientation = this->orientation();
    if (orientation == 0)
      return 0;
    else if (orientation ==2){
      std::swap(vertexes[0],vertexes[2]);
      return 1;
    }
    return 1;
  }

  bool Triangle::counterclockwise() {
    int orientation = this->orientation();
    if (orientation == 0)
      return 0;
    else if (orientation ==1){
      std::swap(vertexes[0],vertexes[2]);
      return 1;
    }
    return 1;
  }

  std::vector<double> Triangle::get_max_kb() const {
    std::vector<double> max_kb(3);
    max_kb = {std::abs(area()), std::abs(area()), std::abs(area())};
    return max_kb;
  }

}
