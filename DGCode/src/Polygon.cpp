/*!
    @file   Polygon.cpp
    @author Nicola Melas
    @brief  Implementation of class Polygon
*/

#include "Polygon.hpp"
#include "Triangle.hpp"
#include "Points_Utilities.hpp"

namespace PolygonDG {



  Polygon::Polygon(Points const & v): AbstractPolygon(v) {}


  //! To compute the area of a polygon the divergence theorem is used.
  /*!
    \f$ int_P d\Omega = 1/2  \int{\partial P} xn_x d\gamma\f$
    The integral is computed by using trapezoidal rule on the polygon sides.
  */

  double Polygon::area() const {

    auto siz=this->size();
    if (siz<3) return 0;
    double result(0);
    Points const & Ver(this->vertexes);

    for (std::size_t i=0; i<siz;++i){

      Point2D const & p1(Ver[i]);
      // Other point
      Point2D const & p2(Ver[(i+1) % siz]);
      Point2D const & p0(Ver[(i-1) % siz]);
      result+=p1.x()*(p2.y()-p0.y());
    }

    return 0.5*result;

  }

  void Polygon::showMe(std::ostream & out)const
  {
    std::cout<<" A Generic Polygon"<<std::endl;
    AbstractPolygon::showMe(out);
  }


  std::vector<double> Polygon::get_max_kb() const {

    std::size_t n_edge = size();
    std::vector<double> max_kb(n_edge,0);

    for (std::size_t j = 0; j < n_edge; j++){
      //Fixed two vertices
      const Point2D & v1 = vertexes[j];;
      const Point2D & v2 = vertexes[(j+1)%n_edge];

      for (std::size_t k = 0; k < n_edge; k++) {
        //I search for the triangle with the biggest area with some area checks
        if (k!=j && k!=j+1){
          const Point2D & v3 = vertexes[k];
          Triangle triangle(Points({v1,v2,v3}));

          bool orient = triangle.clockwise();

          if (orient !=0) {

            double area = triangle.area();

            Polygon elemcoords;

            for (std::size_t bb = 0; bb < n_edge; bb++)
              elemcoords.addvertex(vertexes[n_edge-1-bb]);

            auto x1y1v = elemcoords.polyintersection(std::make_shared<Triangle>(triangle));

            if (x1y1v.size() != 1)
              throw std::runtime_error("There is an error in generating the mesh. The insersection of polygon must return a vector of size 1");
            auto sharedtoraw = x1y1v[0].get();
            Polygon x1y1 = *(dynamic_cast<Polygon *>(sharedtoraw));

            std::vector<double> x1 = x1y1.get_x();
            if (1-(std::find_if(x1.begin(),x1.end(),[](double d) { return std::isnan(d); })!=x1.end()) && std::abs(x1y1.area()-area)<1e-8){
              if (std::abs(area) > max_kb[j])
                max_kb[j] = std::abs(area);
            }
          }
        }
      }
    }
    return max_kb;
  }

}
