/*!
    @file   Points_Utilities.cpp
    @author Nicola Melas
    @brief  Implementation of some utility functions
*/

#include "Points_Utilities.hpp"


namespace PolygonDG {

  std::vector<double> get_x(const Points & p) {

    std::vector<double> xx(p.size());

    for (std::size_t j = 0; j < xx.size(); j++)
      xx[j] = p[j](0);

    return xx;

  }


  std::vector<double> get_y(const Points & p) {

    std::vector<double> xx(p.size());

    for (std::size_t j = 0; j < xx.size(); j++)
      xx[j] = p[j](1);

    return xx;

  }

  Points operator+(const Points & points, const Point2D & p2d)  {
    Points P = points;
    for (std::size_t jj = 0; jj < points.size(); jj++){
      P[jj] += p2d;
    }
    return P;
  }


  double operator*(const Points & p1, const Points & p2) {
    double result = 0;
    for (std::size_t jj = 0; jj < p1.size(); jj++){
      result += (p1[jj]).cwiseProduct(p2[jj]).sum();
    }
    return result;
  }

  Points cwiseprod(const Points & p1, const Points & p2) {
    Points result(p1.size());
    for (std::size_t k = 0; k < p1.size(); k++)
      result[k] = (p1[k].cwiseProduct(p2[k]));

    return result;
  }
}
