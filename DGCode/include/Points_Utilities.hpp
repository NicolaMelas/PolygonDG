/*!
    @file   Points_Utilities.hpp
    @author Nicola Melas
    @brief  Some utilities function fot class Points
*/

#ifndef PTS_UTILITIES
#define PTS_UTILITIES

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "types.hpp"

namespace PolygonDG {

  /*@brief  Some utilities function fot class Points
    They are used to simplify the readability of the code.
    Be carefull to understand the difference between the operator *, used
    to return the sum of all the scalar products of all Point2D
    and the cwise product that perform only the cwise
    product of each Point2D*/

  //! Return the vector of all the first component of the 2d Points
  std::vector<double> get_x(const Points &);

  //! Return the vector of all the second component of the 2d Points
  std::vector<double> get_y(const Points &);

  //! The sum is computed for each point in points with p2d
  Points operator+(const Points & points, const Point2D & p2d);

  //! Scalar product between points
  double operator*(const Points & p1, const Points & p2);

  //! Component-wise product between vector of points
  Points cwiseprod(const Points & p1, const Points & p2);

}




#endif
