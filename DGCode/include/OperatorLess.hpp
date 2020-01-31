/*!
    @file   OperatorLess.hpp
    @author Nicola Melas
    @brief  Class that defines the Less Operator
*/

#ifndef EIGENLESS
#define EIGENLESS

#include "types.hpp"

namespace PolygonDG {

  /*!
      @brief  Class that defines the Less Operator
      This template function represents the < operator for two Eigen Vector.
      Defined to do a sort operation between eigen type vector
   */
  template <typename T>
  bool Eigenless(const Eigen::Matrix<T,2,1> & a, const Eigen::Matrix<T,2,1> & b) {

    for(int size = 0; size < std::min(a.rows(),b.rows()); size++){
      if (a(size) < b(size))
        return 1;
      else if (a(size) > b(size))
        return 0;
    }
    //to return the shorter
    bool shorter = 0;
    if (a.rows() < b.rows())
      shorter = 1;
    return shorter;
  }


}


#endif
