/*!
    @file   Legendre.hpp
    @author Nicola Melas
    @brief  Definition of Legendre polynomials and their derivative
*/


#ifndef LEGENDRE_P
#define LEGENDRE_P

#include "types.hpp"
#include "Polygon.hpp"
#include "Triangle.hpp"

namespace PolygonDG {
  /*!
      @brief Evaulate a Legendre polynomial
      This function evaulates the Legendre polynomial of degree N at the points x
      belonging to [-1, 1].
  */
  VectorD LegendreP(const VectorD & x, std::size_t N);

  /*!
      @brief Evaulate the derivative of a Legendre polynomial
      This function evaulates the derivative of a Legendre polynomial of
      degree N at the points x belonging to [-1, 1].
  */
  VectorD GradLegendreP(const VectorD & x, std::size_t N);


  //! Evaluate 2D polynomial basis
  VectorD Basis2DP(const VectorD & a, const VectorD & beta, std::size_t i, std::size_t j);


  /*! Return the derivatives of the modal basis (id,jd)*/
  std::pair<VectorD,VectorD> GradBasis2DP(const VectorD & a, const VectorD & b, std::size_t id, std::size_t jd);


  /*! Return the values of the Basis2DP at the given Nodes and given the polynomial degree
  @param polygon    The polygon on which the basis has to be computed.
  @param Nodes            The points where the basis and gradients are computed
  @param fem_degree       Finite element degree
  Note that physical points are mapped back to compute values.*/
  MatrixXd BasisNodes(const std::shared_ptr<AbstractPolygon> polygon,const Points & Nodes, std::size_t fem_degree);


  /*! Return the values of the GradBasis2DP at the given Nodes and given the polynomial degree
  @param polygon    The polygon on which the basis has to be computed
  @param Nodes            The points where the basis and gradients are computed @n
  @param fem_degree       Finite element degree
  Note that physical points are mapped back to compute values.*/
  VecMatrixXd GradBasisNodes(const std::shared_ptr<AbstractPolygon> polygon,const Points & Nodes, std::size_t fem_degree);

}
#endif
