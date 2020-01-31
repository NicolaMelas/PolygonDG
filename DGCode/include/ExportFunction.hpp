/*!
    @file   ExportFunction.hpp
    @author Nicola Melas
    @brief  Export function to be read in paraview
*/
#ifndef EXPORTF
#define EXPORTF
#include "fespace.hpp"
#include "types.hpp"

namespace PolygonDG {
  /*!
      @brief  Export function to be read in paraview. The function reconstruct
      the values at the quadrature nodes with corresponding values. Then they
      are written in paraview format.
  */
  void exportfunctionvtk(const FeSpace & femregion, const VectorD & u_h, std::string filename);

  /*! Return the values of the basis at the given Nodes and given the polynomial degree
  @param polygon    The polygon on which the basis has to be computed
  @param Nodes            The points where the basis and gradients are computed @n
  @param fem_degree       Finite element degree
  Note that physical points are mapped back to compute values.*/
  MatrixXd BasisNodes(const std::shared_ptr<AbstractPolygon> polygon, const Points & Nodes, std::size_t fem_degree);

}


#endif
