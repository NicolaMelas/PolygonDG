/*!
    @file   Transformation.hpp
    @author Nicola Melas
    @brief  Functions to describe tha map between the psysical element and
            the reference element.
*/

#ifndef TRANSFORMATION
#define TRANSFORMATION

#include <utility>
#include "types.hpp"
#include "AbstractPolygon.hpp"
#include "Polygon.hpp"
#include "Triangle.hpp"
#include "Edge.hpp"

namespace PolygonDG {

  /*!
      @brief  Functions to describe tha map between the psysical element and
              the reference element.

      To do this, We need the jacobian of the trasformation such That
      physical = Jacobian* reference + translation. The computation is different
      beetwen triangular and polygonal case.
      This map is computed for the 1D and 2D cases, then used in the numerical
      integration.
  */

  //! Returns the Jacobian and translation vector of the transformation from the reference triangle or Boundingbox
  std::pair<MatrixXd,VectorD> get_jacobian_and_translation(const std::shared_ptr<AbstractPolygon>);

  //! Return the corresponding 1D physical point wrt the quadrature nodes
  Points get_physical_points_faces(const Edge & edge, const std::vector<double> & node_1D);

  //! Return the corresponding 2D physical point wrt the 2D quadrature nodes
  Points get_physical_points(const std::shared_ptr<AbstractPolygon>, const Points & nodes_2D);


}

#endif
