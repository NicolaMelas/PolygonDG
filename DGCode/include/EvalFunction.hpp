/*!
    @file   EvalFunction.hpp
    @author Nicola Melas
    @brief  External function to evaluate e function on a fespace
*/

#ifndef EVALFUNCTION
#define EVALFUNCTION

#include "fespace.hpp"
#include "types.hpp"


namespace PolygonDG {
  //! Given the 2D function and the time, it returns the evaluation of the function on the FeSpace dofs
  /*! The function is computed on the quadrature nodes and then in
   the return vector there are the values of the function on the
   degrees of freedom that depends on the fem_degree of the
   given FeSpace.

  @param femregion  The FeSpace
  @param u_x        First component of the function that has to be evaluated
  @param u_y        Second component of the function that has to be evaluated
  @param time       Time, third component of the function
  */

  VectorD evaluate_function(const FeSpace & femregion, MuParserInterface::muParserXInterface<3, std::array<double, 3>> & u_x,  MuParserInterface::muParserXInterface<3, std::array<double, 3>> & u_y, const double & time) ;

}

#endif
