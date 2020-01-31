/*!
    @file   EvalFunction.cpp
    @author Nicola Melas
    @brief  Implementation of function to evaluate a function
*/

#include "EvalFunction.hpp"
#include <chrono>

namespace PolygonDG {

  VectorD evaluate_function(const FeSpace & femregion, MuParserInterface::muParserXInterface<3, std::array<double, 3>> & u_x,  MuParserInterface::muParserXInterface<3, std::array<double, 3>> & u_y, const double & time = 0.0) {

    std::size_t ndof = femregion.get_ndof();
    std::size_t nln = femregion.get_nln();

    VectorD u1 = VectorD::Zero(ndof);
    VectorD u2 = VectorD::Zero(ndof);

    VectorSt partial_index(nln);
    VectorSt index(nln);
    VectorSt tmp(nln);

    for (std::size_t ii = 0; ii < nln; ii++)
        tmp(ii) = ii;

    //Reference to all needed info
    const std::vector<FeElement> & elements = femregion.get_FeElements();
    std::size_t ne = elements.size();

    for (std::size_t ie = 0; ie < ne; ie++) {

      const FeElement & element = elements[ie];

      partial_index = VectorSt::Constant(nln, (ie)*nln);
      index = partial_index + tmp;
      std::size_t n_tria = element.get_number_tria();

      for (std::size_t iTria = 0; iTria < n_tria; iTria++) {

        const Points & phys_tria = element.get_physical_points_tria()[iTria];
        const auto & dx_tria = element.get_dx()[iTria];
        std::size_t dx_size = dx_tria.size();
        const MatrixXd & phi_tria = element.get_phi_tria()[iTria];

        for (std::size_t k = 0; k < dx_size; k++) {

          //scaling quadrature coefficient
          const auto & dx = dx_tria[k];

          const auto & x =  (phys_tria)[k](0);
          const auto & y =  (phys_tria)[k](1);

          //muparser evaluation
          std::array<double,3> xx;
          xx[0] = x;
          xx[1] = y;
          xx[2] = time;

          double U1 = u_x(xx);
          double U2 = u_y(xx);

          //fill the vector with i-th phi at the given quadrature nodes
          for (std::size_t i = 0; i < nln; i++) {
            u1(index(i)) += U1 * phi_tria(k,i) * dx;
            u2(index(i)) += U2 * phi_tria(k,i) * dx;
          }
        }
      }
    }
    VectorD u = VectorD::Zero(2*ndof);
    u.block(0,0,ndof,1) = u1;
    u.block(ndof,0,ndof,1) = u2;


    return u;
  }


}
