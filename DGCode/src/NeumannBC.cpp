/*!
    @file   NeumannBC.cpp
    @author Nicola Melas
    @brief  Implementation of class NeumannBC
*/

#include "NeumannBC.hpp"

namespace PolygonDG {


  void NeumannBC::apply(VectorD &rhs, const FeSpace & femregion, double time) const {

    // Getting by const reference things needed to apply BC
    const Neighbour & neighbours = femregion.get_Neighbour();
    const std::vector<FeElement> & FeElements = femregion.get_FeElements();
    const VecVecXi & ntmp = neighbours.get_neigh();

    std::size_t ndof = femregion.get_ndof();
    std::size_t nln = femregion.get_nln();

    // Vector to save degrees of freedom indexes
    VectorSt index(nln);
    VectorSt tmp(nln);

    for (std::size_t ii = 0; ii < nln; ii++)
      tmp[ii] = ii;

    for (std::size_t ie = 0; ie < ntmp.size(); ie++) {

      const FeElement & element = FeElements[ie];
      VectorSt partial_index = VectorSt::Constant(nln, (ie)*nln);
      index = partial_index + tmp;

      const auto & neigh_tmp = ntmp[ie];
      std::size_t size_tmp = neigh_tmp.size();

      // Saving scaled edge length, basis on the edge and physical points on the edge
      const std::vector<std::vector<double>> & ds_edges = element.get_ds();
      const std::vector<MatrixXd> & Phi = element.get_phi_edge();
      const std::vector<Points> & phys_edges = element.get_physical_points_edge();

      for (std::size_t iedg = 0; iedg < size_tmp ; iedg++) {

        //BC applied only where there is the tag of BC itself
        if (neigh_tmp(iedg) == tag) {

          const std::vector<double> & ds_edge = ds_edges[iedg];
          std::size_t ds_size = ds_edge.size();

          const Points & phys_edge = phys_edges[iedg];
          for (std::size_t k = 0; k < ds_size; k++){
            auto ds = ds_edge[k];

            MatrixXd Bedge(nln,1);

            for (std::size_t kk = 0; kk < nln; kk++) {
              Bedge(kk) = Phi[iedg](k,kk);
            }

            //muparser evaluation
            double x = phys_edge[k](0);
            double y = phys_edge[k](1);

            std::array<double,3> xx;
            xx[0] = x;
            xx[1] = y;
            xx[2] = time;

            double gd1 = g1(xx);
            double gd2 = g2(xx);

            for (std::size_t i=0; i < nln;i++) { // loop over scalar shape functions

              rhs(index(i)) += Bedge(i) * gd1 * ds;
              rhs(ndof+index(i)) += Bedge(i) * gd2 * ds;

            }
          }
        }
      }
    }
  }

}
