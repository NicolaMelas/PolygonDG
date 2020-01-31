/*!
    @file   DirichletBC.cpp
    @author Nicola Melas
    @brief  Implementation of class DirichletBC
*/

#include "DirichletBC.hpp"

namespace PolygonDG {

  void DirichletBC::apply(VectorD &rhs, const FeSpace & femregion, double time) const {

    // Getting by const reference things needed to apply BC
    const Neighbour & neighbours = femregion.get_Neighbour();
    const std::vector<FeElement> & FeElements = femregion.get_FeElements();
    const VecVecXi & ntmp = neighbours.get_neigh();

    // Parameters
    double mu = dati.get_coefficients().at("mu");
    double lam = dati.get_coefficients().at("lam");
    std::size_t ndof = femregion.get_ndof();
    std::size_t fem_degree = femregion.get_fem();
    std::size_t nln = femregion.get_nln();

    // Vector to save degrees of freedom indexes
    VectorSt index(nln);
    VectorSt tmp(nln);

    for (std::size_t ii = 0; ii < nln; ii++)
      tmp[ii] = ii;

    double penalty = 0.0;
    double penalty_scaled = 0.0;

    penalty = (fem_degree == 0) ? dati.get_coefficients().at("penalty") * (lam + 2*mu) : dati.get_coefficients().at("penalty") * fem_degree * fem_degree * (lam + 2*mu);

    for (std::size_t ie = 0; ie < ntmp.size(); ie++) {

      const FeElement & element = FeElements[ie];
      VectorSt partial_index = VectorSt::Constant(nln, (ie)*nln);
      index = partial_index + tmp;

      const auto & neigh_tmp = ntmp[ie];
      std::size_t size_tmp = neigh_tmp.size();

      // Saving normals, scaled edge length, penalty parameter of the edges
      const Points & normals = element.get_normal_edges();
      const std::vector<std::vector<double>> ds_edges = element.get_ds();
      const std::vector<double> & penalties = element.get_penalties();

      //Saving basis function evaluated at quadrature nodes of the edge,
      // gradient, and physical points where quad rule is applied
      const std::vector<MatrixXd> & Phi = element.get_phi_edge();
      const std::vector<VecMatrixXd> & Grad = element.get_Grad_edge();
      const std::vector<Points> & phys_edges = element.get_physical_points_edge();

      for (std::size_t iedg = 0; iedg < size_tmp ; iedg++) {

        //BC applied only where there is the tag of BC itself
        if (neigh_tmp(iedg) == tag) {

          penalty_scaled = penalty * penalties[iedg];
          const std::vector<double> & ds_edge = ds_edges[iedg];
          std::size_t ds_size = ds_edge.size();
          //physical point along this edge
          const Points & phys_edge = phys_edges[iedg];

          for (std::size_t k = 0; k < ds_size; k++){
            auto ds = ds_edge[k];

            //values of basis and gradbasis at this quadrature node
            MatrixXd Bedge(nln,1);
            MatrixXd Gedge_x(nln,1);
            MatrixXd Gedge_y(nln,1);

            for (std::size_t kk = 0; kk < nln; kk++) {
              Bedge(kk) = Phi[iedg](k,kk);
              Gedge_x(kk) = Grad[iedg][kk](k,0);
              Gedge_y(kk) = Grad[iedg][kk](k,1);
            }

            //The point to be used to compute the value of the BC source
            double x = phys_edge[k](0);
            double y = phys_edge[k](1);

            //Muparser evaluation
            std::array<double,3> xx;
            xx[0] = x;
            xx[1] = y;
            xx[2] = time;
            double gd1 = g1(xx);
            double gd2 = g2(xx);

            //just coefficients to be used to impose BC
            double  aa = 0.5 * (lam+2*mu) * normals[iedg](0);
            double  ff = 0.5 * (lam+2*mu) * normals[iedg](1);
            double  bb = 0.5 * lam * normals[iedg](0);
            double  gg = 0.5 * lam * normals[iedg](1);
            double  ee = 0.5 * mu * normals[iedg](0);
            double  cc = 0.5 * mu * normals[iedg](1);

            for (std::size_t i=0; i < nln;i++) { // loop over scalar shape functions
              rhs(index(i)) += penalty_scaled * Bedge(i) * gd1 * ds - 2.0*((aa*Gedge_x(i) + cc*Gedge_y(i))*gd1 + (gg*Gedge_x(i) + ee*Gedge_y(i))*gd2 ) * ds;
              rhs(ndof+index(i)) += penalty_scaled * Bedge(i) * gd2 * ds - 2.0*((bb*Gedge_y(i) + cc*Gedge_x(i))*gd1 + (ff*Gedge_y(i) + ee*Gedge_x(i))*gd2 ) * ds;
            }

          }
        }
      }
    }
  }

}
