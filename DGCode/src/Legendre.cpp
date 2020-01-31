/*!
    @file   Legendre.cpp
    @author Nicola Melas
    @brief  Implementation of LegendreFunctions
*/

#include "Legendre.hpp"
#include "Transformation.hpp"

namespace PolygonDG {

  VectorD LegendreP(const VectorD & x, std::size_t N)  {

    std::size_t dims = x.size();

    MatrixXd PL = MatrixXd::Zero(N+1,dims);

    // Initial values P_0(x) and P_1(x)
    PL.block(0,0,1,dims) = MatrixXd::Ones(1,dims);

    if (N==0)
      return PL.transpose();

    PL.block(1,0,1,dims) = x.transpose();
    if (N==1)
      return PL.block(N,0,1,dims).transpose();
    //computing (in block) the coeff of degree i
    for(std::size_t i = 1; i < N ; i++)
      PL.block(i+1,0,1,dims) = 1.0/(i+1.0)*((2.0*i+1.0)*(x.transpose()).cwiseProduct(PL.block(i,0,1,dims)) - double(i)*PL.block(i-1,0,1,dims));


    return PL.block(N,0,1,dims).transpose();
  }

  VectorD GradLegendreP(const VectorD & x, std::size_t N)  {

    std::size_t dims = x.size();

    MatrixXd dPL = MatrixXd::Zero(N+1,x.size());

    // Initial values P_0(x) and P_1(x)
    dPL.block(0,0,1,dims) = MatrixXd::Zero(1,dims);

    if (N==0)
      return dPL.transpose();

    dPL.block(1,0,1,dims) = MatrixXd::Ones(1,dims);

    if (N==1)
      return dPL.block(N,0,1,dims).transpose();

    MatrixXd dP = MatrixXd::Zero(x.size(),1);

    for (int i=N-1;i >=0; i-=2) {
      auto P = LegendreP(x,i);
      dP += 2.0*P/(2.0/(2.0*i+1.0));
    }

    return dP;
  }

  VectorD Basis2DP(const VectorD & a, const VectorD & b, std::size_t i, std::size_t j) {

    VectorD h1 = LegendreP(a,i);
    VectorD h2 = LegendreP(b,j);
    double c = std::sqrt((2.0*i+1.0)*(2.0*j+1.0)/4.0);
    //2D basis
    VectorD P = c*h1.cwiseProduct(h2);

    return P;

  }


  std::pair<VectorD,VectorD> GradBasis2DP(const VectorD & x, const VectorD & y, std::size_t id, std::size_t jd) {

    VectorD Lx = LegendreP(x,id);
    VectorD dLx = GradLegendreP(x, id);
    VectorD Ly = LegendreP(y,jd);
    VectorD dLy = GradLegendreP(y, jd);

    // x-derivative
    double c = std::sqrt((2.0*id+1.0)*(2.0*jd+1.0)/4.0);
    VectorD dmodedx = (c*dLx).cwiseProduct(Ly);

    // y-derivative
    VectorD dmodedy = (c*Lx).cwiseProduct(dLy);
    //2D gradbasis
    return std::make_pair(dmodedx, dmodedy);
  }

  MatrixXd BasisNodes(const std::shared_ptr<AbstractPolygon> polygon,const Points & Nodes, std::size_t fem_degree) {

    std::size_t N = fem_degree;
    std::size_t n_polys = (N+1)*(N+2)*0.5;
    std::size_t nodsiz = Nodes.size();

    Points Nodes_ref(nodsiz);

    std::pair<MatrixXd,VectorD> BJ_trans = get_jacobian_and_translation(polygon);
    const MatrixXd & BJ = BJ_trans.first;  // Jacobian of elemental map
    const VectorD & translation = BJ_trans.second;                // translation vector

    double D = BJ.determinant();
    MatrixXd BJ_inv = BJ.inverse();

    VectorD trans_inv = VectorD::Zero(2);
    trans_inv << -BJ(1,1)*translation(0)+BJ(0,1)*translation(1), BJ(1,0)*translation(0)-BJ(0,0)*translation(1);
    trans_inv /= D;

    std::vector<double> a_(nodsiz,0.0);
    std::vector<double> b_(nodsiz,0.0);

    //coming back from physical points
    for (std::size_t k=0; k < nodsiz;k++) {
      const VectorD & tmp = Nodes[k];
      Nodes_ref[k] = (BJ_inv * tmp + trans_inv);
    }

    a_ = get_x(Nodes_ref);
    b_ = get_y(Nodes_ref);

    VectorD a = Eigen::Map<VectorD, Eigen::Unaligned>(a_.data(), a_.size());
    VectorD b = Eigen::Map<VectorD, Eigen::Unaligned>(b_.data(), b_.size());

    MatrixXd psi = MatrixXd::Zero(nodsiz,n_polys);

    std::size_t sk = 0;
    for (std::size_t i = 0; i < N+1; i++) {
      for (std::size_t j = 0; j < N + 1 - i ; j++) {
        psi.block(0,sk,psi.rows(),1) = Basis2DP(a,b,i,j);
        sk = sk+1;
      }
    }
    return psi;
  }

  VecMatrixXd GradBasisNodes(const std::shared_ptr<AbstractPolygon> polygon,const Points & Nodes, std::size_t fem_degree) {

    std::size_t N = fem_degree;
    std::size_t n_polys = (N+1)*(N+2)*0.5;
    std::size_t nodsiz = Nodes.size();

    Points Nodes_ref(nodsiz);

    std::pair<MatrixXd,VectorD> BJ_trans = get_jacobian_and_translation(polygon);    // Jacobian of elemental map
    const MatrixXd & BJ = BJ_trans.first;
    const VectorD & translation = BJ_trans.second;                // translation vector

    double D = BJ.determinant();
    MatrixXd BJ_inv = BJ.inverse();

    VectorD trans_inv = VectorD::Zero(2);
    trans_inv << -BJ(1,1)*translation(0)+BJ(0,1)*translation(1), BJ(1,0)*translation(0)-BJ(0,0)*translation(1);
    trans_inv /= D;

    std::vector<double> a_(nodsiz,0.0);
    std::vector<double> b_(nodsiz,0.0);

    for (std::size_t k=0; k < nodsiz;k++) {
      const VectorD & tmp = Nodes[k];
      Nodes_ref[k] = (BJ_inv * tmp + trans_inv);
    }

    a_ = get_x(Nodes_ref);
    b_ = get_y(Nodes_ref);

    VectorD a = Eigen::Map<VectorD, Eigen::Unaligned>(a_.data(), a_.size());
    VectorD b = Eigen::Map<VectorD, Eigen::Unaligned>(b_.data(), b_.size());

    VecMatrixXd dpsi(2,MatrixXd::Zero(nodsiz,n_polys));

    std::size_t sk = 0;
    for (std::size_t i = 0; i < N+1; i++) {
      for (std::size_t j = 0; j < N + 1 - i ; j++) {
        auto dsp1_x_y = GradBasis2DP(a,b,i,j);
        dpsi[0].block(0,sk,dpsi[0].rows(),1) = dsp1_x_y.first;
        dpsi[1].block(0,sk,dpsi[1].rows(),1) = dsp1_x_y.second;
        sk = sk+1;
      }
    }

    VecMatrixXd Grad(n_polys);

    for (std::size_t i = 0; i < n_polys; i++) {
      Grad[i] = MatrixXd::Zero(nodsiz,2);
    }

    for (std::size_t i = 0; i < n_polys; i++) {
      for (int j = 0; j < Grad[i].rows(); j++) {
        VectorD tmp = VectorD::Zero(2);
        tmp << dpsi[0](j,i), dpsi[1](j,i);
        (Grad[i]).block(j,0,1,2) = tmp.transpose() * BJ_inv;
      }
    }
    return Grad;
  }
}
