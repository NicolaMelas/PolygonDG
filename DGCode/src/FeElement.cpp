/*!
    @file   FeElement.cpp
    @author Nicola Melas
    @brief  Implementation of class FeElements
*/

#include "FeElement.hpp"

namespace PolygonDG {

  FeElement::FeElement(std::size_t nqn, std::size_t fem_degree_, std::size_t id_, const std::vector<std::shared_ptr<AbstractPolygon>> & polygons_, const Neighbour & neighbour ) : id(id_) {

    fem_degree = fem_degree_;
    auto polygon = polygons_[id_];
    //defining the quadrature
    Gauss_Quadrature quadrature(nqn);

    auto nodes_weights2D = quadrature.GLRefTria2D();
    const std::vector<double> & w_2D = nodes_weights2D.second;
    const auto & qnodes_2D = nodes_weights2D.first;

    auto nodes_weights1D = quadrature.GauLeg();
    const std::vector<double> & w_1D = nodes_weights1D.second;
    const auto & qnodes_1D = nodes_weights1D.first;

    // Set all the things depending directly on the polygon
    normal_edges = polygon->get_normals();
    triangulation = Triangulation(polygon);
    number_tria = triangulation.size();

    VectorI neigh_ie = neighbour.get_neigh()[id];
    VectorI neighedges_ie = neighbour.get_neighedges()[id];

    Points coords_element = polygon->theVertexes();
    std::vector<double> meshsize = polygon->get_edges_length();


    for (std::size_t iTria = 0; iTria < number_tria; iTria++) {

      const std::shared_ptr<Triangle> v1v2v3 = triangulation[iTria];

      MatrixXd BJ = get_jacobian_and_translation(v1v2v3).first;
      Points pphys2D = get_physical_points(v1v2v3, qnodes_2D);
      double Jdet = BJ.determinant();


      auto dphiq_Grad = evalshape2D(polygon, pphys2D);
      //saving physical points of quadrature rule
      physical_points_tria.push_back(pphys2D);
      //saving basis function at the considered triangle
      phi_trias.push_back(dphiq_Grad.first);
      //saving gradbasis function at the considered triangle
      Grad_trias.push_back(dphiq_Grad.second);

      std::vector<double> tmp_dx = w_2D;
      std::transform(tmp_dx.begin(), tmp_dx.end(), tmp_dx.begin(),[&Jdet](const double & a){return a * std::abs(Jdet);});
      //saving coefficients for quadrature
      dx_trias.push_back(tmp_dx);

    }

    // save info about edges
    std::size_t nedge = coords_element.size();
    auto neight = neighbour.get_neigh()[id];
    auto neighedge = neighbour.get_neighedges()[id];
    double area = polygon->area();

    double penalty_edge;

    for (std::size_t iedg = 0; iedg < nedge; iedg++) {

      auto neigedge = neighedge[iedg];
      auto max_kb = polygon->get_max_kb();
      VectorD eigen_kb = Eigen::Map<VectorD, Eigen::Unaligned>(max_kb.data(), max_kb.size());
      VectorD Cinv = VectorD::Constant(nedge,area).cwiseQuotient(eigen_kb);

      if (neight(iedg) < 0) {
        penalty_edge = Cinv(iedg) * meshsize[iedg] / area;
      }
      else if (neight(iedg) > -1){
        auto polyneigh = polygons_[neight(iedg)];
        double areaneigh = polyneigh-> area();
        auto Cinv_ext = areaneigh / polyneigh->get_max_kb()[neigedge];
        auto s1 = Cinv(iedg) * meshsize[iedg] / area;
        auto s2 = Cinv_ext * meshsize[iedg] / areaneigh;
        penalty_edge = std::max(s1,s2);
      }
      //saving penalties
      penalty_edges.push_back(penalty_edge);
      std::vector<double> tmp_ds = w_1D;
      std::transform(tmp_ds.begin(), tmp_ds.end(), tmp_ds.begin(),[&meshsize,&iedg](const double & a){return a * meshsize[iedg];});
      //saving quadrature coefficients
      ds_edges_poly.push_back(tmp_ds);

      std::vector<Edge> edges = polygon->get_edges();
      Points pphys_1D = get_physical_points_faces(edges[iedg], qnodes_1D);
      //Saving pysical edge point for 1d quadrature rule
      physical_points_edge.push_back(pphys_1D);

      auto B_ed_G_ed = evalshape2D(polygon, pphys_1D);

      //Saving basis and gradbasis on the edge
      phi_edges.push_back(B_ed_G_ed.first);
      Grad_edges.push_back(B_ed_G_ed.second);
    }
  }



  std::pair<MatrixXd,VecMatrixXd> FeElement::evalshape2D(const std::shared_ptr<AbstractPolygon> polygon, const Points & Nodes) const {

    std::size_t N = fem_degree;
    auto psi = BasisNodes(polygon, Nodes, N);
    auto Grad = GradBasisNodes(polygon, Nodes, N);
    return std::make_pair(psi, Grad);
  }

}
