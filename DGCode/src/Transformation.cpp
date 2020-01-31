/*!
    @file   Transformation.cpp
    @author Nicola Melas
    @brief  Implementation of transformations function from ref. element
*/

#include "Transformation.hpp"

namespace PolygonDG {
  std::pair<MatrixXd,VectorD> get_jacobian_and_translation(const std::shared_ptr<AbstractPolygon> polygon)  {

    MatrixXd BJ = MatrixXd::Zero(2,2);
    VectorD translation = VectorD::Zero(2);

    Points coords = polygon->theVertexes();
    //Building Jacobian and translation vector
    if (polygon->size() == 3) {

      double x0=coords[0](0);
      double x1=coords[1](0);
      double x2=coords[2](0);

      double y0=coords[0](1);
      double y1=coords[1](1);
      double y2=coords[2](1);

      BJ << x1-x0, x2-x0, y1-y0, y2-y0;
      translation << x0,y0;

    }
    else {

      double x1B = polygon->x_min(); double x2B = polygon->x_max();
      double y1B = polygon->y_min(); double y2B = polygon->y_max();

      double x0=x1B;   // x-coordinates of vertices
      double x1=x2B;
      double x2=x2B;
      double x3=x1B;

      double y0=y1B;   // y-coordinates of vertices
      double y1=y1B;
      double y2=y2B;
      double y3=y2B;

      BJ << 0.25 * (-x0 + x1 + x2 - x3) , 0.25 * ( + x3 -x0 - x1 + x2   ), 0.25 * (-y0 + y1 + y2 - y3) , 0.25 * (-y0 - y1 + y2 + y3);
      translation << (0.25) * (x0 + x1 + x2 + x3), 0.25*( y0 + y1 + y2 + y3);    // translation vector
    }
    return std::make_pair(BJ,translation);
  }

  Points get_physical_points(const std::shared_ptr<AbstractPolygon> polygon, const Points & nodes_2D)  {

    Points coords = polygon->theVertexes();
    //From points beginning to reference element, trasnform to the corresponding points inside the polygon
    std::pair<MatrixXd,VectorD> BJ_trans = get_jacobian_and_translation(polygon);

    const MatrixXd & BJ = BJ_trans.first;
    const MatrixXd & translation = BJ_trans.second;

    Points pphys_2D(nodes_2D.size(),Point2D::Zero());

    for(std::size_t k = 0; k < nodes_2D.size(); k++) {
      VectorD tmp1 = VectorD::Zero(2);
      tmp1 = nodes_2D[k];
      VectorD tmp2 =  BJ * tmp1;
      tmp2 += translation;
      pphys_2D[k] = tmp2.transpose();
    }
    return pphys_2D;
  }

  Points get_physical_points_faces(const Edge & edge, const std::vector<double> & node1D)  {

    Point2D p1 = edge.first();
    Point2D p2 = edge.second();
    Point2D p3 = 0.5*(p1+p2);

    //From points beginning to 1D reference element, trasnform to the corresponding points of the edge
    auto tria = std::make_shared<Triangle>(p1,p2,p3);

    std::size_t nqn_1D=node1D.size();
    MatrixXd nodes_face = MatrixXd::Zero(nqn_1D,2); // zeros(nqn_1D,2,nfaces);
    const VectorD node_1D = Eigen::Map<const VectorD, Eigen::Unaligned>(node1D.data(), node1D.size());

    nodes_face.block(0,0,nqn_1D,1) = node_1D;

    auto BJ_trans = get_jacobian_and_translation(tria);
    const MatrixXd & BJ_face = BJ_trans.first;
    const VectorD & translation = BJ_trans.second;

    Points pphys_1D(nqn_1D,Point2D::Zero());

    for(std::size_t k = 0; k < nqn_1D; k++){
      VectorD tmp1 = nodes_face.block(k,0,1,2).transpose();
      VectorD tmp2 =  BJ_face * tmp1;
      tmp2 += translation;
      pphys_1D[k]= tmp2;
    }

    return pphys_1D;
  }

}
