/*!
    @file   neighbour.cpp
    @author Nicola Melas
    @brief  Implementation of class Neighbour
*/

#include "neighbour.hpp"
#include <iostream>

namespace PolygonDG {

  void Neighbour::create_neigh(const Mesh & region) {

    auto ne = region.get_ne();
    const auto & connectivity = region.get_connectivity();

    neigh.resize(ne);
    neighedges.resize(ne);

    //! Default inizialization to tag = -1
    for(std::size_t i = 0; i < ne; i++) {
      neigh[i] = -1 * VectorI::Ones(connectivity[i].size());
      neighedges[i] = -1 * VectorI::Ones((connectivity[i]).size());
    }

    for (std::size_t i = 0; i < ne -1; i++) {

      auto n_edges = connectivity[i].size();
      MatrixXi edges(n_edges,2);
      VectorI  v(connectivity[i].size());
      v = connectivity[i];

      for (int j = 0; j < n_edges -1; j++) {
        MatrixXi tmp(1,2);
        tmp << v(j), v(j+1);
        edges.block(j,0,1,2) = tmp;
      }
      edges(n_edges-1,0) = v(n_edges-1);
      edges(n_edges-1,1) = v(0);

      //checking if in the other elements there is an element that shares an edge
      for (std::size_t k = i + 1; k < ne; k++) {

        auto n_edgesn = connectivity[k].size();
        MatrixXi edgesn(n_edgesn,2);
        VectorI  vn(n_edgesn);
        vn = connectivity[k];

        for (int l = 0; l < n_edgesn -1; l++) {
          MatrixXi tmp(1,2);
          tmp << vn(l), vn(l+1);
          edgesn.block(l,0,1,2) = tmp;
        }

        edgesn(n_edgesn-1,0) = vn(n_edgesn-1);
        edgesn(n_edgesn-1,1) = vn(0);

        for (int s = 0; s < n_edges; s++) {
          for (int t = 0; t < n_edgesn; t++) {
            if (edges(s,0) == edgesn(t,1) && edges(s,1) == edgesn(t,0) ) {
              neigh[i](s)=k;
              neigh[k](t)=i;
              neighedges[i](s)=t;
              neighedges[k](t)=s;

            }
          }
        }
      }
    }
  }
}
