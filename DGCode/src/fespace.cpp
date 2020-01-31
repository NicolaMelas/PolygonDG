/*!
    @file   fespace.cpp
    @author Nicola Melas
    @brief  Implementation of class FeSpace
*/

#include "fespace.hpp"

namespace PolygonDG {

  FeSpace::FeSpace(const Mesh & region_,Neighbour & neighbour_, std::size_t fem_degree_)  :  region(region_), neighbour(neighbour_), fem_degree(fem_degree_) {

    nqn = 2 * fem_degree_ + 1;
    nln = 0.5 * (fem_degree_ + 1) * (fem_degree_ + 2);
    ndof = nln*region.get_ne();
    compute_FElements();
  }

  void FeSpace::compute_FElements() {

    const std::vector<std::shared_ptr<AbstractPolygon>> & polygons = region.get_coords_element();
    FeElements.reserve(region.get_ne());
    for (std::size_t ii = 0; ii < region.get_ne(); ii++ )
      FeElements.push_back(FeElement(nqn,fem_degree,ii,polygons,neighbour));

  }

}
