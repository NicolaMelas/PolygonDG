/*!
    @file   neighbour.hpp
    @author Nicola Melas
    @brief  Definition of the adiacency information
*/
#ifndef NEIGHBOUR
#define NEIGHBOUR
#include <Eigen/Dense>
#include "Coefficients.hpp"
#include "Mesh.hpp"
#include "types.hpp"


namespace PolygonDG {

  /*!
      @brief  Definition of the adiacency information. The boundaries
      are set by default to -1. The important thing is that is a
      negative number. The class scope is to stored the information of the
      adiacency elements, used in the computation of the averages and the jumps
      along edges in computing the SIP method to solve the
      ElastodynamicsProblem.
  */

  class Neighbour {

  private:
    /*! It stores the indexes of the neigh elements to each edge of the
    considered element*/
    VecVecXi  neigh;
    /*! It stores the indexes of the edge of the neigh elements to each edge of the
    considered element*/
    VecVecXi  neighedges;

    void create_neigh(const Mesh &);

  public:

    //! Constructor
    /*@param region The Mesh where to build the Adiacency Matrix*/
    Neighbour(const Mesh & region) {create_neigh(region);}

    //! Default Assignment operator
    Neighbour & operator=(Neighbour const&)=default;

    //! Default Copy constructor
    Neighbour(Neighbour const &)=default;

    //! Default Move constructor
    Neighbour(Neighbour&&)=default;

    //! Default Move assignement
    Neighbour & operator=(Neighbour&&)=default;

    //! Return for each element the numbers of the near elements
    inline const VecVecXi & get_neigh() const {return neigh;}

    //! Return for each element the edge number of the near elements
    inline const VecVecXi & get_neighedges() const {return neighedges;}

    /*! Set, given the number of the element and the number of the edge,
     the boundarytag to the given tag. Only for boundary edges.
     It's necessary to handle correclty the boundary conditions*/
    inline void set_boundary_tag(std::size_t elem_index, std::size_t edge_index, int tag) {
      if (tag < 0)  {
        neigh[elem_index](edge_index) = tag;
        neighedges[elem_index](edge_index) = tag;
      }
    };


  };
}

#endif
