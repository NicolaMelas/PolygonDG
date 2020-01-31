/*!
    @file   fespace.hpp
    @author Nicola Melas
    @brief  The class containt the information about FeSpace
*/

#ifndef FEMREGION
#define FEMREGION
#include "Coefficients.hpp"
#include "Mesh.hpp"
#include "Edge.hpp"
#include "FeElement.hpp"

namespace PolygonDG {
  /*!
      @brief  The class containt the information about FeSpace. It contains
      the vector of the FeElements. It can be defined passing the Mesh and
      the fem_degree. The FeElements are the key point of the fespace. The
      Neighbour is needed to describe the connectivity.
  */
  class FeSpace {

  private:

    const Mesh & region;
    Neighbour & neighbour;
    std::vector<FeElement> FeElements;
    std::size_t fem_degree;
    std::size_t nln; // local degree of freedom
    std::size_t ndof; // total  degrees of freedom
    std::size_t nqn; // number of quadrature nodes

    //! Computing the vector of FeElements with all the important information
    void compute_FElements();

  public:
    //! Construct the Finite element space of degree fem_degree on the given Mesh
    /*!
    @param region   The Mesh on which the FeSpace will be built
    @param neighbour The adiacency matrix to describe the neighborhood of each element
    @param fem_degree The polynomial degree of basis functions
    */
    FeSpace(const Mesh & region, Neighbour& neighbour, std::size_t fem_degree);

    //! Default Assignment operator
    FeSpace & operator=(FeSpace const&)=default;

    //! Default Copy constructor
    FeSpace(FeSpace const &)=default;

    //! Default Move constructor
    FeSpace(FeSpace&&)=default;

    //! Default Move assignement
    FeSpace & operator=(FeSpace&&)=default;

    //! Returns the fem degree
    inline std::size_t get_fem() const {return fem_degree;}

    //! Number of dof on each element
    inline std::size_t get_nln() const {return nln;}

    //! The dofs of the whole mesh
    inline std::size_t get_ndof() const {return ndof;}

    //! Returns the number of quadrature nodes used to compute ingration
    inline std::size_t get_nqn() const {return nqn;}

    //! Return the number of element of the Mesh
    inline std::size_t get_ne() const {return region.get_ne();}

    //! Return the Mesh (as const Reference)
    inline const Mesh & get_mesh() const {return region;}

    //! Return the Adiacency (as Reference!!)
    inline Neighbour & get_Neighbour() const {return neighbour;}

    //! Return the vector of AbstractPolygon pointers of elements of the Mesh (as const Reference)
    inline const std::vector<std::shared_ptr<AbstractPolygon>> & get_coords_element() const {return region.get_coords_element();}

    //! Return the vector of FeElements as const Reference
    inline const std::vector<FeElement> & get_FeElements() const {return FeElements;};

    //! Return the number of edges of the element ie of the Mesh
    /* @param ie    ie-th element of the Mesh */
    inline std::size_t get_nedge(std::size_t ie) const {return region.get_nedge(ie);};

  };
}

#endif
