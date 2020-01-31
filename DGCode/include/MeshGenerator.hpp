/*!
    @file   MeshGenerator.hpp
    @author Nicola Melas
    @brief  Definition of MeshGenerator
*/

#ifndef MESHGENERATOR
#define MESHGENERATOR


#include "Mesh.hpp"
#include "Polygon.hpp"
#include "Triangle.hpp"
#include <algorithm>
#include <Eigen/KroneckerProduct>
#include <Eigen/Sparse>
#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/foreach.hpp>
#include <cmath>
#include "types.hpp"

namespace PolygonDG {
  /*!
    @brief  Definition of MeshGenerator. It's an abstract class
    due to the presence of the function generate, declared
    pure virtual. This function acts directly on a Mesh to build it.
    To do it, the generator needs to know the domain and the number
    of elements. How the generation works is done by children classes
  */

  class MeshGenerator {

  protected:
    //! domain to be generated
    std::vector<double> domain;

    //! number of element
    std::size_t ne;

  public:
    //! Constructor
    /*!@param domain_  The domain of the Mesh to be generated
       @param ne_      Number of elements
    */
    MeshGenerator(const std::vector<double> & domain_, std::size_t ne_) : domain(domain_), ne(ne_){};

    //! Default Assignment operator
    MeshGenerator & operator=(MeshGenerator const&)=default;

    //! Default Copy constructor
    MeshGenerator(MeshGenerator const &)=default;

    //! Default Move constructor
    MeshGenerator(MeshGenerator&&)=default;

    //! Default Move assignement
    MeshGenerator & operator=(MeshGenerator&&)=default;

    //! Pure virtual generate function. It must be defined in children classes
    virtual void generate(Mesh &) = 0;

    //! Virtual destructor
    virtual ~MeshGenerator()=default;

  };

}

#endif
