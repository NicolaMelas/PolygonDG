/*!
    @file   TriaGenerator.hpp
    @author Nicola Melas
    @brief  Definition of TriaGenerator as child of MeshGenerator
*/

#ifndef TRIAGENERATOR
#define TRIAGENERATOR

#include "MeshGenerator.hpp"

namespace PolygonDG {
  /*!
    @brief  Definition of TriaGenerator as child of MeshGenerator.
    The algorithm is quite simple. Squares are cycled horizontally,
    divided in 2 triangles and the informations are saved in the right
    variables. The number of triangles is ne_x*ne_y*2. For example
    a limit of this generator is that all the dyagonals are with the
    same orientation. It wouldn't be the best choiche in a convection
    problem.
  */
  class TriaGenerator : public MeshGenerator {

  private:

    std::size_t ne_x;
    std::size_t ne_y;

  public:
    //! Constructor
    /*!
       @param domain_  The domain of the Mesh to be generated
       @param ne_x_      Number of subdivision of the x-axis
       @param ne_y_      Number of subdivision of the y-axis
    */
    TriaGenerator(const std::vector<double> & domain_, std::size_t ne_x_, std::size_t ne_y_) : MeshGenerator(domain_, ne_x_ * ne_y_), ne_x(ne_x_), ne_y(ne_y_) {ne *=2;}

    //! Default Assignment operator
    TriaGenerator & operator=(TriaGenerator const&)=default;

    //! Default Copy constructor
    TriaGenerator(TriaGenerator const &)=default;

    //! Default Move constructor
    TriaGenerator(TriaGenerator&&)=default;

    //! Default Move assignement
    TriaGenerator & operator=(TriaGenerator&&)=default;

    /** Specialised version. Modify by reference all the elements needed
        to construct a Mesh. It creates a triangular mesh*/
    void generate(Mesh &) override;

    //! virtual destructor
    virtual ~TriaGenerator()=default;

  };



}


#endif
