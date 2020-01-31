/*!
    @file   NeumannBC.hpp
    @author Nicola Melas
    @brief  Class that defines the child class of Neumann Condition
*/

#ifndef NEUM_BC
#define NEUM_BC

#include "BoundaryCondition.hpp"

namespace PolygonDG{

  using MuParserInterface::muParserXInterface;

  /*!
      @brief Class that defines a child class of BoundaryCondition
      This class simply defines the apply function to add
      \f$ \int{g*v} \f$ on the Neumann Boundary. Moreover it is defined
      the get_type function to distinguish from the Dirichlet case
  */

  class NeumannBC : public BoundaryCondition {


  public:
    //!Constructor
    /*!
    @param inside    Function that define the boundary where to apply the condition
    @param     g1_    First component of the BC source
    @param     g2_    Second component of the BC source
    @param     tag_   Tag identifier of the BC boundary (less than zero)

    It relies on the BoundaryCondition constructor.
    */
    NeumannBC(std::function<bool(Point2D)> inside, muParserXInterface<3, std::array<double, 3>> & g1_,muParserXInterface<3, std::array<double, 3>> & g2_, int tag_)
    : BoundaryCondition(inside, g1_,g2_,tag_) {};

    //! Default Assignment operator
    NeumannBC & operator=(NeumannBC const&)=default;

    //! Default Copy constructor
    NeumannBC(NeumannBC const &)=default;

    //! Default Move constructor
    NeumannBC(NeumannBC&&)=default;

    //! Default Move assignement
    NeumannBC & operator=(NeumannBC&&)=default;

    //! Specialised version for NeumannBC
    void apply(VectorD &rhs,const FeSpace &, double time = 0.0) const override;

    //! Specialised version for NeumannBC, return 2
    int get_type() const override {return 2;};

    //! Destructor
    virtual ~NeumannBC() = default;

  };
}

#endif
