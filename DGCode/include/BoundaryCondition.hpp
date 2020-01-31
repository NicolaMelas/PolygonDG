/*!
    @file   BoundaryCondition.hpp
    @author Nicola Melas
    @brief  Class that defines the abstract class of Boundary Condition
*/

#ifndef BC_
#define BC_

#include "types.hpp"
#include <iostream>
#include "Coefficients.hpp"
#include "muParserXInterface.hpp"
#include <functional>
#include "Mesh.hpp"
#include "neighbour.hpp"
#include "fespace.hpp"
#include "Points_Utilities.hpp"
#include "Gaussquadrature.hpp"


namespace PolygonDG {

  using MuParserInterface::muParserXInterface;

  /*!
      @brief Class that defines the abstract class of boundaryconditions
      This pure virtual base class defines the properties shared by all BCs.@n
      The key points of the class are the function that define where to impose
      the BC, that allows to understand how to tag the boundary and the function
      apply that, for the Elastodynamics problem, has to act only on the rhs,
      given all the information about the femregion and the naighborhood of each
      element.
  */

  class BoundaryCondition {


  protected:
    //! Tag to identify the boundary
    /*! It must be a negative number */
    int tag;

    //! Function returning true when evaluated on the points where to apply BC
    /*! It's used to modify the adiacency matrix. Note that it has to be used
    to modify only points that are already on the boundary.*/
    std::function<bool(Point2D)> inside_domain;

    //! First component source
    muParserXInterface<3, std::array<double, 3>> & g1;

    //! Second component source
    muParserXInterface<3, std::array<double, 3>> & g2;

  public:
    //! Constructor taking a function, the sources and the tag
    /*!
    @param inside    Function that define the boundary where to apply the condition
    @param     g1_    First component of the BC source
    @param     g2_    Second component of the BC source
    @param     tag_   Tag identifier of the BC boundary (less than zero)
    */
    BoundaryCondition(std::function<bool(Point2D)> inside, muParserXInterface<3, std::array<double, 3>> & g1_,muParserXInterface<3, std::array<double, 3>> & g2_, int tag_)
    : tag(tag_), inside_domain(inside), g1(g1_), g2(g2_) {};

    //! Default Assignment operator
    BoundaryCondition & operator=(BoundaryCondition const&)=default;

    //! Default Copy constructor
    BoundaryCondition(BoundaryCondition const &)=default;

    //! Default Move constructor
    BoundaryCondition(BoundaryCondition&&)=default;

    //! Default Move assignement
    BoundaryCondition & operator=(BoundaryCondition&&)=default;

    //! Returns the tag of the BC (less than zero)
    inline int get_tag() const {return tag;};

    //! True if the point is a point where BC is applied
    inline bool inside(const Point2D & p1) {return inside_domain(p1);};

    //! BC is applied to the rhs vector defined using FeSpace. If the BC source is time dependent, give also the time
    virtual void apply(VectorD & rhs, const FeSpace &, double time = 0.0) const = 0;

    //! Return the BC identifier of the BC type, 1 for Dirichlet, 2 for Neumann
    virtual int get_type() const = 0;

    //! Default destructor
    virtual ~BoundaryCondition() = default;
  };
}

#endif
