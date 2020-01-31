/*!
    @file   DirichletBC.hpp
    @author Nicola Melas
    @brief  Class that defines the Dirichlet Boundary Condition
*/

#ifndef DIRI_BC
#define DIRI_BC

#include "BoundaryCondition.hpp"

namespace PolygonDG{

  using MuParserInterface::muParserXInterface;

  /*!
      @brief Specialised version of the class BoundaryCondition. In particular
      it has 2 Specialised functions, get_type() that return 1 and apply that
      applies using the penalty coefficient the Dirichlet BC on the weak form.
  */

  class DirichletBC : public BoundaryCondition {

  private:
    const Coefficients & dati;

  public:

    //! Constructor taking a function, the sources and the tag
    /*!
    @param inside    Function that define the boundary where to apply the condition
    @param     g1_    First component of the BC source
    @param     g2_    Second component of the BC source
    @param     tag_   Tag identifier of the BC boundary (less than zero)
    @param     dati_  Coefficients to be used only for Dirichlet BC
    It recalls the father class constructor
    */
    DirichletBC(std::function<bool(Point2D)> inside, muParserXInterface<3, std::array<double, 3>> & g1_,muParserXInterface<3, std::array<double, 3>> & g2_, int tag_, const Coefficients & dati_)
    : BoundaryCondition(inside, g1_,g2_,tag_), dati(dati_) {};

    //! Default Assignment operator
    DirichletBC & operator=(DirichletBC const&)=default;

    //! Default Copy constructor
    DirichletBC(DirichletBC const &)=default;

    //! Default Move constructor
    DirichletBC(DirichletBC&&)=default;

    //! Default Move assignement
    DirichletBC & operator=(DirichletBC&&)=default;

    //! Specialised version for Dirichlet Booundary Condition
    /* Note that to impose DirichletBC, the source function act both on the Test
    function itself and its gradient*/
    void apply(VectorD &rhs,const FeSpace &, double time = 0.0) const override;

    //Specialised version for DirichletBC, return 1
    inline int get_type() const override {return 1;};

    //! Destructor
    virtual ~DirichletBC() = default;

  };
}

#endif
