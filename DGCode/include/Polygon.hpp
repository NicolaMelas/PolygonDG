/*!
    @file   Polygon.hpp
    @author Nicola Melas
    @brief  Class that defines the class of generic polygon
*/

#ifndef POLYGON_HH
#define POLYGON_HH

#include "AbstractPolygon.hpp"

namespace PolygonDG {
  /*!
      @brief  Class that defines the class of generic polygon, child
      of father class AbstractPolygon.
      No checks are needed to addvertex. Specialised versions for
      area(), showMe() and get_max_kb() are built. There are no private members in
      this class. Private/Protected Members are all inside the father class.
  */

  class Polygon: public AbstractPolygon
  {
  public:

    //! Default constructor.
    Polygon() : AbstractPolygon() {};

    //! Default Assignment operator
    Polygon & operator=(Polygon const&)=default;

    //! Default Copy constructor
    Polygon(Polygon const &)=default;

    //! Default Move constructor
    Polygon(Polygon&&)=default;

    //! Default Move assignement
    Polygon & operator=(Polygon&&)=default;

    //! Polygon may be constructed giving Points;
    Polygon(Points const & v);

    //! Destructor
    virtual ~Polygon(){};

    //! Area of the polygon, Specialised version.
    /*!
      \f$ \int_P d\Omega = 1/2 \int{\partial P} xn_x d\gamma \f$
      The integral is computed by using trapezoidal rule on the polygon sides.
    */
    virtual double area() const override;

    //! Specialised version for generic polygons.
    virtual void showMe(std::ostream & out=std::cout) const override;

    //! Specialised version for generic polygons
    std::vector<double> get_max_kb() const override;

    //! Specialised version for generic polygons
    void addvertex(const Point2D & vertex) override {vertexes.push_back(vertex);};

  };
}
#endif
