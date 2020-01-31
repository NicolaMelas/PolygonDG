/*!
    @file   Triangle.hpp
    @author Nicola Melas
    @brief  Class that defines the class of generic triangle
*/

#ifndef TRIANGLE_HH
#define TRIANGLE_HH

#include "AbstractPolygon.hpp"

namespace PolygonDG {
  /*!
      @brief  Class that defines the class of generic triangle, child
      of father class AbstractPolygon.
      Checks about number of vertexes are needed to addvertex. Specialised versions for
      area(), showMe() and get_max_kb() are built. There are no private members in
      this class. Private/Protected Members are all inside the father class.
      It has one constructor more than other related classes
  */

  class Triangle final: public AbstractPolygon
  {
  public:

    //! Constructor given a vector of points
    /*! Needs to check how many vertexes there are already*/
    Triangle(Points const &);

    //! Constructor given three points
    Triangle(Point2D const &,Point2D const &,Point2D const &);

    //! Copy Constructor
    Triangle(Triangle const &)=default;

    //! Move Constructor
    Triangle(Triangle&&)=default;

    //! assignement operator
    Triangle & operator=(const Triangle &)=default;

    //! Move operator
    Triangle & operator=(Triangle &&)=default;

    //! Specialised for Triangles
    virtual double area() const override;

    //! Adding manually a vertex to the triangle
    inline virtual void addvertex(const Point2D & vertex) override { if (vertexes.size() < 3)  vertexes.push_back(vertex); else (std::cout << "vertexes has not been added") << std::endl;};

    //! Specialised for Triangles
    virtual void showMe(std::ostream & out=std::cout) const override;

    //! 0 if not orderable, 1 if it has been ordered correctly
    bool clockwise();

    //! 0 if not orderable, 1 if it has been ordered correctly
    bool counterclockwise();

    //! 1 if clockwise, 2 if counterclockwise, 0 otherwise
    int orientation() const;

    //! Specialised version for triangles
    /*! The area is the same for all edges because the
        triangle inside the main triangle is the triangle itself*/
    std::vector<double> get_max_kb() const override;

  };
}

#endif
