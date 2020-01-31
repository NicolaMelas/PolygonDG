/*!
    @file   FeElement.hpp
    @author Nicola Melas
    @brief  The class containt the information of an element of the FeSpace
*/

#ifndef FE_ELEMENT
#define FE_ELEMENT
#include "types.hpp"
#include "neighbour.hpp"
#include "Gaussquadrature.hpp"
#include "Transformation.hpp"
#include "Triangulation.hpp"
#include "Legendre.hpp"

namespace PolygonDG {
  /*!
      @brief The class containt the information of an element of the FeSpace.
      It stores everything needed to assemble the system. For each Polygon
      of the mesh, a triangulation has been done. Then the reference triangle
      quadrature rules are mapped on the triangles of the triangulation. On those
      triangles basis functions are computed and then, all information about basis
      function and quadrature rules are available. the same info are saved for edge
      with 1D quadrature.
  */

  class FeElement {

  private:

    //! Id of the polygon
    std::size_t id;

    //! Fem degree
    std::size_t fem_degree;

    //! Number of triangles inside the Polygon
    /*!  The method requires a triangulation process, so
    the number of triangles is stored */
    std::size_t number_tria;

    //! Info about the triangulation
    /* Each element of the Vector contains info about the vertexes*/
    std::vector<std::shared_ptr<Triangle>> triangulation;

    //! For each quadrature node of the triangle, the integration coefficients are stored
    std::vector<std::vector<double>> dx_trias;

    //! Value of the basis functions on each triangle
    /*! The values are stored for each dof(column) and for each
    quadrature nodes(rows)*/
    std::vector<MatrixXd> phi_trias;

    //! values of the Gradient
    /*For each Gradient..for each degree of freedom of the triangle, the matrix contains for nqn rows,
     grad_x in the first column and grad_y in the second one*/
    std::vector<VecMatrixXd> Grad_trias;

    //! For each edge, the penalty coefficient used for the S Matrix
    std::vector<double> penalty_edges;

    //! For each edge and for each quadrature nodes, integration coefficients are stored
    std::vector<std::vector<double>> ds_edges_poly;

    //! value of the basis functions on edges
    /*! The values are stored for each dof(column) and for each
    quadrature nodes(rows)*/
    std::vector<MatrixXd> phi_edges;

    //! value of the Gradient of the basis functions on edges
    /*For each Gradient..for each degree of freedom of the edge, the matrix contains for nqn rows,
     grad_x in the first column and grad_y in the second one*/
    std::vector<VecMatrixXd> Grad_edges;

    //! the normals to the edges
    Points normal_edges;

    //! For each triangle, return the physical points of the integration rule
    /*! These are the points used to perform 2D quadrature rule*/
    std::vector<Points> physical_points_tria;

    //! For each edge, return the physical points of the integration rule
    /*! These are the points used to perform 1D quadrature rule*/
    std::vector<Points> physical_points_edge;

    /*! Evaluate the basis function and its gradient on the Physical Nodes.
    @param polygon    The polygon on which the basis has to be computed
    @Nodes            The points where the basis and gradients are computed @n
    Note that physical points are mapped back to compute values.*/
    std::pair<MatrixXd,VecMatrixXd> evalshape2D(const std::shared_ptr<AbstractPolygon> polygon, const Points & Nodes) const;

  public:
    //! Constructor
    /*!
    @param nqn      Number of quadrature nodes
    @param fem_degree_ Polynomial degree
    @param id_      Polygon Id
    @param polygons Vector with all the pointers to the Mesh elements
    @param neighbour  Adiacency matrix
    */
    FeElement(std::size_t nqn, std::size_t fem_degree_, std::size_t id_, const std::vector<std::shared_ptr<AbstractPolygon>> & polygons, const Neighbour & neighbour);

    //! Default Assignment operator
    FeElement & operator=(FeElement const&)=default;

    //! Default Copy constructor
    FeElement(FeElement const &)=default;

    //! Default Move constructor
    FeElement(FeElement&&)=default;

    //! Default Move assignement
    FeElement & operator=(FeElement&&)=default;

    //!Return the number of triangles which are used to subdivide the polygon
    inline std::size_t get_number_tria() const {return number_tria;}

    /*! Return the physical integration point used to compute numerical integration
     on the triangles*/
    inline const std::vector<Points> & get_physical_points_tria() const {return physical_points_tria;}

    /*! Return the physical integration point used to compute numerical integration
     on the edges*/
    inline const std::vector<Points> & get_physical_points_edge() const {return physical_points_edge;}

    /*! Return for each triangle of the triangulation the scaled coefficients
    used in the quadrature rules*/
    inline const std::vector<std::vector<double>> & get_dx() const {return dx_trias;}

    //! Return the basis function computed inside the triangles (of the triangulation)
    inline const std::vector<MatrixXd> & get_phi_tria() const {return phi_trias;};

    //! Return the grad of basis function
    /*For each Gradient..for each degree of freedom of the triangle, the matrix contains for nqn rows,
     grad_x in the first column and grad_y in the second one*/
    inline const std::vector<VecMatrixXd> & get_Grad_tria() const {return Grad_trias;}

    //!Return penalty for each edge. It is used to compute correctly the Matrix S
    inline const std::vector<double> & get_penalties() const {return penalty_edges;}

    //! Return size of each scaled edge. Used in quadrature integration.
    inline const std::vector<std::vector<double>> & get_ds() const {return ds_edges_poly;}

    //! Return basis function on the edges boundary
    inline const std::vector<MatrixXd> & get_phi_edge()  const {return phi_edges;}

    //! Return value of the Gradient of the basis functions on edges
    /*For each Gradient..for each degree of freedom of the edge, the matrix contains for nqn rows,
     grad_x in the first column and grad_y in the second one*/
    inline const std::vector<VecMatrixXd> & get_Grad_edge() const {return Grad_edges;}

    //! Return the normals to the edges
    inline const Points & get_normal_edges() const {return normal_edges;}

    //! Return the triangles that compose the polygon
    inline const std::vector<std::shared_ptr<Triangle>> & get_triangulation() const {return triangulation;};

  };
}

#endif
