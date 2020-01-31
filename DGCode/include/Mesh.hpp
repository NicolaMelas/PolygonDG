/*!
    @file   Mesh.hpp
    @author Nicola Melas
    @brief  Definition of Mesh
*/

#ifndef REGION
#define REGION
#include "Coefficients.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/utility.h>
#include "Polygon.hpp"
#include "types.hpp"


namespace PolygonDG {
  /*!
      @brief  Definition of Mesh.
      It represents the region where to study the problem. It is
      composed of polygons, through pointers to them. Then, information about
      the connectivity, the max diameter, and the points are directly stored
      in a Mesh Object. The Mesh class has as friend classes, the generators of
      Mesh.
  */
  class Mesh {


  protected:
    //! domain representing [xmin,xmax,ymin,ymax]
    std::vector<double> domain;
    //!number of elements
    std::size_t  ne;
    //! all the points
    Points  coord;
    //!Pointers to polygons that represents the mesh
    std::vector<std::shared_ptr<AbstractPolygon>>  polygons;
    //! Connectivity matrix
    VecVecXi  connectivity;
    //! maximum diameter
    double diameter;

    //! Function that compute the diameter of the mesh (it's protected)
    void compute_diameter();

  public:
    //! Default constructor
    Mesh() = default;
    /*! Constructor of the mesh. It automatically calls the generator of mesh
      @param p1     Bottom-Left Point of the Rectangle Mesh
      @param p2     Top-Right Point of the Rectangle Mesh
      @param ne_x   Number of subdivision of the x-axis
      @param ne_y   Number of subdivision of the y-axis
      @param mesh_type Polygon or Triangle
      */
    Mesh(const Point2D & p1, const Point2D & p2, std::size_t ne_x, std::size_t ne_y, std::string mesh_type = "Polygon");

    //! Default Assignment operator
    Mesh & operator=(Mesh const&)=default;

    //! Default Copy constructor
    Mesh(Mesh const &)=default;

    //! Default Move constructor
    Mesh(Mesh&&)=default;

    //! Default Move assignement
    Mesh & operator=(Mesh&&)=default;

    //! destructor
    ~Mesh()= default;

    //! Number of edges of the element ie
    /*! @param ie   ie-th element of the Mesh*/
    inline std::size_t get_nedge(std::size_t ie) const {return polygons[ie]->size();}

    //! Return the number of the Mesh elements
    inline std::size_t get_ne() const {return ne;}

    //! Vertices of the Mesh (as const reference!)
    inline const Points & get_coord() const {return coord;}

    //! Returns the BBox of the given element, containing, for each element, x_min, x_max, y_min, y_max
      /*! @param ie ie-th element of the Mesh */
    inline std::vector<double> get_BBox(std::size_t ie) const {return polygons[ie]->get_BBox();}

    //! Returns the pointers to the elements of the mesh
    inline const std::vector<std::shared_ptr<AbstractPolygon>> & get_coords_element() const {return polygons;}

    //! Get the connectivity, the indices of the points contained in each element
    inline const VecVecXi & get_connectivity() const {return connectivity;}

    //! Return area of an element
    /*! @param ie ie-th element of the Mesh */
    inline double get_area(std::size_t ie) const {return polygons[ie]->area();}

    /*! For each edge of each element it contains the max area of the triangle composed
    by the edge and one of the other vertexes of the element
    @param ie ie-th element of the Mesh */
    inline std::vector<double> get_max_kb(std::size_t ie) const {return polygons[ie]->get_max_kb();}

    //! Return the mesh size (maximum diameter)
    inline double get_diameter() const {return diameter;}

    //! Print info about the mesh in an output stream
    void print_mesh(std::ostream & out=std::cout, bool print_polygons = 0) const;

    //! Export as paraview file
    void printvtk(std::string fileName) const;

    friend class MeshGenerator;
    friend class PolyGenerator;
    friend class TriaGenerator;

  };

}


#endif
