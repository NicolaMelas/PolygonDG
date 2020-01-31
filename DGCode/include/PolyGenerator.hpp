/*!
    @file   PolyGenerator.hpp
    @author Nicola Melas
    @brief  Definition of PolyGenerator as child of MeshGenerator
*/

#ifndef POLYGENERATOR
#define POLYGENERATOR

#include "MeshGenerator.hpp"

namespace PolygonDG {

  //! Componentwise-distance
  RowMatrixXd dRectangle(const Points & P, const std::vector<double> & BdBox);

  /*!
    @brief  Definition of PolyGenerator as child of MeshGenerator.
    The algorithm is based on the voronoi diagram performed on a set of points,
    in particular on the given points(centroid of polygons) and their reflection
    with respect to an edge of the boundary. This assure that the mesh is built
    correctly. The Mesh is homogeneized using the LLyod algorithm.The Generator
    uses the Qhull library to perform the voronoi diagram and the boost library
    to perform the perform Cuthill-mkgee ordering to optimize how
    the mesh is stored.
  */

  class PolyGenerator : public MeshGenerator {

    std::function<RowMatrixXd(Points,std::vector<double>)> distance = dRectangle;

  private:

    // Maximum number of iteration in the generation loop
    int MaxIter = 1000;

    //! Compute Distant functions
    RowMatrixXd DistFnc(const Points & P, const std::vector<double> & BdBox) const {return distance(P, BdBox);};

    //! Specify boundary conditions
    VecMatrixXi BndryCnds(const Points & Node, const VectorD & BdBox) const;

    //! Generate Random Pointset
    Points PolyMshr_RndPtSet() const;

    //! Reflect PointSet
    Points PolyMshr_Rflct(const Points & P, double alpha) const;

    //! Compute Centroids and areas of polygon
    std::pair<Points, VectorD> PolyMshr_CntrdPly(const VecVecXi & Element, const Points & Node) const;

    //! Extract Node List
    std::pair<Points,VecVecXi> PolyMshr_ExtrNds(const Points & Node0,const VecVecXi & Element0) const;

    //! Collapse small edges
    std::pair<Points,VecVecXi> PolyMshr_CllpsEdgs(Points & Node0, VecVecXi & Element0, double Tol) const;

    //! Resequence Nodes
    std::pair<Points,VecVecXi> PolyMshr_RsqsNds(const Points & Node0, const VecVecXi & Element0) const;

    //!Rebuild List
    std::pair<Points,VecVecXi> PolyMshr_RbldLists(const Points & Node0, const VecVecXi & Element0, const VectorI & cNode) const;

    //!Construct Voronoi Diagram from a starting set of points (Using Qhull)
    std::pair<Points,VecVecXi> construct_voronoi(const Points & generators) const;

    //!Reordering algorithm
    std::vector<int> reverse_cuthill(const Eigen::SparseMatrix<int> & A) const;


  public:

    //! Constructor
    /*!
       @param domain_  The domain of the Mesh to be generated
       @param ne       Number of elements
    */

    PolyGenerator(const std::vector<double> & domain_, std::size_t ne) : MeshGenerator(domain_,ne) {}

    //! Default Assignment operator
    PolyGenerator & operator=(PolyGenerator const&)=default;

    //! Default Copy constructor
    PolyGenerator(PolyGenerator const &)=default;

    //! Default Move constructor
    PolyGenerator(PolyGenerator&&)=default;

    //! Default Move assignement
    PolyGenerator & operator=(PolyGenerator&&)=default;

    /** Specialised version. Modify by reference all the elements needed
        to construct a Mesh. It creates a polygonal mesh*/
    void generate(Mesh &) override;

    //! set distance function
    inline void set_distance(std::function<RowMatrixXd(Points,std::vector<double>)> dist_) {distance=dist_;}

    //! set number of iterations
    inline void set_max_iterations(int it_){MaxIter = it_;}

    //! virtual destructor
    virtual ~PolyGenerator()=default;



  };

}


#endif
