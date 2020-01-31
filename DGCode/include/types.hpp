/*!
    @file   types.hpp
    @author Nicola Melas
    @brief  File that defines the types used in all files
*/

#ifndef TYPES
#define TYPES

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/utility.h>

//! Generic Point
typedef Eigen::Matrix<double,2,1> Point2D;
//! std::Vector of points
typedef std::vector<Point2D> Points;
//! Dynamic-size Matrix of double
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
//! Dynamic-size Matrix of integers
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixXi;
//! Dynamic-size Vector of integers
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorI;
//! Dynamic-size Vector of std::size_t
typedef Eigen::Matrix<std::size_t, Eigen::Dynamic, 1> VectorSt;
//! Dynamic-size Vector of double
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorD;
//! std::vector of dynamic-size Matrix of double
typedef std::vector<MatrixXd> VecMatrixXd;
//! std::vector of dynamic-size Matrix of integers
typedef std::vector<MatrixXi> VecMatrixXi;
//! std::vector of dynamic-size vector of double
typedef std::vector<VectorD> VecVecXd;
//! std::vector of dynamic-size vector of integers
typedef std::vector<VectorI> VecVecXi;
//! Sparse version for MatrixXd
typedef Eigen::SparseMatrix<double> SparseMatrixXd;
//! Sparse version for VecMatrixXd
typedef std::vector<SparseMatrixXd> VecSMatrixXd;
//! Rowmajor version for MatrixXd
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXd;
//! Rowmajor version for MatrixXi
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXi;
//! Eigen triplet of double
typedef Eigen::Triplet<double> TripletD;

/*! CGAL Kernel*/
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//! 3D point CGAL
typedef Kernel::Point_3 Point_3;
//! Int indexes triplet
typedef CGAL::Triple<int, int, int> Triangle_int;

//! Vector of CGAL triplets
typedef Eigen::Matrix<Triangle_int, Eigen::Dynamic, 1> VectorTri;


#endif
