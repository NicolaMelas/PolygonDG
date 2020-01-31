/*!
    @file   ElastodynamicsProblem.hpp
    @author Nicola Melas
    @brief  Class that defines the problem to solve
*/

#ifndef ASSEMBLEMATR
#define ASSEMBLEMATR

#include "fespace.hpp"
#include <string>
#include <vector>
#include "Gaussquadrature.hpp"
#include <utility>
#include <cmath>
#include <iterator>
#include <array>
#include "neighbour.hpp"
#include "muParserXInterface.hpp"
#include "BoundaryCondition.hpp"
#include "DirichletBC.hpp"
#include "NeumannBC.hpp"
#include "types.hpp"


namespace PolygonDG {
  /*!
      @brief Class that defines the elastodynamics problem to solve on
      polygonal grids. It stores all the matrixes. The two necessary to solve
      the system are M and A. @n
      The vector of BCs are stored to apply differently in time,
      advancing in time. Diri_tags are necessary to distinguish all the dirichlet
      boundaries from the Neumann ones. A solver is stored as a member to avoid
      factorization of the matrix at each temporal step. A boolean is used to
      control the need of this operation.@n
      The system to solve is: M* /ddot(u) + A*u = f
  */
  class ElastodynamicsProblem {

  private:

    //! Mass Matrix \f$ \int_{@Omega} (u . v ) dx \f$
    SparseMatrixXd M;
    //! \f$ \int_{E_h} penalty h_e^{-1} [v]*[u] ds \f$
    SparseMatrixXd S;
    //! \f$ \int_{E_h} {\sigma(v)} * [u]ds \f$
    SparseMatrixXd IT;
    //! Matrix given by \f$  \int_{\Omega} \sigma(u) \epsilon(v) dx \f$
    SparseMatrixXd V;
    //Stiffness Matrix A = M + V - IT - IT^T
    SparseMatrixXd A;
    //rhs
    VectorD f;
    FeSpace & femregion;
    const Coefficients & dati;
    std::vector<std::shared_ptr<BoundaryCondition>> BCs;
    std::vector<int> Diri_tags;
    Eigen::SparseLU<SparseMatrixXd> solver;
    bool solver_not_define = 1;

    //! Function that, given the 4 block matrices, stores them in a matrix M = [M1,M2;M3,M4]
    void finalize_matrix(SparseMatrixXd & M, const SparseMatrixXd & M1, const SparseMatrixXd & M2,
                                                const SparseMatrixXd & M3,const SparseMatrixXd & M4);

    //! Function that assemble the contribution given by the neigh elements
    void assemble_neigh(SparseMatrixXd& M,const VectorSt& row, const VectorI& neight,const VecSMatrixXd& M1,std::size_t nln,std::size_t n_edge);

    //! Function that checks the data
    void checkdata() const;
  public:

    //! Constructor.
    /*!
    @param fespace  Finite Element space on which the problem is defined
    @param coeff    Coefficients used to assemble the system
    @param BCs      The boundary conditions applied to the problem
    */
    ElastodynamicsProblem(FeSpace & fespace, const Coefficients & coeff, std::vector<std::shared_ptr<BoundaryCondition>> BCs);

    //! Default Assignment operator
    ElastodynamicsProblem & operator=(ElastodynamicsProblem const&)=default;

    //! Default Copy constructor
    ElastodynamicsProblem(ElastodynamicsProblem const &)=default;

    //! Default Move constructor
    ElastodynamicsProblem(ElastodynamicsProblem&&)=default;

    //! Default Move assignement
    ElastodynamicsProblem & operator=(ElastodynamicsProblem&&)=default;

    //! Assemble all the matrices and the rhs of the Elastodynamics problem
    void assemble(double time = 0.0);

    //! Update the system at the given time
    /*! It has to apply boundary condition updating the BC source at
    the given time */
    void update_rhs(double time);

    //! Return the Mass Matrix \f$ \int_{\Omega} (u * v ) dx \f$
    inline const SparseMatrixXd & get_M() const {return M;};

    //! Return the matrix given by \f$ \int_{\Omega}\sigma(u) \epsilon(v) dx \f$
    inline const SparseMatrixXd & get_V() const {return V;};

    //! \f$ \int_{E_h} \sigma(v) * [u]ds \f$
    inline const SparseMatrixXd & get_IT() const {return IT;};

    //! \f$ \int_{E_h} penalty  h_e^{-1} [v]*[u] ds \f$
    inline const SparseMatrixXd & get_S() const {return S;};

    //! Stiffness Matrix A = M + V - IT - IT^T
    inline const SparseMatrixXd & get_A() const {return A;};

    //! \f$ \int_{\Omega} (f * v) dx \f$ + B.C.
    inline const VectorD & get_f() const {return f;};

    //! Return the Adiacency property of the femregion of the problem
    inline const Neighbour & get_Neighbour() const {return femregion.get_Neighbour();};

    //! Return the L2 error
    /*! The vector u_h is reconstruct at the quadrature nodes where
    in the meantime the exact solution is computed. Then at those nodes the
    L2 error is computed.*/
    double L2_error(const VectorD & u_h,  MuParserInterface::muParserXInterface<3, std::array<double, 3>> & u_x,  MuParserInterface::muParserXInterface<3, std::array<double, 3>> & u_y, double time = 0.0) const;

    //! Return the H01 error
    /*! The vector u_h is reconstruct at the quadrature nodes where
    in the meantime the exact solution is computed. Then at those nodes the
    H01 error is computed.*/
    double H1_error(const VectorD & u_h,  MuParserInterface::muParserXInterface<3, std::array<double, 3>> & du_x,  MuParserInterface::muParserXInterface<3, std::array<double, 3>> & du_y,
                                          MuParserInterface::muParserXInterface<3, std::array<double, 3>> & dv_x,  MuParserInterface::muParserXInterface<3, std::array<double, 3>> & dv_y,
                                          double time = 0.0) const;

    //! If the user want to modify the right hand side with the time component multiplier
    inline void modify_time_rhs(double multiplier) {f *= multiplier;}

    //! Solve the ElastodynamicsProblem given the two initial guesses
    /*! Leap Frog scheme is used to discretize in time. It's a 2-step method */
    void solve(VectorD & u_h, const VectorD & u1old, const VectorD & u2old);

  };
}

#endif
