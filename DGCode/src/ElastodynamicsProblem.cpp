/*!
    @file   ElastodynamicsProblem.cpp
    @author Nicola Melas
    @brief  Implementation of class ElastodynamicsProblem
*/
#include "ElastodynamicsProblem.hpp"
#include "EvalFunction.hpp"

namespace PolygonDG {

  ElastodynamicsProblem::ElastodynamicsProblem(FeSpace & femregion_, const Coefficients & dati_, std::vector<std::shared_ptr<BoundaryCondition>> Bc_) :
                  femregion(femregion_), dati(dati_), BCs(Bc_) {
    // TO distinguish in the assembling procedure Dirichlet from Neumann
    for (auto bc : BCs) {
      if (bc->get_type() == 1)
        Diri_tags.push_back(bc->get_tag());
    }

    //check if there are all needed coeff to solve the problem
    checkdata();

    //Getting Neighbour (by ref to modify it)
    Neighbour & neighbours = femregion.get_Neighbour();
    const auto & tmp_neigh = neighbours.get_neigh();
    const auto & tmp_conn = femregion.get_coords_element();

    //Cyclying all the element and the edges to set properly the tag of BC
    for (std::size_t jj = 0; jj < tmp_neigh.size(); jj++){

      const auto & tmp_neigh_ie = tmp_neigh[jj];
      const auto & tmp_poly_ie = tmp_conn[jj]->theVertexes();
      std::size_t tmp_size = tmp_neigh_ie.size();

      for (std::size_t kk = 0; kk < tmp_size; kk++) {
        // If positive, i.e. it is an internal element, it's useles to continue
        if (tmp_neigh_ie(kk) < 0) {

          const auto & p1 = tmp_poly_ie[kk];
          const auto & p2 = tmp_poly_ie[(kk+1)%tmp_size];
          for (auto bc : BCs) { //check that the given tag is negative
            if (bc->inside(p1) && bc->inside(p2) && bc->get_tag() < 0)
              neighbours.set_boundary_tag(jj,kk,bc->get_tag());
            else if (bc->inside(p1) && bc->inside(p2) && bc->get_tag() >= 0)
              throw std::runtime_error("The boundary tag cannot be positive. ElastodynamicsProblem has noot been constructed");
          }
        }
      }
    }

  }

  void ElastodynamicsProblem::assemble(double time) {

    //Elements containing all the necessary info
    const std::vector<FeElement> & elements = femregion.get_FeElements();

    double mu = dati.get_coefficients().at("mu");
    double lam = dati.get_coefficients().at("lam");
    double rho = dati.get_coefficients().at("rho");

    std::size_t ndof = femregion.get_ndof();
    std::size_t fem_degree = femregion.get_fem();
    std::size_t nln = femregion.get_nln();

    double penalty = 0.0;
    double penalty_scaled = 0.0;

    penalty = (fem_degree == 0) ? dati.get_coefficients().at("penalty") * (lam + 2*mu): dati.get_coefficients().at("penalty") * fem_degree * fem_degree * (lam + 2*mu);

    //Defining the matrices

    // \f$ \int_{\Omega} (u . v ) dx}
    SparseMatrixXd M1(ndof, ndof), M2(ndof, ndof), M3(ndof, ndof), M4(ndof, ndof);

    // \f$ \int_{\Omega} (sigma(u) eps(v) dx}
    SparseMatrixXd V1(ndof, ndof), V2(ndof, ndof), V3(ndof, ndof), V4(ndof, ndof);

    // \f$ \int_{E_h} {sigma(v)} . [u]ds
    SparseMatrixXd IT1(ndof, ndof), IT2(ndof, ndof), IT3(ndof, ndof), IT4(ndof, ndof);

    //  \f$ \int_{E_h} penalty  h_e^(-1) [v].[u] ds
    SparseMatrixXd S1(ndof, ndof), S2(ndof, ndof), S3(ndof, ndof), S4(ndof, ndof);


    VectorSt partial_index(nln);
    VectorSt index(nln);
    VectorSt tmp(nln);

    for (std::size_t ii = 0; ii < nln; ii++)
      tmp[ii] = ii;

    std::size_t ne = femregion.get_mesh().get_ne();
    const Neighbour & neighbour = femregion.get_Neighbour();

    for (std::size_t ie = 0; ie < ne; ie++) {

      partial_index = VectorSt::Constant(nln, (ie)*nln);
      index = partial_index + tmp;

      //neighbourhood info
      VectorI neigh_ie = neighbour.get_neigh()[ie];
      VectorI neighedges_ie = neighbour.get_neighedges()[ie];

      //Getting by reference all the necessary info
      const FeElement & element_ie = elements[ie];
      std::size_t n_tria = element_ie.get_number_tria();
      const Points & normals = element_ie.get_normal_edges();
      const std::vector<std::vector<double>> & ds_all = element_ie.get_ds();
      const std::vector<std::vector<double>> & dx_all = element_ie.get_dx();
      const std::vector<MatrixXd> & phi_trias = element_ie.get_phi_tria();
      const std::vector<MatrixXd> & phi_edges = element_ie.get_phi_edge();
      const std::vector<VecMatrixXd> & Grad_trias = element_ie.get_Grad_tria();
      const std::vector<VecMatrixXd> & Grad_edge = element_ie.get_Grad_edge();

      //for each triangle of the polygon
      for (std::size_t iTria = 0; iTria < n_tria; iTria++) {
        //Saving the scaled coeff for the integration inside this triangle
        const std::vector<double> & dx_Tria = dx_all[iTria];

        std::size_t dx_size = dx_Tria.size();
        //For each quadrature node
        for (std::size_t k = 0; k < dx_size; k++) {

          auto dx = dx_Tria[k];

          VectorD phi = VectorD::Zero(nln);
          VectorD grad_x = VectorD::Zero(nln);
          VectorD grad_y = VectorD::Zero(nln);

          //fill the values of basis and gradbasis at quadrature nodes
          for (std::size_t kk = 0; kk < nln; kk++) {
            phi(kk) = phi_trias[iTria](k,kk);
            grad_x(kk) = Grad_trias[iTria][kk](k,0);
            grad_y(kk) = Grad_trias[iTria][kk](k,1);
          }
          //For each degree of freedom in both directions(x and y)
          for (std::size_t ii = 0; ii < nln; ii++) {
            for (std::size_t jj = 0; jj < nln; jj++) {

              //Assembling Volume Matrices: phi_x*phi_x, phi_x*phi_y, phi_y*phi_x,phi_y*phi_y
              V1.coeffRef(index(ii), index(jj)) += ((lam+2*mu)*(grad_x(jj) * grad_x(ii)) + mu*(grad_y(jj) * grad_y(ii))) * dx;
              V2.coeffRef(index(ii), index(jj)) += (lam*(grad_y(jj) * grad_x(ii)) + mu*(grad_x(jj) * grad_y(ii)) ) * dx;
              V3.coeffRef(index(ii), index(jj)) += ( mu*(grad_y(jj) * grad_x(ii)) + lam*(grad_x(jj) * grad_y(ii))) * dx;
              V4.coeffRef(index(ii), index(jj)) += ((lam+2*mu)*(grad_y(jj) * grad_y(ii)) + mu*(grad_x(jj) * grad_x(ii))) * dx;

              //Thanks to simmetry and orthogonality, it's sufficient computing only this
              M1.coeffRef(index(ii), index(jj)) += (rho * phi(ii) * phi(jj)) * dx;


            }
          }
        }
      }

      std::size_t nedge_ie = femregion.get_nedge(ie);

      //Computing contributes affected by adiacency elements
      VecSMatrixXd ITN1(nedge_ie,SparseMatrixXd(nln,nln));
      VecSMatrixXd ITN2(nedge_ie,SparseMatrixXd(nln,nln));
      VecSMatrixXd ITN3(nedge_ie,SparseMatrixXd(nln,nln));
      VecSMatrixXd ITN4(nedge_ie,SparseMatrixXd(nln,nln));
      VecSMatrixXd  SN1(nedge_ie,SparseMatrixXd(nln,nln));
      VecSMatrixXd  SN4(nedge_ie,SparseMatrixXd(nln,nln));


      const std::vector<double> & penalties_edge = element_ie.get_penalties();

      for (std::size_t iedg = 0; iedg < nedge_ie; iedg++) {
        double penalty_edge = penalties_edge[iedg];
        penalty_scaled = penalty * penalty_edge;

        const std::vector<double> & ds_edge = ds_all[iedg];
        std::size_t ds_size = ds_edge.size();
        // for each quadrature edge node
        for (std::size_t k = 0; k < ds_size; k++){

          const auto & ds = ds_edge[k];

          MatrixXd Bedge = MatrixXd::Zero(nln,1);
          MatrixXd Gedge_x = MatrixXd::Zero(nln,1);
          MatrixXd Gedge_y = MatrixXd::Zero(nln,1);

          //fill the values of basis and gradbasis at quadrature nodes on the edge
          for (std::size_t kk = 0; kk < nln; kk++) {
            Bedge(kk) = phi_edges[iedg](k,kk);
            Gedge_x(kk) = Grad_edge[iedg][kk](k,0);
            Gedge_y(kk) = Grad_edge[iedg][kk](k,1);
          }

          double  aa = 0.5 * (lam+2*mu) * normals[iedg](0);
          double  ff = 0.5 * (lam+2*mu) * normals[iedg](1);
          double  bb = 0.5 * lam * normals[iedg](0);
          double  gg = 0.5 * lam * normals[iedg](1);
          double  ee = 0.5 * mu * normals[iedg](0);
          double  cc = 0.5 * mu * normals[iedg](1);

          for (std::size_t i=0; i < nln;i++) { // loop over scalar shape functions
            for (std::size_t j=0; j < nln;j++) { // loop over scalar shape functions
              //All edges but not Neumann ones
              if (std::find(Diri_tags.begin(), Diri_tags.end(), neigh_ie(iedg)) != Diri_tags.end() || neigh_ie(iedg) > -1) {
                S1.coeffRef(index(i), index(j)) +=  penalty_scaled * Bedge(j) * Bedge(i) * ds;
                S4.coeffRef(index(i), index(j)) +=  penalty_scaled * Bedge(j) * Bedge(i) * ds;
              }
              //only internal edges
              if (neigh_ie(iedg) > -1) {
                //values of the opposite edge. They have to be reversed due to reverse order of the nodes.
                const MatrixXd & phiedgeneigh = elements[neigh_ie[iedg]].get_phi_edge()[neighedges_ie[iedg]];
                IT1.coeffRef(index(i), index(j)) += ( aa * Gedge_x(i) * Bedge (j) + cc * Gedge_y(i) * Bedge(j)) * ds;
                IT2.coeffRef(index(i), index(j)) += ( ee * Gedge_y(i) * Bedge (j) + gg * Gedge_x(i) * Bedge(j)) * ds;
                IT3.coeffRef(index(i), index(j)) += ( bb * Gedge_y(i) * Bedge (j) + cc * Gedge_x(i) * Bedge(j)) * ds;
                IT4.coeffRef(index(i), index(j)) += ( ee * Gedge_x(i) * Bedge (j) + ff * Gedge_y(i) * Bedge(j)) * ds;

                (ITN1[iedg]).coeffRef((i), (j)) += - (aa * Gedge_x(i) * phiedgeneigh(ds_size-k-1,j) + cc * Gedge_y(i) * phiedgeneigh(ds_size-k-1,j)) * ds;
                (ITN2[iedg]).coeffRef((i), (j)) += - (ee * Gedge_y(i) * phiedgeneigh(ds_size-k-1,j) + gg * Gedge_x(i) * phiedgeneigh(ds_size-k-1,j)) * ds;
                (ITN3[iedg]).coeffRef((i), (j)) += - (bb * Gedge_y(i) * phiedgeneigh(ds_size-k-1,j) + cc * Gedge_x(i) * phiedgeneigh(ds_size-k-1,j)) * ds;
                (ITN4[iedg]).coeffRef((i), (j)) += - (ee * Gedge_x(i) * phiedgeneigh(ds_size-k-1,j) + ff * Gedge_y(i) * phiedgeneigh(ds_size-k-1,j)) * ds;

                (SN1[iedg]).coeffRef((i), (j)) += - penalty_scaled * Bedge(i) * phiedgeneigh(ds_size-k-1,j) * ds;
                (SN4[iedg]).coeffRef((i), (j)) += - penalty_scaled * Bedge(i) * phiedgeneigh(ds_size-k-1,j) * ds;

              }
              //only dirichlet edges
              else if (std::find(Diri_tags.begin(), Diri_tags.end(), neigh_ie(iedg)) != Diri_tags.end()) {

                IT1.coeffRef(index(i), index(j)) += 2.0 * ( aa * Gedge_x(i) * Bedge (j) + cc * Gedge_y(i) * Bedge(j)) * ds;
                IT2.coeffRef(index(i), index(j)) += 2.0 * ( ee * Gedge_y(i) * Bedge (j) + gg * Gedge_x(i) * Bedge(j)) * ds;
                IT3.coeffRef(index(i), index(j)) += 2.0 * ( bb * Gedge_y(i) * Bedge (j) + cc * Gedge_x(i) * Bedge(j)) * ds;
                IT4.coeffRef(index(i), index(j)) += 2.0 * ( ee * Gedge_x(i) * Bedge (j) + ff * Gedge_y(i) * Bedge(j)) * ds;

              }
            }
          }
        }
      }
      //Modyfyng the matrices with the neighborhood contribution
      assemble_neigh(IT1, index, neigh_ie, ITN1, nln, nedge_ie);
      assemble_neigh(IT2, index, neigh_ie, ITN2, nln, nedge_ie);
      assemble_neigh(IT3, index, neigh_ie, ITN3, nln, nedge_ie);
      assemble_neigh(IT4, index, neigh_ie, ITN4, nln, nedge_ie);

      assemble_neigh(S1, index, neigh_ie, SN1, nln, nedge_ie);
      assemble_neigh(S4, index, neigh_ie, SN4, nln, nedge_ie);

    }

    //Thanks to simmetry. Moreover M2 and M3 = 0
    M4 = M1;

    M.makeCompressed();
    V.makeCompressed();
    IT.makeCompressed();
    S.makeCompressed();

    //Setting from triplets M from [M1,M2;M3,M4]
    finalize_matrix(M,M1,M2,M3,M4);
    finalize_matrix(V,V1,V2,V3,V4);
    finalize_matrix(IT,IT1,IT2,IT3,IT4);
    finalize_matrix(S,S1,S2,S3,S4);

    //DEleting numerical zero
    M.prune(double(0),1e-10);
    V.prune(double(0),1e-10);
    IT.prune(double(0),1e-10);
    S.prune(double(0),1e-10);

    //Summing to obtain the whole Stiffness matrix
    A = V + S - IT - SparseMatrixXd(IT.transpose());

    //Assembling finally the right hand side
    //Evaluation of the source
    f = evaluate_function(femregion, dati.get_source_1(), dati.get_source_2(), time);
    //Adding BC contributions
    for (auto bc : BCs)
      bc->apply(f, femregion, time);


  }

  void ElastodynamicsProblem::finalize_matrix(SparseMatrixXd & M, const SparseMatrixXd & M1, const SparseMatrixXd & M2,
                                              const SparseMatrixXd & M3,const SparseMatrixXd & M4){
    M.resize(M1.rows()+M3.rows(),M1.cols()+M2.cols());
    M.makeCompressed();

    std::size_t ndof = M1.rows();
    std::vector<TripletD> tripletsM;

    for (int k=0; k<M1.outerSize(); ++k){
      for (SparseMatrixXd::InnerIterator it(M1,k); it; ++it)
      {
        tripletsM.push_back(TripletD(it.row(), it.col(), it.value()));
      }
    }

    for (int k=0; k<M2.outerSize(); ++k){
      for (SparseMatrixXd::InnerIterator it(M2,k); it; ++it)
      {
        tripletsM.push_back(TripletD(it.row(), it.col()+ndof, it.value()));
      }
    }

    for (int k=0; k<M3.outerSize(); ++k){
      for (SparseMatrixXd::InnerIterator it(M3,k); it; ++it)
      {
        tripletsM.push_back(TripletD(it.row() + ndof, it.col(), it.value()));
      }
    }

    for (int k=0; k<M4.outerSize(); ++k){
      for (SparseMatrixXd::InnerIterator it(M4,k); it; ++it)
      {
        tripletsM.push_back(TripletD(it.row()+ndof, it.col()+ndof, it.value()));
      }
    }

    M.setFromTriplets(tripletsM.begin(), tripletsM.end());

  }

  void ElastodynamicsProblem::update_rhs(double time) {
    //Modify only the rhs, only time dependent component
    f = evaluate_function(femregion, dati.get_source_1(), dati.get_source_2(), time);

    for (auto bc : BCs) {
      bc->apply(f, femregion, time);
    }
  }


  double ElastodynamicsProblem::L2_error(const VectorD & u_h,  MuParserInterface::muParserXInterface<3, std::array<double, 3>> & u_x,  MuParserInterface::muParserXInterface<3, std::array<double, 3>> & u_y, double time) const {

    // PURPOSE:
    //
    //! This routine computes the error in L2 norm, i.e

    //!  \f$ ||u-u_h||_{0,\Omega} \f$

    std::size_t ne = femregion.get_ne();
    std::size_t nln = femregion.get_nln();

    // initialization
    double E_L2 = 0;
    double t = time;

    VectorSt partial_index(nln);
    VectorSt index(nln);
    VectorSt tmp(nln);

    for (std::size_t ii = 0; ii < nln; ii++)
      tmp(ii) = ii;

    for (std::size_t ie = 0; ie < ne; ie++) {

      partial_index = VectorSt::Constant(nln, (ie)*nln);
      index = partial_index + tmp;

      VectorD local_uh1(nln);
      VectorD local_uh2(nln);

      //Reconstruncting the local u_h
      for (std::size_t lll = 0; lll < nln; lll++) {
        local_uh1(lll) = u_h(index(lll));
        local_uh2(lll) = u_h(index(lll)+femregion.get_ndof());

      }

      const FeElement & element = femregion.get_FeElements()[ie];

      std::size_t n_tria = element.get_number_tria();

      for (std::size_t iTria = 0; iTria < n_tria; iTria++) {
        //Getting point where to compute L2-error
        const Points & phys_tria = element.get_physical_points_tria()[iTria];
        const auto & dx_tria = element.get_dx()[iTria];
        std::size_t dx_size = dx_tria.size();
        const MatrixXd & dphiq = element.get_phi_tria()[iTria];

        for (std::size_t k = 0; k < dx_size; k++) {

          auto dx = dx_tria[k];

          auto x =  phys_tria[k](0);
          auto y =  phys_tria[k](1);

          std::array<double,3> xx;
          xx[0] = x;
          xx[1] = y;
          xx[2] = t;

          double local_exact1 = u_x(xx);
          double local_exact2 = u_y(xx);

          double local_aprox1=0;
          double local_aprox2=0;

          for (std::size_t s = 0; s < nln ; s++) {
            // reconstruct the discrete solution at the quadrature nodes
            local_aprox1 = local_aprox1 + dphiq(k,s)*local_uh1(s);
            local_aprox2 = local_aprox2 + dphiq(k,s)*local_uh2(s);
          }

          E_L2 += ((local_aprox1 - local_exact1)*(local_aprox1 - local_exact1))*dx + ((local_aprox2 - local_exact2)*(local_aprox2 - local_exact2))*dx;
        }
      }
    }
    return std::sqrt(E_L2);
  }

  double ElastodynamicsProblem::H1_error(const VectorD & u_h, MuParserInterface::muParserXInterface<3, std::array<double, 3>> & du_x,  MuParserInterface::muParserXInterface<3, std::array<double, 3>> & du_y,
                                                              MuParserInterface::muParserXInterface<3, std::array<double, 3>> & dv_x,  MuParserInterface::muParserXInterface<3, std::array<double, 3>> & dv_y,
                                                              double time) const {
    // PURPOSE:
    //
    //! This routine computes the error in H01 norm, i.e

    //!  \f$ ||u-u_h||_{1,\Omega} \f$

    std::size_t ne = femregion.get_ne();
    std::size_t nln = femregion.get_nln();

    // initialization
    double E_H1 = 0;
    double t = time;

    VectorSt partial_index(nln);
    VectorSt index(nln);
    VectorSt tmp(nln);

    for (std::size_t ii = 0; ii < nln; ii++)
      tmp(ii) = ii;

    for (std::size_t ie = 0; ie < ne; ie++) {

      partial_index = VectorSt::Constant(nln, (ie)*nln);
      index = partial_index + tmp;

      VectorD local_uh1(nln);
      VectorD local_uh2(nln);

      //Reconstruncting the local u_h
      for (std::size_t lll = 0; lll < nln; lll++) {
        local_uh1(lll) = u_h(index(lll));
        local_uh2(lll) = u_h(index(lll)+femregion.get_ndof());

      }

      const FeElement & element = femregion.get_FeElements()[ie];

      std::size_t n_tria = element.get_number_tria();

      for (std::size_t iTria = 0; iTria < n_tria; iTria++) {
        //Getting point where to compute L2-error
        const Points & phys_tria = element.get_physical_points_tria()[iTria];
        const auto & dx_tria = element.get_dx()[iTria];
        std::size_t dx_size = dx_tria.size();
        const auto & Grad = element.get_Grad_tria()[iTria];

        for (std::size_t k = 0; k < dx_size; k++) {

          auto dx = dx_tria[k];

          auto x =  phys_tria[k](0);
          auto y =  phys_tria[k](1);

          std::array<double,3> xx;
          xx[0] = x;
          xx[1] = y;
          xx[2] = t;

          double local_exact_grad1_x = du_x(xx);
          double local_exact_grad1_y = du_y(xx);
          double local_exact_grad2_x = dv_x(xx);
          double local_exact_grad2_y = dv_y(xx);

          double local_aprox_grad1_x = 0;
          double local_aprox_grad1_y = 0;
          double local_aprox_grad2_x = 0;
          double local_aprox_grad2_y = 0;

          for (std::size_t s = 0; s < nln ; s++) {
            // reconstruct the discrete solution at the quadrature nodes
            local_aprox_grad1_x = local_aprox_grad1_x + Grad[s](k,0)*local_uh1(s);
            local_aprox_grad1_y = local_aprox_grad1_y + Grad[s](k,1)*local_uh1(s);

            local_aprox_grad2_x = local_aprox_grad2_x + Grad[s](k,0)*local_uh2(s);
            local_aprox_grad2_y = local_aprox_grad2_y + Grad[s](k,1)*local_uh2(s);
          }

          E_H1 += ((local_aprox_grad1_x - local_exact_grad1_x)*(local_aprox_grad1_x - local_exact_grad1_x))*dx +
                  ((local_aprox_grad1_y - local_exact_grad1_y)*(local_aprox_grad1_y - local_exact_grad1_y))*dx +
                  ((local_aprox_grad2_x - local_exact_grad2_x)*(local_aprox_grad2_x - local_exact_grad2_x))*dx +
                  ((local_aprox_grad2_y - local_exact_grad2_y)*(local_aprox_grad2_y - local_exact_grad2_y))*dx;


        }

      }

    }

    return std::sqrt(E_H1);
  }


  void ElastodynamicsProblem::assemble_neigh(SparseMatrixXd& M,const VectorSt& row, const VectorI& neight,const VecSMatrixXd& M1,std::size_t nln,std::size_t n_edge){

    VectorD ones;
    VectorSt jv;
    VectorSt tmp(nln);

    for (std::size_t ii = 0; ii < nln; ii++)
        tmp(ii) = ii + 1;

    for(std::size_t iedg = 0; iedg < n_edge; iedg++) {
        if(neight(iedg) > -1) {

            jv = VectorSt::Constant(nln, (neight(iedg))*nln);
            jv = jv + tmp;

            for (int i = 0; i < row.rows(); i++) {
              for (int j = 0; j < jv.rows(); j++) {
                double tmpdouble = M.coeff(row(i),jv(j)-1) + M1[iedg].coeff(i,j);
                M.coeffRef(row(i),jv(j)-1) = tmpdouble;
              }
            }
        }
    }

  }

  //! Solve the ElastodynamicsProblem given the two initial guesses
  void ElastodynamicsProblem::solve(VectorD & u_h, const VectorD & u1old, const VectorD & u2old) {

    //to factorize the matrix only once
    if (solver_not_define) {
      M.makeCompressed();
      solver.compute(M);
    }

    double dt = dati.get_coefficients().at("dt");
    VectorD rhs = (2*M - dt*dt*A) * u1old - M*u2old + dt*dt*f;
    u_h = solver.solve(rhs);

  }

  //! Check the given coefficients
  void ElastodynamicsProblem::checkdata() const {

    const auto & coeffmap = dati.get_coefficients();
    if (coeffmap.size()<5)
      throw std::runtime_error("There are some parameters missing in the coefficients of the problem.");
    else if (coeffmap.find("rho")==coeffmap.end())
      throw std::runtime_error("rho coefficient is missing. The problem cannot be assembled.");
    else if (coeffmap.find("lam")==coeffmap.end())
      throw std::runtime_error("lam coefficient is missing. The problem cannot be assembled.");
    else if (coeffmap.find("mu")==coeffmap.end())
      throw std::runtime_error("mu coefficient is missing. The problem cannot be assembled.");
    else if (coeffmap.find("dt")==coeffmap.end())
      throw std::runtime_error("dt coefficient is missing. The problem cannot be assembled.");
    else if (coeffmap.find("penalty")==coeffmap.end())
      throw std::runtime_error("penalty coefficient is missing. The problem cannot be assembled.");

  }
}
