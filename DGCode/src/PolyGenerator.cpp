#include "MeshGenerator.hpp"
#include "PolyGenerator.hpp"
#include "OperatorLess.hpp"
#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include <array>


namespace PolygonDG {


  RowMatrixXd dRectangle(const Points & P, const std::vector<double> & BdBox)  {

    RowMatrixXd d = RowMatrixXd::Zero(P.size(),5);

    for (std::size_t i = 0; i < P.size(); i++){
      const VectorD & pp = P[i];

      // [x1-P(:,0), P(:,0)-x2, y1-P(:,1), P(:,1)-y2];
      d.block(i,0,1,4) <<  BdBox[0]-pp(0),-BdBox[1]+ pp(0),BdBox[2]- pp(1),-BdBox[3]+ pp(1);
    }

    d.block(0,4,P.size(),1) = (d.block(0,0,P.size(),4)).rowwise().maxCoeff();

    return d;
  }

  Points PolyGenerator::PolyMshr_RndPtSet() const {

    Points P(ne,Point2D::Zero());
    Points Y(ne, Point2D::Zero());

    const std::vector<double> & BdBox = domain;
    std::size_t Ctr = 0;

    // creating ne random points
    while (Ctr < ne) { //until ne points are built correctly
      for ( std::size_t jj = 0; jj < Y.size(); jj++) {
        Y[jj] << (BdBox[1]-BdBox[0]) * (VectorD::Random(1)* 0.5 + VectorD::Constant(1,0.5)) + VectorD::Constant(1,BdBox[0]),
                 (BdBox[3]-BdBox[2]) * (VectorD::Random(1)* 0.5 + VectorD::Constant(1,0.5)) + VectorD::Constant(1,BdBox[2]);
      }

      auto d = DistFnc(Y, domain);

      std::vector<int> tmp;
      for (int j = 0; j < d.rows(); j++) {
        //if negative, it means that the point is inside the domain
        if(d(j,d.cols()-1) < 0)
          tmp.push_back(j);
      }

      VectorI I = Eigen::Map<VectorI,Eigen::Unaligned>(tmp.data(), tmp.size());
      std::size_t Numadded = std::min(int(ne-Ctr), int(I.size()));
      for (std::size_t jj = Ctr; jj < Ctr + Numadded;jj++) //adding points to the variable to return
        P[Ctr+jj] = Y[I(jj)];
      Ctr+=Numadded;
    }
    return P;
  }



  Points PolyGenerator::PolyMshr_Rflct(const Points & P, double alpha) const{

    double epsilon = 1e-8;
    double eta = 0.9;
    auto d = DistFnc(P, domain);
    std::size_t NBdrySegs= d.cols()-1;

    //perturbations on first or second coordinate
    Point2D tmp1;
    tmp1 << epsilon,0;
    Point2D tmp2;
    tmp2 << 0,epsilon;

    //Relative distance
    auto dd = DistFnc(P+tmp1, domain);
    auto ddd = DistFnc(P+tmp2, domain);
    RowMatrixXd n1 = (dd-d)/epsilon;
    RowMatrixXd n2 = (ddd-d)/epsilon;

    Eigen::Matrix<bool, Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> I(d.rows(),NBdrySegs);
    for (int i = 0; i < d.rows(); i++) {
      for (std::size_t j = 0; j < NBdrySegs; j++) {
        //points close to the boundary are 1
        I(i,j) = (std::abs(d(i,j))<alpha);
      }
    }

    MatrixXd P1(P.size(),NBdrySegs);
    MatrixXd P2(P.size(),NBdrySegs);

    VectorD P_X(P.size());
    VectorD P_Y(P.size());

    for (auto it = P.begin(); it != P.end();it++){
      //save all x and y separately
      P_X(it - P.begin()) = (*it)(0);
      P_Y(it - P.begin()) = (*it)(1);
    }

    for (std::size_t kk = 0; kk < NBdrySegs; kk++) {
      //fill all the column with same vectors
      P1.block(0,kk,P1.rows(),1) = P_X;
      P2.block(0,kk,P1.rows(),1) = P_Y;
    }

    RowMatrixXi II = I.cast<int>();
    //number of points close to the boundary
    std::size_t tmpsize = II.sum();
    Points R_P(tmpsize,Point2D::Zero());
    int count = 0;

    for (int i = 0; i < I.rows(); i++) {
      for (int j = 0; j < I.cols(); j++) {
        if (I(i,j) == 1) {
          //reflection w.r.t. the boundary
          R_P[count] << P1(i,j)-2.0*n1(i,j)*d(i,j),P2(i,j)-2.0*n2(i,j)*d(i,j);
          count++;
        }
      }
    }

    auto d_R_P = DistFnc(R_P, domain);
    Eigen::Matrix<bool, Eigen::Dynamic,1> J(d_R_P.rows());

    count = 0;

    for (int i = 0; i < I.rows(); i++) {
      for (int j = 0; j < I.cols(); j++) {
        if (I(i,j) == 1) {
          //Points Reflected correctly
          J(count) = std::abs(d_R_P(count,d_R_P.cols()-1)) >= eta* std::abs(d(i,j)) && d_R_P(count,d_R_P.cols()-1) > 0;
          count++;
        }
      }
    }

    std::size_t i = 0;
    //deleting wrong reflection
    for (int count = 0; count < J.rows(); count++){
      if (J(count) != 1) {
        for (std::size_t kk = i; kk < R_P.size()-1; kk++ )
          R_P[i] = R_P[i+1];
        }
        else
          i++;
    }

    MatrixXi JJ = J.cast<int>();
    //deleting useless memory
    R_P.erase(R_P.begin()+JJ.sum(),R_P.end());
    //sorting and eliminating duplicates
    std::sort(R_P.begin(),R_P.end(),Eigenless<double>);
    auto it = std::unique(R_P.begin(), R_P.end());
    R_P.erase(it, R_P.end());

    return R_P;

  }


  std::pair<Points, VectorD> PolyGenerator::PolyMshr_CntrdPly(const VecVecXi & Element, const Points & Node) const {

    Points Pc(ne, Point2D::Zero());
    VectorD A = VectorD::Zero(ne,1);

    for (std::size_t el = 0; el < ne; el++) {

      std::size_t nv = Element[el].size();
      const auto & el_elem = Element[el];
      Polygon polygon;

      for (std::size_t jj = 0; jj < nv; jj++)
        polygon.addvertex(Node[el_elem(jj)]);

      Pc[el] = polygon.compute_centroid();
      A[el] = polygon.area();

    }

    return std::make_pair(Pc,A);

  }


  std::pair<Points,VecVecXi> PolyGenerator::PolyMshr_ExtrNds(const Points & Node0,const VecVecXi & Element0) const {

    // Eliminating duplicates and sorting to extract nodes. Then converting to
    // eigen sintax
    std::vector<int> map;
    for (std::size_t kk = 0; kk < ne; kk++) {
      std::vector<int> tmp(Element0[kk].data(),Element0[kk].data()+Element0[kk].size());
      map.insert(map.end(), tmp.begin(), tmp.end());
    }

    std::sort(map.begin(),map.end());
    auto unique = std::unique(map.begin(), map.end());
    map.erase(unique, map.end());
    VectorI mapp = Eigen::Map<VectorI,Eigen::Unaligned>(map.data(), map.size());

    std::vector<int> Cnode(Node0.size());

    for (std::size_t kk = 0; kk < Cnode.size(); kk++) {
      Cnode[kk] = kk;
    }
    // v are index of nodes not used for elements
    std::vector<int> v(Cnode.size());
    auto it = std::set_difference (Cnode.begin(), Cnode.end(), map.begin(), map.end() , v.begin());
    v.resize(it-v.begin());

    // those points are set to a fixed value
    for (std::size_t kk = 0; kk < v.size(); kk++) {
      Cnode[v[kk]] = map[map.size()-1];
    }
    //converting to eigen sintax
    VectorI cNode = Eigen::Map<VectorI, Eigen::Unaligned>(Cnode.data(), Cnode.size());
    VecVecXi Ee0 = Element0;
    //deleting not necessary elements
    Ee0.resize(ne);

    auto Node_Element = PolyMshr_RbldLists(Node0,Ee0,cNode);

    return std::make_pair(Node_Element.first, Node_Element.second);
  }

  std::pair<Points,VecVecXi> PolyGenerator::PolyMshr_CllpsEdgs(Points & Node0, VecVecXi & Element0, double Tol) const {

    bool boolean = true;

    VecVecXi Element;
    while (boolean) {
      MatrixXi cEdge;
      /* For each element searching for nodes too close*/
      for (std::size_t el = 0; el< Element0.size(); el++) {
        const auto & el_elem = Element0[el];
        auto elsz = el_elem.size();

        //Triangle cannot collapse
        if (elsz < 4)
          continue;

        VectorD vx(elsz);
        VectorD vy(elsz);

        for(int kk = 0; kk < elsz; kk++) {
          vx(kk) = Node0[el_elem(kk)](0);
          vy(kk) = Node0[el_elem(kk)](1);
        }

        VectorD beta(elsz), beta2(elsz);
        for(int kk = 0; kk < elsz; kk++) {
          //function to distinguish how close points are close
          beta(kk) = std::atan2(vy(kk)-vy.sum()/elsz, vx(kk)-vx.sum()/elsz);

          if (kk==0)
            beta2(elsz-1) = beta(kk);
          else
            beta2(kk-1) = beta(kk);
        }

        for(int kk = 0; kk < elsz; kk++) {
          double dd = std::fmod((beta2(kk)-beta(kk)), (2*M_PI));
          beta(kk) = dd + (dd < 0)*2.0*M_PI;
        }

        double betaIdeal = 2.0*M_PI/elsz;
        MatrixXi Edge(elsz,2);
        Edge.block(0,0,elsz,1) = el_elem;
        Edge(elsz-1,1) = el_elem(0);
        Edge.block(0,1,elsz-1,1) = el_elem.block(1,0,elsz-1,1);

        std::vector<int> useindex;
        for(int kk = 0; kk < beta.rows(); kk++){
          if (beta(kk) < Tol*betaIdeal)
            //Saving indexes with points very close
            useindex.push_back(kk);
        }

        std::size_t oldrows = cEdge.rows();
        cEdge.conservativeResize(oldrows + useindex.size(),std::max(cEdge.cols(),Edge.cols()));

        for(std::size_t kk = 0; kk < useindex.size(); kk++) {
          //Edge that has to be throw away
          cEdge.block(oldrows+kk,0,1,2) = Edge.block(useindex[kk],0,1,2);
        }
      }

      if(cEdge.size()==0)
        break;

      MatrixXi tmpEdge(cEdge.rows(), cEdge.cols());
      //smaller index on the left and bigger on the right
      tmpEdge.col(0) = cEdge.rowwise().minCoeff(); //min(cEdge(kk,0),cEdge(kk,1));
      tmpEdge.col(1) = cEdge.rowwise().maxCoeff(); //max(cEdge(kk,0),cEdge(kk,1));

      int count_duplicates = 0;

      //to delete duplicates
      for (int i = 0; i < tmpEdge.rows()-1-count_duplicates; i++){
        for(int j = i+1; j < tmpEdge.rows()-count_duplicates;j++) {
          if (tmpEdge(i,0) == tmpEdge(j,0) && tmpEdge(i,1) == tmpEdge(j,1)) {
            tmpEdge.block(j,0,tmpEdge.rows()-j-1,2) = tmpEdge.block(j+1,0,tmpEdge.rows()-j-1,2);
            count_duplicates++;
            j--;
          }
        }
      }

      tmpEdge.conservativeResize(tmpEdge.rows()-count_duplicates,2);

      std::vector<int> Cnode(Node0.size());
      for (std::size_t kk = 0; kk < Node0.size(); kk++) {
        Cnode[kk] = kk;
      }

      //Saving nodes without duplicates
      for (int i = 0; i < tmpEdge.rows(); i++) {
        //edge with same vertex as both ends, so nodes collapse
        Cnode[tmpEdge(i,1)] = Cnode[tmpEdge(i,0)];
      }
      //Converting to eigen sintax
      VectorI cNode = Eigen::Map<VectorI, Eigen::Unaligned>(Cnode.data(), Cnode.size());

      auto a_b =  PolyMshr_RbldLists(Node0,Element0,cNode);

      Node0 = a_b.first;
      Element0 = a_b.second;
    }

    return std::make_pair(Node0,Element0);

  }

  std::pair<Points,VecVecXi> PolyGenerator::PolyMshr_RsqsNds(const Points & Node0, const VecVecXi & Element0) const {

      //Resequence the nodes in a better order
      std::size_t NNode0 = Node0.size();
      std::size_t NElem0 = Element0.size();

      int nn = 0;
      for (std::size_t kk = 0; kk  < NElem0; kk++) {
        nn += Element0[kk].size()*Element0[kk].size();
      }

      VectorI i = VectorI::Zero(nn);
      VectorI j = VectorI::Zero(nn);
      VectorI s = VectorI::Ones(nn);
      int index = 0;

      for (std::size_t el = 0; el < NElem0; el++ ) {
        //preparing elements for reordering
        VectorI eNode = Element0[el];
        VectorI tmpi = Eigen::KroneckerProduct<VectorI,VectorI>(eNode, VectorI::Constant(Element0[el].size(),1));
        MatrixXi tmpj = Eigen::KroneckerProduct<VectorI,MatrixXi>(eNode, MatrixXi::Constant(1,Element0[el].size(),1));
        //middle step to convert to vector (for eigen)
        std::vector<int> tmpjj(tmpj.data(),tmpj.data()+tmpj.size());
        VectorI tmppp = Eigen::Map<VectorI,Eigen::Unaligned>(tmpjj.data(),tmpjj.size());
        i.block(index,0,eNode.size()*eNode.size(),1) = tmpi;
        j.block(index,0,eNode.size()*eNode.size(),1) = tmppp;
        index += eNode.size()*eNode.size();
      }

      //Reordering based on the matrix
      Eigen::SparseMatrix<int> K(NNode0, NNode0);
      for (int iii = 0; iii < i.size();iii++) {
          K.coeffRef(i(iii),j(iii)) = s(iii);
      }

      auto p = reverse_cuthill(K);


      VectorI cNode(p.size());
      for (std::size_t kk = 0; kk  < p.size(); kk++) {
        cNode(p[kk]) = kk;
      }

      return PolyMshr_RbldLists(Node0,Element0,cNode);
    }

  std::pair<Points,VecVecXi> PolyGenerator::PolyMshr_RbldLists(const Points & Node0, const VecVecXi & Element0, const VectorI & cNode) const {

    VecVecXi Element(Element0.size());
    std::vector<int> Cnode(cNode.data(),cNode.data()+cNode.size());
    std::vector<int> foo, ix, iy;
    foo = Cnode;
    std::sort(foo.begin(), foo.end());
    auto it = std::unique(foo.begin(), foo.end());
    foo.erase(it,foo.end());

    for (auto kk = foo.begin() ; kk != foo.end(); kk++) {
      auto itt = std::find(Cnode.begin(), Cnode.end(), *kk);
      if(itt != Cnode.end()) {
        ix.push_back(itt-Cnode.begin());
      }
    }

    for (auto kk = Cnode.begin() ; kk != Cnode.end(); kk++) {
      auto itt = std::find(foo.begin(), foo.end(), *kk);
      if(itt != foo.end()) {
        iy.push_back(itt-foo.begin());
      }
    }

    if(Node0.size() > ix.size())
      ix[ix.size()-1] = cNode.maxCoeff();


    Points Node(ix.size());
    for (std::size_t kk = 0; kk < ix.size(); kk++) {
      Node[kk] = Node0[ix[kk]];
    }

    for (std::size_t el = 0; el < Element0.size(); el++) {
      VectorI tmp = Element0[el];
      std::vector<int> tmpvec(tmp.size());
      for (int kk = 0; kk < tmp.size(); kk++)
        tmpvec[kk] = iy[tmp(kk)];
      std::sort(tmpvec.begin(), tmpvec.end());
      auto last = std::unique(tmpvec.begin(), tmpvec.end());
      tmpvec.erase(last, tmpvec.end());
      Element[el] = Eigen::Map<VectorI, Eigen::Unaligned>(tmpvec.data(), tmpvec.size());

      std::size_t nv = Element[el].size();
      VectorD vx(nv);
      VectorD vy(nv);
      const auto & el_elem = Element[el];

      for(std::size_t kk = 0; kk < nv; kk++) {
        vx(kk) = Node[el_elem(kk)](0);
        vy(kk) = Node[el_elem(kk)](1);
      }

      std::vector<double> beta(nv);
      for(std::size_t kk = 0; kk < nv; kk++) {
        beta[kk] = std::atan2(vy(kk)-vy.sum()/nv, vx(kk)-vx.sum()/nv);
      }

      std::vector<double> beta2 = beta;
      std::sort(beta2.begin(), beta2.end());
      std::vector<int> iix;
      for (auto kk = beta2.begin() ; kk != beta2.end(); kk++) {
        auto itt = std::find(beta.begin(), beta.end(), *kk);
        if(itt != beta.end()) {
          iix.push_back(itt-beta.begin());
        }
      }

      VectorI tmpp(nv);
      for (std::size_t kk = 0; kk < nv; kk++){
        tmpp(kk) = el_elem(iix[kk]);
      }

      Element[el] = tmpp;
    }
    return std::make_pair(Node, Element);
  }


  void PolyGenerator::generate(Mesh & mesh) {

    mesh.ne = ne;

    Points P = PolyMshr_RndPtSet();

    const std::vector<double> & BdBox = domain;
    // tolerance for collapsing edge
    double Tol = 5e-3;
    int It = 0;
    double Err = 1;
    //constant for the algorithm
    double c = 1.5;
    double Area = (BdBox[1]-BdBox[0])*(BdBox[3]-BdBox[2]);
    Points Pc = P;

    Points Node;
    VecVecXi Element;

    while(It<=MaxIter && (Err>Tol)) {

      double Alpha = c * std::sqrt(Area/ne);
      P = Pc;
      Points R_P = PolyMshr_Rflct(P, Alpha);
      Points tmp(P);
      std::move(R_P.begin(),R_P.end(),std::back_inserter(tmp));

      std::pair<Points,VecVecXi> voro = construct_voronoi(tmp);
      Node = voro.first; //vettore di vertici
      Element = voro.second;

      //used to check if the points of the mesh will be the centroid of the elements
      auto Pc_A = PolyMshr_CntrdPly(Element, Node);

      Pc = Pc_A.first;
      auto A = Pc_A.second;

      double area = (A.cwiseAbs()).sum();
      auto err1 = (A.cwiseProduct(A));

      Points PC_P(P.size());
      for (std::size_t it =0; it < P.size(); it++)
        PC_P[it] =  Pc[it] - P[it];

      auto PCPsq = cwiseprod(PC_P,PC_P);

      VectorD err2(PCPsq.size());

      //sum per row
      for (int k = 0; k < err2.size(); k++)
        err2(k) = PCPsq[k].sum();
      //scaled centroids error
      Err = std::sqrt((err1.cwiseProduct(err2)).sum())*ne/std::pow(area,1.5);
      It++;
    }

    std::pair<Points,VecVecXi> Node_Elem = PolyMshr_ExtrNds(Node,Element); //Extract Node List
    Node = Node_Elem.first;
    Element = Node_Elem.second;

    Node_Elem = PolyMshr_CllpsEdgs(Node,Element,0.1);  //Remove small edges
    Node = Node_Elem.first;
    Element = Node_Elem.second;

    Node_Elem = PolyMshr_RsqsNds(Node,Element); //Reorder nodes
    Node = Node_Elem.first;
    Element = Node_Elem.second;

    //Building mesh
    for (std::size_t i = 0; i < Element.size(); i++) {

      Points tmpcoordsel(Element[i].size());

      for(int j = 0; j < Element[i].size();j++)
        tmpcoordsel[j] = Node[Element[i](j)];

      mesh.polygons.push_back(std::make_shared<Polygon>(tmpcoordsel));
    }
    mesh.coord = Node;
    mesh.connectivity = Element;
  }

  std::pair<Points,VecVecXi> PolyGenerator::construct_voronoi(const Points & generators) const {

    //To convert the vector to be used in Qhull
    std::vector<double> points;
    for (auto it = generators.begin(); it != generators.end(); it++)
      std::copy(it->data(),it->data()+it->size(),std::back_inserter(points));

    orgQhull::Qhull qhull("",2,0.5*points.size(), points.data(),"v Qbb Qz o");

    std::stringstream oss;
    qhull.setOutputStream(&oss);
    if (qhull.hasOutputStream())
      qhull.outputQhull();

    auto a = qhull.facetList();
    auto b = a.toStdVector();
    Points Node(b.size()+1);
    double inf = std::numeric_limits<double>::infinity();
    Point2D first(2);
    first(0) = inf;
    first(1) = inf;
    Node[0] = first;

    for (std::size_t kk = 0; kk < b.size(); kk++) {
      Node[kk+1] << b[kk].getCenter()[0], b[kk].getCenter()[1];
    }

    //Reading anc copying the output to know information about index elements
    std::string line;
    std::getline(oss, line);
    std::getline(oss,line,' ');
    int lines_to_jump = std::stoi(line);

    std::getline(oss,line,' ');
    int n_region = std::stoi(line);
    VecVecXi Element(n_region);

    for (int hh = 0; hh <= lines_to_jump; hh++)
      std::getline(oss, line);

    for (int kk = 0; kk < n_region; kk++){

      std::getline(oss,line,' ');
      int size_el = std::stoi(line);
      VectorI tmpvec(size_el);

      if (size_el != 0) {
        for (int jj = 0; jj < size_el -1; jj++) {
          std::getline(oss,line,' ');
          int index = std::stoi(line);
          tmpvec(jj) = index;
        }
        std::getline(oss,line);
        int index = std::stoi(line);
        tmpvec(size_el-1) = index;
      }
      Element[kk] = tmpvec;
    }
    return std::make_pair(Node, Element);
  }

  std::vector<int> PolyGenerator::reverse_cuthill(const Eigen::SparseMatrix<int> & A) const
  {
    using namespace boost;
    using namespace std;
    typedef adjacency_list<vecS, vecS, undirectedS,
            property<vertex_color_t, default_color_type,
            property<vertex_degree_t,int> > > Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;

    typedef std::pair<std::size_t, std::size_t> Pair;

    //preparing the graph to be used in boost library
    std::vector<Pair> edges(A.nonZeros());
    int count = 0;
    for (int k=0; k<A.outerSize(); ++k)
      for (Eigen::SparseMatrix<int>::InnerIterator it(A,k); it; ++it)
      {
        edges[count]=Pair(it.row(),it.col());
        count++;
      }

    Graph G(A.rows());
    for (int i = 0; i < A.nonZeros(); ++i)
      add_edge(edges[i].first, edges[i].second, G);

    graph_traits<Graph>::vertex_iterator ui, ui_end;

    property_map<Graph,vertex_degree_t>::type deg = get(vertex_degree, G);

    for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
      deg[*ui] = degree(*ui, G);

    property_map<Graph, vertex_index_t>::type index_map = get(vertex_index, G);


    std::vector<Vertex> inv_perm(num_vertices(G));

    //reverse cuthill_mckee_ordering
    cuthill_mckee_ordering(G, inv_perm.rbegin(), get(vertex_color, G),
                           make_degree_map(G));

    //Fillin the result
    std::vector<int> result;
    result.reserve(inv_perm.size());
    for (std::vector<Vertex>::const_iterator i=inv_perm.begin();
       i != inv_perm.end(); ++i) {
      result.push_back(index_map[*i]);
    }
    return result;
  }

}
