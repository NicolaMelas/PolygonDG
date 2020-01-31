/*!
    @file   ExportFunction.cpp
    @author Nicola Melas
    @brief  Implementation of tool to ExportFunction
*/

#include "ExportFunction.hpp"
#include <fstream>
#include "Legendre.hpp"
#include "FeElement.hpp"

namespace PolygonDG {

  void exportfunctionvtk(const FeSpace & femregion, const VectorD & u_h, std::string fileName) {

    // PURPOSE:
    //
    //! Export the vector u_h in a vtk file. Export as solution of the problem

    if (int(2*femregion.get_ndof()) != u_h.size())
      throw std::runtime_error("The vector cannot belong to the given FeSpace");

    //Reconstruction values and points
    std::size_t ne = femregion.get_ne();
    std::size_t nln = femregion.get_nln();

    const Points & coord = femregion.get_mesh().get_coord();

    VecVecXi tria_connectivity;
    VectorSt partial_index(nln);
    VectorSt index(nln);
    VectorSt tmp(nln);

    //values of first and second component
    std::vector<double> u_points(coord.size(),0.0);
    std::vector<double> v_points(coord.size(),0.0);


    for (std::size_t ii = 0; ii < nln; ii++)
      tmp(ii) = ii;

    for (std::size_t ie = 0; ie < ne; ie++) {

      partial_index = VectorSt::Constant(nln, (ie)*nln);
      index = partial_index + tmp;

      VectorD local_uh1(nln);
      VectorD local_uh2(nln);

      // Reconstruction of the local solution
      for (std::size_t lll = 0; lll < nln; lll++) {
        local_uh1(lll) = u_h(index(lll));
        local_uh2(lll) = u_h(index(lll)+femregion.get_ndof());

      }


      const FeElement & element = femregion.get_FeElements()[ie];
      const auto & polygon = femregion.get_coords_element()[ie];
      std::size_t n_tria = element.get_number_tria();
      const auto & triangulations = element.get_triangulation();

      //FOr each element cycling the triangles
      for (std::size_t iTria = 0; iTria < n_tria; iTria++) {

        const auto tria = triangulations[iTria];
        const Points & ptria = tria->theVertexes();
        //find indexes to build the connectivity
        auto it1 = std::find(coord.cbegin(),coord.cend(),ptria[0]);
        auto it2 = std::find(coord.cbegin(),coord.cend(),ptria[1]);
        auto it3 = std::find(coord.cbegin(),coord.cend(),ptria[2]);
        VectorI tmpconn(3);
        tmpconn << it1 - coord.cbegin(),it2 - coord.cbegin(),it3 - coord.cbegin();
        //builting connectivity of triangles to export
        tria_connectivity.push_back(tmpconn);
        std::size_t p_size = 3;

        const MatrixXd & dphiq = BasisNodes(polygon, ptria, femregion.get_fem());

        for (std::size_t k = 0; k < p_size; k++) {

          double local_aprox1=0;
          double local_aprox2=0;

          for (std::size_t s = 0; s < nln ; s++) {  // reconstruct the discrete solution at the nodes
            local_aprox1 = local_aprox1 + dphiq(k,s)*local_uh1(s);
            local_aprox2 = local_aprox2 + dphiq(k,s)*local_uh2(s);
          }

          u_points[tmpconn[k]] = (u_points[tmpconn[k]]!=0) ? (u_points[tmpconn[k]]*0.5 + 0.5*local_aprox1) : local_aprox1 ;
          v_points[tmpconn[k]] = (v_points[tmpconn[k]]!=0) ? (v_points[tmpconn[k]]*0.5 + 0.5*local_aprox2) : local_aprox2 ;

        }
      }
    }

    std::ofstream fout;

    if(fileName.substr(fileName.size() - 4, 4) != ".vtu")
      fout.open(fileName + ".vtu");
    else
      fout.open(fileName);

    // Print the header
    fout << "<?xml version=\"1.0\"?>\n";
    fout << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    fout << "  <UnstructuredGrid>\n";

    fout << "    <Piece NumberOfPoints=\"" << coord.size() << "\" NumberOfCells=\"" << tria_connectivity.size() << "\">\n";

    fout << std::setprecision(10) << std::scientific;

    // Print the nodes coordinates.
    fout << "      <Points>\n";
    fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n         ";
    for(std::size_t ii = 0; ii < coord.size(); ii++)
      fout << ' ' << coord[ii](0) << ' ' << coord[ii](1) << ' ' << 0;
    fout << "\n        </DataArray>\n";
    fout << "      </Points>\n";

    // Print the cells connectivity and type.
    fout << "      <Cells>\n";
    fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n         ";
    for(std::size_t i = 0; i <tria_connectivity.size() ; i++) // Loop over elements
      for(int j = 0; j < tria_connectivity[i].size(); j++) // Loop over vertices
        fout << ' ' << tria_connectivity[i](j);
    fout << "\n        </DataArray>\n";

    fout << "         <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n         ";
    unsigned offset = 0;
    for(std::size_t i = 0; i < tria_connectivity.size(); i++) {
      offset += tria_connectivity[i].size();
      fout << ' ' << offset;
    }
    fout << "\n        </DataArray>\n";

    fout << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n         ";
    for(unsigned i = 0; i < tria_connectivity.size(); i++)
      fout << " 5";
    fout << "\n        </DataArray>\n";
    fout << "      </Cells>\n";

    // Print a value for the Polygon.
    fout << "      <PointData Scalars=\"u\">\n";
    fout << "        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n         ";
    for(std::size_t ii = 0; ii < u_points.size(); ii++)
      fout << ' ' << u_points[ii] << ' ' << v_points[ii] << ' ' << 0;
    fout << "\n        </DataArray>\n";
    fout << "      </PointData>\n";

    fout << "    </Piece>\n";


    fout << "  </UnstructuredGrid>\n";
    fout << "</VTKFile>" << std::endl;

    fout.close();

  }
}
