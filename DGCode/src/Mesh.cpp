/*!
    @file   Mesh.cpp
    @author Nicola Melas
    @brief  Implementation of class Mesh
*/

#include "Mesh.hpp"
#include "MeshGenerator.hpp"
#include "PolyGenerator.hpp"
#include "TriaGenerator.hpp"
#include <fstream>
#include <random>

namespace PolygonDG {

  Mesh::Mesh(const Point2D & p1, const Point2D & p2, std::size_t ne_x, std::size_t ne_y, std::string mesh_type){

    domain.resize(4);

    domain[0] = p1(0);
    domain[1] = p2(0);
    domain[2] = p1(1);
    domain[3] = p2(1);

    //Generation
    if (mesh_type == "Polygon") {
      PolyGenerator pg(domain, ne_x * ne_y);
      pg.generate(*this);
    }
    else if (mesh_type == "Triangle"){
      TriaGenerator tg(domain, ne_x, ne_y);
      tg.generate(*this);
    }
    else
      throw std::runtime_error("A correct mesh type must be selected");
    //computinh h_max
    compute_diameter();
  }

  void Mesh::compute_diameter() {

    auto fun = [](std::vector<double> BB) -> double {return  std::sqrt((BB[1]-BB[0])*(BB[1]-BB[0]) + (BB[3]-BB[2])*(BB[3]-BB[2]));};
    diameter = 0;
    for (std::size_t ii = 0; ii < ne; ii++ )
      diameter = std::max(diameter, fun(get_BBox(ii)));
  }

  void Mesh::print_mesh(std::ostream & out, bool print_polygons) const {

    std::size_t n_of_points = coord.size();
    out << "The Mesh is composed of " << ne << " elements and " << n_of_points << " points." << std::endl;

    out << "           X        Y        " << std::endl;

    for (std::size_t ii = 0; ii < n_of_points; ii++)
      out << "Point " << ii << " " << std::setprecision(4) << coord[ii].transpose() << std::endl;

    if (print_polygons) {
      for (std::size_t ii = 0; ii < ne; ii++)
        polygons[ii]->showMe();
    }

  }

  void Mesh::printvtk(std::string fileName) const {

    std::ofstream fout;

    if(fileName.substr(fileName.size() - 4, 4) != ".vtu")
      fout.open(fileName + ".vtu");
    else
      fout.open(fileName);


      // Create a vector with random integers in order to distinguish elements.
    std::vector<unsigned> elemValues;
    for(unsigned i = 0; i < polygons.size(); i++)
      elemValues.push_back(i);

    std::default_random_engine engine;
    std::shuffle(elemValues.begin(), elemValues.end(), engine);

    // Print the header
    fout << "<?xml version=\"1.0\"?>\n";
    fout << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    fout << "  <UnstructuredGrid>\n";

    fout << "    <Piece NumberOfPoints=\"" << coord.size() << "\" NumberOfCells=\"" << ne << "\">\n";

    fout << std::setprecision(10) << std::scientific;

    // Print the nodes coordinates.
    fout << "      <Points>\n";
    fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n         ";
    for(auto itNod = coord.cbegin(); itNod != coord.cend(); itNod++)
      fout << ' ' << itNod->x() << ' ' << itNod->y() << ' ' << 0;
    fout << "\n        </DataArray>\n";
    fout << "      </Points>\n";

    // Print the cells connectivity and type.
    fout << "      <Cells>\n";
    fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n         ";
    for(std::size_t i = 0; i < ne; i++) // Loop over elements
      for(int j = 0; j < connectivity[i].size(); j++) // Loop over vertices
        fout << ' ' << connectivity[i](j);
    fout << "\n        </DataArray>\n";

    fout << "         <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n         ";
    unsigned offset = 0;
    for(std::size_t i = 0; i < ne; i++) {
      offset += polygons[i]->size();
      fout << ' ' << offset;
    }
    fout << "\n        </DataArray>\n";

    fout << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n         ";
    for(unsigned i = 0; i < ne; i++)
      fout << " 5";
    fout << "\n        </DataArray>\n";
    fout << "      </Cells>\n";

    // Print a value for the Polygon.
    fout << "      <CellData Scalars=\"Mesh\">\n";
    fout << "        <DataArray type=\"UInt32\" Name=\"Mesh\" format=\"ascii\">\n         ";
    for(std::size_t i = 0; i < ne; i++)
      fout << ' ' << elemValues[i];
    fout << "\n        </DataArray>\n";
    fout << "      </CellData>\n";

    fout << "    </Piece>\n";


    fout << "  </UnstructuredGrid>\n";
    fout << "</VTKFile>" << std::endl;

    fout.close();

  }

}
