/*!
    @file   Gaussquadrature.hpp
    @author Nicola Melas
    @brief  The class containt Gauss quadrature rule
*/
#ifndef GAUSSQUAD
#define GAUSSQUAD

#include <array>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <type_traits>
#include <tuple>
#include "Coefficients.hpp"
#include "types.hpp"

namespace PolygonDG {
	/*!
	    @brief  The class containt Gauss quadrature rule. There
			are two main function, one that return the 1D Gauss nodes with default
			 entries(0 and 1) and the other that return the 2D Gauss nodes(and weights)
			 on the reference triangle.
	*/
	class Gauss_Quadrature {

	private:
		//! degree of precision
		std::size_t n;

	public:

		//!Constructor
		//!@param		nn The precision degree
		Gauss_Quadrature(std::size_t nn) : n(nn) {};

		//! 1D Gauss Legendre nodes weigth
		std::pair<std::vector<double>,std::vector<double>> GauLeg(double a = 0, double b = 1);

		//! 2D (triangle) Gauss Legendere nodes weigth
		/**This routine computes the n^2 Gauss-Ledendre nodes and weights on the reference triangle
		 (0,0), (1,0), (0,1)
		*/
		std::pair<Points,std::vector<double>> GLRefTria2D();

	};
}


#endif
