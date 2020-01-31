/*!
    @file   Coefficients.hpp
    @author Nicola Melas
    @brief  Class that defines the coefficients of the problem
*/
#ifndef COEFFICIENTS
#define COEFFICIENTS

#include "Points_Utilities.hpp"
#include <string>
#include <vector>
#include "muParserXInterface.hpp"
#include "types.hpp"
#include <unordered_map>

namespace PolygonDG {
  /*!
      @brief Class that defines the coefficients of the problem.
      \f$ \rho \f$ \f$ \mu \f$ and \f$ \lambda \f$ are properties of the material, dt is the time
      discretization, in particular used in the solving procedure. The penalty
      coeff is used to apply correctly the SIP method.
  */

  class Coefficients {

  private:

    std::unordered_map<std::string,double> coeffs;
    MuParserInterface::muParserXInterface<3, std::array<double, 3>> & source_1;
    MuParserInterface::muParserXInterface<3, std::array<double, 3>> & source_2;

  public:

    //! Constructor
    Coefficients(const std::unordered_map<std::string, double> coeffs_,MuParserInterface::muParserXInterface<3, std::array<double, 3>> & source_1_,MuParserInterface::muParserXInterface<3, std::array<double, 3>> & source_2_ ) :
    coeffs(coeffs_), source_1(source_1_), source_2(source_2_) {};

    //! Default Assignment operator
    Coefficients & operator=(Coefficients const&)=default;

    //! Default Copy constructor
    Coefficients(Coefficients const &)=default;

    //! Default Move constructor
    Coefficients(Coefficients&&)=default;

    //! Default Move assignement
    Coefficients & operator=(Coefficients&&)=default;

    //! Density of the material
    inline const std::unordered_map<std::string,double> & get_coefficients() const {return coeffs;}

    //! First component of the source of the problem
    inline auto & get_source_1() const {return source_1;};

    //! Second component of the source of the problem
    inline auto & get_source_2() const {return source_2;};

  };
}

#endif
