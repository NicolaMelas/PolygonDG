/*!
    @file   Gaussquadrature.cpp
    @author Nicola Melas
    @brief  Implementation of class Gaussquadrature
*/

#include "Gaussquadrature.hpp"

namespace PolygonDG {

  std::pair<std::vector<double>,std::vector<double>> Gauss_Quadrature::GauLeg(double a, double b) {

    std::size_t m = (n+1)/2;
    double xm=0.5*(b+a);
    double xl=0.5*(b-a);
    std::vector<double> xx(n);
    std::vector<double> ww(n);

    double z = 0.0;
    double z1 = 0.0;
    double p1 = 0.0;
    double p2 = 0.0;
    double p3 = 0.0;
    double pp = 0.0;

    for(std::size_t i = 0; i < m; i++) {

      z=std::cos(M_PI*((i+1)-0.25)/(n+0.5));
      z1 = z+1;

      while(std::abs(z-z1) > std::numeric_limits<double>::epsilon()) {
        p1=1.0;
        p2=0.0;
        //Legendre poly computed in z
        for(std::size_t j = 0; j < n;j++) {
          p3=p2;
          p2=p1;
          p1=((2.0*(j+1)-1.0)*z*p2-((j+1)-1.0)*p3)/(j+1);
        }
        pp=n*(z*p1-p2)/(z*z-1.0);
        z1=z;
        z=z1-p1/pp;
      }

      xx[i]=xm-xl*z;
      xx[n-1-i]=xm+xl*z;
      ww[i]=2.0*xl/((1.0-z*z)*pp*pp);
      ww[n-1-i]=ww[i];

    }

    return std::make_pair(xx,ww);

  }

  std::pair<Points,std::vector<double>> Gauss_Quadrature::GLRefTria2D() {

    std::vector<double> w_2D(n*n);
    Points node_2D(n*n);

    std::pair<std::vector<double>,std::vector<double>> x_w_1D = GauLeg(-1,1);
    const std::vector<double> & x = x_w_1D.first;
    const std::vector<double> & w = x_w_1D.second;


    for (std::size_t i = 0;i < n; i++) {
      for (std::size_t j = 0;j < n; j++) {
        node_2D[n*i + j](0)=(1+x[i])/2;
        node_2D[n*i + j](1)=(1-x[i])*(1+x[j])/4;
        w_2D[n*i + j] = (1-x[i])*w[i]*w[j]/8;
      }
    }

    return std::make_pair(node_2D,w_2D);

  }

}
