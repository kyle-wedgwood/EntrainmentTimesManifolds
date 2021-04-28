#include "Kronauer.hpp"
#include "parameters.hpp"

__constant__ double I;
__constant__ double tShift;
__constant__ double taux;
__constant__ double length_scaling;
__constant__ double length_shift;

__device__ void Kronauer::ComputeF( double t, double2 u, double2& f)
{

  bool ft = (length_scaling+sin( 2.0*(pi/24.0)*(t-tShift-length_shift)) > 0);
  double alpha = alpha_0*sqrt(I/I_0);
  double B = G*alpha*ft*beta/(alpha+beta)*(1-0.4*u.y)*(1-0.4*u.x);

  f.x =
    (pi/12.0)*(mu*(u.x-(4.0/3.0)*u.x*u.x*u.x)-u.y*( (24.0/(0.99669*taux))*(24.0/(0.99669*taux))+k*B) );
  f.y = (pi/12.0)*(u.x+B);
}
