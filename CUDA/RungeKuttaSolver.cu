#include <cstdlib>
#include "RungeKuttaSolver.hpp"
#include "NonlinearProblem.hpp"
#include "parameters.hpp"

// HELPER FUNCTIONS FOR RUNGE-KUTTA CLASS
#define _add(a,b) AddDouble2(a,b)
#define _scale(a,b) ScaleDouble2(a,b)

__device__ double2 AddDouble2( double2 a, double2 b)
{
  a.x += b.x;
  a.y += b.y;
  return a;
}

__device__ double2 ScaleDouble2( double a, double2 b)
{
  b.x *= a;
  b.y *= a;
  return b;
}

__device__ RungeKuttaSolver::RungeKuttaSolver( double dt, NonlinearProblem* pProblem)
{
  mDt = dt;
  mpProblem = pProblem;
}

__device__ void RungeKuttaSolver::RungeKuttaStep( double t, double2& u)
{
  double2 f0, f1;

  // Take predictive step
  mpProblem->ComputeF( t, u, f0);

  // Take corrective step
  mpProblem->ComputeF( t+mDt, _add(u,_scale(mDt,f0)), f1);

  // Perform full step
  u = _add( u, _scale(mDt/2.0,_add(f0,f1)));

}
