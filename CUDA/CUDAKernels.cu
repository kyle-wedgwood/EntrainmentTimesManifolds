#include "CUDAKernels.hpp"
#include "RungeKuttaSolver.hpp"
#include "Kronauer.hpp"
#include "parameters.hpp"
#define norm(u,v) ((u.x-v.x)*(u.x-v.x) + (u.y-v.y)*(u.y-v.y))

__global__ void FindEntrainmentTimesKernel( const int2 dim,
                                            const double* pXMeshPts,
                                            const double* pYMeshPts,
                                            const double t_final,
                                            const double2* pRefOrbit,
                                            double* pResult)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index<dim.x*dim.y)
  {
    Kronauer* p_problem = new Kronauer();
    RungeKuttaSolver* p_solver = new RungeKuttaSolver( timestep, p_problem);

    // Initialise system
    int orbit_ind;
    double time = 0.0;
    double2 u;
    double2 ref_pt = pRefOrbit[0];
    u.x = pXMeshPts[index % dim.x];
    u.y = pYMeshPts[index / dim.x];

    do
    {
      time += timestep;

      p_solver->RungeKuttaStep( time, u);
      p_solver->RungeKuttaStep( time, ref_pt);

      //orbit_ind = (int) ( (time+24.0-tShift)/timestep) % (int) (24.0/timestep);
      //ref_pt = pRefOrbit[orbit_ind];

    } while ( (time<t_final) & (norm( u, ref_pt) > entrain_thresh*entrain_thresh));

    pResult[index] = time;

    delete( p_solver);
    delete( p_problem);
  }
}

__global__ void FindEntrainmentTimesPhaseKernel( const int2 dim,
                                                 const double* pXMeshPts,
                                                 const double* pYMeshPts,
                                                 const double t_final,
                                                 double* pResult)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index<dim.x*dim.y)
  {
    Kronauer* p_problem = new Kronauer();
    RungeKuttaSolver* p_solver = new RungeKuttaSolver( timestep, p_problem);

    // Initialise system
    double time = 0.0;
    double2 u;
    double2 u_old;
    u.x = pXMeshPts[index % dim.x];
    u.y = pYMeshPts[index / dim.x];

    do
    {
      time += timestep;

      u_old = u;
      p_solver->RungeKuttaStep( time, u);

    } while ( (time<t_final) & (( abs( fmod( time-tShift-24.0, 24.0) - entrained_phase) > 0.01) || (u.x > 0.0) || (u_old.x < 0.0)));

    pResult[index] = time;

    delete( p_solver);
    delete( p_problem);
  }
}

__global__ void FindInsideFlagKernel( const int2 dim,
                                      const double* pXMeshPts,
                                      const double* pYMeshPts,
                                      const double t_final,
                                      const double* pRefOrbitAmp,
                                      double* pResult)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index<dim.x*dim.y)
  {
    Kronauer* p_problem = new Kronauer();
    RungeKuttaSolver* p_solver = new RungeKuttaSolver( timestep, p_problem);

    // Initialise system
    int orbit_ind;
    double time = 0.0;
    double2 u;
    bool inside_flag = false;
    u.x = pXMeshPts[index % dim.x];
    u.y = pYMeshPts[index / dim.x];

    do
    {
      time += timestep;

      p_solver->RungeKuttaStep( time, u);

      orbit_ind = (int) ( ( atan2( u.y, u.x)+pi)/theta_step);
      inside_flag = ( sqrtf( u.x*u.x+u.y*u.y) < pRefOrbitAmp[orbit_ind]);

    } while ( (time<t_final) && (!inside_flag));

    pResult[index] = (double) inside_flag;

    delete( p_solver);
    delete( p_problem);
  }
}
