#ifndef KRONAUERHEADERDEF
#define KRONAUERHEADERDEF

#include "NonlinearProblem.hpp"

class Kronauer:
  public NonlinearProblem
{

  public:

    __device__ void ComputeF( double t, double2 u, double2& f);

  private:

};

#endif
