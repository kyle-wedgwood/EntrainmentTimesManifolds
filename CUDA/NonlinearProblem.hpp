#ifndef NONLINEARPROBLEMHEADERDEF
#define NONLINEARPROBLEMHEADERDEF

class NonlinearProblem
{

  public:

    __device__ virtual void ComputeF( double t, double2 u, double2& f) = 0;

};

#endif
