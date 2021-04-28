#ifndef RUNGEKUTTAHEADERDER
#define RUNGEKUTTAHEADERDER

#include <iostream>
#include "NonlinearProblem.hpp"

class RungeKuttaSolver
{

  public:

    __device__ RungeKuttaSolver( double dt, NonlinearProblem* pProblem);

    __device__ void RungeKuttaStep( double t, double2& u);

  private:

    double mDt;
    NonlinearProblem* mpProblem;

};

#endif
