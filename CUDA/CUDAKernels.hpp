#ifndef CUDAKERNELSHEADERDEF
#define CUDAKERNELSHEADERDEF

__global__ void FindEntrainmentTimesKernel( const int2 dim,
                                            const double* pXMeshPts,
                                            const double* pYMeshPts,
                                            const double t_final,
                                            const double2* pRefOrbit,
                                            double* pResult);

__global__ void FindEntrainmentTimesPhaseKernel( const int2 dim,
                                                 const double* pXMeshPts,
                                                 const double* pYMeshPts,
                                                 const double t_final,
                                                 double* pResult);

__global__ void FindInsideFlagKernel( const int2 dim,
                                      const double* pXMeshPts,
                                      const double* pYMeshPts,
                                      const double t_final,
                                      const double* pRefOrbitAmp,
                                      double* pResult);
#endif
