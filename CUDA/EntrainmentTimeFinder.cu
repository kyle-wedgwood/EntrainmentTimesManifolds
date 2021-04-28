#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include "cu_error_functions.hpp"
#include "parameters.hpp"
#include "EntrainmentTimeFinder.hpp"
#include "CUDAKernels.hpp"

using namespace std;

EntrainmentTimeFinder::~EntrainmentTimeFinder()
{
  CUDA_CALL( cudaFree( mpDev_result));
  CUDA_CALL( cudaFree( mpDev_xMeshPts));
  CUDA_CALL( cudaFree( mpDev_yMeshPts));
  CUDA_CALL( cudaFree( mpDev_refOrbit));

  free( mpHost_result);
  free( mpHost_xMeshPts);
  free( mpHost_yMeshPts);
}

void EntrainmentTimeFinder::FindEntrainmentTimes()
{

  // Create meshes
  if (CreateMesh() > 0)
  {
    cout << "Simulation aborted." << endl;
    return;
  }

  // Transfer mesh data
  CUDA_CALL( cudaMemcpy( mpDev_xMeshPts, mpHost_xMeshPts,
        mDim.x*sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CALL( cudaMemcpy( mpDev_yMeshPts, mpHost_yMeshPts,
        mDim.y*sizeof(double), cudaMemcpyHostToDevice));

  cout << "Meshes copied to device." << endl;

  // Reset memory
  CUDA_CALL( cudaMemset( mpDev_result, 0.0, mDim.x*mDim.y*sizeof(double)));

  // Actually run the network
  cout << "Starting simulation..." << endl;

  FindEntrainmentTimesKernel<<<mNoBlocks,mNoThreads>>>( mDim,
                                                        mpDev_xMeshPts,
                                                        mpDev_yMeshPts,
                                                        mFinalTime,
                                                        mpDev_refOrbit,
                                                        mpDev_result);
  /*
  FindEntrainmentTimesPhaseKernel<<<mNoBlocks,mNoThreads>>>( mDim,
                                                             mpDev_xMeshPts,
                                                             mpDev_yMeshPts,
                                                             mFinalTime,
                                                             mpDev_result);
                                                             */
  CUDA_CHECK_ERROR();
  CUDA_CALL( cudaDeviceSynchronize());

  cout << "Fourier averages computed successfully." << endl;

  // Copy data back
  CUDA_CALL( cudaMemcpy( mpHost_result, mpDev_result,
        mDim.x*mDim.y*sizeof(double), cudaMemcpyDeviceToHost));

  cout << "Data copied to host." << endl;
}

void EntrainmentTimeFinder::FindInsideFlag()
{
  // Create meshes
  if (CreateMesh() > 0)
  {
    cout << "Simulation aborted." << endl;
    return;
  }

  // Transfer mesh data
  CUDA_CALL( cudaMemcpy( mpDev_xMeshPts, mpHost_xMeshPts,
        mDim.x*sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CALL( cudaMemcpy( mpDev_yMeshPts, mpHost_yMeshPts,
        mDim.y*sizeof(double), cudaMemcpyHostToDevice));

  cout << "Meshes copied to device." << endl;

  // Reset memory
  CUDA_CALL( cudaMemset( mpDev_result, 0.0, mDim.x*mDim.y*sizeof(double)));

  // Actually run the network
  cout << "Starting simulation..." << endl;
  FindInsideFlagKernel<<<mNoBlocks,mNoThreads>>>( mDim,
                                                  mpDev_xMeshPts,
                                                  mpDev_yMeshPts,
                                                  mFinalTime,
                                                  mpDev_refOrbitAmp,
                                                  mpDev_result);

  CUDA_CHECK_ERROR();
  CUDA_CALL( cudaDeviceSynchronize());

  cout << "Fourier averages computed successfully." << endl;

  // Copy data back
  CUDA_CALL( cudaMemcpy( mpHost_result, mpDev_result,
        mDim.x*mDim.y*sizeof(double), cudaMemcpyDeviceToHost));

  cout << "Data copied to host." << endl;
}

void EntrainmentTimeFinder::LoadOrbit( const char* filename)
{
  ifstream file;
  string line;
  std::cout << "Loading orbit from " << filename << std::endl;

  file.open( filename);
  int no_pts = (int)( 24.0/timestep);

  double temp_x;
  double temp_y;

  vector<double> x;
  vector<double> y;

  double2 p_hostRefOrbit[no_pts];

  while( getline( file, line))
  {
    istringstream iss( line);
    iss >> temp_x >> temp_y;
    x.push_back( temp_x);
    y.push_back( temp_y);
  }

  assert( x.size() == no_pts);

  // Put orbit into double2 array for easy copying
  for (int i=0; i<no_pts; i++)
  {
    p_hostRefOrbit[i].x = x[i];
    p_hostRefOrbit[i].y = y[i];
  }

  CUDA_CALL( cudaMalloc( &mpDev_refOrbit, no_pts*sizeof( double2)));
  CUDA_CALL( cudaMemcpy( mpDev_refOrbit, p_hostRefOrbit,
        no_pts*sizeof( double2), cudaMemcpyHostToDevice));

  cout << "Reference orbit loaded." << endl;
}

void EntrainmentTimeFinder::LoadOrbitAmp( const char* filename)
{
  ifstream file;
  string line;
  file.open( filename);

  int no_pts = (int) ( 2*pi/theta_step) + 1;
  double temp_theta;
  double temp_rho;
  double p_hostRefOrbit[no_pts];

  vector<double> rho;

  while( getline( file, line))
  {
    istringstream iss( line);
    iss >> temp_theta >> temp_rho;
    rho.push_back( temp_rho);
  }

  assert( rho.size() == no_pts);

  // Put orbit into double array for easy copying
  for (int i=0; i<no_pts; i++)
  {
    p_hostRefOrbit[i] = rho[i];
  }

  CUDA_CALL( cudaMalloc( &mpDev_refOrbitAmp, no_pts*sizeof( double)));
  CUDA_CALL( cudaMemcpy( mpDev_refOrbitAmp, p_hostRefOrbit,
        no_pts*sizeof( double), cudaMemcpyHostToDevice));

  cout << "Reference orbit loaded." << endl;
}

void EntrainmentTimeFinder::SetDimensions( int2 dim)
{
  mDim = dim;
  mDimFlag = true;
}

void EntrainmentTimeFinder::SetGeometry( double2 xMin, double2 xMax)
{
  mXMin = xMin;
  mXMax = xMax;
  mGeometryFlag = true;
}

void EntrainmentTimeFinder::SetFinalTime( double finalTime)
{
  mFinalTime = finalTime;
  mTimeFlag = true;
}

int EntrainmentTimeFinder::CreateMesh()
{
  if (!mDimFlag)
  {
    cout << "Dimensions not set. Aborting..." << endl;
    return 1;
  }
  if (!mGeometryFlag)
  {
    cout << "Geometry not set. Aborting..." << endl;
    return 2;
  }
  if (!mTimeFlag)
  {
    cout << "Simulation time not set. Aborting..." << endl;
  }

  // Allocate memory
  mpHost_xMeshPts = (double*) malloc( mDim.x*sizeof(double));
  mpHost_yMeshPts = (double*) malloc( mDim.y*sizeof(double));
  mpHost_result = (double*) malloc( mDim.x*mDim.y*sizeof(double));

  CUDA_CALL( cudaMalloc( &mpDev_xMeshPts, mDim.x*sizeof(double)));
  CUDA_CALL( cudaMalloc( &mpDev_yMeshPts, mDim.y*sizeof(double)));
  CUDA_CALL( cudaMalloc( &mpDev_result, mDim.x*mDim.y*sizeof(double)));

  double dx = (mXMax.x-mXMin.x)/(mDim.x-1);
  double dy = (mXMax.y-mXMin.y)/(mDim.y-1);
  for (int i=0;i<mDim.x;++i)
  {
    mpHost_xMeshPts[i] = mXMin.x+i*dx;
  }
  for (int i=0;i<mDim.y;++i)
  {
    mpHost_yMeshPts[i] = mXMin.y+i*dy;
  }

  mNoThreads  = 512;
  mNoBlocks = (mDim.x*mDim.y-1)/mNoThreads;

  cout << "Created mesh object with "
       << mDim.x
       << " x "
       << mDim.y
       << " points."
       << endl;

  return 0;
}

void EntrainmentTimeFinder::SaveData( const char* filename)
{
  ofstream file;
  file.open( filename);
  file << 0.0 << "\t";
  for (int i=0;i<mDim.x;++i)
  {
    file << mpHost_xMeshPts[i] << "\t"; // x mesh points
  }
  file << endl;

  for (int j=0;j<mDim.y;++j)
  {
    file << mpHost_yMeshPts[j] << "\t"; // y mesh points
    for (int i=0;i<mDim.x;++i)
    {
      file << mpHost_result[j*mDim.y+i] << "\t"; // Actual Fourier averages
    }
    file << endl;
  }
  file.close();
  cout << "Entrainment times saved to " << filename << endl;
}
