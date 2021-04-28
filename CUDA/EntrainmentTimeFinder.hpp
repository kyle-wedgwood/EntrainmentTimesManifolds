#ifndef ENTRAINMENTTIMEFINDERHEADERDEF
#define ENTRAINMENTTIMEFINDERHEADERDEF

#include <iostream>

class EntrainmentTimeFinder
{
  public:

    ~EntrainmentTimeFinder();

    void FindEntrainmentTimes();

    void FindInsideFlag();

    void LoadOrbit( const char*);

    void LoadOrbitAmp( const char*);

    void SetDimensions( int2 dim);

    void SetGeometry( double2 xMin, double2 xMax);

    void SetFinalTime( double finalTime);

    void SaveData( const char*);

    // For debugging
    void SetDebugFlag( bool val);

  private:

    int CreateMesh();

    // Flags to show that parameters have been set
    bool mGeometryFlag = false;
    bool mTimeFlag = false;
    bool mDimFlag = false;

    // Arrays to store reference orbit
    double2* mpDev_refOrbit;
    double* mpDev_refOrbitAmp;

    // Arrays to store results
    double* mpHost_result;
    double* mpDev_result;

    // Mesh to compute points
    double* mpHost_xMeshPts;
    double* mpHost_yMeshPts;
    double* mpDev_xMeshPts;
    double* mpDev_yMeshPts;

    // Dimensions of mesh
    int2 mDim;
    double2 mXMin;
    double2 mXMax;

    // Time options
    double mFinalTime;

    // CUDA stuff
    unsigned int mNoThreads;
    unsigned int mNoBlocks;

    // For bookkeeping
    unsigned int* mpNoFinished;

    // Debugging
    bool mDebug;

};

#endif
