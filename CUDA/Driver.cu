/* Driver for use of finder operator to find isochrons
   Currently only works for planar models */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "EntrainmentTimeFinder.hpp"
#include "parameters.hpp"

void SetDayLengthPars( const double day_length,
                       double& host_length_scaling,
                       double& host_length_shift)
{
  host_length_scaling = sqrt(pow(sin(pi/12.0*day_length),
        2)/(2.0*(1.0-cos(pi/12.0*day_length))));

  if (((day_length < 12) & (host_length_scaling > 0)) |
      ((day_length > 12) & (host_length_scaling < 0)))
  {
    host_length_scaling *= -1;
  }
  host_length_shift = 12.0/pi*asin(host_length_scaling);

  std::cout << host_length_scaling << std::endl;
  std::cout << host_length_shift << std::endl;
}

int main(int argc, char* argv[])
{

  double host_I = 50.0;
  double host_tShift = 0.0;
  double host_taux = 24.2;
  double day_length = 4.0; // in hours

  double dI = 10.0;
  double dtShift = 0.0;
  double dtaux = 0.0;
  double dday_length = 4.0;
  double host_length_scaling, host_length_shift;

  int npts = 5;
  char filename[80];

  EntrainmentTimeFinder* p_finder = new EntrainmentTimeFinder();

  p_finder->SetDimensions(make_int2(1024, 1024));
  p_finder->SetGeometry(make_double2(-2.0, -2.0),
                         make_double2(2.0, 2.0));
  p_finder->SetFinalTime(20000.0);

  for (int i=0; i<npts; i++)
  {
    // Load reference orbit
    sprintf(filename, "orbits/stable_orbit_FS_I_%d_N_%d_taux_%.1f.dat", (int)
        (host_I), (int) day_length, host_taux);
    p_finder->LoadOrbit(filename);
    SetDayLengthPars(day_length, host_length_scaling, host_length_shift);

    // Copy parameters to device
    cudaMemcpyToSymbol(I, &host_I, sizeof(double));
    cudaMemcpyToSymbol(taux, &host_taux, sizeof(double));
    cudaMemcpyToSymbol(tShift, &host_tShift, sizeof(double));
    cudaMemcpyToSymbol(length_scaling, &host_length_scaling, sizeof(double));
    cudaMemcpyToSymbol(length_shift, &host_length_shift, sizeof(double));

    // Find the entrainment times
    p_finder->FindEntrainmentTimes();

    // Save the data
    sprintf(filename, "results/EntrainmentTimes_FS_I_%d_N_%d_taux_%.1f.dat",
        (int) host_I, (int) day_length, host_taux);

    p_finder->SaveData(filename);

    // Update for next loop
    host_I += dI;
    host_tShift += dtShift;
    host_taux += dtaux;
    day_length += dday_length;

    std::cout << "Done " << i+1 << " of " << npts << std::endl;
  }

  delete(p_finder);

  return 0;
}
