#ifndef PARAMETERSHEADERDEF
#define PARAMETERSHEADERDEF

// MODEL PARAMETERS
#define I_0 9500
__constant__ extern double I;
#define mu 0.23
__constant__ extern double taux;
#define k 0.55
__constant__ extern double tShift;
#define G 33.75
#define alpha_0 0.05
#define beta 0.0075

#define pi M_PI

// FOR ADJUSTING DAY LENGTH
__constant__ extern double length_scaling;
__constant__ extern double length_shift;

// TIME STEPPER
#define timestep 0.001

// PERIOD OF ORBIT
#define Delta (24.0)

// FOR POLAR REPRESENTATION OF ORBIT
#define theta_step 0.001

// FOR ASSESSING ENTRAINMENT
#define entrain_thresh 0.01
//#define entrained_phase 5.48 // For B_max = 0.1
//#define entrained_phase 4.5379  // For B_max = 0.3
//#define entrained_phase 2.7274 // (B_max,N) = (0.3,8)
#define entrained_phase 6.1517 // (B_max,N) = (0.3,16)

#endif
