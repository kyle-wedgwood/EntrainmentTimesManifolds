cmake_minimum_required( VERSION 2.8)

find_package( CUDA QUIET REQUIRED)

set( CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; --std=c++11 -dc -gencode arch=compute_30,code=sm_30)

if (APPLE)
  link_directories( /usr/local/lib)
else()
  link_directories( /usr/local/share/arrayfire-3/lib
                    /usr/local/lib)
                  set( CMAKE_INCLUDE_PATH
                    "/usr/local/share/arrayfire-3/include/;${CMAKE_INCLUDE_PATH}")
endif()

cuda_add_executable( Driver
                     parameters.hpp
                     StuartLandau.hpp
                     StuartLandau.cu
                     NonlinearProblem.cu
                     NonlinearProblem.hpp
                     KoopmanIsochron.cu
                     KoopmanIsochron.hpp
                     RungeKuttaSolver.cu
                     RungeKuttaSolver.hpp
                     CUDAKernels.cu
                     CUDAKernels.hpp
                     cu_error_functions
                     Driver.cu)
