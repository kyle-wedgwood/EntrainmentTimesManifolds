NVCC = nvcc
RM = rm -f
CUFLAGS  = -x cu --std=c++11 -arch=compute_30 -code=sm_30
LDFLAGS = -arch=compute_30 -code=sm_30 -L/usr/local/lib -L/usr/lib
CUINC  = -I. -dc -I/usr/local/include

OBJS = Driver.o \
       EntrainmentTimeFinder.o \
       Kronauer.o \
       RungeKuttaSolver.o \
       CUDAKernels.o

all: $(OBJS)
		$(NVCC) $(LDFLAGS) $(OBJS) -o Driver

%.o: %.cu
		$(NVCC) $(CUFLAGS) $(CUINC) $< -o $@

clean:
		rm -f *.o Driver
