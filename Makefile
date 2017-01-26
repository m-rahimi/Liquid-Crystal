PROG = chanel-midwayDSB 

OBJS =  module.o main.o input.o allocate_master.o initial.o functions.o visual.o parameter.o convert.o scater.o GL.o free_energy.o ginzburg_landau.o wtraj.o rtraj.o output.o Sbend.o

F90 = mpifort
F90FLAGS = -O3 #-check bounds #-O3 -fast
LDFLAGS = -O3  
LIBS =  
MKLPATH = /software/intel/mkl/lib/intel64
MKLINCLUDE = /software/intel/mkl/include
#MKLPATH = /opt/intel/mkl/lib/intel64
#MKLINCLUDE = /opt/intel/mkl/include

#OBJS = -Wl,--start-group
MKLLIB = -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_intel_thread.a $(MKLPATH)/libmkl_core.a -Wl,--end-group  -liomp5 -lpthread


all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) -L$(MKLPATH) -I$(MKLINCLUDE) $(MKLLIB) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

%.o : %.f90
	$(F90) $(F90FLAGS) -c $<
#-L/opt/mpi/lib/  -L/opt/intel/mkl/lib/intel64 -I/opt/intel/mkl/inlude  -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group  -liomp5 -lpthread 
