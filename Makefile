EXEC   = N-GenIC

OBJS   = main.o power.o allvars.o save.o read_param.o  read_glass.o  \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/polint.o nrsrc/trapzd.o

INCL   = allvars.h proto.h  nrsrc/nrutil.h  Makefile



#OPT   +=  -DPRODUCEGAS   # Set this to automatically produce gas particles 
                          # for a single DM species in the input file by 
			  # interleaved by a half a grid spacing


OPT   +=  -DMULTICOMPONENTGLASSFILE  # set this if the initial glass file
                                      # contains multiple components

OPT   +=  -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a 
                                      # transfer function that depends on
                                      # particle type

#OPT    +=  -DNO64BITID    # switch this on if you want normal 32-bit IDs
#OPT   +=  -DCORRECT_CIC   # only switch this on if particles are homogenously 
	                   # distributed over mesh cells (say glass)

OPT   +=  -DNEUTRINOS   # this will produce a second component as slight 
		        # neutrinos (needs to be in initial glass)
#OPT   +=  -DNEUTRINO_PAIRS  # this will produce an additional partner for
                             # every neutrino with opposite thermal velocities

OPT   += -DOUTPUT_DF     # turn this on to output the linear density field
                         # generated by N-GenIC


OPTIONS =  $(OPT)

#SYSTYPE="hpcf"
#SYSTYPE="cosmos"
#SYSTYPE="som"
#SYSTYPE="odyssey"
#SYSTYPE="Paco"
#SYSTYPE="zefiro"
SYSTYPE="flatiron"
#SYSTYPE="gordon"
#SYSTYPE="GordonS"


CC       =   mpicc        # sets the C-compiler (default)
OPTIMIZE =   -O2 -w -Wall    # optimization and warning flags (default)
MPICHLIB =  -lmpich

ifeq ($(SYSTYPE),"flatiron")
CC       =   mpicc     # sets the C-compiler   
OPT      +=  -DMPICH_IGNORE_CXX_SEEK
OPTIMIZE =   -std=gnu99 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function -march=native
GSL_INCL = -I$(GSL_BASE)/include
GSL_LIBS = -L$(GSL_BASE)/lib
FFTW_INCL= -I$(FFTW2_BASE)/include
FFTW_LIBS= -L$(FFTW2_BASE)/lib
MPICHLIB =  -lmpi
endif

ifeq ($(SYSTYPE),"GordonS")
CC = mpicc -g -O2 #-xW -ipo -Wall
CXX = mpiCC -g -O2 -xW -ipo -Wall
OPTIMIZE =
GMP_INCL = -I/opt/gnu/gmp/include
GMP_LIBS = -L/opt/gnu/gmp/lib
GSL_INCL = -I/opt/gsl/2.1/intel/include
GSL_LIBS = -L/opt/gsl/2.1/intel/lib
FFTW_INCL= -I/opt/fftw/2.1.5/intel/mvapich2_ib/include
FFTW_LIBS= -L/opt/fftw/2.1.5/intel/mvapich2_ib/lib
HDF5INCL = -I/opt/hdf5/intel/mvapich2_ib/include
HDF5LIB = -L/opt/hdf5/intel/mvapich2_ib/lib -lhdf5 -lz
endif

ifeq ($(SYSTYPE),"gordon")
CC = mpicc # sets the C-compiler
OPT += -DMPICH_IGNORE_CXX_SEEK
OPTIMIZE = -O3 -xHOST
GSL_INCL = -I/opt/gsl/2.1/intel/include
GSL_LIBS = -L/opt/gsl/2.1/intel/lib
FFTW_INCL= -I/opt/fftw/2.1.5/intel/openmpi_ib/include
FFTW_LIBS= -L/opt/fftw/2.1.5/intel/openmpi_ib/lib
MPICHLIB = -lmpi
HDF5INCL = -DH5_USE_16_API 
HDF5LIB = -lhdf5
endif

ifeq ($(SYSTYPE),"Paco")
CC       =  gcc     # sets the C-compiler
OPT      +=  
OPTIMIZE =   -O3 -g
GSL_INCL = -I/Users/fvillaescusa/Software/gsl-2.4_lib/include
GSL_LIBS = -L/Users/fvillaescusa/Software/gsl-2.4_lib/lib
FFTW_INCL= -I/Users/fvillaescusa/Software/fftw-2.1.5_lib/include
FFTW_LIBS= -L/Users/fvillaescusa/Software/fftw-2.1.5_lib/lib
MPICHLIB = -L/Users/fvillaescusa/Software/openmpi-3.0.0_lib/lib -lmpi
MPI_INCL = -I/Users/fvillaescusa/Software/openmpi-3.0.0_lib/include
endif

ifeq ($(SYSTYPE),"zefiro")
CC       =  mpicc     # sets the C-compiler
OPT      +=  
OPTIMIZE =  -O3 -g
GSL_INCL = -I/home/users/villaes/Libraries/GSL/gsl-1.16_lib/include 
GSL_LIBS = -L/home/users/villaes/Libraries/GSL/gsl-1.16_lib/lib
FFTW_INCL= -I/home/users/villaes/Libraries/FFTW/fftw-2.1.5_lib/include
FFTW_LIBS= -L/home/users/villaes/Libraries/FFTW/fftw-2.1.5_lib/lib
MPICHLIB = -L/usr/lib64/openmpi/include -lmpi
MPI_INCL = -I/usr/lib64/openmpi/lib
endif

ifeq ($(SYSTYPE),"odyssey")
CC       =  mpicc     # sets the C-compiler
OPT      +=  #-DNOCALLSOFSYSTEM  -DMPICH_IGNORE_CXX_SEEK
OPTIMIZE =   -O3 -g #-m64 -vec_report0 -xhost
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL = #-DH5_USE_16_API
HDF5LIB  = #-lhdf5
endif

ifeq ($(SYSTYPE),"cosmos")
CC       =  icc
OPTIMIZE =  
GSL_INCL =
GSL_LIBS =  
FFTW_INCL= -I/home/cosmos/users/mv249/Paco/Libraries/fftw-2.1.5_lib/include
FFTW_LIBS= -L/home/cosmos/users/mv249/Paco/Libraries/fftw-2.1.5_lib/lib64
MPICHLIB = -lmpi
endif

ifeq ($(SYSTYPE),"som")
CC       =  icc
OPTIMIZE =  -g -O2
MPI_INCL = -I/software/openmpi-1.6.2/intel/include
GSL_INCL =
GSL_LIBS =  
FFTW_INCL= -I/software/fftw-2.1.5/intel/include
FFTW_LIBS= -L/software/fftw-2.1.5/intel/lib
MPICHLIB = -L/software/openmpi-1.6.2/intel/lib -lmpi
endif


ifeq ($(SYSTYPE),"hpcf")
CC       =   mpicc  # sets the C-compiler
#OPT     +=  -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG 
OPTIMIZE =  -O3 -xHOST -ip -ipo
GSL_INCL = -I/usr/local/Cluster-Apps/gsl/1.9/include/
GSL_LIBS = -L/usr/local/Cluster-Apps/gsl/1.9/lib/
FFTW_INCL= -I/usr/local/Cluster-Apps.sandybridge/fftw/intel/2.1.5/double/include/
FFTW_LIBS= -L/usr/local/Cluster-Apps.sandybridge/fftw/intel/2.1.5/double/lib/
#GSL_INCL = -I/usr/local/Cluster-Users/mv249/GSL/include/
#GSL_LIBS = -L/usr/local/Cluster-Users/mv249/GSL/lib/
#FFTW_INCL= -I/usr/local/Cluster-Users/mv249/FFTW/include
#FFTW_LIBS= -L/usr/local/Cluster-Users/mv249/FFTW/lib -Wl,--rpath -Wl,/usr/local/Cluster-Users/mv249/FFTW/lib/
endif



FFTW_LIB =  $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw

LIBS   =   -lm  $(MPILIB)  $(FFTW_LIB)  $(GSL_LIBS)  -lgsl -lgslcblas

ifeq ($(SYSTYPE),"Solaris")
LIBS   =   -R/opt/local/lib/sparcv9 -lm  -lmpi   $(GSL_LIBS) -lgsl -lgslcblas  $(FFTW_LIB)
endif



CFLAGS =   $(OPTIONS)  $(OPTIMIZE)  $(FFTW_INCL) $(GSL_INCL) $(MPI_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS) $(MPICHLIB) -o  $(EXEC)  

$(OBJS): $(INCL)


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)


