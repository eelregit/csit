EXEC   = 2LPTic

OBJS   = main.o power.o allvars.o save.o read_param.o  read_glass.o  \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/polint.o nrsrc/trapzd.o

INCL   = allvars.h proto.h  nrsrc/nrutil.h  Makefile



#OPT   +=  -DPRODUCEGAS   # Set this to automatically produce gas particles 
                         # for a single DM species in the input file by interleaved by a half a grid spacing


#OPT   +=  -DMULTICOMPONENTGLASSFILE  # set this if the initial glass file contains multiple components

#OPT   +=  -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a transfer function that depends on
                                     # particle type

#OPT   +=  -DNO64BITID    # switch this on if you want normal 32-bit IDs
#OPT   +=  -DCORRECT_CIC  # only switch this on if particles start from a glass (as opposed to grid)

#OPT += -DONLY_ZA # swith this on if you want ZA initial conditions (2LPT otherwise)


OPTIONS =  $(OPT)

SYSTYPE="ipmu-gw"
#SYSTYPE="popeye"

ifeq ($(SYSTYPE),"ipmu-gw")
#echo "Computing for ipmu-gw"
FFTW_INCL = -I/home/kazuyuki.akitsu/include
FFTW_LIBS = -L/home/kazuyuki.akitsu/lib
GSL_LIBS =  -L/home/kazuyuki.akitsu/lib/
GSL_INCL =  -I/home/kazuyuki.akitsu/include/gsl
CC       = mpicc
MPICHLIB = -L/opt/openmpi/1.10.7-gcc/lib -lmpi
OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)
endif

ifeq ($(SYSTYPE),"popeye")
#echo "Computing for popeye"
FFTW_INCL = -I/usr/include/
FFTW_LIBS = -L/usr/lib64/
GSL_LIBS = -L/usr/lib64/
GSL_INCL = -I/usr/include/gsl/
CC	 = mpicc
MPICHLIB = -L/lib64/openmpi/lib
#MPICHLIB = -L/usr/mpi/gcc/openmpi-2.1.2-hfi/lib64
endif

FFTW_LIB =  $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw

LIBS   =   -lm  $(MPICHLIB)  $(FFTW_LIB)  $(GSL_LIBS)  -lgsl -lgslcblas

CFLAGS =   $(OPTIONS)  $(OPTIMIZE)  $(FFTW_INCL) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



