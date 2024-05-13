# ######## GNU ########
F90 = mpif90
CC  = mpicc
FFLAGS = -O3 -Wall
#FFLAGS = -g -Wall
# -ffpe-trap=zero,invalid,underflow -fbacktrace
CFLAGS = -O3 -Wall
LDFLAGS = -O3

# FFTW flags
FFLAGS +=  -I$(FFTW_INC_DIR)
LDFLAGS += -L$(FFTW_LIB_DIR) -lfftw3_threads -lfftw3

#LDFLAGS      += -L$(ARMPL_LIBRARIES) -larmpl_mp -fopenmp
#FFLAGS       += -I$(ARMPL_INCLUDES)

SRCDIR = .
OBJ = \
	fftw_test_3d_multithread.o

fftw_test_3d_multithread: $(OBJ)
	$(F90) -o $@ $(OBJ) $(LDFLAGS) 

clean:
	rm -f *.o *.mod fftw_test_3d_multithread


# fortran rule
%.o:    $(SRCDIR)/%.f90
	$(F90) -cpp $(FFLAGS) -c $<

# C rule
%.o:    $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $<

