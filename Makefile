CC=g++

RM=/bin/rm -f

INTEL  = -mcmodel=medium
CFLAGS = $(INTEL) -fopenmp

# HEALPIX
HEALPIX     = $(HOME)/Healpix_3.40/src/cxx/basic_gcc
HEALPIX_INC = -I$(HEALPIX)/include
HEALPIX_LIB = -L$(HEALPIX)/lib -lhealpix_cxx -lcxxsupport -lsharp

INC = $(HEALPIX_INC)
LIB = $(HEALPIX_LIB)

EXE  = 3d_radtransfer.x
OBJ  = main.o

$(EXE): $(OBJ) Makefile
	$(CC) $(CFLAGS) -o $(EXE) $(OBJ) $(LIB) $(INC)

.cc.o: MAKEFILE
	$(CC) $(CFLAGS) $(INC) -o $*.o -c $*.cc
	
clean:
	$(RM) $(EXE) $(OBJ)


