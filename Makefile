#################################################################################################
# Makefile for CHAPSim2, by Wei Wang, July 2021                                                 #
# Usage:                                                                                        #
#       make all        to make all files with -O2                                              #
#       make cfg=gnu    to debug for gfortran compiler                                          #
#       make cfg=intel  to debug for intel compiler                                             #
#       make cfg=cray   to make for cray                                                        #
# For debugging run:                                                                            #
# mpirun -np 4 valgrind --leak-check=full --track-origins=yes \                                 #
#                       --log-file=valgrind_output.txt ./CHAPSIM*                               #
#                          < solver_input > solver_output                                       #
#                                                                                               #
#################################################################################################

.PHONY: debug default clean
.SUFFIXES:

PROGRAM= CHAPSim

ifeq ($(cfg), gnu)
	FOPTS= -g \
		   -fbacktrace \
		   -fbounds-check \
		   -fcheck=all \
		   -ffpe-trap=invalid,zero,overflow \
		   -finit-real=snan -ftrapv \
		   -ffree-line-length-512 \
       -fallow-argument-mismatch -O3 \
		   -Wall
	FFLGS= -DOUBLE_PREC -DDEBUG
else ifeq ($(cfg), intel)
	FOPTS= -g -assume ieee_fpe_flags -check all -check bounds -check uninit -debug all \
	-fp-stack-check fpe0 -fpe3 -fpe-all=3 -ftrapuv -ftz -warn all, nounused
	FFLGS= -DOUBLE_PREC -DDEBUG
	FOPTS= -O3  -march=native  -fimplicit-none  -Wall  -Wline-truncation  -fwhole-file  -std=gnu
else ifeq ($(cfg), cray)
	FOPTS= # -m 3
	FFLGS= # -s default64
else ifeq ($(cfg), pg)
	FOPTS= -O3  -march=native  -Wall -fimplicit-none  -ffree-line-length-512  -fwhole-file  -std=gnu \
	-ffpe-trap=invalid,zero,overflow -fall-intrinsics
	FFLGS= -DOUBLE_PREC -pg
else
	FOPTS= -O3  -march=native  -Wall -fimplicit-none  -ffree-line-length-512  -fwhole-file  -std=gnu \
	-ffpe-trap=invalid,zero,overflow -fall-intrinsics
	FFLGS= -DOUBLE_PREC
endif


include ./lib/2decomp_fft/src/Makefile.inc
INCLUDE = -I ./lib/2decomp_fft/include
LIBS = -L ./lib/2decomp_fft/lib -l2decomp_fft

DIR_SRC= ./src
DIR_BIN= ./bin
DIR_OBJ= ./obj
DIR_MOD= ./mod

OBJS= mpi_mod.o\
      modules.o\
      tools_general.o\
      input_thermo.o\
      boundary_conditions.o\
      input_general.o\
      algorithms.o\
      operations.o\
      tools_solver.o\
      restart.o\
      geometry.o\
      domain_decomposition.o\
      poisson_interface.o\
      poisson.o\
      eq_continuity.o\
      eq_energy.o\
      eq_momentum.o\
      test_algrithms.o\
      flow_initialization.o\
      chapsim.o


default :
	@cd $(DIR_BIN)
	make $(PROGRAM) -f Makefile
	@mv *.mod $(DIR_MOD)
	@mv *.o $(DIR_OBJ)
	@mv $(PROGRAM) $(DIR_BIN)

$(PROGRAM): $(OBJS)
	$(F90) -o $@ $(OBJS) $(FOPTS) $(FFLGS) $(LIBS)

%.o : $(DIR_SRC)/%.f90
	$(F90) $(INCLUDE) $(FOPTS) $(FFLGS) $(F90FLAGS) -c $<

all:
	@make clean
	@cd $(DIR_BIN)
	make $(PROGRAM) -f Makefile
	@mv *.mod $(DIR_MOD)
	@mv *.o $(DIR_OBJ)
	@mv $(PROGRAM) $(DIR_BIN)

clean:
	@rm -f $(DIR_OBJ)/*.o $(DIR_BIN)/$(PROGRAM)
	@rm -f *.mod *.o $(DIR_SRC)/*.mod $(DIR_SRC)/*.o

