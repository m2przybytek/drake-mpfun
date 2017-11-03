COMPILER = /home/mitek/programy/gcc-7.1.0/install
FC       = $(COMPILER)/bin/gfortran
# -Wno-surprising -Wno-compare-reals
#FFLAGS   = -fcheck=all -Wall -Wextra -Wunused -Wunused-parameter -Warray-temporaries -Wrealloc-lhs -Wrealloc-lhs-all -pedantic -std=f2008
FFLAGS   = -O3
OPENMP   = -fopenmp
MPINC    = -Impfun
MPLIB    = -Lmpfun -lmpfun
LFLAGS   = -Wl,-rpath -Wl,$(COMPILER)/lib64

NAME = test
OBJ  = main.o integrals.o diis.o eproblem_mp.o eproblem_real.o general.o

$(NAME) : $(OBJ)
	$(FC) -o $@ $^ $(LFLAGS) $(OPENMP) $(MPLIB)

%.o : %.mod
%.o : %.f90
	$(FC) -c $(FFLAGS) $(OPENMP) $(MPINC) $<

.PHONY : clean
clean :
	rm *.o *.mod $(NAME)

main.o      : integrals.o diis.o eproblem_mp.o eproblem_real.o general.o
integrals.o : general.o
