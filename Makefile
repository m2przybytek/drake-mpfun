COMPILER  = $(HOME)/programy/gcc-7.2.0/install
FC        = $(COMPILER)/bin/gfortran
# -Wno-surprising -Wno-compare-reals
#FFLAGS    = -fcheck=all -Wall -Wextra -Wunused -Wunused-parameter -Warray-temporaries -Wrealloc-lhs -Wrealloc-lhs-all -pedantic -std=f2008
FFLAGS    = -O3
OPENMP    = -fopenmp
LFLAGS    = -Wl,-rpath -Wl,$(COMPILER)/lib64

MPFUNTYPE = mpfr
MPFUNINC  = -Impfun-$(MPFUNTYPE)
MPFUNLIB  = -Lmpfun-$(MPFUNTYPE) -lmpfun

ifeq ($(MPFUNTYPE),fort)
else ifeq ($(MPFUNTYPE),mpfr)
#MPFUNLIB += /home/mitek/programy/mpfr-libs/install/lib/libmpfr.a /home/mitek/programy/mpfr-libs/install/lib/libgmp.a
MPFUNLIB += -lmpfr
else
$(error MPFUNTYPE not defined, possible values: fort,mpfr)
endif

NAME = test
OBJ  = main.o integrals.o diis.o lproblem_mp.o eproblem_mp.o eproblem_real.o time.o inputread.o general.o

$(NAME) : $(OBJ)
	$(FC) -o $@ $^ $(LFLAGS) $(OPENMP) $(MPFUNLIB)

%.o : %.mod
%.o : %.f90
	$(FC) -c -cpp -D$(MPFUNTYPE) $(FFLAGS) $(OPENMP) $(MPFUNINC) $<

.PHONY : clean
clean :
	rm *.o *.mod $(NAME)

main.o      : integrals.o diis.o lproblem_mp.o eproblem_mp.o eproblem_real.o inputread.o general.o
integrals.o : time.o general.o
