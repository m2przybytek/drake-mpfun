FC       = /usr/bin/gfortran  #$(COMPILER)/bin/gfortran
# COMPILER  = /home/mitek/programy/gcc-7.2.0/install
# FC        = $(COMPILER)/bin/gfortran
FFLAGS    = -fcheck=all -Wall -Wextra -Wunused -Wunused-parameter -Warray-temporaries -Wrealloc-lhs -Wrealloc-lhs-all -Wno-target-lifetime -pedantic -std=f2008
#FFLAGS    = -O3
OPENMP    = -fopenmp
LFLAGS    = -Wl,-rpath -Wl,$(COMPILER)/lib64

MPFUNTYPE = mpfr
MPFUNINC  = -Impfun-$(MPFUNTYPE)
MPFUNLIB  = -Lmpfun-$(MPFUNTYPE) -lmpfun

ifeq ($(MPFUNTYPE),fort)
else ifeq ($(MPFUNTYPE),mpfr)
MPFUNLIB += -L/usr/lib -lmpfr -lgmp
else
$(error MPFUNTYPE not defined, possible values: fort,mpfr)
endif

NAME = test
OBJ  = hermite.o driver.o integrals.o integrals_core.o misc.o global.o \
       eproblem.o lproblem.o diis.o time.o
#integrals.o diis.o lproblem_mp.o eproblem_mp.o eproblem_real.o time.o general.o

$(NAME) : $(OBJ)
	$(FC) -o $@ $^ $(LFLAGS) $(OPENMP) $(MPFUNLIB)

%.o : %.mod
%.o : %.f90
	$(FC) -c $(FFLAGS) $(OPENMP) $(MPFUNINC) $<

.PHONY : clean
clean :
	rm *.o *.mod $(NAME)

hermite.o        : global.o driver.o
driver.o         : global.o misc.o integrals.o eproblem.o lproblem.o diis.o time.o
integrals.o      : global.o integrals_core.o misc.o time.o
integrals_core.o : global.o misc.o
misc.o           : global.o
#main.o      : integrals.o diis.o lproblem_mp.o eproblem_mp.o eproblem_real.o general.o
#integrals.o : time.o general.o
