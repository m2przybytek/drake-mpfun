COMPILER = /home/mitek/programy/gcc-7.1.0/install
FC       = $(COMPILER)/bin/gfortran
FFLAGS   = -O3

NAME = libmpfun.a

OBJ  = mpmodule.o mpfuna.o mpfunbq.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfungq1.o

$(NAME) : $(OBJ)
	ar rvs $@ $?

%.o : %.mod
%.o : %.f90
	$(FC) -c $(FFLAGS) $<

.PHONY : clean
clean :
	rm *.o *.mod $(NAME)

mpmodule.o : mpfuna.o mpfunbq.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfungq1.o
mpfunbq.o  : mpfuna.o
mpfunc.o   : mpfuna.o mpfunbq.o
mpfund.o   : mpfuna.o mpfunbq.o mpfunc.o
mpfune.o   : mpfuna.o mpfunbq.o mpfunc.o mpfund.o
mpfunf.o   : mpfuna.o mpfunbq.o mpfunc.o mpfund.o mpfune.o 
mpfungq1.o : mpfuna.o mpfunbq.o mpfunc.o mpfund.o mpfune.o mpfunf.o 
