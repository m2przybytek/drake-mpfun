COMPILER = /home/mitek/programy/gcc-7.2.0/install
FC       = $(COMPILER)/bin/gfortran
FFLAGS   = -O3

NAME = libmpfun.a

OBJ  = mpfuna.o mpfunbq.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfungq1.o mpmodule.o

$(NAME) : $(OBJ)
	ar rvs $@ $?

%.o : %.mod
%.o : %.f90
	$(FC) -c $(FFLAGS) $<

.PHONY : clean
clean :
	rm *.o *.mod $(NAME)

mpfunbq.o  : mpfuna.o
mpfunc.o   : mpfuna.o mpfunbq.o
mpfund.o   : mpfuna.o mpfunbq.o mpfunc.o
mpfune.o   : mpfuna.o mpfunbq.o mpfunc.o mpfund.o
mpfunf.o   : mpfuna.o mpfunbq.o mpfunc.o mpfund.o mpfune.o 
mpfungq1.o : mpfuna.o mpfunbq.o mpfunc.o mpfund.o mpfune.o mpfunf.o
mpmodule.o : mpfuna.o mpfunbq.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfungq1.o
