# COMPILER = /home/tgrining/programy/gcc-7.2.0/install
FC       = /usr/bin/gfortran  #$(COMPILER)/bin/gfortran
CC       = /usr/bin/gcc  #$(COMPILER)/bin/gcc
FFLAGS   = -O3 -fno-underscoring
CFLAGS   = -I/usr/include   #/home/mitek/programy/mpfr-libs/install/include

NAME = libmpfun.a

OBJ  = mpfuna.o mpfunf.o mpfung1.o mpmodule.o mpinterface.o

$(NAME) : $(OBJ)
	ar rvs $@ $?

%.o : %.mod
%.o : %.f90
	$(FC) -c $(FFLAGS) $<
%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean :
	rm *.o *.mod $(NAME)

mpfunf.o   : mpfuna.o
mpfung1.o  : mpfuna.o mpfunf.o
mpmodule.o : mpfuna.o mpfunf.o mpfung1.o
