program hermite
use mpmodule
use global
use driver
implicit none
type(ControlData) :: Control
type(mp_real) :: energSCF
type(mp_real),allocatable :: orbE(:)
type(mp_real),allocatable :: matC(:,:)
integer,allocatable :: parC(:)
type(mp_real),allocatable :: OOOO(:,:,:,:)

call init_Control(Control)

if(Control%nbas>=Control%nocc) then

   allocate(orbE(Control%nbas))
   allocate(matC(Control%nbas,Control%nbas),parC(Control%nbas))
   allocate(OOOO(Control%nocc,Control%nocc,Control%nocc,Control%nocc))

   call calc_SCF(energSCF,orbE,matC,parC,OOOO,Control)
   if(Control%omega>=0) then
      call calc_CC(energSCF,orbE,matC,parC,OOOO,Control)
   else
      write(*,'()')
      write(*,'(3a)') '!!! omega<0 : skipping ',trim(Control%calc_type),' !!!'
   endif

   deallocate(OOOO)
   deallocate(matC,parC)
   deallocate(orbE)

else
   write(*,'()')
   write(*,'(a)') '!!! nbas < nocc !!!'
endif

call free_Control(Control)

end program hermite
