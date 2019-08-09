module global
use iso_fortran_env, only : real64
use mpmodule
implicit none

private
public norm_definition,next_pair
public ControlData,init_Control,free_Control
public write_fform,FLD_DEFAULT,ACC_DEFAULT

integer,parameter :: FLD_DEFAULT = 60
integer,parameter :: ACC_DEFAULT = 40

type ControlData
type(mp_real) :: eps
character(8) :: calc_type
logical :: extended_L = .true.
integer :: nocc,nbas
integer :: omega,npair
type(mp_real) :: gfac,alpha
integer :: nC
type(mp_real),allocatable :: C(:)
integer :: MAX_ETA_POWER,MIN_ETA_POWER
type(mp_real) :: SCF_THR
integer :: SCF_MAXITER,SCF_DIIS_START,SCF_DIIS_SIZE
type(mp_real) :: CC_THR
integer :: CC_MAXITER,CC_DIIS_START,CC_DIIS_SIZE
integer :: IPRINT
end type ControlData

type ClistData
type(mp_real) :: C
type(ClistData),pointer :: next
end type ClistData

character(*),parameter :: digit     = '1234567890'
character(*),parameter :: other     = '+-.eEdD'
character(1),parameter :: separator = '/'

contains

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine norm_definition(n,norm)
implicit none
integer,intent(in) :: n
type(mp_real),intent(out) :: norm(0:n)
integer :: i
type(mp_real) :: tmp

norm(0) = 1/mpnrt(mppi(),4)
do i=1,n
   tmp = 2*i
   norm(i) = norm(i-1)/sqrt(tmp)
enddo
!!$do i=0,n
!!$   norm(i) = 1
!!$enddo

end subroutine norm_definition

function next_pair(i,j,omega) result(not_last)
implicit none
logical :: not_last
integer,intent(inout) :: i,j
integer,intent(in) :: omega

not_last = .true.

if(i<0.or.j<0) then
   i = 0
   j = 0
else
   i = i + 1
!!$   if(i>j.or.i+j>omega) then
   if(i>j) then
      i = 0
      j = j + 1
      if(j>omega) then
         i = -1
         j = -1
         not_last = .false.
      endif
   endif
endif

end function next_pair

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine init_Control(Control)
implicit none
type(ControlData) :: Control
real(real64) :: mp_digits
integer :: iunit,i,j,nC
type(ClistData),pointer :: Clist,Ctmp
character(100) :: string

mp_digits = mpwds*mpdpw

write(*,'(a,t50,i10)') 'Approximate number of digits:',int(mp_digits)
write(*,'(a,t50,i10)') 'Relative size real(mp_real)/real(real64):',mpwds6
write(*,'()')

Control%eps = mpreal(10._real64)**mpreald(-mp_digits)

open(newunit=iunit,file='HERMITE.INP',status='old')

read(iunit,*) Control%calc_type

read(iunit,*) Control%nocc
read(iunit,*) Control%nbas
read(iunit,*) Control%omega

read(iunit,'(a100)') string
Control%gfac  = read_number2(string)
read(iunit,'(a100)') string
Control%alpha = read_number2(string)

nC = 0
nullify(Clist)
do
   read(iunit,'(a100)') string
   if(index(string,'END')/=0) exit

   Ctmp => Clist
   allocate(Clist)
   Clist%next => Ctmp

   nC = nC + 1
   Clist%C = read_number2(string)
enddo

Control%nC = nC - 1
allocate(Control%C(0:Control%nC))
do while(nC>0)
   nC = nC - 1
   Control%C(nC) = Clist%C

   Ctmp => Clist%next
   deallocate(Clist)
   Clist => Ctmp
enddo

read(iunit,*) Control%MAX_ETA_POWER,Control%MIN_ETA_POWER

read(iunit,'(a100)') string
Control%SCF_THR = read_number2(string)
read(iunit,*) Control%SCF_MAXITER
read(iunit,*) Control%SCF_DIIS_START
read(iunit,*) Control%SCF_DIIS_SIZE

read(iunit,'(a100)') string
Control%CC_THR = read_number2(string)
read(iunit,*) Control%CC_MAXITER
read(iunit,*) Control%CC_DIIS_START
read(iunit,*) Control%CC_DIIS_SIZE

read(iunit,*) Control%IPRINT

close(iunit)

Control%npair = 0
i = -1; j = -1
do while(next_pair(i,j,Control%omega))
   Control%npair = Control%npair + 1
enddo

write(*,'(a)') Control%calc_type
write(*,'()')

write(*,'(a,i6)') 'nocc  = ',Control%nocc
write(*,'(a,i6)') 'nbas  = ',Control%nbas
write(*,'(a,i6,2x,a,i10)') 'omega = ',Control%omega,'-> npair = ',Control%npair

call write_fform(Control%gfac,string)
write(*,'(a)',advance='no') 'gfac  = '
write(*,'(a)') trim(string)
call write_fform(Control%alpha,string)
write(*,'(a)',advance='no') 'alpha = '
write(*,'(a)') trim(string)
do i=0,Control%nC
   call write_fform(Control%C(i),string)
   write(*,'(a,i2.2,a)',advance='no') 'C(',i,') = '
   write(*,'(a)') trim(string)
enddo

write(*,'()')

write(*,'(a,i6)') 'max eta power = ',Control%MAX_ETA_POWER
write(*,'(a,i6)') 'min eta power = ',Control%MIN_ETA_POWER

write(*,'()')

write(*,'(a,es12.3)') 'SCF_THR    = ',dble(Control%SCF_THR)
write(*,'(a,i6)') 'SCF_MAXITER    = ',Control%SCF_MAXITER
write(*,'(a,i6)') 'SCF_DIIS_START = ',Control%SCF_DIIS_START
write(*,'(a,i6)') 'SCF_DIIS_SIZE  = ',Control%SCF_DIIS_SIZE
write(*,'()')
write(*,'(a,es12.3)') 'CC_THR     = ',dble(Control%CC_THR)
write(*,'(a,i6)') 'CC_MAXITER     = ',Control%CC_MAXITER
write(*,'(a,i6)') 'CC_DIIS_START  = ',Control%CC_DIIS_START
write(*,'(a,i6)') 'CC_DIIS_SIZE   = ',Control%CC_DIIS_SIZE

write(*,'()')

write(*,'(a,i6)') 'PRINT          = ',Control%IPRINT

!!$block
!!$  type(mp_real) :: v
!!$  character(FLD_DEFAULT) :: s
!!$  integer :: i
!!$  v = -erfc(mpreal(9.d0))
!!$  call write_fform(v,s); write(*,'(a)') s
!!$  v = -v
!!$  call write_fform(v,s); write(*,'(a)') s
!!$  v = log(v)
!!$  call write_fform(v,s); write(*,'(a)') s
!!$  v = exp(-v/5)
!!$  call write_fform(v,s); write(*,'(a)') s
!!$  v = -100
!!$  call write_fform(v,s); write(*,'(a)') s
!!$  v = v + mpreal('4.e-41')
!!$  call write_fform(v,s); write(*,'(a)') s
!!$  v = mpreal('-5.e-41')
!!$  call write_fform(v,s); write(*,'(a)') s
!!$  do i=23,-44,-1
!!$     v = -(mpreal(10.d0)**i - mpreal('5.e-41'))
!!$     call write_fform(v,s)
!!$     write(*,'(i3,2x,a)') i,s
!!$enddo
!!$end block

end subroutine init_Control

subroutine free_Control(Control)
implicit none
type(ControlData) :: Control

deallocate(Control%C)

end subroutine free_Control

function read_number2(string) result(val)
implicit none
type(mp_real) :: val
character(*) :: string
integer :: ipos

string = adjustl(string)

ipos = verify(string,digit//other)
val = mpreal(string(1:ipos-1))

string(1:) = string(ipos:)
string = adjustl(string)

if(string(1:1)==separator) then

   string(1:1) = ' '
   string = adjustl(string)

   ipos = verify(string,digit//other)
   val = val/mpreal(string(1:ipos-1))

endif

end function read_number2

subroutine write_fform(val,sval,n1_IN,n2_IN)
implicit none
type(mp_real),intent(in) :: val
character(*),intent(out) :: sval
integer,intent(in),optional :: n1_IN,n2_IN
integer :: n1,n2
integer :: ipos,ipos_orig,mpos,t1,t2,power,offset,k
logical :: minus,minus_power,visible,carry
character(1),allocatable :: table(:)

if(present(n1_IN).and.present(n2_IN)) then
   n1 = n1_IN
   n2 = n2_IN
else
   n1 = FLD_DEFAULT
   n2 = ACC_DEFAULT
endif

if(n1>len(sval)) then
   write(*,'()') 'write_fform: sval is too short'
   stop
endif

if(n1<n2+3) then
   write(*,'()') 'write_fform: assure n1>=n2+3'
   stop
endif

sval = ''
ipos = n1-n2-1
sval(ipos:ipos) = '0'
ipos = n1-n2
sval(ipos:ipos) = '.'
do ipos=n1-n2+1,n1
   sval(ipos:ipos) = '0'
enddo

t2 = n1 + 5
t1 = t2 + 20
allocate(table(t1))
call mpeform(val,t1,t2,table)

if(table(1)=='-') then
   minus = .true.
   table(1) = ' '
   offset = 1
else
   minus = .false.
   offset = 0
endif
table(offset+2) = table(offset+1)
table(offset+1) = ' '
offset = offset + 1

ipos = offset + 1
do while(ipos<=t1)
   if(table(ipos)=='e') exit
   ipos = ipos + 1
enddo
table(ipos) = ' '
ipos = ipos + 1

if(table(ipos)=='-') then
   minus_power = .true.
   table(ipos) = ' '
   ipos = ipos + 1
else
   minus_power = .false.
endif
power = 0
do while(ipos<=t1)
   if(table(ipos)==' ') exit
   power = 10*power + mod(index(digit,table(ipos)),10)
   table(ipos) = ' '
   ipos = ipos + 1
enddo
if(minus_power) power = -power

ipos = n1-n2-power
if(power>=0) ipos = ipos - 1
ipos_orig = ipos

mpos = min(ipos-1,n1-n2-2)
visible = (merge(mpos,ipos,minus)>0)

if(visible) then
   
   k = offset + 1
   do while(ipos<=n1)
      if(ipos/=(n1-n2)) then
         sval(ipos:ipos) = table(k)
         k = k + 1
      endif
      ipos = ipos + 1
   enddo

   if(ipos==n1+1) then

      k = mod(index(digit,table(k)),10)
      select case(k)
      case(0:4)
         carry = .false.
      case(5)
         carry = (mod(index(digit,sval(n1:n1)),2)==1)
      case(6:9)
         carry = .true.
      end select

      if(carry) then
         ipos = ipos - 1
         do while(ipos>0)
            k = index(digit,sval(ipos:ipos))
            if(k==0) then
               if(sval(ipos:ipos)==' ') then
                  sval(ipos:ipos) = '1'
                  exit
               endif
            else
               k = mod(k,10) + 1
               sval(ipos:ipos) = digit(k:k)
               if(k<10) exit
            endif
            ipos = ipos - 1
         enddo
      endif

      mpos = min(ipos-1,mpos)
      visible = (merge(mpos,ipos,minus)>0)

   endif

endif

if(visible) then
   if(minus.and.min(ipos_orig,ipos)<=n1) sval(mpos:mpos) = '-'
else
   do ipos=1,n1
      sval(ipos:ipos) = '*'
   enddo
endif

deallocate(table)

end subroutine write_fform

end module global
