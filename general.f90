module general
use mpmodule
implicit none

integer,save :: G_nocc
integer,save :: G_nbas
integer,save :: G_nprim
integer,save :: G_npair
type(mp_real),save :: G_gfac
type(mp_real),save :: G_alpha

interface G_mp_to_real
module procedure mp_to_dble,mp_to_quad
end interface G_mp_to_real

interface G_real_to_mp
module procedure dble_to_mp,quad_to_mp
end interface G_real_to_mp

private mp_to_dble,mp_to_quad
private dble_to_mp,quad_to_mp

contains

function G_npair_total() result(cnt)
implicit none
integer :: cnt
integer :: i,j

i = -1
j = -1

cnt = 0
do while(G_next_pair(i,j))
   cnt = cnt + 1
enddo

end function G_npair_total

function G_next_pair(i,j) result(not_last)
implicit none
logical :: not_last
integer,intent(inout) :: i,j

not_last = .true.

if(i<0.or.j<0) then
   i = 0
   j = 0
else
   i = i + 1
   if(i>j.or.i+j>G_nprim) then
      i = 0
      j = j + 1
      if(j>G_nprim) then
         i = -1
         j = -1
         not_last = .false.
      endif
   endif
endif

end function G_next_pair

subroutine G_norm_definition(n,norm)
implicit none
integer,intent(in) :: n
type(mp_real),intent(out) :: norm(0:n)
integer :: i

norm(0) = 1/mpnrt(mppi(),4)
do i=1,n
   norm(i) = norm(i-1)/sqrt(mpreal(2.d0*i))
enddo
!norm=mpreal(1.d0)

end subroutine G_norm_definition

impure elemental subroutine mp_to_dble(val_mp,val)
implicit none
type(mp_real),intent(in) :: val_mp
real(8),intent(out) :: val
val = dble(val_mp)
end subroutine mp_to_dble

impure elemental subroutine dble_to_mp(val,val_mp)
implicit none
real(8),intent(in) :: val
type(mp_real),intent(out) :: val_mp
val_mp = mpreald(val)
end subroutine dble_to_mp

#if defined (mpfr)

impure elemental subroutine mp_to_quad(val_mp,val)
implicit none
type(mp_real),intent(in) :: val_mp
real(16),intent(out) :: val
write(*,*) 'qudruple precision not implemented in mpfun-mpfr'
stop
end subroutine mp_to_quad

impure elemental subroutine quad_to_mp(val,val_mp)
implicit none
real(16),intent(in) :: val
type(mp_real),intent(out) :: val_mp
write(*,*) 'qudruple precision not implemented in mpfun-mpfr'
stop
end subroutine quad_to_mp

#else

impure elemental subroutine mp_to_quad(val_mp,val)
implicit none
type(mp_real),intent(in) :: val_mp
real(16),intent(out) :: val
val = qreal(val_mp)
end subroutine mp_to_quad

impure elemental subroutine quad_to_mp(val,val_mp)
implicit none
real(16),intent(in) :: val
type(mp_real),intent(out) :: val_mp
val_mp = mpreal(val)
end subroutine quad_to_mp

#endif

subroutine my_mpfform(val,i1,i2,sval)
implicit none
type(mp_real),intent(in) :: val
integer,intent(in) :: i1,i2
character(*),intent(out) :: sval
character(1),allocatable :: table(:)
character(*),parameter :: digit = '0123456789'
integer :: ipos,n1,n2,offset,power,k
logical :: minus,carry

if(i1<i2+3) then
   write(*,*) 'Wrong format, assure i1>=i2+3'
   stop
endif

sval = ''
ipos = i1-i2-1
sval(ipos:ipos) = '0'
ipos = i1-i2
sval(ipos:ipos) = '.'
do ipos=i1-i2+1,i1
   sval(ipos:ipos) = '0'
enddo

n2 = i1 + 5
n1 = n2 + 20
allocate(table(n1))
call mpeform(val,n1,n2,table)

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
do while(ipos<=n1)
   if(table(ipos)=='e') exit
   ipos = ipos + 1
enddo
table(ipos) = ' '
if(table(ipos+1)=='-') then
   power = -1
   table(ipos+1) = ' '
   ipos = ipos + 1
else
   power = 1
endif
ipos = ipos + 1
power = power*(index(digit,table(ipos))-1)
table(ipos) = ' '
ipos = ipos + 1
do while(ipos<=n1)
   if(table(ipos)==' ') exit
   power = 10*power + sign(index(digit,table(ipos))-1,power)
   table(ipos) = ' '
enddo

k = 0
if(power>=0) then
   ipos = i1-i2-2-power
   if(ipos>0) then
      if(minus) sval(ipos:ipos) = '-'
      k = offset + 1
      do ipos=i1-i2-1-power,i1-i2-1
         sval(ipos:ipos) = table(k)
         k = k + 1
      enddo
      do ipos=i1-i2+1,i1
         sval(ipos:ipos) = table(k)
         k = k + 1
      enddo
   else
      do ipos=1,i1
         sval(ipos:ipos) = '*'
      enddo
   endif
else
   ipos = i1-i2-power
   if(ipos<=i1) then
      ipos = i1-i2-2
      if(minus) sval(ipos:ipos) = '-'
      k = offset + 1
      do ipos=i1-i2-power,i1
         sval(ipos:ipos) = table(k)
         k = k + 1
      enddo
   endif
endif

if(.not.(k==0.or.ipos>i1+1)) then
   k = index(digit,table(k))-1
   select case(k)
   case(0:4)
      carry = .false.
   case(5)
      carry = (mod(index(digit,sval(i1:i1)),2)==0)
   case(6:9)
      carry = .true.
   end select
   ipos = i1
   do while(carry)
      k = index(digit,sval(ipos:ipos))
      select case(k)
      case(0)
         select case(sval(ipos:ipos))
         case('.')
         case(' ')
            sval(ipos:ipos) = '1'
            carry = .false.
         case('-')
            sval(ipos:ipos) = '1'
            ipos = ipos - 1
            if(ipos>0) then
               sval(ipos:ipos) = '-'
            else
               do ipos=1,i1
                  sval(ipos:ipos) = '*'
               enddo
            endif
            carry = .false.
         end select
      case(1:9)
         k = k + 1
         sval(ipos:ipos) = digit(k:k)
         carry = .false.
      case(10)
         sval(ipos:ipos) = '0'
         carry = .true.
      end select
      ipos = ipos - 1
   enddo
endif

deallocate(table)

end subroutine my_mpfform

end module general
