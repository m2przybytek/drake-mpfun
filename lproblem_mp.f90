module lproblem_mp
use mpmodule
implicit none

private
public symU_linearsolve_mp

contains

subroutine symU_linearsolve_mp(n,A,b)
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:),b(:)
integer,allocatable :: perm(:)

if(n<0) then

   write(*,'(a)') 'Negative matrix dimension in linear equations solver!'
   stop

elseif(n==0) then

   write(*,'(a)') 'WARNING!!! Zero matrix dimension in linear equations solver!'
   return

elseif(n==1) then

   b(1) = b(1)/A(1,1)

else

   allocate(perm(n))

   call LDL_factor(n,A,perm)
   call LDL_solve(n,A,perm,b)

   deallocate(perm)

endif

end subroutine symU_linearsolve_mp

subroutine LDL_factor(n,A,perm)
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:)
integer,intent(out) :: perm(:)
integer :: i,j,k,k_max
type(mp_real) :: D,rtmp
logical :: do_warn

do i=1,n
   perm(i) = i
enddo

do_warn = .true.
do k=n,1,-1

   rtmp  = abs(A(k,k))
   k_max = k
   do i=1,k-1
      if(abs(A(i,i))>rtmp) then
         rtmp  = abs(A(i,i))
         k_max = i
      endif
   enddo
   if(A(k_max,k_max)<=0.and.do_warn) then
      write(*,'(a)') 'WARNING!!! &
           &Matrix in linear solver is not positive definite!'
      do_warn = .false.
   endif

   call exchange(n,A,perm,k_max,k)

   D = 1/A(k,k)

   do j=1,k-1
      rtmp   = A(j,k)
      A(j,k) = A(j,k)*D
      do i=1,j
         A(i,j) = A(i,j) - A(i,k)*rtmp
      enddo
   enddo

   A(k,k) = D

enddo

end subroutine LDL_factor

subroutine LDL_solve(n,A,perm,x)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: A(:,:)
integer,intent(in) :: perm(:)
type(mp_real),intent(inout) :: x(:)
integer :: i,j
type(mp_real) :: tmp

do j=n,1,-1
   tmp = x(perm(j))
   do i=1,j-1
      x(perm(i)) = x(perm(i)) - A(i,j)*tmp
   enddo
   x(perm(j)) = A(j,j)*tmp
enddo

do j=2,n
   tmp = x(perm(j))
   do i=1,j-1
      tmp = tmp - x(perm(i))*A(i,j)
   enddo
   x(perm(j)) = tmp
enddo

end subroutine LDL_solve

subroutine exchange(n,A,perm,j1_IN,j2_IN)
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:)
integer,intent(inout) :: perm(:)
integer,intent(in) :: j1_IN,j2_IN
integer :: j1,j2
integer :: i
type(mp_real) :: rtmp
integer :: itmp

if(j1_IN/=j2_IN) then

   j1 = min(j1_IN,j2_IN)
   j2 = max(j1_IN,j2_IN)

   rtmp     = A(j1,j1)
   A(j1,j1) = A(j2,j2)
   A(j2,j2) = rtmp

   do i=1,j1-1
      rtmp    = A(i,j1)
      A(i,j1) = A(i,j2)
      A(i,j2) = rtmp
   enddo
   do i=j1+1,j2-1
      rtmp    = A(j1,i)
      A(j1,i) = A(i,j2)
      A(i,j2) = rtmp
   enddo
   do i=j2+1,n
      rtmp    = A(j1,i)
      A(j1,i) = A(j2,i)
      A(j2,i) = rtmp
   enddo

   itmp     = perm(j1)
   perm(j1) = perm(j2)
   perm(j2) = itmp

endif

end subroutine exchange

end module lproblem_mp
