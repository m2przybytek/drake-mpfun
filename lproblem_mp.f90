module lproblem_mp
use mpmodule
implicit none

private
public symU_linearsolve_mp

contains

subroutine symU_linearsolve_mp(n,A,x)
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:),x(:)

if(n<0) then

   write(*,'(a)') 'Negative matrix dimension in linear equations solver!'
   stop

elseif(n==0) then

   write(*,'(a)') 'WARNING!!! Zero matrix dimension in linear equations solver!'
   return

elseif(n==1) then

   x(1) = x(1)/A(1,1)

else

   call cholesky_factor(n,A)
   call cholesky_solve(n,A,x)

endif

end subroutine symU_linearsolve_mp

subroutine cholesky_factor(n,A)
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:)
integer :: i,j,k
type(mp_real) :: rtmp

do j=1,n
   do i=1,j-1
      rtmp = A(i,j)
      do k=1,i-1
         rtmp = rtmp - A(k,i)*A(k,j)
      enddo
      A(i,j) = rtmp*A(i,i)
   enddo
   rtmp = A(j,j)
   do k=1,j-1
      rtmp = rtmp - A(k,j)**2
   enddo
   if(rtmp<0) then
      write(*,'(a)') &
           'Matrix in linear equations solver is not positive definite!'
      stop
   endif
   A(j,j) = 1/sqrt(rtmp)
enddo

end subroutine cholesky_factor

subroutine cholesky_solve(n,A,x)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: A(:,:)
type(mp_real),intent(inout) :: x(:)
integer :: i,j
type(mp_real) :: rtmp

do j=1,n
   rtmp = x(j)
   do i=1,j-1
      rtmp = rtmp - A(i,j)*x(i)
   enddo
   x(j) = rtmp*A(j,j)
enddo

do j=n,1,-1
   rtmp = x(j)*A(j,j)
   do i=1,j-1
      x(i) = x(i) - rtmp*A(i,j)
   enddo
   x(j) = rtmp
enddo

end subroutine cholesky_solve

end module lproblem_mp
