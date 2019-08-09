module lproblem
use mpmodule
implicit none

private
public symU_linearsolve,DecompositionData

type DecompositionData
integer :: n
type(mp_real),allocatable :: A(:,:),work(:)
integer,allocatable :: perm(:)
contains
procedure :: init => init_DecompositionData
procedure :: free => free_DecompositionData
procedure :: solve => use_DecompositionData
end type DecompositionData

contains

subroutine symU_linearsolve(n,A,b)
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:),b(:)
type(mp_real),allocatable :: work(:)
integer,allocatable :: perm(:)

if(n<0) then

   write(*,'(a)') 'Negative matrix dimension in linear equations solver!'
   stop

elseif(n==0) then

   return

elseif(n==1) then

   b(1) = b(1)/A(1,1)

else

   allocate(work(n))
   allocate(perm(n))

   call UDUT_factor(n,A,perm,work)
   call UDUT_solve(n,A,perm,b,work)

   deallocate(perm)
   deallocate(work)

endif

end subroutine symU_linearsolve

subroutine init_DecompositionData(LHS,n,A)
implicit none
class(DecompositionData) :: LHS
integer,intent(in) :: n
type(mp_real),intent(in) :: A(:,:)
integer :: i,j

if(n<0) then

   write(*,'(a)') 'Negative matrix dimension in linear equations solver!'
   stop

elseif(n==0) then

   LHS%n = 0

   return

elseif(n==1) then

   LHS%n = 1

   allocate(LHS%A(1,1))

   LHS%A(1,1) = 1/A(1,1)

else

   LHS%n = n

   allocate(LHS%A(LHS%n,LHS%n))
   allocate(LHS%work(LHS%n))
   allocate(LHS%perm(LHS%n))

!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
   do j=1,LHS%n
      do i=1,j
         LHS%A(i,j) = A(i,j)
      enddo
   enddo
!$OMP END PARALLEL DO

   call UDUT_factor(LHS%n,LHS%A,LHS%perm,LHS%work)

endif

end subroutine init_DecompositionData

subroutine free_DecompositionData(LHS)
implicit none
class(DecompositionData) :: LHS

if(LHS%n>0) then
   if(LHS%n>1) then
      deallocate(LHS%perm)
      deallocate(LHS%work)
   endif
   deallocate(LHS%A)
endif

end subroutine free_DecompositionData

subroutine use_DecompositionData(LHS,n,x)
implicit none
class(DecompositionData) :: LHS
integer,intent(in) :: n
type(mp_real),intent(inout) :: x(:)

if(LHS%n/=n) then

   write(*,'(a)') 'RHS vector has improper size in linear equations solver!'
   stop

endif

if(LHS%n==1) then

   x(1) = x(1)*LHS%A(1,1)

elseif(LHS%n>1) then

   call UDUT_solve(LHS%n,LHS%A,LHS%perm,x,LHS%work)

endif

end subroutine use_DecompositionData

subroutine UDUT_factor(n,A,perm,work)
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:)
integer,intent(out) :: perm(:)
type(mp_real) :: work(:)
integer :: i,j,k,k_max,itmp
type(mp_real) :: rtmp
logical :: do_warn

do i=1,n
   perm(i) = i
enddo

do_warn = .true.
do k=n,1,-1

   rtmp = -1
   do i=1,k
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

!$OMP PARALLEL PRIVATE(rtmp)

   if(k_max/=k) then
!$OMP SECTIONS
!$OMP SECTION
      itmp        = perm(k)
      perm(k)     = perm(k_max)
      perm(k_max) = itmp
!$OMP SECTION
      rtmp           = A(k,k)
      A(k,k)         = A(k_max,k_max)
      A(k_max,k_max) = rtmp
!$OMP END SECTIONS NOWAIT
!$OMP DO
      do i=1,n
         if(i<k_max) then
            rtmp       = A(i,k)
            A(i,k)     = A(i,k_max)
            A(i,k_max) = rtmp
         elseif(i>k) then
            rtmp       = A(k,i)
            A(k,i)     = A(k_max,i)
            A(k_max,i) = rtmp
         elseif(i/=k_max.and.i/=k) then
            rtmp       = A(i,k)
            A(i,k)     = A(k_max,i)
            A(k_max,i) = rtmp
         endif
      enddo
!$OMP END DO
   endif

!$OMP SINGLE
   A(k,k) = 1/A(k,k)
!$OMP END SINGLE

   if(k/=1) then
!$OMP DO
      do i=1,k-1
         work(i) = A(i,k)
         A(i,k)  = A(i,k)*A(k,k)
      enddo
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
      do j=1,k-1
         do i=1,j
            A(i,j) = A(i,j) - A(i,k)*work(j)
         enddo
      enddo
!$OMP END DO
   endif

!$OMP END PARALLEL

enddo

!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
do j=1,n-1
   do i=j+1,n
      A(i,j) = A(j,i)
   enddo
enddo
!$OMP END PARALLEL DO

end subroutine UDUT_factor

subroutine UDUT_solve(n,A,perm,x,work)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: A(:,:)
integer,intent(in) :: perm(:)
type(mp_real),intent(inout) :: x(:)
type(mp_real) :: work(:)
integer :: i,j

!$OMP PARALLEL

!$OMP DO
do i=1,n
   work(i) = x(perm(i))
enddo
!$OMP END DO

do j=n,2,-1
!$OMP DO
   do i=1,j-1
      work(i) = work(i) - A(i,j)*work(j)
   enddo
!$OMP END DO
enddo

!$OMP DO
do i=1,n
   work(i) = A(i,i)*work(i)
enddo
!$OMP END DO

do j=1,n-1
!$OMP DO
   do i=j+1,n
      work(i) = work(i) - A(i,j)*work(j)
   enddo
!$OMP END DO
enddo

!$OMP DO
do i=1,n
   x(perm(i)) = work(i)
enddo
!$OMP END DO

!$OMP END PARALLEL

end subroutine UDUT_solve

end module lproblem
