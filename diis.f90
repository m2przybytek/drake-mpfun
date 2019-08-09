module diis
use mpmodule
implicit none

private
public DIISData

type DIISData
type(mp_real) :: eps
integer :: nD1,nD2
integer :: max_n,act_n,act_i
type(mp_real),allocatable :: A(:,:,:),E(:,:,:)
type(mp_real),allocatable :: mat(:,:)
type(mp_real),allocatable :: sol_mat(:,:),sol_vec(:)
integer,allocatable :: sol_perm(:)
contains
procedure :: init => init_DIIS
procedure :: free => free_DIIS
procedure :: step => step_DIIS
end type DIISData

contains

subroutine init_DIIS(DIIS,eps,nD1,nD2,max_n)
implicit none
class(DIISData) :: DIIS
type(mp_real),intent(in) :: eps
integer,intent(in) :: nD1,nD2
integer,intent(in) :: max_n

DIIS%eps   = eps
DIIS%nD1   = nD1
DIIS%nD2   = nD2
DIIS%max_n = max_n
DIIS%act_n = 0
DIIS%act_i = 0

allocate(&
     DIIS%A(DIIS%nD1,DIIS%nD2,DIIS%max_n),&
     DIIS%E(DIIS%nD1,DIIS%nD2,DIIS%max_n),&
     DIIS%mat(DIIS%max_n,DIIS%max_n))

allocate(&
     DIIS%sol_mat(DIIS%max_n+1,DIIS%max_n+1),&
     DIIS%sol_vec(DIIS%max_n+1),&
     DIIS%sol_perm(DIIS%max_n+1))

end subroutine init_DIIS

subroutine free_DIIS(DIIS)
implicit none
class(DIISData) :: DIIS

deallocate(DIIS%sol_mat,DIIS%sol_vec,DIIS%sol_perm)
deallocate(DIIS%A,DIIS%E,DIIS%mat)

end subroutine free_DIIS

subroutine step_DIIS(DIIS,A,E)
implicit none
class(DIISData) :: DIIS
type(mp_real),intent(inout) :: A(:,:)
type(mp_real),intent(in) :: E(:,:)
integer :: i,j,k
type(mp_real) :: val,part

if(DIIS%act_n<DIIS%max_n) DIIS%act_n = DIIS%act_n + 1

DIIS%act_i = DIIS%act_i + 1
if(DIIS%act_i>DIIS%max_n) DIIS%act_i = 1

!$OMP PARALLEL
!$OMP DO COLLAPSE(2)
do j=1,DIIS%nD2
   do i=1,DIIS%nD1
      DIIS%A(i,j,DIIS%act_i) = A(i,j)
   enddo
enddo
!$OMP END DO NOWAIT
!$OMP DO COLLAPSE(2)
do j=1,DIIS%nD2
   do i=1,DIIS%nD1
      DIIS%E(i,j,DIIS%act_i) = E(i,j)
   enddo
enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL

do k=1,DIIS%act_n
   val = 0
!$OMP PARALLEL PRIVATE(part)
   part = 0
!$OMP DO COLLAPSE(2)
   do j=1,DIIS%nD2
      do i=1,DIIS%nD1
         part = part + DIIS%E(i,j,k)*DIIS%E(i,j,DIIS%act_i)
      enddo
   enddo
!$OMP END DO NOWAIT
!$OMP CRITICAL(use_diis_reduction)
   val = val + part
!$OMP END CRITICAL(use_diis_reduction)
!$OMP END PARALLEL
   DIIS%mat(k,DIIS%act_i) = val
   if(k/=DIIS%act_i) DIIS%mat(DIIS%act_i,k) = val
enddo

do j=1,DIIS%act_n
   do i=1,DIIS%act_n
      DIIS%sol_mat(i,j) = DIIS%mat(i,j)
   enddo
   DIIS%sol_mat(DIIS%act_n+1,j) = -1
enddo
do i=1,DIIS%act_n
   DIIS%sol_mat(i,DIIS%act_n+1) = -1
enddo
DIIS%sol_mat(DIIS%act_n+1,DIIS%act_n+1) = 0

do i=1,DIIS%act_n
   DIIS%sol_vec(i) = 0
enddo
DIIS%sol_vec(DIIS%act_n+1) = -1

call gauss_elimination_symmetric(DIIS%eps,&
     DIIS%act_n+1,DIIS%sol_mat,DIIS%sol_vec,DIIS%sol_perm)

do k=1,DIIS%act_n
   val = DIIS%sol_vec(k)
   if(k==1) then
!$OMP PARALLEL DO COLLAPSE(2)
      do j=1,DIIS%nD2
         do i=1,DIIS%nD1
            A(i,j) = val*DIIS%A(i,j,k)
         enddo
      enddo
!$OMP END PARALLEL DO
   else
!$OMP PARALLEL DO
      do j=1,DIIS%nD2
         do i=1,DIIS%nD1
            A(i,j) = A(i,j) + val*DIIS%A(i,j,k)
         enddo
      enddo
!$OMP END PARALLEL DO
   endif
enddo

end subroutine step_DIIS

subroutine gauss_elimination_symmetric(eps,n,A,X,perm)
implicit none
type(mp_real),intent(in) :: eps
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:),X(:)
integer,intent(inout) :: perm(:)
integer :: i,j,k,m,i_max,j_max,itmp
type(mp_real) :: rtmp

do i=1,n
   perm(i) = i
enddo

do k=1,n

   rtmp = -1
   do j=k,n
      do i=k,n
         if(abs(A(i,j))>rtmp) then
            rtmp  = abs(A(i,j))
            i_max = i
            j_max = j
         endif
      enddo
   enddo

   if(i_max/=k) then
      do j=1,n
         rtmp       = A(k,j)
         A(k,j)     = A(i_max,j)
         A(i_max,j) = rtmp
      enddo
      itmp        = perm(k)
      perm(k)     = perm(i_max)
      perm(i_max) = itmp
   endif

   if(j_max/=k) then
      do i=k,n
         rtmp       = A(i,k)
         A(i,k)     = A(i,j_max)
         A(i,j_max) = rtmp
      enddo
      rtmp     = X(k)
      X(k)     = X(j_max)
      X(j_max) = rtmp
   endif

   if(abs(A(k,k))<eps) exit

   A(k,k) = 1/A(k,k)
   do j=k+1,n
      rtmp = A(k,j)*A(k,k)
      A(k,j) = 0
      do i=k+1,n
         A(i,j) = A(i,j) - A(i,k)*rtmp
      enddo
      X(j) = X(j) - X(k)*rtmp
   enddo

enddo
m = k - 1

do j=n,m+1,-1
   A(j,n) = 0
enddo
A(m,n) = X(m)*A(m,m)
do j=m-1,1,-1
   rtmp = X(j)
   do i=j+1,m
      rtmp = rtmp - A(i,j)*A(i,n)
   enddo
   A(j,n) = rtmp*A(j,j)
enddo

do i=1,n
   X(perm(i)) = A(i,n)
enddo

end subroutine gauss_elimination_symmetric

end module diis
