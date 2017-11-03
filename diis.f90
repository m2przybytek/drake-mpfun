module diis
use mpmodule
implicit none

private
public DIISData
public init_DIIS,free_DIIS,use_DIIS

type DIISData
integer :: nD
integer :: max_n,act_n,act_i
type(mp_real),allocatable :: A(:,:),E(:,:)
type(mp_real),allocatable :: mat(:,:)
type(mp_real),allocatable :: sol_mat(:,:),sol_vec(:)
end type DIISData

contains

subroutine init_DIIS(DIIS,nD,max_n)
implicit none
type(DIISData) :: DIIS
integer,intent(in) :: nD
integer,intent(in) :: max_n

DIIS%nD    = nD
DIIS%max_n = max_n
DIIS%act_n = 0
DIIS%act_i = 0

allocate(&
     DIIS%A(DIIS%nD**2,DIIS%max_n),&
     DIIS%E(DIIS%nD**2,DIIS%max_n),&
     DIIS%mat(DIIS%max_n,DIIS%max_n))

allocate(DIIS%sol_mat(DIIS%max_n+1,DIIS%max_n+1),DIIS%sol_vec(DIIS%max_n+1))

end subroutine init_DIIS

subroutine free_DIIS(DIIS)
implicit none
type(DIISData) :: DIIS

deallocate(DIIS%sol_mat,DIIS%sol_vec)
deallocate(DIIS%A,DIIS%E,DIIS%mat)

end subroutine free_DIIS

subroutine use_DIIS(DIIS,A,E)
implicit none
type(DIISData) :: DIIS
type(mp_real),intent(inout) :: A(:,:)
type(mp_real),intent(in) :: E(:,:)
integer :: i,j,ij,k
type(mp_real) :: val

if(DIIS%act_n<DIIS%max_n) DIIS%act_n = DIIS%act_n + 1

DIIS%act_i = DIIS%act_i + 1
if(DIIS%act_i>DIIS%max_n) DIIS%act_i = 1

ij = 0
do j=1,DIIS%nD
   do i=1,DIIS%nD
      ij = ij + 1
      DIIS%A(ij,DIIS%act_i) = A(i,j)
   enddo
enddo
ij = 0
do j=1,DIIS%nD
   do i=1,DIIS%nD
      ij = ij + 1
      DIIS%E(ij,DIIS%act_i) = E(i,j)
   enddo
enddo

do i=1,DIIS%act_n
   val = 0
   do ij=1,DIIS%nD**2
      val = val + DIIS%E(ij,i)*DIIS%E(ij,DIIS%act_i)
   enddo
   DIIS%mat(i,DIIS%act_i) = val
   if(i/=DIIS%act_i) DIIS%mat(DIIS%act_i,i) = val
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

call gauss_elimination_symmetric(DIIS%act_n+1,DIIS%sol_mat,DIIS%sol_vec)

A(1:DIIS%nD,1:DIIS%nD) = mpreal(0.d0)
do k=1,DIIS%act_n
   val = DIIS%sol_vec(k)
   ij = 0
   do j=1,DIIS%nD
      do i=1,DIIS%nD
         ij = ij + 1
         A(i,j) = A(i,j) + val*DIIS%A(ij,k)
      enddo
   enddo
enddo

end subroutine use_DIIS

subroutine gauss_elimination_symmetric(n,A,X)
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:),X(:)
integer :: i,j,k
type(mp_real) :: factor,rtmp

if(n==1) then

   X(1) = X(1)/A(1,1)

else

   do i=1,n-1
      factor = 1/A(i,i)
      do j=i+1,n
         rtmp = A(i,j)*factor
         do k=i,n
            A(k,j) = A(k,j) - A(k,i)*rtmp
         enddo
         X(j) = X(j) - X(i)*rtmp
      enddo
      A(i,i) = factor
   enddo

   X(n) = X(n)/A(n,n)

   do j=n-1,1,-1
      do k=j+1,n
         X(j) = X(j) - A(k,j)*X(k)
      enddo
      X(j) = X(j)*A(j,j)
   enddo

endif

end subroutine gauss_elimination_symmetric

end module diis
