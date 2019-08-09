module eproblem
use mpmodule
implicit none

private
public symU_diagonalize
public test_diagonalize

integer,parameter :: max_QL_iterSingle = 30

contains

subroutine symU_diagonalize(eps,n,diag,V,H,S)
implicit none
type(mp_real),intent(in) :: eps
integer,intent(in) :: n
type(mp_real),intent(out) :: diag(:)
type(mp_real),intent(out) :: V(:,:)
type(mp_real),intent(inout) :: H(:,:)
type(mp_real),intent(inout),optional :: S(:,:)
logical :: generalized
type(mp_real),allocatable :: super(:)
integer,allocatable :: perm(:)
integer :: i,j

generalized = present(S)

if(n<0) then

   write(*,'(a)') 'Negative matrix dimension in diagonalization!'
   stop

elseif(n==0) then

   return

elseif(n==1) then

   if(generalized) then
      if(S(1,1)<0) then
         write(*,'(a)') 'S matrix in diagonalization is not positive definite!'
         stop
      endif
      diag(1) = H(1,1)/S(1,1)
      V(1,1) = 1/sqrt(S(1,1))
   else
      diag(1) = H(1,1)
      V(1,1) = 1
   endif

else

   allocate(super(n),perm(n))

   if(generalized) then
      call UDUT_factor(n,S,perm,diag)
      call UDUT_transform(n,H,S,perm,V)
   endif

!$OMP PARALLEL DO COLLAPSE(2)
   do j=1,n
      do i=1,n
         if(i==j) then
            V(i,j) = 1
         else
            V(i,j) = 0
         endif
      enddo
   enddo
!$OMP END PARALLEL DO

   call tridiagonal_upper(n,H,V)

!$OMP PARALLEL
!$OMP DO
   do i=1,n
      diag(i) = H(i,i)
   enddo
!$OMP END DO NOWAIT
!$OMP DO
   do i=1,n-1
      super(i) = H(i,i+1)
   enddo
!$OMP END DO NOWAIT
!$OMP SINGLE
   super(n) = 0
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL

   call implicitQL(eps,n,diag,super,V)

   if(generalized) then
      call UDUT_solve(n,V,S,perm,H)
   endif

   call CombSort(n,diag,perm)

!$OMP PARALLEL
!$OMP DO COLLAPSE(2)
   do j=1,n
      do i=1,n
         H(i,j) = V(i,perm(j))
      enddo
   enddo
!$OMP END DO
!$OMP DO COLLAPSE(2)
   do j=1,n
      do i=1,n
         V(i,j) = H(i,j)
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
   
   deallocate(super,perm)

endif

end subroutine symU_diagonalize

subroutine UDUT_factor(n,A,perm,work)
! A -> U.D.UT factorization, with U upper-unit-diagonal
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:)
integer,intent(out) :: perm(:)
type(mp_real) :: work(:)
integer :: i,j,k,k_max,itmp
type(mp_real) :: rtmp
logical :: positive

do i=1,n
   perm(i) = i
enddo

do k=n,1,-1

   rtmp = -1
   do i=1,k
      if(abs(A(i,i))>rtmp) then
         rtmp  = abs(A(i,i))
         k_max = i
      endif
   enddo

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

positive = .true.
!$OMP PARALLEL DO
do k=1,n
   if(A(k,k)<0) then
!$OMP ATOMIC WRITE
      positive = .false.
!$OMP END ATOMIC
   else
      A(k,k) = sqrt(A(k,k))
   endif
enddo
!$OMP END PARALLEL DO

if(.not.positive) then
   write(*,'(a)') 'S matrix in diagonalization is not positive definite!'
   stop
endif

end subroutine UDUT_factor

subroutine UDUT_transform(n,A,U,perm,B)
! transformation iDh.U^(-1).A.UT^(-1).iDh -> A
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:)
type(mp_real),intent(in) :: U(:,:)
integer,intent(in) :: perm(:)
type(mp_real) :: B(:,:)
integer :: i,j,k

!$OMP PARALLEL

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      if(perm(i)<=j) then
         B(i,j) = A(perm(i),j)
      else
         B(i,j) = A(j,perm(i))
      endif
   enddo
enddo
!$OMP END DO

!$OMP DO
do k=1,n
   associate(x => B(:,k))
     do j=n,1,-1
        do i=1,j-1
           x(i) = x(i) - U(i,j)*x(j)
        enddo
        x(j) = U(j,j)*x(j)
     enddo
   end associate
enddo
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC)
do j=1,n
   do i=1,j
      A(i,j) = B(i,perm(j))
   enddo
enddo
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC)
do k=1,n
   associate(x => A(k,:))
     do j=n,k,-1
        do i=k,j-1
           x(i) = x(i) - U(i,j)*x(j)
        enddo
        x(j) = U(j,j)*x(j)
     enddo
   end associate
enddo
!$OMP END DO

!$OMP END PARALLEL

end subroutine UDUT_transform

subroutine UDUT_solve(n,A,U,perm,B)
! solution of UT^(-1).iDh.A -> A
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:)
type(mp_real),intent(in) :: U(:,:)
integer,intent(in) :: perm(:)
type(mp_real) :: B(:,:)
integer :: i,j,k

!$OMP PARALLEL

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      B(i,j) = A(i,j)*U(i,i)
   enddo
enddo
!$OMP END DO

!$OMP DO
do k=1,n
   associate(x => B(:,k))
     do j=2,n
        do i=1,j-1
           x(j) = x(j) - U(i,j)*x(i)
        enddo
     enddo
   end associate
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      A(perm(i),j) = B(i,j)
   enddo
enddo
!$OMP END DO

!$OMP END PARALLEL

end subroutine UDUT_solve

subroutine tridiagonal_upper(n,A,V)
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:,:),V(:,:)
integer :: i,j,k
type(mp_real) :: scale,xnorm,factor
type(mp_real) :: rtmp

do k=n,3,-1

   rtmp = 0
   do i=1,k-2
      rtmp = rtmp + abs(A(i,k))
   enddo
   if(rtmp==0) cycle

   scale = rtmp + abs(A(k-1,k))

   factor = 1/scale
   rtmp = 0
   do i=1,k-1
      A(i,k) = factor*A(i,k)
      rtmp = rtmp + A(i,k)**2
   enddo
   
   xnorm = sign(sqrt(rtmp),A(k-1,k))

   A(k-1,k) = A(k-1,k) + xnorm
   factor = 1/sqrt(A(k-1,k)*xnorm)
   do j=1,k-1
      V(j,k) = factor*A(j,k)
      A(j,k) = A(j,j)*V(j,k)
      do i=1,j-1
         A(i,k) = A(i,k) + A(i,j)*V(j,k)
         A(j,k) = A(j,k) + A(i,j)*V(i,k)
      enddo
   enddo

   rtmp = 0
   do i=1,k-1
      rtmp = rtmp + A(i,k)*V(i,k)
   enddo
   factor = rtmp/2
   
!$OMP PARALLEL
!$OMP DO
   do i=1,k-1
      A(i,k) = A(i,k) - factor*V(i,k)
   enddo
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
   do j=1,k-1
      do i=1,j
         A(i,j) = A(i,j) - A(i,k)*V(j,k) - A(j,k)*V(i,k)
      enddo
   enddo
!$OMP END DO
!$OMP DO
   do i=1,k-1
      A(i,k) = 0
   enddo
!$OMP END DO
!$OMP END PARALLEL

   A(k-1,k) = -scale*xnorm

enddo

do k=3,n

   if(abs(V(k-1,k))==0) cycle

!$OMP PARALLEL PRIVATE(rtmp)
!$OMP DO SCHEDULE(DYNAMIC)
   do j=1,k-1
      rtmp = 0
      do i=1,k-1
         rtmp = rtmp + V(i,j)*V(i,k)
      enddo
      do i=1,k-1
         V(i,j) = V(i,j) - rtmp*V(i,k)
      enddo
   enddo
!$OMP END DO
!$OMP DO
   do i=1,k-1
      V(i,k) = 0
   enddo
!$OMP END DO
!$OMP END PARALLEL

enddo

end subroutine tridiagonal_upper

subroutine implicitQL(eps,n,A,B,V)
implicit none
type(mp_real),intent(in) :: eps
integer,intent(in) :: n
type(mp_real),intent(inout) :: A(:),B(:),V(:,:)
type(mp_real) :: thr_EPSILON
type(mp_real) :: thr_ROTATION
integer :: i,j,jstart,jend
integer :: iter,iterSingle,Nrot
type(mp_real) :: theta,t,c,s,cs,shift,rtmp
type(mp_real) :: d1,d2,e,f1,f2,l12

thr_EPSILON = eps
thr_ROTATION = sqrt(1/eps)

iter = 0
iterSingle = 0
Nrot = 0

jstart = 1
do while(jstart<n)

   if(iterSingle<max_QL_iterSingle) then
      jend = jstart
      do while(jend<n)
         if(abs(B(jend))<thr_EPSILON*(abs(A(jend))+abs(A(jend+1)))) exit
         jend = jend + 1
      enddo
   else
      write(*,'(a,i6)') &
           'WARNING!!! Too many iterations to converge eigenvalue: ',jstart
      jend = jstart
   endif

   if(jend==jstart) then

      jstart = jstart + 1
      iterSingle = 0

   else

      iter = iter + 1
      iterSingle = iterSingle + 1

      d1 = A(jstart)
      d2 = A(jstart+1)
      e  = B(jstart)
      theta = (d2 - d1)/(2*e)
      if(abs(theta)<thr_ROTATION) then
         t = 1/(theta + sign(sqrt(theta**2 + 1),theta))
      else
         t = 1/(2*theta)
      endif
      shift = d1 - t*e

      c = 1
      s = 1
      d1 = A(jend)
      f2 = d1 - shift
      do j=jend,jstart+1,-1

         d2 = d1
         d1 = A(j-1)
         e  = c*B(j-1)
         f1 = s*B(j-1)

         cs = hypot(f1,f2)
         c  = f2/cs
         s  = f1/cs

         l12  = s*(d1 - d2) + 2*c*e
         d1   = d1 - s*l12
         A(j) = d2 + s*l12
         B(j) = cs
         f2   = c*l12 - e

         Nrot = Nrot + 1
!$OMP PARALLEL DO PRIVATE(rtmp)
         do i=1,n
            rtmp     = c*V(i,j-1) - s*V(i,j)
            V(i,j)   = s*V(i,j-1) + c*V(i,j)
            V(i,j-1) = rtmp
         enddo
!$OMP END PARALLEL DO

      enddo
      A(jstart) = d1
      B(jstart) = f2
      B(jend) = 0

   endif

enddo

end subroutine implicitQL

subroutine CombSort(n,val,perm)
implicit none
integer,intent(in) :: n
type(mp_real),intent(inout) :: val(:)
integer,intent(out) :: perm(:)
integer :: i,gap,itmp
type(mp_real) :: rtmp
logical :: swapped

do i=1,n
   perm(i) = i
enddo

gap = n
swapped = .true.
do while(gap>1.or.swapped)
   gap = max(1,(gap*10)/13)
   if(gap==9.or.gap==10) gap = 11
   swapped = .false.

   do i=1,n-gap
      if(val(i+gap)<val(i)) then

         rtmp = val(i)
         val(i) = val(i+gap)
         val(i+gap) = rtmp

         itmp = perm(i)
         perm(i) = perm(i+gap)
         perm(i+gap) = itmp

         swapped = .true.
      endif
   enddo

enddo

end subroutine CombSort

subroutine test_diagonalize(n,diag,V,H,S)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: diag(:)
type(mp_real),intent(in) :: V(:,:)
type(mp_real),intent(in) :: H(:,:)
type(mp_real),intent(in),optional :: S(:,:)
logical :: generalized
type(mp_real),allocatable :: HV(:,:),SV(:,:)
type(mp_real) :: max_resid,mean_resid,max_dia,mean_dia,max_off,mean_off
type(mp_real) :: lambda,rtmp
integer :: i,j,k

generalized = present(S)

allocate(HV(n,n),SV(n,n))

call multiply_upper(n,H,V,HV)
if(generalized) then
   call multiply_upper(n,S,V,SV)
else
   SV(1:n,1:n) = V(1:n,1:n)
endif

max_resid  = 0
mean_resid = 0
do j=1,n
   lambda = diag(j)
   rtmp = 0
   do i=1,n
      rtmp = rtmp + (HV(i,j) - lambda*SV(i,j))**2
   enddo
   rtmp = sqrt(rtmp)
   max_resid  = max(max_resid,rtmp)
   mean_resid = mean_resid + rtmp
enddo
mean_resid = mean_resid/n

max_dia  = 0
mean_dia = 0
do j=1,n
   rtmp = 0
   do i=1,n
      rtmp = rtmp + V(i,j)*SV(i,j)
   enddo
   rtmp = -1 + sqrt(rtmp)
   max_dia  = max(max_dia,abs(rtmp))
   mean_dia = mean_dia + rtmp
enddo
mean_dia = mean_dia/n

max_off  = 0
mean_off = 0
do j=1,n
   do i=1,n
      if(i==j) cycle
      rtmp = 0
      do k=1,n
         rtmp = rtmp + V(k,i)*SV(k,j)
      enddo
      max_off  = max(max_off,abs(rtmp))
      mean_off = mean_off + rtmp
   enddo
enddo
mean_off = mean_off/max(n*(n-1),1)

deallocate(HV,SV)

if(generalized) then
   write(*,'(a)') 'Norm of residual vector:          ||HV - e SV|| ='
else
   write(*,'(a)') 'Norm of residual vector:          ||HV -  e V|| ='
endif
write(*,'(a)',advance='no') 'max:  '; call mpwrite(6,30,20,max_resid)
write(*,'(a)',advance='no') 'mean: '; call mpwrite(6,30,20,mean_resid)

write(*,'(a)') 'Norm of eigenvector:               ||i|| - 1 =   '
write(*,'(a)',advance='no') 'max:  '; call mpwrite(6,30,20,max_dia)
write(*,'(a)',advance='no') 'mean: '; call mpwrite(6,30,20,mean_dia)

write(*,'(a,2es11.3)') 'Product of different eigenvectors:    <i|j> =    '
write(*,'(a)',advance='no') 'max:  '; call mpwrite(6,30,20,max_off)
write(*,'(a)',advance='no') 'mean: '; call mpwrite(6,30,20,mean_off)

end subroutine test_diagonalize

subroutine multiply_upper(n,A,B,C)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: A(:,:),B(:,:)
type(mp_real),intent(out) :: C(:,:)
integer :: i,j,k
type(mp_real) :: rtmp1,rtmp2

do j=1,n
   do i=1,n
      rtmp1 = 0
      rtmp2 = B(i,j)
      do k=1,i-1
         rtmp1  = rtmp1  + A(k,i)*B(k,j)
         C(k,j) = C(k,j) + A(k,i)*rtmp2
      enddo
      C(i,j) = rtmp1 + A(i,i)*rtmp2
   enddo
enddo

end subroutine multiply_upper

end module eproblem
