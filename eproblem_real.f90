module eproblem_real
implicit none

private
public eprec
public symU_diagonalize_real
public test_diagonalize_real

integer,parameter :: eprec = 16

integer,parameter :: max_QL_iterSingle = 30
real(eprec),parameter :: thr_ZERO     = tiny(0._eprec)
real(eprec),parameter :: thr_EPSILON  = epsilon(0._eprec)
real(eprec),parameter :: thr_ROTATION = sqrt(1/epsilon(0._eprec))

contains

subroutine symU_diagonalize_real(n,diag,V,H,S)
implicit none
integer,intent(in) :: n
real(eprec),intent(out) :: diag(:)
real(eprec),intent(out) :: V(:,:)
real(eprec),intent(inout) :: H(:,:)
real(eprec),intent(inout),optional :: S(:,:)
logical :: generalized
real(eprec),allocatable :: super(:)
integer,allocatable :: perm(:)
integer :: i,j

generalized = present(S)

if(n<0) then

   write(*,'(a)') 'Negative matrix dimension in diagonalization!'
   stop

elseif(n==0) then

   write(*,'(a)') 'WARNING!!! Zero matrix dimension in diagonalization!'
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
      call LDL_factor(n,S,perm)
      call LDL_transform(n,H,S,perm,V)
   endif

   V(1:n,1:n) = 0
   do i=1,n
      V(i,i) = 1
   enddo

   call tridiagonal_upper(n,H,V)
   do i=1,n
      diag(i) = H(i,i)
   enddo
   do i=1,n-1
      super(i) = H(i,i+1)
   enddo
   super(n) = 0
   call implicitQL(n,diag,super,V)

   if(generalized) then
      call LDL_solve(n,V,S,perm,super)
   endif

   call CombSort(n,diag,perm)
   do j=1,n
      H(1:n,j) = V(1:n,perm(j))
   enddo
   V(1:n,1:n) = H(1:n,1:n)

   deallocate(super,perm)

endif

end subroutine symU_diagonalize_real

subroutine LDL_factor(n,A,perm)
implicit none
integer,intent(in) :: n
real(eprec),intent(inout) :: A(:,:)
integer,intent(out) :: perm(:)
integer :: i,j,k,k_max
real(eprec) :: rtmp
integer :: itmp

do i=1,n
   perm(i) = i
enddo

do k=n,2,-1

   rtmp  = abs(A(k,k))
   k_max = k
   do i=1,k-1
      if(abs(A(i,i))>rtmp) then
         rtmp  = abs(A(i,i))
         k_max = i
      endif
   enddo

   if(k_max<k) then
      rtmp           = A(k,k)
      A(k,k)         = A(k_max,k_max)
      A(k_max,k_max) = rtmp
      do i=1,k_max-1
         rtmp       = A(i,k)
         A(i,k)     = A(i,k_max)
         A(i,k_max) = rtmp
      enddo
      do i=k_max+1,k-1
         rtmp       = A(i,k)
         A(i,k)     = A(k_max,i)
         A(k_max,i) = rtmp
      enddo
      do i=k+1,n
         rtmp       = A(k,i)
         A(k,i)     = A(k_max,i)
         A(k_max,i) = rtmp
      enddo
      itmp        = perm(k)
      perm(k)     = perm(k_max)
      perm(k_max) = itmp
   endif

   A(k,k) = 1/A(k,k)

   do j=1,k-1
      rtmp   = A(j,k)
      A(j,k) = A(j,k)*A(k,k)
      do i=1,j
         A(i,j) = A(i,j) - A(i,k)*rtmp
      enddo
   enddo

enddo

A(1,1) = 1/A(1,1)

do k=1,n
   if(A(k,k)<0) then
      write(*,'(a)') 'S matrix in diagonalization is not positive definite!'
      stop
   endif
   A(k,k) = sqrt(A(k,k))
enddo

end subroutine LDL_factor

subroutine LDL_transform(n,A,U,perm,B)
implicit none
integer,intent(in) :: n
real(eprec),intent(inout) :: A(:,:)
real(eprec),intent(in) :: U(:,:)
integer,intent(in) :: perm(:)
real(eprec) :: B(:,:)
integer :: i,j,k
real(eprec) :: rtmp

do i=1,n
   k = perm(i)
   do j=1,k
      B(i,j) = A(j,k)
   enddo
   do j=k+1,n
      B(i,j) = A(k,j)
   enddo
enddo

do k=1,n
   do j=n,1,-1
      rtmp = B(j,k)
      do i=1,j-1
         B(i,k) = B(i,k) - U(i,j)*rtmp
      enddo
      B(j,k) = U(j,j)*rtmp
   enddo
enddo

do j=1,n
   k = perm(j)
   do i=1,j
      A(i,j) = B(i,k)
   enddo
enddo

do j=n,1,-1
   do i=1,j-1
      rtmp = U(i,j)
      do k=1,i
         A(k,i) = A(k,i) - A(k,j)*rtmp
      enddo
   enddo
   rtmp = U(j,j)
   do k=1,j
      A(k,j) = A(k,j)*rtmp
   enddo
enddo

end subroutine LDL_transform

subroutine LDL_solve(n,A,U,perm,b)
implicit none
integer,intent(in) :: n
real(eprec),intent(inout) :: A(:,:)
real(eprec),intent(in) :: U(:,:)
integer,intent(in) :: perm(:)
real(eprec) :: b(:)
integer :: i,j,k
real(eprec) :: rtmp

do k=1,n
   do j=1,n
      rtmp = A(j,k)*U(j,j)
      do i=1,j-1
         rtmp = rtmp - b(i)*U(i,j)
      enddo
      b(j) = rtmp
   enddo
   do j=1,n
      A(perm(j),k) = b(j)
   enddo
enddo

end subroutine LDL_solve

subroutine tridiagonal_upper(n,A,V)
implicit none
integer,intent(in) :: n
real(eprec),intent(inout) :: A(:,:),V(:,:)
integer :: i,j,k,j_1
real(eprec) :: sumoff,scale,inv_scale,xnorm,xnorm2,kfac
real(eprec) :: rtmp,rtmp1,rtmp2

do j=n,3,-1
   j_1 = j - 1

   sumoff = 0
   do i=1,j-2
      sumoff = sumoff + abs(A(i,j))
   enddo
   if(sumoff<thr_ZERO) cycle

   scale = sumoff + abs(A(j_1,j))
   inv_scale = 1/scale
   do i=1,j_1
      A(i,j) = A(i,j)*inv_scale
   enddo

   xnorm2 = 0
   do i=1,j_1
      xnorm2 = xnorm2 + A(i,j)**2
   enddo
   xnorm = sign(sqrt(xnorm2),A(j_1,j))
   rtmp = 1/sqrt(xnorm2 + xnorm*A(j_1,j))
   do i=1,j-2
      V(i,j) = A(i,j)*rtmp
   enddo
   V(j_1,j) = (A(j_1,j) + xnorm)*rtmp

   do i=1,j_1
      rtmp1 = 0
      rtmp2 = V(i,j)
      do k=1,i-1
         rtmp1  = rtmp1  + A(k,i)*V(k,j)
         A(k,j) = A(k,j) + A(k,i)*rtmp2
      enddo
      A(i,j) = rtmp1 + A(i,i)*rtmp2
   enddo

   rtmp = 0
   do i=1,j_1
      rtmp = rtmp + A(i,j)*V(i,j)
   enddo
   kfac = rtmp*0.5d0
   do i=1,j_1
      A(i,j) = A(i,j) - kfac*V(i,j)
   enddo

   do i=1,j_1
      rtmp1 = V(i,j)
      rtmp2 = A(i,j)
      do k=1,i
         A(k,i) = A(k,i) - rtmp1*A(k,j) - rtmp2*V(k,j)
      enddo
   enddo

   do i=1,j-2
      A(i,j) = 0
   enddo
   A(j_1,j) = -xnorm*scale

enddo

do j=3,n
   j_1 = j - 1

   if(abs(V(j_1,j))<thr_ZERO) cycle

   do i=1,j_1
      rtmp = 0
      do k=1,j_1
         rtmp = rtmp + V(k,i)*V(k,j)
      enddo
      do k=1,j_1
         V(k,i) = V(k,i) - rtmp*V(k,j)
      enddo
   enddo

   do i=1,j_1
      V(i,j) = 0
   enddo

enddo

end subroutine tridiagonal_upper

subroutine implicitQL(n,A,B,V)
implicit none
integer,intent(in) :: n
real(eprec),intent(inout) :: A(:),B(:),V(:,:)
integer :: jstart,jend,j,j_1,i
integer :: iter,iterSingle,Nrot
real(eprec) :: theta,t,c,s,shift,rot
real(eprec) :: d1,d2,e,f1,f2,pythag,l12

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
         j_1 = j - 1

         d2 = d1
         d1 = A(j_1)
         e  = c*B(j_1)
         f1 = s*B(j_1)

         pythag = hypot(f1,f2)
         c = f2/pythag
         s = f1/pythag

         l12 = s*(d1 - d2) + 2*c*e
         d1   = d1 - s*l12
         A(j) = d2 + s*l12
         B(j) = pythag
         f2   = c*l12 - e

         Nrot = Nrot + 1
         do i=1,n
            rot = V(i,j_1)
            V(i,j_1) = c*rot - s*V(i,j)
            V(i,j)   = s*rot + c*V(i,j)
         enddo
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
real(eprec),intent(inout) :: val(:)
integer,intent(out) :: perm(:)
integer :: i,gap,itmp
real(eprec) :: rtmp
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
      if (val(i+gap)<val(i)) then

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

subroutine test_diagonalize_real(n,diag,V,H,S)
implicit none
integer,intent(in) :: n
real(eprec),intent(in) :: diag(:)
real(eprec),intent(in) :: V(:,:)
real(eprec),intent(in) :: H(:,:)
real(eprec),intent(in),optional :: S(:,:)
logical :: generalized
real(eprec),allocatable :: HV(:,:),SV(:,:)
real(eprec) :: max_resid,mean_resid,max_dia,mean_dia,max_off,mean_off
real(eprec) :: lambda,rtmp
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
write(*,'(a,es20.10)') 'max:  ',max_resid
write(*,'(a,es20.10)') 'mean: ',mean_resid

write(*,'(a)') 'Norm of eigenvector:               ||i|| - 1 =   '
write(*,'(a,es20.10)') 'max:  ',max_dia
write(*,'(a,es20.10)') 'mean: ',mean_dia

write(*,'(a,2es11.3)') 'Product of different eigenvectors:    <i|j> =    '
write(*,'(a,es20.10)') 'max:  ',max_off
write(*,'(a,es20.10)') 'mean: ',mean_off

end subroutine test_diagonalize_real

subroutine multiply_upper(n,A,B,C)
implicit none
integer,intent(in) :: n
real(eprec),intent(in) :: A(:,:),B(:,:)
real(eprec),intent(out) :: C(:,:)
integer :: i,j,k
real(eprec) :: rtmp1,rtmp2

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

end module eproblem_real
