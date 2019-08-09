module integrals
use iso_fortran_env, only : real64
!$ use omp_lib
use mpmodule
use time
use global
use misc
use integrals_core
implicit none

private
public integrals_scfSH
public scfJauxData,init_scfJaux,free_scfJaux,integrals_scfJ
public integrals_X,integrals_Hrest,integrals_J,integrals_P
public integrals_L00,integrals_L10,integrals_L11
public integrals_vecX,integrals_vecP

type scfJauxData
integer,allocatable :: idx(:,:)
type(prod2HData) :: prod2H
type(mp_real),allocatable :: ints(:,:)
type(mp_real),allocatable :: sum_0(:,:),sum_1(:,:)
integer :: n_0,n_1
end type scfJauxData

type contData
integer :: n_0,n_1
integer :: par,n
type(mp_real),allocatable :: val(:)
end type contData

contains

subroutine integrals_scfSH(scfS,scfH,Control)
implicit none
type(mp_real),intent(out) :: scfS(0:,0:)
type(mp_real),intent(out) :: scfH(0:,0:)
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu,Twall
integer :: nbas,i,j
type(mp_real),allocatable :: norm(:)
type(mp_real) :: val

if(Control%IPRINT>=1) call timer('START',Tcpu,Twall)

nbas = Control%nbas

allocate(norm(0:nbas-1))
call norm_definition(nbas-1,norm)

val = sqrt(mppi())
norm(0) = val*norm(0)**2
do i=1,nbas-1
   val = val*(2*i)
   norm(i) = val*norm(i)**2
enddo

val = 1; val = val/2

!$OMP PARALLEL

!$OMP DO COLLAPSE(2)
do j=0,nbas-1
   do i=0,nbas-1
      if(i==j) then
         scfS(i,j) = norm(i)
      else
         scfS(i,j) = 0
      endif
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do j=0,nbas-1
   do i=0,nbas-1
      if(i==j) then
         scfH(i,j) = (i+val)*norm(i)
      else
         scfH(i,j) = 0
      endif
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

deallocate(norm)

if(Control%IPRINT>=1) call timer('scfSH',Tcpu,Twall)

end subroutine integrals_scfSH

subroutine init_scfJaux(scfJaux,Control)
implicit none
type(scfJauxData) :: scfJaux
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu,Twall
integer :: nbas,nidx_bas,n1,n2,i,j,ij,v
type(int2Data) :: int2

if(Control%IPRINT>=1) call timer('START',Tcpu,Twall)

nbas = Control%nbas
nidx_bas = (nbas*(nbas+1))/2

allocate(scfJaux%idx(2,nidx_bas))
call fill_idx(scfJaux%idx,nbas-1)

n1 = nbas-1
n2 = nbas-1
call prod2H_create(scfJaux%prod2H,n1,n2)
call prod2H_norm(scfJaux%prod2H)

n1 = nbas
n2 = nidx_bas
allocate(scfJaux%ints(n1,n2))

n1 = nbas
n2 = 1
!$ n2 = omp_get_max_threads()
allocate(scfJaux%sum_0(n1,n2),scfJaux%sum_1(n1,n2))

n1 = 2*(nbas-1)
n2 = 2*(nbas-1)
call int2_create('1 ','1 ','22',int2,n1,n2,&
     Control%alpha,Control%nC,Control%C)

scfJaux%n_0 = int2%n1_0
scfJaux%n_1 = int2%n1_1

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v)
do ij=1,nidx_bas
   i = scfJaux%idx(1,ij)
   j = scfJaux%idx(2,ij)
   v = i + j
   select case(mod(v,2))
   case(0)
      v = v/2+1
      call product2_unit(int2%n1_0,int2%val_0,v,scfJaux%prod2H%unit(i,j),&
           scfJaux%ints(:,ij))
   case(1)
      v = (v+1)/2
      call product2_unit(int2%n1_1,int2%val_1,v,scfJaux%prod2H%unit(i,j),&
           scfJaux%ints(:,ij))
   end select
enddo
!$OMP END PARALLEL DO

call int2_free(int2)

if(Control%IPRINT>=1) call timer('prepare_scfJ',Tcpu,Twall)

end subroutine init_scfJaux

subroutine free_scfJaux(scfJaux)
implicit none
type(scfJauxData) :: scfJaux

deallocate(scfJaux%sum_0,scfJaux%sum_1)
deallocate(scfJaux%ints)
call prod2H_free(scfJaux%prod2H)

deallocate(scfJaux%idx)

end subroutine free_scfJaux

subroutine integrals_scfJ(scfJ,scfD,scfJaux,Control)
implicit none
type(mp_real),intent(out) :: scfJ(0:,0:)
type(mp_real),intent(in) :: scfD(:,:)
type(scfJauxData) :: scfJaux
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu,Twall
integer :: nbas,nidx_bas
integer :: ij,i,j,v
type(mp_real) :: mlt,val

if(Control%IPRINT>=2) call timer('START',Tcpu,Twall)

nbas = Control%nbas
nidx_bas = (nbas*(nbas+1))/2

mlt = Control%gfac

call sums_D(scfJaux%n_0,scfJaux%n_1,scfJaux%sum_0,scfJaux%sum_1,&
     scfD,scfJaux%ints,&
     nidx_bas,scfJaux%idx)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v) &
!$OMP PRIVATE(val)
do ij=1,nidx_bas
   i = scfJaux%idx(1,ij)
   j = scfJaux%idx(2,ij)
   v = i + j
   select case(mod(v,2))
   case(0)
      v = v/2+1
      val = mlt*product1_unit(v,scfJaux%prod2H%unit(i,j),scfJaux%sum_0(:,1))
      scfJ(i,j) = val
      if(i/=j) scfJ(j,i) = val
   case(1)
      v = (v+1)/2
      val = mlt*product1_unit(v,scfJaux%prod2H%unit(i,j),scfJaux%sum_1(:,1))
      scfJ(i,j) = val
      scfJ(j,i) = val
   end select
enddo
!$OMP END PARALLEL DO

if(Control%IPRINT>=2) call timer('use_scfJ',Tcpu,Twall)

end subroutine integrals_scfJ

subroutine integrals_X(typeX,matX_S,matX_T,Control)
implicit none
character(1),intent(in) :: typeX
type(mp_real),intent(out) :: matX_S(:,:),matX_T(:,:)
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu,Twall
integer :: omega,nidx_omg,npair,n1,n2
integer :: i,j,k,l,ij,kl,v,v1,v2
type(mp_real) :: val_A,val_B
integer,allocatable :: idx(:,:),pair(:,:)
type(prod2HData) :: prod2H
type(int2Data) :: int2
type(mp_real),allocatable :: ints(:,:)

if(Control%IPRINT>=1) call timer('START',Tcpu,Twall)

omega = Control%omega
nidx_omg = ((omega+1)*(omega+2))/2
npair = Control%npair

allocate(idx(2,nidx_omg),pair(2,npair))
call fill_idx(idx,omega)
call fill_pair(pair,omega)

n1 = omega
n2 = omega
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = 2*omega
n2 = 2*omega
select case(typeX)
case('S')
   call int2_create('f ','f ','22',int2,n1,n2,&
        Control%alpha,Control%nC,Control%C)
case('H')
   call int2_create('df','df','22',int2,n1,n2,&
        Control%alpha,Control%nC,Control%C)
end select

n1 = omega+1
n2 = nidx_omg
allocate(ints(n1,n2))

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v)
do ij=1,nidx_omg
   i = idx(1,ij)
   j = idx(2,ij)
   v = i + j
   select case(mod(v,2))
   case(0)
      v = v/2+1
      call product2_unit(int2%n1_0,int2%val_0,v,prod2H%unit(i,j),&
           ints(:,ij))
   case(1)
      v = (v+1)/2
      call product2_unit(int2%n1_1,int2%val_1,v,prod2H%unit(i,j),&
           ints(:,ij))
   end select
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(i,j,k,l,v1,v2) &
!$OMP PRIVATE(val_A,val_B)
do kl=1,npair
   do ij=1,npair

      i = pair(1,ij)
      j = pair(2,ij)
      k = pair(1,kl)
      l = pair(2,kl)

      if(mod(i+j+k+l,2)==0) then

         v1 = i + k
         v1 = (v1+2-mod(v1,2))/2
         v2 = max(j,l)
         v2 = 1 + min(j,l) + (v2*(v2+1))/2
         val_A = product1_unit(v1,prod2H%unit(i,k),ints(:,v2))

         v1 = i + l
         v1 = (v1+2-mod(v1,2))/2
         v2 = max(j,k)
         v2 = 1 + min(j,k) + (v2*(v2+1))/2
         val_B = product1_unit(v1,prod2H%unit(i,l),ints(:,v2))

         matX_S(ij,kl) = val_A + val_B
         matX_T(ij,kl) = val_A - val_B

      else

         matX_S(ij,kl) = 0
         matX_T(ij,kl) = 0

      endif

   enddo
enddo
!$OMP END PARALLEL DO

deallocate(ints)

call int2_free(int2)

call prod2H_free(prod2H)

deallocate(idx,pair)

if(Control%IPRINT>=1) call timer(typeX,Tcpu,Twall)

end subroutine integrals_X

subroutine integrals_Hrest(matH_S,matH_T,matS_S,matS_T,Control)
implicit none
type(mp_real),intent(inout) :: matH_S(:,:),matH_T(:,:)
type(mp_real),intent(in) :: matS_S(:,:),matS_T(:,:)
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu,Twall
integer :: omega,npair,i,j,k,l,ij,kl,ijkl
integer,allocatable :: pair(:,:)

if(Control%IPRINT>=1) call timer('START',Tcpu,Twall)

omega = Control%omega
npair = Control%npair

allocate(pair(2,npair))
call fill_pair(pair,omega)

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k,l,ijkl)
do kl=1,npair
   do ij=1,npair

      i = pair(1,ij)
      j = pair(2,ij)
      k = pair(1,kl)
      l = pair(2,kl)
      ijkl = i + j + k + l

      if(mod(ijkl,2)==0) then

         matH_S(ij,kl) = matH_S(ij,kl) + matS_S(ij,kl)*(ijkl/2+1)
         matH_T(ij,kl) = matH_T(ij,kl) + matS_T(ij,kl)*(ijkl/2+1)

      endif

   enddo
enddo
!$OMP END PARALLEL DO

deallocate(pair)

if(Control%IPRINT>=1) call timer('Hrest',Tcpu,Twall)

end subroutine integrals_Hrest

subroutine integrals_J(matJ_S,matJ_T,matM_S,matM_T,matC,parC,Control)
implicit none
type(mp_real),intent(out) :: matJ_S(:,:,:,:),matJ_T(:,:,:,:)
type(mp_real),intent(out) :: matM_S(:,:,:,:),matM_T(:,:,:,:)
type(mp_real),intent(in) :: matC(:,:)
integer,intent(in) :: parC(:)
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu,Twall
integer :: nocc,nbas,omega,nidx_bas,nidx_omg,npair,n1,n2,nC
integer :: iocc,jocc,i,j,k,l,ij,kl,v,v1,v2
type(mp_real) :: mlt,val_A1,val_A2,val_B1,val_B2,val_A,val_B
integer,allocatable :: idx(:,:),pair(:,:)
type(prod2HData) :: prod2H
type(int2Data) :: int2
type(auxData) :: auxC1
type(mp_real),allocatable :: ints(:,:)
type(mp_real),allocatable :: sum_0(:,:),sum_1(:,:)
type(contData),allocatable :: conts(:)

if(Control%IPRINT>=1) call timer('START',Tcpu,Twall)

nocc  = Control%nocc
nbas  = Control%nbas
omega = Control%omega
nidx_bas = (nbas*(nbas+1))/2
nidx_omg = ((omega+1)*(omega+2))/2
npair = Control%npair

allocate(idx(2,max(nidx_omg,nidx_bas)),pair(2,npair))
call fill_idx(idx,max(omega,nbas-1))
call fill_pair(pair,omega)

n1 = max(omega,nbas-1)
n2 = max(omega,nbas-1)
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = omega + nbas
n2 = nidx_omg
allocate(ints(n1,n2))

n1 = 2*(omega + nbas-1)
n2 = 2* omega
call int2_create('f ','f ','42',int2,n1,n2,&
     Control%alpha,Control%nC,Control%C)
mlt = Control%gfac

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v)
do ij=1,nidx_omg
   i = idx(1,ij)
   j = idx(2,ij)
   v = i + j
   select case(mod(v,2))
   case(0)
      v = v/2+1
      call product2_unit(int2%n1_0,int2%val_0,v,prod2H%unit(i,j),&
           ints(:,ij))
   case(1)
      v = (v+1)/2
      call product2_unit(int2%n1_1,int2%val_1,v,prod2H%unit(i,j),&
           ints(:,ij))
   end select
enddo
!$OMP END PARALLEL DO

call int2_free(int2)

n1 = 2*(omega + nbas-1)
n2 = 2* omega
nC = 2*(nbas-1)
call aux_init(auxC1,n1,n2,nC)

n1 = nbas
n2 = 1
!$ n2 = omp_get_max_threads()
allocate(sum_0(n1,n2),sum_1(n1,n2))

nC = nidx_omg
allocate(conts(nC))
do ij=1,nidx_omg
   i = idx(1,ij)
   j = idx(2,ij)
   associate(cont => conts(ij))
     nC = i + j + 2*(nbas-1)
     cont%n_0 = nC/2+1
     cont%n_1 = (nC+1)/2
     allocate(cont%val(max(cont%n_0,cont%n_1)))
   end associate
enddo

do jocc=1,nocc
   do iocc=1,jocc

      call sums_CC(nbas,nbas,sum_0,sum_1,&
           matC(:,iocc),matC(:,jocc),prod2H,&
           nidx_bas,idx)
      call auxC1_create(auxC1,sum_0(:,1),sum_1(:,1))

      select case(mod(parC(iocc)+parC(jocc),2))
      case(0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v)
         do ij=1,nidx_omg
            associate(cont => conts(ij))
              i = idx(1,ij)
              j = idx(2,ij)
              v = i + j
              select case(mod(v,2))
              case(0)
                 v = v/2+1
                 cont%par = 0
                 cont%n   = cont%n_0
                 call product2_se_unit(&
                      cont%n,auxC1%se_00,auxC1%val_00,v,prod2H%unit(i,j),&
                      cont%val)
              case(1)
                 v = (v+1)/2
                 cont%par = 1
                 cont%n   = cont%n_1
                 call product2_se_unit(&
                      cont%n,auxC1%se_11,auxC1%val_11,v,prod2H%unit(i,j),&
                      cont%val)
              end select
            end associate
         enddo
!$OMP END PARALLEL DO
      case(1)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v)
         do ij=1,nidx_omg
            associate(cont => conts(ij))
              i = idx(1,ij)
              j = idx(2,ij)
              v = i + j
              select case(mod(v,2))
              case(0)
                 v = v/2+1
                 cont%par = 1
                 cont%n   = cont%n_1
                 call product2_se_unit(&
                      cont%n,auxC1%se_10,auxC1%val_10,v,prod2H%unit(i,j),&
                      cont%val)
              case(1)
                 v = (v+1)/2
                 cont%par = 0
                 cont%n   = cont%n_0
                 call product2_se_unit(&
                      cont%n,auxC1%se_01,auxC1%val_01,v,prod2H%unit(i,j),&
                      cont%val)
              end select
            end associate
         enddo
!$OMP END PARALLEL DO
      end select

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(i,j,k,l,v1,v2) &
!$OMP PRIVATE(val_A1,val_A2,val_B1,val_B2,val_A,val_B)
      do kl=1,npair
         do ij=1,npair

            i = pair(1,ij)
            j = pair(2,ij)
            k = pair(1,kl)
            l = pair(2,kl)

            v1 = max(i,k)
            v1 = 1 + min(i,k) + (v1*(v1+1))/2
            v2 = j + l
            if(conts(v1)%par==mod(v2,2)) then
               v2 = max(j,l)
               v2 = 1 + min(j,l) + (v2*(v2+1))/2
               val_A1 = mlt*product1(conts(v1)%n,conts(v1)%val,ints(:,v2))
            else
               val_A1 = 0
            endif

            v1 = max(j,l)
            v1 = 1 + min(j,l) + (v1*(v1+1))/2
            v2 = i + k
            if(conts(v1)%par==mod(v2,2)) then
               v2 = max(i,k)
               v2 = 1 + min(i,k) + (v2*(v2+1))/2
               val_A2 = mlt*product1(conts(v1)%n,conts(v1)%val,ints(:,v2))
            else
               val_A2 = 0
            endif

            v1 = max(i,l)
            v1 = 1 + min(i,l) + (v1*(v1+1))/2
            v2 = j + k
            if(conts(v1)%par==mod(v2,2)) then
               v2 = max(j,k)
               v2 = 1 + min(j,k) + (v2*(v2+1))/2
               val_B1 = mlt*product1(conts(v1)%n,conts(v1)%val,ints(:,v2))
            else
               val_B1 = 0
            endif

            v1 = max(j,k)
            v1 = 1 + min(j,k) + (v1*(v1+1))/2
            v2 = i + l
            if(conts(v1)%par==mod(v2,2)) then
               v2 = max(i,l)
               v2 = 1 + min(i,l) + (v2*(v2+1))/2
               val_B2 = mlt*product1(conts(v1)%n,conts(v1)%val,ints(:,v2))
            else
               val_B2 = 0
            endif

            val_A = val_A1 + val_A2
            val_B = val_B1 + val_B2

            matJ_S(ij,kl,iocc,jocc) = val_A + val_B
            matJ_T(ij,kl,iocc,jocc) = val_A - val_B
            if(iocc/=jocc) then
               matJ_S(ij,kl,jocc,iocc) = matJ_S(ij,kl,iocc,jocc)
               matJ_T(ij,kl,jocc,iocc) = matJ_T(ij,kl,iocc,jocc)
            endif

            val_A = val_A1 - val_A2
            val_B = val_B1 - val_B2

            matM_S(ij,kl,iocc,jocc) = val_A - val_B
            matM_T(ij,kl,iocc,jocc) = val_A + val_B
            if(iocc/=jocc) then
               matM_S(ij,kl,jocc,iocc) = matM_S(ij,kl,iocc,jocc)
               matM_T(ij,kl,jocc,iocc) = matM_T(ij,kl,iocc,jocc)
            endif

         enddo
      enddo
!$OMP END PARALLEL DO

   enddo
enddo

do ij=1,nidx_omg
   associate(cont => conts(ij))
     deallocate(cont%val)
   end associate
enddo
deallocate(conts)
deallocate(sum_0,sum_1)
call aux_free(auxC1)

deallocate(ints)
call prod2H_free(prod2H)

deallocate(idx,pair)

if(Control%IPRINT>=1) call timer('J',Tcpu,Twall)

end subroutine integrals_J

subroutine integrals_P(matP1_S,matP1_T,matPe_S,matPe_T,orbE,matC,parC,Control)
implicit none
type(mp_real),intent(out) :: matP1_S(:,:),matP1_T(:,:)
type(mp_real),intent(out) :: matPe_S(:,:),matPe_T(:,:)
type(mp_real),intent(in) :: orbE(:),matC(:,:)
integer,intent(in) :: parC(:)
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu,Twall
integer :: nocc,nbas,omega,nidx_omg,npair,n1,n2,n3,nC
integer :: iocc,i,j,k,l,ij,kl,v1,v2,ip,jp
type(mp_real) :: val_A1,val_A2,val_B1,val_B2,val_A,val_B,val_S,val_T
integer,allocatable :: idx(:,:),pair(:,:)
type(prod2HData) :: prod2H
type(int3Data) :: int3
type(auxData) :: aux2H
type(mp_real),allocatable :: ints(:,:)

if(Control%IPRINT>=1) call timer('START',Tcpu,Twall)

nocc  = Control%nocc
nbas  = Control%nbas
omega = Control%omega
nidx_omg = ((omega+1)*(omega+2))/2
npair = Control%npair

allocate(idx(2,nidx_omg),pair(2,npair))
call fill_idx(idx,omega)
call fill_pair(pair,omega)

n1 = 2*omega
n2 = omega + nbas-1
n3 = omega + nbas-1
call int3_create('222',int3,n1,n2,n3,&
     Control%alpha,Control%nC,Control%C)

n1 = max(omega,nbas-1)
n2 = omega
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = omega + nbas-1
n2 = omega
nC = nbas-1
call aux_init(aux2H,n1,n2,nC)

n1 = max(int3%n1_0,int3%n1_1)
n2 = nidx_omg
allocate(ints(n1,n2))

!$OMP PARALLEL

!$OMP DO COLLAPSE(2)
do kl=1,npair
   do ij=1,npair
      matP1_S(ij,kl) = 0
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do kl=1,npair
   do ij=1,npair
      matP1_T(ij,kl) = 0
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do kl=1,npair
   do ij=1,npair
      matPe_S(ij,kl) = 0
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do kl=1,npair
   do ij=1,npair
      matPe_T(ij,kl) = 0
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

do iocc=1,nocc

   call aux2H_create(aux2H,matC(:,iocc),prod2H)

   select case(parC(iocc))
   case(0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,ip,jp)
      do ij=1,nidx_omg
         i = idx(1,ij)
         j = idx(2,ij)
         select case(10*mod(i,2) + mod(j,2))
         case(00)
            ip = i/2+1
            jp = j/2+1
            call product23_se(int3%n1_0,int3%val_000,&
                 aux2H%se_00(:,ip),aux2H%val_00(:,ip),&
                 aux2H%se_00(:,jp),aux2H%val_00(:,jp),&
                 ints(:,ij))
         case(10)
            ip = (i+1)/2
            jp = j/2+1
            call product23_se(int3%n1_1,int3%val_110,&
                 aux2H%se_11(:,ip),aux2H%val_11(:,ip),&
                 aux2H%se_00(:,jp),aux2H%val_00(:,jp),&
                 ints(:,ij))
         case(01)
            ip = i/2+1
            jp = (j+1)/2
            call product23_se(int3%n1_1,int3%val_101,&
                 aux2H%se_00(:,ip),aux2H%val_00(:,ip),&
                 aux2H%se_11(:,jp),aux2H%val_11(:,jp),&
                 ints(:,ij))
         case(11)
            ip = (i+1)/2
            jp = (j+1)/2
            call product23_se(int3%n1_0,int3%val_011,&
                 aux2H%se_11(:,ip),aux2H%val_11(:,ip),&
                 aux2H%se_11(:,jp),aux2H%val_11(:,jp),&
                 ints(:,ij))
         end select
      enddo
!$OMP END PARALLEL DO
   case(1)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,ip,jp)
      do ij=1,nidx_omg
         i = idx(1,ij)
         j = idx(2,ij)
         select case(10*mod(i,2) + mod(j,2))
         case(00)
            ip = i/2+1
            jp = j/2+1
            call product23_se(int3%n1_0,int3%val_011,&
                 aux2H%se_10(:,ip),aux2H%val_10(:,ip),&
                 aux2H%se_10(:,jp),aux2H%val_10(:,jp),&
                 ints(:,ij))
         case(10)
            ip = (i+1)/2
            jp = j/2+1
            call product23_se(int3%n1_1,int3%val_101,&
                 aux2H%se_01(:,ip),aux2H%val_01(:,ip),&
                 aux2H%se_10(:,jp),aux2H%val_10(:,jp),&
                 ints(:,ij))
         case(01)
            ip = i/2+1
            jp = (j+1)/2
            call product23_se(int3%n1_1,int3%val_110,&
                 aux2H%se_10(:,ip),aux2H%val_10(:,ip),&
                 aux2H%se_01(:,jp),aux2H%val_01(:,jp),&
                 ints(:,ij))
         case(11)
            ip = (i+1)/2
            jp = (j+1)/2
            call product23_se(int3%n1_0,int3%val_000,&
                 aux2H%se_01(:,ip),aux2H%val_01(:,ip),&
                 aux2H%se_01(:,jp),aux2H%val_01(:,jp),&
                 ints(:,ij))
         end select
      enddo
!$OMP END PARALLEL DO
   end select

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(i,j,k,l,v1,v2) &
!$OMP PRIVATE(val_A1,val_A2,val_B1,val_B2,val_A,val_B,val_S,val_T)
   do kl=1,npair
      do ij=1,npair

         i = pair(1,ij)
         j = pair(2,ij)
         k = pair(1,kl)
         l = pair(2,kl)

         if(mod(i+j+k+l,2)==0) then

            v1 = j + l
            v1 = (v1+2-mod(v1,2))/2
            v2 = max(i,k)
            v2 = 1 + min(i,k) + (v2*(v2+1))/2
            val_A1 = product1_unit(v1,prod2H%unit(j,l),ints(:,v2))

            v1 = i + k
            v1 = (v1+2-mod(v1,2))/2
            v2 = max(j,l)
            v2 = 1 + min(j,l) + (v2*(v2+1))/2
            val_A2 = product1_unit(v1,prod2H%unit(i,k),ints(:,v2))

            v1 = j + k
            v1 = (v1+2-mod(v1,2))/2
            v2 = max(i,l)
            v2 = 1 + min(i,l) + (v2*(v2+1))/2
            val_B1 = product1_unit(v1,prod2H%unit(j,k),ints(:,v2))

            v1 = i + l
            v1 = (v1+2-mod(v1,2))/2
            v2 = max(j,k)
            v2 = 1 + min(j,k) + (v2*(v2+1))/2
            val_B2 = product1_unit(v1,prod2H%unit(i,l),ints(:,v2))

            val_A = val_A1 + val_A2
            val_B = val_B1 + val_B2

            val_S = val_A + val_B
            val_T = val_A - val_B

            matP1_S(ij,kl) = matP1_S(ij,kl) + val_S
            matP1_T(ij,kl) = matP1_T(ij,kl) + val_T
            matPe_S(ij,kl) = matPe_S(ij,kl) + val_S*orbE(iocc)
            matPe_T(ij,kl) = matPe_T(ij,kl) + val_T*orbE(iocc)

         endif

      enddo
   enddo
!$OMP END PARALLEL DO

enddo

deallocate(ints)
call aux_free(aux2H)
call prod2H_free(prod2H)

call int3_free(int3)

deallocate(idx,pair)

if(Control%IPRINT>=1) call timer('P',Tcpu,Twall)

end subroutine integrals_P

subroutine integrals_L00(matL,Control)
implicit none
type(mp_real),intent(out) :: matL(:,:)
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu,Twall
integer :: omega,nidx_omg,npair,n1,n2
integer :: i,j,k,l,ij,kl,v,v1,v2
type(mp_real) :: mlt,val
integer,allocatable :: idx(:,:),pair(:,:)
type(prod2HData) :: prod2H
type(int2Data) :: int2
type(mp_real),allocatable :: ints(:,:)

if(Control%IPRINT>=1) call timer('START',Tcpu,Twall)

omega = Control%omega
nidx_omg = ((omega+1)*(omega+2))/2
npair = Control%npair

allocate(idx(2,nidx_omg),pair(2,npair))
call fill_idx(idx,omega)
call fill_pair(pair,omega)

n1 = omega
n2 = omega
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = 2*omega
n2 = 2*omega
call int2_create('f0','f0','22',int2,n1,n2,&
     Control%alpha,Control%nC,Control%C)
mlt = Control%gfac

n1 = omega+1
n2 = nidx_omg
allocate(ints(n1,n2))

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v)
do ij=1,nidx_omg
   i = idx(1,ij)
   j = idx(2,ij)
   v = i + j
   select case(mod(v,2))
   case(0)
      v = v/2+1
      call product2_unit(int2%n1_0,int2%val_0,v,prod2H%unit(i,j),&
           ints(:,ij))
   case(1)
      v = (v+1)/2
      call product2_unit(int2%n1_1,int2%val_1,v,prod2H%unit(i,j),&
           ints(:,ij))
   end select
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(i,j,k,l,v1,v2) &
!$OMP PRIVATE(val)
do kl=1,npair
   do ij=1,npair

      i = pair(1,ij)
      j = pair(2,ij)
      k = pair(1,kl)
      l = pair(2,kl)

      if(mod(i+j+k+l,2)==0) then

         v1 = i + k
         v1 = (v1+2-mod(v1,2))/2
         v2 = max(j,l)
         v2 = 1 + min(j,l) + (v2*(v2+1))/2
         val = mlt*product1_unit(v1,prod2H%unit(i,k),ints(:,v2))

         matL(ij,kl) = val

      else

         matL(ij,kl) = 0

      endif

   enddo
enddo
!$OMP END PARALLEL DO

deallocate(ints)

call int2_free(int2)

call prod2H_free(prod2H)

deallocate(idx,pair)

if(Control%IPRINT>=1) call timer('L00',Tcpu,Twall)

end subroutine integrals_L00

subroutine integrals_L10(matL,matC,parC,Control)
implicit none
type(mp_real),intent(out) :: matL(:,:)
type(mp_real),intent(in) :: matC(:,:)
integer,intent(in) :: parC(:)
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu,Twall
integer :: nocc,nbas,omega,npair,n1,n2,n3,nC
integer :: iocc,i,j,k,l,ij,kl,v,ip,jp
type(mp_real) :: mlt,val_1,val_2,val
integer,allocatable :: pair(:,:)
type(prod2HData) :: prod2H
type(int3Data) :: int3
type(auxData) :: aux2H
type(mp_real),allocatable :: ints(:,:,:)

if(Control%IPRINT>=1) call timer('START',Tcpu,Twall)

nocc  = Control%nocc
nbas  = Control%nbas
omega = Control%omega
npair = Control%npair

allocate(pair(2,npair))
call fill_pair(pair,omega)

n1 = 2*omega
n2 = omega + nbas-1
n3 = omega + nbas-1
call int3_from2_create('f ','f0',int3,n1,n2,n3,&
     Control%alpha,Control%nC,Control%C)
mlt = Control%gfac

n1 = max(omega,nbas-1)
n2 = omega
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = omega + nbas-1
n2 = omega
nC = nbas-1
call aux_init(aux2H,n1,n2,nC)

n1 = max(int3%n1_0,int3%n1_1)
n2 = omega
n3 = omega
allocate(ints(n1,0:n2,0:n3))

!$OMP PARALLEL DO COLLAPSE(2)
do kl=1,npair
   do ij=1,npair
      matL(ij,kl) = 0
   enddo
enddo
!$OMP END PARALLEL DO

do iocc=1,nocc

   call aux2H_create(aux2H,matC(:,iocc),prod2H)

   select case(parC(iocc))
   case(0)
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(ip,jp)
      do j=0,omega
         do i=0,omega
            select case(10*mod(i,2) + mod(j,2))
            case(00)
               ip = i/2+1
               jp = j/2+1
               call product23_se(int3%n1_0,int3%val_000,&
                    aux2H%se_00(:,ip),aux2H%val_00(:,ip),&
                    aux2H%se_00(:,jp),aux2H%val_00(:,jp),&
                    ints(:,i,j))
            case(10)
               ip = (i+1)/2
               jp = j/2+1
               call product23_se(int3%n1_1,int3%val_110,&
                    aux2H%se_11(:,ip),aux2H%val_11(:,ip),&
                    aux2H%se_00(:,jp),aux2H%val_00(:,jp),&
                    ints(:,i,j))
            case(01)
               ip = i/2+1
               jp = (j+1)/2
               call product23_se(int3%n1_1,int3%val_101,&
                    aux2H%se_00(:,ip),aux2H%val_00(:,ip),&
                    aux2H%se_11(:,jp),aux2H%val_11(:,jp),&
                    ints(:,i,j))
            case(11)
               ip = (i+1)/2
               jp = (j+1)/2
               call product23_se(int3%n1_0,int3%val_011,&
                    aux2H%se_11(:,ip),aux2H%val_11(:,ip),&
                    aux2H%se_11(:,jp),aux2H%val_11(:,jp),&
                    ints(:,i,j))
            end select
         enddo
      enddo
!$OMP END PARALLEL DO
   case(1)
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(ip,jp)
      do j=0,omega
         do i=0,omega
            select case(10*mod(i,2) + mod(j,2))
            case(00)
               ip = i/2+1
               jp = j/2+1
               call product23_se(int3%n1_0,int3%val_011,&
                    aux2H%se_10(:,ip),aux2H%val_10(:,ip),&
                    aux2H%se_10(:,jp),aux2H%val_10(:,jp),&
                    ints(:,i,j))
            case(10)
               ip = (i+1)/2
               jp = j/2+1
               call product23_se(int3%n1_1,int3%val_101,&
                    aux2H%se_01(:,ip),aux2H%val_01(:,ip),&
                    aux2H%se_10(:,jp),aux2H%val_10(:,jp),&
                    ints(:,i,j))
            case(01)
               ip = i/2+1
               jp = (j+1)/2
               call product23_se(int3%n1_1,int3%val_110,&
                    aux2H%se_10(:,ip),aux2H%val_10(:,ip),&
                    aux2H%se_01(:,jp),aux2H%val_01(:,jp),&
                    ints(:,i,j))
            case(11)
               ip = (i+1)/2
               jp = (j+1)/2
               call product23_se(int3%n1_0,int3%val_000,&
                    aux2H%se_01(:,ip),aux2H%val_01(:,ip),&
                    aux2H%se_01(:,jp),aux2H%val_01(:,jp),&
                    ints(:,i,j))
            end select
         enddo
      enddo
!$OMP END PARALLEL DO
   end select

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(i,j,k,l,v) &
!$OMP PRIVATE(val_1,val_2,val)   
   do kl=1,npair
      do ij=1,npair

         i = pair(1,ij)
         j = pair(2,ij)
         k = pair(1,kl)
         l = pair(2,kl)

         if(mod(i+j+k+l,2)==0) then

            v = k + l
            v = (v+2-mod(v,2))/2
            val_1 = mlt*product1_unit(v,prod2H%unit(k,l),ints(:,j,i))
            val_2 = mlt*product1_unit(v,prod2H%unit(k,l),ints(:,i,j))

            val = val_1 + val_2

            matL(ij,kl) = matL(ij,kl) + val

         endif

      enddo
   enddo
!$OMP END PARALLEL DO

enddo

deallocate(ints)
call aux_free(aux2H)
call prod2H_free(prod2H)

call int3_free(int3)

deallocate(pair)

if(Control%IPRINT>=1) call timer('L10',Tcpu,Twall)

end subroutine integrals_L10

subroutine integrals_L11(matL,matC,parC,Control)
implicit none
type(mp_real),intent(out) :: matL(:,:)
type(mp_real),intent(in) :: matC(:,:)
integer,intent(in) :: parC(:)
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu,Twall
integer :: nocc,nbas,omega,nidx_bas,nidx_omg,npair,ithr,n1,n2,n3,nC
integer :: iocc,jocc,i,j,k,l,ij,kl,v,v1,v2,ip,jp,n
type(mp_real) :: mlt,val,val_11,val_21,val_12,val_22
integer,allocatable :: idx(:,:),pair(:,:)
type(prod2HData) :: prod2H
type(int3Data) :: int3
type(auxData),allocatable :: aux2H(:)
type(auxData) :: auxC1
type(mp_real),allocatable :: ints(:,:)
type(mp_real),allocatable :: sum_0(:,:),sum_1(:,:)
type(mp_real),allocatable :: conts(:,:)

if(Control%IPRINT>=1) call timer('START',Tcpu,Twall)

nocc  = Control%nocc
nbas  = Control%nbas
omega = Control%omega
nidx_bas = (nbas*(nbas+1))/2
nidx_omg = ((omega+1)*(omega+2))/2
npair = Control%npair

allocate(idx(2,max(nidx_omg,nidx_bas)),pair(2,npair))
call fill_idx(idx,max(omega,nbas-1))
call fill_pair(pair,omega)

n1 = 2*(omega + nbas-1)
n2 = omega + nbas-1
n3 = omega + nbas-1
call int3_create('422',int3,n1,n2,n3,&
     Control%alpha,Control%nC,Control%C)
mlt = Control%gfac

n1 = max(omega,nbas-1)
n2 = max(omega,nbas-1)
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = omega + nbas-1
n2 = omega
nC = nbas-1
allocate(aux2H(nocc))
do iocc=1,nocc
   call aux_init(aux2H(iocc),n1,n2,nC)
enddo

n1 = 2*(omega + nbas-1)
n2 = 2* omega
nC = 2*(nbas-1)
call aux_init(auxC1,n1,n2,nC)

n1 = max(int3%n1_0,int3%n1_1)
n2 = nidx_omg
allocate(ints(n1,n2))

n1 = nbas
n2 = 1
!$ n2 = omp_get_max_threads()
allocate(sum_0(n1,n2),sum_1(n1,n2))

n1 = omega + 1
n2 = 1
!$ n2 = omp_get_max_threads()
allocate(conts(n1,n2))

!$OMP PARALLEL DO COLLAPSE(2)
do kl=1,npair
   do ij=1,npair
      matL(ij,kl) = 0
   enddo
enddo
!$OMP END PARALLEL DO

do iocc=1,nocc
   call aux2H_create(aux2H(iocc),matC(:,iocc),prod2H)
enddo

do jocc=1,nocc
   do iocc=1,nocc

      call sums_CC(nbas,nbas,sum_0,sum_1,&
           matC(:,iocc),matC(:,jocc),prod2H,&
           nidx_bas,idx)
      call auxC1_create(auxC1,sum_0(:,1),sum_1(:,1))

      select case(10*parC(iocc)+parC(jocc))
      case(00)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,ip,jp)
         do ij=1,nidx_omg
            i = idx(1,ij)
            j = idx(2,ij)
            select case(10*mod(i,2) + mod(j,2))
            case(00)
               ip = i/2+1
               jp = j/2+1
               call product23_se(int3%n1_0,int3%val_000,&
                    aux2H(iocc)%se_00(:,ip),aux2H(iocc)%val_00(:,ip),&
                    aux2H(jocc)%se_00(:,jp),aux2H(jocc)%val_00(:,jp),&
                    ints(:,ij))
            case(10)
               ip = (i+1)/2
               jp = j/2+1
               call product23_se(int3%n1_1,int3%val_110,&
                    aux2H(iocc)%se_11(:,ip),aux2H(iocc)%val_11(:,ip),&
                    aux2H(jocc)%se_00(:,jp),aux2H(jocc)%val_00(:,jp),&
                    ints(:,ij))
            case(01)
               ip = i/2+1
               jp = (j+1)/2
               call product23_se(int3%n1_1,int3%val_101,&
                    aux2H(iocc)%se_00(:,ip),aux2H(iocc)%val_00(:,ip),&
                    aux2H(jocc)%se_11(:,jp),aux2H(jocc)%val_11(:,jp),&
                    ints(:,ij))
            case(11)
               ip = (i+1)/2
               jp = (j+1)/2
               call product23_se(int3%n1_0,int3%val_011,&
                    aux2H(iocc)%se_11(:,ip),aux2H(iocc)%val_11(:,ip),&
                    aux2H(jocc)%se_11(:,jp),aux2H(jocc)%val_11(:,jp),&
                    ints(:,ij))
            end select
         enddo
!$OMP END PARALLEL DO
      case(10)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,ip,jp)
         do ij=1,nidx_omg
            i = idx(1,ij)
            j = idx(2,ij)
            select case(10*mod(i,2) + mod(j,2))
            case(00)
               ip = i/2+1
               jp = j/2+1
               call product23_se(int3%n1_1,int3%val_110,&
                    aux2H(iocc)%se_10(:,ip),aux2H(iocc)%val_10(:,ip),&
                    aux2H(jocc)%se_00(:,jp),aux2H(jocc)%val_00(:,jp),&
                    ints(:,ij))
            case(10)
               ip = (i+1)/2
               jp = j/2+1
               call product23_se(int3%n1_0,int3%val_000,&
                    aux2H(iocc)%se_01(:,ip),aux2H(iocc)%val_01(:,ip),&
                    aux2H(jocc)%se_00(:,jp),aux2H(jocc)%val_00(:,jp),&
                    ints(:,ij))
            case(01)
               ip = i/2+1
               jp = (j+1)/2
               call product23_se(int3%n1_0,int3%val_011,&
                    aux2H(iocc)%se_10(:,ip),aux2H(iocc)%val_10(:,ip),&
                    aux2H(jocc)%se_11(:,jp),aux2H(jocc)%val_11(:,jp),&
                    ints(:,ij))
            case(11)
               ip = (i+1)/2
               jp = (j+1)/2
               call product23_se(int3%n1_1,int3%val_101,&
                    aux2H(iocc)%se_01(:,ip),aux2H(iocc)%val_01(:,ip),&
                    aux2H(jocc)%se_11(:,jp),aux2H(jocc)%val_11(:,jp),&
                    ints(:,ij))
            end select
         enddo
!$OMP END PARALLEL DO
      case(01)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,ip,jp)
         do ij=1,nidx_omg
            i = idx(1,ij)
            j = idx(2,ij)
            select case(10*mod(i,2) + mod(j,2))
            case(00)
               ip = i/2+1
               jp = j/2+1
               call product23_se(int3%n1_1,int3%val_101,&
                    aux2H(iocc)%se_00(:,ip),aux2H(iocc)%val_00(:,ip),&
                    aux2H(jocc)%se_10(:,jp),aux2H(jocc)%val_10(:,jp),&
                    ints(:,ij))
            case(10)
               ip = (i+1)/2
               jp = j/2+1
               call product23_se(int3%n1_0,int3%val_011,&
                    aux2H(iocc)%se_11(:,ip),aux2H(iocc)%val_11(:,ip),&
                    aux2H(jocc)%se_10(:,jp),aux2H(jocc)%val_10(:,jp),&
                    ints(:,ij))
            case(01)
               ip = i/2+1
               jp = (j+1)/2
               call product23_se(int3%n1_0,int3%val_000,&
                    aux2H(iocc)%se_00(:,ip),aux2H(iocc)%val_00(:,ip),&
                    aux2H(jocc)%se_01(:,jp),aux2H(jocc)%val_01(:,jp),&
                    ints(:,ij))
            case(11)
               ip = (i+1)/2
               jp = (j+1)/2
               call product23_se(int3%n1_1,int3%val_110,&
                    aux2H(iocc)%se_11(:,ip),aux2H(iocc)%val_11(:,ip),&
                    aux2H(jocc)%se_01(:,jp),aux2H(jocc)%val_01(:,jp),&
                    ints(:,ij))
            end select
         enddo
!$OMP END PARALLEL DO
      case(11)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,ip,jp)
         do ij=1,nidx_omg
            i = idx(1,ij)
            j = idx(2,ij)
            select case(10*mod(i,2) + mod(j,2))
            case(00)
               ip = i/2+1
               jp = j/2+1
               call product23_se(int3%n1_0,int3%val_011,&
                    aux2H(iocc)%se_10(:,ip),aux2H(iocc)%val_10(:,ip),&
                    aux2H(jocc)%se_10(:,jp),aux2H(jocc)%val_10(:,jp),&
                    ints(:,ij))
            case(10)
               ip = (i+1)/2
               jp = j/2+1
               call product23_se(int3%n1_1,int3%val_101,&
                    aux2H(iocc)%se_01(:,ip),aux2H(iocc)%val_01(:,ip),&
                    aux2H(jocc)%se_10(:,jp),aux2H(jocc)%val_10(:,jp),&
                    ints(:,ij))
            case(01)
               ip = i/2+1
               jp = (j+1)/2
               call product23_se(int3%n1_1,int3%val_110,&
                    aux2H(iocc)%se_10(:,ip),aux2H(iocc)%val_10(:,ip),&
                    aux2H(jocc)%se_01(:,jp),aux2H(jocc)%val_01(:,jp),&
                    ints(:,ij))
            case(11)
               ip = (i+1)/2
               jp = (j+1)/2
               call product23_se(int3%n1_0,int3%val_000,&
                    aux2H(iocc)%se_01(:,ip),aux2H(iocc)%val_01(:,ip),&
                    aux2H(jocc)%se_01(:,jp),aux2H(jocc)%val_01(:,jp),&
                    ints(:,ij))
            end select
         enddo
!$OMP END PARALLEL DO
      end select

      select case(mod(parC(iocc)+parC(jocc),2))
      case(0)
!$OMP PARALLEL PRIVATE(ithr,i,j,v,n)
         ithr = 0
!$ ithr = omp_get_thread_num()
         associate(cont => conts(:,ithr+1))
!$OMP DO SCHEDULE(DYNAMIC)
           do ij=1,nidx_omg
              i = idx(1,ij)
              j = idx(2,ij)
              v = i + j
              select case(mod(v,2))
              case(0)
                 n = auxC1%n2_0
                 do k=1,n
                    cont(k) = product1_se(&
                         auxC1%se_00(:,k),auxC1%val_00(:,k),ints(:,ij))
                 enddo
              case(1)
                 n = auxC1%n2_1
                 do k=1,n
                    cont(k) = product1_se(&
                         auxC1%se_11(:,k),auxC1%val_11(:,k),ints(:,ij))
                 enddo
              end select
              do k=1,n
                 ints(k,ij) = cont(k)
              enddo
           enddo
!$OMP END DO
         end associate
!$OMP END PARALLEL
      case(1)
!$OMP PARALLEL PRIVATE(ithr,i,j,v,n)
         ithr = 0
!$ ithr = omp_get_thread_num()
         associate(cont => conts(:,ithr+1))
!$OMP DO SCHEDULE(DYNAMIC)
           do ij=1,nidx_omg
              i = idx(1,ij)
              j = idx(2,ij)
              v = i + j
              select case(mod(v,2))
              case(0)
                 n = auxC1%n2_0
                 do k=1,n
                    cont(k) = product1_se(&
                         auxC1%se_10(:,k),auxC1%val_10(:,k),ints(:,ij))
                 enddo
              case(1)
                 n = auxC1%n2_1
                 do k=1,n
                    cont(k) = product1_se(&
                         auxC1%se_01(:,k),auxC1%val_01(:,k),ints(:,ij))
                 enddo
              end select
              do k=1,n
                 ints(k,ij) = cont(k)
              enddo
           enddo
!$OMP END DO
         end associate
!$OMP END PARALLEL
      end select

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(i,j,k,l,v1,v2) &
!$OMP PRIVATE(val_11,val_21,val_12,val_22,val)
      do kl=1,npair
         do ij=1,npair

            i = pair(1,ij)
            j = pair(2,ij)
            k = pair(1,kl)
            l = pair(2,kl)

            if(mod(i+j+k+l,2)==0) then

               v1 = j + l
               v1 = (v1+2-mod(v1,2))/2
               v2 = max(i,k)
               v2 = 1 + min(i,k) + (v2*(v2+1))/2
               val_11 = mlt*product1_unit(v1,prod2H%unit(j,l),ints(:,v2))

               v1 = i + l
               v1 = (v1+2-mod(v1,2))/2
               v2 = max(j,k)
               v2 = 1 + min(j,k) + (v2*(v2+1))/2
               val_21 = mlt*product1_unit(v1,prod2H%unit(i,l),ints(:,v2))

               v1 = j + k
               v1 = (v1+2-mod(v1,2))/2
               v2 = max(i,l)
               v2 = 1 + min(i,l) + (v2*(v2+1))/2
               val_12 = mlt*product1_unit(v1,prod2H%unit(j,k),ints(:,v2))

               v1 = i + k
               v1 = (v1+2-mod(v1,2))/2
               v2 = max(j,l)
               v2 = 1 + min(j,l) + (v2*(v2+1))/2
               val_22 = mlt*product1_unit(v1,prod2H%unit(i,k),ints(:,v2))

               val = val_11 + val_21 + val_12 + val_22

               matL(ij,kl) = matL(ij,kl) + val

            endif

         enddo
      enddo
!$OMP END PARALLEL DO

   enddo
enddo

deallocate(conts)
deallocate(sum_0,sum_1)
deallocate(ints)

call aux_free(auxC1)
do iocc=1,nocc
   call aux_free(aux2H(iocc))
enddo
deallocate(aux2H)
call prod2H_free(prod2H)

call int3_free(int3)

deallocate(idx,pair)

if(Control%IPRINT>=1) call timer('L11',Tcpu,Twall)

end subroutine integrals_L11

subroutine integrals_vecX(typeX,vecX,matC,parC,Control)
implicit none
character(1),intent(in) :: typeX
type(mp_real),intent(out) :: vecX(:,:,:)
type(mp_real),intent(in) :: matC(:,:)
integer,intent(in) :: parC(:)
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu,Twall
integer :: nocc,nbas,omega,npair,n1,n2,nC
integer :: iocc,jocc,ij,i,j,v
type(mp_real) :: mlt,val
integer,allocatable :: pair(:,:)
type(prod2HData) :: prod2H
type(int2Data) :: int2
type(auxData),allocatable,target :: aux2H(:)
type(mp_real),allocatable :: ints(:,:)
integer,allocatable :: par_ints(:)
type(mp_real),pointer :: cont_0(:,:),cont_1(:,:)
integer,pointer :: se_0(:,:),se_1(:,:)
integer,allocatable :: par_cont(:)

if(Control%IPRINT>=1) call timer('START',Tcpu,Twall)

nocc  = Control%nocc
nbas  = Control%nbas
omega = Control%omega
npair = Control%npair

allocate(pair(2,npair))
call fill_pair(pair,omega)

n1 = nbas-1
n2 = omega
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = omega + nbas-1
n2 = omega + nbas-1
select case(typeX)
case('1')
   call int2_create('f0','1 ','22',int2,n1,n2,&
        Control%alpha,Control%nC,Control%C)
   mlt = Control%gfac
case('0')
   call int2_create('f ','1 ','22',int2,n1,n2,&
        Control%alpha,Control%nC,Control%C)
   mlt = 1
end select

n1 = omega + nbas-1
n2 = omega
nC = nbas-1
allocate(aux2H(nocc))
do iocc=1,nocc
   call aux_init(aux2H(iocc),n1,n2,nC)
   call aux2H_create(aux2H(iocc),matC(:,iocc),prod2H)
enddo

n1 = max(int2%n1_0,int2%n1_1)
n2 = omega
allocate(ints(n1,0:n2),par_ints(0:n2))
allocate(par_cont(0:n2))

do jocc=1,nocc
   associate(jaux => aux2H(jocc))

     select case(parC(jocc))
     case(0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(v)
        do i=0,omega
           select case(mod(i,2))
           case(0)
              v = i/2+1
              par_ints(i) = 0
              call product2_se(int2%n1_0,int2%val_0,&
                   jaux%se_00(:,v),jaux%val_00(:,v),ints(:,i))
           case(1)
              v = (i+1)/2
              par_ints(i) = 1
              call product2_se(int2%n1_1,int2%val_1,&
                   jaux%se_11(:,v),jaux%val_11(:,v),ints(:,i))
           end select
        enddo
!$OMP END PARALLEL DO
     case(1)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(v)
        do i=0,omega
           select case(mod(i,2))
           case(0)
              v = i/2+1
              par_ints(i) = 1
              call product2_se(int2%n1_1,int2%val_1,&
                   jaux%se_10(:,v),jaux%val_10(:,v),ints(:,i))
           case(1)
              v = (i+1)/2
              par_ints(i) = 0
              call product2_se(int2%n1_0,int2%val_0,&
                   jaux%se_01(:,v),jaux%val_01(:,v),ints(:,i))
           end select
        enddo
!$OMP END PARALLEL DO
     end select

     do iocc=1,nocc
        associate(iaux => aux2H(iocc))

          select case(parC(iocc))
          case(0)
             se_0   => iaux%se_00
             cont_0 => iaux%val_00
             se_1   => iaux%se_11
             cont_1 => iaux%val_11
             par_cont(0:omega:2) = 0
             par_cont(1:omega:2) = 1
          case(1)
             se_0   => iaux%se_10
             cont_0 => iaux%val_10
             se_1   => iaux%se_01
             cont_1 => iaux%val_01
             par_cont(0:omega:2) = 1
             par_cont(1:omega:2) = 0
          end select

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v) &
!$OMP PRIVATE(val)
          do ij=1,npair

             i = pair(1,ij)
             j = pair(2,ij)

             if(par_cont(i)==par_ints(j)) then

                select case(mod(i,2))
                case(0)
                   v = i/2+1
                   val = mlt*product1_se(se_0(:,v),cont_0(:,v),ints(:,j))
                case(1)
                   v = (i+1)/2
                   val = mlt*product1_se(se_1(:,v),cont_1(:,v),ints(:,j))
                end select

                vecX(ij,iocc,jocc) = val

             else

                vecX(ij,iocc,jocc) = 0

             endif

          enddo
!$OMP END PARALLEL DO

        end associate
     enddo

   end associate
enddo

deallocate(par_cont)
deallocate(ints,par_ints)

do iocc=1,nocc
   call aux_free(aux2H(iocc))
enddo
deallocate(aux2H)

call int2_free(int2)

call prod2H_free(prod2H)

deallocate(pair)

if(Control%IPRINT>=1) call timer('vec'//typeX,Tcpu,Twall)

end subroutine integrals_vecX

subroutine integrals_vecP(vecP,matC,parC,Control)
implicit none
type(mp_real),intent(out) :: vecP(:,:,:)
type(mp_real),intent(in) :: matC(:,:)
integer,intent(in) :: parC(:)
type(ControlData),intent(in) :: Control
real(8) :: Tcpu,Twall
integer :: nocc,nbas,omega,npair,nocc3,nidx_bas,n1,n2,nC
integer :: iocc,jocc,kocc,ijkocc,i,j,ij,v
type(mp_real) :: mlt,val,val_1,val_2
integer,allocatable :: idx(:,:),pair(:,:)
type(prod2HData) :: prod2H
type(int2Data) :: int2
type(auxData),allocatable :: aux2H(:)
type(auxData),allocatable,target :: auxC2(:)
type(mp_real),allocatable :: sum_0(:,:),sum_1(:,:)
type(mp_real),allocatable :: ints(:,:)
integer,allocatable :: par_ints(:)
type(mp_real),pointer :: cont_0(:,:),cont_1(:,:)
integer,pointer :: se_0(:,:),se_1(:,:)
integer,allocatable :: par_cont(:)

if(Control%IPRINT>=1) call timer('START',Tcpu,Twall)

nocc  = Control%nocc
nbas  = Control%nbas
omega = Control%omega
npair = Control%npair
nocc3    = (nocc*(nocc+1)*(nocc+2))/6
nidx_bas = (nbas*(nbas+1))/2

allocate(idx(2,nidx_bas),pair(2,npair))
call fill_idx(idx,nbas-1)
call fill_pair(pair,omega)

n1 = omega + nbas-1
n2 = omega
nC = nbas-1
allocate(aux2H(nocc))
do iocc=1,nocc
   call aux_init(aux2H(iocc),n1,n2,nC)
enddo

n1 = omega + 3*(nbas-1)
n2 = omega
nC = 2*(nbas-1)
allocate(auxC2(nocc3))
do ijkocc=1,nocc3
   call aux_init(auxC2(ijkocc),n1,n2,nC)
enddo

n1 = nbas-1
n2 = max(omega,nbas-1)
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = nbas
n2 = 1
!$ n2 = omp_get_max_threads()
allocate(sum_0(n1,n2),sum_1(n1,n2))

do iocc=1,nocc
   call aux2H_create(aux2H(iocc),matC(:,iocc),prod2H)
enddo

ijkocc = 0
do jocc=1,nocc
   do iocc=1,jocc

      call sums_CC(nbas,nbas,sum_0,sum_1,&
           matC(:,iocc),matC(:,jocc),prod2H,&
           nidx_bas,idx)

      do kocc=1,iocc
         ijkocc = ijkocc + 1
         call auxC2_create(auxC2(ijkocc),sum_0(:,1),sum_1(:,1),aux2H(kocc))
      enddo

   enddo
enddo

deallocate(sum_0,sum_1)

call prod2H_free(prod2H)

n1 = omega + 3*(nbas-1)
n2 = omega + nbas-1
call int2_create('f ','1 ','42',int2,n1,n2,&
     Control%alpha,Control%nC,Control%C)
mlt = Control%gfac

n1 = max(int2%n1_0,int2%n1_1)
n2 = omega
allocate(ints(n1,0:n2),par_ints(0:n2))
allocate(par_cont(0:n2))

!$OMP PARALLEL DO COLLAPSE(3)
do jocc=1,nocc
   do iocc=1,nocc
      do ij=1,npair
         vecP(ij,iocc,jocc) = 0
      enddo
   enddo
enddo
!$OMP END PARALLEL DO

do kocc=1,nocc
   associate(kaux => aux2H(kocc))

     select case(parC(kocc))
     case(0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(v)
        do i=0,omega
           select case(mod(i,2))
           case(0)
              v = i/2+1
              par_ints(i) = 0
              call product2_se(int2%n1_0,int2%val_0,&
                   kaux%se_00(:,v),kaux%val_00(:,v),ints(:,i))
           case(1)
              v = (i+1)/2
              par_ints(i) = 1
              call product2_se(int2%n1_1,int2%val_1,&
                   kaux%se_11(:,v),kaux%val_11(:,v),ints(:,i))
           end select
        enddo
!$OMP END PARALLEL DO
     case(1)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(v)
        do i=0,omega
           select case(mod(i,2))
           case(0)
              v = i/2+1
              par_ints(i) = 1
              call product2_se(int2%n1_1,int2%val_1,&
                   kaux%se_10(:,v),kaux%val_10(:,v),ints(:,i))
           case(1)
              v = (i+1)/2
              par_ints(i) = 0
              call product2_se(int2%n1_0,int2%val_0,&
                   kaux%se_01(:,v),kaux%val_01(:,v),ints(:,i))
           end select
        enddo
!$OMP END PARALLEL DO
     end select

     do jocc=1,nocc
        do iocc=1,jocc
           if(kocc<=iocc) then
              ijkocc = kocc + ((iocc-1)*iocc)/2 + ((jocc-1)*jocc*(jocc+1))/6
           elseif(kocc<=jocc) then
              ijkocc = iocc + ((kocc-1)*kocc)/2 + ((jocc-1)*jocc*(jocc+1))/6
           else
              ijkocc = iocc + ((jocc-1)*jocc)/2 + ((kocc-1)*kocc*(kocc+1))/6
           endif
           associate(ijkaux => auxC2(ijkocc))

             select case(mod(parC(iocc)+parC(jocc)+parC(kocc),2))
             case(0)
                se_0   => ijkaux%se_00
                cont_0 => ijkaux%val_00
                se_1   => ijkaux%se_11
                cont_1 => ijkaux%val_11
                par_cont(0:omega:2) = 0
                par_cont(1:omega:2) = 1
             case(1)
                se_0   => ijkaux%se_10
                cont_0 => ijkaux%val_10
                se_1   => ijkaux%se_01
                cont_1 => ijkaux%val_01
                par_cont(0:omega:2) = 1
                par_cont(1:omega:2) = 0
             end select

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v) &
!$OMP PRIVATE(val_1,val_2,val)
             do ij=1,npair

                i = pair(1,ij)
                j = pair(2,ij)

                if(par_cont(i)==par_ints(j)) then

                   select case(mod(j,2))
                   case(0)
                      v = j/2+1
                      val_1 = mlt*product1_se(se_0(:,v),cont_0(:,v),ints(:,i))
                   case(1)
                      v = (j+1)/2
                      val_1 = mlt*product1_se(se_1(:,v),cont_1(:,v),ints(:,i))
                   end select

                   select case(mod(i,2))
                   case(0)
                      v = i/2+1
                      val_2 = mlt*product1_se(se_0(:,v),cont_0(:,v),ints(:,j))
                   case(1)
                      v = (i+1)/2
                      val_2 = mlt*product1_se(se_1(:,v),cont_1(:,v),ints(:,j))
                   end select

                   val = val_1 + val_2

                   vecP(ij,iocc,jocc) = vecP(ij,iocc,jocc) + val
                   if(iocc/=jocc) vecP(ij,jocc,iocc) = vecP(ij,jocc,iocc) + val

                endif

             enddo
!$OMP END PARALLEL DO

           end associate
        enddo
     enddo

   end associate
enddo

deallocate(par_cont)
deallocate(ints,par_ints)

call int2_free(int2)

do ijkocc=1,nocc3
   call aux_free(auxC2(ijkocc))
enddo
deallocate(auxC2)
do iocc=1,nocc
   call aux_free(aux2H(iocc))
enddo
deallocate(aux2H)

deallocate(idx,pair)

if(Control%IPRINT>=1) call timer('vecP',Tcpu,Twall)

end subroutine integrals_vecP

subroutine sums_D(n_0,n_1,sum_0,sum_1,D,ints,nidx,idx)
implicit none
integer,intent(in) :: n_0,n_1
type(mp_real),intent(out) :: sum_0(:,0:),sum_1(:,0:)
type(mp_real),intent(in) :: D(0:,0:)
type(mp_real),intent(in) :: ints(:,:)
integer,intent(in) :: nidx,idx(:,:)
integer :: nthr,ithr,ij,i,j,v
type(mp_real) :: val

!$OMP PARALLEL PRIVATE(ithr,i,j,v,val)

!$OMP SINGLE
nthr = 1
!$ nthr = omp_get_num_threads()
!$OMP END SINGLE

!$OMP DO COLLAPSE(2)
do j=0,nthr-1
   do i=1,n_0
      sum_0(i,j) = 0
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do j=0,nthr-1
   do i=1,n_1
      sum_1(i,j) = 0
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP BARRIER
      
!$OMP DO SCHEDULE(DYNAMIC)
do ij=1,nidx
   ithr = 0
!$ ithr = omp_get_thread_num()
   i = idx(1,ij)
   j = idx(2,ij)
   v = i + j
   select case(mod(v,2))
   case(0)
      val = D(i,j)
      if(i/=j) val = val + D(j,i)
      call accumulate(val,n_0,ints(:,ij),sum_0(:,ithr))
   case(1)
      val = D(i,j) + D(j,i)
      call accumulate(val,n_1,ints(:,ij),sum_1(:,ithr))
   end select
enddo
!$OMP END DO

!$OMP SECTIONS
!$OMP SECTION
do ithr=1,nthr-1
   call addition(n_0,sum_0(:,ithr),sum_0(:,0))
enddo
!$OMP SECTION
do ithr=1,nthr-1
   call addition(n_1,sum_1(:,ithr),sum_1(:,0))
enddo
!$OMP END SECTIONS

!$OMP END PARALLEL

end subroutine sums_D

subroutine sums_CC(n_0,n_1,sum_0,sum_1,C1,C2,prod2H,nidx,idx)
implicit none
integer,intent(in) :: n_0,n_1
type(mp_real),intent(out) :: sum_0(:,0:),sum_1(:,0:)
type(mp_real),intent(in) :: C1(0:),C2(0:)
type(prod2HData),intent(in) :: prod2H
integer,intent(in) :: nidx,idx(:,:)
integer :: nthr,ithr,ij,i,j,v
type(mp_real) :: val

!$OMP PARALLEL PRIVATE(ithr,i,j,v,val)

!$OMP SINGLE
nthr = 1
!$ nthr = omp_get_num_threads()
!$OMP END SINGLE

!$OMP DO COLLAPSE(2)
do j=0,nthr-1
   do i=1,n_0
      sum_0(i,j) = 0
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do j=0,nthr-1
   do i=1,n_1
      sum_1(i,j) = 0
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP BARRIER
      
!$OMP DO SCHEDULE(DYNAMIC)
do ij=1,nidx
   ithr = 0
!$ ithr = omp_get_thread_num()
   i = idx(1,ij)
   j = idx(2,ij)
   v = i + j
   select case(mod(v,2))
   case(0)
      v = v/2+1
      val = C1(i)*C2(j)
      if(i/=j) val = val + C1(j)*C2(i)
      call accumulate_unit(val,v,prod2H%unit(i,j),sum_0(:,ithr))
   case(1)
      v = (v+1)/2
      val = C1(i)*C2(j) + C1(j)*C2(i)
      call accumulate_unit(val,v,prod2H%unit(i,j),sum_1(:,ithr))
   end select
enddo
!$OMP END DO

!$OMP SECTIONS
!$OMP SECTION
do ithr=1,nthr-1
   call addition(n_0,sum_0(:,ithr),sum_0(:,0))
enddo
!$OMP SECTION
do ithr=1,nthr-1
   call addition(n_1,sum_1(:,ithr),sum_1(:,0))
enddo
!$OMP END SECTIONS

!$OMP END PARALLEL

end subroutine sums_CC

subroutine fill_idx(idx,n)
implicit none
integer,intent(out) :: idx(:,:)
integer,intent(in) :: n
integer :: i,j,ij

ij = 0
do j=0,n
   do i=0,j
      ij = ij + 1
      idx(1,ij) = i
      idx(2,ij) = j
   enddo
enddo

end subroutine fill_idx

subroutine fill_pair(pair,omega)
implicit none
integer,intent(out) :: pair(:,:)
integer,intent(in) :: omega
integer :: i,j,ij

ij = 0; i = -1; j = -1
do while(next_pair(i,j,omega))
   ij = ij + 1
   pair(1,ij) = i
   pair(2,ij) = j
enddo

end subroutine fill_pair

end module integrals
