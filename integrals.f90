module integrals
!$ use omp_lib
use mpmodule
use general
use time
implicit none

private
public integrals_scfSH
public scfJauxData,prepare_scfJaux,free_scfJaux,integrals_scfJ
public integrals_SH,integrals_J
public integrals_vec0,integrals_vecP
public integrals_f12,integrals_matP
public outer_vecproduct,fullMO_transform,scf_FDS_SDF

type table2Data
integer :: n1_0,n2_0
integer :: n1_1,n2_1
type(mp_real),allocatable :: val_0(:,:),val_1(:,:)
end type table2Data

type table4Data
integer :: n1_0,n2_0
integer :: n1_1,n2_1
integer,allocatable :: se_00(:,:),se_10(:,:),se_01(:,:),se_11(:,:)
type(mp_real),allocatable :: val_00(:,:),val_10(:,:),val_01(:,:),val_11(:,:)
end type table4Data

type contData
integer :: n_0,n_1
integer :: par,n
type(mp_real),allocatable :: val(:)
end type contData

type prodUnitData
integer :: nu
type(mp_real),allocatable :: val(:)
end type prodUnitData

type prod2HData
integer :: n1,n2
type(prodUnitData),allocatable :: unit(:,:)
end type prod2HData

type scfJauxData
type(prod2HData) :: prod2H
type(mp_real),allocatable :: ints(:,:)
type(mp_real),allocatable :: sum_0(:,:),sum_1(:,:)
integer :: n_0,n_1
end type scfJauxData

contains

include 'integrals/scf.f90'
include 'integrals/sh.f90'
include 'integrals/j.f90'

subroutine integrals_matP(IPRINT)
  implicit none
  integer,intent(in) :: IPRINT
  real(8) :: TCPU,Twall

  if(IPRINT>=1) call timer('START',Tcpu,Twall)


  if(IPRINT>=1) call timer('matP',Tcpu,Twall)
  
end subroutine integrals_matP

subroutine integrals_vec0(vec0,matC,parC,basicF12,IPRINT)
implicit none
type(mp_real),intent(out) :: vec0(:,:,:)
type(mp_real),intent(in) :: matC(:,:)
integer,intent(in) :: parC(:)
integer,intent(in) :: IPRINT
real(8) :: TCPU,Twall
integer :: n1,n2,nC,iocc,jocc,i,j,ij,v1,v2
type(mp_real) :: a1,a2,a12,int0,c0,mult
type(mp_real) :: val
type(prod2HData) :: prod2H
type(table2Data) :: aux,auxOP
type(table4Data),allocatable,target :: aux2C(:)
integer,allocatable :: pair_list(:,:)
type(mp_real),allocatable :: ints(:,:)
integer,allocatable :: par_ints(:)
type(mp_real),pointer :: cont_0(:,:),cont_1(:,:)
integer,pointer :: se_0(:,:),se_1(:,:)
integer,allocatable :: par_cont(:)
type(mp_real) :: basicF12(0:)

if(IPRINT>=1) call timer('START',Tcpu,Twall)

n1 = G_nbas-1
n2 = G_nprim
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = G_nprim + (G_nbas-1)
n2 = G_nprim + (G_nbas-1)
! a1 =  mpreal(-0.5d0)
! a2 =  mpreal(-0.5d0)
! a12 = mpreal(1.d0)
! call table2_create(aux,n1,n2,a1,a2,a12)

! int0 = sqrt(2*mppi())/(1 + 2*G_alpha)
! c0   = 1/(1 + 2*G_alpha)
! mult = mpreal(-1.d0)
! call auxXXF1_create(auxOP,int0,c0,mult,aux)

! call table2_fact(auxOP)
! call table2_free(aux)

call table2_easy(auxOP,basicF12,n1,n2)

n1 = G_nprim + (G_nbas-1)
n2 = G_nprim
nC = G_nbas-1
allocate(aux2C(G_nocc))
do iocc=1,G_nocc
   call aux2C_create(aux2C(iocc),n1,n2,nC,matC(:,iocc),prod2H)
enddo

n1 = 2
n2 = G_npair
allocate(pair_list(n1,n2))
i = -1
j = -1
ij = 0
do while(G_next_pair(i,j))
   ij = ij + 1
   pair_list(1,ij) = i
   pair_list(2,ij) = j
enddo

n1 = max(auxOP%n1_0,auxOP%n1_1)
n2 = G_nprim
allocate(ints(n1,0:n2),par_ints(0:n2),par_cont(0:n2))

do jocc=1,G_nocc
   associate(jaux => aux2C(jocc))

     select case(parC(jocc))
     case(0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(v2)
        do i=0,G_nprim
           select case(mod(i,2))
           case(0)
              v2 = i/2+1
              par_ints(i) = 0
              call product2_se(auxOP%n1_0,auxOP%val_0,&
                   jaux%se_00(:,v2),jaux%val_00(:,v2),ints(:,i))
           case(1)
              v2 = (i+1)/2
              par_ints(i) = 1
              call product2_se(auxOP%n1_1,auxOP%val_1,&
                   jaux%se_11(:,v2),jaux%val_11(:,v2),ints(:,i))
           end select
        enddo
!$OMP END PARALLEL DO
     case(1)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(v2)
        do i=0,G_nprim
           select case(mod(i,2))
           case(0)
              v2 = i/2+1
              par_ints(i) = 1
              call product2_se(auxOP%n1_1,auxOP%val_1,&
                   jaux%se_10(:,v2),jaux%val_10(:,v2),ints(:,i))
           case(1)
              v2 = (i+1)/2
              par_ints(i) = 0
              call product2_se(auxOP%n1_0,auxOP%val_0,&
                   jaux%se_01(:,v2),jaux%val_01(:,v2),ints(:,i))
           end select
        enddo
!$OMP END PARALLEL DO
     end select

     do iocc=1,G_nocc
        associate(iaux => aux2C(iocc))

          select case(parC(iocc))
          case(0)
             se_0   => iaux%se_00
             cont_0 => iaux%val_00
             se_1   => iaux%se_11
             cont_1 => iaux%val_11
             par_cont(0:G_nprim:2) = 0
             par_cont(1:G_nprim:2) = 1
          case(1)
             se_0   => iaux%se_10
             cont_0 => iaux%val_10
             se_1   => iaux%se_01
             cont_1 => iaux%val_01
             par_cont(0:G_nprim:2) = 1
             par_cont(1:G_nprim:2) = 0
          end select

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v1,val)
          do ij=1,G_npair

             i = pair_list(1,ij)
             j = pair_list(2,ij)

             if(par_cont(i)==par_ints(j)) then
                select case(mod(i,2))
                case(0)
                   v1 = i/2+1
                   val = product1_se(se_0(:,v1),cont_0(:,v1),ints(:,j))
                case(1)
                   v1 = (i+1)/2
                   val = product1_se(se_1(:,v1),cont_1(:,v1),ints(:,j))
                end select
             else
                val = 0
             endif

             vec0(ij,iocc,jocc) = val

          enddo
!$OMP END PARALLEL DO

        end associate
     enddo

   end associate
enddo

deallocate(ints,par_ints,par_cont)
deallocate(pair_list)

do iocc=1,G_nocc
   call table4_free(aux2C(iocc))
enddo
deallocate(aux2C)

call table2_free(auxOP)
call prod2H_free(prod2H)

if(IPRINT>=1) call timer('vec0',Tcpu,Twall)

end subroutine integrals_vec0

subroutine integrals_vecP(vecP,matC,parC,IPRINT)
implicit none
type(mp_real),intent(out) :: vecP(:,:,:)
type(mp_real),intent(in) :: matC(:,:)
integer,intent(in) :: parC(:)
integer,intent(in) :: IPRINT
real(8) :: TCPU,Twall
integer :: nthr,ithr,n1,n2,nC,iocc,jocc,kocc,ijkocc,i,j,ij,v1,v2
real(8) :: safe_nbas
type(mp_real) :: a1,a2,a12,int0,c0,mult
type(mp_real) :: val,val_1,val_2
type(prod2HData) :: prod2H
type(table2Data) :: aux,auxOP
type(table4Data),allocatable :: aux2C(:)
type(table4Data),allocatable,target :: aux22C(:)
type(mp_real),allocatable :: CC(:,:),sum_0(:,:),sum_1(:,:)
integer,allocatable :: pair_list(:,:)
type(mp_real),allocatable :: ints(:,:)
integer,allocatable :: par_ints(:)
type(mp_real),pointer :: cont_0(:,:),cont_1(:,:)
integer,pointer :: se_0(:,:),se_1(:,:)
integer,allocatable :: par_cont(:)

if(IPRINT>=1) call timer('START',Tcpu,Twall)

safe_nbas = 0.5d0 - 0.5d0/G_nbas

n1 = G_nbas-1
n2 = max(G_nprim,G_nbas-1)
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = G_nprim + 3*(G_nbas-1)
n2 = G_nprim + (G_nbas-1)
a1 =  mpreal(-2.d0)/3
a2 =  mpreal(-2.d0)/3
a12 = mpreal(2.d0)/3
call table2_create(aux,n1,n2,a1,a2,a12)

int0 = sqrt(3*mppi())/(2 + 3*G_alpha)
c0   = (mpreal(2.d0)/3)/(2 + 3*G_alpha)
mult = mpreal(-2.d0)
call auxXXF1_create(auxOP,int0,c0,mult,aux)
call table2_fact(auxOP)

call table2_free(aux)

n1 = G_nprim + (G_nbas-1)
n2 = G_nprim
nC = G_nbas-1
allocate(aux2C(G_nocc))
do iocc=1,G_nocc
   call aux2C_create(aux2C(iocc),n1,n2,nC,matC(:,iocc),prod2H)
enddo

allocate(aux22C((G_nocc*(G_nocc+1)*(G_nocc+2))/6))

allocate(CC(0:G_nbas-1,0:G_nbas-1))

n1 = G_nbas
n2 = 1
!$ n2 = omp_get_max_threads()
allocate(sum_0(n1,0:n2-1),sum_1(n1,0:n2-1))

ijkocc = 0
do jocc=1,G_nocc
   do iocc=1,jocc

      call outer_vecproduct(G_nbas,matC(:,iocc),matC(:,jocc),CC)

!$OMP PARALLEL PRIVATE(ithr,i,j,v2,val)

!$OMP SINGLE
      nthr = 1
!$ nthr = omp_get_num_threads()
!$OMP END SINGLE

!$OMP SECTIONS
!$OMP SECTION
      sum_0(1:G_nbas,0:nthr-1) = mpreal(0.d0)
!$OMP SECTION
      sum_1(1:G_nbas,0:nthr-1) = mpreal(0.d0)
!$OMP END SECTIONS

!$OMP DO SCHEDULE(DYNAMIC)
      do ij=1,(G_nbas*(G_nbas+1))/2
         ithr = 0
!$ ithr = omp_get_thread_num()
         j = int(sqrt(0.25d0 + 2.d0*(ij-1)) - safe_nbas)
         i = (ij-1) - (j*(j+1))/2
         v2 = i + j
         select case(mod(v2,2))
         case(0)
            v2 = v2/2+1
            val = CC(i,j)
            if(i/=j) val = val + CC(j,i)
            call accumulate_unit(val,v2,prod2H%unit(i,j),sum_0(:,ithr))
         case(1)
            v2 = (v2+1)/2
            val = CC(i,j) + CC(j,i)
            call accumulate_unit(val,v2,prod2H%unit(i,j),sum_1(:,ithr))
         end select
      enddo
!$OMP END DO

!$OMP SECTIONS
!$OMP SECTION
      do ithr=1,nthr-1
         call addition(G_nbas,sum_0(:,ithr),sum_0(:,0))
      enddo
!$OMP SECTION
      do ithr=1,nthr-1
         call addition(G_nbas,sum_1(:,ithr),sum_1(:,0))
      enddo
!$OMP END SECTIONS

!$OMP END PARALLEL

      do kocc=1,iocc
         ijkocc = ijkocc + 1
         n1 = G_nprim + 3*(G_nbas-1)
         n2 = G_nprim
         nC = 2*(G_nbas-1)
         call aux22C_create(aux22C(ijkocc),&
              n1,n2,nC,sum_0(:,0),sum_1(:,0),aux2C(kocc))
      enddo

   enddo
enddo

deallocate(sum_0,sum_1)
deallocate(CC)

n1 = 2
n2 = G_npair
allocate(pair_list(n1,n2))
i = -1
j = -1
ij = 0
do while(G_next_pair(i,j))
   ij = ij + 1
   pair_list(1,ij) = i
   pair_list(2,ij) = j
enddo

n1 = max(auxOP%n1_0,auxOP%n1_1)
n2 = G_nprim
allocate(ints(n1,0:n2),par_ints(0:n2),par_cont(0:n2))

vecP(1:G_npair,1:G_nocc,1:G_nocc) = mpreal(0.d0)

do kocc=1,G_nocc
   associate(kaux => aux2C(kocc))

     select case(parC(kocc))
     case(0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(v2)
        do i=0,G_nprim
           select case(mod(i,2))
           case(0)
              v2 = i/2+1
              par_ints(i) = 0
              call product2_se(auxOP%n1_0,auxOP%val_0,&
                   kaux%se_00(:,v2),kaux%val_00(:,v2),ints(:,i))
           case(1)
              v2 = (i+1)/2
              par_ints(i) = 1
              call product2_se(auxOP%n1_1,auxOP%val_1,&
                   kaux%se_11(:,v2),kaux%val_11(:,v2),ints(:,i))
           end select
        enddo
!$OMP END PARALLEL DO
     case(1)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(v2)
        do i=0,G_nprim
           select case(mod(i,2))
           case(0)
              v2 = i/2+1
              par_ints(i) = 1
              call product2_se(auxOP%n1_1,auxOP%val_1,&
                   kaux%se_10(:,v2),kaux%val_10(:,v2),ints(:,i))
           case(1)
              v2 = (i+1)/2
              par_ints(i) = 0
              call product2_se(auxOP%n1_0,auxOP%val_0,&
                   kaux%se_01(:,v2),kaux%val_01(:,v2),ints(:,i))
           end select
        enddo
!$OMP END PARALLEL DO
     end select

     do jocc=1,G_nocc
        do iocc=1,jocc
           if(kocc<=iocc) then
              ijkocc = kocc + ((iocc-1)*iocc)/2 + ((jocc-1)*jocc*(jocc+1))/6
           elseif(kocc<=jocc) then
              ijkocc = iocc + ((kocc-1)*kocc)/2 + ((jocc-1)*jocc*(jocc+1))/6
           else
              ijkocc = iocc + ((jocc-1)*jocc)/2 + ((kocc-1)*kocc*(kocc+1))/6
           endif
           associate(ijkaux => aux22C(ijkocc))

             select case(mod(parC(iocc)+parC(jocc)+parC(kocc),2))
             case(0)
                se_0   => ijkaux%se_00
                cont_0 => ijkaux%val_00
                se_1   => ijkaux%se_11
                cont_1 => ijkaux%val_11
                par_cont(0:G_nprim:2) = 0
                par_cont(1:G_nprim:2) = 1
             case(1)
                se_0   => ijkaux%se_10
                cont_0 => ijkaux%val_10
                se_1   => ijkaux%se_01
                cont_1 => ijkaux%val_01
                par_cont(0:G_nprim:2) = 1
                par_cont(1:G_nprim:2) = 0
             end select

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v1,val_1,val_2,val)
             do ij=1,G_npair

                i = pair_list(1,ij)
                j = pair_list(2,ij)

                if(par_cont(i)==par_ints(j)) then

                   select case(mod(j,2))
                   case(0)
                      v1 = j/2+1
                      val_1 = product1_se(se_0(:,v1),cont_0(:,v1),ints(:,i))
                   case(1)
                      v1 = (j+1)/2
                      val_1 = product1_se(se_1(:,v1),cont_1(:,v1),ints(:,i))
                   end select

                   select case(mod(i,2))
                   case(0)
                      v1 = i/2+1
                      val_2 = product1_se(se_0(:,v1),cont_0(:,v1),ints(:,j))
                   case(1)
                      v1 = (i+1)/2
                      val_2 = product1_se(se_1(:,v1),cont_1(:,v1),ints(:,j))
                   end select

                   val = G_gfac*(val_1 + val_2)

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

deallocate(ints,par_ints,par_cont)
deallocate(pair_list)

do ijkocc=1,(G_nocc*(G_nocc+1)*(G_nocc+2))/6
   call table4_free(aux22C(ijkocc))
enddo
deallocate(aux22C)
do iocc=1,G_nocc
   call table4_free(aux2C(iocc))
enddo
deallocate(aux2C)

call table2_free(auxOP)
call prod2H_free(prod2H)

if(IPRINT>=1) call timer('vecP',Tcpu,Twall)

end subroutine integrals_vecP


subroutine integrals_f12(matp_basicint,IPRINT)
implicit none
integer,intent(in) :: IPRINT
real(8) :: TCPU,Twall

type(mp_real),allocatable :: matp_binoms(:,:)
type(mp_real),intent(out) :: matp_basicint(0:)
type(mp_real) :: val,sumval
integer :: i,l,k

if(IPRINT>=1) call timer('START',Tcpu,Twall)

!!! Binomial section
allocate(matp_binoms(0:G_naux,0:G_naux))
do i = 1, G_naux
   matp_binoms(0,i) = mpreal(1.d0)
end do
do i = 1, G_naux
   do l = 1, G_naux
      matp_binoms(l,i) = matp_binoms(l-1,i)*(i-l+0.5d0)/l
      ! Calculates binomials of the form i-1/2 over l
   end do
end do

!!! Summing of binomials and basic integrals calculation
matp_basicint(0) = sqrt(2.d0*mppi())
do i = 1, G_naux
   matp_basicint(i) = matp_basicint(i-1) * (-1.d0) * 2.d0 * i
end do
do i = 1, G_naux
   sumval = 0d0
   do k = 0, i-1
      val = 1.0d0
      val = val * (k+0.5d0)/(i-k) * matp_binoms(i-k-1,i) * (-(2.0d0*G_alpha+1.d0))**(-k-1)
      sumval = sumval + val
   end do
   sumval = sumval + mpreal(-3.d0)**(-i-1)
   matp_basicint(i) = (-1.d0) * matp_basicint(i) * sumval
   ! i is the compund index (m+n)/2 for even m and even n. Multiply by -1 for odd values
   ! and multiply by 0 for other possibilities.
end do

matp_basicint(0) = matp_basicint(0)/3

if(IPRINT>=1) call timer('f12',Tcpu,Twall)
end subroutine integrals_f12

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

function product1(n,vec1,vec2) result(val)
implicit none
type(mp_real) :: val
integer,intent(in) :: n
type(mp_real),intent(in) :: vec1(:),vec2(:)
integer :: i

val = 0
do i=1,n
   val = val + vec1(i)*vec2(i)
enddo

end function product1

function product1_se(se,vec1,vec2) result(val)
implicit none
type(mp_real) :: val
integer,intent(in) :: se(:)
type(mp_real),intent(in) :: vec1(:),vec2(:)
integer :: i

val = 0
do i=se(1),se(2)
   val = val + vec1(i)*vec2(i)
enddo

end function product1_se

function product1_unit(v,unit,vec) result(val)
implicit none
type(mp_real) :: val
integer,intent(in) :: v
type(prodUnitData),intent(in) :: unit
type(mp_real),intent(in) :: vec(:)
integer :: off,i

off = v - unit%nu

val = unit%val(1)*vec(off)
do i=1,unit%nu
   val = val + unit%val(1+i)*vec(off+i)
enddo

end function product1_unit

subroutine product2_se(n,mat,se,cont,vec)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: mat(:,:)
integer,intent(in) :: se(:)
type(mp_real),intent(in) :: cont(:)
type(mp_real),intent(out) :: vec(:)
integer :: i1,i2

vec(1:n) = mpreal(0.d0)

do i2=se(1),se(2)
   do i1=1,n
      vec(i1) = vec(i1) + mat(i1,i2)*cont(i2)
   enddo
enddo

end subroutine product2_se

subroutine product2_unit(n,mat,v,unit,vec)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: mat(:,:)
integer,intent(in) :: v
type(prodUnitData),intent(in) :: unit
type(mp_real),intent(out) :: vec(:)
integer :: off2,i1,i2

off2 = v - unit%nu

do i1=1,n
   vec(i1) = mat(i1,off2)*unit%val(1)
enddo
do i2=1,unit%nu
   do i1=1,n
      vec(i1) = vec(i1) + mat(i1,off2+i2)*unit%val(1+i2)
   enddo
enddo

end subroutine product2_unit

subroutine product2_se_unit(n,se,mat,v,unit,vec)
implicit none
integer,intent(in) :: n
integer,intent(in) :: se(:,:)
type(mp_real),intent(in) :: mat(:,:)
integer,intent(in) :: v
type(prodUnitData),intent(in) :: unit
type(mp_real),intent(out) :: vec(:)
integer :: off2,i1,i2

vec(1:n) = mpreal(0.d0)

off2 = v - unit%nu

do i2=0,unit%nu
   do i1=se(1,off2+i2),se(2,off2+i2)
      vec(i1) = vec(i1) + mat(i1,off2+i2)*unit%val(1+i2)
   enddo
enddo

end subroutine product2_se_unit

subroutine accumulate(alpha,n,vec,acc)
implicit none
type(mp_real),intent(in) :: alpha
integer,intent(in) :: n
type(mp_real),intent(in) :: vec(:)
type(mp_real),intent(inout) :: acc(:)
integer :: i

do i=1,n
   acc(i) = acc(i) + alpha*vec(i)
enddo

end subroutine accumulate

subroutine accumulate_unit(alpha,v,unit,acc)
implicit none
type(mp_real),intent(in) :: alpha
integer,intent(in) :: v
type(prodUnitData),intent(in) :: unit
type(mp_real),intent(inout) :: acc(:)
integer :: off,i

off = v - unit%nu

do i=0,unit%nu
   acc(off+i) = acc(off+i) + alpha*unit%val(1+i)
enddo

end subroutine accumulate_unit

subroutine addition(n,vec,acc)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: vec(:)
type(mp_real),intent(inout) :: acc(:)
integer :: i

do i=1,n
   acc(i) = acc(i) + vec(i)
enddo

end subroutine addition

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine aux22_create(aux22,n1,n2)
implicit none
type(table2Data) :: aux22
integer,intent(in) :: n1,n2
integer :: i1,i2,i12,j
type(mp_real),allocatable :: dfact(:)
type(mp_real) :: tmp,val

aux22%n1_0 = n1/2+1
aux22%n2_0 = n2/2+1

aux22%n1_1 = (n1+1)/2
aux22%n2_1 = (n2+1)/2

allocate(&
     aux22%val_0(aux22%n1_0,aux22%n2_0),&
     aux22%val_1(aux22%n1_1,aux22%n2_1))

allocate(dfact(0:(n1+n2)/2))

dfact(0) = sqrt(0.5d0*mppi())
do j=1,(n1+n2)/2
   dfact(j) = dfact(j-1)*(2*j-1)
enddo

!$OMP PARALLEL PRIVATE(i12,tmp,val)

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
do i2=0,n2,2
   do i1=0,n1,2
      i12 = (i1+i2)/2
      tmp = (1 - 2*mod(i12,2))
      val = dfact(i12)*tmp
      do j=0,min(i1,i2)-1
         tmp = -(tmp/(j+1))*(2*(i1-j)*(i2-j))
         val = val + dfact(i12-j-1)*tmp
      enddo
      aux22%val_0(i1/2+1,i2/2+1) = val
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
do i2=1,n2,2
   do i1=1,n1,2
      i12 = (i1+i2)/2
      tmp = (1 - 2*mod(i12,2))
      val = dfact(i12)*tmp
      do j=0,min(i1,i2)-1
         tmp = -(tmp/(j+1))*(2*(i1-j)*(i2-j))
         val = val + dfact(i12-j-1)*tmp
      enddo
      aux22%val_1((i1+1)/2,(i2+1)/2) = val
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

deallocate(dfact)

end subroutine aux22_create

subroutine aux22F2_create(aux22F2,c0,c2,aux)
implicit none
type(table2Data) :: aux22F2
type(mp_real),intent(in) :: c0,c2
type(table2Data),intent(in) :: aux
integer :: i1,i2

aux22F2%n1_0 = aux%n1_0
aux22F2%n2_0 = aux%n2_0

aux22F2%n1_1 = aux%n1_1
aux22F2%n2_1 = aux%n2_1

allocate(&
     aux22F2%val_0(aux22F2%n1_0,aux22F2%n2_0),&
     aux22F2%val_1(aux22F2%n1_1,aux22F2%n2_1))

!$OMP PARALLEL

!$OMP SINGLE
aux22F2%val_0(1,1) = c0*aux%val_0(1,1)
!$OMP END SINGLE NOWAIT

!$OMP DO
do i1=2,aux22F2%n1_0
   aux22F2%val_0(i1,1) = c0*aux%val_0(i1,1) &
        + c2*aux%val_0(i1-1,1)
enddo
!$OMP END DO NOWAIT

!$OMP DO
do i2=2,aux22F2%n2_0
   aux22F2%val_0(1,i2) = c0*aux%val_0(1,i2) &
        + c2*aux%val_0(1,i2-1)
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do i2=2,aux22F2%n2_0
   do i1=2,aux22F2%n1_0
      aux22F2%val_0(i1,i2) = c0*aux%val_0(i1,i2) &
           + c2*(aux%val_0(i1-1,i2) + aux%val_0(i1,i2-1) &
           - 2*aux%val_1(i1-1,i2-1))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP SINGLE
aux22F2%val_1(1,1) = c0*aux%val_1(1,1) &
     - 2*c2*aux%val_0(1,1)
!$OMP END SINGLE NOWAIT

!$OMP DO
do i1=2,aux22F2%n1_1
   aux22F2%val_1(i1,1) = c0*aux%val_1(i1,1) &
        + c2*(aux%val_1(i1-1,1) - 2*aux%val_0(i1,1))
enddo
!$OMP END DO NOWAIT

!$OMP DO
do i2=2,aux22F2%n2_1
   aux22F2%val_1(1,i2) = c0*aux%val_1(1,i2) &
        + c2*(aux%val_1(1,i2-1) - 2*aux%val_0(1,i2))
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do i2=2,aux22F2%n2_1
   do i1=2,aux22F2%n1_1
      aux22F2%val_1(i1,i2) = c0*aux%val_1(i1,i2) &
           + c2*(aux%val_1(i1-1,i2) + aux%val_1(i1,i2-1) &
           - 2*aux%val_0(i1,i2))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

end subroutine aux22F2_create

subroutine aux22DF2_create(aux22DF2,c0,c2,c4,aux)
implicit none
type(table2Data) :: aux22DF2
type(mp_real),intent(in) :: c0,c2,c4
type(table2Data),intent(in) :: aux
integer :: i1,i2

aux22DF2%n1_0 = aux%n1_0
aux22DF2%n2_0 = aux%n2_0

aux22DF2%n1_1 = aux%n1_1
aux22DF2%n2_1 = aux%n2_1

allocate(&
     aux22DF2%val_0(aux22DF2%n1_0,aux22DF2%n2_0),&
     aux22DF2%val_1(aux22DF2%n1_1,aux22DF2%n2_1))

!$OMP PARALLEL

!$OMP SECTIONS
!$OMP SECTION
aux22DF2%val_0(1,1) = c0*aux%val_0(1,1)
!$OMP SECTION
aux22DF2%val_0(2,1) = c0*aux%val_0(2,1) &
     + c2*aux%val_0(1,1)
!$OMP SECTION
aux22DF2%val_0(1,2) = c0*aux%val_0(1,2) &
     + c2*aux%val_0(1,1)
!$OMP SECTION
aux22DF2%val_0(2,2) = c0*aux%val_0(2,2) &
     + c2*(aux%val_0(1,2) + aux%val_0(2,1) - 2*aux%val_1(1,1)) &
     + 6*c4*aux%val_0(1,1)
!$OMP END SECTIONS NOWAIT

!$OMP DO
do i1=3,aux22DF2%n1_0
   aux22DF2%val_0(i1,1) = c0*aux%val_0(i1,1) &
        + c2*aux%val_0(i1-1,1) &
        + c4*aux%val_0(i1-2,1)
enddo
!$OMP END DO NOWAIT

!$OMP DO
do i1=3,aux22DF2%n1_0
   aux22DF2%val_0(i1,2) = c0*aux%val_0(i1,2) &
        + c2*(aux%val_0(i1-1,2) + aux%val_0(i1,1) - 2*aux%val_1(i1-1,1))&
        + c4*(aux%val_0(i1-2,2) + 6*aux%val_0(i1-1,1) - 4*aux%val_1(i1-2,1))
enddo
!$OMP END DO NOWAIT

!$OMP DO
do i2=3,aux22DF2%n2_0
   aux22DF2%val_0(1,i2) = c0*aux%val_0(1,i2) &
        + c2*aux%val_0(1,i2-1) &
        + c4*aux%val_0(1,i2-2)
enddo
!$OMP END DO NOWAIT

!$OMP DO
do i2=3,aux22DF2%n2_0
   aux22DF2%val_0(2,i2) = c0*aux%val_0(2,i2) &
        + c2*(aux%val_0(2,i2-1) + aux%val_0(1,i2) - 2*aux%val_1(1,i2-1))&
        + c4*(aux%val_0(2,i2-2) + 6*aux%val_0(1,i2-1) - 4*aux%val_1(1,i2-2))
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do i2=3,aux22DF2%n2_0
   do i1=3,aux22DF2%n1_0
      aux22DF2%val_0(i1,i2) = c0*aux%val_0(i1,i2) &
           + c2*(aux%val_0(i1-1,i2) + aux%val_0(i1,i2-1) &
           - 2*aux%val_1(i1-1,i2-1)) &
           + c4*(aux%val_0(i1-2,i2) + aux%val_0(i1,i2-2) &
           + 6*aux%val_0(i1-1,i2-1) &
           - 4*aux%val_1(i1-2,i2-1) - 4*aux%val_1(i1-1,i2-2))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP SECTIONS
!$OMP SECTION
aux22DF2%val_1(1,1) = c0*aux%val_1(1,1) &
     - 2*c2*aux%val_0(1,1)
!$OMP SECTION
aux22DF2%val_1(2,1) = c0*aux%val_1(2,1) &
     + c2*(aux%val_1(1,1) - 2*aux%val_0(2,1)) &
     - 4*c4*aux%val_0(1,1)
!$OMP SECTION
aux22DF2%val_1(1,2) = c0*aux%val_1(1,2) &
     + c2*(aux%val_1(1,1) - 2*aux%val_0(1,2)) &
     - 4*c4*aux%val_0(1,1)
!$OMP SECTION
aux22DF2%val_1(2,2) = c0*aux%val_1(2,2) &
     + c2*(aux%val_1(1,2) + aux%val_1(2,1) - 2*aux%val_0(2,2)) &
     + c4*(6*aux%val_1(1,1) - 4*aux%val_0(1,2) - 4*aux%val_0(2,1))
!$OMP END SECTIONS NOWAIT

!$OMP DO
do i1=3,aux22DF2%n1_1
   aux22DF2%val_1(i1,1) = c0*aux%val_1(i1,1) &
        + c2*(aux%val_1(i1-1,1) - 2*aux%val_0(i1,1)) &
        + c4*(aux%val_1(i1-2,1) - 4*aux%val_0(i1-1,1))
enddo
!$OMP END DO NOWAIT

!$OMP DO
do i1=3,aux22DF2%n1_1
   aux22DF2%val_1(i1,2) = c0*aux%val_1(i1,2) &
        + c2*(aux%val_1(i1-1,2) + aux%val_1(i1,1) - 2*aux%val_0(i1,2)) &
        + c4*(aux%val_1(i1-2,2) + 6*aux%val_1(i1-1,1) &
        - 4*aux%val_0(i1-1,2) - 4*aux%val_0(i1,1))
enddo
!$OMP END DO NOWAIT

!$OMP DO
do i2=3,aux22DF2%n2_1
   aux22DF2%val_1(1,i2) = c0*aux%val_1(1,i2) &
        + c2*(aux%val_1(1,i2-1) - 2*aux%val_0(1,i2)) &
        + c4*(aux%val_1(1,i2-2) - 4*aux%val_0(1,i2-1))
enddo
!$OMP END DO NOWAIT

!$OMP DO
do i2=3,aux22DF2%n2_1
   aux22DF2%val_1(2,i2) = c0*aux%val_1(2,i2) &
        + c2*(aux%val_1(2,i2-1) + aux%val_1(1,i2) - 2*aux%val_0(2,i2)) &
        + c4*(aux%val_1(2,i2-2) + 6*aux%val_1(1,i2-1) &
        - 4*aux%val_0(2,i2-1) - 4*aux%val_0(1,i2))
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do i2=3,aux22DF2%n2_1
   do i1=3,aux22DF2%n1_1
      aux22DF2%val_1(i1,i2) = c0*aux%val_1(i1,i2) &
           + c2*(aux%val_1(i1-1,i2) + aux%val_1(i1,i2-1) &
           - 2*aux%val_0(i1,i2)) &
           + c4*(aux%val_1(i1-2,i2) + aux%val_1(i1,i2-2) &
           + 6*aux%val_1(i1-1,i2-1) &
           - 4*aux%val_0(i1-1,i2) - 4*aux%val_0(i1,i2-1))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

end subroutine aux22DF2_create

subroutine aux42F2_create(aux42F2,c0,c2,aux)
implicit none
type(table2Data) :: aux42F2
type(mp_real),intent(in) :: c0,c2
type(table2Data),intent(in) :: aux
integer :: i1,i2

aux42F2%n1_0 = aux%n1_0
aux42F2%n2_0 = aux%n2_0

aux42F2%n1_1 = aux%n1_1
aux42F2%n2_1 = aux%n2_1

allocate(&
     aux42F2%val_0(aux42F2%n1_0,aux42F2%n2_0),&
     aux42F2%val_1(aux42F2%n1_1,aux42F2%n2_1))

!$OMP PARALLEL

!$OMP SINGLE
aux42F2%val_0(1,1) = c0*aux%val_0(1,1)
!$OMP END SINGLE NOWAIT

!$OMP DO
do i1=2,aux42F2%n1_0
   aux42F2%val_0(i1,1) = c0*aux%val_0(i1,1) &
        + c2*aux%val_0(i1-1,1)
enddo
!$OMP END DO NOWAIT

!$OMP DO
do i2=2,aux42F2%n2_0
   aux42F2%val_0(1,i2) = c0*aux%val_0(1,i2) &
        + 4*c2*aux%val_0(1,i2-1)
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do i2=2,aux42F2%n2_0
   do i1=2,aux42F2%n1_0
      aux42F2%val_0(i1,i2) = c0*aux%val_0(i1,i2) &
           + c2*(aux%val_0(i1-1,i2) + 4*aux%val_0(i1,i2-1) &
           - 4*aux%val_1(i1-1,i2-1))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP SINGLE
aux42F2%val_1(1,1) = c0*aux%val_1(1,1) &
     - 4*c2*aux%val_0(1,1)
!$OMP END SINGLE NOWAIT

!$OMP DO
do i1=2,aux42F2%n1_1
   aux42F2%val_1(i1,1) = c0*aux%val_1(i1,1) &
        + c2*(aux%val_1(i1-1,1) - 4*aux%val_0(i1,1))
enddo
!$OMP END DO NOWAIT

!$OMP DO
do i2=2,aux42F2%n2_1
   aux42F2%val_1(1,i2) = c0*aux%val_1(1,i2) &
        + c2*(4*aux%val_1(1,i2-1) - 4*aux%val_0(1,i2))
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do i2=2,aux42F2%n2_1
   do i1=2,aux42F2%n1_1
      aux42F2%val_1(i1,i2) = c0*aux%val_1(i1,i2) &
           + c2*(aux%val_1(i1-1,i2) + 4*aux%val_1(i1,i2-1) &
           - 4*aux%val_0(i1,i2))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

end subroutine aux42F2_create

subroutine auxXXF1_create(auxXXF1,int0,c0,mult,aux)
implicit none
type(table2Data) :: auxXXF1
type(mp_real),intent(in) :: int0,c0,mult
type(table2Data),intent(in) :: aux
integer :: nc,n,k,i1,i2
type(mp_real),allocatable :: coeff(:,:),tmp(:)
type(mp_real) :: val

auxXXF1%n1_0 = aux%n1_0
auxXXF1%n2_0 = aux%n2_0

auxXXF1%n1_1 = aux%n1_1
auxXXF1%n2_1 = aux%n2_1

allocate(&
     auxXXF1%val_0(auxXXF1%n1_0,auxXXF1%n2_0),&
     auxXXF1%val_1(auxXXF1%n1_1,auxXXF1%n2_1))

nc = max(1,auxXXF1%n1_0+auxXXF1%n2_0,auxXXF1%n1_0+auxXXF1%n2_0)
allocate(coeff(0:2*nc,0:nc),tmp(0:nc))

!$OMP PARALLEL PRIVATE(val)

!$OMP SECTIONS
!$OMP SECTION
tmp(0) = int0
do n=1,nc
   tmp(n) = tmp(n-1)*c0/(2*n-1)
enddo
!$OMP SECTION
coeff = mpreal(0.d0)
coeff(0,0) = 1
do n=1,nc
   coeff(0,n) = 1
   coeff(1,n) = 2*n
   do k=2,2*n
      coeff(k,n) = coeff(k-2,n-1) + 2*coeff(k-1,n-1) + coeff(k,n-1)
   enddo
enddo
!$OMP END SECTIONS

!$OMP DO SCHEDULE(DYNAMIC)
do n=0,nc
   val = tmp(n)
   coeff(0,n) = val*coeff(0,n)
   do k=1,2*n
      val = mult*val
      coeff(k,n) = val*coeff(k,n)
   enddo
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
do i2=1,auxXXF1%n2_0
   do i1=1,auxXXF1%n1_0
      auxXXF1%val_0(i1,i2) = coeff(0,0)*aux%val_0(i1,i2)
      do n=1,i1+i2-2
         val = 0
         do k=n-min(n,i1-1),min(n,i2-1)
            val = val + coeff(2*k,n)*aux%val_0(i1-(n-k),i2-k)
         enddo
         auxXXF1%val_0(i1,i2) = auxXXF1%val_0(i1,i2) + val
         val = 0
         do k=n-min(n,i1-1),min(n,i2-1)-1
            val = val + coeff(2*k+1,n)*aux%val_1(i1-(n-k),i2-k-1)
         enddo
         auxXXF1%val_0(i1,i2) = auxXXF1%val_0(i1,i2) + val
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
do i2=1,auxXXF1%n2_1
   do i1=1,auxXXF1%n1_1
      auxXXF1%val_1(i1,i2) = coeff(0,0)*aux%val_1(i1,i2)
      do n=1,i1+i2-1
         val = 0
         do k=n-min(n,i1-1),min(n,i2-1)
            val = val + coeff(2*k,n)*aux%val_1(i1-(n-k),i2-k)
         enddo
         auxXXF1%val_1(i1,i2) = auxXXF1%val_1(i1,i2) + val
         val = 0
         do k=n-min(n,i1),min(n,i2)-1
            val = val + coeff(2*k+1,n)*aux%val_0(i1-(n-k)+1,i2-k)
         enddo
         auxXXF1%val_1(i1,i2) = auxXXF1%val_1(i1,i2) + val
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

deallocate(coeff,tmp)

end subroutine auxXXF1_create

subroutine table2_easy(table2,basicF12,n1,n2)
  implicit none
  type(table2Data) :: table2
  type(mp_real) :: basicF12(0:)
  integer :: n1,n2
  integer :: i,j,ij

  table2%n1_0 = max(2,n1/2+1)
  table2%n2_0 = max(2,n2/2+1)

  table2%n1_1 = max(2,(n1+1)/2)
  table2%n2_1 = max(2,(n2+1)/2)

  allocate(&
       table2%val_0(table2%n1_0,table2%n2_0),&
       table2%val_1(table2%n1_1,table2%n2_1))

  do i = 1, table2%n1_0
     do j = 1, table2%n2_0
        ij=(2*(i-1)+2*(j-1))/2
        table2%val_0(i,j) = basicF12(ij)
     end do
  end do

  do i = 1, table2%n1_1
     do j = 1, table2%n2_1
        ij=(2*(i-1)+2*(j-1))/2+1
        table2%val_1(i,j) = -1.d0*basicF12(ij)
     end do
  end do

end subroutine table2_easy

subroutine table2_create(table2,n1,n2,a1,a2,a12)
implicit none
type(table2Data) :: table2
integer,intent(in) :: n1,n2
type(mp_real),intent(in) :: a1,a2,a12
type(mp_real) :: a12_2
integer :: i1,i2,j
type(mp_real),allocatable :: tmp1(:),tmp2(:),tmp12(:)
type(mp_real) :: val

a12_2 = a12**2

table2%n1_0 = max(2,n1/2+1)
table2%n2_0 = max(2,n2/2+1)

table2%n1_1 = max(2,(n1+1)/2)
table2%n2_1 = max(2,(n2+1)/2)

allocate(&
     table2%val_0(table2%n1_0,table2%n2_0),&
     table2%val_1(table2%n1_1,table2%n2_1))

allocate(&
     tmp1(table2%n1_0),&
     tmp2(table2%n2_0),&
     tmp12(min(table2%n1_0,table2%n2_0)))

!$OMP PARALLEL PRIVATE(val)

!$OMP SECTIONS
!$OMP SECTION
tmp1(1) = 1
do i1=1,table2%n1_0-1
   tmp1(i1+1) = tmp1(i1)*(a1/i1)
enddo
!$OMP SECTION
tmp2(1) = 1
do i2=1,table2%n2_0-1
   tmp2(i2+1) = tmp2(i2)*(a2/i2)
enddo
!$OMP SECTION
tmp12(1) = 1
do j=1,min(table2%n1_0,table2%n2_0)-1
   tmp12(j+1) = tmp12(j)*(a12_2/((2*j-1)*(2*j)))
enddo
!$OMP END SECTIONS

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
do i2=1,table2%n2_0
   do i1=1,table2%n1_0

      val = 0
      do j=0,min(i1,i2)-1
         val = val + tmp1(i1-j)*tmp2(i2-j)*tmp12(j+1)
      enddo

      table2%val_0(i1,i2) = val

   enddo
enddo
!$OMP END DO

!$OMP DO
do j=1,min(table2%n1_0,table2%n2_0)
   tmp12(j) = tmp12(j)*(a12/(2*j-1))
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
do i2=1,table2%n2_1
   do i1=1,table2%n1_1

      val = 0
      do j=0,min(i1,i2)-1
         val = val + tmp1(i1-j)*tmp2(i2-j)*tmp12(j+1)
      enddo

      table2%val_1(i1,i2) = val

   enddo
enddo
!$OMP END DO

!$OMP END PARALLEL

deallocate(tmp1,tmp2,tmp12)

end subroutine table2_create

subroutine table2_fact(table2)
implicit none
type(table2Data) :: table2
type(mp_real),allocatable :: fact(:)
integer :: i1,i2,j

allocate(fact(max(table2%n1_0,table2%n2_0)))

fact(1) = 1
do j=1,max(table2%n1_0,table2%n2_0)-1
   fact(j+1) = fact(j)*((2*j)*(2*j-1))
enddo

!$OMP PARALLEL

!$OMP DO COLLAPSE(2)
do i2=1,table2%n2_0
   do i1=1,table2%n1_0
      table2%val_0(i1,i2) = table2%val_0(i1,i2)*fact(i1)*fact(i2)
   enddo
enddo
!$OMP END DO

!$OMP DO
do j=1,max(table2%n1_0,table2%n2_0)
   fact(j) = fact(j)*(2*j-1)
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
do i2=1,table2%n2_1
   do i1=1,table2%n1_1
      table2%val_1(i1,i2) = table2%val_1(i1,i2)*fact(i1)*fact(i2)
   enddo
enddo
!$OMP END DO

!$OMP END PARALLEL

deallocate(fact)

end subroutine table2_fact

subroutine table2_free(table2)
implicit none
type(table2Data) :: table2

deallocate(table2%val_0,table2%val_1)

end subroutine table2_free

subroutine aux2C_create(aux2C,n1,n2,nC,vec,prod2H)
implicit none
type(table4Data) :: aux2C
integer,intent(in) :: n1,n2,nC
type(mp_real),intent(in) :: vec(0:)
type(prod2HData),intent(in) :: prod2H
integer :: nC_0,nC_1
integer :: i2,iC,ip1,ip2

nC_0 = nC - mod(nC,2)
nC_1 = nC - mod(nC+1,2)

aux2C%n1_0 = n1/2+1
aux2C%n2_0 = n2/2+1

aux2C%n1_1 = (n1+1)/2
aux2C%n2_1 = (n2+1)/2

allocate(&
     aux2C%se_00(2,aux2C%n2_0),&
     aux2C%se_10(2,aux2C%n2_0),&
     aux2C%se_01(2,aux2C%n2_1),&
     aux2C%se_11(2,aux2C%n2_1))

allocate(&
     aux2C%val_00(aux2C%n1_0,aux2C%n2_0),&
     aux2C%val_10(aux2C%n1_1,aux2C%n2_0),&
     aux2C%val_01(aux2C%n1_0,aux2C%n2_1),&
     aux2C%val_11(aux2C%n1_1,aux2C%n2_1))

!$OMP PARALLEL PRIVATE(ip1,ip2)

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   aux2C%se_00(1,ip2) = max(i2-nC_0,0)/2+1
   aux2C%se_00(2,ip2) = (i2+nC_0)/2+1
   aux2C%val_00(aux2C%se_00(1,ip2):aux2C%se_00(2,ip2),ip2) = mpreal(0.d0)
   do iC=0,nC,2
      ip1 = (i2+iC)/2+1
      call accumulate_unit(vec(iC),ip1,prod2H%unit(iC,i2),aux2C%val_00(:,ip2))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   aux2C%se_10(1,ip2) = (max(i2-nC_1,1)+1)/2
   aux2C%se_10(2,ip2) = ((i2+nC_1)+1)/2
   aux2C%val_10(aux2C%se_10(1,ip2):aux2C%se_10(2,ip2),ip2) = mpreal(0.d0)
   do iC=1,nC,2
      ip1 = ((i2+iC)+1)/2
      call accumulate_unit(vec(iC),ip1,prod2H%unit(iC,i2),aux2C%val_10(:,ip2))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=1,n2,2
   ip2 = (i2+1)/2
   aux2C%se_01(1,ip2) = max(i2-nC_1,0)/2+1
   aux2C%se_01(2,ip2) = (i2+nC_1)/2+1
   aux2C%val_01(aux2C%se_01(1,ip2):aux2C%se_01(2,ip2),ip2) = mpreal(0.d0)
   do iC=1,nC,2
      ip1 = (i2+iC)/2+1
      call accumulate_unit(vec(iC),ip1,prod2H%unit(iC,i2),aux2C%val_01(:,ip2))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=1,n2,2
   ip2 = (i2+1)/2
   aux2C%se_11(1,ip2) = (max(i2-nC_0,1)+1)/2
   aux2C%se_11(2,ip2) = ((i2+nC_0)+1)/2
   aux2C%val_11(aux2C%se_11(1,ip2):aux2C%se_11(2,ip2),ip2) = mpreal(0.d0)
   do iC=0,nC,2
      ip1 = ((i2+iC)+1)/2
      call accumulate_unit(vec(iC),ip1,prod2H%unit(iC,i2),aux2C%val_11(:,ip2))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

end subroutine aux2C_create

subroutine aux22C_create(aux22C,n1,n2,nC,sum_0,sum_1,aux)
implicit none
type(table4Data) :: aux22C
integer,intent(in) :: n1,n2,nC
type(mp_real),intent(in) :: sum_0(:),sum_1(:)
type(table4Data),intent(in) :: aux
integer :: nC_0,nC_1
integer :: i2,iC1,iC2,j,ip1,ip2
type(mp_real) :: tmp

nC_0 = nC - mod(nC,2)
nC_1 = nC - mod(nC+1,2)

aux22C%n1_0 = n1/2+1
aux22C%n2_0 = n2/2+1

aux22C%n1_1 = (n1+1)/2
aux22C%n2_1 = (n2+1)/2

allocate(&
     aux22C%se_00(2,aux22C%n2_0),&
     aux22C%se_10(2,aux22C%n2_0),&
     aux22C%se_01(2,aux22C%n2_1),&
     aux22C%se_11(2,aux22C%n2_1))

allocate(&
     aux22C%val_00(aux22C%n1_0,aux22C%n2_0),&
     aux22C%val_10(aux22C%n1_1,aux22C%n2_0),&
     aux22C%val_01(aux22C%n1_0,aux22C%n2_1),&
     aux22C%val_11(aux22C%n1_1,aux22C%n2_1))

!$OMP PARALLEL PRIVATE(ip1,ip2,tmp)

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   aux22C%se_00(1,ip2) = 1
   aux22C%se_00(2,ip2) = &
        max(nC_0 + 2*(aux%se_00(2,ip2)-1),nC_1 + 2*aux%se_10(2,ip2)-1)/2+1
   aux22C%val_00(aux22C%se_00(1,ip2):aux22C%se_00(2,ip2),ip2) = mpreal(0.d0)
   do iC2=2*(aux%se_00(1,ip2)-1),2*(aux%se_00(2,ip2)-1),2
      do iC1=0,nC,2
         ip1 = (iC1+iC2)/2+1
         tmp = sum_0(iC1/2+1)*aux%val_00(iC2/2+1,ip2)
         aux22C%val_00(ip1,ip2) = aux22C%val_00(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            aux22C%val_00(ip1-j-1,ip2) = aux22C%val_00(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
   do iC2=2*aux%se_10(1,ip2)-1,2*aux%se_10(2,ip2)-1,2
      do iC1=1,nC,2
         ip1 = (iC1+iC2)/2+1
         tmp = sum_1((iC1+1)/2)*aux%val_10((iC2+1)/2,ip2)
         aux22C%val_00(ip1,ip2) = aux22C%val_00(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            aux22C%val_00(ip1-j-1,ip2) = aux22C%val_00(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   aux22C%se_10(1,ip2) = 1
   aux22C%se_10(2,ip2) = &
        (max(nC_1 + 2*(aux%se_00(2,ip2)-1),nC_0 + 2*aux%se_10(2,ip2)-1)+1)/2
   aux22C%val_10(aux22C%se_10(1,ip2):aux22C%se_10(2,ip2),ip2) = mpreal(0.d0)
   do iC2=2*(aux%se_00(1,ip2)-1),2*(aux%se_00(2,ip2)-1),2
      do iC1=1,nC,2
         ip1 = (iC1+iC2+1)/2
         tmp = sum_1((iC1+1)/2)*aux%val_00(iC2/2+1,ip2)
         aux22C%val_10(ip1,ip2) = aux22C%val_10(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            aux22C%val_10(ip1-j-1,ip2) = aux22C%val_10(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
   do iC2=2*aux%se_10(1,ip2)-1,2*aux%se_10(2,ip2)-1,2
      do iC1=0,nC,2
         ip1 = (iC1+iC2+1)/2
         tmp = sum_0(iC1/2+1)*aux%val_10((iC2+1)/2,ip2)
         aux22C%val_10(ip1,ip2) = aux22C%val_10(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            aux22C%val_10(ip1-j-1,ip2) = aux22C%val_10(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=1,n2,2
   ip2 = (i2+1)/2
   aux22C%se_01(1,ip2) = 1
   aux22C%se_01(2,ip2) = &
        max(nC_0 + 2*(aux%se_01(2,ip2)-1),nC_1 + 2*aux%se_11(2,ip2)-1)/2+1
   aux22C%val_01(aux22C%se_01(1,ip2):aux22C%se_01(2,ip2),ip2) = mpreal(0.d0)
   do iC2=2*(aux%se_01(1,ip2)-1),2*(aux%se_01(2,ip2)-1),2
      do iC1=0,nC,2
         ip1 = (iC1+iC2)/2+1
         tmp = sum_0(iC1/2+1)*aux%val_01(iC2/2+1,ip2)
         aux22C%val_01(ip1,ip2) = aux22C%val_01(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            aux22C%val_01(ip1-j-1,ip2) = aux22C%val_01(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
   do iC2=2*aux%se_11(1,ip2)-1,2*aux%se_11(2,ip2)-1,2
      do iC1=1,nC,2
         ip1 = (iC1+iC2)/2+1
         tmp = sum_1((iC1+1)/2)*aux%val_11((iC2+1)/2,ip2)
         aux22C%val_01(ip1,ip2) = aux22C%val_01(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            aux22C%val_01(ip1-j-1,ip2) = aux22C%val_01(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=1,n2,2
   ip2 = (i2+1)/2
   aux22C%se_11(1,ip2) = 1
   aux22C%se_11(2,ip2) = &
        (max(nC_1 + 2*(aux%se_01(2,ip2)-1),nC_0 + 2*aux%se_11(2,ip2)-1)+1)/2
   aux22C%val_11(aux22C%se_11(1,ip2):aux22C%se_11(2,ip2),ip2) = mpreal(0.d0)
   do iC2=2*(aux%se_01(1,ip2)-1),2*(aux%se_01(2,ip2)-1),2
      do iC1=1,nC,2
         ip1 = (iC1+iC2+1)/2
         tmp = sum_1((iC1+1)/2)*aux%val_01(iC2/2+1,ip2)
         aux22C%val_11(ip1,ip2) = aux22C%val_11(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            aux22C%val_11(ip1-j-1,ip2) = aux22C%val_11(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
   do iC2=2*aux%se_11(1,ip2)-1,2*aux%se_11(2,ip2)-1,2
      do iC1=0,nC,2
         ip1 = (iC1+iC2+1)/2
         tmp = sum_0(iC1/2+1)*aux%val_11((iC2+1)/2,ip2)
         aux22C%val_11(ip1,ip2) = aux22C%val_11(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            aux22C%val_11(ip1-j-1,ip2) = aux22C%val_11(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

end subroutine aux22C_create

subroutine aux4C_create(aux4C,n1,n2,nC,sum_0,sum_1)
implicit none
type(table4Data) :: aux4C
integer,intent(in) :: n1,n2,nC
type(mp_real),intent(in) :: sum_0(:),sum_1(:)
integer :: nC_0,nC_1
integer :: i2,iC,j,ip1,ip2
type(mp_real) :: tmp

nC_0 = nC - mod(nC,2)
nC_1 = nC - mod(nC+1,2)

aux4C%n1_0 = n1/2+1
aux4C%n2_0 = n2/2+1

aux4C%n1_1 = (n1+1)/2
aux4C%n2_1 = (n2+1)/2

allocate(&
     aux4C%se_00(2,aux4C%n2_0),&
     aux4C%se_10(2,aux4C%n2_0),&
     aux4C%se_01(2,aux4C%n2_1),&
     aux4C%se_11(2,aux4C%n2_1))

allocate(&
     aux4C%val_00(aux4C%n1_0,aux4C%n2_0),&
     aux4C%val_10(aux4C%n1_1,aux4C%n2_0),&
     aux4C%val_01(aux4C%n1_0,aux4C%n2_1),&
     aux4C%val_11(aux4C%n1_1,aux4C%n2_1))

!$OMP PARALLEL PRIVATE(ip1,ip2,tmp)

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   aux4C%se_00(1,ip2) = max(i2-nC_0,0)/2+1
   aux4C%se_00(2,ip2) = (i2+nC_0)/2+1
   aux4C%val_00(aux4C%se_00(1,ip2):aux4C%se_00(2,ip2),ip2) = mpreal(0.d0)
   do iC=0,nC,2
      ip1 = (i2+iC)/2+1
      tmp = sum_0(iC/2+1)
      aux4C%val_00(ip1,ip2) = aux4C%val_00(ip1,ip2) + tmp
      do j=0,min(i2,iC)-1
         tmp = (tmp/(j+1))*(2*(i2-j)*(iC-j))
         aux4C%val_00(ip1-j-1,ip2) = aux4C%val_00(ip1-j-1,ip2) + tmp
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   aux4C%se_10(1,ip2) = (max(i2-nC_1,1)+1)/2
   aux4C%se_10(2,ip2) = ((i2+nC_1)+1)/2
   aux4C%val_10(aux4C%se_10(1,ip2):aux4C%se_10(2,ip2),ip2) = mpreal(0.d0)
   do iC=1,nC,2
      ip1 = ((i2+iC)+1)/2
      tmp = sum_1((iC+1)/2)
      aux4C%val_10(ip1,ip2) = aux4C%val_10(ip1,ip2) + tmp
      do j=0,min(i2,iC)-1
         tmp = (tmp/(j+1))*(2*(i2-j)*(iC-j))
         aux4C%val_10(ip1-j-1,ip2) = aux4C%val_10(ip1-j-1,ip2) + tmp
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=1,n2,2
   ip2 = (i2+1)/2
   aux4C%se_01(1,ip2) = max(i2-nC_1,0)/2+1
   aux4C%se_01(2,ip2) = (i2+nC_1)/2+1
   aux4C%val_01(aux4C%se_01(1,ip2):aux4C%se_01(2,ip2),ip2) = mpreal(0.d0)
   do iC=1,nC,2
      ip1 = (i2+iC)/2+1
      tmp = sum_1((iC+1)/2)
      aux4C%val_01(ip1,ip2) = aux4C%val_01(ip1,ip2) + tmp
      do j=0,min(i2,iC)-1
         tmp = (tmp/(j+1))*(2*(i2-j)*(iC-j))
         aux4C%val_01(ip1-j-1,ip2) = aux4C%val_01(ip1-j-1,ip2) + tmp
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=1,n2,2
   ip2 = (i2+1)/2
   aux4C%se_11(1,ip2) = (max(i2-nC_0,1)+1)/2
   aux4C%se_11(2,ip2) = ((i2+nC_0)+1)/2
   aux4C%val_11(aux4C%se_11(1,ip2):aux4C%se_11(2,ip2),ip2) = mpreal(0.d0)
   do iC=0,nC,2
      ip1 = ((i2+iC)+1)/2
      tmp = sum_0(iC/2+1)
      aux4C%val_11(ip1,ip2) = aux4C%val_11(ip1,ip2) + tmp
      do j=0,min(i2,iC)-1
         tmp = (tmp/(j+1))*(2*(i2-j)*(iC-j))
         aux4C%val_11(ip1-j-1,ip2) = aux4C%val_11(ip1-j-1,ip2) + tmp
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

end subroutine aux4C_create

subroutine table4_free(table4)
implicit none
type(table4Data) :: table4

deallocate(table4%val_00,table4%val_10,table4%val_01,table4%val_11)
deallocate(table4%se_00,table4%se_10,table4%se_01,table4%se_11)

end subroutine table4_free

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine prod2H_create(prod2H,n1,n2)
implicit none
type(prod2HData) :: prod2H
integer,intent(in) :: n1,n2
integer :: i1,i2,j

prod2H%n1 = n1
prod2H%n2 = n2

allocate(prod2H%unit(0:prod2H%n1,0:prod2H%n2))

do i2=0,prod2H%n2
   do i1=0,prod2H%n1
      associate(unit => prod2H%unit(i1,i2))
        unit%nu = min(i1,i2)
        allocate(unit%val(unit%nu+1))
      end associate
   enddo
enddo

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC)
do i2=0,prod2H%n2
   do i1=0,prod2H%n1
      associate(unit => prod2H%unit(i1,i2))
        unit%val(unit%nu+1) = 1
        do j=0,unit%nu-1
           unit%val(unit%nu-j) = (unit%val(unit%nu-j+1)/(j+1))*(2*(i1-j)*(i2-j))
        enddo
      end associate
   enddo
enddo
!$OMP END PARALLEL DO

end subroutine prod2H_create

subroutine prod2H_norm(prod2H)
implicit none
type(prod2HData) :: prod2H
integer :: i1,i2,j
type(mp_real),allocatable :: norm(:)
type(mp_real) :: val

allocate(norm(0:max(prod2H%n1,prod2H%n2)))
call G_norm_definition(max(prod2H%n1,prod2H%n2),norm)

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(val)
do i2=0,prod2H%n2
   do i1=0,prod2H%n1
      associate(unit => prod2H%unit(i1,i2))
        val = norm(i1)*norm(i2)
        do j=1,unit%nu+1
           unit%val(j) = unit%val(j)*val
        enddo
      end associate
   enddo
enddo
!$OMP END PARALLEL DO

deallocate(norm)

end subroutine prod2H_norm

subroutine prod2H_free(prod2H)
implicit none
type(prod2HData) :: prod2H
integer :: i1,i2

do i2=0,prod2H%n2
   do i1=0,prod2H%n1
      associate(unit => prod2H%unit(i1,i2))
        deallocate(unit%val)
      end associate
   enddo
enddo

deallocate(prod2H%unit)

end subroutine prod2H_free

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine outer_vecproduct(n,orb1,orb2,mat,set_clean)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: orb1(:),orb2(:)
type(mp_real),intent(inout) :: mat(:,:)
logical,intent(in),optional :: set_clean
logical :: clean
integer :: i,j

if(present(set_clean)) then
   clean = set_clean
else
   clean = .true.
endif

if(clean) then
   do j=1,n
      do i=1,n
         mat(i,j) = orb1(i)*orb2(j)
      enddo
   enddo
else
   do j=1,n
      do i=1,n
         mat(i,j) = mat(i,j) + orb1(i)*orb2(j)
      enddo
   enddo
endif

end subroutine outer_vecproduct

subroutine fullMO_transform(n,C,A,work)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: C(:,:)
type(mp_real),intent(inout) :: a(:,:)
type(mp_real) :: work(:,:)
integer :: i,j,k
type(mp_real) :: val

!$OMP PARALLEL PRIVATE(val)

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      val = 0
      do k=1,n
         val = val + C(k,i)*A(k,j)
      enddo
      work(i,j) = val
   enddo
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      val = 0
      do k=1,n
         val = val + work(i,k)*C(k,j)
      enddo
      A(i,j) = val
   enddo
enddo
!$OMP END DO

!$OMP END PARALLEL

end subroutine fullMO_transform

subroutine scf_FDS_SDF(n,error,P,F,D,S,work)
implicit none
integer,intent(in) :: n
type(mp_real),intent(out) :: error
type(mp_real),intent(out) :: P(:,:)
type(mp_real),intent(in) :: F(:,:),D(:,:),S(:,:)
type(mp_real) :: work(:,:)
integer :: i,j,k
type(mp_real) :: val

! FDS = ( D^ F )^ S  for F^ = F

error = 0

!$OMP PARALLEL PRIVATE(val)

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      val = 0
      do k=1,n
         val = val + D(k,i)*F(k,j)
      enddo
      work(i,j) = val
   enddo
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      val = 0
      do k=1,n
         val = val + work(k,i)*S(k,j)
      enddo
      P(i,j) = val
   enddo
enddo
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC)
do j=1,n
   do i=1,j-1
      val = P(i,j)
      P(i,j) = P(i,j) - P(j,i)
      P(j,i) = P(j,i) - val
   enddo
   P(j,j) = 0
enddo
!$OMP END DO

val = 0

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      val = val + P(i,j)**2
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP CRITICAL(FDS_SDF_error)
error = error + val
!$OMP END CRITICAL(FDS_SDF_error)

!$OMP END PARALLEL

error = sqrt(error)

end subroutine scf_FDS_SDF

end module integrals
