module integrals
use mpmodule
use general
implicit none

private
public integrals_scfSH
public scfJauxData,prepare_scfJaux,free_scfJaux,integrals_scfJ
public integrals_SH,integrals_J
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

type cont4Data
integer :: n_0,n_1
type(mp_real),allocatable :: val_0(:),val_1(:)
end type cont4Data

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
type(mp_real),allocatable :: sum_0(:),sum_1(:)
integer :: n_0,n_1
end type scfJauxData

contains

subroutine integrals_scfSH(scfS,scfH)
implicit none
type(mp_real),intent(out) :: scfS(0:,0:)
type(mp_real),intent(out) :: scfH(0:,0:)
type(mp_real),allocatable :: norm(:)
type(mp_real) :: val
integer :: i

allocate(norm(0:G_nbas-1))
call G_norm_definition(G_nbas-1,norm)

scfS(0:G_nbas-1,0:G_nbas-1) = mpreal(0.d0)

val = sqrt(mppi())
scfS(0,0) = val*norm(0)**2
do i=1,G_nbas-1
   val = val*(2*i)
   scfS(i,i) = val*norm(i)**2
enddo

scfH(0:G_nbas-1,0:G_nbas-1) = mpreal(0.d0)

do i=0,G_nbas-1
   scfH(i,i) = (i+0.5d0)*scfS(i,i)
enddo

deallocate(norm)

end subroutine integrals_scfSH

subroutine prepare_scfJaux(scfJaux)
implicit none
type(scfJauxData) :: scfJaux
integer :: n1,n2,i,j,ij,v2
real(8) :: safe_half
type(table2Data) :: aux22

safe_half = 0.5d0 - 0.5d0/G_nbas

n1 = G_nbas-1
n2 = G_nbas-1
call prod2H_create(scfJaux%prod2H,n1,n2)
call prod2H_norm(scfJaux%prod2H)

n1 = G_nbas
n2 = (G_nbas*(G_nbas+1))/2
allocate(scfJaux%ints(n1,n2))
allocate(scfJaux%sum_0(n1),scfJaux%sum_1(n1))

n1 = 2*(G_nbas-1)
n2 = 2*(G_nbas-1)
call aux22_create(aux22,n1,n2)

scfJaux%n_0 = aux22%n1_0
scfJaux%n_1 = aux22%n1_1

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v2)
do ij=1,(G_nbas*(G_nbas+1))/2
   j = int(sqrt(0.25d0 + 2.d0*(ij-1)) - safe_half)
   i = (ij-1) - (j*(j+1))/2
   v2 = i + j
   if(mod(v2,2)==0) then
      v2 = v2/2+1
      call ints_product2(aux22%n1_0,v2,aux22%val_0,scfJaux%prod2H%unit(i,j),&
           scfJaux%ints(:,ij))
   else
      v2 = (v2+1)/2
      call ints_product2(aux22%n1_1,v2,aux22%val_1,scfJaux%prod2H%unit(i,j),&
           scfJaux%ints(:,ij))
   endif
enddo
!$OMP END PARALLEL DO

call table2_free(aux22)

end subroutine prepare_scfJaux

subroutine integrals_scfJ(scfJ,scfD,scfJaux)
implicit none
type(mp_real),intent(out) :: scfJ(0:,0:)
type(mp_real),intent(in) :: scfD(0:,0:)
type(scfJauxData) :: scfJaux
real(8) :: safe_half
integer :: i,j,ij,v1
type(mp_real) :: val

safe_half = 0.5d0 - 0.5d0/G_nbas

!$OMP PARALLEL

!$OMP SECTIONS PRIVATE(ij,val)
!$OMP SECTION
scfJaux%sum_0(1:scfJaux%n_0) = mpreal(0.d0)
do j=0,G_nbas-1
   do i=mod(j,2),j,2
      ij = 1 + i + (j*(j+1))/2
      val = scfD(i,j)
      if(i/=j) val = val + scfD(j,i)
      call ints_accumulate(scfJaux%n_0,val,scfJaux%ints(:,ij),scfJaux%sum_0)
   enddo
enddo
!$OMP SECTION
scfJaux%sum_1(1:scfJaux%n_1) = mpreal(0.d0)
do j=0,G_nbas-1
   do i=mod(j+1,2),j,2
      ij = 1 + i + (j*(j+1))/2
      val = scfD(i,j) + scfD(j,i)
      call ints_accumulate(scfJaux%n_1,val,scfJaux%ints(:,ij),scfJaux%sum_1)
   enddo
enddo
!$OMP END SECTIONS

!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v1,val)
do ij=1,(G_nbas*(G_nbas+1))/2
   j = int(sqrt(0.25d0 + 2.d0*(ij-1)) - safe_half)
   i = (ij-1) - (j*(j+1))/2
   v1 = i + j
   if(mod(v1,2)==0) then
      v1 = v1/2+1
      val = ints_product1(v1,scfJaux%prod2H%unit(i,j),scfJaux%sum_0)
   else
      v1 = (v1+1)/2
      val = ints_product1(v1,scfJaux%prod2H%unit(i,j),scfJaux%sum_1)
   endif
   val = G_gfac*val
   scfJ(i,j) = val
   if(i/=j) scfJ(j,i) = val
enddo
!$OMP END DO

!$OMP END PARALLEL

end subroutine integrals_scfJ

subroutine free_scfJaux(scfJaux)
implicit none
type(scfJauxData) :: scfJaux

deallocate(scfJaux%sum_0,scfJaux%sum_1)
deallocate(scfJaux%ints)
call prod2H_free(scfJaux%prod2H)

end subroutine free_scfJaux

subroutine integrals_SH(matS_S,matS_T,matH_S,matH_T)
implicit none
type(mp_real),intent(out) :: matS_S(:,:),matS_T(:,:)
type(mp_real),intent(out) :: matH_S(:,:),matH_T(:,:)
integer :: n1,n2,iOP,i,j,k,l,ij,kl,ijkl,v1,v2
real(8) :: safe_half
type(mp_real) :: alpha_invsqrt,a1,a2,a12,c0,c2,c4,val_A,val_B
type(prod2HData) :: prod2H
type(table2Data) :: aux,auxOP
type(mp_real),allocatable :: ints(:,:)
integer,allocatable :: pair_list(:,:)

safe_half = 0.5d0 - 0.5d0/(G_nprim+1)
alpha_invsqrt = 1/sqrt(1 + 4*G_alpha)

n1 = G_nprim
n2 = G_nprim
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = 2*G_nprim
n2 = 2*G_nprim
a1  = -2*G_alpha/(1 + 4*G_alpha)
a2  = -2*G_alpha/(1 + 4*G_alpha)
a12 =  4*G_alpha/(1 + 4*G_alpha)
call table2_create(aux,n1,n2,a1,a2,a12)

n1 = max(aux%n1_0,aux%n1_1)
n2 = ((G_nprim+1)*(G_nprim+2))/2
allocate(ints(n1,n2))

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

do iOP=1,2

   select case(iOP)
   case(1)
      c0 = mppi()*alpha_invsqrt**3
      c2 = mppi()*alpha_invsqrt**5
      call aux22F2_create(auxOP,c0,c2,aux)
   case(2)
      c0 = mppi()*alpha_invsqrt**5*(1 + 4*G_alpha*(1 + 3*G_alpha))
      c2 = mppi()*alpha_invsqrt**7*(4*G_alpha*(2*G_alpha - 1))
      c4 = mppi()*alpha_invsqrt**9*(4*G_alpha**2)
      call aux22DF2_create(auxOP,c0,c2,c4,aux)
   end select

   call table2_fact(auxOP)

!$OMP PARALLEL

!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v2)
   do ij=1,((G_nprim+1)*(G_nprim+2))/2
      j = int(sqrt(0.25d0 + 2.d0*(ij-1)) - safe_half)
      i = (ij-1) - (j*(j+1))/2
      v2 = i + j
      if(mod(v2,2)==0) then
         v2 = v2/2+1
         call ints_product2(auxOP%n1_0,v2,auxOP%val_0,prod2H%unit(i,j),&
              ints(:,ij))
      else
         v2 = (v2+1)/2
         call ints_product2(auxOP%n1_1,v2,auxOP%val_1,prod2H%unit(i,j),&
              ints(:,ij))
      endif
   enddo
!$OMP END DO

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(i,j,k,l,ijkl,v1,v2,val_A,val_B)
   do kl=1,G_npair
      do ij=1,G_npair

         i = pair_list(1,ij)
         j = pair_list(2,ij)
         k = pair_list(1,kl)
         l = pair_list(2,kl)
         ijkl = i + j + k + l

         if(mod(ijkl,2)==0) then

            v1 = i + k
            if(mod(v1,2)==0) then
               v1 = v1/2+1
            else
               v1 = (v1+1)/2
            endif
            v2 = max(j,l)
            v2 = 1 + min(j,l) + (v2*(v2+1))/2
            val_A = ints_product1(v1,prod2H%unit(i,k),ints(:,v2))

            v1 = i + l
            if(mod(v1,2)==0) then
               v1 = v1/2 + 1
            else
               v1 = (v1+1)/2
            endif
            v2 = max(j,k)
            v2 = 1 + min(j,k) + (v2*(v2+1))/2
            val_B = ints_product1(v1,prod2H%unit(i,l),ints(:,v2))

            select case(iOP)
            case(1)
               matS_S(ij,kl) = (val_A + val_B)
               matS_T(ij,kl) = (val_A - val_B)
            case(2)
               matH_S(ij,kl) = (val_A + val_B) + matS_S(ij,kl)*((ijkl+2)*0.5d0)
               matH_T(ij,kl) = (val_A - val_B) + matS_T(ij,kl)*((ijkl+2)*0.5d0)
            end select

         else

            select case(iOP)
            case(1)
               matS_S(ij,kl) = 0
               matS_T(ij,kl) = 0
            case(2)
               matH_S(ij,kl) = 0
               matH_T(ij,kl) = 0
            end select

         endif

      enddo
   enddo
!$OMP END DO

!$OMP END PARALLEL

   call table2_free(auxOP)

enddo

deallocate(pair_list)
deallocate(ints)
call table2_free(aux)
call prod2H_free(prod2H)

!!$block
!!$  integer :: idat,ios,i1,i2
!!$  character(210) :: sval
!!$  type(mp_real) :: thr,math,calc,tmp
!!$  thr = '1.e-199'
!!$  open(newunit=idat,file='matHS.dat')
!!$  do
!!$     read(idat,*,iostat=ios) i1,i2,sval
!!$     if(ios/=0) exit
!!$     if(i1<=G_npair.and.i2<=G_npair) then
!!$        math = mpreal(sval)
!!$        calc = matH_S(i1,i2)
!!$        if(math/=0) then
!!$           tmp = (calc - math)/math
!!$        else
!!$           tmp = calc - math
!!$        endif
!!$        if(abs(tmp)>thr) then
!!$           write(*,'(2i5,a)',advance='no') i1,i2,'   '
!!$           call mpwrite(6,60,40,tmp)
!!$        endif
!!$     endif
!!$  enddo
!!$  close(idat)
!!$end block

end subroutine integrals_SH

subroutine integrals_J(matJ_S,matJ_T,matC)
implicit none
type(mp_real),intent(out) :: matJ_S(:,:,:,:),matJ_T(:,:,:,:)
type(mp_real),intent(in) :: matC(:,:)
integer :: n1,n2,nC,iocc,jocc,i,j,k,l,ij,kl,v2
real(8) :: safe_half
type(mp_real) :: alpha_invsqrt,a1,a2,a12,c0,c2
type(mp_real) :: val
type(prod2HData) :: prod2H
type(table2Data) :: aux,auxOP
type(table4Data) :: aux2C
type(mp_real),allocatable :: ints(:,:)
type(mp_real),allocatable :: matD(:,:),sum_0(:),sum_1(:)
type(cont4Data),allocatable :: cont4(:)

safe_half = 0.5d0 - 0.5d0/(G_nprim+1)
alpha_invsqrt = 1/sqrt(1 + 3*G_alpha)

n1 = max(G_nprim,G_nbas-1)
n2 = max(G_nprim,G_nbas-1)
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = G_nprim + G_nbas
n2 = ((G_nprim+1)*(G_nprim+2))/2
allocate(ints(n1,n2))

n1 = 2*G_nprim + 2*(G_nbas-1)
n2 = 2*G_nprim
a1  = -(1 + 4*G_alpha)/(2*(1 + 3*G_alpha))
a2  = -2*G_alpha/(1 + 3*G_alpha)
a12 =  2*G_alpha/(1 + 3*G_alpha)
call table2_create(aux,n1,n2,a1,a2,a12)

c0 = mppi()*(3/(4*sqrt(mpreal(2.d0))))*alpha_invsqrt**3
c2 = mppi()*(1/(4*sqrt(mpreal(2.d0))))*alpha_invsqrt**5
call aux42F2_create(auxOP,c0,c2,aux)
call table2_fact(auxOP)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v2)
do ij=1,((G_nprim+1)*(G_nprim+2))/2
   j = int(sqrt(0.25d0 + 2.d0*(ij-1)) - safe_half)
   i = (ij-1) - (j*(j+1))/2
   v2 = i + j
   if(mod(v2,2)==0) then
      v2 = v2/2+1
      call ints_product2(auxOP%n1_0,v2,auxOP%val_0,prod2H%unit(i,j),&
           ints(:,ij))
   else
      v2 = (v2+1)/2
      call ints_product2(auxOP%n1_1,v2,auxOP%val_1,prod2H%unit(i,j),&
           ints(:,ij))
   endif
enddo
!$OMP END PARALLEL DO

call table2_free(auxOP)
call table2_free(aux)

allocate(matD(0:G_nbas-1,0:G_nbas-1),sum_0(G_nbas),sum_1(G_nbas))

nC = ((G_nprim+1)*(G_nprim+2))/2
allocate(cont4(nC))
ij = 0
do j=0,G_nprim
   do i=0,j
      ij = ij + 1
      associate(this => cont4(ij))
        nC = i + j + 2*(G_nbas-1)
        this%n_0 = nC/2+1
        this%n_1 = (nC+1)/2
        allocate(this%val_0(this%n_0),this%val_1(this%n_1))
      end associate
   enddo
enddo

do jocc=1,G_nocc
   do iocc=1,jocc

      call outer_vecproduct(G_nbas,matC(:,iocc),matC(:,jocc),matD)

!$OMP PARALLEL SECTIONS PRIVATE(v2,val)
!$OMP SECTION
      sum_0(1:G_nbas) = mpreal(0.d0)
      do j=0,G_nbas-1
         do i=mod(j,2),j,2
            v2 = (i+j)/2+1
            val = matD(i,j)
            if(i/=j) val = val + matD(j,i)
            call prod2H_accumulate(v2,val,prod2H%unit(i,j),sum_0)
         enddo
      enddo
!$OMP SECTION
      sum_1(1:G_nbas) = mpreal(0.d0)
      do j=0,G_nbas-1
         do i=mod(j+1,2),j,2
            v2 = (i+j+1)/2
            val = matD(i,j) + matD(j,i)
            call prod2H_accumulate(v2,val,prod2H%unit(i,j),sum_1)
         enddo
      enddo
!$OMP END PARALLEL SECTIONS

      n1 = 2*G_nprim + 2*(G_nbas-1)
      n2 = 2*G_nprim
      nC = 2*(G_nbas-1)
      call aux2C_create(aux2C,n1,n2,nC,sum_0,sum_1)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v2)
      do ij=1,((G_nprim+1)*(G_nprim+2))/2
         j = int(sqrt(0.25d0 + 2.d0*(ij-1)) - safe_half)
         i = (ij-1) - (j*(j+1))/2
         v2 = i + j
         if(mod(v2,2)==0) then
            v2 = v2/2+1
            call cont_product2(&
                 aux2C%se_00,aux2C%val_00,v2,prod2H%unit(i,j),&
                 cont4(ij)%n_0,cont4(ij)%val_0)
            call cont_product2(&
                 aux2C%se_10,aux2C%val_10,v2,prod2H%unit(i,j),&
                 cont4(ij)%n_1,cont4(ij)%val_1)
         else
            v2 = (v2+1)/2
            call cont_product2(&
                 aux2C%se_01,aux2C%val_01,v2,prod2H%unit(i,j),&
                 cont4(ij)%n_0,cont4(ij)%val_0)
            call cont_product2(&
                 aux2C%se_11,aux2C%val_11,v2,prod2H%unit(i,j),&
                 cont4(ij)%n_1,cont4(ij)%val_1)
         endif
      enddo
!$OMP END PARALLEL DO

      call table4_free(aux2C)

      kl = 0
      do l=0,G_nprim
         do k=0,l
            kl = kl + 1

            v2 = mod(k+l,2)

            ij = 0
            do j=0,G_nprim
               do i=0,j
                  ij = ij + 1

!!$                  select case(v2)
!!$                  case(0)
!!$                     val = cont_product1(cont4(ij)%n_0,cont4(ij)%val_0,&
!!$                          ints(:,kl))
!!$                  case(1)
!!$                     val = cont_product1(cont4(ij)%n_1,cont4(ij)%val_1,&
!!$                          ints(:,kl))
!!$                  end select
!!$
!!$                  write(*,*) i,j,k,l,dble(val)

!!$                  associate(this => cont4(ij))
!!$              if(mod(i+j,2)==1) then
!!$                 write(*,*) i,j,this%n_1
!!$                 do v2=1,this%n_1
!!$                    write(*,*)  v2,dble(this%val_1(v2))
!!$                 enddo
!!$              endif
!!$            end associate

               enddo
            enddo

         enddo
      enddo

   enddo
enddo

ij = 0
do j=0,G_nprim
   do i=0,j
      ij = ij + 1
      associate(this => cont4(ij))
        deallocate(this%val_0,this%val_1)
      end associate
   enddo
enddo
deallocate(cont4)

deallocate(matD,sum_0,sum_1)

deallocate(ints)
call prod2H_free(prod2H)

!!$block
!!$  integer :: i,j,k,l,p,q,v1,v2,v
!!$  type(mp_real),allocatable :: tabs(:)
!!$  type(mp_real) :: val
!!$  do l=0,G_nprim
!!$     do k=0,G_nprim
!!$        do q=0,G_nbas-1
!!$           do p=0,G_nbas-1
!!$              do j=0,G_nprim
!!$                 do i=0,G_nprim
!!$  v1 = i + j + p + q
!!$  v2 = k + l
!!$  if(mod(v1+v2,2)==0) then
!!$     if(mod(v1,2)==0) then
!!$        v1 = v1/2+1
!!$     else
!!$        v1 = (v1+1)/2
!!$     endif
!!$     v2 = max(k,l)
!!$     v2 = 1 + min(k,l) + (v2*(v2+1))/2
!!$     allocate(tabs(v1))
!!$     call combine_units(v1,tabs,prod2H%unit(i,j),prod2H%unit(p,q),&
!!$          i+j,p+q,prod2)
!!$     val = 0
!!$     do v=1,v1
!!$        val = val + tabs(v)*ints(v,v2)
!!$     enddo
!!$     deallocate(tabs)
!!$  else
!!$     val = 0
!!$  endif
!!$  write(*,'(6i4,3x)',advance='no') i,j,p,q,k,l
!!$  call mpwrite(6,50,40,val)
!!$enddo
!!$enddo
!!$enddo
!!$enddo
!!$enddo
!!$enddo
!!$end block

!!$block
!!$  integer :: idat,ios,i1,i2
!!$  character(210) :: sval
!!$  type(mp_real) :: thr,math,calc,tmp
!!$  thr = '1.e-199'
!!$  open(newunit=idat,file='e.dat')
!!$  do
!!$     read(idat,*,iostat=ios) i1,i2,sval
!!$     if(ios/=0) exit
!!$     if(i1<=n1.and.i2<=n2) then
!!$        math = mpreal(sval)
!!$        calc = table2_take(auxOP,i1,i2)
!!$        if(math/=0) then
!!$           tmp = (calc - math)/math
!!$        else
!!$           tmp = calc - math
!!$        endif
!!$        if(abs(tmp)>thr) then
!!$           write(*,'(2i5,a)',advance='no') i1,i2,'   '
!!$           call mpwrite(6,60,40,tmp)
!!$        endif
!!$     endif
!!$  enddo
!!$  close(idat)
!!$end block

end subroutine integrals_J

function ints_product1(v1,unit1,ints) result(val)
implicit none
type(mp_real) :: val
integer,intent(in) :: v1
type(prodUnitData),intent(in) :: unit1
type(mp_real),intent(in) :: ints(:)
integer :: off1,i1

off1 = v1 - unit1%nu

val = unit1%val(1)*ints(off1)
do i1=1,unit1%nu
   val = val + unit1%val(1+i1)*ints(off1+i1)
enddo

end function ints_product1

subroutine ints_product2(n1,v2,aux,unit2,ints)
implicit none
integer,intent(in) :: n1,v2
type(mp_real),intent(in) :: aux(:,:)
type(prodUnitData),intent(in) :: unit2
type(mp_real),intent(out) :: ints(:)
integer :: off2,i1,i2

off2 = v2 - unit2%nu

do i1=1,n1
   ints(i1) = aux(i1,off2)*unit2%val(1)
enddo
do i2=1,unit2%nu
   do i1=1,n1
      ints(i1) = ints(i1) + aux(i1,off2+i2)*unit2%val(1+i2)
   enddo
enddo

end subroutine ints_product2

subroutine ints_accumulate(n,alpha,ints,acc)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: alpha,ints(:)
type(mp_real),intent(inout) :: acc(:)
integer :: i

do i=1,n
   acc(i) = acc(i) + alpha*ints(i)
enddo

end subroutine ints_accumulate

subroutine prod2H_accumulate(v,alpha,unit,acc)
implicit none
integer,intent(in) :: v
type(mp_real),intent(in) :: alpha
type(prodUnitData),intent(in) :: unit
type(mp_real),intent(inout) :: acc(:)
integer :: off,i

off = v - unit%nu

do i=0,unit%nu
   acc(off+i) = acc(off+i) + alpha*unit%val(1+i)
enddo

end subroutine prod2H_accumulate

function cont_product1(n,vec1,vec2) result(val)
implicit none
type(mp_real) :: val
integer,intent(in) :: n
type(mp_real),intent(in) :: vec1(:),vec2(:)
integer :: i

val = 0
do i=1,n
   val = val + vec1(i)*vec2(i)
enddo

end function cont_product1

subroutine cont_product2(se,mat,v,unit,n,cont)
implicit none
integer,intent(in) :: se(:,:)
type(mp_real),intent(in) :: mat(:,:)
integer,intent(in) :: v
type(prodUnitData),intent(in) :: unit
integer,intent(in) :: n
type(mp_real),intent(out) :: cont(:)
integer :: off,i1,i2

cont(1:n) = mpreal(0.d0)

off = v - unit%nu

do i2=0,unit%nu
   do i1=se(1,off+i2),se(2,off+i2)
      cont(i1) = cont(i1) + mat(i1,off+i2)*unit%val(1+i2)
   enddo
enddo

end subroutine cont_product2

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

!$OMP PARALLEL

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(i12,tmp,val)
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
!$OMP END DO

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(i12,tmp,val)
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
!$OMP END DO

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
!$OMP END SINGLE

!$OMP DO
do i1=2,aux22F2%n1_0
   aux22F2%val_0(i1,1) = c0*aux%val_0(i1,1) &
        + c2*aux%val_0(i1-1,1)
enddo
!$OMP END DO

!$OMP DO
do i2=2,aux22F2%n2_0
   aux22F2%val_0(1,i2) = c0*aux%val_0(1,i2) &
        + c2*aux%val_0(1,i2-1)
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
do i2=2,aux22F2%n2_0
   do i1=2,aux22F2%n1_0
      aux22F2%val_0(i1,i2) = c0*aux%val_0(i1,i2) &
           + c2*(aux%val_0(i1-1,i2) + aux%val_0(i1,i2-1) &
           - 2*aux%val_1(i1-1,i2-1))
   enddo
enddo
!$OMP END DO

!$OMP SINGLE
aux22F2%val_1(1,1) = c0*aux%val_1(1,1) &
     - 2*c2*aux%val_0(1,1)
!$OMP END SINGLE

!$OMP DO
do i1=2,aux22F2%n1_1
   aux22F2%val_1(i1,1) = c0*aux%val_1(i1,1) &
        + c2*(aux%val_1(i1-1,1) - 2*aux%val_0(i1,1))
enddo
!$OMP END DO

!$OMP DO
do i2=2,aux22F2%n2_1
   aux22F2%val_1(1,i2) = c0*aux%val_1(1,i2) &
        + c2*(aux%val_1(1,i2-1) - 2*aux%val_0(1,i2))
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
do i2=2,aux22F2%n2_1
   do i1=2,aux22F2%n1_1
      aux22F2%val_1(i1,i2) = c0*aux%val_1(i1,i2) &
           + c2*(aux%val_1(i1-1,i2) + aux%val_1(i1,i2-1) &
           - 2*aux%val_0(i1,i2))
   enddo
enddo
!$OMP END DO

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
!$OMP END SECTIONS

!$OMP DO
do i1=3,aux22DF2%n1_0
   aux22DF2%val_0(i1,1) = c0*aux%val_0(i1,1) &
        + c2*aux%val_0(i1-1,1) &
        + c4*aux%val_0(i1-2,1)
enddo
!$OMP END DO

!$OMP DO
do i1=3,aux22DF2%n1_0
   aux22DF2%val_0(i1,2) = c0*aux%val_0(i1,2) &
        + c2*(aux%val_0(i1-1,2) + aux%val_0(i1,1) - 2*aux%val_1(i1-1,1))&
        + c4*(aux%val_0(i1-2,2) + 6*aux%val_0(i1-1,1) - 4*aux%val_1(i1-2,1))
enddo
!$OMP END DO

!$OMP DO
do i2=3,aux22DF2%n2_0
   aux22DF2%val_0(1,i2) = c0*aux%val_0(1,i2) &
        + c2*aux%val_0(1,i2-1) &
        + c4*aux%val_0(1,i2-2)
enddo
!$OMP END DO

!$OMP DO
do i2=3,aux22DF2%n2_0
   aux22DF2%val_0(2,i2) = c0*aux%val_0(2,i2) &
        + c2*(aux%val_0(2,i2-1) + aux%val_0(1,i2) - 2*aux%val_1(1,i2-1))&
        + c4*(aux%val_0(2,i2-2) + 6*aux%val_0(1,i2-1) - 4*aux%val_1(1,i2-2))
enddo
!$OMP END DO

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
!$OMP END DO

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
!$OMP END SECTIONS

!$OMP DO
do i1=3,aux22DF2%n1_1
   aux22DF2%val_1(i1,1) = c0*aux%val_1(i1,1) &
        + c2*(aux%val_1(i1-1,1) - 2*aux%val_0(i1,1)) &
        + c4*(aux%val_1(i1-2,1) - 4*aux%val_0(i1-1,1))
enddo
!$OMP END DO

!$OMP DO
do i1=3,aux22DF2%n1_1
   aux22DF2%val_1(i1,2) = c0*aux%val_1(i1,2) &
        + c2*(aux%val_1(i1-1,2) + aux%val_1(i1,1) - 2*aux%val_0(i1,2)) &
        + c4*(aux%val_1(i1-2,2) + 6*aux%val_1(i1-1,1) &
        - 4*aux%val_0(i1-1,2) - 4*aux%val_0(i1,1))
enddo
!$OMP END DO

!$OMP DO
do i2=3,aux22DF2%n2_1
   aux22DF2%val_1(1,i2) = c0*aux%val_1(1,i2) &
        + c2*(aux%val_1(1,i2-1) - 2*aux%val_0(1,i2)) &
        + c4*(aux%val_1(1,i2-2) - 4*aux%val_0(1,i2-1))
enddo
!$OMP END DO

!$OMP DO
do i2=3,aux22DF2%n2_1
   aux22DF2%val_1(2,i2) = c0*aux%val_1(2,i2) &
        + c2*(aux%val_1(2,i2-1) + aux%val_1(1,i2) - 2*aux%val_0(2,i2)) &
        + c4*(aux%val_1(2,i2-2) + 6*aux%val_1(1,i2-1) &
        - 4*aux%val_0(2,i2-1) - 4*aux%val_0(1,i2))
enddo
!$OMP END DO

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
!$OMP END DO

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
!$OMP END SINGLE

!$OMP DO
do i1=2,aux42F2%n1_0
   aux42F2%val_0(i1,1) = c0*aux%val_0(i1,1) &
        + c2*aux%val_0(i1-1,1)
enddo
!$OMP END DO

!$OMP DO
do i2=2,aux42F2%n2_0
   aux42F2%val_0(1,i2) = c0*aux%val_0(1,i2) &
        + 4*c2*aux%val_0(1,i2-1)
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
do i2=2,aux42F2%n2_0
   do i1=2,aux42F2%n1_0
      aux42F2%val_0(i1,i2) = c0*aux%val_0(i1,i2) &
           + c2*(aux%val_0(i1-1,i2) + 4*aux%val_0(i1,i2-1) &
           - 4*aux%val_1(i1-1,i2-1))
   enddo
enddo
!$OMP END DO

!$OMP SINGLE
aux42F2%val_1(1,1) = c0*aux%val_1(1,1) &
     - 4*c2*aux%val_0(1,1)
!$OMP END SINGLE

!$OMP DO
do i1=2,aux42F2%n1_1
   aux42F2%val_1(i1,1) = c0*aux%val_1(i1,1) &
        + c2*(aux%val_1(i1-1,1) - 4*aux%val_0(i1,1))
enddo
!$OMP END DO

!$OMP DO
do i2=2,aux42F2%n2_1
   aux42F2%val_1(1,i2) = c0*aux%val_1(1,i2) &
        + c2*(4*aux%val_1(1,i2-1) - 4*aux%val_0(1,i2))
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
do i2=2,aux42F2%n2_1
   do i1=2,aux42F2%n1_1
      aux42F2%val_1(i1,i2) = c0*aux%val_1(i1,i2) &
           + c2*(aux%val_1(i1-1,i2) + 4*aux%val_1(i1,i2-1) &
           - 4*aux%val_0(i1,i2))
   enddo
enddo
!$OMP END DO

!$OMP END PARALLEL

end subroutine aux42F2_create

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

!$OMP PARALLEL

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

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(val)
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

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(val)
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

subroutine aux2C_create(aux2C,n1,n2,nC,sum_0,sum_1)
implicit none
type(table4Data) :: aux2C
integer,intent(in) :: n1,n2,nC
type(mp_real),intent(in) :: sum_0(:),sum_1(:)
integer :: nC_0,nC_1
integer :: i2,iC,j,ip1,ip2
type(mp_real) :: tmp

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

!$OMP PARALLEL PRIVATE(ip1,ip2,tmp)

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   aux2C%se_00(1,ip2) = max(i2-nC_0,0)/2+1
   aux2C%se_00(2,ip2) = (i2+nC_0)/2+1
   aux2C%val_00(aux2C%se_00(1,ip2):aux2C%se_00(2,ip2),ip2) = mpreal(0.d0)
   do iC=0,nC,2
      ip1 = (i2+iC)/2+1
      tmp = sum_0(iC/2+1)
      aux2C%val_00(ip1,ip2) = aux2C%val_00(ip1,ip2) + tmp
      do j=0,min(i2,iC)-1
         tmp = (tmp/(j+1))*(2*(i2-j)*(iC-j))
         aux2C%val_00(ip1-j-1,ip2) = aux2C%val_00(ip1-j-1,ip2) + tmp
      enddo
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
      tmp = sum_1((iC+1)/2)
      aux2C%val_10(ip1,ip2) = aux2C%val_10(ip1,ip2) + tmp
      do j=0,min(i2,iC)-1
         tmp = (tmp/(j+1))*(2*(i2-j)*(iC-j))
         aux2C%val_10(ip1-j-1,ip2) = aux2C%val_10(ip1-j-1,ip2) + tmp
      enddo
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
      tmp = sum_1((iC+1)/2)
      aux2C%val_01(ip1,ip2) = aux2C%val_01(ip1,ip2) + tmp
      do j=0,min(i2,iC)-1
         tmp = (tmp/(j+1))*(2*(i2-j)*(iC-j))
         aux2C%val_01(ip1-j-1,ip2) = aux2C%val_01(ip1-j-1,ip2) + tmp
      enddo
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
      tmp = sum_0(iC/2+1)
      aux2C%val_11(ip1,ip2) = aux2C%val_11(ip1,ip2) + tmp
      do j=0,min(i2,iC)-1
         tmp = (tmp/(j+1))*(2*(i2-j)*(iC-j))
         aux2C%val_11(ip1-j-1,ip2) = aux2C%val_11(ip1-j-1,ip2) + tmp
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

end subroutine aux2C_create

subroutine table4_free(table4)
implicit none
type(table4Data) :: table4

deallocate(table4%val_00,table4%val_10,table4%val_01,table4%val_11)
deallocate(table4%se_00,table4%se_10,table4%se_01,table4%se_11)

end subroutine table4_free

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

!$OMP PARALLEL

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      work(i,j) = 0
      do k=1,n
         work(i,j) = work(i,j) + C(k,i)*A(k,j)
      enddo
   enddo
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      A(i,j) = 0
      do k=1,n
         A(i,j) = A(i,j) + work(i,k)*C(k,j)
      enddo
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

!$OMP PARALLEL

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      work(i,j) = 0
      do k=1,n
         work(i,j) = work(i,j) + D(k,i)*F(k,j)
      enddo
   enddo
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      P(i,j) = 0
      do k=1,n
         P(i,j) = P(i,j) + work(k,i)*S(k,j)
      enddo
   enddo
enddo
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(val)
do j=1,n
   do i=1,j-1
      val = P(i,j)
      P(i,j) = P(i,j) - P(j,i)
      P(j,i) = P(j,i) - val
   enddo
   P(j,j) = 0
enddo
!$OMP END DO

!$OMP END PARALLEL

error = 0
do j=1,n
   do i=1,n
      error = error + P(i,j)**2
   enddo
enddo
error = sqrt(error)

end subroutine scf_FDS_SDF

!-------------------------------------------------------------------------------

!!$block
!!$  integer :: i,j,k,l,ijkl,v1,v2
!!$  type(mp_real) :: val
!!$  i = 4
!!$  j = 7
!!$  k = 6
!!$  l = 9
!!$  ijkl = i + j + k + l
!!$  if(mod(ijkl,2)==0) then
!!$     v1 = i + j
!!$     v2 = k + l
!!$     if(mod(v2,2)==0) then
!!$        v1 = v1/2 + 1
!!$        v2 = v2/2 + 1
!!$        val = ints_product(v1,v2,aux22%val_0,&
!!$             scfJaux%prod2H%unit(i,j),&
!!$             scfJaux%prod2H%unit(k,l))
!!$     else
!!$        v1 = (v1+1)/2
!!$        v2 = (v2+1)/2
!!$        val = ints_product(v1,v2,aux22%val_1,&
!!$             scfJaux%prod2H%unit(i,j),&
!!$             scfJaux%prod2H%unit(k,l))
!!$     endif
!!$  else
!!$     val = 0
!!$  endif
!!$  call mpwrite(6,60,40,val)
!!$end block

!!$block
!!$  integer :: i,j,k,l,ijkl,v1,v2
!!$  type(mp_real) :: val
!!$  i = 4
!!$  j = 5
!!$  k = 2
!!$  l = 3
!!$  ijkl = i + j + k + l
!!$  if(mod(ijkl,2)==0) then
!!$     v1 = i + j
!!$     if(mod(v1,2)==0) then
!!$        v1 = v1/2+1
!!$     else
!!$        v1 = (v1+1)/2
!!$     endif
!!$     v2 = max(k,l)
!!$     v2 = 1 + min(k,l) + (v2*(v2+1))/2
!!$     val = ints_product1(v1,scfJaux%prod2H%unit(i,j),scfJaux%ints(:,v2))
!!$  else
!!$     val = 0
!!$  endif
!!$  call mpwrite(6,60,40,val)
!!$end block

!!$function ints_product(v1,v2,aux,unit1,unit2) result(val)
!!$implicit none
!!$type(mp_real) :: val
!!$integer,intent(in) :: v1,v2
!!$type(mp_real),intent(in) :: aux(:,:)
!!$type(prodUnitData),intent(in) :: unit1,unit2
!!$integer :: off1,off2,i1,i2
!!$
!!$off1 = v1 - unit1%nu
!!$off2 = v2 - unit2%nu
!!$
!!$val = 0
!!$do i2=0,unit2%nu
!!$   do i1=0,unit1%nu
!!$      val = val + unit1%val(1+i1)*aux(off1+i1,off2+i2)*unit2%val(1+i2)
!!$   enddo
!!$enddo
!!$
!!$end function ints_product

!!$function table2_take(table2,i1,i2) result(val)
!!$implicit none
!!$type(mp_real) :: val
!!$type(table2Data),intent(in) :: table2
!!$integer,intent(in) :: i1,i2
!!$
!!$if(mod(i2,2)==0) then
!!$   if(mod(i1,2)==0) then
!!$      val = table2%val_0(i1/2+1,i2/2+1)
!!$      return
!!$   endif
!!$else
!!$   if(mod(i1,2)==1) then
!!$      val = table2%val_1((i1+1)/2,(i2+1)/2)
!!$      return
!!$   endif
!!$endif
!!$
!!$val = 0
!!$
!!$end function table2_take

!!$subroutine combine_units(n,tabs,unit1,unit2,v1,v2,prod2)
!!$implicit none
!!$integer,intent(in) :: n
!!$type(mp_real),intent(out) :: tabs(:)
!!$type(prodUnitData),intent(in) :: unit1,unit2
!!$integer,intent(in) :: v1,v2
!!$type(prod2HData),intent(in) :: prod2
!!$integer :: i1,i2,off,j
!!$type(mp_real) :: val
!!$
!!$tabs = mpreal(0.d0)
!!$
!!$do i2=0,unit2%nu
!!$   do i1=0,unit1%nu
!!$
!!$      val = unit1%val(unit1%nu+1-i1)*unit2%val(unit2%nu+1-i2)
!!$
!!$      associate(unit => prod2%unit(v1-2*i1,v2-2*i2))
!!$        off = n - (i1 + i2) - unit%nu
!!$        do j=0,unit%nu
!!$           tabs(off+j) = tabs(off+j) + val*unit%val(1+j)
!!$        enddo
!!$      end associate
!!$
!!$   enddo
!!$enddo
!!$
!!$end subroutine combine_units

end module integrals
