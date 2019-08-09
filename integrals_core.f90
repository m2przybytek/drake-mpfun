module integrals_core
use mpmodule
use global
use misc
implicit none

private
public prod2HData,prod2H_create,prod2H_norm,prod2H_free
public int2Data,int2_create,int2_free
public int3Data,int3_create,int3_from2_create,int3_free
public auxData,aux_init,aux_free,aux2H_create,auxC1_create,auxC2_create

type prod2HData
integer :: n1,n2
type(prodUnitData),allocatable :: unit(:,:)
end type prod2HData

type int2Data
integer :: n1_0,n2_0
integer :: n1_1,n2_1
type(mp_real),allocatable :: val_0(:,:),val_1(:,:)
end type int2Data

type int3Data
integer :: n1_0,n2_0,n3_0
integer :: n1_1,n2_1,n3_1
type(mp_real),allocatable :: val_000(:,:,:)
type(mp_real),allocatable :: val_110(:,:,:)
type(mp_real),allocatable :: val_101(:,:,:)
type(mp_real),allocatable :: val_011(:,:,:)
end type int3Data

type auxData
integer :: n1  ,n2  ,nC
integer :: n1_0,n2_0,nC_0
integer :: n1_1,n2_1,nC_1
integer,allocatable :: se_00(:,:),se_10(:,:),se_01(:,:),se_11(:,:)
type(mp_real),allocatable :: val_00(:,:),val_10(:,:),val_01(:,:),val_11(:,:)
end type auxData

contains

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

!!$block
!!$  character(FLD_DEFAULT) :: s
!!$  write(*,*) '-----------------------------------------------------------------'
!!$  do i2=0,prod2H%n2
!!$     do i1=0,prod2H%n1
!!$        associate(unit => prod2H%unit(i1,i2))
!!$          do j=unit%nu+1,1,-1
!!$             call write_fform(unit%val(j),s)
!!$             write(*,*) i1,i2,j,trim(s)
!!$          enddo
!!$        end associate
!!$        write(*,*) '----------'
!!$     enddo
!!$  enddo
!!$end block
!!$stop

end subroutine prod2H_create

subroutine prod2H_norm(prod2H)
implicit none
type(prod2HData) :: prod2H
integer :: i1,i2,j
type(mp_real) :: val
type(mp_real),allocatable :: norm(:)

allocate(norm(0:max(prod2H%n1,prod2H%n2)))
call norm_definition(max(prod2H%n1,prod2H%n2),norm)

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

subroutine int2_create(f1,f2,int2_type,int2,n1,n2,alpha_f,nC_f,C_f)
implicit none
character(2),intent(in) :: f1,f2,int2_type
type(int2Data) :: int2
integer,intent(in) :: n1,n2
type(mp_real),intent(in) :: alpha_f
integer,intent(in) :: nC_f
type(mp_real),intent(in) :: C_f(0:)
integer :: alpha_mult,nP,nC1,nC2,nC,i1,i2,i,j
type(mp_real) :: val,a0,a1,ai,bi
type(mp_real),allocatable :: P(:),P0(:),PC(:)
type(mp_real),allocatable :: C1(:),C2(:),C(:)
type(mp_real),allocatable :: comb1(:,:),power2(:)

int2%n1_0 = n1/2+1
int2%n2_0 = n2/2+1

int2%n1_1 = (n1+1)/2
int2%n2_1 = (n2+1)/2

allocate(&
     int2%val_0(int2%n1_0,int2%n2_0),&
     int2%val_1(int2%n1_1,int2%n2_1))

nP = (n1+n2)/2
allocate(P(0:nP))

alpha_mult = 0
call int2_choose(f1,nC1,C1,alpha_mult,alpha_f,nC_f,C_f)
call int2_choose(f2,nC2,C2,alpha_mult,alpha_f,nC_f,C_f)

nC = nC1 + nC2
allocate(C(0:nC))

!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
do i=0,nC
   C(i) = 0
   do i1=max(0,i-nC2),min(nC1,i)
      C(i) = C(i) + C1(i1)*C2(i-i1) 
   enddo
enddo
!$OMP END PARALLEL DO

if(alpha_mult==0) then

   P(0) = C(0)*sqrt(mppi()/2)
   do i=1,nP
      P(i) = P(i-1)*(1-2*i)
   enddo

else

   select case(int2_type)
   case('22')

      val = 1/(1 + 2*alpha_mult*alpha_f)
      a0 = mppi()*sqrt(val)
      a1 = sqrt(2*mppi())*val
      ai = val

      bi = -1

      do i=0,nP
         P(i) = 1
      enddo

   case('42')

      val = 1/(2 + 3*alpha_mult*alpha_f)
      a0 = mppi()*sqrt(val)
      a1 = sqrt(3*mppi())*val
      ai = (val*3)/2

      bi = -3; bi = bi/4

      val = 4; val = val/9
      P(0) = 1
      do i=1,nP
         P(i) = P(i-1)*val
      enddo

   case default
      write(*,'(2a)') 'Unknown int2_type in int2_create: ',int2_type
      stop
   end select

   allocate(P0(0:nP),PC(0:(nP+nC/2)))

!$OMP PARALLEL PRIVATE(val)

!$OMP DO
   do i=0,nP
      P0(i) = 0
   enddo
!$OMP END DO

!$OMP SINGLE
   PC(0) = a0
   do i=1,nP+nC/2
      PC(i) = PC(i-1)*ai*(2*i-1)
   enddo
!$OMP END SINGLE

!$OMP DO
   do i=0,nP
      do j=0,nC,2
         P0(i) = P0(i) + PC(i+j/2)*C(j)
      enddo
   enddo
!$OMP END DO
   
!$OMP SINGLE
   PC(0) = a1
   do i=1,nP+(nC-1)/2
      PC(i) = PC(i-1)*ai*(2*i)
   enddo
!$OMP END SINGLE

!$OMP DO
   do i=0,nP
      do j=1,nC,2
         P0(i) = P0(i) + PC(i+(j-1)/2)*C(j)
      enddo
   enddo
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC)   
   do i=0,nP
      val = P(i)
      P(i) = P0(i)*val
      do j=1,i
         val = (val*bi/j)*((i-j+1)*(2*(i-j)+1))
         P(i) = P(i) + P0(i-j)*val
      enddo
   enddo
!$OMP END DO

!$OMP END PARALLEL   

   deallocate(P0,PC)

endif

deallocate(C1,C2,C)

select case(int2_type)
case('22')

!$OMP PARALLEL

!$OMP DO COLLAPSE(2)
   do i2=1,int2%n2_0
      do i1=1,int2%n1_0
         int2%val_0(i1,i2) = P(i1+i2-2)
      enddo
   enddo
!$OMP END DO

!$OMP DO
   do i=0,nP
      P(i) = -P(i)
   enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
   do i2=1,int2%n2_1
      do i1=1,int2%n1_1
         int2%val_1(i1,i2) = P(i1+i2-1)
      enddo
   enddo
!$OMP END DO

!$OMP END PARALLEL

case('42')

   allocate(comb1(int2%n1_0,int2%n1_0),power2(int2%n2_0))

!$OMP PARALLEL PRIVATE(val)

!$OMP DO SCHEDULE(DYNAMIC)
   do i1=1,int2%n1_0
      comb1(1,i1) = 1
      do j=1,i1-1
         comb1(j+1,i1) = (comb1(j,i1)/j)*((j-i1)*(2*(i1-j)-1))
      enddo
   enddo
!$OMP END DO NOWAIT

!$OMP SINGLE
   power2(1) = 1
   do i2=2,int2%n2_0
      power2(i2) = 4*power2(i2-1)
   enddo
!$OMP END SINGLE

!$OMP DO COLLAPSE(2)
   do i2=1,int2%n2_0
      do i1=1,int2%n1_0
         val = 0
         do j=1,i1
            val = val + comb1(j,i1)*P(i1+i2-j-1)
         enddo
         int2%val_0(i1,i2) = val*power2(i2)
      enddo
   enddo
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC)
   do i1=1,int2%n1_1
      comb1(1,i1) = 1
      do j=1,i1-1
         comb1(j+1,i1) = (comb1(j,i1)/j)*((j-i1)*(2*(i1-j)+1))
      enddo
   enddo
!$OMP END DO NOWAIT

!$OMP SINGLE
   power2(1) = -2
   do i2=2,int2%n2_1
      power2(i2) = 4*power2(i2-1)
   enddo
!$OMP END SINGLE

!$OMP DO COLLAPSE(2)
   do i2=1,int2%n2_1
      do i1=1,int2%n1_1
         val = 0
         do j=1,i1
            val = val + comb1(j,i1)*P(i1+i2-j)
         enddo
         int2%val_1(i1,i2) = val*power2(i2)
      enddo
   enddo
!$OMP END DO

!$OMP END PARALLEL

   deallocate(comb1,power2)

end select

deallocate(P)

!!$block
!!$  character(FLD_DEFAULT) :: s
!!$  write(*,*) '-----------------------------------------------------------------'
!!$  do i2=0,n2,2
!!$     do i1=0,n1,2
!!$        call write_fform(int2%val_0(i1/2+1,i2/2+1),s)
!!$        write(*,*) i1,i2,trim(s)
!!$     enddo
!!$  enddo
!!$  write(*,*) '-----------------------------------------------------------------'
!!$  do i2=1,n2,2
!!$     do i1=1,n1,2
!!$        call write_fform(int2%val_1((i1+1)/2,(i2+1)/2),s)
!!$        write(*,*) i1,i2,trim(s)
!!$     enddo
!!$  enddo
!!$end block
!!$stop

end subroutine int2_create

subroutine int2_choose(f,nC,C,alpha_mult,alpha_f,nC_f,C_f)
implicit none
character(2),intent(in) :: f
integer,intent(out) :: nC
type(mp_real),intent(out),allocatable :: C(:)
integer,intent(inout) :: alpha_mult
type(mp_real),intent(in) :: alpha_f
integer,intent(in) :: nC_f
type(mp_real),intent(in) :: C_f(0:)
integer :: i
type(mp_real) :: val

select case(f)
case('1 ')

   nC = 0
   allocate(C(0:nC))

   C(0) = 1

case('f0')

   nC = 0
   allocate(C(0:nC))

   C(0) = C_f(0)

case('f ')

   nC = nC_f
   allocate(C(0:nC))

   do i=0,nC
      C(i) = C_f(i)
   enddo

   alpha_mult = alpha_mult + 1

case('df')

   nC = nC_f + 1
   allocate(C(0:nC))
   do i=0,nC
      C(i) = 0
   enddo

   do i=0,nC-2
      C(i) = C(i) + (i+1)*C_f(i+1)
   enddo

   val = -2*alpha_f
   do i=1,nC
      C(i) = C(i) + val*C_f(i-1)
   enddo

   alpha_mult = alpha_mult + 1

case default
   write(*,'(2a)') 'Unknown f in int2_choose: ',f
   stop
end select

end subroutine int2_choose

subroutine int2_free(int2)
implicit none
type(int2Data) :: int2

deallocate(int2%val_0,int2%val_1)

end subroutine int2_free

subroutine int3_create(int3_type,int3,n1,n2,n3,alpha,nC,C)
implicit none
character(3),intent(in) :: int3_type
type(int3Data) :: int3
integer,intent(in) :: n1,n2,n3
type(mp_real),intent(in) :: alpha
integer,intent(in) :: nC
type(mp_real),intent(in) :: C(0:)
integer :: nTOT,nHALF,nUV,nMAX,nC2
integer :: mTOT,mHALF,mUV,m1,m2,m3,j1,j2,j3
integer :: ai,aj,ai2,aj2,aij,ad,i,j,ij,iC,jC
integer :: m1p,m2p,m3p,j1p,j2p,j3p
integer :: W
type(mp_real) :: Zi,Z1,Z2,Zs,Zm1,Zm2,Zms
type(mp_real) :: A_00,A_10,A_01,A_11,B_00,B_10,B_01,B_11
type(mp_real) :: val,tmp
type(mp_real),allocatable :: val_000(:,:,:),mul_000(:,:,:)
type(mp_real),allocatable :: val_110(:,:,:),mul_110(:,:,:)
type(mp_real),allocatable :: val_101(:,:,:),mul_101(:,:,:)
type(mp_real),allocatable :: val_011(:,:,:),mul_011(:,:,:)
type(mp_real),allocatable :: F(:,:),F_0(:,:),F_1(:,:)
type(mp_real),allocatable :: CC(:,:),FC_0(:,:),FC_1(:,:)
type(mp_real),allocatable :: powIJ(:),powIk(:),powUV(:),binom(:,:),vec(:)
type(mp_real),allocatable :: fact_0(:),fact_1(:)

int3%n1_0 = n1/2+1
int3%n2_0 = n2/2+1
int3%n3_0 = n3/2+1

int3%n1_1 = (n1+1)/2
int3%n2_1 = (n2+1)/2
int3%n3_1 = (n3+1)/2

allocate(&
     int3%val_000(int3%n1_0,int3%n2_0,int3%n3_0),&
     int3%val_110(int3%n1_1,int3%n2_1,int3%n3_0),&
     int3%val_101(int3%n1_1,int3%n2_0,int3%n3_1),&
     int3%val_011(int3%n1_0,int3%n2_1,int3%n3_1))

allocate(&
     val_000(int3%n1_0,int3%n2_0,int3%n3_0),&
     val_110(int3%n1_1,int3%n2_1,int3%n3_0),&
     val_101(int3%n1_1,int3%n2_0,int3%n3_1),&
     val_011(int3%n1_0,int3%n2_1,int3%n3_1))

allocate(&
     mul_000(int3%n1_0,int3%n2_0,int3%n3_0),&
     mul_110(int3%n1_1,int3%n2_1,int3%n3_0),&
     mul_101(int3%n1_1,int3%n2_0,int3%n3_1),&
     mul_011(int3%n1_0,int3%n2_1,int3%n3_1))

nTOT  = n1 + n2 + n3
nHALF = nTOT/2
nUV   = n2 + n3
nMAX  = max(n1,n2,n3)
nC2   = nC/2

select case(int3_type)
case('222')
   W = 3
case('422')
   W = 4
case default
   write(*,'(2a)') 'Unknown int3_type in int3_create: ',int3_type
   stop
end select

allocate(F(0:nTOT,0:nHALF))

allocate(CC(0:nC,0:nC))
allocate(&
     F_0(0:nHALF,0:nHALF),&
     F_1(0:nHALF,0:nHALF))
allocate(&
     FC_0(0:nHALF+nC2  ,0:nHALF+nC2),&
     FC_1(0:nHALF+nC2+1,0:nHALF+nC2))

!$OMP PARALLEL

!$OMP DO COLLAPSE(2)
do jC=0,nC
   do iC=0,nC
      CC(iC,jC) = C(iC)*C(jC)
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do j=0,nHALF
   do i=0,nHALF
      F_0(i,j) = 0
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2)
do j=0,nHALF
   do i=0,nHALF
      F_1(i,j) = 0
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

Zi = (W-1) + W*alpha
Z1 = 1/Zi
Z2 = Z1**2
Zs = sqrt(Z1)
tmp = 1/(1-Z2)
Zm1 = Z1*tmp
Zm2 = Z2*tmp
Zms = Z1*sqrt(tmp)

B_00 = Zms*Z1
B_01 = Zm2*Zs

tmp = atan(Zms)
A_00 = tmp
A_11 = tmp*Zms
B_11 = (A_11 + Z2)*Zm1

tmp = asinh(Zms)
A_01 = tmp*Zs
A_10 = A_01
B_10 = A_01*Z1 + B_01

val = mppi()
tmp = val*2
A_00 = A_00*tmp
A_10 = A_10*tmp
A_01 = A_01*tmp
tmp = val*W
B_00 = B_00*tmp
B_10 = B_10*tmp
B_01 = B_01*tmp
tmp = sqrt(val/W)
A_00 = A_00*tmp
A_11 = A_11*tmp*(2*W)
B_00 = B_00*tmp
B_11 = B_11*tmp*W**2

Zm1 = (Zm1*W)/2
Zm2 = (Zm2*W**2)/2
Zms = W; Zms = Zms/2

do aj=0,1
   do ai=0,1
      ai2 = ai+ai-1
      aj2 = aj+aj-1
      aij = ai+aj-1
      ad  = ai-aj

      select case(10*ai+aj)
      case(00)
         FC_1(0,0) = A_00
         FC_1(1,0) = B_00
         val       = B_00
      case(10)
         FC_1(0,0) = A_10
         FC_1(1,0) = B_10
         val       = B_01
      case(01)
         FC_1(0,0) = A_01
         FC_1(1,0) = B_01
         val       = B_10
      case(11)
         FC_1(0,0) = A_11
         FC_1(1,0) = B_11
         val       = B_11
      end select

      do j=1,nHALF+nC2
         tmp = FC_1(0,j-1)*Zm2*((1-2*j)*(j+aj-1))
         FC_1(0,j) = val
         val       = val*Zm1
         FC_1(1,j) = tmp + val*(2*j+aij)
         val       = tmp + val*(4*j+aj2 - (2*j-ad)*Z2)
      enddo

!$OMP PARALLEL PRIVATE(i,j)

!$OMP DO SCHEDULE(DYNAMIC)
      do j=0,nHALF+nC2
         do i=1,min(nHALF+2*nC2-j,nHALF+nC2)
            FC_1(i+1,j) = &
                 FC_1(i-1,j)*Zm2*((1-2*i)*(i+ai-1)) + &
                 FC_1(i  ,j)*Zm1*(4*i+ai2 - (2*(i-j)+ad)*Z2)
         enddo
      enddo
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC)
      do j=0,nHALF+nC2
         do i=0,min(nHALF+2*nC2-j,nHALF+nC2)
            FC_0(i,j) = FC_1(i+1,j)*Zi - FC_1(i,j)*Zms*(2*i+ai)
         enddo
      enddo
!$OMP END DO

!$OMP DO
      do ij=0,((nHALF+1)*(nHALF+2))/2-1
         i = mod(ij,nHALF+1)
         j = ij/(nHALF+1)
         if(i+j>nHALF) then
            i = nHALF - i
            j = nHALF - j + 1
         endif
         do jC=aj,nC,2
            do iC=ai,nC,2
               F_0(i,j) = F_0(i,j) &
                    + FC_0(i+(iC-ai)/2,j+(jC-aj)/2)*CC(iC,jC)
            enddo
         enddo
      enddo
!$OMP END DO NOWAIT

!$OMP DO
      do ij=0,(nHALF*(nHALF+1))/2-1
         i = mod(ij,nHALF)+1
         j = ij/nHALF+1
         if(i+j>nHALF+1) then
            i = nHALF - i + 1
            j = nHALF - j + 2
         endif
         do jC=aj,nC,2
            do iC=ai,nC,2
               F_1(i,j) = F_1(i,j) &
                    + FC_1(i+(iC-ai)/2,j+(jC-aj)/2)*CC(iC,jC)
            enddo
         enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

   enddo
enddo

!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
do j=nHALF,0,-1
   do i=0,2*j
      if(mod(i,2)==0) then
         F(i,j) = F_0(i/2,j-i/2)
      else
         F(i,j) = F_1((i+1)/2,j-(i-1)/2)
      endif
   enddo
enddo
!$OMP END PARALLEL DO

deallocate(FC_0,FC_1)
deallocate(F_0,F_1)
deallocate(CC)

allocate(powIJ(0:nMAX),powIk(0:nMAX),powUV(0:nUV),binom(0:nUV+1,0:nUV+1))

!$OMP PARALLEL PRIVATE(vec,val) &
!$OMP PRIVATE(m1p,m2p,m3p,mTOT,mHALF,mUV) &
!$OMP PRIVATE(j1p,j2p,j3p)

!$OMP CRITICAL(int3_vec)
allocate(vec(0:nTOT))
!$OMP END CRITICAL(int3_vec)

!$OMP SECTIONS
!$OMP SECTION
val = 2; val = val/W
powIJ(0) = 1
do i=1,nMAX
   powIJ(i) = powIJ(i-1)*val/i
enddo
!$OMP SECTION
val = 1-W; val = val/W
powIk(0) = 1
do i=1,nMAX
   powIk(i) = powIk(i-1)*val/i
enddo
!$OMP SECTION
val = -W
powUV(0) = 1
do i=1,nUV
   powUV(i) = powUV(i-1)*val
enddo
!$OMP SECTION
binom(0,0) = 1
!$OMP END SECTIONS NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do j=1,nUV+1
   binom(0,j) = 1
   do i=1,j/2
      binom(i,j) = (binom(i-1,j)*(j-i+1))/i
   enddo
   do i=0,(j-1)/2
      binom(j-i,j) = binom(i,j)
   enddo
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
do m3=0,n3,2
   do m2=0,n2,2
      m2p = m2/2+1
      m3p = m3/2+1

      mUV = m2 + m3

      val = powIJ(m2)*powIJ(m3)

      do m1=0,n1,2
         m1p = m1/2+1

         mTOT  = mUV + m1
         mHALF = mTOT/2

         if(m1==0) then
            call int3_firstUV(vec,mTOT,m2,m3,powUV,binom)
         else
            call int3_nextUV(vec,mTOT)
         endif

         val_000(m1p,m2p,m3p) = product1(mTOT+1,vec,F(:,mHALF))*powIJ(m1)*val

      enddo

   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
do m3=0,n3,2
   do m2=1,n2,2
      m2p = (m2+1)/2
      m3p = m3/2+1

      mUV = m2 + m3

      val = powIJ(m2)*powIJ(m3)

      do m1=1,n1,2
         m1p = (m1+1)/2

         mTOT  = mUV + m1
         mHALF = mTOT/2

         if(m1==1) then
            call int3_firstUV(vec,mTOT,m2,m3,powUV,binom)
         else
            call int3_nextUV(vec,mTOT)
         endif

         val_110(m1p,m2p,m3p) = product1(mTOT+1,vec,F(:,mHALF))*powIJ(m1)*val

      enddo

   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
do m3=1,n3,2
   do m2=0,n2,2
      m2p = m2/2+1
      m3p = (m3+1)/2

      mUV = m2 + m3

      val = powIJ(m2)*powIJ(m3)

      do m1=1,n1,2
         m1p = (m1+1)/2

         mTOT  = mUV + m1
         mHALF = mTOT/2

         if(m1==1) then
            call int3_firstUV(vec,mTOT,m2,m3,powUV,binom)
         else
            call int3_nextUV(vec,mTOT)
         endif

         val_101(m1p,m2p,m3p) = product1(mTOT+1,vec,F(:,mHALF))*powIJ(m1)*val

      enddo

   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
do m3=1,n3,2
   do m2=1,n2,2
      m2p = (m2+1)/2
      m3p = (m3+1)/2

      mUV = m2 + m3

      val = powIJ(m2)*powIJ(m3)

      do m1=0,n1,2
         m1p = m1/2+1

         mTOT  = mUV + m1
         mHALF = mTOT/2

         if(m1==0) then
            call int3_firstUV(vec,mTOT,m2,m3,powUV,binom)
         else
            call int3_nextUV(vec,mTOT)
         endif

         val_011(m1p,m2p,m3p) = product1(mTOT+1,vec,F(:,mHALF))*powIJ(m1)*val

      enddo

   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
do j3=0,n3,2
   do j2=0,n2,2
      do j1=0,n1,2
         j1p = j1/2+1
         j2p = j2/2+1
         j3p = j3/2+1

         mul_000(j1p,j2p,j3p) = int3_mulJ(j1,j2,j3,powIJ,powIk)

      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
do j3=0,n3,2
   do j2=1,n2,2
      do j1=1,n1,2
         j1p = (j1+1)/2
         j2p = (j2+1)/2
         j3p = j3/2+1

         mul_110(j1p,j2p,j3p) = int3_mulJ(j1,j2,j3,powIJ,powIk)

      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
do j3=1,n3,2
   do j2=0,n2,2
      do j1=1,n1,2
         j1p = (j1+1)/2
         j2p = j2/2+1
         j3p = (j3+1)/2

         mul_101(j1p,j2p,j3p) = int3_mulJ(j1,j2,j3,powIJ,powIk)

      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
do j3=1,n3,2
   do j2=1,n2,2
      do j1=0,n1,2
         j1p = j1/2+1
         j2p = (j2+1)/2
         j3p = (j3+1)/2

         mul_011(j1p,j2p,j3p) = int3_mulJ(j1,j2,j3,powIJ,powIk)

      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP CRITICAL(int3_vec)
deallocate(vec)
!$OMP END CRITICAL(int3_vec)

!$OMP END PARALLEL

deallocate(powIJ,powIk,powUV,binom)

deallocate(F)

allocate(fact_0(nMAX/2+1),fact_1((nMAX+1)/2))

val = 1
fact_0(1) = 1
do i=1,nMAX
   val = val*i
   if(mod(i,2)==0) then
      fact_0(i/2+1) = val
   else
      fact_1((i+1)/2) = val
   endif
enddo

!$OMP PARALLEL PRIVATE(val)

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
do m3p=1,int3%n3_0
   do m2p=1,int3%n2_0
      do m1p=1,int3%n1_0

         val = 0

         do j3p=1,m3p
            do j2p=1,m2p
               do j1p=1,m1p
                  val = val &
                       + val_000(m1p-j1p+1,m2p-j2p+1,m3p-j3p+1) &
                       * mul_000(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         do j3p=1,m3p
            do j2p=1,m2p-1
               do j1p=1,m1p-1
                  val = val &
                       + val_110(m1p-j1p  ,m2p-j2p  ,m3p-j3p+1) &
                       * mul_110(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         do j3p=1,m3p-1
            do j2p=1,m2p
               do j1p=1,m1p-1
                  val = val &
                       + val_101(m1p-j1p  ,m2p-j2p+1,m3p-j3p  ) &
                       * mul_101(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         do j3p=1,m3p-1
            do j2p=1,m2p-1
               do j1p=1,m1p
                  val = val &
                       + val_011(m1p-j1p+1,m2p-j2p  ,m3p-j3p  ) &
                       * mul_011(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         int3%val_000(m1p,m2p,m3p) = val*fact_0(m1p)*fact_0(m2p)*fact_0(m3p)

      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
do m3p=1,int3%n3_0
   do m2p=1,int3%n2_1
      do m1p=1,int3%n1_1

         val = 0

         do j3p=1,m3p
            do j2p=1,m2p
               do j1p=1,m1p
                  val = val &
                       + val_000(m1p-j1p+1,m2p-j2p+1,m3p-j3p+1) &
                       * mul_110(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         do j3p=1,m3p
            do j2p=1,m2p
               do j1p=1,m1p
                  val = val &
                       + val_110(m1p-j1p+1,m2p-j2p+1,m3p-j3p+1) &
                       * mul_000(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         do j3p=1,m3p-1
            do j2p=1,m2p
               do j1p=1,m1p
                  val = val &
                       + val_101(m1p-j1p+1,m2p-j2p+1,m3p-j3p  ) &
                       * mul_011(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         do j3p=1,m3p-1
            do j2p=1,m2p
               do j1p=1,m1p
                  val = val &
                       + val_011(m1p-j1p+1,m2p-j2p+1,m3p-j3p  ) &
                       * mul_101(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         int3%val_110(m1p,m2p,m3p) = val*fact_1(m1p)*fact_1(m2p)*fact_0(m3p)

      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
do m3p=1,int3%n3_1
   do m2p=1,int3%n2_0
      do m1p=1,int3%n1_1

         val = 0

         do j3p=1,m3p
            do j2p=1,m2p
               do j1p=1,m1p
                  val = val &
                       + val_000(m1p-j1p+1,m2p-j2p+1,m3p-j3p+1) &
                       * mul_101(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         do j3p=1,m3p
            do j2p=1,m2p-1
               do j1p=1,m1p
                  val = val &
                       + val_110(m1p-j1p+1,m2p-j2p  ,m3p-j3p+1) &
                       * mul_011(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         do j3p=1,m3p
            do j2p=1,m2p
               do j1p=1,m1p
                  val = val &
                       + val_101(m1p-j1p+1,m2p-j2p+1,m3p-j3p+1) &
                       * mul_000(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         do j3p=1,m3p
            do j2p=1,m2p-1
               do j1p=1,m1p
                  val = val &
                       + val_011(m1p-j1p+1,m2p-j2p  ,m3p-j3p+1) &
                       * mul_110(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         int3%val_101(m1p,m2p,m3p) = val*fact_1(m1p)*fact_0(m2p)*fact_1(m3p)

      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
do m3p=1,int3%n3_1
   do m2p=1,int3%n2_1
      do m1p=1,int3%n1_0

         val = 0

         do j3p=1,m3p
            do j2p=1,m2p
               do j1p=1,m1p
                  val = val &
                       + val_000(m1p-j1p+1,m2p-j2p+1,m3p-j3p+1) &
                       * mul_011(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         do j3p=1,m3p
            do j2p=1,m2p
               do j1p=1,m1p-1
                  val = val &
                       + val_110(m1p-j1p  ,m2p-j2p+1,m3p-j3p+1) &
                       * mul_101(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         do j3p=1,m3p
            do j2p=1,m2p
               do j1p=1,m1p-1
                  val = val &
                       + val_101(m1p-j1p  ,m2p-j2p+1,m3p-j3p+1) &
                       * mul_110(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         do j3p=1,m3p
            do j2p=1,m2p
               do j1p=1,m1p
                  val = val &
                       + val_011(m1p-j1p+1,m2p-j2p+1,m3p-j3p+1) &
                       * mul_000(j1p,j2p,j3p)
               enddo
            enddo
         enddo

         int3%val_011(m1p,m2p,m3p) = val*fact_0(m1p)*fact_1(m2p)*fact_1(m3p)

      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

deallocate(fact_0,fact_1)

deallocate(mul_000,mul_110,mul_101,mul_011)
deallocate(val_000,val_110,val_101,val_011)

!!$block
!!$  integer,parameter :: ss = 70
!!$  character(ss) :: s
!!$  do m3=0,n3
!!$     do m2=0,n2
!!$        do m1=0,n1
!!$           if(mod(m1+m2+m3,2)==0) then
!!$              if(mod(m2+m3,2)==0) then
!!$                 if(mod(m3,2)==0) then
!!$                    call write_fform(&
!!$                         int3%val_000(m1/2+1,m2/2+1,m3/2+1),s,ss,40)
!!$                 else
!!$                    call write_fform(&
!!$                         int3%val_011(m1/2+1,(m2+1)/2,(m3+1)/2),s,ss,40)
!!$                 endif
!!$              else
!!$                 if(mod(m3,2)==0) then
!!$                    call write_fform(&
!!$                         int3%val_110((m1+1)/2,(m2+1)/2,m3/2+1),s,ss,40)
!!$                 else
!!$                    call write_fform(&
!!$                         int3%val_101((m1+1)/2,m2/2+1,(m3+1)/2),s,ss,40)
!!$                 endif
!!$              endif
!!$              write(*,'(3i5,a)') m1,m2,m3,trim(s)
!!$           endif
!!$        enddo
!!$     enddo
!!$  enddo
!!$end block

end subroutine int3_create

subroutine int3_from2_create(f1,f2,int3,n1,n2,n3,alpha,nC,C)
implicit none
character(2),intent(in) :: f1,f2
type(int3Data) :: int3
integer,intent(in) :: n1,n2,n3
type(mp_real),intent(in) :: alpha
integer,intent(in) :: nC
type(mp_real),intent(in) :: C(:)
integer :: i1,i2,i3,i1p,i2p,i3p,v
type(int2Data) :: int2
type(prod2HData) :: prod2H

int3%n1_0 = n1/2+1
int3%n2_0 = n2/2+1
int3%n3_0 = n3/2+1

int3%n1_1 = (n1+1)/2
int3%n2_1 = (n2+1)/2
int3%n3_1 = (n3+1)/2

allocate(&
     int3%val_000(int3%n1_0,int3%n2_0,int3%n3_0),&
     int3%val_110(int3%n1_1,int3%n2_1,int3%n3_0),&
     int3%val_101(int3%n1_1,int3%n2_0,int3%n3_1),&
     int3%val_011(int3%n1_0,int3%n2_1,int3%n3_1))

call int2_create(f1,f2,'42',int2,n1+n2,n3,alpha,nC,C)
call prod2H_create(prod2H,n1,n2)

!$OMP PARALLEL PRIVATE(i1p,i2p,i3p,v)

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
do i3=0,n3,2
   do i2=0,n2,2
      do i1=0,n1,2
         i1p = i1/2+1
         i2p = i2/2+1
         i3p = i3/2+1

         v = (i1+i2)/2+1

         int3%val_000(i1p,i2p,i3p) = &
              product1_unit(v,prod2H%unit(i1,i2),int2%val_0(:,i3p))

      enddo
   enddo
enddo
!$OMP END DO

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
do i3=0,n3,2
   do i2=1,n2,2
      do i1=1,n1,2
         i1p = (i1+1)/2
         i2p = (i2+1)/2
         i3p = i3/2+1

         v = (i1+i2)/2+1

         int3%val_110(i1p,i2p,i3p) = &
              product1_unit(v,prod2H%unit(i1,i2),int2%val_0(:,i3p))

      enddo
   enddo
enddo
!$OMP END DO

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
do i3=1,n3,2
   do i2=0,n2,2
      do i1=1,n1,2
         i1p = (i1+1)/2
         i2p = i2/2+1
         i3p = (i3+1)/2

         v = (i1+i2+1)/2

         int3%val_101(i1p,i2p,i3p) = &
              product1_unit(v,prod2H%unit(i1,i2),int2%val_1(:,i3p))

      enddo
   enddo
enddo
!$OMP END DO

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
do i3=1,n3,2
   do i2=1,n2,2
      do i1=0,n1,2
         i1p = i1/2+1
         i2p = (i2+1)/2
         i3p = (i3+1)/2

         v = (i1+i2+1)/2

         int3%val_011(i1p,i2p,i3p) = &
              product1_unit(v,prod2H%unit(i1,i2),int2%val_1(:,i3p))

      enddo
   enddo
enddo
!$OMP END DO

!$OMP END PARALLEL

call prod2H_free(prod2H)
call int2_free(int2)

!!$block
!!$  integer,parameter :: ss = 70
!!$  character(ss) :: s
!!$  do i3=0,n3
!!$     do i2=0,n2
!!$        do i1=0,n1
!!$           if(mod(i1+i2+i3,2)==0) then
!!$              if(mod(i2+i3,2)==0) then
!!$                 if(mod(i3,2)==0) then
!!$                    call write_fform(&
!!$                         int3%val_000(i1/2+1,i2/2+1,i3/2+1),s,ss,40)
!!$                 else
!!$                    call write_fform(&
!!$                         int3%val_011(i1/2+1,(i2+1)/2,(i3+1)/2),s,ss,40)
!!$                 endif
!!$              else
!!$                 if(mod(i3,2)==0) then
!!$                    call write_fform(&
!!$                         int3%val_110((i1+1)/2,(i2+1)/2,i3/2+1),s,ss,40)
!!$                 else
!!$                    call write_fform(&
!!$                         int3%val_101((i1+1)/2,i2/2+1,(i3+1)/2),s,ss,40)
!!$                 endif
!!$              endif
!!$              write(*,'(3i5,a)') i1,i2,i3,trim(s)
!!$           endif
!!$        enddo
!!$     enddo
!!$  enddo
!!$end block

end subroutine int3_from2_create

subroutine int3_free(int3)
implicit none
type(int3Data) :: int3

deallocate(int3%val_000,int3%val_110,int3%val_101,int3%val_011)

end subroutine int3_free

subroutine int3_firstUV(vec,n,n1,n2,power,binom)
implicit none
type(mp_real),intent(out) :: vec(0:)
integer,intent(in) :: n,n1,n2
type(mp_real),intent(in) :: power(0:),binom(0:,0:)
integer :: k,k1,k2,i
type(mp_real) :: val

do i=0,n
   vec(i) = 0
enddo

do k2=0,n2
   do k1=0,n1
      k = k1 + k2
      val = power(k)*binom(k1,n1)*binom(k2,n2)
      do i=0,n-k
         vec(k1+i) = vec(k1+i) + val*binom(i,n-k)
      enddo
   enddo
enddo

end subroutine int3_firstUV

subroutine int3_nextUV(vec,n)
implicit none
type(mp_real),intent(inout) :: vec(0:)
integer,intent(in) :: n
integer :: i

vec(n-1) = 0
do i=n-1,1,-1
   vec(i) = vec(i) + vec(i-1)
enddo

vec(n) = 0
do i=n,1,-1
   vec(i) = vec(i) + vec(i-1)
enddo

end subroutine int3_nextUV

function int3_mulJ(j1,j2,j3,powIJ,powIk) result(val)
implicit none
type(mp_real) :: val
integer,intent(in) :: j1,j2,j3
type(mp_real),intent(in) :: powIJ(0:),powIk(0:)
integer :: i12,i23,i31,i1k,i2k,i3k
type(mp_real) :: term

val = 0

do i12=0,min(j1,j2)

   do i23=mod(j2-i12,2),min(j2-i12,j3),2

      i2k = (j2-i12-i23)/2

      term = powIJ(i12)*powIJ(i23)*powIk(i2k)

      do i31=mod(j3-i23,2),min(j3-i23,j1-i12),2

         i1k = (j1-i31-i12)/2
         i3k = (j3-i23-i31)/2
         
         val = val + powIJ(i31)*powIk(i1k)*powIk(i3k) * term

      enddo

   enddo

enddo

end function int3_mulJ

subroutine aux_init(aux,n1,n2,nC)
implicit none
type(auxData) :: aux
integer,intent(in) :: n1,n2,nC

aux%n1 = n1
aux%n2 = n2

aux%n1_0 = n1/2+1
aux%n2_0 = n2/2+1

aux%n1_1 = (n1+1)/2
aux%n2_1 = (n2+1)/2

aux%nC   = nC
aux%nC_0 = nC - mod(nC  ,2)
aux%nC_1 = nC - mod(nC+1,2)

allocate(&
     aux%se_00(2,aux%n2_0),&
     aux%se_10(2,aux%n2_0),&
     aux%se_01(2,aux%n2_1),&
     aux%se_11(2,aux%n2_1))

allocate(&
     aux%val_00(aux%n1_0,aux%n2_0),&
     aux%val_10(aux%n1_1,aux%n2_0),&
     aux%val_01(aux%n1_0,aux%n2_1),&
     aux%val_11(aux%n1_1,aux%n2_1))

end subroutine aux_init

subroutine aux_free(aux)
implicit none
type(auxData) :: aux

deallocate(aux%val_00,aux%val_10,aux%val_01,aux%val_11)
deallocate(aux%se_00,aux%se_10,aux%se_01,aux%se_11)

end subroutine aux_free

subroutine aux2H_create(aux2H,vec,prod2H)
implicit none
type(auxData) :: aux2H
type(mp_real),intent(in) :: vec(0:)
type(prod2HData),intent(in) :: prod2H
integer :: n2,nC,nC_0,nC_1
integer :: i2,iC,ip1,ip2

n2 = aux2H%n2

nC   = aux2H%nC
nC_0 = aux2H%nC_0
nC_1 = aux2H%nC_1

!$OMP PARALLEL PRIVATE(ip1,ip2)

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   aux2H%se_00(1,ip2) = max(i2-nC_0,0)/2+1
   aux2H%se_00(2,ip2) = (i2+nC_0)/2+1
   do ip1=aux2H%se_00(1,ip2),aux2H%se_00(2,ip2)
      aux2H%val_00(ip1,ip2) = 0
   enddo
   do iC=0,nC,2
      ip1 = (i2+iC)/2+1
      call accumulate_unit(vec(iC),ip1,prod2H%unit(iC,i2),aux2H%val_00(:,ip2))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   aux2H%se_10(1,ip2) = (max(i2-nC_1,1)+1)/2
   aux2H%se_10(2,ip2) = ((i2+nC_1)+1)/2
   do ip1=aux2H%se_10(1,ip2),aux2H%se_10(2,ip2)
      aux2H%val_10(ip1,ip2) = 0
   enddo
   do iC=1,nC,2
      ip1 = ((i2+iC)+1)/2
      call accumulate_unit(vec(iC),ip1,prod2H%unit(iC,i2),aux2H%val_10(:,ip2))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=1,n2,2
   ip2 = (i2+1)/2
   aux2H%se_01(1,ip2) = max(i2-nC_1,0)/2+1
   aux2H%se_01(2,ip2) = (i2+nC_1)/2+1
   do ip1=aux2H%se_01(1,ip2),aux2H%se_01(2,ip2)
      aux2H%val_01(ip1,ip2) = 0
   enddo
   do iC=1,nC,2
      ip1 = (i2+iC)/2+1
      call accumulate_unit(vec(iC),ip1,prod2H%unit(iC,i2),aux2H%val_01(:,ip2))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=1,n2,2
   ip2 = (i2+1)/2
   aux2H%se_11(1,ip2) = (max(i2-nC_0,1)+1)/2
   aux2H%se_11(2,ip2) = ((i2+nC_0)+1)/2
   do ip1=aux2H%se_11(1,ip2),aux2H%se_11(2,ip2)
      aux2H%val_11(ip1,ip2) = 0
   enddo
   do iC=0,nC,2
      ip1 = ((i2+iC)+1)/2
      call accumulate_unit(vec(iC),ip1,prod2H%unit(iC,i2),aux2H%val_11(:,ip2))
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

end subroutine aux2H_create

subroutine auxC1_create(auxC1,sum_0,sum_1)
implicit none
type(auxData) :: auxC1
type(mp_real),intent(in) :: sum_0(:),sum_1(:)
integer :: n2,nC,nC_0,nC_1
integer :: i2,iC,j,ip1,ip2
type(mp_real) :: tmp

n2 = auxC1%n2

nC   = auxC1%nC
nC_0 = auxC1%nC_0
nC_1 = auxC1%nC_1

!$OMP PARALLEL PRIVATE(ip1,ip2,tmp)

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   auxC1%se_00(1,ip2) = max(i2-nC_0,0)/2+1
   auxC1%se_00(2,ip2) = (i2+nC_0)/2+1
   do ip1=auxC1%se_00(1,ip2),auxC1%se_00(2,ip2)
      auxC1%val_00(ip1,ip2) = 0
   enddo
   do iC=0,nC,2
      ip1 = (i2+iC)/2+1
      tmp = sum_0(iC/2+1)
      auxC1%val_00(ip1,ip2) = auxC1%val_00(ip1,ip2) + tmp
      do j=0,min(i2,iC)-1
         tmp = (tmp/(j+1))*(2*(i2-j)*(iC-j))
         auxC1%val_00(ip1-j-1,ip2) = auxC1%val_00(ip1-j-1,ip2) + tmp
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   auxC1%se_10(1,ip2) = (max(i2-nC_1,1)+1)/2
   auxC1%se_10(2,ip2) = ((i2+nC_1)+1)/2
   do ip1=auxC1%se_10(1,ip2),auxC1%se_10(2,ip2)
      auxC1%val_10(ip1,ip2) = 0
   enddo
   do iC=1,nC,2
      ip1 = ((i2+iC)+1)/2
      tmp = sum_1((iC+1)/2)
      auxC1%val_10(ip1,ip2) = auxC1%val_10(ip1,ip2) + tmp
      do j=0,min(i2,iC)-1
         tmp = (tmp/(j+1))*(2*(i2-j)*(iC-j))
         auxC1%val_10(ip1-j-1,ip2) = auxC1%val_10(ip1-j-1,ip2) + tmp
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=1,n2,2
   ip2 = (i2+1)/2
   auxC1%se_01(1,ip2) = max(i2-nC_1,0)/2+1
   auxC1%se_01(2,ip2) = (i2+nC_1)/2+1
   do ip1=auxC1%se_01(1,ip2),auxC1%se_01(2,ip2)
      auxC1%val_01(ip1,ip2) = 0
   enddo
   do iC=1,nC,2
      ip1 = (i2+iC)/2+1
      tmp = sum_1((iC+1)/2)
      auxC1%val_01(ip1,ip2) = auxC1%val_01(ip1,ip2) + tmp
      do j=0,min(i2,iC)-1
         tmp = (tmp/(j+1))*(2*(i2-j)*(iC-j))
         auxC1%val_01(ip1-j-1,ip2) = auxC1%val_01(ip1-j-1,ip2) + tmp
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=1,n2,2
   ip2 = (i2+1)/2
   auxC1%se_11(1,ip2) = (max(i2-nC_0,1)+1)/2
   auxC1%se_11(2,ip2) = ((i2+nC_0)+1)/2
   do ip1=auxC1%se_11(1,ip2),auxC1%se_11(2,ip2)
      auxC1%val_11(ip1,ip2) = 0
   enddo
   do iC=0,nC,2
      ip1 = ((i2+iC)+1)/2
      tmp = sum_0(iC/2+1)
      auxC1%val_11(ip1,ip2) = auxC1%val_11(ip1,ip2) + tmp
      do j=0,min(i2,iC)-1
         tmp = (tmp/(j+1))*(2*(i2-j)*(iC-j))
         auxC1%val_11(ip1-j-1,ip2) = auxC1%val_11(ip1-j-1,ip2) + tmp
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

end subroutine auxC1_create

subroutine auxC2_create(auxC2,sum_0,sum_1,aux)
implicit none
type(auxData) :: auxC2
type(mp_real),intent(in) :: sum_0(:),sum_1(:)
type(auxData),intent(in) :: aux
integer :: n2,nC,nC_0,nC_1
integer :: i2,iC1,iC2,j,ip1,ip2
type(mp_real) :: tmp

n2 = auxC2%n2

nC   = auxC2%nC
nC_0 = auxC2%nC_0
nC_1 = auxC2%nC_1

!$OMP PARALLEL PRIVATE(ip1,ip2,tmp)

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   auxC2%se_00(1,ip2) = 1
   auxC2%se_00(2,ip2) = &
        max(nC_0 + 2*(aux%se_00(2,ip2)-1),nC_1 + 2*aux%se_10(2,ip2)-1)/2+1
   do ip1=auxC2%se_00(1,ip2),auxC2%se_00(2,ip2)
      auxC2%val_00(ip1,ip2) = 0
   enddo
   do iC2=2*(aux%se_00(1,ip2)-1),2*(aux%se_00(2,ip2)-1),2
      do iC1=0,nC,2
         ip1 = (iC1+iC2)/2+1
         tmp = sum_0(iC1/2+1)*aux%val_00(iC2/2+1,ip2)
         auxC2%val_00(ip1,ip2) = auxC2%val_00(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            auxC2%val_00(ip1-j-1,ip2) = auxC2%val_00(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
   do iC2=2*aux%se_10(1,ip2)-1,2*aux%se_10(2,ip2)-1,2
      do iC1=1,nC,2
         ip1 = (iC1+iC2)/2+1
         tmp = sum_1((iC1+1)/2)*aux%val_10((iC2+1)/2,ip2)
         auxC2%val_00(ip1,ip2) = auxC2%val_00(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            auxC2%val_00(ip1-j-1,ip2) = auxC2%val_00(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=0,n2,2
   ip2 = i2/2+1
   auxC2%se_10(1,ip2) = 1
   auxC2%se_10(2,ip2) = &
        (max(nC_1 + 2*(aux%se_00(2,ip2)-1),nC_0 + 2*aux%se_10(2,ip2)-1)+1)/2
   do ip1=auxC2%se_10(1,ip2),auxC2%se_10(2,ip2)
      auxC2%val_10(ip1,ip2) = 0
   enddo
   do iC2=2*(aux%se_00(1,ip2)-1),2*(aux%se_00(2,ip2)-1),2
      do iC1=1,nC,2
         ip1 = (iC1+iC2+1)/2
         tmp = sum_1((iC1+1)/2)*aux%val_00(iC2/2+1,ip2)
         auxC2%val_10(ip1,ip2) = auxC2%val_10(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            auxC2%val_10(ip1-j-1,ip2) = auxC2%val_10(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
   do iC2=2*aux%se_10(1,ip2)-1,2*aux%se_10(2,ip2)-1,2
      do iC1=0,nC,2
         ip1 = (iC1+iC2+1)/2
         tmp = sum_0(iC1/2+1)*aux%val_10((iC2+1)/2,ip2)
         auxC2%val_10(ip1,ip2) = auxC2%val_10(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            auxC2%val_10(ip1-j-1,ip2) = auxC2%val_10(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=1,n2,2
   ip2 = (i2+1)/2
   auxC2%se_01(1,ip2) = 1
   auxC2%se_01(2,ip2) = &
        max(nC_0 + 2*(aux%se_01(2,ip2)-1),nC_1 + 2*aux%se_11(2,ip2)-1)/2+1
   do ip1=auxC2%se_01(1,ip2),auxC2%se_01(2,ip2)
      auxC2%val_01(ip1,ip2) = 0
   enddo
   do iC2=2*(aux%se_01(1,ip2)-1),2*(aux%se_01(2,ip2)-1),2
      do iC1=0,nC,2
         ip1 = (iC1+iC2)/2+1
         tmp = sum_0(iC1/2+1)*aux%val_01(iC2/2+1,ip2)
         auxC2%val_01(ip1,ip2) = auxC2%val_01(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            auxC2%val_01(ip1-j-1,ip2) = auxC2%val_01(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
   do iC2=2*aux%se_11(1,ip2)-1,2*aux%se_11(2,ip2)-1,2
      do iC1=1,nC,2
         ip1 = (iC1+iC2)/2+1
         tmp = sum_1((iC1+1)/2)*aux%val_11((iC2+1)/2,ip2)
         auxC2%val_01(ip1,ip2) = auxC2%val_01(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            auxC2%val_01(ip1-j-1,ip2) = auxC2%val_01(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(DYNAMIC)
do i2=1,n2,2
   ip2 = (i2+1)/2
   auxC2%se_11(1,ip2) = 1
   auxC2%se_11(2,ip2) = &
        (max(nC_1 + 2*(aux%se_01(2,ip2)-1),nC_0 + 2*aux%se_11(2,ip2)-1)+1)/2
   do ip1=auxC2%se_11(1,ip2),auxC2%se_11(2,ip2)
      auxC2%val_11(ip1,ip2) = 0
   enddo
   do iC2=2*(aux%se_01(1,ip2)-1),2*(aux%se_01(2,ip2)-1),2
      do iC1=1,nC,2
         ip1 = (iC1+iC2+1)/2
         tmp = sum_1((iC1+1)/2)*aux%val_01(iC2/2+1,ip2)
         auxC2%val_11(ip1,ip2) = auxC2%val_11(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            auxC2%val_11(ip1-j-1,ip2) = auxC2%val_11(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
   do iC2=2*aux%se_11(1,ip2)-1,2*aux%se_11(2,ip2)-1,2
      do iC1=0,nC,2
         ip1 = (iC1+iC2+1)/2
         tmp = sum_0(iC1/2+1)*aux%val_11((iC2+1)/2,ip2)
         auxC2%val_11(ip1,ip2) = auxC2%val_11(ip1,ip2) + tmp
         do j=0,min(iC1,iC2)-1
            tmp = (tmp/(j+1))*(2*(iC1-j)*(iC2-j))
            auxC2%val_11(ip1-j-1,ip2) = auxC2%val_11(ip1-j-1,ip2) + tmp
         enddo
      enddo
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

end subroutine auxC2_create

end module integrals_core
