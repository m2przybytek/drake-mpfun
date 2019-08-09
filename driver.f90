module driver
use iso_fortran_env, only : real64
use mpmodule
use time
use diis
use eproblem
use lproblem
use global
use misc
use integrals
implicit none

private
public calc_SCF,calc_CC

character(1) :: par_sign(0:1) = ['+','-']

type InfoData
type(mp_real) :: energy
type(mp_real) :: orthog(2)
end type InfoData

type,extends(InfoData) :: pairInfoData
logical :: nonzero(2)
end type pairInfoData

type CCpairData
type(mp_real) :: mult
type(mp_real),pointer :: tau(:),Qtau(:),R(:),E(:)
type(pairInfoData),allocatable :: Info(:,:)
type(DecompositionData),allocatable :: LHS
end type CCpairData

contains

subroutine calc_SCF(energSCF,orbE,matC,parC,OOOO,Control)
implicit none
type(mp_real),intent(out) :: energSCF
type(mp_real),intent(out) :: orbE(:)
type(mp_real),intent(out) :: matC(:,:)
integer,intent(out) :: parC(:)
type(mp_real),intent(out) :: OOOO(:,:,:,:)
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu,Twall
integer :: nocc,nbas,n_0,n_1,i_0,i_1,ipar,iter,i,j,k,l
integer :: DIIS_START,DIIS_SIZE
type(DIISData) :: DIIS
logical :: do_DIIS
character(FLD_DEFAULT) :: sval
type(mp_real) :: energSCF_prev,energSCF_diff,energSCF_check,error
type(mp_real) :: energMP2,epair
type(mp_real) :: fact
type(mp_real),allocatable :: scfS(:,:),scfH(:,:)
type(mp_real),allocatable :: scfF(:,:),scfD(:,:),scfP(:,:)
type(mp_real),allocatable :: work(:,:)
type(scfJauxData) :: scfJaux
type(mp_real),allocatable :: eval_0(:),evec_0(:,:),tmpF_0(:,:),tmpS_0(:,:)
type(mp_real),allocatable :: eval_1(:),evec_1(:,:),tmpF_1(:,:),tmpS_1(:,:)
type(mp_real),allocatable :: PPOO(:,:,:,:)

!!$block
!!$  integer,parameter :: ns = 5
!!$  allocate(scfS(ns,ns),scfH(ns,ns),eval_0(ns),evec_0(ns,ns))
!!$  do j=1,ns
!!$     do i=1,ns
!!$        fact = i+j
!!$        scfS(i,j) = 1/fact
!!$        fact = i*j
!!$        scfH(i,j) = cos(fact)
!!$     enddo
!!$  enddo
!!$  call symU_diagonalize(Control%eps,ns,eval_0,evec_0,scfH,scfS)
!!$  do i=1,ns
!!$     write(*,*) i,dble(eval_0(i))
!!$  enddo
!!$  write(*,*)
!!$  do i=1,ns
!!$     write(*,*) i,dble(evec_0(i,min(6,ns)))
!!$  enddo
!!$end block
!!$stop

call timer('START',Tcpu,Twall)

write(*,'()')
write(*,'(a)') repeat('-',80)
write(*,'(a)') repeat('*',80)
write(*,'(a)') repeat('-',80)
write(*,'()')

nocc = Control%nocc
nbas = Control%nbas
n_0  = (nbas-1)/2+1
n_1  = nbas/2
DIIS_START = Control%SCF_DIIS_START
DIIS_SIZE  = min(nbas,Control%SCF_DIIS_SIZE)
do_DIIS    = (DIIS_SIZE>1)

allocate(scfS(nbas,nbas),scfH(nbas,nbas))
allocate(scfF(nbas,nbas),scfD(nbas,nbas),scfP(nbas,nbas))

call integrals_scfSH(scfS,scfH,Control)
call init_scfJaux(scfJaux,Control)

allocate(work(nbas,nbas))

if(do_DIIS) call DIIS%init(Control%eps,nbas,nbas,DIIS_SIZE)
allocate(eval_0(n_0),evec_0(n_0,n_0),tmpF_0(n_0,n_0),tmpS_0(n_0,n_0))
allocate(eval_1(n_1),evec_1(n_1,n_1),tmpF_1(n_1,n_1),tmpS_1(n_1,n_1))

write(*,'()')
write(*,'(2x,a,34x,a,25x,a,6x,a)') 'iter','energy','dE','error'

energSCF = 0
call copy_mat(nbas,nbas,scfH,scfF)

do iter=0,Control%SCF_MAXITER

   tmpF_0(:,:) = scfF(1:nbas:2,1:nbas:2)
   tmpS_0(:,:) = scfS(1:nbas:2,1:nbas:2)
   call symU_diagonalize(Control%eps,n_0,eval_0,evec_0,tmpF_0,tmpS_0)
   tmpF_1(:,:) = scfF(2:nbas:2,2:nbas:2)
   tmpS_1(:,:) = scfS(2:nbas:2,2:nbas:2)
   call symU_diagonalize(Control%eps,n_1,eval_1,evec_1,tmpF_1,tmpS_1)

   call zero_mat(nbas,nbas,matC)

   i_0 = 0
   i_1 = 0
   do i=1,nbas
      if(i_0<n_0) then
         if(i_1<n_1) then
            if(eval_0(i_0+1)<=eval_1(i_1+1)) then
               ipar = 0
            else
               ipar = 1
            endif
         else
            ipar = 0
         endif
      else
         ipar = 1
      endif
      select case(ipar)
      case(0)
         i_0 = i_0 + 1
         orbE(i) = eval_0(i_0)
         matC(1:nbas:2,i) = evec_0(:,i_0)
         parC(i) = 0
      case(1)
         i_1 = i_1 + 1
         orbE(i) = eval_1(i_1)
         matC(2:nbas:2,i) = evec_1(:,i_1)
         parC(i) = 1
      end select
   enddo

   call zero_mat(nbas,nbas,scfD)
   do i=1,nocc
      fact = 1
      call add_outer(fact,nbas,nbas,matC(:,i),matC(:,i),scfD)
   enddo
   call integrals_scfJ(work,scfD,scfJaux,Control)

   energSCF_prev = energSCF
   energSCF = 2*dot_mat(nbas,nbas,scfH,scfD) + dot_mat(nbas,nbas,work,scfD)
   energSCF_diff = energSCF - energSCF_prev

   call copy_mat(nbas,nbas,scfH,scfF)
   fact = 1
   call add_mat(fact,nbas,nbas,work,scfF)

   call FDS_SDF(nbas,error,scfP,scfF,scfD,scfS,work)

   call write_fform(energSCF,sval)
   write(*,'(i5,a,2x,2es9.1)') iter,&
        trim(sval),&
        dble(energSCF_diff),&
        dble(error)

   if(iter>0.and.abs(error)<Control%SCF_THR) exit

   if(iter>DIIS_START.and.do_DIIS) call DIIS%step(scfF,scfP)

enddo

write(*,'()')

deallocate(eval_1,evec_1,tmpF_1,tmpS_1)
deallocate(eval_0,evec_0,tmpF_0,tmpS_0)
if(do_DIIS) call DIIS%free

allocate(PPOO(nbas,nbas,nocc,nocc))

do j=1,nocc
   do i=1,nocc
      call outer(nbas,nbas,matC(:,i),matC(:,j),scfD)
      call integrals_scfJ(PPOO(:,:,i,j),scfD,scfJaux,Control)
      call fullMO_transform(nbas,matC,PPOO(:,:,i,j),work)
   enddo
enddo

energSCF_check = 0
do i=1,nocc
   energSCF_check = energSCF_check + 2*orbE(i)
enddo
do j=1,nocc
   do i=1,nocc
      energSCF_check = energSCF_check - PPOO(i,i,j,j)
   enddo
enddo
call write_fform(energSCF,sval)
write(*,'(2a)') 'SCF energy: ',trim(sval)
call write_fform(energSCF_check,sval)
write(*,'(2a)') '     check: ',trim(sval)
write(*,'()')

write(*,'(a)') 'Orbital energies'
do i=1,nocc
   call write_fform(orbE(i),sval)
   write(*,'(i5,2x,a,4x,a)') i,par_sign(parC(i)),trim(sval)
enddo
write(*,'()')

write(*,'(a)') 'Pair energies (standard)'
energMP2 = 0
do j=1,nocc
   do i=1,j
      epair = 0
      do l=nocc+1,nbas
         do k=nocc+1,nbas
            epair = epair &
                 - PPOO(k,l,i,j)**2 / (orbE(k)+orbE(l)-orbE(i)-orbE(j))
         enddo
      enddo
      if(i/=j) epair = 2*epair
      energMP2 = energMP2 + epair
      call write_fform(epair,sval)
      write(*,'(2i5,2x,a)') i,j,trim(sval)
   enddo
enddo

write(*,'()')
call write_fform(energMP2,sval)
write(*,'(a)') 'MP2 energy (standard)'
write(*,'(a)') trim(sval)
write(*,'()')
call write_fform(energSCF+energMP2,sval)
write(*,'(a)') 'Total MP2 energy (standard)'
write(*,'(a)') trim(sval)
write(*,'()')

do j=1,nocc
   do i=1,nocc
      do l=1,nocc
         do k=1,nocc
            OOOO(k,l,i,j) = PPOO(k,i,l,j)
         enddo
      enddo
   enddo
enddo

deallocate(PPOO)

deallocate(work)

call free_scfJaux(scfJaux)

deallocate(scfF,scfD,scfP)
deallocate(scfS,scfH)

call timer('SCF',Tcpu,Twall)

end subroutine calc_SCF

subroutine calc_CC(energSCF,orbE,matC,parC,OOOO,Control)
implicit none
type(mp_real),intent(in) :: energSCF
type(mp_real),intent(in) :: orbE(:)
type(mp_real),intent(in) :: matC(:,:)
integer,intent(in) :: parC(:)
type(mp_real),intent(in) :: OOOO(:,:,:,:)
type(ControlData),intent(in) :: Control
real(real64) :: Tcpu_TOT,Twall_TOT,Tcpu,Twall
integer :: stage
integer :: norb,npair
integer :: iorb,jorb,korb
integer :: ipair,jpair
integer :: neta,ieta,iter,i
type(mp_real) :: energy_this,energy_prev,edelta_this,edelta_prev,error
type(mp_real) :: inv_log10,fact
character(FLD_DEFAULT) :: sval
integer :: DIIS_START,DIIS_SIZE
logical :: do_DIIS
type(DIISData) :: DIIS
type(mp_real),allocatable :: eta(:)
type(mp_real),allocatable :: matS_S(:,:),matS_T(:,:)
type(mp_real),allocatable :: matF_S(:,:),matF_T(:,:)
type(mp_real),allocatable :: matL_S(:,:),matL_T(:,:)
type(mp_real),allocatable :: matJ_S(:,:,:,:),matJ_T(:,:,:,:)
type(mp_real),allocatable :: matM_S(:,:,:,:),matM_T(:,:,:,:)
type(mp_real),allocatable :: matP1_S(:,:),matP1_T(:,:)
type(mp_real),allocatable :: matPe_S(:,:),matPe_T(:,:)
type(mp_real),allocatable :: matP2_S(:,:),matP2_T(:,:)
type(mp_real),allocatable :: vecG(:,:,:)
type(mp_real),allocatable :: vecP(:,:,:),vec1(:,:,:),vec0(:,:,:)
type(mp_real),allocatable :: pair_energy(:,:,:,:)
type(mp_real),allocatable :: TMP(:,:)
type(mp_real),allocatable,target :: DATA_tau(:,:),DATA_Qtau(:,:)
type(mp_real),allocatable,target :: DATA_R(:,:),DATA_E(:,:)
type(CCpairData),allocatable :: CCpairs(:,:)
type(TripletData) :: Triplet
type(DecompositionData) :: project_S,project_T
type(InfoData),allocatable :: fullInfo(:,:)
integer,allocatable :: final_iter(:)

call timer('START',Tcpu_TOT,Twall_TOT)
Tcpu  = Tcpu_TOT
Twall = Twall_TOT

write(*,'()')
write(*,'(a)') repeat('-',80)
write(*,'(a)') repeat('*',80)
write(*,'(a)') repeat('-',80)
write(*,'()')

select case(trim(Control%calc_type))
case('FCCD')
   stage = 3
case('LCCD')
   stage = 2
case('CID')
   stage = 1
case default
   stage = 0
end select

norb  = Control%nocc
npair = Control%npair
DIIS_START = max(1,Control%CC_DIIS_START)
DIIS_SIZE  = min(npair*norb**2,Control%CC_DIIS_SIZE)
do_DIIS    = (stage>0.and.DIIS_SIZE>1)

inv_log10 = 10
inv_log10 = 1/log(inv_log10)

neta = max(0,Control%MAX_ETA_POWER - Control%MIN_ETA_POWER + 1) + 1
allocate(eta(neta))
i = 1
do ieta=Control%MAX_ETA_POWER,Control%MIN_ETA_POWER,-1
   eta(i) = 10
   eta(i) = eta(i)**ieta
   i = i + 1
enddo
eta(i) = 0

allocate(TMP(npair,npair))

!-------------------------------------------------------------------------------

if(stage>0) then

   allocate(matL_S(npair,npair),matL_T(npair,npair))

   call integrals_L00(matL_S,Control)

   call integrals_L10(TMP,matC,parC,Control)
   fact = -1
   call add_mat(fact,npair,npair,TMP,matL_S)

   if(Control%extended_L) then

      fact = -1
      call add_matT(fact,npair,npair,TMP,matL_S)

      call integrals_L11(TMP,matC,parC,Control)
      fact = 1
      call add_mat(fact,npair,npair,TMP,matL_S)

   endif

   call zero_mat(npair,npair,matL_T)

endif

allocate(matP1_S(npair,npair),matP1_T(npair,npair))
allocate(matPe_S(npair,npair),matPe_T(npair,npair))
call integrals_P(matP1_S,matP1_T,matPe_S,matPe_T,orbE,matC,parC,Control)

allocate(&
     matJ_S(npair,npair,norb,norb),&
     matJ_T(npair,npair,norb,norb))
allocate(&
     matM_S(npair,npair,norb,norb),&
     matM_T(npair,npair,norb,norb))
call integrals_J(matJ_S,matJ_T,matM_S,matM_T,matC,parC,Control)

allocate(matS_S(npair,npair),matS_T(npair,npair))
allocate(matF_S(npair,npair),matF_T(npair,npair))
call integrals_X('S',matS_S,matS_T,Control)
call integrals_X('H',matF_S,matF_T,Control)

call integrals_Hrest(matF_S,matF_T,matS_S,matS_T,Control)
do iorb=1,norb
   fact = 2
   fact = fact - 1
   call add_mat(fact,npair,npair,matJ_S(:,:,iorb,iorb),matF_S)
enddo
do iorb=1,norb
   fact = 2
   fact = fact - 1
   call add_mat(fact,npair,npair,matJ_T(:,:,iorb,iorb),matF_T)
enddo

allocate(matP2_S(npair,npair),matP2_T(npair,npair))
allocate(vecG(npair,norb,norb))

allocate(vecP(npair,norb,norb))
allocate(vec1(npair,norb,norb))
allocate(vec0(npair,norb,norb))
call integrals_vecP(    vecP,matC,parC,Control)
call integrals_vecX('1',vec1,matC,parC,Control)
call integrals_vecX('0',vec0,matC,parC,Control)

call zero_mat(npair,npair,matP2_S)
call zero_mat(npair,npair,matP2_T)
do jorb=1,norb
   do iorb=1,norb

      call outer(npair,npair,vec0(:,iorb,jorb),vec0(:,iorb,jorb),TMP)
      fact = 1
      call add_mat(fact,npair,npair,TMP,matP2_S)
      if(iorb/=jorb) then
         fact = 1
         call add_mat(fact,npair,npair,TMP,matP2_T)
      endif

      call outer(npair,npair,vec0(:,iorb,jorb),vec0(:,jorb,iorb),TMP)
      fact = 1
      call add_mat(fact,npair,npair,TMP,matP2_S)
      if(iorb/=jorb) then
         fact = -1
         call add_mat(fact,npair,npair,TMP,matP2_T)
      endif

   enddo
enddo

do jpair=1,norb
   do ipair=1,norb

      call copy_vec(npair,vec1(:,ipair,jpair),vecG(:,ipair,jpair))

      fact = -1
      call add_vec(fact,npair,vecP(:,ipair,jpair),vecG(:,ipair,jpair))

      do jorb=1,norb
         do iorb=1,norb
            fact = OOOO(iorb,jorb,ipair,jpair)
            call add_vec(fact,npair,vec0(:,iorb,jorb),vecG(:,ipair,jpair))
         enddo
      enddo

   enddo
enddo

if(stage>0) then

   do jorb=1,norb
      do iorb=1,norb
         fact = 1
         call add_outer(fact,npair,npair,&
              vec0(:,iorb,jorb),vec1(:,iorb,jorb),matL_S)
      enddo
   enddo

   if(Control%extended_L) then
      do jorb=1,norb
         do iorb=1,norb
            fact = -1
            call add_outer(fact,npair,npair,&
                 vec0(:,iorb,jorb),vecP(:,iorb,jorb),matL_S)
            fact = 1
            call add_outer(fact,npair,npair,&
                 vecG(:,iorb,jorb),vec0(:,iorb,jorb),matL_S)
         enddo
      enddo
   endif

   fact = 2
   call scal_mat(fact,npair,npair,matL_S)

endif

deallocate(vec0)
deallocate(vec1)
deallocate(vecP)

write(*,'()')
call timer('integrals',Tcpu,Twall)

!-------------------------------------------------------------------------------

allocate(DATA_tau(npair,norb**2))
allocate(DATA_Qtau(npair,norb**2))
allocate(DATA_R(npair,norb**2))
allocate(DATA_E(npair,norb**2))
if(do_DIIS) call DIIS%init(Control%eps,npair,norb**2,DIIS_SIZE)

allocate(CCpairs(norb,norb))
i = 0
do jpair=1,norb
   do ipair=1,norb
      associate(CCpair => CCpairs(ipair,jpair))

        if(ipair<jpair) then
           CCpair%mult = 3
        else
           CCpair%mult = 1
           if(ipair==jpair) CCpair%mult = CCpair%mult/2
        endif

        i = i + 1
        CCpair%tau  => DATA_tau(:,i)
        CCpair%Qtau => DATA_Qtau(:,i)
        CCpair%R    => DATA_R(:,i)
        CCpair%E    => DATA_E(:,i)

        allocate(CCpair%Info(0:Control%CC_MAXITER,neta))

        if(stage>1) allocate(CCpair%LHS)

      end associate
   enddo
enddo

allocate(fullInfo(0:Control%CC_MAXITER,neta))
allocate(final_iter(neta))

call Triplet%init(npair,Control%omega)

call copy_mat(npair,npair,matS_S,TMP)
call project_S%init(npair,TMP)
call copy_mat(npair,npair,matS_T,TMP)
call Triplet%shrink(TMP)
call project_T%init(Triplet%npair,TMP)

if(stage>2) allocate(pair_energy(norb,norb,norb,norb))

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

do ieta=1,neta

   write(*,'()')
   write(*,'(a)') repeat('-',80)
   write(*,'()')
   write(*,'(5x,a,es9.1)') 'ETA = ',dble(eta(ieta))
   write(*,'()')

   if(stage>1) then
      do jpair=1,norb
         do ipair=1,jpair

            associate(CCpair => CCpairs(jpair,ipair))

              call MWO(npair,TMP,norb,orbE,jpair,ipair,eta(ieta),&
                   matF_S,matS_S,matP1_S,matPe_S,matP2_S)
              call CCpair%LHS%init(npair,TMP)

            end associate

            if(ipair==jpair) cycle

            associate(CCpair => CCpairs(ipair,jpair))

              call MWO(npair,TMP,norb,orbE,ipair,jpair,eta(ieta),&
                   matF_T,matS_T,matP1_T,matPe_T,matP2_T)
              call Triplet%shrink(TMP)
              call CCpair%LHS%init(Triplet%npair,TMP)

            end associate

         enddo
      enddo
   endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

   do iter=0,Control%CC_MAXITER

!-------------------------------------------------------------------------------

      if(stage>2.and.iter>0) then
         do jpair=1,norb
            do ipair=1,jpair

               associate(CCpair => CCpairs(jpair,ipair),&
                    elms => pair_energy(:,:,jpair,ipair))
                 do jorb=1,norb
                    do iorb=1,jorb
                       fact = &
                            dot_vec(npair,vecG(:,iorb,jorb),CCpair%tau) + &
                            dot_vec(npair,vecG(:,jorb,iorb),CCpair%tau)
                       elms(iorb,jorb) = fact
                       elms(jorb,iorb) = fact
                    enddo
                 enddo
               end associate

               if(ipair==jpair) cycle

               associate(CCpair => CCpairs(ipair,jpair),&
                    elms => pair_energy(:,:,ipair,jpair))
                 do iorb=1,norb
                    elms(iorb,iorb) = 0
                 enddo
                 do jorb=2,norb
                    do iorb=1,jorb-1
                       fact = &
                            dot_vec(npair,vecG(:,iorb,jorb),CCpair%tau) - &
                            dot_vec(npair,vecG(:,jorb,iorb),CCpair%tau)
                       elms(iorb,jorb) =  fact
                       elms(jorb,iorb) = -fact
                    enddo
                 enddo
               end associate

            enddo
         enddo
      endif

!-------------------------------------------------------------------------------

      do jpair=1,norb
         do ipair=1,jpair

            associate(RHS => CCpairs(jpair,ipair)%R)

              call zero_vec(npair,RHS)

              fact = -1
              call add_vec(fact,npair,vecG(:,ipair,jpair),RHS)
              fact = -1
              call add_vec(fact,npair,vecG(:,jpair,ipair),RHS)

              if(stage>0.and.iter>0) then

                 if(Control%extended_L) then
                    fact = -1
                    call add_matvec(fact,npair,npair,&
                         matL_S,CCpairs(jpair,ipair)%tau,RHS)
                 else
                    fact = -1
                    call add_matvec(fact,npair,npair,&
                         matL_S,CCpairs(jpair,ipair)%Qtau,RHS)
                 endif

                 do jorb=1,norb
                    do iorb=1,jorb
                       fact = &
                            - OOOO(iorb,jorb,ipair,jpair) &
                            - OOOO(jorb,iorb,ipair,jpair)
                       if(iorb==jorb) fact = fact/2
                       call add_matvec(fact,npair,npair,&
                            matS_S,CCpairs(jorb,iorb)%Qtau,RHS)
                    enddo
                 enddo

                 do korb=1,norb
                    fact = 2 - 1
                    fact = fact/2
                    call add_matvec(fact,npair,npair,&
                         matJ_S(:,:,korb,jpair),&
                         CCpairs(max(ipair,korb),min(ipair,korb))%Qtau,RHS)
                    if(korb/=ipair) then
                       fact = 3
                       fact = fact/2
                       if(ipair>korb) fact = -fact
                       call add_matvec(fact,npair,npair,&
                            matM_S(:,:,korb,jpair),&
                            CCpairs(min(ipair,korb),max(ipair,korb))%Qtau,RHS)
                    endif
                    fact = 2 - 1
                    fact = fact/2
                    call add_matvec(fact,npair,npair,&
                         matJ_S(:,:,korb,ipair),&
                         CCpairs(max(korb,jpair),min(korb,jpair))%Qtau,RHS)
                    if(korb/=jpair) then
                       fact = -3
                       fact = fact/2
                       if(korb>jpair) fact = -fact
                       call add_matvec(fact,npair,npair,&
                            matM_S(:,:,korb,ipair),&
                            CCpairs(min(korb,jpair),max(korb,jpair))%Qtau,RHS)
                    endif
                 enddo

                 if(stage>2) then

                    do jorb=1,norb
                       do iorb=1,jorb
                          fact = -pair_energy(iorb,jorb,jpair,ipair)
                          if(iorb==jorb) fact = fact/2
                          call add_matvec(fact,npair,npair,&
                               matS_S,CCpairs(jorb,iorb)%Qtau,RHS)
                       enddo
                    enddo

                    do jorb=1,norb
                       do iorb=1,norb
                          fact = pair_energy(iorb,jorb,&
                               max(jpair,jorb),min(jpair,jorb))/2
                          call add_matvec(fact,npair,npair,&
                               matS_S,CCpairs(&
                               max(ipair,iorb),min(ipair,iorb))%Qtau,RHS)
                          fact = pair_energy(iorb,jorb,&
                               max(ipair,jorb),min(ipair,jorb))/2
                          call add_matvec(fact,npair,npair,&
                               matS_S,CCpairs(&
                               max(iorb,jpair),min(iorb,jpair))%Qtau,RHS)
                          if(iorb==jorb) cycle
                          if(jpair/=jorb) then
                             fact = (pair_energy(iorb,jorb,&
                                  min(jpair,jorb),max(jpair,jorb))*3)/2
                             if(jpair>jorb) fact = -fact
                             call add_matvec(fact,npair,npair,&
                                  matS_S,CCpairs(&
                                  max(ipair,iorb),min(ipair,iorb))%Qtau,RHS)
                          endif
                          if(ipair/=jorb) then
                             fact = (pair_energy(iorb,jorb,&
                                  min(ipair,jorb),max(ipair,jorb))*3)/2
                             if(ipair>jorb) fact = -fact
                             call add_matvec(fact,npair,npair,&
                                  matS_S,CCpairs(&
                                  max(iorb,jpair),min(iorb,jpair))%Qtau,RHS)
                          endif
                       enddo
                    enddo

                 endif

              endif

            end associate

            if(ipair==jpair) cycle

            associate(RHS => CCpairs(ipair,jpair)%R)

              call zero_vec(npair,RHS)

              fact = -1
              call add_vec(fact,npair,vecG(:,ipair,jpair),RHS)
              fact = 1
              call add_vec(fact,npair,vecG(:,jpair,ipair),RHS)

              if(stage>0.and.iter>0) then

                 if(Control%extended_L) then
                    fact = -1
                    call add_matvec(fact,npair,npair,&
                         matL_T,CCpairs(ipair,jpair)%tau,RHS)
                 else
                    fact = -1
                    call add_matvec(fact,npair,npair,&
                         matL_T,CCpairs(ipair,jpair)%Qtau,RHS)
                 endif

                 do jorb=2,norb
                    do iorb=1,jorb-1
                       fact = &
                            - OOOO(iorb,jorb,ipair,jpair) &
                            + OOOO(jorb,iorb,ipair,jpair)
                       call add_matvec(fact,npair,npair,&
                            matS_T,CCpairs(iorb,jorb)%Qtau,RHS)
                    enddo
                 enddo

                 do korb=1,norb
                    if(korb/=ipair) then
                       fact = 2 - 3
                       fact = fact/2
                       if(ipair>korb) fact = -fact
                       call add_matvec(fact,npair,npair,&
                            matJ_T(:,:,korb,jpair),&
                            CCpairs(min(ipair,korb),max(ipair,korb))%Qtau,RHS)
                    endif
                    fact = 1
                    fact = fact/2
                    call add_matvec(fact,npair,npair,&
                         matM_T(:,:,korb,jpair),&
                         CCpairs(max(ipair,korb),min(ipair,korb))%Qtau,RHS)
                    if(korb/=jpair) then
                       fact = 2 - 3
                       fact = fact/2
                       if(korb>jpair) fact = -fact
                       call add_matvec(fact,npair,npair,&
                            matJ_T(:,:,korb,ipair),&
                            CCpairs(min(korb,jpair),max(korb,jpair))%Qtau,RHS)
                    endif
                    fact = -1
                    fact = fact/2
                    call add_matvec(fact,npair,npair,&
                         matM_T(:,:,korb,ipair),&
                         CCpairs(max(korb,jpair),min(korb,jpair))%Qtau,RHS)
                 enddo

                 if(stage>2) then

                    do jorb=2,norb
                       do iorb=1,jorb-1
                          fact = -pair_energy(iorb,jorb,ipair,jpair)
                          call add_matvec(fact,npair,npair,&
                               matS_T,CCpairs(iorb,jorb)%Qtau,RHS)
                       enddo
                    enddo

                    do jorb=1,norb
                       do iorb=1,norb
                          if(ipair/=iorb) then
                             fact = pair_energy(iorb,jorb,&
                                  max(jpair,jorb),min(jpair,jorb))/2
                             if(ipair>iorb) fact = -fact
                             call add_matvec(fact,npair,npair,&
                                  matS_T,CCpairs(&
                                  min(ipair,iorb),max(ipair,iorb))%Qtau,RHS)
                          endif
                          if(iorb/=jpair) then
                             fact = pair_energy(iorb,jorb,&
                                  max(ipair,jorb),min(ipair,jorb))/2
                             if(iorb>jpair) fact = -fact
                             call add_matvec(fact,npair,npair,&
                                  matS_T,CCpairs(&
                                  min(iorb,jpair),max(iorb,jpair))%Qtau,RHS)
                          endif
                          if(iorb==jorb) cycle
                          if(ipair/=iorb.and.jpair/=jorb) then
                             fact = (pair_energy(iorb,jorb,&
                                  min(jpair,jorb),max(jpair,jorb))*3)/2
                             if(ipair>iorb) fact = -fact
                             if(jpair>jorb) fact = -fact
                             call add_matvec(fact,npair,npair,&
                                  matS_T,CCpairs(&
                                  min(ipair,iorb),max(ipair,iorb))%Qtau,RHS)
                          endif
                          if(iorb/=jpair.and.ipair/=jorb) then
                             fact = (pair_energy(iorb,jorb,&
                                  min(ipair,jorb),max(ipair,jorb))*3)/2
                             if(iorb>jpair) fact = -fact
                             if(ipair>jorb) fact = -fact
                             call add_matvec(fact,npair,npair,&
                                  matS_T,CCpairs(&
                                  min(iorb,jpair),max(iorb,jpair))%Qtau,RHS)
                          endif
                       enddo
                    enddo

                 endif

              endif

            end associate

         enddo
      enddo

!-------------------------------------------------------------------------------

      do jpair=1,norb
         do ipair=1,jpair

            associate(CCpair => CCpairs(jpair,ipair))

              call copy_vec(npair,CCpair%R,CCpair%tau)
              if(stage>1) then
                 call CCpair%LHS%solve(npair,CCpair%tau)
              else
                 call MWO(npair,TMP,norb,orbE,jpair,ipair,eta(ieta),&
                      matF_S,matS_S,matP1_S,matPe_S,matP2_S)
                 if(stage==1.and.iter>0) then
                    fact = -energy_prev
                    call add_mat(fact,npair,npair,matS_S,TMP)
                 endif
                 call symU_linearsolve(npair,TMP,CCpair%tau)
              endif

              call TQT(npair,TMP,matS_S,matP1_S,matP2_S)
              call matvec(npair,npair,TMP,CCpair%tau,CCpair%Qtau)
              call project_S%solve(npair,CCpair%Qtau)

              associate(Info => CCpair%Info(iter,ieta))
                  
                Info%energy = CCpair%mult*( &
                     dot_vec(npair,CCpair%tau,vecG(:,ipair,jpair)) + &
                     dot_vec(npair,CCpair%tau,vecG(:,jpair,ipair)) )
                call SOcheck(Info%orthog(1),Info%nonzero(1),&
                     npair,CCpair%tau,matS_S,matP1_S,Control%eps,inv_log10)
                call SOcheck(Info%orthog(2),Info%nonzero(2),&
                     npair,CCpair%Qtau,matS_S,matP1_S,Control%eps,inv_log10)

              end associate

            end associate

            if(ipair==jpair) cycle

            associate(CCpair => CCpairs(ipair,jpair))

              call copy_vec(npair,CCpair%R,CCpair%tau)
              call Triplet%shrink(CCpair%tau)
              if(stage>1) then
                 call CCpair%LHS%solve(Triplet%npair,CCpair%tau)
              else
                 call MWO(npair,TMP,norb,orbE,ipair,jpair,eta(ieta),&
                      matF_T,matS_T,matP1_T,matPe_T,matP2_T)
                 if(stage==1.and.iter>0) then
                    fact = -energy_prev
                    call add_mat(fact,npair,npair,matS_T,TMP)
                 endif
                 call Triplet%shrink(TMP)
                 call symU_linearsolve(Triplet%npair,TMP,CCpair%tau)
              endif
              call Triplet%expand(CCpair%tau)

              call TQT(npair,TMP,matS_T,matP1_T,matP2_T)
              call matvec(npair,npair,TMP,CCpair%tau,CCpair%Qtau)
              call Triplet%shrink(CCpair%Qtau)
              call project_T%solve(Triplet%npair,CCpair%Qtau)
              call Triplet%expand(CCpair%Qtau)

              associate(Info => CCpair%Info(iter,ieta))

                Info%energy = CCpair%mult*( &
                     dot_vec(npair,CCpair%tau,vecG(:,ipair,jpair)) - &
                     dot_vec(npair,CCpair%tau,vecG(:,jpair,ipair)) )
                call SOcheck(Info%orthog(1),Info%nonzero(1),&
                     npair,CCpair%tau,matS_T,matP1_T,Control%eps,inv_log10)
                call SOcheck(Info%orthog(2),Info%nonzero(2),&
                     npair,CCpair%Qtau,matS_T,matP1_T,Control%eps,inv_log10)

              end associate

            end associate

         enddo
      enddo

!-------------------------------------------------------------------------------

      if(iter>0) then
         fact = -1
         call add_mat(fact,npair,norb**2,DATA_R,DATA_E)
         error = sqrt(dot_mat(npair,norb**2,DATA_E,DATA_E))
      else
         error = 0
      endif
      
      associate(full => fullInfo(iter,ieta))

        full%energy = 0
        do i=1,2
           full%orthog(i) = -huge(0)
        enddo

        do jpair=1,norb
           do ipair=1,norb
              associate(Info => CCpairs(ipair,jpair)%Info(iter,ieta))
                full%energy = full%energy + Info%energy
                do i=1,2
                   if(Info%nonzero(i)) &
                        full%orthog(i) = max(full%orthog(i),Info%orthog(i))
                enddo
              end associate
           enddo
        enddo

        energy_this = full%energy
        if(iter>0) then
           edelta_this = energy_this - energy_prev
        else
           edelta_this = energy_this
        endif
        energy_prev = energy_this
        edelta_prev = edelta_this

        if(iter==0) write(*,'(2x,a,34x,a,25x,a,6x,a,5x,a,1x,a)') &
             'iter','energy','dE','error','P(X)','P(QX)'
        call write_fform(energy_this,sval)
        write(*,'(i5,a,2x,2es9.1,2x,2f6.1)') iter,&
             trim(sval),&
             dble(edelta_this),&
             dble(error),&
             (dble(full%orthog(i)),i=1,2)

      end associate

      if(stage>0) then

         if(iter>0) then
            if(abs(error)<Control%CC_THR) then
               write(*,'(a)') 'CONVERGENCE'
               exit
            endif
            if(abs(edelta_this/edelta_prev)>10) then
               write(*,'(a)') 'DIVERGENCE'
               exit
            endif
         endif

      else

         exit

      endif

      if(iter>DIIS_START.and.do_DIIS) then

         call DIIS%step(DATA_tau,DATA_E)

         do jpair=1,norb
            do ipair=1,jpair

               associate(CCpair => CCpairs(jpair,ipair))

                 call TQT(npair,TMP,matS_S,matP1_S,matP2_S)
                 call matvec(npair,npair,TMP,CCpair%tau,CCpair%Qtau)
                 call project_S%solve(npair,CCpair%Qtau)

                 call MWO(npair,TMP,norb,orbE,jpair,ipair,eta(ieta),&
                      matF_S,matS_S,matP1_S,matPe_S,matP2_S)
                 call matvec(npair,npair,TMP,CCpair%tau,CCpair%E)

               end associate

               if(ipair==jpair) cycle

               associate(CCpair => CCpairs(ipair,jpair))

                 call TQT(npair,TMP,matS_T,matP1_T,matP2_T)
                 call matvec(npair,npair,TMP,CCpair%tau,CCpair%Qtau)
                 call Triplet%shrink(CCpair%Qtau)
                 call project_T%solve(Triplet%npair,CCpair%Qtau)
                 call Triplet%expand(CCpair%Qtau)

                 call MWO(npair,TMP,norb,orbE,ipair,jpair,eta(ieta),&
                      matF_T,matS_T,matP1_T,matPe_T,matP2_T)
                 call matvec(npair,npair,TMP,CCpair%tau,CCpair%E)

               end associate

            enddo
         enddo

      else

         call copy_mat(npair,norb**2,DATA_R,DATA_E)

      endif

!-------------------------------------------------------------------------------

   enddo
   final_iter(ieta) = iter

   if(final_iter(ieta)>Control%CC_MAXITER) then
      final_iter(ieta) = Control%CC_MAXITER
      write(*,'(a)') 'NO CONVERGENCE'
   endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

   if(stage>1) then
      do jpair=1,norb
         do ipair=1,norb
            call CCpairs(ipair,jpair)%LHS%free
         enddo
      enddo
   endif

   if(Control%IPRINT>=1) then
      write(*,'()')
      write(*,'(a,t48,a,19x,2(2x,a,1x,a))') 'MP2',&
           'energy','P(S)','P(QS)','P(T)','P(QT)'
      do jpair=1,norb
         do ipair=1,jpair
            associate(Info => &
                 CCpairs(jpair,ipair)%Info(0,ieta))
              call write_fform(Info%energy,sval)
              write(*,'(2i5,2x,a,2f6.1)',advance='no') ipair,jpair,sval,&
                   (dble(Info%orthog(i)),i=1,2)
            end associate
            if(ipair/=jpair) then
               associate(Info => &
                    CCpairs(ipair,jpair)%Info(0,ieta))
                 write(*,'(2f6.1)') &
                      (dble(Info%orthog(i)),i=1,2)
               end associate
            else
               write(*,*)
            endif
         enddo
      enddo
      if(stage>0) then
         write(*,'()')
         write(*,'(a,t48,a,19x,2(2x,a,1x,a))') trim(Control%calc_type),&
              'energy','P(S)','P(QS)','P(T)','P(QT)'
         do jpair=1,norb
            do ipair=1,jpair
               associate(Info => &
                    CCpairs(jpair,ipair)%Info(final_iter(ieta),ieta))
                 call write_fform(Info%energy,sval)
                 write(*,'(2i5,2x,a,2f6.1)',advance='no') ipair,jpair,sval,&
                      (dble(Info%orthog(i)),i=1,2)
               end associate
               if(ipair/=jpair) then
                  associate(Info => &
                       CCpairs(ipair,jpair)%Info(final_iter(ieta),ieta))
                    write(*,'(2f6.1)') &
                         (dble(Info%orthog(i)),i=1,2)
                  end associate
               else
                  write(*,*)
               endif
            enddo
         enddo
      endif
   endif

enddo

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

write(*,'()')
write(*,'(a)') repeat('-',80)
write(*,'()')

write(*,'(5x,a)') 'FINAL RESULTS : MP2'
write(*,'()')
write(*,'(5x,a,39x,a,19x,2x,a,1x,a)') 'eta','energy','P(X)','P(QX)'
do ieta=1,neta
   associate(full => fullInfo(0,ieta))
     write(*,'(es10.1,2x)',advance='no') dble(eta(ieta))
     call write_fform(full%energy,sval)
     write(*,'(a,2f6.1)') sval,(dble(full%orthog(i)),i=1,2)
   end associate
enddo
if(stage>0) then
   write(*,'()')
   write(*,'(5x,2a)') 'FINAL RESULTS : ',trim(Control%calc_type)
   write(*,'()')
   write(*,'(5x,a,39x,a,19x,2x,a,1x,a)') 'eta','energy','P(X)','P(QX)'
   do ieta=1,neta
      associate(full => fullInfo(final_iter(ieta),ieta))
        write(*,'(es10.1,2x)',advance='no') dble(eta(ieta))
        call write_fform(full%energy,sval)
        write(*,'(a,2f6.1)') sval,(dble(full%orthog(i)),i=1,2)
      end associate
   enddo
endif

write(*,'()')
write(*,'(5x,2a)') 'FINAL RESULTS : TOTAL MP2'
write(*,'()')
write(*,'(5x,a,39x,a)') 'eta','energy'
do ieta=1,neta
   associate(full => fullInfo(0,ieta))
     write(*,'(es10.1,2x)',advance='no') dble(eta(ieta))
     call write_fform(energSCF+full%energy,sval)
     write(*,'(a)') sval
   end associate
enddo
if(stage>0) then
   write(*,'()')
   write(*,'(5x,2a)') 'FINAL RESULTS : TOTAL ',trim(Control%calc_type)
   write(*,'()')
   write(*,'(5x,a,39x,a)') 'eta','energy'
   do ieta=1,neta
      associate(full => fullInfo(final_iter(ieta),ieta))
        write(*,'(es10.1,2x)',advance='no') dble(eta(ieta))
        call write_fform(energSCF+full%energy,sval)
        write(*,'(a)') sval
      end associate
   enddo
endif

if(stage>2) deallocate(pair_energy)

call project_T%free
call project_S%free

call Triplet%free

deallocate(final_iter)
deallocate(fullInfo)

do jpair=1,norb
   do ipair=1,norb
      associate(CCpair => CCpairs(ipair,jpair))

        if(stage>1) deallocate(CCpair%LHS)

        deallocate(CCpair%Info)

        nullify(CCpair%E)
        nullify(CCpair%R)
        nullify(CCpair%Qtau)
        nullify(CCpair%tau)

      end associate
   enddo
enddo
deallocate(CCpairs)

if(do_DIIS) call DIIS%free
deallocate(DATA_E)
deallocate(DATA_R)
deallocate(DATA_Qtau)
deallocate(DATA_tau)

!-------------------------------------------------------------------------------

deallocate(vecG)
deallocate(matP2_S,matP2_T)
deallocate(matF_S,matF_T)
deallocate(matS_S,matS_T)
deallocate(matM_S,matM_T)
deallocate(matJ_S,matJ_T)
deallocate(matPe_S,matPe_T)
deallocate(matP1_S,matP1_T)
if(stage>0) deallocate(matL_S,matL_T)

!-------------------------------------------------------------------------------

deallocate(TMP)

deallocate(eta)

write(*,'()')
call timer('solve',Tcpu,Twall)
call timer('CC',Tcpu_TOT,Twall_TOT)

end subroutine calc_CC

subroutine MWO(n,TOT,norb,orbE,iorb,jorb,eta,matF,matS,matP1,matPe,matP2)
implicit none
integer,intent(in) :: n
type(mp_real),intent(out) :: TOT(:,:)
integer,intent(in) :: norb
type(mp_real),intent(in) :: orbE(:)
integer,intent(in) :: iorb,jorb
type(mp_real),intent(in) :: eta
type(mp_real),intent(in) :: matF(:,:),matS(:,:),matP1(:,:),matPe(:,:),matP2(:,:)
type(mp_real) :: e_ij,e_HO,fact

e_ij = orbE(iorb) + orbE(jorb)
e_HO = orbE(norb)

call copy_mat(n,n,matF,TOT)
fact = -e_ij
call add_mat(fact,n,n,matS,TOT)

fact = e_ij - e_HO + eta
call add_mat(fact,n,n,matP1,TOT)
fact = -1
call add_mat(fact,n,n,matPe,TOT)
fact = 2*e_HO - e_ij
call add_mat(fact,n,n,matP2,TOT)

end subroutine MWO

subroutine TQT(n,TOT,matS,matP1,matP2)
implicit none
integer,intent(in) :: n
type(mp_real),intent(OUT) :: TOT(:,:)
type(mp_real),intent(in) :: matS(:,:),matP1(:,:),matP2(:,:)
type(mp_real) :: fact

call copy_mat(n,n,matS,TOT)
fact = -1
call add_mat(fact,n,n,matP1,TOT)
fact = 1
call add_mat(fact,n,n,matP2,TOT)

end subroutine TQT

subroutine SOcheck(orthog,nonzero,n,vec,matS,matP1,eps,inv_log10)
implicit none
type(mp_real),intent(out) :: orthog
logical,intent(out) :: nonzero
integer,intent(in) :: n
type(mp_real),intent(in) :: vec(:)
type(mp_real),intent(in) :: matS(:,:),matP1(:,:)
type(mp_real),intent(in) :: eps,inv_log10
type(mp_real) :: val

val = vecmatvec(n,n,matS,vec,vec)

nonzero = abs(val)>eps

if(nonzero) then
   val = vecmatvec(n,n,matP1,vec,vec)/val
   orthog = log(abs(val))*inv_log10
else
   orthog = 0
endif

end subroutine SOcheck

end module driver
