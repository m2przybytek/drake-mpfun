program testy
use mpmodule
use general
use eproblem_real
use eproblem_mp
use diis
use integrals
implicit none
real(8) :: mp_digits
type(mp_real) :: mp_eps
integer :: DIIS_size,MAXITER
logical :: mp_SCF
type(mp_real) :: thr_SCF
integer :: iter,i,j,k,l,s
character(50) :: sval
type(mp_real) :: energSCF,energSCF_prev,energSCF_diff,error
type(mp_real) :: energMP2
type(mp_real),allocatable :: orbE(:),matC(:,:),PPOO(:,:,:,:)
type(mp_real),allocatable :: scfS(:,:),scfH(:,:)
type(mp_real),allocatable :: scfF(:,:),scfD(:,:),scfP(:,:)
type(scfJauxData) :: scfJaux
type(DIISData) :: DIIS
type(mp_real),allocatable :: matS_S(:,:),matS_T(:,:)
type(mp_real),allocatable :: matF_S(:,:),matF_T(:,:)
type(mp_real),allocatable :: matJ_S(:,:,:,:),matJ_T(:,:,:,:)
type(mp_real),allocatable :: tmpF(:,:),tmpS(:,:)
real(eprec),allocatable :: Eeval(:),Eevec(:,:),EtmpF(:,:),EtmpS(:,:)
type(mp_real),allocatable :: work(:,:)
type(mp_real) :: val

G_nocc  = 2
G_nbas  = 6
G_nprim = 10
G_npair = G_npair_total()
G_gfac  = mpreal(2.d0)
G_alpha = mpreal(1.d0)

mp_digits = (mpwds-1)*mpdpw
mp_eps = mpreal(10.d0)**mpreald(-mp_digits)
write(*,*) 'Approximate number of digits: ',mp_digits

DIIS_size = 6
MAXITER   = 100
!mp_SCF    = .true.
mp_SCF    = .false.
if(mp_SCF) then
   thr_SCF = mpreald(1.d-30)
else
   thr_SCF = mpreal(1.d3*epsilon(0._eprec))
endif

allocate(orbE(G_nbas),matC(G_nbas,G_nbas),PPOO(G_nbas,G_nbas,G_nocc,G_nocc))

!-------------------------------------------------------------------------------

allocate(scfS(G_nbas,G_nbas),scfH(G_nbas,G_nbas))
allocate(scfF(G_nbas,G_nbas),scfD(G_nbas,G_nbas),scfP(G_nbas,G_nbas))

call integrals_scfSH(scfS,scfH)
call prepare_scfJaux(scfJaux)

allocate(work(G_nbas,G_nbas))

if(mp_SCF) then
   allocate(tmpF(G_nbas,G_nbas),tmpS(G_nbas,G_nbas))
else
   allocate(&
        Eeval(G_nbas),Eevec(G_nbas,G_nbas),&
        EtmpF(G_nbas,G_nbas),EtmpS(G_nbas,G_nbas))
endif
call init_DIIS(DIIS,G_nbas,DIIS_size)

energSCF = 0
scfF(:,:) = scfH
write(*,*)

do iter=1,MAXITER

   if(mp_SCF) then
      tmpF(:,:) = scfF
      tmpS(:,:) = scfS
      call symU_diagonalize_mp(mp_eps,G_nbas,orbE,matC,tmpF,tmpS)
   else
      call G_mp_to_real(scfF,EtmpF)
      call G_mp_to_real(scfS,EtmpS)
      call symU_diagonalize_real(G_nbas,Eeval,Eevec,EtmpF,EtmpS)
      call G_real_to_mp(Eeval,orbE)
      call G_real_to_mp(Eevec,matC)
   endif

   do i=1,G_nocc
      call outer_vecproduct(G_nbas,matC(:,i),matC(:,i),scfD,(i==1))
   enddo
   call integrals_scfJ(work,scfD,scfJaux)

   energSCF_prev = energSCF
   energSCF = 0
   do j=1,G_nbas
      do i=1,G_nbas
         scfF(i,j) = scfH(i,j) + work(i,j)
         energSCF = energSCF + (scfH(i,j) + scfF(i,j))*scfD(i,j)
      enddo
   enddo
   energSCF_diff = energSCF - energSCF_prev
   call scf_FDS_SDF(G_nbas,error,scfP,scfF,scfD,scfS,work)

   call my_mpfform(energSCF,50,40,sval)
   write(*,'(i5,a,2x,2es12.3)') iter,&
        trim(sval),&
        dble(energSCF_diff),&
        dble(error)

   if(abs(error)<thr_SCF) exit

   call use_DIIS(DIIS,scfF,scfP)

enddo

call free_DIIS(DIIS)
if(mp_SCF) then
   deallocate(tmpF,tmpS)
else
   deallocate(Eeval,Eevec,EtmpF,EtmpS)
endif

do j=1,G_nocc
   do i=1,G_nocc
      call outer_vecproduct(G_nbas,matC(:,i),matC(:,j),scfD)
      call integrals_scfJ(PPOO(:,:,i,j),scfD,scfJaux)
      call fullMO_transform(G_nbas,matC,PPOO(:,:,i,j),work)
   enddo
enddo

energSCF_prev = 0
do i=1,G_nocc
   energSCF_prev = energSCF_prev + 2*orbE(i)
enddo
do j=1,G_nocc
   do i=1,G_nocc
      energSCF_prev = energSCF_prev - PPOO(i,i,j,j)
   enddo
enddo
write(*,*)
call my_mpfform(energSCF,50,40,sval)
write(*,'(2a)') 'SCF energy: ',trim(sval)
call my_mpfform(energSCF_prev,50,40,sval)
write(*,'(2a)') '     check: ',trim(sval)

write(*,*)
energMP2 = 0
do j=1,G_nocc
   do i=1,G_nocc
      if(i>=j) then
         s = 1
      else
         s = 3
      endif
      val = 0
      do k=G_nocc+1,G_nbas
         do l=G_nocc+1,G_nbas
            val = val + &
                 (PPOO(k,l,i,j) + (2-s)*PPOO(l,k,i,j))**2/ &
                 (orbE(k)+orbE(l)-orbE(i)-orbE(j))
         enddo
      enddo
      val = (-0.5d0*s)*val
      if(i==j) val = 0.5d0*val
      energMP2 = energMP2 + val
      call my_mpfform(val,50,40,sval)
      write(*,'(2i5,a)') i,j,trim(sval)
   enddo
enddo
write(*,*)
call my_mpfform(energMP2,50,40,sval)
write(*,'(2a)') 'MP2 energy: ',trim(sval)

deallocate(work)

call free_scfJaux(scfJaux)

deallocate(scfF,scfD,scfP)
deallocate(scfS,scfH)

!-------------------------------------------------------------------------------

allocate(matS_S(G_npair,G_npair),matS_T(G_npair,G_npair))
allocate(matF_S(G_npair,G_npair),matF_T(G_npair,G_npair))
allocate(&
     matJ_S(G_npair,G_npair,G_nocc,G_nocc),&
     matJ_T(G_npair,G_npair,G_nocc,G_nocc))

call integrals_SH(matS_S,matS_T,matF_S,matF_T)
call integrals_J(matJ_S,matJ_T,matC)

block
  type(mp_real),allocatable :: eval(:),evec(:,:)
  allocate(eval(G_npair),evec(G_npair,G_npair))
  allocate(tmpF(G_npair,G_npair),tmpS(G_npair,G_npair))
  tmpF(:,:) = matF_S
  tmpS(:,:) = matS_S
  call symU_diagonalize_mp(mp_eps,G_npair,eval,evec,tmpF,tmpS)
  call mpwrite(6,60,40,eval(1))
  tmpS(:,:) = matS_S
  call symU_diagonalize_mp(mp_eps,G_npair,eval,evec,tmpS)
  call mpwrite(6,60,40,eval(G_npair)/eval(1))
  deallocate(tmpF,tmpS)
  deallocate(eval,evec)
end block

deallocate(matJ_S,matJ_T)
deallocate(matF_S,matF_T)
deallocate(matS_S,matS_T)

!-------------------------------------------------------------------------------

deallocate(orbE,matC,PPOO)

end program testy
