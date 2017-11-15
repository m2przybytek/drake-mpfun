program testy
use mpmodule
use general
use inputread
use eproblem_real
use eproblem_mp
use lproblem_mp
use diis
use integrals
implicit none
integer :: IPRINT
real(8) :: mp_digits
type(mp_real) :: mp_eps
integer :: DIIS_size,MAXITER
logical :: mp_SCF
type(mp_real) :: thr_SCF
integer :: n_0,n_1,iter,i,j,k,l,m
character(50) :: sval
type(mp_real) :: energSCF,energSCF_prev,energSCF_diff,error
type(mp_real) :: energMP2
integer,allocatable :: parC(:)
type(mp_real),allocatable :: orbE(:),matC(:,:)
type(mp_real),allocatable :: PPOO(:,:,:,:)
type(mp_real),allocatable :: scfS(:,:),scfH(:,:)
type(mp_real),allocatable :: scfF(:,:),scfD(:,:),scfP(:,:)
type(scfJauxData) :: scfJaux
type(DIISData) :: DIIS
type(mp_real),allocatable :: matS_S(:,:),matS_T(:,:)
type(mp_real),allocatable :: matF_S(:,:),matF_T(:,:)
type(mp_real),allocatable :: matJ_S(:,:,:,:),matJ_T(:,:,:,:)
type(mp_real),allocatable :: vec0(:,:,:),vec1(:,:,:),vecP(:,:,:)
type(CCpairData),allocatable :: CCpairs(:,:)
type(mp_real) :: val
type(mp_real),allocatable :: work(:,:)
type(mp_real),allocatable :: eval_0(:),evec_0(:,:),tmpF_0(:,:),tmpS_0(:,:)
type(mp_real),allocatable :: eval_1(:),evec_1(:,:),tmpF_1(:,:),tmpS_1(:,:)
real(eprec),allocatable :: Eeval_0(:),Eevec_0(:,:),EtmpF_0(:,:),EtmpS_0(:,:)
real(eprec),allocatable :: Eeval_1(:),Eevec_1(:,:),EtmpF_1(:,:),EtmpS_1(:,:)
type(mp_real),allocatable :: LHS(:,:),RHS(:)
type(mp_real),allocatable :: eval(:),evec(:,:),tmpF(:,:),tmpS(:,:)
integer,allocatable :: list_shrink(:)

IPRINT = 0

! G_nocc  = 1
! G_nbas  = 100
! G_nprim = 50
! call G_set_npair
! G_gfac  = mpreal(2.d0)
! G_alpha = mpreal(1.d0)

call read_Input()
call G_set_npair

mp_digits = mpwds*mpdpw
mp_eps = mpreal(10.d0)**mpreald(-mp_digits)
write(*,*) 'Approximate number of digits: ',mp_digits

DIIS_size = 6
MAXITER   = 100
mp_SCF    = .true.
!mp_SCF    = .false.
if(mp_SCF) then
   thr_SCF = mpreal('1.e-32')
else
   call G_real_to_mp(1.d3*epsilon(0._eprec),thr_SCF)
endif

allocate(orbE(G_nbas),matC(G_nbas,G_nbas),parC(G_nbas))
allocate(PPOO(G_nbas,G_nbas,G_nocc,G_nocc))

!-------------------------------------------------------------------------------

n_0 = (G_nbas-1)/2+1
n_1 = G_nbas/2

allocate(scfS(G_nbas,G_nbas),scfH(G_nbas,G_nbas))
allocate(scfF(G_nbas,G_nbas),scfD(G_nbas,G_nbas),scfP(G_nbas,G_nbas))

call integrals_scfSH(scfS,scfH,IPRINT)
call prepare_scfJaux(scfJaux,IPRINT)

allocate(work(G_nbas,G_nbas))

call init_DIIS(DIIS,G_nbas,DIIS_size)
if(mp_SCF) then
   allocate(eval_0(n_0),evec_0(n_0,n_0),tmpF_0(n_0,n_0),tmpS_0(n_0,n_0))
   allocate(eval_1(n_1),evec_1(n_1,n_1),tmpF_1(n_1,n_1),tmpS_1(n_1,n_1))
else
   allocate(Eeval_0(n_0),Eevec_0(n_0,n_0),EtmpF_0(n_0,n_0),EtmpS_0(n_0,n_0))
   allocate(Eeval_1(n_1),Eevec_1(n_1,n_1),EtmpF_1(n_1,n_1),EtmpS_1(n_1,n_1))
endif

energSCF = 0
scfF(:,:) = scfH
write(*,*)

do iter=1,MAXITER
   
   orbE = mpreal(0.d0)
   matC = mpreal(0.d0)
   if(mp_SCF) then
      tmpF_0(:,:) = scfF(1:G_nbas:2,1:G_nbas:2)
      tmpS_0(:,:) = scfS(1:G_nbas:2,1:G_nbas:2)
      call symU_diagonalize_mp(mp_eps,n_0,eval_0,evec_0,tmpF_0,tmpS_0)
      tmpF_1(:,:) = scfF(2:G_nbas:2,2:G_nbas:2)
      tmpS_1(:,:) = scfS(2:G_nbas:2,2:G_nbas:2)
      call symU_diagonalize_mp(mp_eps,n_1,eval_1,evec_1,tmpF_1,tmpS_1)
      i = 0
      j = 0
      do k=1,G_nbas
         if(i<n_0) then
            if(j<n_1) then
               if(eval_0(i+1)<=eval_1(j+1)) then
                  l = 0
               else
                  l = 1
               endif
            else
               l = 0
            endif
         else
            l = 1
         endif
         select case(l)
         case(0)
            i = i + 1
            orbE(k) = eval_0(i)
            matC(1:G_nbas:2,k) = evec_0(:,i)
            parC(k) = 0
         case(1)
            j = j + 1
            orbE(k) = eval_1(j)
            matC(2:G_nbas:2,k) = evec_1(:,j)
            parC(k) = 1
         end select
      enddo
   else
      call G_mp_to_real(scfF(1:G_nbas:2,1:G_nbas:2),EtmpF_0)
      call G_mp_to_real(scfS(1:G_nbas:2,1:G_nbas:2),EtmpS_0)
      call symU_diagonalize_real(n_0,Eeval_0,Eevec_0,EtmpF_0,EtmpS_0)
      call G_mp_to_real(scfF(2:G_nbas:2,2:G_nbas:2),EtmpF_1)
      call G_mp_to_real(scfS(2:G_nbas:2,2:G_nbas:2),EtmpS_1)
      call symU_diagonalize_real(n_1,Eeval_1,Eevec_1,EtmpF_1,EtmpS_1)
      i = 0
      j = 0
      do k=1,G_nbas
         if(i<n_0) then
            if(j<n_1) then
               if(Eeval_0(i+1)<=Eeval_1(j+1)) then
                  l = 0
               else
                  l = 1
               endif
            else
               l = 0
            endif
         else
            l = 1
         endif
         select case(l)
         case(0)
            i = i + 1
            call G_real_to_mp(Eeval_0(i),orbE(k))
            call G_real_to_mp(Eevec_0(:,i),matC(1:G_nbas:2,k))
            parC(k) = 0
         case(1)
            j = j + 1
            call G_real_to_mp(Eeval_1(j),orbE(k))
            call G_real_to_mp(Eevec_1(:,j),matC(2:G_nbas:2,k))
            parC(k) = 1
         end select
      enddo
   endif

   do i=1,G_nocc
      call outer_vecproduct(G_nbas,matC(:,i),matC(:,i),scfD,(i==1))
   enddo
   call integrals_scfJ(work,scfD,scfJaux,IPRINT)

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

if(mp_SCF) then
   deallocate(eval_1,evec_1,tmpF_1,tmpS_1)
   deallocate(eval_0,evec_0,tmpF_0,tmpS_0)
else
   deallocate(Eeval_1,Eevec_1,EtmpF_1,EtmpS_1)
   deallocate(Eeval_0,Eevec_0,EtmpF_0,EtmpS_0)
endif
call free_DIIS(DIIS)

do j=1,G_nocc
   do i=1,G_nocc
      call outer_vecproduct(G_nbas,matC(:,i),matC(:,j),scfD)
      call integrals_scfJ(PPOO(:,:,i,j),scfD,scfJaux,IPRINT)
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
write(*,'(a)') 'Orbital energies'
do i=1,G_nocc
   call my_mpfform(orbE(i),50,40,sval)
   write(*,'(i5,a)') i,trim(sval)
enddo

write(*,*)
write(*,'(a)') 'Pair energies (standard)'
energMP2 = 0
do j=1,G_nocc
   do i=1,G_nocc
      val = 0
      if(i>=j) then
         do k=G_nocc+1,G_nbas
            do l=G_nocc+1,G_nbas
               val = val + &
                    (PPOO(k,l,i,j) + PPOO(l,k,i,j))**2/ &
                    (orbE(k)+orbE(l)-orbE(i)-orbE(j))
            enddo
         enddo
         if(i==j) then
            val = -0.25d0*val
         else
            val = -0.5d0*val
         endif
      else
         do k=G_nocc+1,G_nbas
            do l=G_nocc+1,G_nbas
               val = val + &
                    (PPOO(k,l,i,j) - PPOO(l,k,i,j))**2/ &
                    (orbE(k)+orbE(l)-orbE(i)-orbE(j))
            enddo
         enddo
         val = -1.5d0*val
      endif
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

write(*,*)
write(*,'(a)') 'stage 0'

allocate(matS_S(G_npair,G_npair),matS_T(G_npair,G_npair))
allocate(matF_S(G_npair,G_npair),matF_T(G_npair,G_npair))
allocate(&
     matJ_S(G_npair,G_npair,G_nocc,G_nocc),&
     matJ_T(G_npair,G_npair,G_nocc,G_nocc))
allocate(&
     vec0(G_npair,G_nocc,G_nocc),&
     vec1(G_npair,G_nocc,G_nocc),&
     vecP(G_npair,G_nocc,G_nocc))

call integrals_SH(matS_S,matS_T,matF_S,matF_T,IPRINT)
call integrals_J(matJ_S,matJ_T,matC,parC,IPRINT)

call integrals_vec0(vec0,matC,parC,IPRINT)
vec1(:,:,:) = mpreal(0.d0)
call integrals_vecP(vecP,matC,parC,IPRINT)

do k=1,G_nocc
   do j=1,G_npair
      do i=1,G_npair
         matF_S(i,j) = matF_S(i,j) + matJ_S(i,j,k,k)
      enddo
   enddo
enddo
do k=1,G_nocc
   do j=1,G_npair
      do i=1,G_npair
         matF_T(i,j) = matF_T(i,j) + matJ_T(i,j,k,k)
      enddo
   enddo
enddo

allocate(CCpairs(G_nocc,G_nocc))
do j=1,G_nocc
   do i=1,G_nocc
      associate(CCpair => CCpairs(i,j))

        call init_CCpair(CCpair)

        if(i>=j) then

           if(i==j) then
              CCpair%mult = mpreal(-0.5d0)
           else
              CCpair%mult = mpreal(-1.d0)
           endif

           do m=1,G_npair
              CCpair%vec(m) = &
                   - (vec1(m,i,j) + vec1(m,j,i))
           enddo
           do m=1,G_npair
              CCpair%vec(m) = CCpair%vec(m) &
                   + (vecP(m,i,j) + vecP(m,j,i))
           enddo
           do l=1,G_nocc
              do k=1,G_nocc
                 val = PPOO(k,l,i,j) + PPOO(k,l,j,i)
                 do m=1,G_npair
                    CCpair%vec(m) = CCpair%vec(m) &
                         - vec0(m,k,l)*val
                 enddo
              enddo
           enddo

        else

           CCpair%mult = mpreal(-3.d0)

           do m=1,G_npair
              CCpair%vec(m) = &
                   - (vec1(m,i,j) - vec1(m,j,i))
           enddo
           do m=1,G_npair
              CCpair%vec(m) = CCpair%vec(m) &
                   + (vecP(m,i,j) - vecP(m,j,i))
           enddo
           do l=1,G_nocc
              do k=1,G_nocc
                 val = PPOO(k,l,i,j) - PPOO(k,l,j,i)
                 do m=1,G_npair
                    CCpair%vec(m) = CCpair%vec(m) &
                         - vec0(m,k,l)*val
                 enddo
              enddo
           enddo

        endif

      end associate
   enddo
enddo

allocate(list_shrink(G_npair_shrink))
i = -1
j = -1
k = 0
l = 0
do while(G_next_pair(i,j))
   k = k + 1
   if(i/=j) then
      l = l + 1
      list_shrink(l) = k
   endif
enddo

allocate(LHS(G_npair,G_npair),RHS(G_npair))
do j=1,G_nocc
   do i=1,G_nocc
      associate(CCpair => CCpairs(i,j))

        if(i>=j) then

           val = orbE(i) + orbE(j)
           do l=1,G_npair
              do k=1,G_npair
                 LHS(k,l) = matF_S(k,l) - val*matS_S(k,l)
              enddo
           enddo
           do m=1,G_npair
              RHS(m) = CCpair%vec(m)
           enddo
           call symU_linearsolve_mp(G_npair,LHS,RHS)
           do m=1,G_npair
              CCpair%tau(m) = RHS(m)
           enddo

        else

           CCpair%tau = mpreal(0.d0)

        endif

        val = 0
        do m=1,G_npair
           val = val + CCpair%vec(m)*CCpair%tau(m)
        enddo
        val = CCpair%mult*val
        call my_mpfform(val,50,40,sval)
        write(*,'(2i5,a)') i,j,trim(sval)

      end associate
   enddo
enddo
deallocate(LHS,RHS)

allocate(eval(G_npair),evec(G_npair,G_npair))
allocate(tmpF(G_npair,G_npair),tmpS(G_npair,G_npair))
write(*,*)
write(*,'(a)') 'Problem S: energia H_S+J_S, stosunek wartosci wlasnych S_S'
tmpF(:,:) = matF_S
tmpS(:,:) = matS_S
call symU_diagonalize_mp(mp_eps,G_npair,eval,evec,tmpF,tmpS)
call mpwrite(6,60,40,eval(1))
tmpS(:,:) = matS_S
call symU_diagonalize_mp(mp_eps,G_npair,eval,evec,tmpS)
call mpwrite(6,60,40,eval(G_npair)/eval(1))
write(*,*)
write(*,'(a)') 'Problem S: energia H_T+J_T, stosunek wartosci wlasnych S_T'
tmpF(1:G_npair_shrink,1:G_npair_shrink) = matF_T(list_shrink,list_shrink)
tmpS(1:G_npair_shrink,1:G_npair_shrink) = matS_T(list_shrink,list_shrink)
call symU_diagonalize_mp(mp_eps,G_npair_shrink,eval,evec,tmpF,tmpS)
call mpwrite(6,60,40,eval(1))
tmpS(1:G_npair_shrink,1:G_npair_shrink) = matS_T(list_shrink,list_shrink)
call symU_diagonalize_mp(mp_eps,G_npair_shrink,eval,evec,tmpS)
call mpwrite(6,60,40,eval(G_npair_shrink)/eval(1))
deallocate(tmpF,tmpS)
deallocate(eval,evec)

deallocate(list_shrink)

do j=1,G_nocc
   do i=1,G_nocc
      associate(CCpair => CCpairs(i,j))
        call free_CCpair(CCpair)
      end associate
   enddo
enddo
deallocate(CCpairs)

deallocate(vec0,vec1,vecP)
deallocate(matJ_S,matJ_T)
deallocate(matF_S,matF_T)
deallocate(matS_S,matS_T)

!-------------------------------------------------------------------------------

deallocate(PPOO)
deallocate(orbE,matC,parC)

end program testy
