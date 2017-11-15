! subroutine integrals_matP(matP,matC,parC,IPRINT)
! implicit none
! type(mp_real),intent(out) :: matP(:,:)
! type(mp_real),intent(in) :: matC(:,:)
! integer,intent(in) :: parC(:)
! integer,intent(in) :: IPRINT
! real(8) :: TCPU,Twall

! integer :: n1,n2,nC,iocc,jocc,i,j,ij,v1,v2
! type(mp_real) :: a1,a2,a12,int0,c0,mult
! type(mp_real) :: val
! type(prod2HData) :: prod2H
! type(table2Data) :: aux,auxOP
! type(table4Data),allocatable,target :: aux2C(:)
! integer,allocatable :: pair_list(:,:)
! type(mp_real),allocatable :: ints(:,:)
! integer,allocatable :: par_ints(:)
! type(mp_real),pointer :: cont_0(:,:),cont_1(:,:)
! integer,pointer :: se_0(:,:),se_1(:,:)
! integer,allocatable :: par_cont(:)

! if(IPRINT>=1) call timer('START',Tcpu,Twall)

! n1 = G_nbas-1
! n2 = G_nprim
! call prod2H_create(prod2H,n1,n2)
! call prod2H_norm(prod2H)

! n1 = G_nprim + (G_nbas-1)
! n2 = G_nprim + (G_nbas-1)
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

! n1 = G_nprim + (G_nbas-1)
! n2 = G_nprim
! nC = G_nbas-1
! allocate(aux2C(G_nocc))
! do iocc=1,G_nocc
!    call aux2C_create(aux2C(iocc),n1,n2,nC,matC(:,iocc),prod2H)
! enddo


! if(IPRINT>=1) call timer('vec0',Tcpu,Twall)
! end subroutine integrals_matP
