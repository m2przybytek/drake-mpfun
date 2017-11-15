subroutine integrals_scfSH(scfS,scfH,IPRINT)
  implicit none
  type(mp_real),intent(out) :: scfS(0:,0:)
  type(mp_real),intent(out) :: scfH(0:,0:)
  integer,intent(in) :: IPRINT
  real(8) :: Tcpu,Twall
  type(mp_real),allocatable :: norm(:)
  type(mp_real) :: val
  integer :: i

  if(IPRINT>=1) call timer('START',Tcpu,Twall)

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

  if(IPRINT>=1) call timer('scfSH',Tcpu,Twall)

end subroutine integrals_scfSH

subroutine prepare_scfJaux(scfJaux,IPRINT)
implicit none
type(scfJauxData) :: scfJaux
integer,intent(in) :: IPRINT
real(8) :: Tcpu,Twall
integer :: n1,n2,i,j,ij,v2
real(8) :: safe_nbas
type(table2Data) :: aux22

if(IPRINT>=1) call timer('START',Tcpu,Twall)

safe_nbas = 0.5d0 - 0.5d0/G_nbas

n1 = G_nbas-1
n2 = G_nbas-1
call prod2H_create(scfJaux%prod2H,n1,n2)
call prod2H_norm(scfJaux%prod2H)

n1 = G_nbas
n2 = (G_nbas*(G_nbas+1))/2
allocate(scfJaux%ints(n1,n2))

n1 = G_nbas
n2 = 1
!$ n2 = omp_get_max_threads()
allocate(scfJaux%sum_0(n1,0:n2-1),scfJaux%sum_1(n1,0:n2-1))

n1 = 2*(G_nbas-1)
n2 = 2*(G_nbas-1)
call aux22_create(aux22,n1,n2)

scfJaux%n_0 = aux22%n1_0
scfJaux%n_1 = aux22%n1_1

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v2)
do ij=1,(G_nbas*(G_nbas+1))/2
   j = int(sqrt(0.25d0 + 2.d0*(ij-1)) - safe_nbas)
   i = (ij-1) - (j*(j+1))/2
   v2 = i + j
   select case(mod(v2,2))
   case(0)
      v2 = v2/2+1
      call product2_unit(aux22%n1_0,aux22%val_0,v2,scfJaux%prod2H%unit(i,j),&
           scfJaux%ints(:,ij))
   case(1)
      v2 = (v2+1)/2
      call product2_unit(aux22%n1_1,aux22%val_1,v2,scfJaux%prod2H%unit(i,j),&
           scfJaux%ints(:,ij))
   end select
enddo
!$OMP END PARALLEL DO

call table2_free(aux22)

if(IPRINT>=1) call timer('prepare_scfJ',Tcpu,Twall)

end subroutine prepare_scfJaux

subroutine integrals_scfJ(scfJ,scfD,scfJaux,IPRINT)
implicit none
type(mp_real),intent(out) :: scfJ(0:,0:)
type(mp_real),intent(in) :: scfD(0:,0:)
type(scfJauxData) :: scfJaux
integer,intent(in) :: IPRINT
real(8) :: Tcpu,Twall
real(8) :: safe_nbas
integer :: nthr,ithr,i,j,ij,v1
type(mp_real) :: val

if(IPRINT>=2) call timer('START',Tcpu,Twall)

safe_nbas = 0.5d0 - 0.5d0/G_nbas

!$OMP PARALLEL PRIVATE(ithr,i,j,v1,val)

!$OMP SINGLE
nthr = 1
!$ nthr = omp_get_num_threads()
!$OMP END SINGLE

!$OMP SECTIONS
!$OMP SECTION
scfJaux%sum_0(1:scfJaux%n_0,0:nthr-1) = mpreal(0.d0)
!$OMP SECTION
scfJaux%sum_1(1:scfJaux%n_1,0:nthr-1) = mpreal(0.d0)
!$OMP END SECTIONS

!$OMP DO
do ij=1,(G_nbas*(G_nbas+1))/2
   ithr = 0
!$ ithr = omp_get_thread_num()
   j = int(sqrt(0.25d0 + 2.d0*(ij-1)) - safe_nbas)
   i = (ij-1) - (j*(j+1))/2
   v1 = i + j
   select case(mod(v1,2))
   case(0)
      val = scfD(i,j)
      if(i/=j) val = val + scfD(j,i)
      call accumulate(val,scfJaux%n_0,scfJaux%ints(:,ij),scfJaux%sum_0(:,ithr))
   case(1)
      val = scfD(i,j) + scfD(j,i)
      call accumulate(val,scfJaux%n_1,scfJaux%ints(:,ij),scfJaux%sum_1(:,ithr))
   end select
enddo
!$OMP END DO

!$OMP SECTIONS
!$OMP SECTION
do ithr=1,nthr-1
   call addition(scfJaux%n_0,scfJaux%sum_0(:,ithr),scfJaux%sum_0(:,0))
enddo
!$OMP SECTION
do ithr=1,nthr-1
   call addition(scfJaux%n_1,scfJaux%sum_1(:,ithr),scfJaux%sum_1(:,0))
enddo
!$OMP END SECTIONS

!$OMP DO SCHEDULE(DYNAMIC)
do ij=1,(G_nbas*(G_nbas+1))/2
   j = int(sqrt(0.25d0 + 2.d0*(ij-1)) - safe_nbas)
   i = (ij-1) - (j*(j+1))/2
   v1 = i + j
   select case(mod(v1,2))
   case(0)
      v1 = v1/2+1
      val = G_gfac*product1_unit(v1,scfJaux%prod2H%unit(i,j),scfJaux%sum_0(:,0))
      scfJ(i,j) = val
      if(i/=j) scfJ(j,i) = val
   case(1)
      v1 = (v1+1)/2
      val = G_gfac*product1_unit(v1,scfJaux%prod2H%unit(i,j),scfJaux%sum_1(:,0))
      scfJ(i,j) = val
      scfJ(j,i) = val
   end select
enddo
!$OMP END DO

!$OMP END PARALLEL

if(IPRINT>=2) call timer('use_scfJ',Tcpu,Twall)

end subroutine integrals_scfJ

subroutine free_scfJaux(scfJaux)
  implicit none
  type(scfJauxData) :: scfJaux

  deallocate(scfJaux%sum_0,scfJaux%sum_1)
  deallocate(scfJaux%ints)
  call prod2H_free(scfJaux%prod2H)

end subroutine free_scfJaux
