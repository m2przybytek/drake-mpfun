subroutine integrals_J(matJ_S,matJ_T,matC,parC,IPRINT)
implicit none
type(mp_real),intent(out) :: matJ_S(:,:,:,:),matJ_T(:,:,:,:)
type(mp_real),intent(in) :: matC(:,:)
integer,intent(in) :: parC(:)
integer,intent(in) :: IPRINT
real(8) :: Tcpu,Twall
integer :: nthr,ithr,n1,n2,nC,iocc,jocc,parCC,i,j,k,l,ij,kl,v1,v2
real(8) :: safe_nbas,safe_nprim
type(mp_real) :: alpha_invsqrt,a1,a2,a12,c0,c2
type(mp_real) :: val,val_A1,val_A2,val_B1,val_B2,val_A,val_B
type(prod2HData) :: prod2H
type(table2Data) :: aux,auxOP
type(table4Data) :: aux4C
type(mp_real),allocatable :: ints(:,:)
integer,allocatable :: pair_list(:,:)
type(mp_real),allocatable :: CC(:,:),sum_0(:,:),sum_1(:,:)
type(contData),allocatable :: cont(:)

if(IPRINT>=1) call timer('START',Tcpu,Twall)

safe_nbas  = 0.5d0 - 0.5d0/G_nbas
safe_nprim = 0.5d0 - 0.5d0/(G_nprim+1)
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
   j = int(sqrt(0.25d0 + 2.d0*(ij-1)) - safe_nprim)
   i = (ij-1) - (j*(j+1))/2
   v2 = i + j
   select case(mod(v2,2))
   case(0)
      v2 = v2/2+1
      call product2_unit(auxOP%n1_0,auxOP%val_0,v2,prod2H%unit(i,j),&
           ints(:,ij))
   case(1)
      v2 = (v2+1)/2
      call product2_unit(auxOP%n1_1,auxOP%val_1,v2,prod2H%unit(i,j),&
           ints(:,ij))
   end select
enddo
!$OMP END PARALLEL DO

call table2_free(auxOP)
call table2_free(aux)

if(IPRINT>=1) call timer('prepare_J',Tcpu,Twall)

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

allocate(CC(0:G_nbas-1,0:G_nbas-1))

n1 = G_nbas
n2 = 1
!$ n2 = omp_get_max_threads()
allocate(sum_0(n1,0:n2-1),sum_1(n1,0:n2-1))

nC = ((G_nprim+1)*(G_nprim+2))/2
allocate(cont(nC))
ij = 0
do j=0,G_nprim
   do i=0,j
      ij = ij + 1
      associate(this => cont(ij))
        nC = i + j + 2*(G_nbas-1)
        this%n_0 = nC/2+1
        this%n_1 = (nC+1)/2
        allocate(this%val(max(this%n_0,this%n_1)))
      end associate
   enddo
enddo

do jocc=1,G_nocc
   do iocc=1,jocc

      if(IPRINT>=1) call timer('START',Tcpu,Twall)

      call outer_vecproduct(G_nbas,matC(:,iocc),matC(:,jocc),CC)
      parCC = mod(parC(iocc)+parC(jocc),2)

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

      n1 = 2*G_nprim + 2*(G_nbas-1)
      n2 = 2*G_nprim
      nC = 2*(G_nbas-1)
      call aux4C_create(aux4C,n1,n2,nC,sum_0(:,0),sum_1(:,0))

      select case(parCC)
      case(0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v2)
         do ij=1,((G_nprim+1)*(G_nprim+2))/2
            associate(this => cont(ij))
              j = int(sqrt(0.25d0 + 2.d0*(ij-1)) - safe_nprim)
              i = (ij-1) - (j*(j+1))/2
              v2 = i + j
              select case(mod(v2,2))
              case(0)
                 v2 = v2/2+1
                 this%par = 0
                 this%n   = this%n_0
                 call product2_se_unit(&
                      this%n,aux4C%se_00,aux4C%val_00,v2,prod2H%unit(i,j),&
                      this%val)
              case(1)
                 v2 = (v2+1)/2
                 this%par = 1
                 this%n   = this%n_1
                 call product2_se_unit(&
                      this%n,aux4C%se_11,aux4C%val_11,v2,prod2H%unit(i,j),&
                      this%val)
              end select
            end associate
         enddo
!$OMP END PARALLEL DO
      case(1)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(i,j,v2)
         do ij=1,((G_nprim+1)*(G_nprim+2))/2
            associate(this => cont(ij))
              j = int(sqrt(0.25d0 + 2.d0*(ij-1)) - safe_nprim)
              i = (ij-1) - (j*(j+1))/2
              v2 = i + j
              select case(mod(v2,2))
              case(0)
                 v2 = v2/2+1
                 this%par = 1
                 this%n   = this%n_1
                 call product2_se_unit(&
                      this%n,aux4C%se_10,aux4C%val_10,v2,prod2H%unit(i,j),&
                      this%val)
              case(1)
                 v2 = (v2+1)/2
                 this%par = 0
                 this%n   = this%n_0
                 call product2_se_unit(&
                      this%n,aux4C%se_01,aux4C%val_01,v2,prod2H%unit(i,j),&
                      this%val)
              end select
            end associate
         enddo
!$OMP END PARALLEL DO
      end select

      call table4_free(aux4C)

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(i,j,k,l,v1,v2) &
!$OMP PRIVATE(val_A1,val_A2,val_B1,val_B2,val_A,val_B)
      do kl=1,G_npair
         do ij=1,G_npair

            i = pair_list(1,ij)
            j = pair_list(2,ij)
            k = pair_list(1,kl)
            l = pair_list(2,kl)

            v1 = max(i,k)
            v1 = 1 + min(i,k) + (v1*(v1+1))/2
            v2 = j + l
            if(cont(v1)%par==mod(v2,2)) then
               v2 = max(j,l)
               v2 = 1 + min(j,l) + (v2*(v2+1))/2
               val_A1 = product1(cont(v1)%n,cont(v1)%val,ints(:,v2))
            else
               val_A1 = 0
            endif

            v1 = max(j,l)
            v1 = 1 + min(j,l) + (v1*(v1+1))/2
            v2 = i + k
            if(cont(v1)%par==mod(v2,2)) then
               v2 = max(i,k)
               v2 = 1 + min(i,k) + (v2*(v2+1))/2
               val_A2 = product1(cont(v1)%n,cont(v1)%val,ints(:,v2))
            else
               val_A2 = 0
            endif

            v1 = max(i,l)
            v1 = 1 + min(i,l) + (v1*(v1+1))/2
            v2 = j + k
            if(cont(v1)%par==mod(v2,2)) then
               v2 = max(j,k)
               v2 = 1 + min(j,k) + (v2*(v2+1))/2
               val_B1 = product1(cont(v1)%n,cont(v1)%val,ints(:,v2))
            else
               val_B1 = 0
            endif

            v1 = max(j,k)
            v1 = 1 + min(j,k) + (v1*(v1+1))/2
            v2 = i + l
            if(cont(v1)%par==mod(v2,2)) then
               v2 = max(i,l)
               v2 = 1 + min(i,l) + (v2*(v2+1))/2
               val_B2 = product1(cont(v1)%n,cont(v1)%val,ints(:,v2))
            else
               val_B2 = 0
            endif

            val_A = val_A1 + val_A2
            val_B = val_B1 + val_B2
            matJ_S(ij,kl,iocc,jocc) = G_gfac*(val_A + val_B)
            matJ_T(ij,kl,iocc,jocc) = G_gfac*(val_A - val_B)
            if(iocc/=jocc) then
               matJ_S(ij,kl,jocc,iocc) = matJ_S(ij,kl,iocc,jocc)
               matJ_T(ij,kl,jocc,iocc) = matJ_T(ij,kl,iocc,jocc)
            endif

         enddo
      enddo
!$OMP END PARALLEL DO

      if(IPRINT>=1) call timer('J',Tcpu,Twall)

   enddo
enddo

ij = 0
do j=0,G_nprim
   do i=0,j
      ij = ij + 1
      associate(this => cont(ij))
        deallocate(this%val)
      end associate
   enddo
enddo
deallocate(cont)

deallocate(sum_0,sum_1)
deallocate(CC)
deallocate(pair_list)

deallocate(ints)
call prod2H_free(prod2H)

end subroutine integrals_J
