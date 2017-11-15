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
