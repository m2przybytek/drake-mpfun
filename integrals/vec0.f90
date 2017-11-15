subroutine integrals_vec0(vec0,matC,parC,IPRINT)
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

if(IPRINT>=1) call timer('START',Tcpu,Twall)

n1 = G_nbas-1
n2 = G_nprim
call prod2H_create(prod2H,n1,n2)
call prod2H_norm(prod2H)

n1 = G_nprim + (G_nbas-1)
n2 = G_nprim + (G_nbas-1)
a1 =  mpreal(-0.5d0)
a2 =  mpreal(-0.5d0)
a12 = mpreal(1.d0)
call table2_create(aux,n1,n2,a1,a2,a12)

int0 = sqrt(2*mppi())/(1 + 2*G_alpha)
c0   = 1/(1 + 2*G_alpha)
mult = mpreal(-1.d0)
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
