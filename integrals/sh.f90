subroutine integrals_SH(matS_S,matS_T,matH_S,matH_T,IPRINT)
implicit none
type(mp_real),intent(out) :: matS_S(:,:),matS_T(:,:)
type(mp_real),intent(out) :: matH_S(:,:),matH_T(:,:)
integer,intent(in) :: IPRINT
real(8) :: Tcpu,Twall
integer :: n1,n2,iOP,i,j,k,l,ij,kl,ijkl,v1,v2
real(8) :: safe_nprim
type(mp_real) :: alpha_invsqrt,a1,a2,a12,c0,c2,c4,val_A,val_B
type(prod2HData) :: prod2H
type(table2Data) :: aux,auxOP
type(mp_real),allocatable :: ints(:,:)
integer,allocatable :: pair_list(:,:)

if(IPRINT>=1) call timer('START',Tcpu,Twall)

safe_nprim = 0.5d0 - 0.5d0/(G_nprim+1)
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

!$OMP PARALLEL PRIVATE(i,j,k,l,ijkl,v1,v2,val_A,val_B)

!$OMP DO SCHEDULE(DYNAMIC)
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
!$OMP END DO

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
   do kl=1,G_npair
      do ij=1,G_npair

         i = pair_list(1,ij)
         j = pair_list(2,ij)
         k = pair_list(1,kl)
         l = pair_list(2,kl)
         ijkl = i + j + k + l

         if(mod(ijkl,2)==0) then

            v1 = i + k
            v1 = (v1+2-mod(v1,2))/2
            v2 = max(j,l)
            v2 = 1 + min(j,l) + (v2*(v2+1))/2
            val_A = product1_unit(v1,prod2H%unit(i,k),ints(:,v2))

            v1 = i + l
            v1 = (v1+2-mod(v1,2))/2
            v2 = max(j,k)
            v2 = 1 + min(j,k) + (v2*(v2+1))/2
            val_B = product1_unit(v1,prod2H%unit(i,l),ints(:,v2))

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

if(IPRINT>=1) call timer('SH',Tcpu,Twall)

end subroutine integrals_SH
