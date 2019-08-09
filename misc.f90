module misc
use mpmodule
use global
implicit none

private
public prodUnitData
public product1,product1_se,product1_unit
public product2_se,product2_unit,product2_se_unit
public product23_se
public accumulate,accumulate_unit,addition
public FDS_SDF,fullMO_transform
public zero_vec,scal_vec,copy_vec,add_vec,dot_vec
public zero_mat,scal_mat,copy_mat,add_mat,add_matT,dot_mat
public outer,add_outer,matvec,add_matvec,vecmatvec
public TripletData

type prodUnitData
integer :: nu
type(mp_real),allocatable :: val(:)
end type prodUnitData

type TripletData
integer :: npair,npair_orig
integer,allocatable :: idxS(:),idxE(:)
contains
procedure :: shrink_Triplet_mat
procedure :: shrink_Triplet_mat_dim
procedure :: shrink_Triplet_vec
procedure :: expand_Triplet_mat
procedure :: expand_Triplet_vec
procedure :: init => init_Triplet
procedure :: free => free_Triplet
generic :: shrink => &
     shrink_Triplet_mat,shrink_Triplet_mat_dim,shrink_Triplet_vec
generic :: expand => &
     expand_Triplet_mat,expand_Triplet_vec
end type TripletData

contains

function product1(n,vec1,vec2) result(val)
implicit none
type(mp_real) :: val
integer,intent(in) :: n
type(mp_real),intent(in) :: vec1(:),vec2(:)
integer :: i

val = 0
do i=1,n
   val = val + vec1(i)*vec2(i)
enddo

end function product1

function product1_se(se,vec1,vec2) result(val)
implicit none
type(mp_real) :: val
integer,intent(in) :: se(:)
type(mp_real),intent(in) :: vec1(:),vec2(:)
integer :: i

val = 0
do i=se(1),se(2)
   val = val + vec1(i)*vec2(i)
enddo

end function product1_se

function product1_unit(v,unit,vec) result(val)
implicit none
type(mp_real) :: val
integer,intent(in) :: v
type(prodUnitData),intent(in) :: unit
type(mp_real),intent(in) :: vec(:)
integer :: off,i

off = v - unit%nu

val = unit%val(1)*vec(off)
do i=1,unit%nu
   val = val + unit%val(1+i)*vec(off+i)
enddo

end function product1_unit

subroutine product2_se(n,mat,se,cont,vec)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: mat(:,:)
integer,intent(in) :: se(:)
type(mp_real),intent(in) :: cont(:)
type(mp_real),intent(out) :: vec(:)
integer :: i1,i2

do i1=1,n
   vec(i1) = 0
enddo

do i2=se(1),se(2)
   do i1=1,n
      vec(i1) = vec(i1) + mat(i1,i2)*cont(i2)
   enddo
enddo

end subroutine product2_se

subroutine product2_unit(n,mat,v,unit,vec)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: mat(:,:)
integer,intent(in) :: v
type(prodUnitData),intent(in) :: unit
type(mp_real),intent(out) :: vec(:)
integer :: off2,i1,i2

off2 = v - unit%nu

do i1=1,n
   vec(i1) = mat(i1,off2)*unit%val(1)
enddo
do i2=1,unit%nu
   do i1=1,n
      vec(i1) = vec(i1) + mat(i1,off2+i2)*unit%val(1+i2)
   enddo
enddo

end subroutine product2_unit

subroutine product2_se_unit(n,se,mat,v,unit,vec)
implicit none
integer,intent(in) :: n
integer,intent(in) :: se(:,:)
type(mp_real),intent(in) :: mat(:,:)
integer,intent(in) :: v
type(prodUnitData),intent(in) :: unit
type(mp_real),intent(out) :: vec(:)
integer :: off2,i1,i2

do i1=1,n
   vec(i1) = 0
enddo

off2 = v - unit%nu

do i2=0,unit%nu
   do i1=se(1,off2+i2),se(2,off2+i2)
      vec(i1) = vec(i1) + mat(i1,off2+i2)*unit%val(1+i2)
   enddo
enddo

end subroutine product2_se_unit

subroutine product23_se(n,mat,se2,cont2,se3,cont3,vec)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: mat(:,:,:)
integer,intent(in) :: se2(:)
type(mp_real),intent(in) :: cont2(:)
integer,intent(in) :: se3(:)
type(mp_real),intent(in) :: cont3(:)
type(mp_real),intent(out) :: vec(:)
integer :: i1,i2,i3
type(mp_real) :: val

do i1=1,n
   vec(i1) = 0
enddo

do i3=se3(1),se3(2)
   do i2=se2(1),se2(2)
      val = cont2(i2)*cont3(i3)
      do i1=1,n
         vec(i1) = vec(i1) + mat(i1,i2,i3)*val
      enddo
   enddo
enddo

end subroutine product23_se

!-------------------------------------------------------------------------------

subroutine accumulate(alpha,n,vec,acc)
implicit none
type(mp_real),intent(in) :: alpha
integer,intent(in) :: n
type(mp_real),intent(in) :: vec(:)
type(mp_real),intent(inout) :: acc(:)
integer :: i

do i=1,n
   acc(i) = acc(i) + alpha*vec(i)
enddo

end subroutine accumulate

subroutine accumulate_unit(alpha,v,unit,acc)
implicit none
type(mp_real),intent(in) :: alpha
integer,intent(in) :: v
type(prodUnitData),intent(in) :: unit
type(mp_real),intent(inout) :: acc(:)
integer :: off,i

off = v - unit%nu

do i=0,unit%nu
   acc(off+i) = acc(off+i) + alpha*unit%val(1+i)
enddo

end subroutine accumulate_unit

subroutine addition(n,vec,acc)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: vec(:)
type(mp_real),intent(inout) :: acc(:)
integer :: i

do i=1,n
   acc(i) = acc(i) + vec(i)
enddo

end subroutine addition

!-------------------------------------------------------------------------------

subroutine FDS_SDF(n,error,P,F,D,S,work)
implicit none
integer,intent(in) :: n
type(mp_real),intent(out) :: error
type(mp_real),intent(out) :: P(:,:)
type(mp_real),intent(in) :: F(:,:),D(:,:),S(:,:)
type(mp_real) :: work(:,:)
integer :: i,j,k
type(mp_real) :: val

! FDS = ( D^ F )^ S  for F^ = F

error = 0

!$OMP PARALLEL PRIVATE(val)

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      val = 0
      do k=1,n
         val = val + D(k,i)*F(k,j)
      enddo
      work(i,j) = val
   enddo
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      val = 0
      do k=1,n
         val = val + work(k,i)*S(k,j)
      enddo
      P(i,j) = val
   enddo
enddo
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC)
do j=1,n
   do i=1,j-1
      val = P(i,j)
      P(i,j) = P(i,j) - P(j,i)
      P(j,i) = P(j,i) - val
   enddo
   P(j,j) = 0
enddo
!$OMP END DO

val = 0

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      val = val + P(i,j)**2
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP CRITICAL(FDS_SDF_error)
error = error + val
!$OMP END CRITICAL(FDS_SDF_error)

!$OMP END PARALLEL

error = sqrt(error)

end subroutine FDS_SDF

subroutine fullMO_transform(n,C,A,work)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: C(:,:)
type(mp_real),intent(inout) :: a(:,:)
type(mp_real) :: work(:,:)
integer :: i,j,k
type(mp_real) :: val

!$OMP PARALLEL PRIVATE(val)

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      val = 0
      do k=1,n
         val = val + C(k,i)*A(k,j)
      enddo
      work(i,j) = val
   enddo
enddo
!$OMP END DO

!$OMP DO COLLAPSE(2)
do j=1,n
   do i=1,n
      val = 0
      do k=1,n
         val = val + work(i,k)*C(k,j)
      enddo
      A(i,j) = val
   enddo
enddo
!$OMP END DO

!$OMP END PARALLEL

end subroutine fullMO_transform

subroutine zero_vec(n,vec)
implicit none
integer,intent(in) :: n
type(mp_real),intent(out) :: vec(:)
integer :: i

!$OMP PARALLEL DO
do i=1,n
   vec(i) = 0
enddo
!$OMP END PARALLEL DO

end subroutine zero_vec

subroutine scal_vec(alpha,n,vec)
implicit none
type(mp_real),intent(in) :: alpha
integer,intent(in) :: n
type(mp_real),intent(inout) :: vec(:)
integer :: i

!$OMP PARALLEL DO
do i=1,n
   vec(i) = alpha*vec(i)
enddo
!$OMP END PARALLEL DO

end subroutine scal_vec

subroutine copy_vec(n,IN,OUT)
implicit none
integer,intent(in) :: n
type(mp_real),intent(in) :: IN(:)
type(mp_real),intent(inout) :: OUT(:)
integer :: i

!$OMP PARALLEL DO
do i=1,n
   OUT(i) = IN(i)
enddo
!$OMP END PARALLEL DO

end subroutine copy_vec

subroutine add_vec(alpha,n,IN,OUT)
implicit none
type(mp_real),intent(in) :: alpha
integer,intent(in) :: n
type(mp_real),intent(in) :: IN(:)
type(mp_real),intent(inout) :: OUT(:)
integer :: i

!$OMP PARALLEL DO
do i=1,n
   OUT(i) = OUT(i) + alpha*IN(i)
enddo
!$OMP END PARALLEL DO

end subroutine add_vec

function dot_vec(n,vec1,vec2) result(val)
implicit none
type(mp_real) :: val
integer,intent(in) :: n
type(mp_real),intent(in) :: vec1(:),vec2(:)
type(mp_real) :: part
integer :: i

val = 0

!$OMP PARALLEL PRIVATE(part)

part = 0

!$OMP DO
do i=1,n
   part = part + vec1(i)*vec2(i)
enddo
!$OMP END DO NOWAIT

!$OMP CRITICAL(dot_vec_reduction)
val = val + part
!$OMP END CRITICAL(dot_vec_reduction)

!$OMP END PARALLEL

end function dot_vec

subroutine zero_mat(n1,n2,mat)
implicit none
integer,intent(in) :: n1,n2
type(mp_real),intent(out) :: mat(:,:)
integer :: i1,i2

!$OMP PARALLEL DO COLLAPSE(2)
do i2=1,n2
   do i1=1,n1
      mat(i1,i2) = 0
   enddo
enddo
!$OMP END PARALLEL DO

end subroutine zero_mat

subroutine scal_mat(alpha,n1,n2,mat)
implicit none
type(mp_real),intent(in) :: alpha
integer,intent(in) :: n1,n2
type(mp_real),intent(inout) :: mat(:,:)
integer :: i1,i2

!$OMP PARALLEL DO COLLAPSE(2)
do i2=1,n2
   do i1=1,n1
      mat(i1,i2) = alpha*mat(i1,i2)
   enddo
enddo
!$OMP END PARALLEL DO

end subroutine scal_mat

subroutine copy_mat(n1,n2,IN,OUT)
implicit none
integer,intent(in) :: n1,n2
type(mp_real),intent(in) :: IN(:,:)
type(mp_real),intent(inout) :: OUT(:,:)
integer :: i1,i2

!$OMP PARALLEL DO COLLAPSE(2)
do i2=1,n2
   do i1=1,n1
      OUT(i1,i2) = IN(i1,i2)
   enddo
enddo
!$OMP END PARALLEL DO

end subroutine copy_mat

subroutine add_mat(alpha,n1,n2,IN,OUT)
implicit none
type(mp_real),intent(in) :: alpha
integer,intent(in) :: n1,n2
type(mp_real),intent(in) :: IN(:,:)
type(mp_real),intent(inout) :: OUT(:,:)
integer :: i1,i2

!$OMP PARALLEL DO COLLAPSE(2)
do i2=1,n2
   do i1=1,n1
      OUT(i1,i2) = OUT(i1,i2) + alpha*IN(i1,i2)
   enddo
enddo
!$OMP END PARALLEL DO

end subroutine add_mat

subroutine add_matT(alpha,n1,n2,IN,OUT)
implicit none
type(mp_real),intent(in) :: alpha
integer,intent(in) :: n1,n2
type(mp_real),intent(in) :: IN(:,:)
type(mp_real),intent(inout) :: OUT(:,:)
integer :: i1,i2

!$OMP PARALLEL DO COLLAPSE(2)
do i2=1,n2
   do i1=1,n1
      OUT(i1,i2) = OUT(i1,i2) + alpha*IN(i2,i1)
   enddo
enddo
!$OMP END PARALLEL DO

end subroutine add_matT

function dot_mat(n1,n2,mat1,mat2) result(val)
implicit none
type(mp_real) :: val
integer,intent(in) :: n1,n2
type(mp_real),intent(in) :: mat1(:,:),mat2(:,:)
type(mp_real) :: part
integer :: i1,i2

val = 0

!$OMP PARALLEL PRIVATE(part)

part = 0

!$OMP DO COLLAPSE(2)
do i2=1,n2
   do i1=1,n1
      part = part + mat1(i1,i2)*mat2(i1,i2)
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP CRITICAL(dot_mat_reduction)
val = val + part
!$OMP END CRITICAL(dot_mat_reduction)

!$OMP END PARALLEL

end function dot_mat

!-------------------------------------------------------------------------------

subroutine outer(n1,n2,vec1,vec2,mat)
implicit none
integer,intent(in) :: n1,n2
type(mp_real),intent(in) :: vec1(:),vec2(:)
type(mp_real),intent(out) :: mat(:,:)
integer :: i1,i2

!$OMP PARALLEL DO COLLAPSE(2)
do i2=1,n2
   do i1=1,n1
      mat(i1,i2) = vec1(i1)*vec2(i2)
   enddo
enddo
!$OMP END PARALLEL DO

end subroutine outer

subroutine add_outer(alpha,n1,n2,vec1,vec2,mat)
implicit none
type(mp_real),intent(in) :: alpha
integer,intent(in) :: n1,n2
type(mp_real),intent(in) :: vec1(:),vec2(:)
type(mp_real),intent(inout) :: mat(:,:)
type(mp_real) :: val
integer :: i1,i2

!$OMP PARALLEL DO PRIVATE(val)
do i2=1,n2
   val = vec2(i2)*alpha
   do i1=1,n1
      mat(i1,i2) = mat(i1,i2) + vec1(i1)*val
   enddo
enddo
!$OMP END PARALLEL DO

end subroutine add_outer

subroutine matvec(n1,n2,mat,vec,res)
implicit none
integer,intent(in) :: n1,n2
type(mp_real),intent(in) :: mat(:,:),vec(:)
type(mp_real),intent(inout) :: res(:)
integer :: i1,i2
type(mp_real) :: val

!$OMP PARALLEL DO PRIVATE(val)
do i1=1,n1
   val = 0
   do i2=1,n2
      val = val + mat(i1,i2)*vec(i2)
   enddo
   res(i1) = val
enddo
!$OMP END PARALLEL DO

end subroutine matvec

subroutine add_matvec(alpha,n1,n2,mat,vec,res)
implicit none
type(mp_real),intent(in) :: alpha
integer,intent(in) :: n1,n2
type(mp_real),intent(in) :: mat(:,:),vec(:)
type(mp_real),intent(inout) :: res(:)
integer :: i1,i2
type(mp_real) :: val

!$OMP PARALLEL DO PRIVATE(val)
do i1=1,n1
   val = 0
   do i2=1,n2
      val = val + mat(i1,i2)*vec(i2)
   enddo
   res(i1) = res(i1) + alpha*val
enddo
!$OMP END PARALLEL DO

end subroutine add_matvec

function vecmatvec(n1,n2,mat,vec1,vec2) result(val)
implicit none
type(mp_real) :: val
integer,intent(in) :: n1,n2
type(mp_real),intent(in) :: mat(:,:),vec1(:),vec2(:)
type(mp_real) :: part
integer :: i1,i2

val = 0

!$OMP PARALLEL PRIVATE(part)

part = 0

!$OMP DO COLLAPSE(2)
do i2=1,n2
   do i1=1,n1
      part = part + vec1(i1)*mat(i1,i2)*vec2(i2)
   enddo
enddo
!$OMP END DO NOWAIT

!$OMP CRITICAL(vecmatvec_reduction)
val = val + part
!$OMP END CRITICAL(vecmatvec_reduction)

!$OMP END PARALLEL

end function vecmatvec

!-------------------------------------------------------------------------------

subroutine init_Triplet(Triplet,npair,omega)
implicit none
class(TripletData) :: Triplet
integer,intent(in) :: npair
integer,intent(in) :: omega
integer :: i,j,ipair,ipair_orig

Triplet%npair      = npair
Triplet%npair_orig = npair

allocate(Triplet%idxS(Triplet%npair_orig))
allocate(Triplet%idxE(Triplet%npair_orig))

Triplet%idxS = 0
Triplet%idxE = 0

ipair      = 0
ipair_orig = 0

i = -1; j = -1
do while(next_pair(i,j,omega))
   ipair_orig = ipair_orig + 1
   if(i/=j) then
      ipair = ipair + 1
      Triplet%idxS(ipair)      = ipair_orig
      Triplet%idxE(ipair_orig) = ipair
   endif
enddo

Triplet%npair = ipair

end subroutine init_Triplet

subroutine free_Triplet(Triplet)
implicit none
class(TripletData) :: Triplet

deallocate(Triplet%idxE)
deallocate(Triplet%idxS)

end subroutine free_Triplet

subroutine shrink_Triplet_mat(Triplet,mat)
implicit none
class(TripletData),intent(in) :: Triplet
type(mp_real),intent(inout) :: mat(:,:)
integer :: i,j,j_idx

do j=1,Triplet%npair

   j_idx = Triplet%idxS(j)
   do i=1,Triplet%npair
      mat(i,j) = mat(Triplet%idxS(i),j_idx)
   enddo
   do i=Triplet%npair+1,Triplet%npair_orig
      mat(i,j) = 0
   enddo

enddo

do j=Triplet%npair+1,Triplet%npair_orig
   do i=1,Triplet%npair_orig
      mat(i,j) = 0
   enddo
enddo

end subroutine shrink_Triplet_mat

subroutine shrink_Triplet_mat_dim(Triplet,mat,dim)
implicit none
class(TripletData),intent(in) :: Triplet
type(mp_real),intent(inout) :: mat(:,:)
integer,intent(in) :: dim
integer :: i,j,j_idx

select case(dim)
case(1)

   do j=1,Triplet%npair_orig
      do i=1,Triplet%npair
         mat(i,j) = mat(Triplet%idxS(i),j)
      enddo
      do i=Triplet%npair+1,Triplet%npair_orig
         mat(i,j) = 0
      enddo
   enddo

case(2)

   do j=1,Triplet%npair
      j_idx = Triplet%idxS(j)
      do i=1,Triplet%npair_orig
         mat(i,j) = mat(i,j_idx)
      enddo
   enddo
   do j=Triplet%npair+1,Triplet%npair_orig
      do i=1,Triplet%npair_orig
         mat(i,j) = 0
      enddo
   enddo

case default
   write(*,'(a)') 'ERROR!!! Incorrect dimension in Triplet%shrink!'
   stop
end select
   
end subroutine shrink_Triplet_mat_dim

subroutine shrink_Triplet_vec(Triplet,vec)
implicit none
class(TripletData),intent(in) :: Triplet
type(mp_real),intent(inout) :: vec(:)
integer :: i

do i=1,Triplet%npair
   vec(i) = vec(Triplet%idxS(i))
enddo
do i=Triplet%npair+1,Triplet%npair_orig
   vec(i) = 0
enddo

end subroutine shrink_Triplet_vec

subroutine expand_Triplet_mat(Triplet,mat)
implicit none
class(TripletData),intent(in) :: Triplet
type(mp_real),intent(inout) :: mat(:,:)
integer :: i,j,j_idx

do j=Triplet%npair_orig,1,-1
   if(Triplet%idxE(j)>0) then

      j_idx = Triplet%idxE(j)
      do i=Triplet%npair_orig,1,-1
         if(Triplet%idxE(i)>0) then
            mat(i,j) = mat(Triplet%idxE(i),j_idx)
         else
            mat(i,j) = 0
         endif
      enddo

   else

      do i=1,Triplet%npair_orig
         mat(i,j) = 0
      enddo

   endif
enddo

end subroutine expand_Triplet_mat

subroutine expand_Triplet_vec(Triplet,vec)
implicit none
class(TripletData),intent(in) :: Triplet
type(mp_real),intent(inout) :: vec(:)
integer :: i

do i=Triplet%npair_orig,1,-1
   if(Triplet%idxE(i)>0) then
      vec(i) = vec(Triplet%idxE(i))
   else
      vec(i) = 0
   endif
enddo

end subroutine expand_Triplet_vec

end module misc
