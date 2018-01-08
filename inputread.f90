module inputread
  use mpmodule
  use general
  implicit none

contains

  subroutine read_Input()
    implicit none
    double precision :: pre_g,pre_a
    open(10,file="input.in",status="old")
    read(10,*) G_nocc
    read(10,*) G_nbas
    read(10,*) G_nprim
    read(10,*) G_naux
    read(10,*) pre_g
    read(10,*) pre_a
    G_gfac = mpreal(pre_g)
    G_alpha = mpreal(pre_a)
    close(10)
    write(*,*) "---------------------------------------------------------"
    write(*,*) "| Some input reading has happened                       |"
    write(*,*) "---------------------------------------------------------"
    write(*,'(A11,F15.2,A32)') " | g     = ",pre_g, "|"
    write(*,'(A11,I15,A32)') "| nbas  = ",G_nbas,"|"
    write(*,'(A11,I15,A32)') "| nocc  = ",G_nocc,"|"
    write(*,'(A11,I15,A32)') "| naux  = ",G_naux,"|"
    write(*,'(A11,I15,A32)') "| nprim = ",G_nprim,"|"
    write(*,'(A11,F15.2,A32)') "| alpha = ",pre_a,"|"
    write(*,*) "---------------------------------------------------------"
  end subroutine read_Input

end module inputread
