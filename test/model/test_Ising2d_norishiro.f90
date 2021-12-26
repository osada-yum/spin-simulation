program test_Ising2d_norishiro
  use, intrinsic :: iso_fortran_env
  use builtin_rand_m
  use ising2d_m
  implicit none
  type(Ising2d)              :: system
  type(builtin_rand_wrapper) :: rng

  system = Ising2d(2.1_rkind, 201, 200)
  call rng%set_seed(42)

  call system%set_updater("Metropolis")
  call system%set_random_spin(rng)

  block
     integer, parameter         :: n = 100
     integer                    :: i
     do i = 1, n
        call system%update_one_mcs(rng)
        call check_norishiro(system)
     end do
     do i = 1, n
        call system%update_one_mcs(rng)
        call check_norishiro(system)
     end do
   end block

contains

  !> check norishiro's equal correspond parts.
  subroutine check_norishiro(sys)
    type(Ising2d), intent(in) :: sys
    integer                    :: i
    integer                    :: bot, top, b_nori, t_nori
    bot    = 0
    b_nori = 0 - sys%offset()
    top    = sys%particles() - sys%offset()
    t_nori = sys%particles()
    do i = 1, sys%x()-1
       if (sys%spin(i+top) /= sys%spin(i+b_nori)) then
          write(error_unit, '(a, i0, a, i0, a)') "spin(", i+top, ") /= spin(", i+b_nori, ")"
          write(error_unit, '(i0, a, i0)') sys%spin(i+top), " /= ", sys%spin(i+b_nori)
          call util_error_stop("Bottom norishiro does not work"&
               , __LINE__, __FILE__)
       end if
       if (sys%spin(i+bot) /= sys%spin(i+t_nori)) then
          write(error_unit, '(a, i0, a, i0, a)') "spin(", i+bot, ") /= spin(", i+t_nori, ")"
          write(error_unit, '(i0, a, i0)') sys%spin(i+bot), " /= ", sys%spin(i+t_nori)
          call util_error_stop("Top norishiro does not work"&
               , __LINE__, __FILE__)
       end if
    end do
  end subroutine check_norishiro

end program test_Ising2d_norishiro
