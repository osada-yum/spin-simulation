program test_Ising2d_element
  use, intrinsic :: iso_fortran_env
  use builtin_rand_m
  use ising2d_m
  implicit none
  type(Ising2d)              :: system
  type(builtin_rand_wrapper) :: rng
  real(real64)               :: magne, energy
  real(real64)               :: magne_calc, energy_calc

  system = Ising2d(0.1_rkind, 101, 100)
  call rng%set_seed(42)

  call system%set_updater("Metropolis")
  call system%set_random_spin(rng)
  call is_all_element_one(system)

  call system%set_order_spin()
  call is_all_element_one(system)

contains

  !> check all elements is 1 or -1.
  subroutine is_all_element_one(sys)
    type(Ising2d), intent(in) :: sys
    integer                   :: i
    do i = 1, sys%particles()
       if (abs(sys%spin(i)) /= 1) then
          write(error_unit, '(2(es30.15), a)') i, sys%spin(i)
          call util_error_stop("sys%spin(i) /= 1"&
               , __LINE__, __FILE__)
       end if
    end do
  end subroutine is_all_element_one

end program test_Ising2d_element
