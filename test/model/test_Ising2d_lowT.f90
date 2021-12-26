program test_Ising2d_lowT
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
  call system%set_order_spin()
  call system%set_kbt(0.1_rkind)

  !> check magne and energy is ground state.
  if (abs(system%magne()-1.0_rkind) > epsilon) then
     write(error_unit, '(2(es30.15), a)') magne, system%magne(), "/= 1.0"
     call util_error_stop("magne /= 1"&
          , __LINE__, __FILE__)
  end if
  if (abs(system%energy()-(-2.0_rkind)) > epsilon) then
     write(error_unit, '(2(es30.15), a)') energy, system%energy(), "/= -2.0"
     call util_error_stop("energy /= -2"&
          , __LINE__, __FILE__)
  end if

  block
    integer, parameter         :: n = 500
    integer                    :: i
    magne  = 0.0_real64
    energy = 0.0_real64
    do i = 1, n
       call system%update_one_mcs(rng)
    end do
    do i = 1, n
       call system%update_one_mcs(rng)
       magne  = magne  + system%magne()
       energy = energy + system%energy()
    end do
    magne  = magne / n
    energy = energy / n
  end block

  call calc_magne_and_energy(system, magne_calc, energy_calc)
  !> check calc_magne_and_energy == system%mange()m system%energy().
  if (abs(magne_calc-system%magne()) > epsilon) then
     write(error_unit, '(es30.15, a, es30.15)') magne, "/=", system%magne()
     call util_error_stop("magne /= s%m()"&
          , __LINE__, __FILE__)
  end if
  if (abs(energy_calc-system%energy()) > epsilon) then
     write(error_unit, '(es30.15, a, es30.15)') energy_calc, "/=", system%energy()
     call util_error_stop("energy /= s%e()"&
          , __LINE__, __FILE__)
  end if
  !> check magne and energy is near ground state.
  if (abs(magne-1.0_rkind) > epsilon) then
     write(error_unit, '(2(es30.15), a)') magne, system%magne(), "/= 1.0"
     call util_error_stop("magne /= 1"&
          , __LINE__, __FILE__)
  end if
  if (abs(energy-(-2.0_rkind)) > epsilon) then
     write(error_unit, '(2(es30.15), a)') energy, system%energy(), "/= -2.0"
     call util_error_stop("energy /= -2"&
          , __LINE__, __FILE__)
  end if

contains

end program test_Ising2d_lowT
