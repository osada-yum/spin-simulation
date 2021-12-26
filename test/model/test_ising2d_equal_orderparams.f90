program test_ising2d_equal_orderparams
  use, intrinsic :: iso_fortran_env
  use builtin_rand_m
  use ising2d_m
  implicit none
  type(Ising2d)              :: system
  type(builtin_rand_wrapper) :: rng
  real(real64)               :: magne, energy
  real(real64)               :: magne_calc, energy_calc

  system = Ising2d(2.269_rkind, 101, 100)
  call rng%set_seed(42)

  call system%set_updater("Metropolis")
  call system%set_order_spin()
  call system%set_kbt(0.1_rkind)

  call system%update_one_mcs(rng)

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

contains

end program test_ising2d_equal_orderparams
