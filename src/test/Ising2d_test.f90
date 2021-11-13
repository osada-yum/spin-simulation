program test_Ising2d
  use, intrinsic :: iso_fortran_env
  use ising2d_m
  implicit none
  type(Ising2d) :: system
  real(real64)   :: magne, energy

  system = Ising2d(1.0_rkind, 201, 200)

  call system%set_order_spin()
  call calc_magne_and_energy(system, magne, energy)
  if (abs(magne-1.0_rkind) > epsilon) error stop 1
  if (abs(energy-(-2.0_rkind)) > epsilon) error stop 1

contains

end program test_Ising2d
