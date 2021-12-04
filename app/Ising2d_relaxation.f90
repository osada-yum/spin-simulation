#ifndef RAND_GEN_TYPE
#define RAND_GEN_TYPE builtin_rand_wrapper
#endif

program Ising2d_relaxation
  use, intrinsic :: iso_fortran_env
  use xorshift_m
  use builtin_rand_m
  use ising2d_m
  use benchmark_m
  implicit none
  type(Ising2d)                   :: system
  type(RAND_GEN_TYPE)             :: gen
  real(rkind), allocatable        :: magne(:), energy(:)
  integer, parameter              :: mcs = 1000, sample = 10
  integer                         :: i, j
  type(benchmark_t)               :: bm
  bm = benchmark_t()

  system = Ising2d(2.269_rkind, 2001, 2000)

  call gen%set_seed(42)

  allocate( magne(mcs), energy(mcs), source = 0.0_rkind)
  write(output_unit, '(a,i0)'      ) "# N = ", system%particles()
  write(output_unit, '(a,f20.10)'  ) "# kbt =", system%kbt()
  write(error_unit , '(3(a, i0))'  ) "nx", system%x(), " ny", system%y(), " MCS: ", mcs
  write(error_unit , '(a, f30.18)' ) "温度: ", system%kbt()
  write(error_unit , '(a, f30.18)' ) "method: METROPOLIS"

  call bm%stamp("start update")
  do j = 1, sample
     !! 初期配置.
     call system%set_order_spin()
     write(error_unit, *) "sample: ", j
     do i = 1, mcs
        call system%update_with_Metropolis_one_mcs(gen)
        magne(i)  = magne(i)  + calc_magne(system)
        energy(i) = energy(i) + calc_energy(system)
     end do
  end do
  call bm%stamp("end update")

  !! 平均とって出力.
  magne(:)  =  magne(:) / sample
  energy(:) = energy(:) / sample
  do i = 1, mcs
     write(output_unit,'(i0,a,i0, " ", i0, " ", *(es20.12))') &
          system%x(), "x", system%y(),&
          i, &
          magne(i), energy(i)
  end do

  deallocate( magne, energy )

  stop "END PROGRAM"

contains

end program Ising2d_relaxation
