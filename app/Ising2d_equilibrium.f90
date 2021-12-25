#ifndef RAND_GEN_TYPE
#define RAND_GEN_TYPE builtin_rand_wrapper
#endif

program Ising2d_equilibrium
  use, intrinsic :: iso_fortran_env
  use utility_m
  use xorshift_m
  use builtin_rand_m
  use ising2d_m
  use benchmark_m
  implicit none
  type(Ising2d)                   :: system
  type(RAND_GEN_TYPE)             :: gen
  real(rkind), allocatable        :: temperature(:)
  real(rkind), parameter          :: temperature_begin = 1.7_rkind, temperature_end = 2.4_rkind
  integer    , parameter          :: num_temperature = 100
  real(rkind), allocatable        :: magne(:), energy(:)
  integer, parameter              :: relx_mcs = 1000, sample_mcs = 1000
  integer                         :: i, j
  type(benchmark_t)               :: bm
  bm = benchmark_t()

  system = Ising2d(1.0_rkind, 201, 200)

  call gen%set_seed(42)

  allocate( magne(num_temperature) , source = 0.0_rkind)
  allocate( energy(num_temperature), source = 0.0_rkind)
  allocate( temperature(num_temperature) )
  temperature(:) = util_linspace(temperature_begin, temperature_end, num_temperature)

  write(output_unit, '(a,i0)'     ) "# N = ", system%particles()
  write(output_unit, '(2(a,i0))'  ) "# MCS(relaxation) = ", relx_mcs, " MCS(sampling) = ", sample_mcs
  write(error_unit , '(3(a, i0))' ) "nx", system%x(), " ny", system%y()
  write(error_unit , '(a, f30.18)') "method: METROPOLIS"

  !! 初期配置(ランダム).
  call system%set_random_spin(gen)

  call system%set_updater("Metropolis")

  call bm%stamp("start update")
  do j = 1, num_temperature
     call system%set_kbt(temperature(j))
     write(error_unit, '(a, i7, es23.15)') "iterate: ", j, temperature(j)
     !! 空回し
     do i = 1, relx_mcs
        call system%update_one_mcs(gen)
     end do
     do i = 1, sample_mcs
        call system%update_one_mcs(gen)
        calc_order_parameters: block
          real(rkind) :: m_tmp, e_tmp
          call calc_magne_and_energy(system, m_tmp, e_tmp)
          magne(j)  = magne(j)  + abs(m_tmp)
          energy(j) = energy(j) + e_tmp
        end block calc_order_parameters
     end do
  end do
  call bm%stamp("end update")

  print_order_parameter: block
    real(rkind) :: magne_mean, energy_mean
    do i = 1, num_temperature
       magne_mean  = magne(i)  / sample_mcs
       energy_mean = energy(i) / sample_mcs
       write(output_unit,'(i0,a,i0, " ", *(es20.12))') &
            system%x(), "x", system%y(),&
            temperature(i), &
            magne_mean, energy_mean
    end do
  end block print_order_parameter
  call bm%dump()
  call destroy_benchmark_t(bm)

  deallocate( magne, energy )

  stop "END PROGRAM"

contains

end program Ising2d_equilibrium
