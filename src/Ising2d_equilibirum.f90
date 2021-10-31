program Ising2d_equilibrium
  use,intrinsic :: iso_fortran_env
  use ising2d_m
  implicit none
  integer, parameter              :: rkind = real64
  type(Ising2d)                   :: system
  real(rkind), allocatable        :: temperature(:)
  real(rkind), parameter          :: temperature_begin = 1.7, temperature_end = 2.4
  integer    , parameter          :: num_temperature = 100
  real(rkind), allocatable        :: magne(:), energy(:)
  integer, parameter              :: relx_mcs = 1000, sample_mcs = 1000
  integer                         :: i, j

  system = Ising2d(1.0_rkind, 201, 200)

  block
    integer              :: seedsize
    integer, allocatable :: seed(:)
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    do i = 1, seedsize
       seed(i) = i
    end do
    call random_seed(put=seed)
  end block

  allocate( magne(num_temperature), energy(num_temperature), source = 0.0_rkind)
  allocate( temperature(num_temperature) )
  do i = 1, num_temperature
     temperature(i) = &
          ( (i-1)                  *temperature_begin&
          + (num_temperature-1-i+1)*temperature_end    ) / (num_temperature-1) ! 内分点.
  end do
  write(output_unit, '(a,i0)'     ) "# N = ", system%particles()
  write(output_unit, '(2(a,i0))'  ) "# MCS(relaxation) = ", relx_mcs, " MCS(sampling) = ", sample_mcs
  write(error_unit , '(3(a, i0))' ) "nx", system%x(), " ny", system%y()
  write(error_unit , '(a, f30.18)') "method: METROPOLIS"

  !! 初期配置(ランダム).
  call system%set_random_spin()

  do j = 1, num_temperature
     call system%set_kbt(temperature(j))
     !! 空回し
     do i = 1, relx_mcs
        call system%update_with_Metropolis_one_mcs()
     end do
     do i = 1, sample_mcs
        call system%update_with_Metropolis_one_mcs()
        block
          real(rkind) :: magne_tmp, energy_tmp
          call calc_magne_and_energy(system, magne_tmp, energy_tmp)
          magne(j)  = magne(j)  + magne_tmp
          energy(j) = energy(j) + energy_tmp
        end block
     end do
  end do

  block
    real(rkind) :: magne_mean, energy_mean
    do i = 1, num_temperature
       magne_mean  = magne(i)  / sample_mcs
       energy_mean = energy(i) / sample_mcs
       write(output_unit,'(i0,a,i0, " ", *(es20.12))') &
            system%x(), "x", system%y(),&
            temperature(i), &
            magne_mean, energy_mean
    end do
  end block

  deallocate( magne, energy )

  stop "END PROGRAM"

contains

end program Ising2d_equilibrium
