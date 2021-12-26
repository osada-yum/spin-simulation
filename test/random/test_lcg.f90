program test_random
  use, intrinsic :: iso_fortran_env
  use utility_m
  implicit none
  integer       , parameter :: num = 100000
  integer(int64), parameter :: maxint32 = huge(0_int32)
  real(real64)  , parameter :: to_real64 = 1 / (real(maxint32, real64) + 1)
  integer                 :: myseed, i
  real(real64)            :: myrnd, summation, summation2

  summation  = 0.0_real64
  summation2 = 0.0_real64
  myseed = 1
  do i = 1, num
     call gen_LCG(myseed, myrnd)
     call test_rangeof_LCG(myrnd)
     summation  = summation  + myrnd
     summation2 = summation2 + myrnd*myrnd
  end do

  !> check average and variance.
  block
    real(real64) :: average, average2, variance
    average  = summation  / num
    average2 = summation2 / num
    variance = average2 - average*average
    if (abs(average - 1/2.0_real64) > 0.1_real64) then
       write(error_unit, '(a, es30.15)') "LCG, average: ", average
       call util_error_stop("Average(LCG) is not near 1/2"&
            , __LINE__, __FILE__)
    end if
    if (abs(average2 - 1/3.0_real64) > 0.1_real64) then
       write(error_unit, '(a, es30.15)') "LCG, average^2: ", average2
       call util_error_stop("Average^2(LCG) is not near 1/3"&
            , __LINE__, __FILE__)
    end if
    if (abs(variance - 1/12.0_real64) > 0.1_real64) then
       write(error_unit, '(a, es30.15)') "LCG, variance: ", variance
       call util_error_stop("Variance(LCG) is not near 1/12"&
            , __LINE__, __FILE__)
    end if
  end block

contains

  !> check reange of random number by LCG.
  subroutine test_rangeof_LCG(rnd)
    real(real64), intent(in) :: rnd
    if (rnd < 0.0_rkind .or. rnd >= 1.0_rkind) then
       write(error_unit, '(a, es30.15)') "LCG: ", rnd
       call util_error_stop("Range of LCG is out of [0,1)"&
            , __LINE__, __FILE__)
    end if
  end subroutine test_rangeof_LCG

  !> generage random number by LCG.
  subroutine gen_LCG(seed, rnd)
    integer     , intent(inout) :: seed
    real(real64), intent(out)   :: rnd
    integer(int64)              :: longseed
    longseed = int(seed*48271, int64)
    seed     = int(iand(longseed, z'FFFFFFFFFFFFFFFF'), int32)
    rnd      = (seed*to_real64+1)/2 ! [0, 1)
  end subroutine gen_LCG

end program test_random
