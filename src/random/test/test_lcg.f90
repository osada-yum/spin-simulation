program test_random
  use, intrinsic :: iso_fortran_env
  implicit none
  integer       , parameter :: num = 10000
  integer(int64), parameter :: maxint32 = huge(1_int32)
  real(real64)  , parameter :: to_real64 = 1 / (real(maxint32, real64) + 1)
  integer                 :: myseed, i
  real(real64)            :: myrnd, summation

  print *, maxint32
  summation = 0.0_real64
  myseed = 1
  do i = 1, num
     call lcg(myseed, myrnd)
     summation = summation + myrnd
     print*, myrnd
  end do
  print*, "average: ", summation / num

contains
  subroutine lcg(seed, rnd)
    integer     , intent(inout) :: seed
    real(real64), intent(out)   :: rnd
    integer(int64)              :: longseed
    longseed = int(seed*48271, int64)
    seed     = int(iand(longseed, maxint32), int32)
    rnd      = (seed*to_real64+1)/2 ! [0, 1)
  end subroutine lcg

end program test_random
