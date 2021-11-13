module builtin_rand_m
  use, intrinsic :: iso_fortran_env
  use random_base_m
  implicit none
  private
  integer       , parameter :: rkind = real64
  integer(int64), parameter :: maxint32 = huge(1_int32)
  real(real64)  , parameter :: to_real64 = 1 / (real(huge(1_int32), real64) + 1)

  public :: builtin_rand_wrapper
  type, extends(random_base_t) :: builtin_rand_wrapper ! builtin_rand_wrapper for fortran builtin random number generator
     private
   contains
     procedure, pass :: set_seed   => set_seed_wrapper
     procedure, pass :: get_rand   => get_rand_wrapper
     procedure, pass :: random_arr => random_arr_wrapper
  end type builtin_rand_wrapper

contains

  subroutine set_seed_wrapper(this, iseed)
    class(builtin_rand_wrapper), intent(inout) :: this
    integer                    , intent(in)    :: iseed
    integer                                    :: i, myseed, seedsize
    integer, allocatable                       :: seed(:)
    real(rkind)                                :: dummy_rnd
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    myseed = iseed
    do i = 1, seedsize
       seed(i) = myseed
       call lcg(myseed, dummy_rnd)
    end do
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine set_seed_wrapper

  real(real64) function get_rand_wrapper(this)
    class(builtin_rand_wrapper), intent(inout) :: this
    real(real64)                               :: rnd
    call random_number(rnd)
    get_rand_wrapper = rnd
    return
  end function get_rand_wrapper

  subroutine random_arr_wrapper(this, arr)
    class(builtin_rand_wrapper), intent(inout) :: this
    real(real64)               , intent(out)   :: arr(:)
    integer                                    :: i
    call random_number(arr)
  end subroutine random_arr_wrapper

  subroutine lcg(seed, rnd)
    integer     , intent(inout) :: seed
    real(rkind), intent(out)    :: rnd
    integer(int64)              :: longseed
    longseed = int(seed*48271, int64)
    seed     = int(iand(longseed, maxint32), int32)
    rnd      = (seed*to_real64+1)/2 ! [0, 1)
  end subroutine lcg

end module builtin_rand_m
