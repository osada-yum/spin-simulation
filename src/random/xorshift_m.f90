module xorshift_m
  use, intrinsic :: iso_fortran_env
  use random_base_m
  implicit none
  private
  real(real64), parameter :: to_real64 = 1 / (real(huge(1_int32), real64) + 1)

  public :: xor96
  !> `xor96`: From [https://ja.wikipedia.org/wiki/Xorshift].
  type, extends(random_base_t) :: xor96
     private
     integer :: x, y, z
   contains
     procedure, pass :: set_seed   => set_seed_xor96
     procedure, pass :: get_rand   => get_rand_xor96
     procedure, pass :: random_arr => random_arr_xor96
  end type xor96

contains

  !> `set_seed_xor96`: Set seed of xorshift.
  subroutine set_seed_xor96(this, iseed)
    class(xor96), intent(inout) :: this
    integer     , intent(in)    :: iseed
    this%x = 123456789
    this%y = 362436069
    this%z = ieor(521288629, iseed)
  end subroutine set_seed_xor96

  !> `get_rand_xor96`: Return random number of xorshift.
  real(real64) function get_rand_xor96(this)
    class(xor96), intent(inout) :: this
    integer                     :: tx, ty, tz
    tx = ieor(this%x, ishft(this%x, 3))
    ty = ieor(this%y, ishft(this%y, -19))
    tz = ieor(this%z, ishft(this%z, 6))
    this%x = this%y
    this%y = this%z
    this%z = ieor(ieor(tx, ty), tz)
    get_rand_xor96 = (this%z*to_real64 + 1.0_real64)/2.0_real64 ! [0, 1)
    return
  end function get_rand_xor96

  !> `random_arr_xor96`: Return random numbers of xorshift by argument `arr`.
  subroutine random_arr_xor96(this, arr)
    class(xor96), intent(inout) :: this
    real(real64), intent(out)   :: arr(:)
    integer                     :: i
    do i = 1, size(arr)
       arr(i) = this%get_rand()
    end do
  end subroutine random_arr_xor96

end module xorshift_m
