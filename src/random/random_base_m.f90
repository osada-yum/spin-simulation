module random_base_m
  use, intrinsic :: iso_fortran_env
  implicit none
  private

  public :: random_base_t
  !> `random_base_t`: Provide interfaces of generating random number.
  type, abstract :: random_base_t
   contains
     procedure(set_seed_i), deferred :: set_seed
     procedure(get_rand_i), deferred :: get_rand
     procedure(arr_rand_i), deferred :: random_arr
  end type random_base_t

  interface
     !> `set_seed_i`: Initialize seed.
     subroutine set_seed_i(this, iseed)
       import random_base_t
       class(random_base_t), intent(inout) :: this
       integer             , intent(in)    :: iseed
     end subroutine set_seed_i
     !> `get_rand_i`: Return one random number.
     real(real64) function get_rand_i(this)
       import random_base_t, real64
       class(random_base_t), intent(inout) :: this
     end function get_rand_i
     !> `arr_rand_i`: Set size(arr) random numbers in arr(:)
     subroutine arr_rand_i(this, arr)
       import random_base_t, real64
       class(random_base_t), intent(inout) :: this
       real(real64)        , intent(out)   :: arr(:)
     end subroutine arr_rand_i
  end interface
end module random_base_m
