module utility_m
  use, intrinsic :: iso_fortran_env, only : int32, real64, error_unit
  implicit none
  integer    , parameter :: ikind = int32, rkind = real64
  real(rkind), parameter :: epsilon = 1e-7_rkind
contains

  !> `util_linspace`: calculate Internal division point.
  pure function util_linspace(from, to, num) result(arr)
    real(rkind), intent(in) :: from, to
    integer    , intent(in) :: num
    real(rkind)             :: arr(num)
    integer                 :: i
    do i = 1, num
       arr(i) = ( (num-1-i+1)*from + (i-1)*to ) / (num-1)
    end do
  end function util_linspace

  !> `util_debug_print`: print mesage if compiler flag is -DDEBUG.
  subroutine util_debug_print(message, linum, filename)
    integer, intent(in) :: linum
    character(len=*)    :: message, filename
#ifdef DEBUG
    write(error_unit, '(a, i0)') "debug message in "//filename//":", linum
    write(error_unit, '(a)'    ) message
#endif
  end subroutine util_debug_print

  !> `util_warning`: print warning. error stop if -DDEBUG.
  subroutine util_warning(message, linum, filename)
    integer, intent(in) :: linum
    character(len=*)    :: message, filename
    write(error_unit, '(a, i0)') "warning in "//filename//":", linum
    write(error_unit, '(a)'    ) message
#ifdef DEBUG
    error stop 1
#endif
  end subroutine util_warning

  !> `util_error_stop`: print error. error stop.
  subroutine util_error_stop(message, linum, filename)
    integer, intent(in) :: linum
    character(len=*)    :: message, filename
    write(error_unit, '(a, i0)') "error in "//filename//":", linum
    write(error_unit, '(a)'    ) message
    error stop 2
  end subroutine util_error_stop

end module utility_m
