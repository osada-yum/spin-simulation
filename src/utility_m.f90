module utility_m
  use, intrinsic :: iso_fortran_env
  implicit none
  integer    , parameter :: ikind = int32, rkind = real64
  real(rkind), parameter :: epsilon = 1e-7_rkind
contains

  subroutine util_error_stop(message, linum, filename)
    integer, intent(in) :: linum
    character(len=*)    :: message, filename
    write(error_unit, '(a, i0)') "error in "//filename//":", linum
    write(error_unit, '(a)'    ) message
#ifdef DEBUG
    error stop 1
#endif
  end subroutine util_error_stop

end module utility_m
