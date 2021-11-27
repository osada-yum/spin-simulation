program test
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: iso_c_binding
    implicit none
    integer(int32) :: a
    real(real64), parameter :: rinvm = 1/ (real(huge(1_int32), real64)+1)
    print *, -2147483648, huge(0_int32), rinvm
    print *, -2147483648*rinvm, huge(0_int32)*rinvm
    print *, -2147483647, huge(0_int32)-1, rinvm
    print *, -2147483647*rinvm, (huge(0_int32)-1)*rinvm
    print *, (-2147483648*rinvm+1.0_real64)/2.0_real64&
        , (huge(0_int32)*rinvm+1.0_real64)/2.0_real64
end program
