module libfortinterpolate

use iso_c_binding
use ISO_FORTRAN_ENV

integer, parameter :: dp = REAL64

contains

subroutine foo(Nx,Ny,input,output) bind(c,name='foo')

    integer(c_int), intent(in)  :: Nx,Ny
    real(c_double), intent(in)  :: input(Nx,Ny)
    real(c_double), intent(out) :: output(Nx,Ny)

    output = 42.0

end subroutine

end module
