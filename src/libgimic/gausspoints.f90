!
! GIMIC Fortran to C interface
!
! The routinese are to be used from C/C++ to access GIMIC functionality
!
module gausspoints
    use iso_c_binding
    use gaussint_module
    implicit none

contains
    subroutine gimic_get_gauss_points(a, b, npts, order, pts, wgts) bind(c)
        real(C_DOUBLE), intent(in) :: a, b
        integer(C_INT), intent(in) :: npts, order
        real(C_DOUBLE), dimension(1), intent(out) :: pts, wgts

        integer :: i
        type(gdata_t) :: g

        allocate(g%pts(npts))
        allocate(g%wgt(npts))
        call setup_gauss_data(a, b, order, g, 'gauss')
        do i=1, npts
            pts(i) = g%pts(i)
            wgts(i) = g%wgt(i)
        end do
        deallocate(g%pts)
        deallocate(g%wgt)
    end subroutine
end module

