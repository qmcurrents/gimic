!
! Calulate the electronic density on a grid.
! Coded by Jonas Juselius <jonas@iki.fi> 2003 & 2011
!

module edens_field_class
    use globals_module
    use settings_module
    use basis_class
    use bfeval_class
    use dens_class
    use grid_class
    use edens_class
    use teletype_module
    use parallel_module
    use cubeplot_module
    implicit none

    public new_edens_field, del_edens_field, edens_field_t
    public edens_plot, edens_field

    type edens_field_t
        type(edens_t) :: edens
        type(grid_t), pointer :: grid
        real(DP), dimension(:,:,:), pointer :: buf
    end type

    private
contains
    ! set up memory (once) for the different components
    subroutine new_edens_field(this, mol, dens, grid, densfile)
        type(edens_field_t) :: this
        type(molecule_t) :: mol
        type(dens_t), target :: dens
        type(grid_t), target :: grid
        character(*), intent(in) :: densfile

        integer(I4) :: n
        integer(I4) :: p1, p2, p3

        call new_edens(this%edens, mol, dens)
        call get_grid_size(grid, p1, p2, p3)
        allocate(this%buf(p1,p2,p3))
        this%grid=>grid
    end subroutine

    subroutine del_edens_field(this)
        type(edens_field_t) :: this

        if (associated(this%buf)) deallocate(this%buf)
        call del_edens(this%edens)
        nullify(this%grid)
    end subroutine

    subroutine edens_plot(this, basename)
        type(edens_field_t), intent(inout) :: this
        character(*), intent(in) :: basename

        integer(I4) :: i,j,p1,p2,p3
        real(DP) :: amax
        real(DP), dimension(3) :: rr
        real(DP), dimension(:,:,:), pointer :: buf
        character(BUFLEN) :: fname

        fname = trim(basename) // '.txt'

        call get_grid_size(this%grid, p1, p2,p3)
        buf=>this%buf

        amax=D0
        open(EDPFD, file=trim(fname))
        do j=1,p2
            do i=1,p1
                rr=gridpoint(this%grid, i, j, 1)
                write(EDPFD, '(4f19.8)') rr, buf(i,j,1)
                if (abs(buf(i,j,1)) > amax) amax=abs(buf(i,j,1))
            end do
            write(EDPFD, *)
        end do
        close(EDPFD)
        write(str_g, '(a,e19.12)') 'Max electronic density:', amax
        call msg_info(str_g)
        call edens_cube(this, basename)
    end subroutine

    subroutine edens_field(this)
        type(edens_field_t) :: this

        integer(I4) :: i, j, k, p1, p2, p3
        real(DP), dimension(3) :: rr
        real(DP), dimension(:,:,:), pointer :: buf
        real(DP), dimension(:), pointer :: bfvec

        call get_grid_size(this%grid, p1, p2, p3)

        buf=>this%buf

        do k=1,p3
            do j=1,p2
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    this%buf(i,j,k)=edens(this%edens, rr)
                end do
            end do
        end do
    end subroutine

    subroutine edens_cube(this, basename)
        type(edens_field_t) :: this
        character(*), intent(in) :: basename

        integer(I4) :: p1, p2, p3
        integer(I4) ::  k
        real(DP), dimension(:,:,:), allocatable :: buf
        character(BUFLEN) :: fname

        if (trim(basename) == '') return
        fname = trim(basename) // '.cube'

        call get_grid_size(this%grid, p1, p2, p3)
        if (p3 < 2) return

        allocate(buf(p1,p2,p3))
        do k=1,p3
            read(EDFD, rec=k) buf(:,:,k)
        end do
        call write_cubeplot(trim(fname), this%grid, buf)
        deallocate(buf)

    end subroutine
end module

! vim:et:sw=4:ts=4
