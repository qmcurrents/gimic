!
! Calulate the electronic density on a grid.
! Coded by Jonas Juselius <jonas@iki.fi> 2003
! 

module edens_class
    use globals_m
    use basis_class
    use bfeval_class
    use dens_class
    use grid_class
    use teletype_m
    use parallel_m
    use cubeplot_m
    implicit none
    
    public new_edens, del_edens, edens_t
    public edens_direct, edens, edens_plot
    
    type edens_t
        real(DP), dimension(:), pointer :: tmp
        real(DP), dimension(:,:), pointer :: aodens
        real(DP), dimension(:,:), pointer :: buf
        type(bfeval_t) :: bf
        type(grid_t), pointer :: grid
    end type

    private

contains
    ! set up memory (once) for the different components
    subroutine new_edens(this, mol, dens, grid, densfile)
        type(edens_t) :: this
        type(molecule_t) :: mol
        type(dens_t), target :: dens
        type(grid_t), target :: grid
        character(*), intent(in) :: densfile

        integer(I4) :: n
        integer(I4) :: p1, p2

        call get_grid_size(grid, p1, p2)
        call get_dens(dens, this%aodens)
        n=get_ncgto(mol)
        allocate(this%tmp(n))
        allocate(this%buf(p1,p2))
        call new_bfeval(this%bf, mol)
        this%grid=>grid

        if (master_p) then
            open(EDFD, file=trim(densfile), access='direct', recl=p1*p2*DP)
        end if
    end subroutine

    subroutine del_edens(this)
        type(edens_t) :: this

        if (associated(this%tmp)) deallocate(this%tmp)
        if (associated(this%buf)) deallocate(this%buf)
        call del_bfeval(this%bf)
        nullify(this%aodens,this%grid)
        if (master_p) then
            close(EDFD)
        end if
    end subroutine

    subroutine set_edens(this, k)
        type(edens_t), intent(in) :: this
        integer(I4), intent(in) :: k

        write(EDFD, rec=k) this%buf
    end subroutine

    subroutine edens_plot(this, basename)
        type(edens_t), intent(inout) :: this
        character(*), intent(in) :: basename
        
        integer(I4) :: i,j,p1,p2,p3
        real(DP) :: amax
        real(DP), dimension(3) :: rr
        real(DP), dimension(:,:), pointer :: buf
        character(BUFLEN) :: fname

        fname = trim(basename) // '.txt'

        call get_grid_size(this%grid, p1, p2,p3)
        buf=>this%buf

        amax=D0
        read(EDFD, rec=1) this%buf
        open(EDPFD, file=trim(fname))
        do j=1,p2
            do i=1,p1
                rr=gridpoint(this%grid, i, j, 1)
                write(EDPFD, '(4f19.8)') rr, buf(i,j)
                if (abs(buf(i,j)) > amax) amax=abs(buf(i,j))
            end do
            write(EDPFD, *) 
        end do
        close(EDPFD)
        write(str_g, '(a,e19.12)') 'Max electronic density:', amax
        call msg_info(str_g)
        call edens_cube(this, basename)
    end subroutine

    subroutine edens_direct(this, k)
        type(edens_t) :: this
        integer(I4), intent(in) :: k

        integer(I4) :: i, j, p1, p2
        integer(I4) :: lo, hi
        real(DP), dimension(3) :: rr
        real(DP), dimension(:,:), pointer :: buf
        real(DP), dimension(:), pointer :: bfvec

        call get_grid_size(this%grid, p1, p2)
        call schedule(p2, lo, hi)

        buf=>this%buf

        do j=lo,hi
            do i=1,p1
                rr=gridpoint(this%grid, i, j, k)
                call bfeval(this%bf,rr, bfvec)
                this%tmp=matmul(this%aodens, bfvec)
                this%buf(i,j)=dot_product(this%tmp, bfvec)
            end do
        end do
        call gather_data(this%buf, this%buf(:,lo:hi))
    end subroutine

    subroutine edens(this)
        type(edens_t) :: this

        integer(I4) :: i, j, k, p1, p2, p3
        integer(I4) :: lo, hi
        real(DP), dimension(3) :: rr
        real(DP), dimension(:,:), pointer :: buf
        real(DP), dimension(:), pointer :: bfvec

        call get_grid_size(this%grid, p1, p2, p3)
        call schedule(p2, lo, hi)

        buf=>this%buf
    
        do k=1,p3
            do j=lo, hi
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    call bfeval(this%bf,rr, bfvec)
                    this%tmp=matmul(this%aodens, bfvec)
                    this%buf(i,j)=dot_product(this%tmp, bfvec)
                end do
            end do
            call gather_data(this%buf, this%buf(:,lo:hi))
            if (master_p) write(EDFD, rec=k) this%buf
        end do
    end subroutine

    subroutine edens_cube(this, basename)
        type(edens_t) :: this
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
