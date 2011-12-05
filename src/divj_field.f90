!
!  Calculate the divergence of a vector field.
!  Useful for checking the convergce of the actual gauge invariance,
!  in contrast to the basis set convergence and gauge independence.
!

module divj_field_class
    use globals_m
    use grid_class
    use divj_m
    use jtensor_class
    use jfield_class
    use teletype_m
    use parallel_m
    implicit none

    type divj_field_t
        real(DP), dimension(:,:), pointer :: buf
        type(grid_t), pointer :: grid
        real(DP), dimension(3) :: bb
    end type

    public divj_direct, divj_direct_plt, new_divj_field, del_divj_field
    public divj_plot, set_divj, divj_field_t, divj_field
    private
    
    real(DP), dimension(5), parameter :: wgt=(/2.d0,-16.d0,0.d0,16.d0,-2.d0/)
    real(DP), parameter :: step=1.d-3
    real(DP), parameter :: hx=1.d0/(24.d0*step)
    real(DP), parameter :: hy=1.d0/(24.d0*step)
    real(DP), parameter :: hz=1.d0/(24.d0*step)

contains
    subroutine new_divj_field(this, grid, magnet)
        type(divj_field_t), intent(inout) :: this
        type(grid_t), target :: grid
        real(DP), dimension(3), intent(in) :: magnet

        integer(I4) :: djrl, p1, p2
        logical :: foo_p

        this%bb=magnet
        p1=0

        call get_grid_size(grid, p1, p2)
        djrl=p1*p2*DP
        if (master_p) then
            open(DIVJFD, file='DIVJ', access='direct', recl=djrl)
        end if

        this%grid=>grid
        allocate(this%buf(p1,p2))
    end subroutine

    subroutine del_divj_field(this)
        type(divj_field_t), intent(inout) :: this

        deallocate(this%buf)
        if (master_p) then
            close(DIVJFD)
        end if
    end subroutine

    subroutine set_divj(this, k)
        type(divj_field_t), intent(in) :: this
        integer(I4), intent(in) :: k

        write(DIVJFD, rec=k) this%buf
    end subroutine

    subroutine divj_plot(this, gopen_file)
        type(divj_field_t), intent(inout) :: this
        character(*), intent(in) :: gopen_file
        
        integer(I4) :: i,j,p1,p2,p3
        real(DP) :: amax
        real(DP), dimension(3) :: rr
        real(DP), dimension(:,:), pointer :: buf

        call get_grid_size(this%grid, p1, p2, p3)

        buf=>this%buf

        amax=D0
        read(DIVJFD, rec=1) this%buf
        str_g='DIVJPLT'
        open(DJPFD, file=trim(str_g))
        do j=1,p2
            do i=1,p1
                rr=gridpoint(this%grid, i, j, 1)
                write(DJPFD, '(4f19.8)') rr, buf(i,j)
                if (abs(buf(i,j)) > amax) amax=abs(buf(i,j))
            end do
            write(DJPFD, *) 
        end do
        close(DJPFD)
        write(str_g, '(a,e19.12)') 'Max divergence:', amax
        call msg_info(str_g)
        call divj_gopenmol(this, gopen_file)
    end subroutine

    subroutine divj_direct_plt(this)
        type(divj_field_t), intent(inout) :: this

        type(jtensor_t) :: jt
        !type(jfield_t) :: jf
        integer(I4) :: i, j, k, p1, p2, p3
        real(DP) :: div, amax
        real(DP), dimension(3) :: rr 
        
        if (mpirun_p) then
            call msg_error('divj_direct_plt(): does not work in parallel')
            stop
        end if
        !jf%jt=>jt
        !jf%grid=>this%grid
        !call jfield_eta(jf, 12.d0)
        call get_grid_size(this%grid, p1, p2, p3)
    
        open(DJPFD, file='DIVJPLT')
        amax=D0
        do k=1,p3
            do j=1,p2
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    div = divj(rr, this%bb)
                    write(DJPFD, '(4f19.8)') rr, div
                    if (abs(div) > amax) amax=abs(div)
                end do
                write(DJPFD, *) 
            end do
        end do
        write(str_g, '(a,e19.12)')'Max divergence:', amax
        call msg_note(str_g)
        call nl
        close(DJPFD)
    end subroutine

    subroutine divj_direct(this, k)
        type(divj_field_t) :: this
        integer(I4), intent(in) :: k

        integer(I4) :: i, j, p1, p2
        integer(I4) :: lo, hi
        real(DP), dimension(3) :: rr
        real(DP), dimension(:,:), pointer :: buf
        real(DP) :: div

        call get_grid_size(this%grid, p1, p2)
        call schedule(p2, lo, hi)

        buf=>this%buf

        do j=1,p2
            do i=1,p1
                rr=gridpoint(this%grid, i, j, k)
                div = divj(rr, this%bb)
            end do
        end do
        call gather_data(buf,buf(:,lo:hi))
    end subroutine

    subroutine divj_field(this)
        type(divj_field_t) :: this

        integer(I4) :: i, j, k, p, p1, p2, p3
        integer(I4) :: lo, hi
        real(DP), dimension(3) :: rr
        real(DP), dimension(:,:), pointer :: buf

        call get_grid_size(this%grid, p1, p2, p3)

        call schedule(p2, lo, hi)

        buf=>this%buf
    
        do k=1,p3
            do j=lo, hi
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    buf(i,j) = divj(rr, this%bb)
                end do
            end do
            call gather_data(buf,buf(:,lo:hi))
            if (master_p) write(DIVJFD, rec=k) this%buf
        end do
    end subroutine

    subroutine divj_gopenmol(this, gopen_file)
        type(divj_field_t) :: this
        character(*), intent(in) :: gopen_file

        integer(I4) :: surface, rank, p1, p2, p3
        integer(I4) :: i, j, k, l
        real(SP), dimension(3) :: qmin, qmax
        real(DP), dimension(:,:), pointer :: buf

        buf=>this%buf
        if (trim(gopen_file) == '') return
        open(GOPFD,file=trim(gopen_file),access='direct',recl=4)

        surface=200
        rank=3

        call get_grid_size(this%grid, p1, p2, p3)
        qmin=real(gridpoint(this%grid,1,1,1)*AU2A)
        qmax=real(gridpoint(this%grid,p1,p2,p3)*AU2A)

        write(GOPFD,rec=1) rank
        write(GOPFD,rec=2) surface
        write(GOPFD,rec=3) p3
        write(GOPFD,rec=4) p2
        write(GOPFD,rec=5) p1
        write(GOPFD,rec=6) qmin(3)
        write(GOPFD,rec=7) qmax(3)
        write(GOPFD,rec=8) qmin(2)
        write(GOPFD,rec=9) qmax(2)
        write(GOPFD,rec=10) qmin(1)
        write(GOPFD,rec=11) qmax(1)

!        write(100,*) rank
!        write(100,*) surface
!        write(100,*) p3
!        write(100,*) p2
!        write(100,*) p1
!        write(100,*) qmin(3)
!        write(100,*) qmax(3)
!        write(100,*) qmin(2)
!        write(100,*) qmax(2)
!        write(100,*) qmin(1)
!        write(100,*) qmax(1)

        l=12
        do k=1,p3
            read(DIVJFD, rec=k) buf
            do j=1,p2
                do i=1,p1
                    write(GOPFD,rec=l) real(buf(i,j))
!                    write(100,*) real(buf(i,j))
                    l=l+1
                end do
            end do
        end do

        close(GOPFD)
    end subroutine
end module

! vim:et:sw=4:ts=4