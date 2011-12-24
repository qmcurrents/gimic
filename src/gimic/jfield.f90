! This module calulates the current field by looping over gridpoints and 
! calling the appropriate routines in jtensor_class
!
! Ok, this is really a mess...
! This modules needs a seriuos clean up. I have been adding, and adding, and
! adding stuff and it's not pretty.
!

module jfield_class
    use globals_m
    use settings_m
    use jtensor_class
    use dens_class
    use grid_class
    use basis_class
    use teletype_m
    use parallel_m
    use tensor_m
    implicit none

    type jfield_t
        !type(jtensor_t), pointer :: jt
        real(DP), dimension(3) :: b
        type(tensor_t), dimension(:,:,:), pointer :: jj
        type(vector_t), dimension(:,:,:), pointer :: vv
        type(grid_t), pointer :: grid
    end type

    public new_jfield, del_jfield, jfield, jvectors
    public jfield_t, jvector_plot, jfield_eta
    private

    real(DP) :: vec_scale=D1
    character(BUFLEN) :: jmod_plt, jvec_plt
contains
    subroutine new_jfield(this, g, magnet)
        type(jfield_t) :: this
!        type(jtensor_t), target :: jt
        type(grid_t), target :: g
        real(DP), dimension(3), intent(in) :: magnet

        integer(I4) :: i, j, k

        this%b=magnet
!        this%jt=>jt
        this%grid=>g

        call get_grid_size(g,i,j,k)
        allocate(this%jj(i,j,k))
        allocate(this%vv(i,j,k))
    end subroutine

    subroutine del_jfield(this)
        type(jfield_t), intent(inout) :: this
        
        this%b=0.d0
        deallocate(this%jj)
        deallocate(this%vv)
        nullify(this%jj)
        nullify(this%vv)
    end subroutine

    subroutine jfield(this, mol, xdens, spin, z)
        type(jfield_t) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens
        character(*), optional :: spin
        integer, optional :: z

        integer(I4) :: i, j, k, p1, p2, p3
        integer :: a, b
        real(DP), dimension(3) :: rr
        character(8) :: spincase
        integer(I4) :: lo, hi, npts
        type(jtensor_t) :: jt

        if (present(spin)) then
            spincase=spin
        else
            spincase='total'
        end if
        call jfield_eta(this, mol, xdens)
        call get_grid_size(this%grid, p1, p2, p3)
        if (present(z)) then
            a=z
            b=z
        else
            a=1
            b=p3
        end if

!$OMP PARALLEL PRIVATE(i,j,k,rr,jt) &
!$OMP SHARED(p1,p2,p3,a,b,settings,this,mol,xdens)
        call new_jtensor(jt, mol, xdens)
        do k=a,b
            !$OMP DO SCHEDULE(STATIC) 
            do j=1,p2
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    call ctensor(jt, rr, this%jj(i,j,k), spincase)
                end do
            end do
            !$OMP END DO
        end do
        call del_jtensor(jt)
!$OMP END PARALLEL
    end subroutine 
    
    subroutine jvectors(this, mol, xdens, spin, z)
        type(jfield_t) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens
        character(*), optional :: spin
        integer, optional :: z

        character(8) :: spincase
        integer(I4) :: i, j, k, p1, p2, p3, a, b

        if (present(spin)) then
            spincase=spin
        else
            spincase='total'
        end if

        call get_grid_size(this%grid, p1, p2, p3)
        if (present(z)) then
            a=z
            b=z
        else
            a=1
            b=p3
        end if
        call msg_note('Contracting j-tensors with magnetic field.')
        call nl
        do k=a,b
            call jfield(this, mol, xdens, spincase, k)
            do j=1,p2
                do i=1,p1
                    this%vv(i,j,k)%v=matmul(this%jj(i,j,k)%t, this%b)
                end do
            end do
!            call plot()
        end do
    end subroutine

    subroutine jfield_eta(this, mol, xdens, fac)
        type(jfield_t) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens
        real(DP), intent(in), optional :: fac

        integer(I4) :: i, j, p1, p2, p3
        type(tensor_t) :: foo
        type(jtensor_t) :: jt
        real(DP) :: delta_t
        real(4) :: tim1, tim2
        real(4), dimension(2) :: times
        real(DP), dimension(3) :: bar=(/D1,D1,D1/)
        real(DP), dimension(3) :: foobar
        real(DP), parameter :: SC=0.25d0
        
        call get_grid_size(this%grid, p1, p2, p3)
        
        call new_jtensor(jt, mol, xdens)
        call etime(times, tim1)
        tim1=times(1)
        do i=1,100
            call jtensor(jt, (/i*SC, i*SC, i*SC/), foo, spin_a)
            foobar=matmul(bar,foo%t)
        end do
        call etime(times, tim2)
        tim2=times(1)
        call del_jtensor(jt)
        
        delta_t=tim2-tim1
        if ( present(fac) ) delta_t=delta_t*fac
        write(str_g, '(a,f9.2,a)') 'Estimated CPU time for &
            &calculation: ', delta_t*real(p1*p2*p3)/100.d0, ' sec'
        call msg_info(str_g)
        call nl
    end subroutine	

    function open_plot(basename, spin) result(fd)
        character(*), intent(in) :: basename
        integer(I4), intent(in) :: spin

        integer(I4) :: fd
        character(132) :: fname

        fd=0
        if (mpi_rank > 0) return
        if (trim(basename) == '') return

        call getfd(fd)
        if (fd == 0) then
            stop 1
        end if

        select case(spin)
            case(1)
                fname=trim(basename)//'.txt'
            case(2)
                fname=trim(basename)//'_a.txt'
            case(3)
                fname=trim(basename)//'_b.txt'
            case(4)
                fname=trim(basename)//'_sd.txt'
        end select

        write(str_g, *) 'Writing ', trim(basename), ' in ', trim(fname)
        call msg_note(str_g)
        open(fd,file=trim(fname),status='unknown')

        return 
    end function

    subroutine jvector_plot(this)
        type(jfield_t), intent(inout) :: this

        integer(I4) :: i, j, k, p1, p2, p3, spin, ispin
        integer(I4) :: fd1, fd2, fd3, fd4
        real(DP), dimension(3) :: v, rr
        type(vector_t), dimension(:,:,:), pointer :: jv

        spin=1
        if (settings%is_uhf) then
            spin=4
        end if

        call get_grid_size(this%grid, p1, p2, p3)

        do ispin=1,spin
            fd1= open_plot('jvec',ispin)
            fd2= open_plot('jmod',ispin)

            jv=>this%vv
            do k=1,p3
                do j=1,p2
                    do i=1,p1
                        rr=gridpoint(this%grid, i,j,k)*AU2A
                        v=jv(i,j,k)%v*AU2A
                        call wrt_jvec(rr,v,fd1)
                        call wrt_jmod(rr,v,fd2)
                    end do
                    if (fd1 /= 0) write(fd1, *)
                    if (fd2 /= 0) write(fd2, *)
                end do
            end do

            call closefd(fd1)
            call closefd(fd2)
            if (p3 > 1) call jmod_cubeplot(this, 'jmod', ispin)
        end do
    end subroutine

    subroutine write_jvectors(this, fd)
        type(jfield_t), intent(inout) :: this
        integer(I4) :: fd

        integer(I4) :: i, j, k, p1, p2, spin, ispin
        real(DP), dimension(3) :: v, rr
        type(vector_t), dimension(:,:,:), pointer :: jv

        call get_grid_size(this%grid, p1, p2)

        k=1
        jv=>this%vv
        do j=1,p2
            do i=1,p1
                rr=gridpoint(this%grid, i,j,k)*AU2A
                v=jv(i,j,k)%v*AU2A
                call wrt_jvec(rr,v,fd)
            end do
            if (fd /= 0) write(fd, *)
        end do
    end subroutine

    subroutine wrt_jvec(rr,v,fd)
        real(DP), dimension(:), intent(in) :: rr, v
        integer(I4), intent(in) :: fd

        if (fd == 0) return

        write(fd, '(6f11.7)')  rr, v
    end subroutine	

    subroutine wrt_njvec(rr,v,fd)
        real(DP), dimension(:), intent(in) :: rr, v
        integer(I4), intent(in) :: fd

        real(DP) :: nfac

        if (fd == 0) return

        nfac=sqrt(sum(v(:)**2))
        if (nfac < 1.d-15) then
            write(fd, '(6f11.7)')  rr, 0.d0, 0.d0, 0.d0
        else 
            write(fd, '(6f11.7)')  rr, v/nfac
        end if
    end subroutine	

    subroutine wrt_jproj(rr,v,grid,fd)
        real(DP), dimension(:), intent(in) :: rr, v
        type(grid_t) :: grid
        integer(I4), intent(in) :: fd

        real(DP) :: jprj
        real(DP), dimension(3) :: norm

        if (fd == 0) return
        norm=get_grid_normal(grid)

        jprj=dot_product(norm,v)
        write(fd, '(3e19.12)') rr, jprj
    end subroutine	

    subroutine wrt_jmod(rr,v,fd)
        real(DP), dimension(:), intent(in) :: rr, v
        integer(I4), intent(in) :: fd

        real(DP) :: jmod

        if (fd == 0) return

        jmod=sqrt(sum(v**2))
        write(fd, '(6f11.7)')  rr, jmod
    end subroutine	

    subroutine print_jt(rr, jt)
        real(DP), dimension(3), intent(in) :: rr
        real(DP), dimension(3,3), intent(in) :: jt

        integer(I4) :: l
        print *
        print '(a,3f12.8)', 'Current tensor at ', rr
        print *, '======================================='
        print '(3e19.12)', (jt(l,:), l=1,3)
        print *
        print '(a,e19.12)', ' Trace:', jt(1,1)+jt(2,2)+jt(3,3)
    end subroutine

    subroutine jmod_cubeplot(this, fname, spin)
        type(jfield_t) :: this
        character(*), intent(in) :: fname
        integer(I4), intent(in) :: spin

        integer(I4) :: p1, p2, p3, fd1, fd2
        integer(I4) :: i, j, k, l
        real(DP), dimension(3) :: qmin, qmax
        real(DP), dimension(3) :: norm, step, mag, v, rr
        real(DP) :: maxi, mini, val, sgn 
        integer(I4), dimension(3) :: npts
        type(vector_t), dimension(:,:,:), pointer :: buf

        call get_grid_size(this%grid, p1, p2, p3)
        npts=(/p1,p2,p3/)
        norm=get_grid_normal(this%grid)
        qmin=gridpoint(this%grid,1,1,1)
        qmax=gridpoint(this%grid,p1,p2,p3)

        step=(qmax-qmin)/(npts-1)

        fd1=opencube(trim(fname), spin, qmin, step, npts)
        fd2=opencube(trim(fname)//'_quasi', spin, qmin, step, npts)

        mag = this%b

        buf => this%vv
        maxi=0.d0
        mini=0.d0
        l=0
        do i=1,p1
            do j=1,p2
                do k=1,p3
                    v=buf(i,j,k)%v
                    rr=gridpoint(this%grid,i,j,k)
                    val=(sqrt(sum(v**2)))
                    rr=rr-dot_product(mag,rr)*mag
                    norm=cross_product(mag,rr)
                    sgn=dot_product(norm,v)
                    if (val > maxi) maxi=val
                    if (val < mini) mini=val
                    if (fd1 /= 0) then 
                        write(fd1,'(f12.6)',advance='no') val
                        if (mod(l,6) == 5) write(fd1,*)
                    end if
                    if (fd2 /= 0) then 
                        if (sgn >= 0.d0 ) then
                            write(fd2,'(f12.6)',advance='no') val
                        else
                            write(fd2,'(f12.6)',advance='no') -1.d0*val
                        end if
                        if (mod(l,6) == 5) write(fd2,*)
                    end if
                    l=l+1
                end do
            end do
        end do
        print *, 'maximini:', maxi, mini
    end subroutine

    function opencube(kname, spin, origin, step, npts) result(fd)
        character(*), intent(in) :: kname
        integer(I4), intent(in) :: spin
        real(DP), dimension(3), intent(in) :: origin, step
        integer(I4), dimension(3), intent(in) :: npts

        integer(I4) :: fd
        character(132) :: fname

        fd=0
        fname=kname
        if (mpi_rank > 0) return
        if (trim(fname) == '') return

        call getfd(fd)
        if (fd == 0) then
            stop 1
        end if

        select case(spin)
            case(1)
                fname=trim(fname)//'.cube'
            case(2)
                fname=trim(fname)//'a'//'.cube'
            case(3)
                fname=trim(fname)//'b'//'.cube'
            case(4)
                fname=trim(fname)//'sd'//'.cube'
        end select

        write(str_g, *) 'Writing ', trim(kname), ' in ', trim(fname)
        call msg_note(str_g)

        open(fd,file=trim(fname),form='formatted',status='unknown')
        write(fd,*) 'Gaussian cube data, generated by genpot'
        write(fd,*) 
        write(fd, '(i5,3f12.6)') 0, origin
        write(fd, '(i5,3f12.6)') npts(1), step(1), 0.d0, 0.d0
        write(fd, '(i5,3f12.6)') npts(2), 0.d0, step(2), 0.d0
        write(fd, '(i5,3f12.6)') npts(3), 0.d0, 0.d0, step(3)

        return 

    end function

end module

! vim:et:sw=4:ts=4
