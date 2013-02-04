! This module calulates the current field by looping over gridpoints and 
! calling the appropriate routines in jtensor_class
!
! Ok, this is really a mess...
! This modules needs a seriuos clean up. I have been adding, and adding, and
! adding stuff and it's not pretty.
!

module jfield_class
    use globals_module
    use settings_module
    use jtensor_class
    use dens_class
    use grid_class
    use basis_class
    use teletype_module
    use parallel_module
    use tensor_module
    ! ACID stuff
    use acid_module
    !
    implicit none

    type jfield_t
        real(DP), dimension(3) :: b
        real(DP), dimension(:,:), pointer :: tens
        real(DP), dimension(:,:), pointer :: vec
        type(grid_t), pointer :: grid
    end type

    public new_jfield, del_jfield, calc_jtensors, calc_jvectors
    public jfield_t, jvector_plots, jfield_eta, compute_jvectors
    public jmod_cubeplot
    private

    real(DP) :: vec_scale=D1
    character(BUFLEN) :: jmod_plt, jvec_plt
contains
    subroutine new_jfield(this, g, magnet)
        type(jfield_t) :: this
        type(grid_t), target :: g
        real(DP), dimension(3), intent(in) :: magnet

        integer(I4) :: i, j, k

        this%b=magnet
        this%grid=>g

        call get_grid_size(g,i,j,k)
        allocate(this%tens(9,i*j*k))
        allocate(this%vec(3,i*j*k))
    end subroutine

    subroutine del_jfield(this)
        type(jfield_t), intent(inout) :: this
        
        this%b=0.d0
        deallocate(this%tens)
        deallocate(this%vec)
    end subroutine

    subroutine calc_jtensors(this, mol, xdens, spin, z)
        type(jfield_t) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens
        character(*), optional :: spin
        integer, optional :: z

        integer(I4) :: i, j, k, p1, p2, p3
        integer :: n, m
        real(DP), dimension(3) :: rr
        real(DP), dimension(:,:), pointer :: tens
        character(8) :: spincase
        integer(I4) :: lo, hi, npts, first, last, fd2
        type(jtensor_t) :: jt

        call get_grid_size(this%grid, p1, p2, p3)
        if (present(z)) then
            n=z
            m=z
        else
            n=1
            m=p3
        end if

        lo = (n-1)*p1*p2 + 1
        hi = m*p1*p2
        npts = hi - lo + 1

        call schedule(npts, lo, hi)
        npts = hi - lo + 1

        first = 1
        last = npts
        ! rank 0 handles any leftover nodes, located at the beginning
        if (mpi_world_size > 1 .and. mpi_rank == 0) then
            first=lo
            last = hi
            lo = 1
            npts = hi - lo + 1
            allocate(tens(9,npts))
        else
            tens=>this%tens
        end if

        if (present(spin)) then
            spincase=spin
        else
            spincase='total'
        end if

        call jfield_eta(this, mol, xdens)
        !$OMP PARALLEL PRIVATE(jt,rr,n,i,j,k) &
        !$OMP SHARED(mol,xdens,spincase,tens,lo,hi)
        call new_jtensor(jt, mol, xdens)
        !$OMP DO SCHEDULE(STATIC) 
        do n=lo,hi
            call get_grid_index(this%grid, n, i, j, k)
            rr = gridpoint(this%grid, i, j, k)
            ! here the ACID T tensor is calculated and put on tens
            call ctensor(jt, rr, tens(:,n-lo+1), spincase)
        end do
        !$OMP END DO
        call del_jtensor(jt)

        !$OMP END PARALLEL
        if (mpi_world_size > 1) then
            call gather_data(tens(:,first:last), this%tens(:,first:))
            if (mpi_rank == 0) then
                this%tens(:,1:first) = tens(:,1:first)
                deallocate(tens)
            end if
        end if
    end subroutine 
    
    subroutine calc_jvectors(this, mol, xdens, spin, z)
        type(jfield_t) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens
        character(*), optional :: spin
        integer, optional :: z

        character(8) :: spincase

        if (present(spin)) then
            spincase=spin
        else
            spincase='total'
        end if

        if (present(z)) then
            call calc_jtensors(this, mol, xdens, spincase, z)
        else
            call calc_jtensors(this, mol, xdens, spincase)
        end if
        ! until this point our J tensor should be the ACID T tensor !!
        if (mpi_rank == 0) then
            ! here the final contraction with B-field is done!!
            call compute_jvectors(this)
        end if
    end subroutine

    subroutine compute_jvectors(this)
        type(jfield_t) :: this

        integer(I4) :: k 
        integer(I4), dimension(2) :: dims
        dims = shape(this%vec)
        !$OMP PARALLEL DO PRIVATE(k) SHARED(this,dims)
        do k=1,dims(2)
            this%vec(:,k)=matmul(reshape(this%tens(:,k),(/3,3/)), this%b)
        end do
        !$OMP END PARALLEL DO
    end subroutine

! Make a rough estimate of how long the calculation will take
! (using one core) 
    subroutine jfield_eta(this, mol, xdens, fac)
        type(jfield_t) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens
        real(DP), intent(in), optional :: fac

        integer(I4) :: i, j, p1, p2, p3
        real(DP), dimension(9) :: foo
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
            foobar=matmul(bar,reshape(foo,(/3,3/)))
        end do
        call etime(times, tim2)
        tim2=times(1)
        call del_jtensor(jt)
        
        delta_t=tim2-tim1
        if ( present(fac) ) delta_t=delta_t*fac
        write(str_g, '(a,f9.2,a)') 'Estimated CPU time for single core &
            &calculation: ', delta_t*real(p1*p2*p3)/100.d0, ' sec'
        call msg_info(str_g)
        call nl
    end subroutine	

    function open_plot(fname) result(fd)
        character(*), intent(in) :: fname

        integer(I4) :: fd

        fd=0
        if (mpi_rank > 0) return
        if (trim(fname) == '') return

        call getfd(fd)
        if (fd == 0) then
            stop 1
        end if
        open(fd,file=trim(fname),status='unknown')
        return 
    end function

    subroutine jvector_plots(this, tag)
        type(jfield_t), intent(inout) :: this
        character(*), optional :: tag

        integer(I4) :: i, j, k, p1, p2, p3
        integer(I4) :: fd1, fd2, fd3, fd4
        real(DP), dimension(3) :: v, rr
        real(DP), dimension(:,:), pointer :: jv

        if (mpi_rank > 0) return

        if (present(tag)) then
            fd1 = open_plot('jvec_' // tag // '.txt')
            fd2 = open_plot('jmod_' // tag // '.txt')
        else
            fd1 = open_plot('jvec.txt')
            fd2 = open_plot('jmod.txt')
        end if

        call get_grid_size(this%grid, p1, p2, p3)
        jv=>this%vec
        do k=1,p3
            do j=1,p2
                do i=1,p1
                    rr=gridpoint(this%grid, i,j,k)*AU2A
                    v=jv(:,i+(j-1)*p1+(k-1)*p1*p2)*AU2A
                    call wrt_jvec(rr,v,fd1)
                    call wrt_jmod(rr,v,fd2)
                end do
                if (fd1 /= 0) write(fd1, *)
                if (fd2 /= 0) write(fd2, *)
            end do
        end do
        call closefd(fd1)
        call closefd(fd2)
        if (grid_is_3d(this%grid)) then 
            if (settings%acid) then
              ! add here ACID plot stuff ! 
              call acid_cube_plot(this)
            end if
            if (present(tag)) then
                call jmod_cubeplot(this, tag)
            else
                call jmod_cubeplot(this)
            end if
        end if
    end subroutine

    subroutine jmod_vtkplot(this)
        use vtkplot_module
        type(jfield_t) :: this

        integer(I4) :: p1, p2, p3, fd1, fd2
        integer(I4) :: i, j, k, l
        real(DP), dimension(3) :: qmin, qmax
        real(DP), dimension(3) :: norm, step, mag, v, rr
        real(DP) :: maxi, mini, sgn 
        integer(I4), dimension(3) :: npts
        real(DP), dimension(:,:), pointer :: buf
        real(DP), dimension(:,:,:), allocatable :: val

        if (mpi_rank > 0) return

        call get_grid_size(this%grid, p1, p2, p3)
        allocate(val(p1,p2,p3))
        npts=(/p1,p2,p3/)
        norm=get_grid_normal(this%grid)
        qmin=gridpoint(this%grid,1,1,1)
        qmax=gridpoint(this%grid,p1,p2,p3)

        step=(qmax-qmin)/(npts-1)

        mag = this%b

        buf => this%vec
        maxi=0.d0
        mini=0.d0
        l=0
        do i=1,p1
            do j=1,p2
                do k=1,p3
                    v=buf(:,i+(j-1)*p1+(k-1)*p1*p2)
                    rr=gridpoint(this%grid,i,j,k)
                    val(i,j,k)=(sqrt(sum(v**2)))
                end do
            end do
        end do
        call write_vtk_imagedata('jmod.vti', this%grid, val)
        deallocate(val)
    end subroutine

    subroutine jmod_cubeplot(this, tag)
        type(jfield_t) :: this
        character(*), optional :: tag

        integer(I4) :: p1, p2, p3, fd1, fd2
        integer(I4) :: i, j, k, l
        real(DP), dimension(3) :: qmin, qmax
        real(DP), dimension(3) :: norm, step, mag, v, rr
        real(DP) :: maxi, mini, val, sgn 
        integer(I4), dimension(3) :: npts
        real(DP), dimension(:,:), pointer :: buf

        if (mpi_rank > 0) return

        call jmod_vtkplot(this)
        call get_grid_size(this%grid, p1, p2, p3)
        npts=(/p1,p2,p3/)
        norm=get_grid_normal(this%grid)
        qmin=gridpoint(this%grid,1,1,1)
        qmax=gridpoint(this%grid,p1,p2,p3)

        step=(qmax-qmin)/(npts-1)

        if (present(tag)) then
            fd1=opencube('jmod_'// tag // '.cube', qmin, step, npts)
            fd2=opencube('jmod_quasi'// tag // '.cube', qmin, step, npts)
        else
            fd1=opencube('jmod.cube', qmin, step, npts)
            fd2=opencube('jmod_quasi.cube', qmin, step, npts)
        end if

        mag = this%b
        buf => this%vec
        maxi=0.d0
        mini=0.d0
        l=0
        do i=1,p1
            do j=1,p2
                do k=1,p3
                    v=buf(:,i+(j-1)*p1+(k-1)*p1*p2)
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
        print *, 'maxi, mini:', maxi, mini
        call closefd(fd1)
        call closefd(fd2)
    end subroutine

    function opencube(fname, origin, step, npts) result(fd)
        character(*), intent(in) :: fname
        real(DP), dimension(3), intent(in) :: origin, step
        integer(I4), dimension(3), intent(in) :: npts

        integer(I4) :: fd

        fd=0
        if (mpi_rank > 0) return
        if (trim(fname) == '') return

        call getfd(fd)
        if (fd == 0) then
            stop 1
        end if

        open(fd,file=trim(fname),form='formatted',status='unknown')
        write(fd,*) 'Gaussian cube data, generated by genpot'
        write(fd,*) 
        write(fd, '(i5,3f12.6)') 0, origin
        write(fd, '(i5,3f12.6)') npts(1), step(1), 0.d0, 0.d0
        write(fd, '(i5,3f12.6)') npts(2), 0.d0, step(2), 0.d0
        write(fd, '(i5,3f12.6)') npts(3), 0.d0, 0.d0, step(3)
    end function

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

    subroutine acid_cube_plot(this)
        type(jfield_t), intent(inout) :: this
        integer(I4), dimension(3) :: npts
        integer(I4) :: fd1, i, j, k, l, m, idx
        real(DP), dimension(:,:), pointer :: jtens
        real(DP), dimension(3) :: qmin, qmax, step, r 
        real(DP) :: val, maxi, mini
        ! collect grid information
        call get_grid_size(this%grid, npts(1), npts(2), npts(3))
        qmin = gridpoint(this%grid, 1, 1, 1)
        qmax = gridpoint(this%grid, npts(1), npts(2), npts(3))
        step = (qmax - qmin)/(npts - 1)
        ! get T tensor assume that this information is kept
        ! not sure because tens gets deallocated at one point...
        jtens => this%tens
        ! open cube file
        fd1 = opencube('acid.cube', qmin, step, npts)
        maxi = 0.0d0
        mini = 0.0d0
        l = 0
        m = 0
        idx = 0 
        do i = 1, npts(1)
           do j = 1, npts(2)
               do k = 1, npts(3)
                   m = m + 1
                   r = gridpoint(this%grid, i, j, k)
                   ! maybe better to call ctens here ?? check code !
                   ! which index for jtens? --> copy from cube plot file 
                   idx = i+(j-1)*npts(1) + (k-1)*npts(1)*npts(2) 
                   val = get_acid(r, jtens(:,idx))  
                   if (val > maxi) maxi = val 
                   if (val > mini) mini = val 
                   if (fd1 /= 0) then
                       write(fd1,'(f12.6)',advance='no') val
                       if (mod(l,6) == 5) write(fd1,*)
                   end if
                   l = l + 1
               end do
           end do
        end do
        print *, 'ACID: maxi, mini', maxi, mini

        call closefd(fd1)
    end subroutine
end module

! vim:et:sw=4:ts=4
