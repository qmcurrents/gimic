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
    use acid_module

    implicit none

    type jfield_t
        real(DP), dimension(3) :: b
        real(DP), dimension(:,:), pointer :: tens
        real(DP), dimension(:,:), pointer :: vec
        type(grid_t), pointer :: grid
    end type

    public new_jfield, del_jfield, calc_jtensors, calc_jvectors
    public jfield_t, jvector_plots, jfield_eta, compute_jvectors
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

        !$omp parallel default(none) &
        !$omp private(jt,rr,n,i,j,k) &
        !$omp shared(this,mol,xdens,spincase,tens,lo,hi)
        call new_jtensor(jt, mol, xdens)

        !$omp do schedule(static)
        do n=lo,hi
            call get_grid_index(this%grid, n, i, j, k)
            rr = gridpoint(this%grid, i, j, k)
            ! here the ACID T tensor is calculated and put on tens
            call ctensor(jt, rr, tens(:,n-lo+1), spincase)
        end do
        !$omp end do

        call del_jtensor(jt)
        !$omp end parallel

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
        !$omp parallel default(none) &
        !$omp shared(this,dims) &
        !$omp private(k)

        !$omp do
        do k=1,dims(2)
            this%vec(:,k)=matmul(reshape(this%tens(:,k),(/3,3/)), this%b)
        end do
        !$omp end do

        !$omp end parallel
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
        real(4) :: hours_per_sec = 1.0/3600.0
        real(4) :: days_per_hour = 1.0/24.0
        real(4) :: hours_per_day = 24.0
        integer :: no_tests = 100

        call get_grid_size(this%grid, p1, p2, p3)

        call new_jtensor(jt, mol, xdens)
        call etime(times, tim1)
        tim1=times(1)
        do i=1,no_tests
            call jtensor(jt, (/i*SC, i*SC, i*SC/), foo, spin_a)
            foobar=matmul(bar,reshape(foo,(/3,3/)))
        end do
        call etime(times, tim2)
        tim2=times(1)
        call del_jtensor(jt)

        delta_t=tim2-tim1
        if ( present(fac) ) delta_t=delta_t*fac
        write(str_g, '(a,f11.2,a,a,f6.1,a)') 'Estimated CPU time for single core &
            calculation: ', delta_t*real(p1*p2*p3)/real(no_tests), ' sec', &
            ' (',  delta_t*real(p1*p2*p3)/real(no_tests) * hours_per_sec, ' h )'
        if ( delta_t*real(p1*p2*p3)/real(no_tests) * hours_per_sec .gt. 2*hours_per_day ) then ! will take more than two days?
            write(str_g, '(a,f4.1,a)') '(', delta_t*real(p1*p2*p3)/real(no_tests) * hours_per_sec*days_per_hour, ' days )' 
        end if
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

    subroutine jvector_plots(this,tag)
        use vtkplot_module

        type(jfield_t), intent(inout) :: this
        character(*), optional :: tag
        logical :: circle_log
        logical :: debug 

        integer(I4) :: i, j, k, p1, p2, p3
        integer(I4) :: fd1, fd2, fd4, fd5
        integer(I4) :: idx, ptf
        real(DP), dimension(3) :: v, rr
        real(DP), dimension(3) :: center, normal
        real(DP), dimension(:,:), pointer :: jv
        real(DP), dimension(:,:), pointer :: jtens
        real(DP), dimension(:,:,:,:), allocatable :: jval
        real(DP) :: val
        real(DP) :: bound, r

        if (mpi_rank > 0) return
        ! get rid of print out of txt files but keep them for debugging
        ! keep also jmod.txt for 2D case to check integration planes

        debug = .false.
        if (debug) then
          if (present(tag)) then
              fd1 = open_plot('jvec' // tag // '.txt')
              fd2 = open_plot('jmod' // tag // '.txt')
          else
              fd1 = open_plot('jvec.txt')
              fd2 = open_plot('jmod.txt')
          end if
        end if
        ! print only jmod.txt for gauss type grid
        if (this%grid%gauss) then
          if (present(tag)) then
              fd2 = open_plot('jmod' // tag // '.txt')
          else
              fd2 = open_plot('jmod.txt')
          end if
        end if

        ! print out magnetic field
        write(*,*) "magnetic field"
        write(*,*) this%b
        write(*,*) ""
        ! 
        call get_grid_size(this%grid, p1, p2, p3)
        ! this is for cdens visualization when radius option is used
        call grid_center(this%grid,center)
        bound=1.d+10
        bound=this%grid%radius
        circle_log = .false.
        ! this only applies for grid base
        if (grid_is_3d(this%grid)) then
            circle_log = .false.
        else
            if (trim(this%grid%mode) .eq. 'bond') then
                if (this%grid%radius.gt.0.1d0) then
                    circle_log = .true.
                end if
            end if
        end if
        ! get grid normal vector n
        normal=get_grid_normal(this%grid)
        normal = normal*AU2A

        allocate(jval(p1,p2,p3,3))

        jv=>this%vec
        do k=1,p3
            do j=1,p2
                do i=1,p1
                    rr=gridpoint(this%grid, i,j,k)*AU2A
                    if (circle_log) then
                       ! for radius option
                       r = sqrt(sum((rr-center)**2))
                       if (r > bound) then
                           v = 0.0d0
                       else
                           v = jv(:,i+(j-1)*p1+(k-1)*p1*p2)
                       end if
                    else
                        v = jv(:,i+(j-1)*p1+(k-1)*p1*p2)
                    end if
                    ! collect jv vec information to put it on vti file
                    jval(i,j,k,1:3) = v
                    if (debug) then
                      call wrt_jvec(rr,v,fd1)
                      call wrt_jmod(rr,v,fd2)
                    end if
                    if (this%grid%gauss) then
                      call wrt_jmod(rr,v,fd2)
                    end if
                end do
                    if (debug) then
                      if (fd1 /= 0) write(fd1, *)
                      if (fd2 /= 0) write(fd2, *)
                    end if
                    if (this%grid%gauss) then
                      if (fd2 /= 0) write(fd2, *)
                    end if
            end do
        end do
        if (debug) then
          call closefd(fd1)
          call closefd(fd2)
        end if
        if (this%grid%gauss) then
          call closefd(fd2)
        end if
        ! case 3D grid
        if (grid_is_3d(this%grid)) then
          if (settings%acid) then
            call acid_vtkplot(this)
          end if
          ! put modulus info on file
          if (settings%jmod) then
            if (present(tag)) then
              call jmod2_vtkplot(this, tag)
            else
              call jmod2_vtkplot(this)
            end if
          end if
        end if
        ! case external 3D grid - numgrid
        if (settings%prop) then
            call get_property(this)
        end if

        ! put jvec information on vti file
        if (this%grid%gauss .or. trim(this%grid%mode)=='file') then
          write(*,*) "VTK files are only written for even grids"
        else 
          if (present(tag)) then
            call write_vtk_vector_imagedata('jvec'// tag // '.vti', this%grid, jval)
          else
            call write_vtk_vector_imagedata("jvec.vti", this%grid, jval)
          end if
        end if
        deallocate(jval)

    end subroutine

    subroutine jmod2_vtkplot(this, tag)
    ! based on jmod_vtkplot
    ! purpose: write all information on one file
        use vtkplot_module
        type(jfield_t) :: this
        character(*), optional :: tag

        integer(I4) :: p1, p2, p3
        integer(I4) :: i, j, k
        real(DP), dimension(3) :: norm, mag, v, rr
        real(DP) :: sgn
        real(DP), dimension(:,:), pointer :: buf
        real(DP), dimension(:,:,:), allocatable :: val

        if (mpi_rank > 0) return

        call get_grid_size(this%grid, p1, p2, p3)
        allocate(val(p1,p2,p3))

        mag = this%b

        buf => this%vec
        do k=1,p3
            do j=1,p2
                do i=1,p1
                    v = buf(:,i+(j-1)*p1+(k-1)*p1*p2)
                    rr=gridpoint(this%grid,i,j,k)
                    val(i,j,k)=(sqrt(sum(v**2)))
                    rr=rr-dot_product(mag,rr)*mag
                    norm=cross_product(mag,rr)
                    sgn=dot_product(norm,v)
                    if (sgn.lt.0.0d0) then
                        val(i,j,k) = -1.0d0*val(i,j,k)
                    end if
                end do
            end do
        end do
        if (present(tag)) then
          call write_vtk_imagedata('jmod'// tag // '.vti', this%grid, val)
        else
          call write_vtk_imagedata('jmod.vti', this%grid, val)
        end if
        deallocate(val)
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

    subroutine acid_vtkplot(this)
        use vtkplot_module
        type(jfield_t) :: this

        integer(I4) :: p1, p2, p3, idx
        integer(I4) :: i, j, k
        real(DP), dimension(:,:), pointer :: jtens
        real(DP), dimension(:,:,:), allocatable :: val

        if (mpi_rank > 0) return

        call get_grid_size(this%grid, p1, p2, p3)
        allocate(val(p1,p2,p3))

        jtens => this%tens
        idx = 0
        do k=1,p3
            do j=1,p2
                do i=1,p1
                    idx = i+(j-1)*p1 + (k-1)*p1*p2
                    val(i,j,k)= get_acid(jtens(:,idx))
                end do
            end do
        end do
        call write_vtk_imagedata('acid.vti', this%grid, val)
        deallocate(val)
    end subroutine

    subroutine get_property(this)
        ! assume external grid is read in correctly via read grid 
        ! assume grid_w.grd, coord.au, gridfile.grd are in the same directory
        type(jfield_t) :: this
        real(DP), dimension(:,:), pointer :: jtens
        real(DP), allocatable :: wg(:), coord(:,:), grd(:,:) 
        real(DP), dimension(3) :: d, bb, jvec, sigma, chi
        real(DP) :: f, tmp
        integer :: i, j, k, npts, natoms

        jtens => this%tens

        ! get atom coordinates in bohr
        open(unit=15, file="coord.au")
        natoms = getnlines(15)
        allocate(coord(natoms,3))
        do i=1, natoms
          read(15,*) coord(i,1:3)
        end do
        close(15)

        ! get weights from numgrid files (external)
        open(unit=15, file="grid_w.grd")
        open(unit=16, file="gridfile.grd")
        npts = getnlines(15)
        write(*,*) "npts", npts
        allocate(wg(npts),grd(npts,3))

        do i=1, npts
          read(15,*) wg(i)
          read(16,*) grd(i, 1:3)
        end do
        close(15)
        close(16)

        ! assume jtens has been calculated in the same order
        ! loope atoms
        do k=1, natoms
          ! loop grid points
          sigma = 0.0d0
          chi = 0.0d0
          do i=1, npts
            ! loop xyz
            do j = 1, 3
              d(j) = grd(i,j) - coord(k,j)   
            end do
            f = -1.0d0/((d(1)*d(1) + d(2)*d(2) + d(3)*d(3))**1.5d0)/(137.0359998d0**2.0d0)
            ! contract with Bx
            bb(1) = -1.0d0
            bb(2) = 0.0d0
            bb(3) = 0.0d0
            jvec = matmul(reshape(jtens(:,i),(/3,3/)),bb) 
            ! sigma_xx
            sigma(1) = sigma(1) + 1.0d6*wg(i)*f*(d(2)*jvec(3) - d(3)*jvec(2))
            ! chi_xx
            chi(1) = chi(1) + wg(i)*0.5d0*(grd(i,2)*jvec(3) - grd(i,3)*jvec(2))
            ! contract with By
            bb(1) = 0.0d0
            bb(2) = -1.0d0
            bb(3) = 0.0d0
            jvec = matmul(reshape(jtens(:,i),(/3,3/)),bb) 
            ! sigma_yy
            sigma(2) = sigma(2) + 1.0d6*wg(i)*f*(d(3)*jvec(1) - d(1)*jvec(3))
            ! chi_yy
            chi(2) = chi(2) + wg(i)*0.5d0*(grd(i,3)*jvec(1) - grd(i,1)*jvec(3))
            ! contract with Bz
            bb(1) = 0.0d0
            bb(2) = 0.0d0
            bb(3) = -1.0d0
            jvec = matmul(reshape(jtens(:,i),(/3,3/)),bb) 
            ! sigma_zz
            sigma(3) = sigma(3) + 1.0d6*wg(i)*f*(d(1)*jvec(2) - d(2)*jvec(1))
            ! chi_zz
            chi(3) = chi(3) + wg(i)*0.5d0*(grd(i,1)*jvec(2) - grd(i,2)*jvec(1))
          end do
          write(*,*) "atom ",k 
          write(*,*) "sigma_xx ", sigma(1) 
          write(*,*) "sigma_yy ", sigma(2) 
          write(*,*) "sigma_zz ", sigma(3) 
          tmp = (sigma(1) + sigma(2) + sigma(3))/3.0d0  
          write(*,*) "shielding constant sigma = ", tmp 

        end do
        write(*,*) "" 
        write(*,*) "chi_xx ", chi(1) 
        write(*,*) "chi_yy ", chi(2) 
        write(*,*) "chi_zz ", chi(3) 
        tmp = (chi(1) + chi(2) + chi(3))/3.0d0  
        write(*,*) "isotropic magnetizability chi = ", tmp 

        ! clean up 
        deallocate(wg, coord)

    end subroutine
end module

! vim:et:sw=4:ts=4
