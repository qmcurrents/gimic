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
        real(DP), dimension(3) :: coord
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
        !$omp private(jt,coord,n,i,j,k) &
        !$omp shared(this,mol,xdens,spincase,tens,lo,hi)
        call new_jtensor(jt, mol, xdens)

        !$omp do schedule(static)
        do n=lo,hi
            call get_grid_index(this%grid, n, i, j, k)
            coord = gridpoint(this%grid, i, j, k)
            ! here the ACID T tensor is calculated and put on tens
            call ctensor(jt, coord, tens(:,n-lo+1), spincase)
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

        implicit none

        type(jfield_t), intent(inout) :: this
        character(*), optional :: tag
        logical :: circle_log
        logical :: debug

        integer(I4) :: i, j, k, p1, p2, p3
        integer(I4) :: fd1, fd2, fd4, fd5
        integer(I4) :: idx, ptf
        real(DP), dimension(3)             :: v, coord
        real(DP), dimension(3)             :: center, normal
        real(DP), dimension(:,:), pointer  :: jv
        real(DP), dimension(:,:), pointer  :: jtens
        real(DP), allocatable              :: jval_regular(:, :, :, :)
        real(real64), allocatable          :: jval_unstructured(:, :)
        integer(int32), allocatable        :: cells(:, :) ! ncells x 4

        real(DP) :: val
        real(DP) :: bound, r
        integer(int32) :: dummy, ncells
        logical                            :: fd2_isopen = .false.
        logical                            :: elements_exists

        if (mpi_rank > 0) return
        ! get rid of print out of txt files but keep them for debugging
        ! keep also jmod.txt for 2D case to check integration planes

        debug = .false.
! lnw: outcommented because inactive anyway
!        if (debug) then
!          if (present(tag)) then
!              fd1 = open_plot('jvec' // tag // '.txt')
!              fd2 = open_plot('jmod' // tag // '.txt')
!          else
!              fd1 = open_plot('jvec.txt')
!              fd2 = open_plot('jmod.txt')
!          end if
!        end if

        ! print only jmod.txt for gauss type grid
        if (this%grid%gauss) then
          if (present(tag)) then
              fd2 = open_plot('jmod' // tag // '.txt')
          else
              fd2 = open_plot('jmod.txt')
          end if
          fd2_isopen = .true.
        end if

        ! print out magnetic field
        write(*,*) "magnetic field"
        write(*,*) this%b
        write(*,*) ""

        call get_grid_size(this%grid, p1, p2, p3)

        if(trim(this%grid%mode) == 'bond' .or. trim(this%grid%mode) == 'base' .or. trim(this%grid%mode) == 'std') then ! bond+base only, not file
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
          normal = normal*AU2A  ! FIXME: normal is never used
        end if

        jv=>this%vec
        if ( ((trim(this%grid%mode)=='std' .or. trim(this%grid%mode)=='base') .and. this%grid%gtype=='even') &
             .or. trim(this%grid%mode)=='bond') then ! (base && even) || bond
          allocate(jval_regular(p1,p2,p3,3))
          do k=1,p3
            do j=1,p2
              do i=1,p1
                coord=gridpoint(this%grid, i,j,k)*AU2A
                if (circle_log) then
                  ! for radius option
                  r = sqrt(sum((coord-center)**2))
                  if (r > bound) then
                    v = 0.0d0
                  else
                    v = jv(:,i+(j-1)*p1+(k-1)*p1*p2)
                  end if
                else
                  v = jv(:,i+(j-1)*p1+(k-1)*p1*p2)
                end if
                ! collect jv vec information to put it on vti file
                jval_regular(i,j,k,1:3) = v
!                if (debug) then
!                  call write_jvec(coord,v,fd1)
!                  call write_jmod(coord,v,fd2)
!                end if
                if (this%grid%gauss) then
                  call write_jmod(coord,v,fd2)
                end if
              end do
!              if (debug) then
!                if (fd1 /= 0) write(fd1, *)
!                if (fd2 /= 0) write(fd2, *)
!              end if
              if (fd2_isopen) write(fd2, *)
            end do
          end do
        else if ( ((trim(this%grid%mode)=='base' .or. trim(this%grid%mode)=='std') .and. this%grid%gauss) &
              .or. trim(this%grid%mode)=='file' ) then ! (base && gauss) || file
          allocate(jval_unstructured(this%grid%npts(1), 3))
          dummy = 0
          do i=1, p1 ! which is the total number of points
            v = jv(:, i)
            if ( (trim(this%grid%mode)=='base' .or. trim(this%grid%mode)=='std') .and. this%grid%gauss) then ! write jmod file for base && gauss
              coord=gridpoint(this%grid, i, dummy, dummy)*AU2A
              call write_jmod(coord,v,fd2)
            end if
            jval_unstructured(i, 1:3) = v
          end do
        end if

!        if (debug) then
!          call closefd(fd1)
!        end if
        if (fd2_isopen) then
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

        ! write jvec information to vti file
        if ( (trim(this%grid%mode)=='base' .or. trim(this%grid%mode)=='std' .or. trim(this%grid%mode)=='bond') &
                  .and. this%grid%gtype=='even' ) then
          if (present(tag)) then
            call write_vtk_vector_imagedata('jvec'// tag // '.vti', this%grid, jval_regular)
          else
            call write_vtk_vector_imagedata("jvec.vti", this%grid, jval_regular)
          end if
          deallocate(jval_regular)
        ! write jvec information to vtu file
        else if ( ((trim(this%grid%mode)=='base' .or. trim(this%grid%mode)=='std') .and. this%grid%gauss) &
              .or. trim(this%grid%mode)=='file' ) then

          ! read element file: first number in first line is number of lines
          inquire(FILE='grid.1.ele', EXIST=elements_exists)
          if(elements_exists) then
            open(GRIDELE, file='grid.1.ele')
            ! read(GRIDELE, '(3i4)') ncells, dummy, dummy
            read(GRIDELE, *) ncells, dummy, dummy
            ! write(*,*) "ncells ", ncells
            allocate(cells(4,ncells))
            do i=1,ncells
              ! read(GRIDELE,'(5i7)') dummy, cells(1:4, i)
              read(GRIDELE, *) dummy, cells(1:4, i)
            end do
            close(GRIDELE)

            call write_vtk_vector_unstructuredgrid("jvec.vtu", this%grid%xdata, jval_unstructured, cells)
            deallocate(cells)
          else
            write(*,*) 'not writing a vtu file, because the file grid.1.ele was not found.'
          endif
          deallocate(jval_unstructured)
        end if

    end subroutine


    subroutine jmod2_vtkplot(this, tag)
    ! based on jmod_vtkplot
    ! purpose: write all information on one file
        use vtkplot_module
        type(jfield_t) :: this
        character(*), optional :: tag

        integer(I4) :: p1, p2, p3
        integer(I4) :: i, j, k
        real(DP), dimension(3) :: norm, mag, v, coord
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
                    coord=gridpoint(this%grid,i,j,k)
                    val(i,j,k)=(sqrt(sum(v**2)))
                    coord=coord-dot_product(mag,coord)*mag
                    norm=cross_product(mag,coord)
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

    subroutine write_jvec(coord, v, fd)
        real(DP), dimension(:), intent(in) :: coord, v
        integer(I4), intent(in)            :: fd

        if (fd == 0) return

        write(fd, '(6f11.7)')  coord, v
    end subroutine

!    subroutine write_njvec(coord,v,fd)
!        real(DP), dimension(:), intent(in) :: coord, v
!        integer(I4), intent(in) :: fd
!
!        real(DP) :: nfac
!
!        if (fd == 0) return
!
!        nfac=sqrt(sum(v(:)**2))
!        if (nfac < 1.d-15) then
!            write(fd, '(6f11.7)')  coord, 0.d0, 0.d0, 0.d0
!        else
!            write(fd, '(6f11.7)')  coord, v/nfac
!        end if
!    end subroutine

!    subroutine write_jproj(coord,v,grid,fd)
!        real(DP), dimension(:), intent(in) :: coord, v
!        type(grid_t) :: grid
!        integer(I4), intent(in) :: fd
!
!        real(DP) :: jprj
!        real(DP), dimension(3) :: norm
!
!        if (fd == 0) return
!        norm=get_grid_normal(grid)
!
!        jprj=dot_product(norm,v)
!        write(fd, '(3e19.12)') coord, jprj
!    end subroutine

    subroutine write_jmod(coord,v,fd)
        real(DP), dimension(:), intent(in) :: coord, v
        integer(I4), intent(in) :: fd

        real(DP) :: jmod

        if (fd == 0) return

        jmod=sqrt(sum(v**2))
        write(fd, '(6f11.7)')  coord, jmod
    end subroutine

    subroutine print_jt(coord, jt)
        real(DP), dimension(3), intent(in)   :: coord
        real(DP), dimension(3,3), intent(in) :: jt

        integer(I4) :: l
        print *
        print '(a,3f12.8)', 'Current tensor at ', coord
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
        real(DP), dimension(:,:), pointer       :: jtens
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
      use vtkplot_module
        ! assume external grid is read in correctly via read grid
        ! assume grid_w.grd, coord.au, gridfile.grd are in the same directory
        ! for visualization of unstructured grid file grid.1.ele is needed
        ! assume this has been generated with tetgen and is in the same directory
        type(jfield_t) :: this
        real(DP), dimension(:,:), pointer :: jtens
        real(DP), allocatable :: wg(:), coord(:,:), grd(:,:), grd_tmp(:,:)
        real(DP), allocatable :: intchi(:,:), pdata(:), intsigma(:,:)  
        real(DP), allocatable :: scont(:,:) 
        real(DP), dimension(3) :: d, bb, jvec, sigma, chi, csum
        real(DP), dimension(3) :: tmpscont, R0
        real(DP) :: f, tmp, si_tmp, spos, sneg
        integer :: i, j, k, npts, natoms, dummy, ncells
        integer :: iatom, iatom_npts, l
        integer(int32), allocatable :: cells(:,:) ! for vtu file
        integer(int32), allocatable :: nelpts(:,:) 
        logical :: coords_exists, points_exists, weights_exists
        logical :: elements_exists, nelpts_exists
        logical :: testflag 
        real(DP), parameter :: fac_au2simag = 7.89104d-29
        character(len=70) :: filename
        character(len=70) :: formatstring

        ! testflag to check if int J d^3 is indeed zero
        testflag = .false.

        inquire(FILE='coord.au',     EXIST=coords_exists)
        inquire(FILE='gridfile.grd', EXIST=points_exists)
        inquire(FILE='grid_w.grd',   EXIST=weights_exists)
        inquire(FILE='nelpts.info',  EXIST=nelpts_exists)
        if(.not. (coords_exists .and. points_exists .and. weights_exists)) then
          write(*,*) 'at least one of the files coord.au, gridfile.grd, and grid_w.grd is missing.',&
                     'Therefore any property calculation is skipped.'
          return
        endif

        jtens => this%tens

        ! get atom number X and generated grid points for X
        ! needed for atom decomposition of contributions
        open(unit=15, file="nelpts.info")
        natoms = getnlines(15)
        allocate(nelpts(natoms,2)) 
        do i = 1, natoms
          read(15,*) nelpts(i,1:2)
          ! write(*,*) nelpts(i,1:2)
        end do
        close(15)

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

        ! test if int J d^3 eq. zero
        if (testflag) then
          allocate(grd_tmp(npts,3))
          ! hardcode here gauge origin
          R0(1) = 100.0d0 
          R0(2) = 100.0d0 
          R0(3) = 100.0d0 
          do i=1, npts
            grd_tmp(:,1) = grd(:,1) - R0(1) 
            grd_tmp(:,2) = grd(:,2) - R0(2) 
            grd_tmp(:,3) = grd(:,3) - R0(3) 
          end do
          ! check if int j d^3 eq. zero
          !grd_tmp = 1.0d0
          grd = grd_tmp
          deallocate(grd_tmp) 
        end if

        ! get cell information for printing vtu file
        inquire(FILE='grid.1.ele', EXIST=elements_exists)
        if(elements_exists) then
          open(GRIDELE, file='grid.1.ele')
          read(GRIDELE, *) ncells, dummy, dummy
          allocate(cells(4,ncells))
          do i=1,ncells
            read(GRIDELE, *) dummy, cells(1:4, i)
          end do
          close(GRIDELE)
        end if 

        ! assume jtens has been calculated in the same order
        ! calculate shielding tensor diagonal elements
        allocate(intsigma(npts,3),pdata(npts))
        allocate(scont(natoms,3))
        ! loop  atoms
        do k=1, natoms
          ! loop grid points
          sigma = 0.0d0
          spos  = 0.0d0
          sneg  = 0.0d0
          iatom = 1
          scont = 0.0d0
          iatom_npts = nelpts(iatom,2)
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
            intsigma(i,1) =  1.0d6*f*(d(2)*jvec(3) - d(3)*jvec(2))
            sigma(1) = sigma(1) + 1.0d6*wg(i)*f*(d(2)*jvec(3) - d(3)*jvec(2))
            ! contract with By
            bb(1) = 0.0d0
            bb(2) = -1.0d0
            bb(3) = 0.0d0
            jvec = matmul(reshape(jtens(:,i),(/3,3/)),bb)
            ! sigma_yy
            intsigma(i,2) =  1.0d6*f*(d(3)*jvec(1) - d(1)*jvec(3))
            sigma(2) = sigma(2) + 1.0d6*wg(i)*f*(d(3)*jvec(1) - d(1)*jvec(3))
            ! contract with Bz
            bb(1) = 0.0d0
            bb(2) = 0.0d0
            bb(3) = -1.0d0
            jvec = matmul(reshape(jtens(:,i),(/3,3/)),bb)
            ! sigma_zz
            intsigma(i,3) = 1.0d6*f*(d(1)*jvec(2) - d(2)*jvec(1))
            sigma(3) = sigma(3) + 1.0d6*wg(i)*f*(d(1)*jvec(2) - d(2)*jvec(1))
            pdata(i) = (intsigma(i,1) + intsigma(i,2) + intsigma(i,3)) 
            ! collect positive and negative contributions
            if (pdata(i).ge.0.0d0) then
                spos = spos + pdata(i)*wg(i)
            else 
                sneg = sneg + pdata(i)*wg(i)
            end if
            ! collect atom contributions 
            if (i == iatom_npts) then
                scont(iatom,1) = (sigma(1)+sigma(2)+sigma(3))/3.0d0  
                scont(iatom,2) = (spos)/3.0d0  
                scont(iatom,3) = (sneg)/3.0d0  
                !write(*,*) scont(iatom, 1), scont(iatom, 2), scont(iatom,3)
                iatom = iatom + 1 
                if (iatom.le.natoms) then
                  iatom_npts = iatom_npts + nelpts(iatom,2) 
                end if 
            end if
            !
          end do ! loop npts
          !
          write(*,*) "atom ",k
          write(*,*) "in ppm" 
          write(*,"(X,A10,2X,F14.6)") "sigma_xx ", sigma(1)
          write(*,"(X,A10,2X,F14.6)") "sigma_yy ", sigma(2)
          write(*,"(X,A10,2X,F14.6)") "sigma_zz ", sigma(3)
          tmp = (sigma(1) + sigma(2) + sigma(3))/3.0d0
          write(*,"(A30,2X,F14.6)") "shielding constant    = ", tmp
          write(*,"(A30,2X,F14.6)") "positive contribution = ", spos/3.0d0
          write(*,"(A30,2X,F14.6)") "negative contribution = ", sneg/3.0d0
          write(*,"(A30,2X,F14.6)") "sum = ", (spos + sneg)/3.0d0
          write(*,*) " "
          write(*,*) "atom contributions, total, positive, negative" 
          csum = 0.0d0
          tmpscont = 0.0d0
          do l=1, natoms
            if (l==1) then
              write(*,"(A5,I5,3F14.6)") "atom ", l, scont(l,1),scont(l,2),scont(l,3)
              csum(1) = csum(1) + scont(l,1) 
              csum(2) = csum(2) + scont(l,2) 
              csum(3) = csum(3) + scont(l,3) 
            else 
              tmpscont(1) = scont(l,1) - scont(l-1,1) 
              tmpscont(2) = scont(l,2) - scont(l-1,2) 
              tmpscont(3) = scont(l,3) - scont(l-1,3) 
              write(*,"(A5,I5,3F14.6)") "atom ", l, tmpscont(1), tmpscont(2), tmpscont(3)  
              csum(1) = csum(1) + tmpscont(1) 
              csum(2) = csum(2) + tmpscont(2) 
              csum(3) = csum(3) + tmpscont(3) 
            end if
          end do
          write(*,"(A10, 3F14.6)") "sum ", csum(1), csum(2), csum(3)
          write(*,*) "****************************************************"

          ! plot integrand intsigma
          if(elements_exists) then
            if (k < 10) then
              formatstring = "(A5,I1,A4)"
            else
              formatstring = "(A5,I2,A4)"
            endif
            write (filename,formatstring) "sigma", k, ".vtu"
            ! print *, trim(filename)
            call write_vtk_scalar_unstructuredgrid(trim(filename), this%grid%xdata, pdata, cells)
            if (k < 10) then
              formatstring = "(A8,I1,A4)"
            else
              formatstring = "(A8,I2,A4)"
            endif
            write (filename,formatstring) "sigma_xx", k, ".vtu"
            print *, filename 
            call write_vtk_scalar_unstructuredgrid(trim(filename), this%grid%xdata, intsigma(:,1), cells)
            write (filename,formatstring) "sigma_yy", k, ".vtu"
            print *, filename 
            call write_vtk_scalar_unstructuredgrid(trim(filename), this%grid%xdata, intsigma(:,2), cells)
            write (filename,formatstring) "sigma_zz", k, ".vtu"
            print *, filename 
            call write_vtk_scalar_unstructuredgrid(trim(filename), this%grid%xdata, intsigma(:,3), cells)
          end if
        end do !loop atoms
        write(*,*) ""
        
        ! calculate diagonal elements of magnetizabilitz tensor
        ! for magnetizability no loop over coords is needed 
        allocate(intchi(npts,3))
        chi = 0.0d0
        spos = 0.0d0
        sneg = 0.0d0
        iatom = 1
        scont = 0.0d0
        iatom_npts = nelpts(iatom,2)
        do i=1, npts
          ! loop xyz
          ! contract with Bx
          bb(1) = -1.0d0
          bb(2) = 0.0d0
          bb(3) = 0.0d0
          jvec = matmul(reshape(jtens(:,i),(/3,3/)),bb)
          ! chi_xx
          intchi(i,1) = 0.5d0*(grd(i,2)*jvec(3) - grd(i,3)*jvec(2))
          chi(1) = chi(1) + wg(i)*0.5d0*(grd(i,2)*jvec(3) - grd(i,3)*jvec(2))
          ! contract with By
          bb(1) = 0.0d0
          bb(2) = -1.0d0
          bb(3) = 0.0d0
          jvec = matmul(reshape(jtens(:,i),(/3,3/)),bb)
          ! chi_yy
          intchi(i,2) = 0.5d0*(grd(i,3)*jvec(1) - grd(i,1)*jvec(3))
          chi(2) = chi(2) + wg(i)*0.5d0*(grd(i,3)*jvec(1) - grd(i,1)*jvec(3))
          ! contract with Bz
          bb(1) = 0.0d0
          bb(2) = 0.0d0
          bb(3) = -1.0d0
          jvec = matmul(reshape(jtens(:,i),(/3,3/)),bb)
          ! chi_zz
          intchi(i,3) = 0.5d0*(grd(i,1)*jvec(2) - grd(i,2)*jvec(1))
          chi(3) = chi(3) + wg(i)*0.5d0*(grd(i,1)*jvec(2) - grd(i,2)*jvec(1))
          ! calculate sum of integrands
          pdata(i) = (intchi(i,1) + intchi(i,2) + intchi(i,3)) 
          ! collect positive and negative contributions
          if (pdata(i).ge.0.0d0) then
              spos = spos + pdata(i)*wg(i)
          else 
              sneg = sneg + pdata(i)*wg(i)
          end if
          ! collect atom contributions 
          if (i == iatom_npts) then
              scont(iatom,1) = (chi(1)+chi(2)+chi(3))/3.0d0  
              scont(iatom,2) = (spos)/3.0d0  
              scont(iatom,3) = (sneg)/3.0d0  
              !write(*,*) scont(iatom, 1), scont(iatom, 2), scont(iatom,3)
              iatom = iatom + 1 
              if (iatom.le.natoms) then
                iatom_npts = iatom_npts + nelpts(iatom,2) 
              end if 
          end if
          !
        end do !npts
        ! 
        write(*,*) ""
        write(*,"(X,A7,2X,F14.8)") "chi_xx ", chi(1)
        write(*,"(X,A7,2X,F14.8)") "chi_yy ", chi(2)
        write(*,"(X,A7,2X,F14.8)") "chi_zz ", chi(3)
        tmp = (chi(1) + chi(2) + chi(3))/3.0d0
        write(*,*) "in au"
        write(*,"(X,A30,2X,F14.6)") "isotropic magnetizability chi = ", tmp
        write(*,"(X,A30,2X,F14.6)") "positive contribution         = ", spos/3.0d0
        write(*,"(X,A30,2X,F14.6)") "negative contribution         = ", sneg/3.0d0
        write(*,"(X,A30,2X,F14.6)") "sum ", (spos + sneg)/3.0d0
        write(*,*) ""
        si_tmp = tmp*fac_au2simag 
        write(*,*) "in SI units J/T^2 " 
        write(*,*) "conversion factor: 7.89104*10^-29 J/T^2 " 
        write(*,*) ""
        write(*,"(A30,2X,E14.6)") "isotropic magnetizability = ", si_tmp 
        write(*,"(A30,2X,E14.6)") "positive contribution     = ", (spos/3.0d0)*fac_au2simag
        write(*,"(A30,2X,E14.6)") "negative contribution     = ", (sneg/3.0d0)*fac_au2simag
        write(*,"(A30,2X,E14.6)") "sum ", ((spos + sneg)/3.0d0)*fac_au2simag
        write(*,*) "****************************************************"
        write(*,*) ""
        write(*,*) "atom contributions, total, positive, negative" 
        csum = 0.0d0
        tmpscont = 0.0d0
        do l=1, natoms
          if (l==1) then
            write(*,"(A5,I5,3F14.6)") "atom ", l, scont(l,1),scont(l,2),scont(l,3)
            csum(1) = csum(1) + scont(l,1) 
            csum(2) = csum(2) + scont(l,2) 
            csum(3) = csum(3) + scont(l,3) 
          else 
            tmpscont(1) = scont(l,1) - scont(l-1,1) 
            tmpscont(2) = scont(l,2) - scont(l-1,2) 
            tmpscont(3) = scont(l,3) - scont(l-1,3) 
            write(*,"(A5,I5,3F14.6)") "atom ", l, tmpscont(1), tmpscont(2), tmpscont(3)  
            csum(1) = csum(1) + tmpscont(1) 
            csum(2) = csum(2) + tmpscont(2) 
            csum(3) = csum(3) + tmpscont(3) 
          end if
        end do
        write(*,"(A10, 3F14.6)") "sum ", csum(1), csum(2), csum(3)
        write(*,*) "****************************************************"

        ! plot integrand intchi 
        if(elements_exists) then
          call write_vtk_scalar_unstructuredgrid("intchi.vtu", this%grid%xdata, pdata, cells)
          call write_vtk_scalar_unstructuredgrid("intchi_xx.vtu", this%grid%xdata, intchi(:,1), cells)
          call write_vtk_scalar_unstructuredgrid("intchi_yy.vtu", this%grid%xdata, intchi(:,2), cells)
          call write_vtk_scalar_unstructuredgrid("intchi_zz.vtu", this%grid%xdata, intchi(:,3), cells)
        end if 


        ! clean up
        deallocate(wg, coord, nelpts)
        deallocate(intchi, pdata, intsigma)
        if(elements_exists) then
          deallocate(cells)
        end if 

    end subroutine
end module

! vim:et:sw=4:ts=4
