! Setup grids for all possible purposes.
! Written by Barnabas Plomlund.
!
! This module is bit of a mess...
! 

module grid_class
    use globals_m
    use settings_m
    use tensor_m
    use gaussint_m
    use basis_class
    use teletype_m
    implicit none

    type grid_t
        logical :: gauss  ! integration grid, 1,2 or 3-D
        real(DP), dimension(3,3) :: basv   ! grid basis vectors
        real(DP), dimension(3) :: l        ! |v|
        real(DP), dimension(3) :: origin, ortho
        real(DP), dimension(3) :: step
        integer(I4), dimension(3) :: npts
        type(gdata_t), dimension(3) :: gdata
        character(BUFLEN) :: mode, gtype
        real(DP), dimension(:,:), pointer :: xdata
    end type grid_t

    public grid_t
    public new_grid, del_grid, gridpoint, get_grid_normal
    public get_grid_size, get_weight, write_grid, read_grid
    public get_grid_length, is_gauss_grid, realpoint, copy_grid
    public grid_center, plot_grid_xyz, get_basvec, get_ortho
    public get_grid_range

    interface get_basvec
        module procedure get_basv1
        module procedure get_basv3
    end interface

    private
    real(DP), parameter :: DPTOL=1.d-10
    type(getkw_t), pointer :: input
contains

    subroutine new_grid(this, inp, mol)
        type(grid_t) :: this
        type(getkw_t), target :: inp
        type(molecule_t) :: mol

        real(DP) :: ll 
        real(DP), dimension(3) :: angle
        integer(I4), dimension(3) :: ngp
        integer(I4) :: i, j

        input=>inp
        
        this%step=1.d0
        this%gtype='even'
        this%basv=0.0
        this%l=0.0
        this%ortho=0.d0
        call getkw(input, 'grid', this%mode)

        ! first figure out where and how to place the grid
        select case (trim(this%mode))
            case ('file')
                call extgrid(this)
                return
            case ('std','base')
                call setup_std_grid(this)
            case ('bond')
                call setup_bond_grid(this,mol)
            case default
                call msg_error('Unknown grid type: ' // trim(this%mode))
                stop
        end select

        call msg_out('Grid mode = ' // trim(this%mode))

        call normalise(this%basv)
        call ortho_coordsys(this)
        
        ! rotate basis vectors if needed
        if (keyword_is_set(input, 'grid.rotation')) then
            call getkw(input, 'grid.rotation', angle)
            call rotate(this, angle)
        end if

        ! calculate distibution of grid points
        call getkw(input, 'grid.type', this%gtype)
        i=len(trim(this%gtype))
        select case (this%gtype(1:i))
            case ('even')
                call setup_even_grid(this)
            case ('gauss')
                call setup_gauss_grid(this, 'gauss')
            case ('lobato')
                call setup_gauss_grid(this, 'lobato')
            case default
                call msg_error('Unknown grid type: ' // trim(this%gtype))
                stop
        end select

        write(str_g, '(2x,a,3i5)') 'Number of grid points <v1,v2>:', &
            this%npts(1), this%npts(2), this%npts(3)
        call msg_out(str_g)
        write(str_g, '(2x,a,i10)') 'Total number of grid points  :', &
            product(this%npts)
        call msg_out(str_g)
        call nl
!        call write_grid(this, 50)

    end subroutine

    subroutine setup_std_grid(this)
        type(grid_t) :: this

        integer(I4) :: i
        real(DP), dimension(3) :: normv

        call getkw(input, 'grid.origin', this%origin)
        call getkw(input, 'grid.ivec', this%basv(:,1))
        call getkw(input, 'grid.jvec', this%basv(:,2))
        call getkw(input, 'grid.lengths', this%l)

        if (keyword_is_set(input, 'grid.spacing')) then
            call getkw(input, 'grid.spacing', this%step)
        else
            call getkw(input, 'grid.grid_points', this%npts)
            this%step=this%l/(this%npts-1)
        end if

        this%basv(:,3)=cross_product(this%basv(:,1),this%basv(:,2)) 
    end subroutine

    subroutine setup_bond_grid(this,mol)
        type(grid_t) :: this
        type(molecule_t) :: mol

        type(atom_t), pointer :: atom

        integer(I4) :: i
        real(DP), dimension(3) :: normv
        integer(I4), dimension(3) :: atoms
        real(DP), dimension(3) :: v1, v2, v3, oo
        real(DP), dimension(2) :: hgt, wdt
        real(DP) :: l3

        if (keyword_is_set(input, 'grid.bond')) then
            call getkw(input, 'grid.bond', atoms(1:2))
            call get_atom(mol, atoms(1), atom)
            call get_coord(atom, this%basv(:,1))
            call get_atom(mol, atoms(2), atom)
            call get_coord(atom, this%basv(:,2))
        else
            call getkw(input, 'grid.coord1', this%basv(:,1))
            call getkw(input, 'grid.coord2', this%basv(:,2))
        end if

        if (keyword_is_set(input,'grid.fixpoint')) then
            call getkw(input, 'grid.fixpoint', atoms(3))
            call get_atom(mol, atoms(3), atom)
            call get_coord(atom, this%origin)
        else
            call getkw(input, 'grid.fixcoord', this%origin)
        end if

        !defaults, etc.
        l3=-1.d0
        wdt=-1.d0
        hgt=-1.d0
        call getkw(input, 'grid.distance', l3)
        if (keyword_is_set(input, 'grid.height')) then
            call getkw(input, 'grid.height', hgt)
            call getkw(input, 'grid.width', wdt)
            hgt(1)=-hgt(1)
            wdt(1)=-wdt(1)
        else
            call getkw(input, 'grid.in', wdt(1))
            call getkw(input, 'grid.out', wdt(2))
            call getkw(input, 'grid.up', hgt(1))
            call getkw(input, 'grid.down', hgt(2))
        end if
        this%l=(/sum(hgt), sum(wdt), 0.d0/)

        if ( l3 < 0.d0 ) then
            call msg_warn('grid.distance < 0!')
        end if
        if ( sum(wdt) < 0.d0 ) then
            call msg_critical('Grid width < 0!')
        end if
        if ( sum(hgt) < 0.d0 ) then
            call msg_critical('Grid heigth < 0!')
        end if

        ! figure out "orthogonal" axis for magnetic field
        v1=this%basv(:,1)-this%origin
        v2=this%basv(:,2)-this%origin
        this%ortho=cross_product(v1,v2)
        if (vcmp(this%ortho, NILL_VECTOR)) then
            call msg_error('Basis vectors are linearly dependent, field &
            &direction undefined!')
            stop
        end if
        this%ortho=norm(this%ortho)
    
        v3=norm(v2-v1)
        v1=-this%ortho
        v2=norm(cross_product(v3, v1))
        oo=this%basv(:,1)+l3*v3
        this%origin=oo-wdt(2)*v2-hgt(2)*v1

        this%basv(:,1)=v1
        this%basv(:,2)=v2
        this%basv(:,3)=v3

        call nl
        call msg_out('Integration grid data')
        call msg_out('------------------------------------------------')
        write(str_g, '(a,3f12.6)') 'center ', oo
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') 'origin ', this%origin
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') 'basv1  ', this%basv(:,1)
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') 'basv2  ', this%basv(:,2)
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') 'basv3  ', this%basv(:,3)
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') 'lenghts', this%l
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') 'magnet ', this%ortho
        call msg_out(str_g)
        call nl
    end subroutine
    
    subroutine normalise(v)
        real(DP), dimension(:,:), intent(inout) :: v

        integer :: i
        real(DP) :: norm

        do i=1,3
            norm=sqrt(sum(v(:,i)**2))
            if (norm > 0.d0) v(:,i)=v(:,i)/norm
        end do
    end subroutine


    subroutine setup_gauss_grid(this, quadr)
        type(grid_t) :: this
        character(*), intent(in) :: quadr

        integer(I4) ::  i, rem, order
        real(DP), dimension(3) :: spc
        logical :: flag

        flag=.false.

        call msg_info('Integration grid selected.')
        call getkw(input, 'grid.gauss_order', order)
        if (keyword_is_set(input,'grid.grid_points')) then
            call getkw(input, 'grid.grid_points', this%npts)
        else
            call getkw(input, 'grid.spacing', spc)
            do i=1,3
                if (abs(spc(i)) < 1.d-10 .or.spc(i) < 0.d0) then
                    this%npts(i)=0
                else
                    this%npts(i)=nint(this%l(i)/spc(i))
                end if
            end do
        end if
        
        do i=1,3
            if (.not.this%npts(i) > 1) then
                this%npts(i)=0
            end if
            rem=mod(this%npts(i),order)
            if (rem /= 0) then
                this%npts(i)=this%npts(i)-rem+order
                flag=.true.
            end if
        end do
        
        if (flag) then
            write(str_g, '(a,3i5)'), &
            'Adjusted number of grid points for quadrature: ', this%npts
            call msg_info(str_g)
        end if

        this%gauss=.true.
        do i=1,3
            if (this%npts(i) > 0) then
                allocate(this%gdata(i)%pts(this%npts(i)))
                allocate(this%gdata(i)%wgt(this%npts(i)))
                call setup_gauss_data(0.d0, this%l(i), order, &
                    this%gdata(i), quadr)
            else
                this%npts(i)=1
                allocate(this%gdata(i)%pts(1))
                allocate(this%gdata(i)%wgt(1))
                call setup_gauss_data(0.d0, this%l(i), 1, &
                    this%gdata(i), quadr)
            end if
        end do
    end subroutine

    subroutine setup_even_grid(this)
        type(grid_t), intent(inout) :: this

        integer(I4) :: i, n

        this%gauss=.false.
        this%npts(1)=nint(this%l(1)/this%step(1))+1
        this%npts(2)=nint(this%l(2)/this%step(2))+1
        if (this%l(3) == 0.d0 .or. this%step(3) == 0.d0) then
            this%npts(3)=1
        else
            this%npts(3)=nint(this%l(3)/this%step(3))+1
        end if

        do n=1,3
            allocate(this%gdata(n)%pts(this%npts(n)))
            allocate(this%gdata(n)%wgt(this%npts(n)))
            do i=1,this%npts(n)
                this%gdata(n)%pts(i)=real(i-1)*this%step(n)
                this%gdata(n)%wgt(i)=1.d0
            end do
        end do
    end subroutine

    subroutine copy_grid(this, g2)
        type(grid_t), intent(in) :: this
        type(grid_t), intent(inout) :: g2

        integer(I4) :: i, n

        g2%gauss=this%gauss
        g2%basv=this%basv
        g2%l=this%l
        g2%origin=this%origin
        g2%step=this%step
        g2%npts=this%npts

        do n=1,3
            allocate(g2%gdata(n)%pts(g2%npts(n)))
            allocate(g2%gdata(n)%wgt(g2%npts(n)))
            g2%gdata(n)%pts=this%gdata(n)%pts
            g2%gdata(n)%wgt=this%gdata(n)%wgt
        end do
    end subroutine

    subroutine del_gdata(this)
        type(gdata_t), intent(inout) :: this

        deallocate(this%pts)
        deallocate(this%wgt)
    end subroutine

    subroutine ortho_coordsys(this)
        type(grid_t), intent(inout) :: this

        integer(I4) :: i
        real(DP), dimension(3) :: tvec
        real(DP) :: dpr
        
        dpr=dot_product(this%basv(:,1), this%basv(:,2))
        if (abs(dpr) > DPTOL ) then
            tvec=cross_product(this%basv(:,1), this%basv(:,3))
            tvec=tvec/sqrt(sum(tvec**2))
            this%basv(:,2)=tvec
            call normalise(this%basv)
            call msg_info( 'new_grid():&
                & You specified a nonorthogonal coordinate system.' )
            call nl
            call msg_out('    New unit coordinate system is:')
            call msg_out('    -------------------------------')
            write(str_g, 99) '     v1 = (', this%basv(:,1), ' )'
            call msg_out(str_g)
            write(str_g, 99) '     v2 = (', this%basv(:,2), ' )'
            call msg_out(str_g)
            write(str_g, 99) '     v3 = (', this%basv(:,3), ' )'
            call msg_out(str_g)
            call nl
            return
        end if
99		format(a,3f12.8,a)
    end subroutine
        
    subroutine get_grid_size(this, i, j, k)
        type(grid_t), intent(in) :: this
        integer(I4), intent(out) :: i
        integer(I4), intent(out), optional :: j, k

        i=this%npts(1)
        if (present(j)) j=this%npts(2)
        if (present(k)) k=this%npts(3)
    end subroutine

    subroutine del_grid(this)
        type(grid_t) :: this

        if (this%gtype == 'file') then
            deallocate(this%xdata)
        else
            call del_gdata(this%gdata(1))
            call del_gdata(this%gdata(2))
            call del_gdata(this%gdata(3))
        end if
        call msg_note('Deallocated grid data')
    end subroutine 

    function get_weight(this, i, d) result(w)
        integer(I4), intent(in) :: i, d
        type(grid_t), intent(in) :: this
        real(DP) :: w

        w=this%gdata(d)%wgt(i)
    end function

    function is_gauss_grid(this) result(r)
        type(grid_t), intent(in) :: this
        logical :: r
        
        r=this%gauss
    end function

    function get_grid_length(this) result(l)
        type(grid_t), intent(in) :: this
        real(DP), dimension(3) :: l
        
        l=this%l
    end function

    function gridpoint(this, i, j, k) result(r)
        type(grid_t), intent(in) :: this
        integer(I4), intent(in) :: i, j, k
        real(DP), dimension(3) :: r
        
        if (this%gtype == 'file') then
            r=this%xdata(:,i)
        else
            r=this%origin+&
              this%gdata(1)%pts(i)*this%basv(:,1)+&
              this%gdata(2)%pts(j)*this%basv(:,2)+&
              this%gdata(3)%pts(k)*this%basv(:,3)
        end if
    end function 

    function realpoint(this, i, j) result(r)
        integer(I4), intent(in) :: i, j
        type(grid_t), intent(in) :: this
        real(DP), dimension(3) :: r

        r=this%origin+real(i)*this%step(1)*this%basv(:,1)+&
            real(j)*this%step(2)*this%basv(:,2)
    end function 

    function get_grid_normal(this) result(n)
        type(grid_t), intent(in) :: this
        real(DP), dimension(3) :: n
        
        n=this%basv(:,3)
    end function

    subroutine grid_center(this, center)
        type(grid_t), intent(in) :: this
        real(DP), dimension(3), intent(out) :: center

        real(DP), dimension(3) :: v1, v2
        
        v1=gridpoint(this, this%npts(1), 1, 1)
        v2=gridpoint(this, 1, this%npts(2), 1)
        center=(v1+v2)*0.5d0
!        write(str_g, '(a,3f10.5)') 'Grid center:', center
!        call msg_note(str_g)
!        call nl
    end subroutine

    subroutine extgrid(this)
        type(grid_t), intent(inout) :: this

        integer(I4) :: nlines, i

        this%basv=0.d0
        this%l=0.d0
        this%origin=0.d0
        this%step=0.d0
        this%gtype='file'

        if (is_mpirun) then
            call msg_error('grid type ''file'' does not work with the &
            &parallel version (yet)!')
            call exit(1)
        end if

        if (keyword_is_set(input,'grid.file')) then
            call getkw(input, 'grid.file', str_g)
            open(GRIDFD, file=trim(str_g))
        else
            open(GRIDFD, file='GRIDDATA')
        end if

        nlines=getnlines(GRIDFD)
        allocate(this%xdata(3,nlines))
        this%npts=(/nlines,1,1/)
        read(GRIDFD,*) this%xdata
        close(GRIDFD)
        write(str_g, '(2x,a,i10)') 'Total number of grid points  :', &
            product(this%npts)
        call msg_out(str_g)
        call nl
    end subroutine

    function norm(v) result(n)
        real(DP), dimension(:), intent(in) :: v
        real(DP), dimension(3) :: n
        
        real(DP) :: l

        l=sqrt(sum(v**2))
        n=v/l
    end function

    subroutine plot_grid_xyz(this, fname, mol)
        type(grid_t), intent(inout) :: this
        character(*), intent(in) :: fname
        type(molecule_t) :: mol

        integer(I4) :: natoms, i
        integer(I4) :: p1, p2, p3
        real(DP), dimension(3) :: r, coord
        character(2) :: symbol
        type(atom_t), pointer :: atom
        logical :: show_axis


        natoms=get_natoms(mol)
        
        call get_grid_size(this,p1,p2,p3)
        call getkw(input, 'show_up_axis', show_axis)
        i=0
        if (show_axis) i=1


        write(str_g, '(2a)') 'Grid plot in ', trim(fname)
        call msg_note(str_g)
        open(77,file=trim(fname))

        if (p3 > 1) then
            write(77,*) natoms+8+i
            write(77,*)
        else if (p3 == 1) then
            write(77,*) natoms+4+i
            write(77,*)
        else 
            write(77,*) natoms
            write(77,*)
        end if

        do i=1,natoms
            call get_atom(mol, i, atom)
            call get_symbol(atom, symbol)
            call get_coord(atom, coord)
            write(77,'(a, 3f16.10)') symbol, coord*au2a
        end do

        if (p3 > 1) then
            r=gridpoint(this,1,1,1)
            write(77,'(a,3f16.10)') 'X ', r*au2a
            r=gridpoint(this,p1,1,1)
            write(77,'(a,3f16.10)') 'X ', r*au2a
            r=gridpoint(this,1,p2,1)
            write(77,'(a,3f16.10)') 'X ', r*au2a
            r=gridpoint(this,1,1,p3)
            write(77,'(a,3f16.10)') 'X ', r*au2a
            r=gridpoint(this,p1,p2,1)
            write(77,'(a,3f16.10)') 'X ', r*au2a
            r=gridpoint(this,p1,1,p3)
            write(77,'(a,3f16.10)') 'X ', r*au2a
            r=gridpoint(this,1,p2,p3)
            write(77,'(a,3f16.10)') 'X ', r*au2a
            r=gridpoint(this,p1,p2,p3)
            write(77,'(a,3f16.10)') 'X ', r*au2a
        else if (p3 == 1) then
            r=gridpoint(this,1,1,1)
            write(77,'(a,3f16.10)') 'X ', r*au2a
            r=gridpoint(this,p1,1,1)
            write(77,'(a,3f16.10)') 'X ', r*au2a
            r=gridpoint(this,1,p2,1)
            write(77,'(a,3f16.10)') 'X ', r*au2a
            r=gridpoint(this,p1,p2,1)
            write(77,'(a,3f16.10)') 'X ', r*au2a
        end if

        if (show_axis) then
            select case (trim(this%mode))
                case ('std','base')
                    write(77,'(a,3f16.10)') 'X ', &
                        (this%origin+this%basv(:,3))*au2a
                case ('bond')
                    write(77,'(a,3f16.10)') 'X ', &
                        (this%origin+this%ortho)*au2a
            end select
        end if
        
        close(77)
    end subroutine

    subroutine plot_grid_xyz_old(fname, this, mol, np)
        character(*), intent(in) :: fname
        type(grid_t), intent(inout) :: this
        type(molecule_t) :: mol
        integer(I4), intent(in) :: np

        integer(I4) :: natoms, i, j, k, ir, jr, kr
        integer(I4) :: p1, p2, p3, d1, d2, d3
        real(DP), dimension(3) :: r, coord
        character(2) :: symbol
        type(atom_t), pointer :: atom
        
        natoms=get_natoms(mol)
        
        if (np > 0) then
            call get_grid_size(this,p1,p2,p3)
            if (mod(p1,2) > 0) then
                d1=pltdiv(p1-1,np)
            else
                d1=pltdiv(p1,np)
            end if
            if (mod(p2,2) > 0) then
                d2=pltdiv(p2-1,np)
            else
                d2=pltdiv(p2,np)
            end if
            if (mod(p3,2) > 0) then
                d3=pltdiv(p3-1,np)
            else
                d3=pltdiv(p3,np)
            end if
            i=p1/d1+1
            j=p2/d2+1
            k=p3/d3+1
        else
            i=0; j=0; k=0
            p1=-1; p2=-1; p3=-1;
            d1=1; d2=1; d3=1;
        end if

        write(str_g, '(3a,3i6)') 'Grid plot in ', trim(fname),':', i,j,k
        call msg_note(str_g)
        open(77,file=trim(fname))
        write(77,*) natoms+i*j*k
        write(77,*)
        
        do k=0,p3,d3
            do j=0,p2,d2
                do i=0,p1,d1
                    if (i == 0) then
                        ir=1
                    else 
                        ir=i
                    end if

                    if (j == 0) then
                        jr=1
                    else 
                        jr=j
                    end if

                    if (k == 0) then
                        kr=1
                    else 
                        kr=k
                    end if
                    
                    r=gridpoint(this,jr,ir,kr)
                    write(77,'(a,3f16.10)') 'X ', r*au2a
                end do
            end do
        end do

        do i=1,natoms
            call get_atom(mol, i, atom)
            call get_symbol(atom, symbol)
            call get_coord(atom, coord)
            write(77,'(a, 3f16.10)') symbol, coord*au2a
        end do
        close(77)
    end subroutine

    function pltdiv(p,np) result(d)
        integer(I4), intent(in) :: p, np
        integer(I4) :: d

        integer(I4) :: i, j

        d=1
        if ( p > np) then
            i=1
            j=2
            do while (i > 0)
                i=mod(p,j)
                if ( i == 0) then
                    if (p/j > np) then
                        i=1
                    else
                        d=j
                    end if
                end if
                j=j+1
            end do
        end if
    end function

    subroutine write_grid(this, fd)
        type(grid_t), intent(in) :: this
        integer(I4), intent(in) :: fd

        integer(I4) :: i

        write(fd, *) this%mode
        write(fd, *) this%gauss
        write(fd, *) this%basv
        write(fd, *) this%l
        write(fd, *) this%origin
        write(fd, *) this%step
        write(fd, *) this%npts
        do i=1,3
            write(fd, *) this%gdata(i)%pts
            write(fd, *) this%gdata(i)%wgt
        end do
        write(fd, *)

    end subroutine

    subroutine read_grid(this, fd)
        type(grid_t), intent(inout) :: this
        integer(I4), intent(in) :: fd

        integer(I4) :: i

        read(fd, *) this%gauss
        read(fd, *) this%basv
        read(fd, *) this%l
        read(fd, *) this%origin
        read(fd, *) this%step
        read(fd, *) this%npts
        do i=1,3
            allocate(this%gdata(i)%pts(this%npts(i)))
            allocate(this%gdata(i)%wgt(this%npts(i)))
            read(fd, *) this%gdata(i)%pts
            read(fd, *) this%gdata(i)%wgt
        end do
        write(str_g, '(a,3i5)') 'Number of grid points <v1,v2>:', &
            this%npts(1), this%npts(2), this%npts(3)
        call msg_out(str_g)
        write(str_g, '(a,i7)') 'Total number of grid points  :', &
            product(this%npts)
        call msg_out(str_g)
        call nl
    end subroutine
    
    ! rotate basis vectors (prior to grid setup)
    subroutine rotate(this, angle)
        type(grid_t) :: this
        real(DP), dimension(3), intent(in) :: angle

        real(DP) :: x
        real(DP), dimension(3) :: rad
        real(DP), dimension(3,3) :: euler, rot

        rad=0.d0
        euler=0.d0
        rot=0.d0
        rad=angle/180.d0*PII
        write(str_g, '(a,3f9.5)') 'Rotation is: ', rad
        call msg_info(str_g)

!ol rotations according to right hand rule

! z-mat
        x=rad(3)
        rot(1,1)=cos(x)
        rot(2,2)=cos(x)
        rot(3,3)=1.d0
        !jj rot(1,2)=sin(x)
        !jj rot(2,1)=-sin(x)

        !ol begin
        rot(1,2)=-sin(x)
    rot(2,1)=sin(x)
        !ol end

! y-mat
        x=rad(2)
        euler(1,1)=cos(x)
        euler(2,2)=1.d0
        euler(3,3)=cos(x)
        !jj euler(1,3)=-sin(x)
        !jj euler(3,1)=sin(x)
    
        !ol begin
        euler(1,3)=sin(x)
        euler(3,1)=-sin(x)
        !ol end

        euler=matmul(euler,rot)

! x-mat
        x=rad(1)
        !jj rot(1,1)=1.d0
        !jj rot(2,2)=cos(x)
        !jj rot(3,3)=cos(x)
        !jj rot(1,2)=sin(x)
        !jj rot(2,1)=-sin(x)

        !ol begin
        rot=0.d0
        rot(1,1)=1.d0
        rot(2,2)=cos(x)
        rot(3,3)=cos(x)
        rot(2,3)=-sin(x)
        rot(3,2)=sin(x)
        !ol end


        euler=matmul(rot,euler)

        this%basv=matmul(euler,this%basv)
        this%origin=matmul(euler,this%origin)
    end subroutine

!
! Obsolete stuff
! 
    subroutine hcsmbd(this)
        type(grid_t), intent(inout) :: this

        integer(I4) :: i
        real(DP), dimension(3) :: v1, v2, v3, oo
        real(DP), dimension(2) :: lh, ht
        real(DP) :: l3, r1, r2

        call msg_critical("Sorry, you just stumbled upon a piece of dead &
            &code. Better luck next time.")
        l3=-1.d0
        lh=-1.d0
        ht=-1.d0
        call getkw(input, 'grid.distance', l3)
        call getkw(input, 'grid.width', lh)
        call getkw(input, 'grid.height', ht)
        if ( l3 < 0.d0 ) then
            call msg_critical('grid.distace < 0!')
            stop 
        end if
        if ( sum(lh) < 0.d0 ) then
            call msg_critical('grid.width < 0!')
            stop 
        end if
        if ( sum(ht) < 0.d0 ) then
            call msg_critical('grid.height < 0!')
            stop 
        end if

        i=len(trim(this%mode))
        if (this%mode(1:i) == 'foo') then
            v1=norm(this%basv(:,1)-this%basv(:,2))
            oo=this%basv(:,1)-l3*v1
            v2=norm(oo-this%origin)
            v3=norm(cross_product(v1, v2))
            v1=cross_product(v2, v3)
            this%origin=oo-lh(1)*v2+ht(2)*v3
        else if (this%mode(1:i) == 'bond') then
            v1=norm(this%basv(:,1)-this%basv(:,2))
            oo=this%basv(:,1)-l3*v1
            v3=norm(cross_product(v1, this%origin-oo))
            v2=cross_product(v1, v3)
            this%origin=oo-lh(1)*v2+ht(2)*v3
        else if (this%mode(1:i) == 'bar') then
            v3=norm(this%basv(:,1)-this%basv(:,2))
            oo=this%basv(:,1)-sqrt(sum((this%basv(:,1)-&
            this%basv(:,2))**2))*v3*0.5d0
            v1=norm(cross_product(v3, this%origin-oo))
            v2=norm(cross_product(v1, v3))
            oo=oo+2.d0*v2 !uhh...
            this%origin=oo-lh(1)*v2+ht(2)*v3+l3*v1
        else
            call msg_error('Unknown grid mode: ' // trim(this%mode))
            call exit(1)
        end if

        this%basv(:,1)=this%origin+sum(lh)*v2
        this%basv(:,2)=this%origin-sum(ht)*v3

        r1=sqrt(sum((this%origin+v2)**2))
        r2=sqrt(sum((this%origin-v2)**2))

        call nl
        call msg_out('Grid data')
        call msg_out('------------------------------------------------')
        write(str_g, '(a,3f12.6)') 'v1     ', v1
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') 'v2     ', v2 
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') 'v3     ', v3
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') 'center ', oo
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') 'origin ', this%origin
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') 'basv1  ', this%basv(:,1)
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') 'basv2  ', this%basv(:,2)
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') '|basv1|', &
            sqrt(sum((this%basv(:,1)-this%origin)**2))
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') '|basv2|', &
            sqrt(sum((this%basv(:,2)-this%origin)**2))
        call msg_out(str_g)
        call nl

    end subroutine

    subroutine bondage(this)
        type(grid_t), intent(inout) :: this

        real(DP), dimension(3) :: v1, v2, v3, oo
        real(DP) :: l3, r

        l3=-1.d0
        r=-1.d0
        call getkw(input, 'grid.l3', l3)
        call getkw(input, 'grid.radius', r)
        if ( l3 < 0.d0 ) stop 'grid.l3 < 0!'
        if ( r < 0.d0 ) stop 'grid.radius < 0!'

        v1=norm(this%basv(:,2)-this%basv(:,1))
        v2(1)=-v1(2)-v1(3)
        v2(2)=v1(1)-v1(3)
        v2(3)=v1(1)+v1(2)
        v2=norm(v2)
        v3=cross_product(v1, v2)

        oo=this%basv(:,1)+l3*v1
        this%origin=oo+r*(v2+v3)
        this%basv(:,1)=oo-r*(v2-v3)
        this%basv(:,2)=oo-r*(-v2+v3)
!        print *
!        print *, 'center', oo
!        print *, 'origin', this%origin
!        print *, 'basv1', this%basv(:,1)
!        print *, 'basv2', this%basv(:,2)
!        print *, '|basv1|', sqrt(sum((this%basv(:,1)-this%origin)**2))
!        print *, '|basv2|', sqrt(sum((this%basv(:,2)-this%origin)**2))
!        print *
        
    end subroutine

    subroutine get_basv3(this, i, j, k)
        type(grid_t) :: this
        real(DP), dimension(3) :: i, j, k
        i=this%basv(:,1)
        j=this%basv(:,2)
        k=this%basv(:,3)
    end subroutine

    subroutine get_basv1(this, n, v)
        type(grid_t) :: this
        integer(I4), intent(in) :: n
        real(DP), dimension(3), intent(out) :: v
        v=this%basv(:,n)
    end subroutine

    subroutine get_ortho(this, v)
        type(grid_t) :: this
        real(DP), dimension(3), intent(out) :: v
        v=this%ortho
    end subroutine
    
    subroutine quadrant()
        integer(I4) :: i, j, k
        real(DP), dimension(3) :: v1, v2, v3
        integer, dimension(3,3) :: vtab2, vtab1 = reshape((/ &
        1,0,0, &
        0,1,0, &
        0,0,1/), (/3,3/))

        vtab2=vtab1
        do i=1,3
            do j=i+1,3
                v3=cross_product(dble(vtab1(:,i)), dble(vtab2(:,j)))
                print '(3i4,a,3i4,a,3f6.2)', vtab1(:,i), ':', vtab2(:,j), &
                ' => ', v3
            end do
        end do
        print *
        vtab2=-vtab1
        do i=1,3
            do j=i+1,3
                v3=cross_product(dble(vtab1(:,i)), dble(vtab2(:,j)))
                print '(3i4,a,3i4,a,3f6.2)', vtab1(:,i), ':', vtab2(:,j), &
                ' => ', v3
            end do
        end do
        print *
        vtab1=-vtab1
        do i=1,3
            do j=i+1,3
                v3=cross_product(dble(vtab1(:,i)), dble(vtab2(:,j)))
                print '(3i4,a,3i4,a,3f6.2)', vtab1(:,i), ':', vtab2(:,j), &
                ' => ', v3
            end do
        end do
!        stop

    end subroutine

    function get_grid_range(this,crd) result(r)
        type(grid_t) :: this
        integer(SP), intent(in) :: crd
        real(DP), dimension(2) :: r

        r(1)=this%gdata(crd)%pts(1)
        r(2)=this%gdata(crd)%pts(this%npts(crd))
    end function

end module

!    function gridmap(this, i, j) result(r)
!        type(grid_t), intent(in) :: this
!        integer(I4), intent(in) :: i, j
!        real(DP), dimension(2) :: r

!        real(DP), dimension(2) :: m1, m2
!        real(DP) :: q1, q2, w1, w2

!        if (this%gtype == 'file') then
!            r=0.d0
!            return
!        end if

!        m1=(/0.0, 1.0/)
!        m2=(/1.0, 0.0/)

!        q1=this%gdata(1)%pts(i)
!        q2=this%gdata(2)%pts(j)

!        r=this%map+q1*m1+q2*m2 
!    end function 


! vim:et:sw=4:ts=4
