! Setup grids for all possible purposes.
! Written by Barnabas Plomlund.
!
! This module is bit of a mess...
! 

module grid_class
	use globals_m
	use tensor_m
	use gaussint_m
	use basis_class
	use teletype_m
	implicit none

	type grid_t
		logical :: lobato  ! integration grid, 1,2 or 3-D
		real(DP), dimension(3,3) :: basv   ! grid basis vectors
		real(DP), dimension(3) :: l        ! |v|
		real(DP), dimension(3) :: origin
		real(DP), dimension(3) :: step
		real(DP), dimension(2) :: map
		integer(I4), dimension(3) :: npts
		type(gdata_t), dimension(3) :: gdata
		character(BUFLEN) :: mode, gtype
		real(DP), dimension(:,:), pointer :: xdata
	end type grid_t
	
	public init_grid, del_grid, gridpoint, gridmap, get_grid_normal
	public grid_t, get_grid_size, get_weight, write_grid, read_grid
	public get_grid_length, is_lobo_grid, realpoint, copy_grid
	public grid_center, plot_grid_xyz

	private
	real(DP), parameter :: DPTOL=1.d-10
contains

	subroutine init_grid(g, mol)
		type(grid_t) :: g
		type(molecule_t) :: mol

		real(DP) :: ll 
		real(DP), dimension(3) :: angle
		integer(I4), dimension(3) :: ngp
		integer(I4) :: i, j
		
		g%step=1.d0
		g%gtype='even'
		g%map=0.0
		g%basv=0.0
		g%l=0.0
		call getkw(input, 'grid', g%mode)

		! first figure out where and how to place the grid
		select case (trim(g%mode))
			case ('file')
				call extgrid(g)
				return
			case ('std','base')
				call setup_std_grid(g)
			case ('bond')
				call setup_bond_grid(g,mol)
			case default
				call msg_error('Unknown grid type: ' // trim(g%mode))
				stop
		end select

		call msg_out('Grid mode = ' // trim(g%mode))

		if (keyword_is_set(input, 'grid.map')) then
			call getkw(input, 'grid.map', g%map)
		end if

		call normalise(g%basv)
		call ortho_coordsys(g)
		call check_handedness(g)
		
		! rotate basis vectors if needed
		if (keyword_is_set(input, 'grid.angle')) then
			call getkw(input, 'grid.angle', angle)
			call rotate(g, angle)
		end if

		! calculate distibution of grid points
		call getkw(input, 'grid.type', g%gtype)
		i=len(trim(g%gtype))
		select case (g%gtype(1:i))
			case ('even')
				call setup_even_gdata(g)
			case ('gauss')
				call setup_gauss_gdata(g)
			case default
				call msg_error('Unknown grid type: ' // trim(g%gtype))
				stop
		end select

		write(str_g, '(2x,a,3i5)') 'Number of grid points <v1,v2>:', &
			g%npts(1), g%npts(2), g%npts(3)
		call msg_out(str_g)
		write(str_g, '(2x,a,i10)') 'Total number of grid points  :', &
			product(g%npts)
		call msg_out(str_g)
		call nl

	end subroutine

	subroutine setup_std_grid(g)
		type(grid_t) :: g

		integer(I4) :: i
		real(DP), dimension(3) :: normv

		call getkw(input, 'grid.origin', g%origin)
		call getkw(input, 'grid.ivec', g%basv(:,1))
		call getkw(input, 'grid.jvec', g%basv(:,2))
		call getkw(input, 'grid.kvec', g%basv(:,3))
		call getkw(input, 'grid.spacing', g%step)
		call getkw(input, 'grid.lengths', g%l)

		g%basv(:,1)=g%basv(:,1)
		g%basv(:,2)=g%basv(:,2)
		g%basv(:,3)=cross_product(g%basv(:,1),g%basv(:,2)) 
	end subroutine

	subroutine setup_bond_grid(g,mol)
		type(grid_t) :: g
		type(molecule_t) :: mol

		type(atom_t), pointer :: atom

		integer(I4) :: i
		real(DP), dimension(3) :: normv
		integer(I4), dimension(2) :: atoms
		real(DP), dimension(3) :: v1, v2, v3, oo
		real(DP), dimension(2) :: lh, ht
		real(DP) :: l3

		call getkw(input, 'grid.origin', g%origin)

		if (keyword_is_set(input, 'grid.atoms')) then
			call getkw(input, 'grid.atoms', atoms)
			call get_atom(mol, atoms(1), atom)
			call get_coord(atom, g%basv(:,1))
			call get_atom(mol, atoms(2), atom)
			call get_coord(atom, g%basv(:,2))
		else
			call getkw(input, 'grid.coord1', g%basv(:,1))
			call getkw(input, 'grid.coord2', g%basv(:,2))
		end if
 		!defaults, etc.
		call getkw(input, 'grid.spacing', g%step)

		l3=-1.d0
		lh=-1.d0
		ht=-1.d0
		call getkw(input, 'grid.distance', l3)
		call getkw(input, 'grid.width', lh)
		call getkw(input, 'grid.height', ht)
		g%l=(/sum(lh), sum(ht), 0.d0/)

		if ( l3 < 0.d0 ) then
			call msg_critical('grid.distance < 0!')
			stop 
		end if
		if ( sum(lh) < 0.d0 ) then
			call msg_critical('grid.width < 0!')
			stop 
		end if
		if ( sum(ht) < 0.d0 ) then
			call msg_critical('grid.heigth < 0!')
			stop 
		end if

		v1=norm(g%basv(:,1)-g%basv(:,2))
		oo=g%basv(:,1)-l3*v1
		v3=norm(cross_product(v1, g%origin-oo))
		v2=cross_product(v1, v3)
		g%origin=oo-lh(1)*v2-ht(2)*v3

        g%basv(:,1)=v2
        g%basv(:,2)=v3
		g%basv(:,3)=cross_product(v2,v3) 

		call nl
		call msg_out('Grid data')
		call msg_out('------------------------------------------------')
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'center ', oo
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'origin ', g%origin
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'basv1  ', g%basv(:,1)
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'basv2  ', g%basv(:,2)
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'basv3  ', g%basv(:,3)
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

	subroutine setup_gauss_gdata(g)
		type(grid_t) :: g

		integer(I4), dimension(3) :: ngp
		integer(I4) :: i

		call msg_info('Integration grid selected.')
		call nl
		call getkw(input, 'grid.gauss_points', ngp)
		call getkw(input, 'grid.grid_points', g%npts)

		g%lobato=.true.
		do i=1,3
			if (g%npts(i) > 0) then
				allocate(g%gdata(i)%pts(g%npts(i)))
				allocate(g%gdata(i)%wgt(g%npts(i)))
			else
				ngp(i)=1
				g%npts(i)=1
				allocate(g%gdata(i)%pts(1))
				allocate(g%gdata(i)%wgt(1))
			end if
			call setup_lobby(0.d0, g%l(i), ngp(i), g%gdata(i))
		end do
	end subroutine

	subroutine check_handedness(g) 
		type(grid_t) :: g

		real(DP), dimension(3) :: magnet
		real(DP) :: x

		! bugger... we can be in "wrong" sect, so we need to ensure we are at
		! the top
		call push_section(input, '.')
		call getkw(input, 'cdens.magnet', magnet) 
		call pop_section(input)

		x=dot_product(g%basv(:,3), magnet)
		if (x < 0.d0) then
			call nl
			call msg_warn('*****************************')
			call msg_warn('Left handed coordinate system')
			call msg_warn('*****************************')
			call nl
		end if
		if (abs(x) /= 1.d0 .and. abs(x) > 1.d-10) then
			call msg_info('Magnetic field not orthogonal to grid')
			print *, x
			call nl
		end if
	end subroutine

	subroutine setup_even_gdata(g)
		type(grid_t), intent(inout) :: g

		integer(I4) :: i, n

		g%lobato=.false.
		g%npts(1)=nint(g%l(1)/g%step(1))+1
		g%npts(2)=nint(g%l(2)/g%step(2))+1
		if (g%l(3) == 0.d0 .or. g%step(3) == 0.d0) then
			g%npts(3)=1
		else
			g%npts(3)=nint(g%l(3)/g%step(3))+1
		end if

		do n=1,3
			allocate(g%gdata(n)%pts(g%npts(n)))
			allocate(g%gdata(n)%wgt(g%npts(n)))
			do i=1,g%npts(n)
				g%gdata(n)%pts(i)=real(i-1)*g%step(n)
				g%gdata(n)%wgt(i)=1.d0
			end do
		end do
	end subroutine

	subroutine copy_grid(g, g2)
		type(grid_t), intent(in) :: g
		type(grid_t), intent(inout) :: g2

		integer(I4) :: i, n

		g2%lobato=g%lobato
		g2%basv=g%basv
		g2%l=g%l
		g2%origin=g%origin
		g2%step=g%step
		g2%npts=g%npts

		do n=1,3
			allocate(g2%gdata(n)%pts(g2%npts(n)))
			allocate(g2%gdata(n)%wgt(g2%npts(n)))
			g2%gdata(n)%pts=g%gdata(n)%pts
			g2%gdata(n)%wgt=g%gdata(n)%wgt
		end do
	end subroutine

	subroutine del_gdata(g)
		type(gdata_t), intent(inout) :: g

		deallocate(g%pts)
		deallocate(g%wgt)
	end subroutine

	subroutine ortho_coordsys(g)
		type(grid_t), intent(inout) :: g

		integer(I4) :: i
		real(DP), dimension(3) :: tvec
		real(DP) :: dpr
		
		dpr=dot_product(g%basv(:,1), g%basv(:,2))
		if (abs(dpr) > DPTOL ) then
			tvec=cross_product(g%basv(:,1), g%basv(:,3))
			tvec=tvec/sqrt(sum(tvec**2))
			g%basv(:,2)=tvec
			call normalise(g%basv)
			call msg_info( 'init_grid():&
				& You specified a nonorthogonal coordinate system.' )
			call nl
			call msg_out('    New unit coordinate system is:')
			call msg_out('    -------------------------------')
			write(str_g, 99) '     v1 = (', g%basv(:,1), ' )'
			call msg_out(str_g)
			write(str_g, 99) '     v2 = (', g%basv(:,2), ' )'
			call msg_out(str_g)
			write(str_g, 99) '     v3 = (', g%basv(:,3), ' )'
			call msg_out(str_g)
			call nl
			return
		end if
99		format(a,3f12.8,a)
	end subroutine
		
	subroutine get_grid_size(g, i, j, k)
		type(grid_t), intent(in) :: g
		integer(I4), intent(out) :: i
		integer(I4), intent(out), optional :: j, k

		i=g%npts(1)
		if (present(j)) j=g%npts(2)
		if (present(k)) k=g%npts(3)
	end subroutine

	subroutine del_grid(g)
		type(grid_t) :: g

		if (g%gtype == 'file') then
			deallocate(g%xdata)
		else
			call del_gdata(g%gdata(1))
			call del_gdata(g%gdata(2))
			call del_gdata(g%gdata(3))
		end if
		call msg_note('Deallocated grid data')
	end subroutine 

	function get_weight(g, i, d) result(w)
		integer(I4), intent(in) :: i, d
		type(grid_t), intent(in) :: g
		real(DP) :: w

		w=g%gdata(d)%wgt(i)
	end function

	function is_lobo_grid(g) result(r)
		type(grid_t), intent(in) :: g
		logical :: r
		
		r=g%lobato
	end function

	function get_grid_length(g) result(l)
		type(grid_t), intent(in) :: g
		real(DP), dimension(3) :: l
		
		l=g%l
	end function

	function gridpoint(g, i, j, k) result(r)
		type(grid_t), intent(in) :: g
		integer(I4), intent(in) :: i, j, k
		real(DP), dimension(3) :: r
		
		real(DP) :: q1, q2, q3

		if (g%gtype == 'file') then
			r=g%xdata(:,i)
		else
			r=g%origin+&
			  g%gdata(1)%pts(i)*g%basv(:,1)+&
			  g%gdata(2)%pts(j)*g%basv(:,2)+&
			  g%gdata(3)%pts(k)*g%basv(:,3)
		end if
	end function 

	function realpoint(g, i, j) result(r)
		integer(I4), intent(in) :: i, j
		type(grid_t), intent(in) :: g
		real(DP), dimension(3) :: r

		r=g%origin+real(i)*g%step(1)*g%basv(:,1)+real(j)*g%step(2)*g%basv(:,2)
	end function 

	function gridmap(g, i, j) result(r)
		type(grid_t), intent(in) :: g
		integer(I4), intent(in) :: i, j
		real(DP), dimension(2) :: r

		real(DP), dimension(2) :: m1, m2
		real(DP) :: q1, q2, w1, w2

		if (g%gtype == 'file') then
			r=0.d0
			return
		end if

		m1=(/0.0, 1.0/)
		m2=(/1.0, 0.0/)

		q1=g%gdata(1)%pts(i)
		q2=g%gdata(2)%pts(j)

		r=g%map+q1*m1+q2*m2 
	end function 

	function get_grid_normal(g) result(n)
		type(grid_t), intent(in) :: g
		real(DP), dimension(3) :: n
		
		n=g%basv(:,3)
	end function

	subroutine grid_center(grid, center)
		type(grid_t), intent(in) :: grid
		real(DP), dimension(3), intent(out) :: center

		real(DP), dimension(3) :: v1, v2
		
		v1=gridpoint(grid, grid%npts(1), 1, 1)
		v2=gridpoint(grid, 1, grid%npts(2), 1)
		center=(v1+v2)*0.5d0
!        write(str_g, '(a,3f10.5)') 'Grid center:', center
!        call msg_note(str_g)
!        call nl
	end subroutine

	subroutine extgrid(g)
		type(grid_t), intent(inout) :: g

		integer(I4) :: nlines, i

		g%basv=0.d0
		g%l=0.d0
		g%origin=0.d0
		g%step=0.d0
		g%map=0.d0
		g%gtype='file'

		if (mpirun_p) then
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
		allocate(g%xdata(3,nlines))
		g%npts=(/nlines,1,1/)
		read(GRIDFD,*) g%xdata
		close(GRIDFD)
		write(str_g, '(2x,a,i10)') 'Total number of grid points  :', &
			product(g%npts)
		call msg_out(str_g)
		call nl
	end subroutine

	function norm(v) result(n)
		real(DP), dimension(3), intent(in) :: v
		real(DP), dimension(3) :: n
		
		real(DP) :: l

		l=sqrt(sum(v**2))
		n=v/l
	end function

	subroutine plot_grid_xyz(fname, g, mol, np)
		character(*), intent(in) :: fname
		type(grid_t), intent(inout) :: g
		type(molecule_t) :: mol
		integer(I4), intent(in) :: np

		integer(I4) :: natoms, i
		integer(I4) :: p1, p2, p3
		real(DP), dimension(3) :: r, coord
		character(2) :: symbol
		type(atom_t), pointer :: atom
		
		natoms=get_natoms(mol)
		
		p3=0
		if (np > 0) then
			call get_grid_size(g,p1,p2,p3)
		end if

		write(str_g, '(2a)') 'Grid plot in ', trim(fname)
		call msg_note(str_g)
		open(77,file=trim(fname))

		if (p3 > 1) then
			write(77,*) natoms+8
			write(77,*)
			r=gridpoint(g,1,1,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(g,p1,1,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(g,1,p2,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(g,1,1,p3)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(g,p1,p2,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(g,p1,1,p3)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(g,1,p2,p3)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(g,p1,p2,p3)
			write(77,'(a,3f16.10)') 'X ', r*au2a
		else if (p3 == 1) then
			write(77,*) natoms+4
			write(77,*)
			r=gridpoint(g,1,1,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(g,p1,1,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(g,1,p2,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(g,p1,p2,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
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
		close(77)
	end subroutine

	subroutine plot_grid_xyz_old(fname, g, mol, np)
		character(*), intent(in) :: fname
		type(grid_t), intent(inout) :: g
		type(molecule_t) :: mol
		integer(I4), intent(in) :: np

		integer(I4) :: natoms, i, j, k, ir, jr, kr
		integer(I4) :: p1, p2, p3, d1, d2, d3
		real(DP), dimension(3) :: r, coord
		character(2) :: symbol
		type(atom_t), pointer :: atom
		
		natoms=get_natoms(mol)
		
		if (np > 0) then
			call get_grid_size(g,p1,p2,p3)
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
					
					r=gridpoint(g,jr,ir,kr)
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
		else 
			d=1
		end if
	end function

	subroutine write_grid(g, fd)
		type(grid_t), intent(in) :: g
		integer(I4), intent(in) :: fd

		integer(I4) :: i

		write(fd, *) g%lobato
		write(fd, *) g%basv
		write(fd, *) g%l
		write(fd, *) g%origin
		write(fd, *) g%step
		write(fd, *) g%npts
        do i=1,3
			write(fd, *) g%gdata(i)%pts
			write(fd, *) g%gdata(i)%wgt
        end do
		write(fd, *)
	end subroutine

	subroutine read_grid(g, fd)
		type(grid_t), intent(inout) :: g
		integer(I4), intent(in) :: fd

		integer(I4) :: i

		call getkw(input, 'grid.map', g%map)

		read(fd, *) g%lobato
		read(fd, *) g%basv
		read(fd, *) g%l
		read(fd, *) g%origin
		read(fd, *) g%step
		read(fd, *) g%npts
		do i=1,3
			allocate(g%gdata(i)%pts(g%npts(i)))
			allocate(g%gdata(i)%wgt(g%npts(i)))
			read(fd, *) g%gdata(i)%pts
			read(fd, *) g%gdata(i)%wgt
		end do
		write(str_g, '(a,3i5)') 'Number of grid points <v1,v2>:', &
			g%npts(1), g%npts(2), g%npts(3)
		call msg_out(str_g)
		write(str_g, '(a,i7)') 'Total number of grid points  :', &
			product(g%npts)
		call msg_out(str_g)
		call nl
	end subroutine
	
	! rotate basis vectors (prior to grid setup)
	subroutine rotate(g, angle)
		type(grid_t) :: g
		real(DP), dimension(3), intent(in) :: angle

		real(DP), dimension(3) :: rad
		real(DP), dimension(3,3) :: euler, rot

		rad=angle/180.d0*PII
		write(str_g, '(a,3f9.5)') 'Rotation is: ', rad
		call msg_info(str_g)

! z-mat
		rot(1,1)=cos(rad(3))
		rot(2,2)=cos(rad(3))
		rot(3,3)=1
		rot(1,2)=sin(rad(3))
		rot(2,1)=-sin(rad(3))
! y-mat
		euler(1,1)=cos(rad(2))
		euler(2,2)=1
		euler(3,3)=cos(rad(2))
		euler(1,3)=-sin(rad(2))
		euler(3,1)=sin(rad(2))

		euler=matmul(euler,rot)
! x-mat
		rot(1,1)=1
		rot(2,2)=cos(rad(1))
		rot(3,3)=cos(rad(1))
		rot(1,2)=sin(rad(1))
		rot(2,1)=-sin(rad(1))

		euler=matmul(rot,euler)

		g%basv=matmul(euler,g%basv)
		g%origin=matmul(euler,g%origin)
	end subroutine

!
! Obsolete stuff
! 
	subroutine hcsmbd(g)
		type(grid_t), intent(inout) :: g

		integer(I4) :: i
		real(DP), dimension(3) :: v1, v2, v3, oo
		real(DP), dimension(2) :: lh, ht
		real(DP) :: l3, r1, r2

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

		i=len(trim(g%mode))
		if (g%mode(1:i) == 'foo') then
			v1=norm(g%basv(:,1)-g%basv(:,2))
			oo=g%basv(:,1)-l3*v1
			v2=norm(oo-g%origin)
			v3=norm(cross_product(v1, v2))
			v1=cross_product(v2, v3)
			g%origin=oo-lh(1)*v2+ht(2)*v3
		else if (g%mode(1:i) == 'bond') then
			v1=norm(g%basv(:,1)-g%basv(:,2))
			oo=g%basv(:,1)-l3*v1
			v3=norm(cross_product(v1, g%origin-oo))
			v2=cross_product(v1, v3)
			g%origin=oo-lh(1)*v2+ht(2)*v3
		else if (g%mode(1:i) == 'bar') then
			v3=norm(g%basv(:,1)-g%basv(:,2))
			oo=g%basv(:,1)-sqrt(sum((g%basv(:,1)-&
			g%basv(:,2))**2))*v3*0.5d0
			v1=norm(cross_product(v3, g%origin-oo))
			v2=norm(cross_product(v1, v3))
			oo=oo+2.d0*v2 !uhh...
			g%origin=oo-lh(1)*v2+ht(2)*v3+l3*v1
		else
			call msg_error('Unknown grid mode: ' // trim(g%mode))
			call exit(1)
		end if

        g%basv(:,1)=g%origin+sum(lh)*v2
        g%basv(:,2)=g%origin-sum(ht)*v3

		r1=sqrt(sum((g%origin+v2)**2))
		r2=sqrt(sum((g%origin-v2)**2))

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
		write(str_g, '(a,3f12.6)') 'origin ', g%origin
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'basv1  ', g%basv(:,1)
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'basv2  ', g%basv(:,2)
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') '|basv1|', &
			sqrt(sum((g%basv(:,1)-g%origin)**2))
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') '|basv2|', &
			sqrt(sum((g%basv(:,2)-g%origin)**2))
		call msg_out(str_g)
		call nl

	end subroutine

	subroutine bondage(g)
		type(grid_t), intent(inout) :: g

		real(DP), dimension(3) :: v1, v2, v3, oo
		real(DP) :: l3, r

		l3=-1.d0
		r=-1.d0
		call getkw(input, 'grid.l3', l3)
		call getkw(input, 'grid.radius', r)
		if ( l3 < 0.d0 ) stop 'grid.l3 < 0!'
		if ( r < 0.d0 ) stop 'grid.radius < 0!'

		v1=norm(g%basv(:,2)-g%basv(:,1))
		v2(1)=-v1(2)-v1(3)
		v2(2)=v1(1)-v1(3)
		v2(3)=v1(1)+v1(2)
		v2=norm(v2)
		v3=cross_product(v1, v2)

		oo=g%basv(:,1)+l3*v1
		g%origin=oo+r*(v2+v3)
        g%basv(:,1)=oo-r*(v2-v3)
        g%basv(:,2)=oo-r*(-v2+v3)
!        print *
!        print *, 'center', oo
!        print *, 'origin', g%origin
!        print *, 'basv1', g%basv(:,1)
!        print *, 'basv2', g%basv(:,2)
!        print *, '|basv1|', sqrt(sum((g%basv(:,1)-g%origin)**2))
!        print *, '|basv2|', sqrt(sum((g%basv(:,2)-g%origin)**2))
!        print *
		
	end subroutine
end module

