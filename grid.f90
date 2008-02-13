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
		real(DP), dimension(3) :: origin, ortho
		real(DP), dimension(3) :: step
		integer(I4), dimension(3) :: npts
		type(gdata_t), dimension(3) :: gdata
		character(BUFLEN) :: mode, gtype
		real(DP), dimension(:,:), pointer :: xdata
	end type grid_t

	public grid_t
	public init_grid, del_grid, gridpoint, get_grid_normal
	public get_grid_size, get_weight, write_grid, read_grid
	public get_grid_length, is_lobo_grid, realpoint, copy_grid
	public grid_center, plot_grid_xyz, get_basvec, get_ortho

	interface get_basvec
		module procedure get_basv1
		module procedure get_basv3
	end interface

	private
	real(DP), parameter :: DPTOL=1.d-10
contains

	subroutine init_grid(self, mol)
		type(grid_t) :: self
		type(molecule_t) :: mol

		real(DP) :: ll 
		real(DP), dimension(3) :: angle
		integer(I4), dimension(3) :: ngp
		integer(I4) :: i, j
		
		self%step=1.d0
		self%gtype='even'
		self%basv=0.0
		self%l=0.0
		self%ortho=0.d0
		call getkw(input, 'grid', self%mode)

		! first figure out where and how to place the grid
		select case (trim(self%mode))
			case ('file')
				call extgrid(self)
				return
			case ('std','base')
				call setup_std_grid(self)
			case ('bond')
				call setup_bond_grid(self,mol)
			case default
				call msg_error('Unknown grid type: ' // trim(self%mode))
				stop
		end select

		call msg_out('Grid mode = ' // trim(self%mode))

		call normalise(self%basv)
		call ortho_coordsys(self)
		
		! rotate basis vectors if needed
		if (keyword_is_set(input, 'grid.rotation')) then
			call getkw(input, 'grid.rotation', angle)
			call rotate(self, angle)
		end if

		! calculate distibution of grid points
		call getkw(input, 'grid.type', self%gtype)
		i=len(trim(self%gtype))
		select case (self%gtype(1:i))
			case ('even')
				call setup_even_gdata(self)
			case ('gauss')
				call setup_gauss_gdata(self)
			case default
				call msg_error('Unknown grid type: ' // trim(self%gtype))
				stop
		end select

		write(str_g, '(2x,a,3i5)') 'Number of grid points <v1,v2>:', &
			self%npts(1), self%npts(2), self%npts(3)
		call msg_out(str_g)
		write(str_g, '(2x,a,i10)') 'Total number of grid points  :', &
			product(self%npts)
		call msg_out(str_g)
		call nl
!        call write_grid(self, 50)

	end subroutine

	subroutine setup_std_grid(self)
		type(grid_t) :: self

		integer(I4) :: i
		real(DP), dimension(3) :: normv

		call getkw(input, 'grid.origin', self%origin)
		call getkw(input, 'grid.ivec', self%basv(:,1))
		call getkw(input, 'grid.jvec', self%basv(:,2))
		call getkw(input, 'grid.spacing', self%step)
		call getkw(input, 'grid.lengths', self%l)
! grid_points

		self%basv(:,3)=cross_product(self%basv(:,1),self%basv(:,2)) 
	end subroutine

	subroutine setup_bond_grid(self,mol)
		type(grid_t) :: self
		type(molecule_t) :: mol

		type(atom_t), pointer :: atom

		integer(I4) :: i
		real(DP), dimension(3) :: normv
		integer(I4), dimension(3) :: atoms
		real(DP), dimension(3) :: v1, v2, v3, oo
		real(DP), dimension(2) :: lh, ht
		real(DP) :: l3

		if (keyword_is_set(input, 'grid.atoms')) then
			call getkw(input, 'grid.atoms', atoms)
			call get_atom(mol, atoms(1), atom)
			call get_coord(atom, self%basv(:,1))
			call get_atom(mol, atoms(2), atom)
			call get_coord(atom, self%basv(:,2))
			call get_atom(mol, atoms(3), atom)
			call get_coord(atom, self%origin)
		else
			call getkw(input, 'grid.coord1', self%basv(:,1))
			call getkw(input, 'grid.coord2', self%basv(:,2))
        	call getkw(input, 'grid.coord3', self%origin)
		end if
 		!defaults, etc.

		l3=-1.d0
		lh=-1.d0
		ht=-1.d0
		call getkw(input, 'grid.distance', l3)
		call getkw(input, 'grid.width', lh)
		call getkw(input, 'grid.height', ht)
		self%l=(/sum(lh), sum(ht), 0.d0/)

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

		! figure out "orthogonal" axis for magnetic field
		v1=self%basv(:,1)-self%origin
		v2=self%basv(:,2)-self%origin
		self%ortho=cross_product(v1,v2)
		if (vcmp(self%ortho, NILL_VECTOR)) then
			call msg_error('Basis vectors are linearly dependent, field &
			&direction undefined!')
			stop
		end if
		self%ortho=norm(self%ortho)

		! set up basis vectors for actual grid
		v1=norm(self%basv(:,1)-self%basv(:,2))
		oo=self%basv(:,1)-l3*v1
		v3=norm(cross_product(v1, self%origin-oo))
		v2=cross_product(v1, v3)
		self%origin=oo-lh(1)*v2-ht(2)*v3

        self%basv(:,1)=v2
        self%basv(:,2)=v3
		self%basv(:,3)=cross_product(v2,v3) 

		call nl
		call msg_out('Integration grid data')
		call msg_out('------------------------------------------------')
		write(str_g, '(a,3f12.6)') 'center ', oo
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'origin ', self%origin
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'basv1  ', self%basv(:,1)
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'basv2  ', self%basv(:,2)
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'basv3  ', self%basv(:,3)
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'magnet ', self%ortho
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


	subroutine setup_gauss_gdata(self)
		type(grid_t) :: self

		integer(I4) ::  i, rem, order
		real(DP), dimension(3) :: spc
		logical :: flag

		flag=.false.

		call msg_info('Integration grid selected.')
		call getkw(input, 'grid.gauss_order', order)
		if (keyword_is_set(input,'grid.grid_points')) then
			call getkw(input, 'grid.grid_points', self%npts)
		else
			call getkw(input, 'grid.spacing', spc)
			do i=1,3
				if (abs(spc(i)) < 1.d-10 .or.spc(i) < 0.d0) then
					self%npts(i)=0
				else
					self%npts(i)=nint(self%l(i)/spc(i))
				end if
			end do
		end if
		
		do i=1,3
			if (.not.self%npts(i) > 1) then
				self%npts(i)=0
			end if
			rem=mod(self%npts(i),order)
			if (rem /= 0) then
				self%npts(i)=self%npts(i)-rem+order
				flag=.true.
			end if
		end do
		
		if (flag) then
			write(str_g, '(a,3i5)'), &
			'Adjusted number of grid points for quadrature: ', self%npts
			call msg_info(str_g)
		end if

		self%lobato=.true.
		do i=1,3
			if (self%npts(i) > 0) then
				allocate(self%gdata(i)%pts(self%npts(i)))
				allocate(self%gdata(i)%wgt(self%npts(i)))
				call setup_lobby(0.d0, self%l(i), order, self%gdata(i))
			else
				self%npts(i)=1
				allocate(self%gdata(i)%pts(1))
				allocate(self%gdata(i)%wgt(1))
				call setup_lobby(0.d0, self%l(i), 1, self%gdata(i))
			end if
		end do
	end subroutine

	subroutine setup_even_gdata(self)
		type(grid_t), intent(inout) :: self

		integer(I4) :: i, n

		self%lobato=.false.
		self%npts(1)=nint(self%l(1)/self%step(1))+1
		self%npts(2)=nint(self%l(2)/self%step(2))+1
		if (self%l(3) == 0.d0 .or. self%step(3) == 0.d0) then
			self%npts(3)=1
		else
			self%npts(3)=nint(self%l(3)/self%step(3))+1
		end if

		do n=1,3
			allocate(self%gdata(n)%pts(self%npts(n)))
			allocate(self%gdata(n)%wgt(self%npts(n)))
			do i=1,self%npts(n)
				self%gdata(n)%pts(i)=real(i-1)*self%step(n)
				self%gdata(n)%wgt(i)=1.d0
			end do
		end do
	end subroutine

	subroutine copy_grid(self, g2)
		type(grid_t), intent(in) :: self
		type(grid_t), intent(inout) :: g2

		integer(I4) :: i, n

		g2%lobato=self%lobato
		g2%basv=self%basv
		g2%l=self%l
		g2%origin=self%origin
		g2%step=self%step
		g2%npts=self%npts

		do n=1,3
			allocate(g2%gdata(n)%pts(g2%npts(n)))
			allocate(g2%gdata(n)%wgt(g2%npts(n)))
			g2%gdata(n)%pts=self%gdata(n)%pts
			g2%gdata(n)%wgt=self%gdata(n)%wgt
		end do
	end subroutine

	subroutine del_gdata(self)
		type(gdata_t), intent(inout) :: self

		deallocate(self%pts)
		deallocate(self%wgt)
	end subroutine

	subroutine ortho_coordsys(self)
		type(grid_t), intent(inout) :: self

		integer(I4) :: i
		real(DP), dimension(3) :: tvec
		real(DP) :: dpr
		
		dpr=dot_product(self%basv(:,1), self%basv(:,2))
		if (abs(dpr) > DPTOL ) then
			tvec=cross_product(self%basv(:,1), self%basv(:,3))
			tvec=tvec/sqrt(sum(tvec**2))
			self%basv(:,2)=tvec
			call normalise(self%basv)
			call msg_info( 'init_grid():&
				& You specified a nonorthogonal coordinate system.' )
			call nl
			call msg_out('    New unit coordinate system is:')
			call msg_out('    -------------------------------')
			write(str_g, 99) '     v1 = (', self%basv(:,1), ' )'
			call msg_out(str_g)
			write(str_g, 99) '     v2 = (', self%basv(:,2), ' )'
			call msg_out(str_g)
			write(str_g, 99) '     v3 = (', self%basv(:,3), ' )'
			call msg_out(str_g)
			call nl
			return
		end if
99		format(a,3f12.8,a)
	end subroutine
		
	subroutine get_grid_size(self, i, j, k)
		type(grid_t), intent(in) :: self
		integer(I4), intent(out) :: i
		integer(I4), intent(out), optional :: j, k

		i=self%npts(1)
		if (present(j)) j=self%npts(2)
		if (present(k)) k=self%npts(3)
	end subroutine

	subroutine del_grid(self)
		type(grid_t) :: self

		if (self%gtype == 'file') then
			deallocate(self%xdata)
		else
			call del_gdata(self%gdata(1))
			call del_gdata(self%gdata(2))
			call del_gdata(self%gdata(3))
		end if
		call msg_note('Deallocated grid data')
	end subroutine 

	function get_weight(self, i, d) result(w)
		integer(I4), intent(in) :: i, d
		type(grid_t), intent(in) :: self
		real(DP) :: w

		w=self%gdata(d)%wgt(i)
	end function

	function is_lobo_grid(self) result(r)
		type(grid_t), intent(in) :: self
		logical :: r
		
		r=self%lobato
	end function

	function get_grid_length(self) result(l)
		type(grid_t), intent(in) :: self
		real(DP), dimension(3) :: l
		
		l=self%l
	end function

	function gridpoint(self, i, j, k) result(r)
		type(grid_t), intent(in) :: self
		integer(I4), intent(in) :: i, j, k
		real(DP), dimension(3) :: r
		
		if (self%gtype == 'file') then
			r=self%xdata(:,i)
		else
			r=self%origin+&
			  self%gdata(1)%pts(i)*self%basv(:,1)+&
			  self%gdata(2)%pts(j)*self%basv(:,2)+&
			  self%gdata(3)%pts(k)*self%basv(:,3)
		end if
	end function 

	function realpoint(self, i, j) result(r)
		integer(I4), intent(in) :: i, j
		type(grid_t), intent(in) :: self
		real(DP), dimension(3) :: r

		r=self%origin+real(i)*self%step(1)*self%basv(:,1)+&
			real(j)*self%step(2)*self%basv(:,2)
	end function 

	function get_grid_normal(self) result(n)
		type(grid_t), intent(in) :: self
		real(DP), dimension(3) :: n
		
		n=self%basv(:,3)
	end function

	subroutine grid_center(self, center)
		type(grid_t), intent(in) :: self
		real(DP), dimension(3), intent(out) :: center

		real(DP), dimension(3) :: v1, v2
		
		v1=gridpoint(self, self%npts(1), 1, 1)
		v2=gridpoint(self, 1, self%npts(2), 1)
		center=(v1+v2)*0.5d0
!        write(str_g, '(a,3f10.5)') 'Grid center:', center
!        call msg_note(str_g)
!        call nl
	end subroutine

	subroutine extgrid(self)
		type(grid_t), intent(inout) :: self

		integer(I4) :: nlines, i

		self%basv=0.d0
		self%l=0.d0
		self%origin=0.d0
		self%step=0.d0
		self%gtype='file'

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
		allocate(self%xdata(3,nlines))
		self%npts=(/nlines,1,1/)
		read(GRIDFD,*) self%xdata
		close(GRIDFD)
		write(str_g, '(2x,a,i10)') 'Total number of grid points  :', &
			product(self%npts)
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

	subroutine plot_grid_xyz(self, fname, mol)
		type(grid_t), intent(inout) :: self
		character(*), intent(in) :: fname
		type(molecule_t) :: mol

		integer(I4) :: natoms, i
		integer(I4) :: p1, p2, p3
		real(DP), dimension(3) :: r, coord
		character(2) :: symbol
		type(atom_t), pointer :: atom
		
		natoms=get_natoms(mol)
		
		call get_grid_size(self,p1,p2,p3)

		write(str_g, '(2a)') 'Grid plot in ', trim(fname)
		call msg_note(str_g)
		open(77,file=trim(fname))

		if (p3 > 1) then
			write(77,*) natoms+8
			write(77,*)
			r=gridpoint(self,1,1,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(self,p1,1,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(self,1,p2,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(self,1,1,p3)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(self,p1,p2,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(self,p1,1,p3)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(self,1,p2,p3)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(self,p1,p2,p3)
			write(77,'(a,3f16.10)') 'X ', r*au2a
		else if (p3 == 1) then
			write(77,*) natoms+4
			write(77,*)
			r=gridpoint(self,1,1,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(self,p1,1,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(self,1,p2,1)
			write(77,'(a,3f16.10)') 'X ', r*au2a
			r=gridpoint(self,p1,p2,1)
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

	subroutine plot_grid_xyz_old(fname, self, mol, np)
		character(*), intent(in) :: fname
		type(grid_t), intent(inout) :: self
		type(molecule_t) :: mol
		integer(I4), intent(in) :: np

		integer(I4) :: natoms, i, j, k, ir, jr, kr
		integer(I4) :: p1, p2, p3, d1, d2, d3
		real(DP), dimension(3) :: r, coord
		character(2) :: symbol
		type(atom_t), pointer :: atom
		
		natoms=get_natoms(mol)
		
		if (np > 0) then
			call get_grid_size(self,p1,p2,p3)
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
					
					r=gridpoint(self,jr,ir,kr)
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

	subroutine write_grid(self, fd)
		type(grid_t), intent(in) :: self
		integer(I4), intent(in) :: fd

		integer(I4) :: i

		write(fd, *) self%mode
		write(fd, *) self%lobato
		write(fd, *) self%basv
		write(fd, *) self%l
		write(fd, *) self%origin
		write(fd, *) self%step
		write(fd, *) self%npts
        do i=1,3
			write(fd, *) self%gdata(i)%pts
			write(fd, *) self%gdata(i)%wgt
        end do
		write(fd, *)

	end subroutine

	subroutine read_grid(self, fd)
		type(grid_t), intent(inout) :: self
		integer(I4), intent(in) :: fd

		integer(I4) :: i

		read(fd, *) self%lobato
		read(fd, *) self%basv
		read(fd, *) self%l
		read(fd, *) self%origin
		read(fd, *) self%step
		read(fd, *) self%npts
		do i=1,3
			allocate(self%gdata(i)%pts(self%npts(i)))
			allocate(self%gdata(i)%wgt(self%npts(i)))
			read(fd, *) self%gdata(i)%pts
			read(fd, *) self%gdata(i)%wgt
		end do
		write(str_g, '(a,3i5)') 'Number of grid points <v1,v2>:', &
			self%npts(1), self%npts(2), self%npts(3)
		call msg_out(str_g)
		write(str_g, '(a,i7)') 'Total number of grid points  :', &
			product(self%npts)
		call msg_out(str_g)
		call nl
	end subroutine
	
	! rotate basis vectors (prior to grid setup)
	subroutine rotate(self, angle)
		type(grid_t) :: self
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

! z-mat
		x=rad(3)
		rot(1,1)=cos(x)
		rot(2,2)=cos(x)
		rot(3,3)=1.d0
		rot(1,2)=sin(x)
		rot(2,1)=-sin(x)
! y-mat
		x=rad(2)
		euler(1,1)=cos(x)
		euler(2,2)=1.d0
		euler(3,3)=cos(x)
		euler(1,3)=-sin(x)
		euler(3,1)=sin(x)

		euler=matmul(euler,rot)
! x-mat
		x=rad(1)
		rot(1,1)=1.d0
		rot(2,2)=cos(x)
		rot(3,3)=cos(x)
		rot(1,2)=sin(x)
		rot(2,1)=-sin(x)

		euler=matmul(rot,euler)

		self%basv=matmul(euler,self%basv)
		self%origin=matmul(euler,self%origin)
	end subroutine

!
! Obsolete stuff
! 
	subroutine hcsmbd(self)
		type(grid_t), intent(inout) :: self

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

		i=len(trim(self%mode))
		if (self%mode(1:i) == 'foo') then
			v1=norm(self%basv(:,1)-self%basv(:,2))
			oo=self%basv(:,1)-l3*v1
			v2=norm(oo-self%origin)
			v3=norm(cross_product(v1, v2))
			v1=cross_product(v2, v3)
			self%origin=oo-lh(1)*v2+ht(2)*v3
		else if (self%mode(1:i) == 'bond') then
			v1=norm(self%basv(:,1)-self%basv(:,2))
			oo=self%basv(:,1)-l3*v1
			v3=norm(cross_product(v1, self%origin-oo))
			v2=cross_product(v1, v3)
			self%origin=oo-lh(1)*v2+ht(2)*v3
		else if (self%mode(1:i) == 'bar') then
			v3=norm(self%basv(:,1)-self%basv(:,2))
			oo=self%basv(:,1)-sqrt(sum((self%basv(:,1)-&
			self%basv(:,2))**2))*v3*0.5d0
			v1=norm(cross_product(v3, self%origin-oo))
			v2=norm(cross_product(v1, v3))
			oo=oo+2.d0*v2 !uhh...
			self%origin=oo-lh(1)*v2+ht(2)*v3+l3*v1
		else
			call msg_error('Unknown grid mode: ' // trim(self%mode))
			call exit(1)
		end if

        self%basv(:,1)=self%origin+sum(lh)*v2
        self%basv(:,2)=self%origin-sum(ht)*v3

		r1=sqrt(sum((self%origin+v2)**2))
		r2=sqrt(sum((self%origin-v2)**2))

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
		write(str_g, '(a,3f12.6)') 'origin ', self%origin
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'basv1  ', self%basv(:,1)
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') 'basv2  ', self%basv(:,2)
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') '|basv1|', &
			sqrt(sum((self%basv(:,1)-self%origin)**2))
		call msg_out(str_g)
		write(str_g, '(a,3f12.6)') '|basv2|', &
			sqrt(sum((self%basv(:,2)-self%origin)**2))
		call msg_out(str_g)
		call nl

	end subroutine

	subroutine bondage(self)
		type(grid_t), intent(inout) :: self

		real(DP), dimension(3) :: v1, v2, v3, oo
		real(DP) :: l3, r

		l3=-1.d0
		r=-1.d0
		call getkw(input, 'grid.l3', l3)
		call getkw(input, 'grid.radius', r)
		if ( l3 < 0.d0 ) stop 'grid.l3 < 0!'
		if ( r < 0.d0 ) stop 'grid.radius < 0!'

		v1=norm(self%basv(:,2)-self%basv(:,1))
		v2(1)=-v1(2)-v1(3)
		v2(2)=v1(1)-v1(3)
		v2(3)=v1(1)+v1(2)
		v2=norm(v2)
		v3=cross_product(v1, v2)

		oo=self%basv(:,1)+l3*v1
		self%origin=oo+r*(v2+v3)
        self%basv(:,1)=oo-r*(v2-v3)
        self%basv(:,2)=oo-r*(-v2+v3)
!        print *
!        print *, 'center', oo
!        print *, 'origin', self%origin
!        print *, 'basv1', self%basv(:,1)
!        print *, 'basv2', self%basv(:,2)
!        print *, '|basv1|', sqrt(sum((self%basv(:,1)-self%origin)**2))
!        print *, '|basv2|', sqrt(sum((self%basv(:,2)-self%origin)**2))
!        print *
		
	end subroutine

	subroutine get_basv3(self, i, j, k)
		type(grid_t) :: self
		real(DP), dimension(3) :: i, j, k
		i=self%basv(:,1)
		j=self%basv(:,2)
		k=self%basv(:,3)
	end subroutine

	subroutine get_basv1(self, n, v)
		type(grid_t) :: self
		integer(I4), intent(in) :: n
		real(DP), dimension(3), intent(out) :: v
		v=self%basv(:,n)
	end subroutine

	subroutine get_ortho(self, v)
		type(grid_t) :: self
		real(DP), dimension(3), intent(out) :: v
		v=self%ortho
	end subroutine
end module

!    function gridmap(self, i, j) result(r)
!        type(grid_t), intent(in) :: self
!        integer(I4), intent(in) :: i, j
!        real(DP), dimension(2) :: r

!        real(DP), dimension(2) :: m1, m2
!        real(DP) :: q1, q2, w1, w2

!        if (self%gtype == 'file') then
!            r=0.d0
!            return
!        end if

!        m1=(/0.0, 1.0/)
!        m2=(/1.0, 0.0/)

!        q1=self%gdata(1)%pts(i)
!        q2=self%gdata(2)%pts(j)

!        r=self%map+q1*m1+q2*m2 
!    end function 

