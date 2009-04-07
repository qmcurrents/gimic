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
	
	public init_edens, del_edens, edens_t
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
	subroutine init_edens(self, mol, dens, grid)
		type(edens_t) :: self
		type(molecule_t) :: mol
		type(dens_t), target :: dens
		type(grid_t), target :: grid

		integer(I4) :: n
		integer(I4) :: p1, p2

		call get_grid_size(grid, p1, p2)
		call get_dens(dens, self%aodens)
		n=get_ncgto(mol)
		allocate(self%tmp(n))
		allocate(self%buf(p1,p2))
		call init_bfeval(self%bf, mol)
		self%grid=>grid

		if (master_p) then
			call getkw(input, 'edens.density', str_g)
			open(EDFD, file=trim(str_g), access='direct', recl=p1*p2*DP)
		end if
	end subroutine

	subroutine del_edens(self)
		type(edens_t) :: self

		if (associated(self%tmp)) deallocate(self%tmp)
		if (associated(self%buf)) deallocate(self%buf)
		call del_bfeval(self%bf)
		nullify(self%aodens,self%grid)
		if (master_p) then
			close(EDFD)
		end if
	end subroutine

	subroutine set_edens(self, k)
		type(edens_t), intent(in) :: self
		integer(I4), intent(in) :: k

		write(EDFD, rec=k) self%buf
	end subroutine

	subroutine edens_plot(self)
		type(edens_t), intent(inout) :: self
		
		integer(I4) :: i,j,p1,p2,p3
		real(DP) :: amax
		real(DP), dimension(3) :: rr
		real(DP), dimension(:,:), pointer :: buf

		call get_grid_size(self%grid, p1, p2,p3)
		buf=>self%buf

		amax=D0
		read(EDFD, rec=1) self%buf
		call getkw(input, 'edens.density_plot', str_g)
		open(EDPFD, file=trim(str_g))
		do j=1,p2
			do i=1,p1
				rr=gridpoint(self%grid, i, j, 1)
				write(EDPFD, '(4f19.8)') rr, buf(i,j)
				if (abs(buf(i,j)) > amax) amax=abs(buf(i,j))
			end do
			write(EDPFD, *) 
		end do
		close(EDPFD)
		write(str_g, '(a,e19.12)') 'Max electronic density:', amax
		call msg_info(str_g)
		call edens_gopenmol(self)
		call edens_cube(self)
	end subroutine

	subroutine edens_direct(self, k)
		type(edens_t) :: self
		integer(I4), intent(in) :: k

		integer(I4) :: i, j, p1, p2
		integer(I4) :: lo, hi
		real(DP), dimension(3) :: rr
		real(DP), dimension(:,:), pointer :: buf
		real(DP), dimension(:), pointer :: bfvec

		call get_grid_size(self%grid, p1, p2)
		call schedule(p2, lo, hi)

		buf=>self%buf

		do j=lo,hi
			do i=1,p1
				rr=gridpoint(self%grid, i, j, k)
				call bfeval(self%bf,rr, bfvec)
				self%tmp=matmul(self%aodens, bfvec)
				self%buf(i,j)=dot_product(self%tmp, bfvec)
			end do
		end do
		call gather_data(self%buf, self%buf(:,lo:hi))
	end subroutine

	subroutine edens(self)
		type(edens_t) :: self

		integer(I4) :: i, j, k, p1, p2, p3
		integer(I4) :: lo, hi
		real(DP), dimension(3) :: rr
		real(DP), dimension(:,:), pointer :: buf
		real(DP), dimension(:), pointer :: bfvec

		call get_grid_size(self%grid, p1, p2, p3)
		call schedule(p2, lo, hi)

		buf=>self%buf
	
		do k=1,p3
			do j=lo, hi
				do i=1,p1
					rr=gridpoint(self%grid, i, j, k)
					call bfeval(self%bf,rr, bfvec)
					self%tmp=matmul(self%aodens, bfvec)
					self%buf(i,j)=dot_product(self%tmp, bfvec)
				end do
			end do
			call gather_data(self%buf, self%buf(:,lo:hi))
			if (master_p) write(EDFD, rec=k) self%buf
		end do
	end subroutine

	subroutine edens_gopenmol(self)
		type(edens_t) :: self

		integer(I4) :: surface, rank, p1, p2, p3
		integer(I4) :: i, j, k, l
		real(SP), dimension(3) :: qmin, qmax
		real(DP), dimension(:,:), pointer :: buf
		character(BUFLEN) :: gopen_file

		buf=>self%buf
		gopen_file=''
		call getkw(input, 'edens.gopenmol', gopen_file)
		if (trim(gopen_file) == '') return
		open(GOPFD,file=trim(gopen_file),access='direct',recl=I4)

		surface=200
		rank=3

		call get_grid_size(self%grid, p1, p2, p3)
		qmin=real(gridpoint(self%grid,1,1,1)*AU2A)
		qmax=real(gridpoint(self%grid,p1,p2,p3)*AU2A)

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

		l=12
		do k=1,p3
			read(EDFD, rec=k) buf
			do j=1,p2
				do i=1,p1
					write(GOPFD,rec=l) real(buf(i,j))
					l=l+1
				end do
			end do
		end do

		close(GOPFD)
	end subroutine

	subroutine edens_cube(self)
		type(edens_t) :: self

		integer(I4) :: p1, p2, p3
		integer(I4) ::  k
		real(DP), dimension(:,:,:), allocatable :: buf
		character(BUFLEN) :: fname

		fname=''
		call getkw(input, 'edens.cube',fname)
		if (trim(fname) == '') return

		call get_grid_size(self%grid, p1, p2, p3)

		allocate(buf(p1,p2,p3))
		do k=1,p3
			read(EDFD, rec=k) buf(:,:,k)
		end do
		call write_cubeplot(fname,self%grid,buf)
		deallocate(buf)

	end subroutine
end module
