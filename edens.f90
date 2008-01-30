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
	implicit none
	
	public init_edens, del_edens, edens_t
	public edens_direct, edens, edens_plot
	public edens_classaster, edens_slave
	
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
	subroutine init_edens(ed, mol, dens, grid)
		type(edens_t) :: ed
		type(molecule_t) :: mol
		type(dens_t), target :: dens
		type(grid_t), target :: grid

		integer(I4) :: n
		integer(I4) :: p1, p2

		call get_grid_size(grid, p1, p2)
		call get_dens(dens, ed%aodens)
		n=get_ncgto(mol)
		allocate(ed%tmp(n))
		allocate(ed%buf(p1,p2))
		call init_bfeval(ed%bf, mol)
		ed%grid=>grid

		if (master_p) then
			open(EDFD, file='EDENS', access='direct', recl=p1*p2*DP)
		end if
	end subroutine

	subroutine del_edens(ed)
		type(edens_t) :: ed

		if (associated(ed%tmp)) deallocate(ed%tmp)
		if (associated(ed%buf)) deallocate(ed%buf)
		call del_bfeval(ed%bf)
		nullify(ed%aodens,ed%grid)
		if (master_p) then
			close(EDFD)
		end if
	end subroutine

	subroutine set_edens(ed, k)
		type(edens_t), intent(in) :: ed
		integer(I4), intent(in) :: k

		write(EDFD, rec=k) ed%buf
	end subroutine

	subroutine edens_plot(ed)
		type(edens_t), intent(inout) :: ed
		
		integer(I4), dimension(:), pointer :: z
		integer(I4) :: i,j,k,p1,p2,p3
		real(DP) :: amax
		real(DP), dimension(3) :: rr
		real(DP), dimension(:,:), pointer :: buf

		call get_grid_size(ed%grid, p1, p2,p3)
		buf=>ed%buf
		if (keyword_is_set(input,'edens.plots')) then
!            call get_kw_size('edens.plots', i)
!            allocate(z(i))
			call getkw(input, 'edens.plots', z)
		else
			return
!            allocate(z(1))
!            z=1
		end if

		do k=1,size(z)
			amax=D0
			if (z(k) < 1 .or. z(k) > p3 ) then
				write(str_g, '(a,i6)') 'edens_plot(): &
					&invalid value for polts:', z(k)
				call msg_error(str_g)
				cycle
			end if
			read(EDFD, rec=z(k)) ed%buf
			str_g=enumfile('EDENSPLT', k)
			open(EDPFD, file=trim(str_g))
			do j=1,p2
				do i=1,p1
					rr=gridpoint(ed%grid, i, j, z(k))
					write(EDPFD, '(4f19.8)') rr, buf(i,j)
					if (abs(buf(i,j)) > amax) amax=abs(buf(i,j))
				end do
				write(EDPFD, *) 
			end do
			close(EDPFD)
			write(str_g, '(a,e19.12)') 'Max electronic density:', amax
			call msg_info(str_g)
		end do
		call edens_gopenmol(ed)
	end subroutine

	subroutine edens_direct(ed, k)
		type(edens_t) :: ed
		integer(I4), intent(in) :: k

		integer(I4) :: i, j, p1, p2
		real(DP), dimension(3) :: rr
		real(DP), dimension(:,:), pointer :: buf
		real(DP), dimension(:), pointer :: bfvec

		call get_grid_size(ed%grid, p1, p2)

		buf=>ed%buf

		do j=1,p2
			do i=1,p1
				rr=gridpoint(ed%grid, i, j, k)
				call bfeval(ed%bf,rr, bfvec)
				ed%tmp=matmul(ed%aodens, bfvec)
				ed%buf(i,j)=dot_product(ed%tmp, bfvec)
			end do
		end do
	end subroutine

	subroutine edens(ed)
		type(edens_t) :: ed

		integer(I4) :: i, j, k, p1, p2, p3
		real(DP), dimension(3) :: rr
		real(DP), dimension(:,:), pointer :: buf
		real(DP), dimension(:), pointer :: bfvec

		call get_grid_size(ed%grid, p1, p2, p3)

		buf=>ed%buf
	
		do k=1,p3
			do j=1,p2
				do i=1,p1
					rr=gridpoint(ed%grid, i, j, k)
					call bfeval(ed%bf,rr, bfvec)
					ed%tmp=matmul(ed%aodens, bfvec)
					ed%buf(i,j)=dot_product(ed%tmp, bfvec)
				end do
			end do
			write(EDFD, rec=k) ed%buf
		end do
	end subroutine

	subroutine edens_classaster(ed, who)
		use mpi_m
		type(edens_t) :: ed
		integer(I4), intent(in) :: who

#ifdef HAVE_MPI
		integer(I4) :: p1, p2, jnum
		integer(I4) :: bufsz, ierr
		integer(I4), dimension(MPI_STATUS_SIZE) :: stat
		real(DP), dimension(:,:), allocatable :: buf
		real(DP) :: amax

		call get_grid_size(ed%grid, p1, p2)

		bufsz=size(ed%buf)

		call mpi_recv(jnum,1,MPI_INTEGER,who,EDENS_TAG, &
			MPI_COMM_WORLD, stat, ierr)
		call mpi_recv(ed%buf,bufsz,MPI_DOUBLE_PRECISION,who, &
			EDENS_TAG, MPI_COMM_WORLD, stat, ierr)

		call set_edens(ed, jnum)

#endif
	end subroutine

	subroutine edens_slave(ed, jnum)
		use mpi_m
		type(edens_t) :: ed
		type(grid_t) :: grid
		integer(I4), intent(in) :: jnum

#ifdef HAVE_MPI
		integer(I4) :: bufsz, ierr
		integer(I4), dimension(MPI_STATUS_SIZE) :: stat

		bufsz=size(ed%buf)
		ed%buf=D0

		call edens_direct(ed, jnum)
		call mpi_send(EDENS_TAG, 1, MPI_INTEGER, 0, JOB_TAG,&
			MPI_COMM_WORLD, ierr)
		call mpi_send(jnum, 1, MPI_INTEGER, 0, EDENS_TAG,&
			MPI_COMM_WORLD, ierr)
		call mpi_send(ed%buf, bufsz, MPI_DOUBLE_PRECISION, 0, EDENS_TAG, &
			MPI_COMM_WORLD, ierr)
#endif
	end subroutine
		
	subroutine edens_gopenmol(ed)
		type(edens_t) :: ed

		integer(I4) :: surface, rank, p1, p2, p3
		integer(I4) :: i, j, k, l
		real(SP), dimension(3) :: qmin, qmax
		real(DP), dimension(:,:), pointer :: buf
		character(BUFLEN) :: gopen_file

		buf=>ed%buf
		gopen_file=''
		call getkw(input, 'edens.gopenmol', gopen_file)
		if (trim(gopen_file) == '') return
		open(GOPFD,file=trim(gopen_file),access='direct',recl=I4)

		surface=200
		rank=3

		call get_grid_size(ed%grid, p1, p2, p3)
		qmin=gridpoint(ed%grid,1,1,1)*AU2A
		qmax=gridpoint(ed%grid,p1,p2,p3)*AU2A

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
end module
