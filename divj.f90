!
!  Calculate the divergence of a vector field.
!  Useful for checking the convergce of the actual gauge invariance,
!  in contrast to the basis set convergence and gauge independence.
!

module divj_class
	use globals_m
	use grid_class
	use jfield_class
	use jtensor_class
	use teletype_m
	implicit none

	type divj_t
		real(DP), dimension(:,:), pointer :: buf
		type(grid_t), pointer :: grid
		type(jtensor_t), pointer :: jt
		real(DP), dimension(3) :: bb
	end type

	public divj_direct, divj_direct_plt, init_divj, del_divj
	public divj_plot, set_divj, divj_t, divj
	public divj_classaster, divj_slave
	private
	
	real(DP), dimension(5), parameter :: wgt=(/2.d0,-16.d0,0.d0,16.d0,-2.d0/)
	real(DP), parameter :: step=1.d-3
	real(DP), parameter :: hx=1.d0/(24.d0*step)
	real(DP), parameter :: hy=1.d0/(24.d0*step)
	real(DP), parameter :: hz=1.d0/(24.d0*step)

contains
	subroutine init_divj(dj, grid, jt)
		type(divj_t), intent(inout) :: dj
		type(grid_t), target :: grid
		type(jtensor_t), target :: jt

		integer(I4) :: djrl, p1, p2
		logical :: foo_p

		dj%bb=D0
		call getkw(input, 'cdens.magnet', dj%bb)
		call getkw(input, 'cdens.orthogonal_magnet', foo_p)
		if (foo_p) then
			call msg_note('init_divj(): &
				&Magnetic field defined to be orthogonal to the grid')
			dj%bb=get_grid_normal(dj%grid)
		end if

		call get_grid_size(grid, p1, p2)
		djrl=p1*p2*DP
		if (master_p) then
			open(DIVJFD, file='DIVJ', access='direct', recl=djrl)
		end if

		dj%grid=>grid
		dj%jt=>jt
		allocate(dj%buf(p1,p2))
	end subroutine

	subroutine del_divj(dj)
		type(divj_t), intent(inout) :: dj

		deallocate(dj%buf)
		if (master_p) then
			close(DIVJFD)
		end if
	end subroutine

	subroutine set_divj(dj, k)
		type(divj_t), intent(in) :: dj
		integer(I4), intent(in) :: k

		write(DIVJFD, rec=k) dj%buf
	end subroutine

	subroutine divj_plot(dj)
		type(divj_t), intent(inout) :: dj
		
		integer(I4), dimension(:), pointer :: z
		integer(I4) :: i,j,k,p1,p2,p3
		real(DP) :: amax
		real(DP), dimension(3) :: rr
		real(DP), dimension(:,:), pointer :: buf

		call get_grid_size(dj%grid, p1, p2, p3)

		buf=>dj%buf
		if (keyword_is_set(input, 'divj.plots')) then
!            call get_kw_size('divj.plots', i)
!            allocate(z(i))
			call getkw_ptr(input, 'divj.plots', z)
		else
!            allocate(z(1))
!            z=1
		return
		end if

		do k=1,size(z)
			amax=D0
			if (z(k) < 1 .or. z(k) > p3 ) then
				write(str_g, '(a,i6)') 'divj_plot(): &
					&invalid value for polts:', z(k)
				call msg_error(str_g)
				cycle
			end if
			read(DIVJFD, rec=z(k)) dj%buf
			open(DJPFD, file='DIVJPLT' // '.' // trim(xchar(k)))
			do j=1,p2
				do i=1,p1
					rr=gridpoint(dj%grid, i, j, z(k))
					write(DJPFD, '(4f)') rr, buf(i,j)
					if (abs(buf(i,j)) > amax) amax=abs(buf(i,j))
				end do
				write(DJPFD, *) 
			end do
			close(DJPFD)
			write(str_g, '(a,f)') 'Max divergence:', amax
			call msg_info(str_g)
		end do
		call divj_gopenmol(dj)
	end subroutine

	subroutine divj_direct_plt(dj)
		type(divj_t), intent(inout) :: dj

		integer(I4) :: i, j, k, p1, p2, p3
		real(DP) :: div, amax
		real(DP), dimension(3) :: rr 
		
		call eta(dj%jt, dj%grid, 12.d0)
		call get_grid_size(dj%grid, p1, p2, p3)
	
		open(DJPFD, file='DIVJPLT')
		amax=D0
		do k=1,p3
			do j=1,p2
				do i=1,p1
					rr=gridpoint(dj%grid, i, j, k)
					call divergence(dj, rr, div)
!                    call divergence2(dj, rr, div)
					write(DJPFD, '(4f)') rr, div
					if (abs(div) > amax) amax=abs(div)
				end do
				write(DJPFD, *) 
			end do
		end do
		write(str_g, '(a,f)')'Max divergence:', amax
		call msg_note(str_g)
		call nl
		close(DJPFD)
	end subroutine

	subroutine divj_direct(dj, k)
		type(divj_t) :: dj
		integer(I4), intent(in) :: k

		integer(I4) :: i, j, p1, p2
		real(DP), dimension(3) :: rr
		real(DP), dimension(:,:), pointer :: buf

		call get_grid_size(dj%grid, p1, p2)

		buf=>dj%buf

		do j=1,p2
			do i=1,p1
				rr=gridpoint(dj%grid, i, j, k)
				call divergence(dj,rr,buf(i,j))
!                call divergence2(dj,rr,buf(i,j))
			end do
		end do
	end subroutine

	subroutine divj(dj)
		type(divj_t) :: dj

		integer(I4) :: i, j, k, p, p1, p2, p3
		real(DP), dimension(3) :: rr
		real(DP), dimension(:,:), pointer :: buf

		call get_grid_size(dj%grid, p1, p2, p3)

		buf=>dj%buf
	
		do k=1,p3
			do j=1,p2
				do i=1,p1
					rr=gridpoint(dj%grid, i, j, k)
					call divergence(dj,rr,buf(i,j))
!                    call divergence2(dj,rr,buf(i,j))
				end do
			end do
			write(DIVJFD, rec=k) dj%buf
		end do
	end subroutine

	subroutine divergence2(dj, rr, div)
		type(divj_t), intent(in) :: dj
		real(DP), dimension(3), intent(in) :: rr
		real(DP), intent(out) :: div

		integer(I4) :: q
		type(tensor_t) :: jtxp, jtyp, jtzp
		type(tensor_t) :: jtxd, jtyd, jtzd
		real(DP) :: djx,djy,djz
		real(DP), dimension(5) :: jx,jy,jz
		real(DP), dimension(3) :: tvec

		jx=D0; jy=D0; jz=D0
		do q=-2,2
			call ctensor2(dj%jt, rr+(/step*real(q),D0,D0/), jtxp,jtxd, 'total')
			call ctensor2(dj%jt, rr+(/D0,step*real(q),D0/), jtyp,jtyd, 'total')
			call ctensor2(dj%jt, rr+(/D0,D0,step*real(q)/), jtzp,jtzd, 'total')
			call jvector(jtxp%t, jtxd%t, dj%bb, tvec)
			jx(q+3)=tvec(1)
			call jvector(jtyp%t, jtyd%t, dj%bb, tvec)
			jy(q+3)=tvec(2)
			call jvector(jtzp%t, jtzd%t, dj%bb, tvec)
			jz(q+3)=tvec(3)
		end do
		djx=hx*sum(wgt*jx)
		djy=hy*sum(wgt*jy)
		djz=hz*sum(wgt*jz)
		div=djx+djy+djz
	end subroutine

	subroutine divergence(dj, rr, div)
		type(divj_t), intent(in) :: dj
		real(DP), dimension(3), intent(in) :: rr
		real(DP), intent(out) :: div

		integer(I4) :: q
		type(tensor_t) :: jtx, jty, jtz
		real(DP) :: djx,djy,djz
		real(DP), dimension(5) :: jx,jy,jz
		real(DP), dimension(3) :: tvec

		jx=D0; jy=D0; jz=D0
		do q=-2,2
			if ( q == 0 ) cycle
			call ctensor(dj%jt, rr+(/step*real(q),D0,D0/), jtx, 'total')
			call ctensor(dj%jt, rr+(/D0,step*real(q),D0/), jty, 'total')
			call ctensor(dj%jt, rr+(/D0,D0,step*real(q)/), jtz, 'total')
			tvec=matmul(jtx%t, dj%bb)
			jx(q+3)=tvec(1)
			tvec=matmul(jty%t, dj%bb)
			jy(q+3)=tvec(2)
			tvec=matmul(jtz%t, dj%bb)
			jz(q+3)=tvec(3)
		end do
		djx=hx*sum(wgt*jx)
		djy=hy*sum(wgt*jy)
		djz=hz*sum(wgt*jz)
		div=djx+djy+djz
	end subroutine

	subroutine divj_classaster(dj, who)
		use mpi_m
		type(divj_t) :: dj
		type(grid_t) :: grid
		integer(I4), intent(in) :: who

#ifdef HAVE_MPI
		integer(I4) :: p1, p2, jnum
		integer(I4) :: bufsz, ierr
		integer(I4), dimension(MPI_STATUS_SIZE) :: stat
		real(DP), dimension(:,:), allocatable :: buf
		real(DP) :: amax

		call get_grid_size(grid, p1, p2)

		bufsz=size(dj%buf)

		call mpi_recv(jnum,1,MPI_INTEGER,who,DIVJ_TAG, &
			MPI_COMM_WORLD, stat, ierr)
		call mpi_recv(dj%buf,bufsz,MPI_DOUBLE_PRECISION,who, &
			DIVJ_TAG, MPI_COMM_WORLD, stat, ierr)

		call set_divj(dj, jnum)

#endif
	end subroutine

	subroutine divj_slave(dj, jnum)
		use mpi_m
		type(divj_t) :: dj
		type(grid_t) :: grid
		integer(I4), intent(in) :: jnum

#ifdef HAVE_MPI
		integer(I4) :: bufsz, ierr
		integer(I4), dimension(MPI_STATUS_SIZE) :: stat

		bufsz=size(dj%buf)
		dj%buf=D0

		call divj_direct(dj, jnum)
		call mpi_send(DIVJ_TAG, 1, MPI_INTEGER, 0, JOB_TAG,&
			MPI_COMM_WORLD, ierr)
		call mpi_send(jnum, 1, MPI_INTEGER, 0, DIVJ_TAG,&
			MPI_COMM_WORLD, ierr)
		call mpi_send(dj%buf, bufsz, MPI_DOUBLE_PRECISION, 0, DIVJ_TAG, &
			MPI_COMM_WORLD, ierr)
#endif
	end subroutine

	subroutine divj_gopenmol(dj)
		type(divj_t) :: dj

		integer(I4) :: surface, rank, p1, p2, p3
		integer(I4) :: i, j, k, l
		real(SP), dimension(3) :: qmin, qmax
		real(DP), dimension(:,:), pointer :: buf
		character(BUFLEN) :: gopen_file

		buf=>dj%buf
		gopen_file=''
		call getkw(input, 'divj.gopenmol', gopen_file)
		if (trim(gopen_file) == '') return
		open(GOPFD,file=trim(gopen_file),access='direct',recl=4)

		surface=200
		rank=3

		call get_grid_size(dj%grid, p1, p2, p3)
		qmin=gridpoint(dj%grid,1,1,1)*AU2A
		qmax=gridpoint(dj%grid,p1,p2,p3)*AU2A

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
			read(DIVJFD, rec=k) buf
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
