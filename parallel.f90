module parallel_m
	use kinds_m
	use globals_m
	use teletype_m
	implicit none
#ifdef HAVE_MPI
	include 'mpif.h'
	integer(I4), parameter :: mpi_compiled=1
	integer(I4), parameter :: GATHER_DATA_TAG=1
	integer(I4), parameter :: DATA_SIZE_TAG=2
#else
	integer(I4), parameter :: mpi_compiled=0
#endif

	public start_mpi, stop_mpi, get_mpi_rank, rankname, schedule
	public gather_data, collect_sum
	private

	interface gather_data
		module procedure gather_data1d
		module procedure gather_data2d
		module procedure gather_data1d_tensor
		module procedure gather_data2d_tensor
	end interface

	integer(I4) :: ierr
	integer(I4) :: rank
	character(80) :: sys
contains
	function start_mpi() result(mpirank)
		integer(I4) :: mpirank
		integer(I4) :: hostnm
		
		if (mpi_compiled == 0) then
			rank=-1
			mpirank=-1
			return
		end if
#ifdef HAVE_MPI
		ierr=hostnm(sys)
		call msg_note('Initializing MPI on ' // trim(sys))
		call nl
		call mpi_init(ierr)
		call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
		if (rank == 0) then
			master_p=.true.
		end if
		mpirank=rank
#endif
	end function
		
	subroutine stop_mpi()
#ifdef HAVE_MPI
		if ( mpirun_p ) then
			call mpi_finalize(ierr)
		end if
#endif
	end subroutine

	function get_mpi_rank() result(r)
		integer(I4) :: r
		r=rank
	end function


!
! Generate names based on the execution rank
!
	subroutine rankname(string)
		character(*), intent(inout) :: string

		write(string,'(a1,i3.3)') 'r', rank

	end subroutine

	subroutine schedule(npts, lo, hi)
		integer(I4), intent(in) :: npts
		integer(I4), intent(out) :: lo, hi

		integer(I4) :: sz, num, reminder, i

		if (.not.mpirun_p) then
			lo=1
			hi=npts
			return
		end if
#ifdef HAVE_MPI
		call mpi_comm_size(MPI_COMM_WORLD, sz,ierr)
		reminder=mod(npts,sz)
		num=(npts-reminder)/sz
		hi=0
		do i=0,rank
			lo=hi+1
			hi=hi+num
			if (rank < reminder) hi=hi+1
		end do
#endif
	end subroutine

	subroutine gather_data1d(dest, source)
		real(DP), dimension(:), intent(out) :: dest
		real(DP), dimension(:), intent(in) :: source
#ifdef HAVE_MPI
		integer(I4) :: nproc, client, npts
		integer(I4) :: lo, hi 
		npts=size(source)

		if (.not.mpirun_p) then
			dest=source
			return
		end if
!        call mpi_barrier(MPI_COMM_WORLD, ierr)

		if (rank == 0) then
			call mpi_comm_size(MPI_COMM_WORLD, nproc,ierr)
			dest(1:npts)=source
			lo=npts+1
			do client=1,nproc-1
				call mpi_recv(npts,1,MPI_INTEGER,client,&
					DATA_SIZE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				hi=lo+npts-1
				call mpi_recv(dest(lo:hi),npts,MPI_DOUBLE_PRECISION,client,&
					GATHER_DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				lo=hi+1
			end do
		else
			call mpi_send(npts,1,MPI_INTEGER,0, DATA_SIZE_TAG, &
				MPI_COMM_WORLD, ierr)
			call mpi_send(source,npts,MPI_DOUBLE_PRECISION,0,GATHER_DATA_TAG, &
				MPI_COMM_WORLD, ierr)
		end if
#endif
	end subroutine

	subroutine gather_data2d(dest, source)
		real(DP), dimension(:,:), intent(out) :: dest
		real(DP), dimension(:,:), intent(in) :: source

		integer(I4) :: nproc, client, ni, nj
		integer(I4) :: lo, hi,i 

		if (.not.mpirun_p) then
			dest=source
			return
		end if
#ifdef HAVE_MPI
		ni=size(source(:,1))
		nj=size(source(1,:))
		if (size(dest(:,1)) /= ni) then
			call msg_error('gather_data2d(): dimension mismatch!')
			stop
		end if

!        call mpi_barrier(MPI_COMM_WORLD, ierr)

		if (rank == 0) then
			call mpi_comm_size(MPI_COMM_WORLD, nproc,ierr)
			dest(:,1:nj)=source(:,1:nj)
			lo=nj+1
			do client=1,nproc-1
				call mpi_recv(nj,1,MPI_INTEGER,client,&
					DATA_SIZE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				hi=lo+nj-1
				do i=lo,lo+nj-1
					call mpi_recv(dest(:,i),ni,MPI_DOUBLE_PRECISION,client,&
					GATHER_DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				end do
				lo=hi+1
			end do
		else
			call mpi_send(nj,1,MPI_INTEGER,0, DATA_SIZE_TAG, &
				MPI_COMM_WORLD, ierr)
			do i=1,nj
				call mpi_send(source(:,i),ni,MPI_DOUBLE_PRECISION,0, &
					GATHER_DATA_TAG, MPI_COMM_WORLD, ierr)
			end do
		end if
#endif
	end subroutine

	subroutine gather_data1d_tensor(dest, source)
		type(tensor_t), dimension(:), intent(out) :: dest
		type(tensor_t), dimension(:), intent(in) :: source

		integer(I4) :: nproc, client, ni
		integer(I4) :: lo, hi,i 
		real(8), dimension(:), allocatable :: tmp

		if (.not.mpirun_p) then
			call copy_tensor(source,dest)
			return
		end if
#ifdef HAVE_MPI

		if (master_p) call msg_warn('gather_1d_tensor(): not tested, pls vrfy')

		ni=size(source)

		allocate(tmp(ni*9))

		if (rank == 0) then
			call mpi_comm_size(MPI_COMM_WORLD, nproc,ierr)
			call copy_tensor(source(1:ni),dest(1:ni))
			lo=ni+1
			do client=1,nproc-1
				call mpi_recv(ni,1,MPI_INTEGER,client,&
					DATA_SIZE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				hi=lo+ni-1
				call mpi_recv(tmp,ni*9,MPI_DOUBLE_PRECISION,client,&
				GATHER_DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				call pack_tensors(tmp,dest(lo:hi)) 
				lo=hi+1
			end do
		else
			call mpi_send(ni,1,MPI_INTEGER,0, DATA_SIZE_TAG, &
				MPI_COMM_WORLD, ierr)
			call unpack_tensors(source, tmp)
			call mpi_send(tmp,ni*9,MPI_DOUBLE_PRECISION,0, &
				GATHER_DATA_TAG, MPI_COMM_WORLD, ierr)
		end if
		deallocate(tmp)
#endif
	end subroutine

	subroutine gather_data2d_tensor(dest, source)
		type(tensor_t), dimension(:,:), intent(out) :: dest
		type(tensor_t), dimension(:,:), intent(in) :: source

		integer(I4) :: nproc, client, ni, nj
		integer(I4) :: lo, hi,i 
		real(8), dimension(:), allocatable :: tmp

		if (.not.mpirun_p) then
			call copy_tensor(source,dest)
			return
		end if
#ifdef HAVE_MPI

		ni=size(source(:,1))
		nj=size(source(1,:))
		if (size(dest(:,1)) /= ni) then
			call msg_error('gather_data2d(): dimension mismatch!')
			stop
		end if

		allocate(tmp(ni*9))

		if (rank == 0) then
			call mpi_comm_size(MPI_COMM_WORLD, nproc,ierr)
			call copy_tensor(source(:,1:nj),dest(:,1:nj))
			lo=nj+1
			do client=1,nproc-1
				call mpi_recv(nj,1,MPI_INTEGER,client,&
					DATA_SIZE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				hi=lo+nj-1
				do i=lo,lo+nj-1
					call mpi_recv(tmp,ni*9,MPI_DOUBLE_PRECISION,client,&
					GATHER_DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
					call pack_tensors(tmp,dest(:,i)) 
				end do
				lo=hi+1
			end do
		else
			call mpi_send(nj,1,MPI_INTEGER,0, DATA_SIZE_TAG, &
				MPI_COMM_WORLD, ierr)
			do i=1,nj
				call unpack_tensors(source(:,i), tmp)
				call mpi_send(tmp,ni*9,MPI_DOUBLE_PRECISION,0, &
					GATHER_DATA_TAG, MPI_COMM_WORLD, ierr)
			end do
		end if
		deallocate(tmp)
#endif
	end subroutine

#ifdef HAVE_MPI
	subroutine pack_tensors(tmp, foo)
		real(DP), dimension(:) :: tmp
		type(tensor_t), dimension(:) :: foo
		
		integer(I4) :: i,j,k,l

		k=1
		do l=1,size(foo)
			do j=1,3
				do i=1,3
					foo(l)%t(i,j)=tmp(k)
					k=k+1
				end do
			end do
		end do
	end subroutine

	subroutine unpack_tensors(foo,tmp)
		type(tensor_t), dimension(:) :: foo
		real(DP), dimension(:) :: tmp
		
		integer(I4) :: i,j,k,l

		k=1
		do l=1,size(foo)
			do j=1,3
				do i=1,3
					tmp(k)=foo(l)%t(i,j)
					k=k+1
				end do
			end do
		end do

	end subroutine
#endif
	
	subroutine collect_sum(source, dest)
		real(DP), intent(in) :: source
		real(DP), intent(out) :: dest

		if (.not.mpirun_p) then
			dest=source
			return
		end if
#ifdef HAVE_MPI
		call mpi_reduce(source,dest,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
			MPI_COMM_WORLD,ierr)
#endif
	end subroutine

end module
