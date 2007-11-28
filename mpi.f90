module mpi_m
	use kinds_m
	implicit none
#ifdef HAVE_MPI
	include 'mpif.h'
	integer(I4), parameter :: mpi_compiled=1
#else
	integer(I4), parameter :: mpi_compiled=0
#endif
	integer(I4), parameter :: DIE_TAG=-1
	integer(I4), parameter :: JOB_TAG=1
	integer(I4), parameter :: JNUM_TAG=2

	public start_mpi, stop_mpi, get_mpi_rank, rankname
	private

	integer(I4) :: ierr
	integer(I4) :: rank
	character(80) :: sys
contains
	function start_mpi() result(mpirank)
		integer(I4) :: mpirank
		external hostnm
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

end module
