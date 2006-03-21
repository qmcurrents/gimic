!
! $Id$
!

module mpi_m
	use precision_m
	implicit none
#ifdef HAVE_MPI
	include 'mpif.h'
	integer(I4), parameter :: mpi_compiled=1

	integer(I4), parameter :: DIE_TAG=-1
	integer(I4), parameter :: JOB_TAG=1
	integer(I4), parameter :: JNUM_TAG=2
#else
	integer(I4), parameter :: mpi_compiled=0
#endif
end module
