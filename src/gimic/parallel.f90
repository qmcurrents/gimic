module parallel_m
    use kinds_m
    use globals_m
    use settings_m
    use teletype_m
    implicit none
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    public init_mpi, stop_mpi, get_mpi_rank, rankname, schedule
    public gather_data, collect_sum
    private

    integer(I4) :: ierr
    integer(I4) :: rank
    character(80) :: sys
contains
    function init_mpi(is_mpirun) result(mpirank)
        logical, intent(out) :: is_mpirun
        integer(I4) :: mpirank
        integer(I4) :: hostnm

        mpi_world_size = 1
#ifdef HAVE_MPI
        ierr=hostnm(sys)
        call msg_note('Initializing MPI on ' // trim(sys))
        call nl
        call mpi_init(ierr)
        call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
        call mpi_comm_size(MPI_COMM_WORLD, mpi_world_size, ierr)
        mpirank=rank
#else
        rank = 0
        mpirank = 0
#endif
        if (mpi_world_size == 1) then
            settings%is_mpirun = .false.
        else
            settings%is_mpirun = .true.
        end if
    end function
        
    subroutine stop_mpi()
#ifdef HAVE_MPI
        call mpi_finalize(ierr)
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

        integer(I4) :: nPerHost, nLeftOver 

        if (mpi_world_size == 1) then
            lo=1
            hi=npts 
            return
        end if
#ifdef HAVE_MPI

        nPerHost = npts / mpi_world_size
        nLeftOver = mod(npts, mpi_world_size)
        lo = nLeftOver + rank * nPerHost + 1
        hi = lo + nPerHost - 1
#endif
    end subroutine

    subroutine gather_data(source, dest)
        real(DP), dimension(:,:), intent(in) :: source
        real(DP), dimension(:,:), intent(out) :: dest

        if (mpi_world_size == 1) then
            dest=source
            return
        end if
#ifdef HAVE_MPI
!        call mpi_barrier(MPI_COMM_WORLD, ierr)

        call mpi_gather(source, size(source), MPI_DOUBLE_PRECISION, &
            dest, size(source), MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr)
!        call mpi_gather(dest, n, MPI_DOUBLE_PRECISION, &
!            MPI_IN_PLACE, n, MPI_DOUBLE_PRECISION, 0, &
!            MPI_COMM_WORLD, ierr)
        if (ierr /= 0) then
            stop 'MPI error in gather_data()'
        end if
#endif
    end subroutine

    subroutine collect_sum(source, dest)
        real(DP), intent(in) :: source
        real(DP), intent(out) :: dest

        if (mpi_world_size == 1) then
            dest=source
            return
        end if
#ifdef HAVE_MPI
        call mpi_reduce(source, dest, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
            MPI_COMM_WORLD, ierr)
#endif
    end subroutine

end module

! vim:et:sw=4:ts=4
