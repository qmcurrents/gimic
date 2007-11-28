!
! This module contains routines to handle communication between MPI nodes.
!

module parallel_class
	use globals_m
	use mpi_m
	use grid_class 
	use jfield_class 
	use divj_class
	use dens_class
	use edens_class
	use integral_class
	use teletype_m
	implicit none

	! deprecated
	type parallel_t
		logical :: it_p, dj_p, jf_p, ed_p
		type(divj_t), pointer :: dj
		type(edens_t), pointer :: ed
		type(jfield_t), pointer :: jf
		type(integral_t), pointer :: it
		integer(I4) :: njobs
		integer(I4), dimension(:), pointer :: work, jnum
	end type

	public start_mpi, stop_mpi, get_mpi_rank
!    public init_parallel, del_parallel, parallel_t
!    public bcast_dens, bcast_molbuf 
!    public bcast_inpbuf, scheduler, ceo
	private

	integer(I4) :: ierr
	integer(I4) :: rank
	character(80) :: sys

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Deprecated routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine init_parallel(pt, calc, jf, it, dj, ed)
		type(parallel_t) :: pt
		integer(I4), dimension(:), intent(in) :: calc
		type(jfield_t), target :: jf
		type(integral_t), target :: it
		type(divj_t), target :: dj
		type(edens_t), target :: ed

#ifdef HAVE_MPI
!        pt%it_p=.false.; pt%jf_p=.false.; 
!        pt%dj_p=.false.; pt%ed_p=.false.;
!        call getkw(input, 'cdens', pt%jf_p)
!        call getkw(input, 'integral', pt%it_p)
!        call getkw(input, 'divj', pt%dj_p)
!        call getkw(input, 'edens', pt%ed_p)

		pt%jf=>jf
		pt%dj=>dj
		pt%it=>it
		pt%ed=>ed

		call setup_work_vectors(pt, calc)
#endif
	end subroutine

	subroutine del_parallel(pt)
		type(parallel_t) :: pt

		nullify(pt%dj)
		nullify(pt%jf)
		nullify(pt%ed)
		deallocate(pt%work)
		deallocate(pt%jnum)
	end subroutine

	subroutine scheduler(pt)
		type(parallel_t) :: pt

#ifdef HAVE_MPI
		integer(I4) :: who, ncpu, nodes, job, gotitm, gotit, q, tun
		integer(I4), dimension(MPI_STATUS_SIZE) :: stat

		call mpi_comm_size(MPI_COMM_WORLD, nodes, ierr)

		! inital distribution of work
		if (pt%njobs <= nodes) then 
			ncpu=pt%njobs
		else
			ncpu=nodes-1
		end if

		do q=1,ncpu
           	call mpi_send(pt%work(q),1,MPI_INTEGER,q,JOB_TAG,&
				MPI_COMM_WORLD, ierr)
           	call mpi_send(pt%jnum(q),1,MPI_INTEGER,q,JNUM_TAG,&
				MPI_COMM_WORLD, ierr)
		end do

		gotit=0
		do while (gotit < pt%njobs) !??
			call mpi_recv(job,1,MPI_INTEGER,MPI_ANY_SOURCE,JOB_TAG, &
				MPI_COMM_WORLD, stat, ierr)
			who=stat(MPI_SOURCE)
			gotit=gotit+1
			select case(job)
				case(INTGRL_TAG)
					call integral_classaster(pt%it, who)
				case(CDENS_TAG)
					call jfield_classaster(pt%jf, who)
				case(DIVJ_TAG)
					call divj_classaster(pt%dj, who)
				case(EDENS_TAG)
					call edens_classaster(pt%ed, who)
				case default
					call msg_error('sheduler(): Unknown command')
					call exit(1)
			end select
			
			if ( q <= pt%njobs ) then
				call mpi_send(pt%work(q), 1, MPI_INTEGER, who, JOB_TAG,&
					MPI_COMM_WORLD, ierr)
				call mpi_send(pt%jnum(q), 1, MPI_INTEGER, who, JNUM_TAG,&
					MPI_COMM_WORLD, ierr)
				q=q+1
			end if
		end do

		! tell the puppets to quit
		do q=1,nodes-1
           	call mpi_send(DIE_TAG,1,MPI_INTEGER,q,JOB_TAG,&
				MPI_COMM_WORLD, ierr)
		end do

#endif
	end subroutine

	subroutine setup_work_vectors(pt, calc)
		type(parallel_t) :: pt
		integer(I4), dimension(:), intent(in) :: calc

		integer(I4) :: i,j,k,l
		integer(I4) :: p1, p2, jf, dj, it, ed
		integer(I4), dimension(:), pointer :: work, jnum

		jf=0; dj=0; it=0; ed=0;
		! ok, we break the oo rules here, but what the heck...
		do i=1,size(calc)
			select case(calc(i))
			case (CDENS_TAG) 
				call get_grid_size(pt%jf%grid, p1, p2, jf)
			case (INTGRL_TAG) 
				call get_grid_size(pt%it%grid, p1, p2, it)
			case (DIVJ_TAG) 
				call get_grid_size(pt%dj%grid, p1, p2, dj)
			case (EDENS_TAG) 
				call get_grid_size(pt%ed%grid, p1, p2, ed)
			end select
		end do

		pt%njobs=jf+it+dj+ed
		allocate(pt%work(pt%njobs))
		allocate(pt%jnum(pt%njobs))
		work=>pt%work
		jnum=>pt%jnum

		j=1; k=1
		do l=1,size(calc)
			select case(calc(l))
			case (CDENS_TAG) 
				do i=j,j+jf-1
					work(i)=CDENS_TAG
					jnum(i)=k
					k=k+1
				end do
				j=i; k=1
			case (INTGRL_TAG) 
				do i=j,j+it-1
					work(i)=INTGRL_TAG
					jnum(i)=k
					k=k+1
				end do
				j=i; k=1
			case (DIVJ_TAG) 
				do i=j,j+dj-1
					work(i)=DIVJ_TAG
					jnum(i)=k
					k=k+1
				end do
				j=i; k=1
			case (EDENS_TAG) 
				do i=j,j+ed-1
					work(i)=EDENS_TAG
					jnum(i)=k
					k=k+1
				end do
				j=i; k=1
			end select
		end do
	end subroutine

	subroutine ceo(pt)
		type(parallel_t) :: pt

#ifdef HAVE_MPI
		integer(I4) :: job, jnum
		integer(I4), dimension(MPI_STATUS_SIZE) :: stat

		call msg_debug('The CEO is in charge on ' // trim(sys), 2)
		do while (1)
			call mpi_recv(job, 1, MPI_INTEGER, 0, JOB_TAG, &
				MPI_COMM_WORLD, stat, ierr)
			if ( job < 0 ) then
				call msg_debug('got quit... ' // trim(sys), 2)
				exit
			end if
			call mpi_recv(jnum, 1, MPI_INTEGER, 0, JNUM_TAG, &
				MPI_COMM_WORLD, stat, ierr)
			select case(job)
				case(INTGRL_TAG)
					call msg_debug('integral_slave on ' // trim(sys) // &
						': ' // xchar(jnum), 2)
					call integral_slave(pt%it, jnum)
				case(CDENS_TAG)
					call msg_debug('cotton pick on ' // trim(sys) // &
						': ' // xchar(jnum), 2)
					call jfield_slave(pt%jf, jnum)
				case(DIVJ_TAG)
					call msg_debug('divj_slave on ' // trim(sys) // &
						': ' // xchar(jnum), 2)
					call divj_slave(pt%dj, jnum)
				case(EDENS_TAG)
					call msg_debug('edens_slave on ' // trim(sys) // &
						': ' // xchar(jnum), 2)
					call edens_slave(pt%ed, jnum)
				case default
					call msg_error('CEO(): Unknown command')
			end select
		end do
		call msg_debug('The CEO resigns on ' // trim(sys), 2)
#endif
	end subroutine

	subroutine send_dens(xdens)
		type(dens_t), intent(in) :: xdens

#ifdef HAVE_MPI
		real(DP), dimension(:,:), pointer :: dd
		integer(I4) :: i, n
		integer(I4), dimension(2) :: nn

		call mpi_barrier(MPI_COMM_WORLD, ierr)
		call get_dens(xdens, dd)
		n=size(dd)
		nn=shape(dd)
		call mpi_bcast(dd,n,MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, ierr)
		call mpi_barrier(MPI_COMM_WORLD, ierr)
!        print *, 'dens sent', dd(1,1), dd(nn(1),nn(2))
		do i=1,3
			call get_pdens(xdens, i, dd)
			n=size(dd)
			call mpi_bcast(dd,n,MPI_DOUBLE_PRECISION, rank, &
				MPI_COMM_WORLD, ierr)
			call mpi_barrier(MPI_COMM_WORLD, ierr)
!            print *, 'pdens sent', i, dd(1,1), dd(nn(1),nn(2))
		end do
#endif
	end subroutine

	subroutine bcast_dens(xdens)
		type(dens_t), intent(inout) :: xdens
#ifdef HAVE_MPI
		real(DP), dimension(:,:), pointer :: dd
		integer(I4) :: i, n
		integer(I4), dimension(2) :: nn

		call mpi_barrier(MPI_COMM_WORLD, ierr)
		call get_dens(xdens, dd)
		n=size(dd)
		nn=shape(dd)
		call mpi_bcast(dd,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
		call mpi_barrier(MPI_COMM_WORLD, ierr)
		call set_dens(xdens, dd)
        write(str_g, *) 'dens recvd: ', trim(sys), rank
		call msg_debug(str_g, 1)
		if (xdens%pdens_p) then
			do i=1,3
				call get_pdens(xdens, i, dd)
				call mpi_bcast(dd,n,MPI_DOUBLE_PRECISION, 0, &
					MPI_COMM_WORLD, ierr)
				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call set_pdens(xdens, i, dd)
				write(str_g, *) 'pdens recvd: ', trim(sys), rank
				call msg_debug(str_g, 1)
			end do
		end if
#endif
	end subroutine

	subroutine bcast_molbuf(molbuf)
		character(MAX_LINE_LEN), dimension(:), pointer :: molbuf
#ifdef HAVE_MPI
		integer(I4) ::  n
		
		call mpi_barrier(MPI_COMM_WORLD, ierr)
		if (associated(molbuf)) then
			n=size(molbuf)
		end if
		call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
		if (.not.associated(molbuf)) then
			allocate(molbuf(n))
		end if
		n=size(molbuf)*len(molbuf(1))
		call mpi_bcast(molbuf,n,MPI_CHARACTER,0,MPI_COMM_WORLD, ierr)
		call mpi_barrier(MPI_COMM_WORLD, ierr)
        write(str_g, *) 'molbuf recvd: ', trim(sys), rank
		call msg_debug(str_g, 5)
#endif
	end subroutine

	subroutine bcast_inpbuf(inpbuf)
		character, dimension(:), pointer :: inpbuf
#ifdef HAVE_MPI
		integer(I4) :: n

		call mpi_barrier(MPI_COMM_WORLD, ierr)
		if (associated(inpbuf)) then
			n=size(inpbuf)
		end if
		call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
		if (.not.associated(inpbuf)) then
			allocate(inpbuf(n))
		end if
		n=size(inpbuf)
		call mpi_bcast(inpbuf,n,MPI_CHARACTER,0,MPI_COMM_WORLD, ierr)
		call mpi_barrier(MPI_COMM_WORLD, ierr)
        write(str_g, *) 'inbuf recvd: ', trim(sys), rank
		call msg_debug(str_g, 5)
#endif
	end subroutine

end module
