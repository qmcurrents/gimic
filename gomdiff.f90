program gomdiff
	implicit none

	type gom_t
		integer(4) :: surface, rank, p1, p2, p3
		real(4), dimension(3) :: qmin, qmax
		real(4), dimension(:,:,:), pointer :: buf
	end type
	
	character(80) :: in1, in2, outf
	integer(4) :: argc, iargc
	type(gom_t) :: gt1, gt2

	argc=iargc()
	if (argc < 3) then
		print *, 'usage: gomdiff plt1 plt2 outfile'
		stop 
	end if
	call getarg(1, in1)
	call getarg(2, in2)
	call getarg(3, outf)

	call read_gomplt(in1, gt1)
	call read_gomplt(in2, gt2)
	call write_gomplt(outf, gt2)

	deallocate(gt1%buf,gt2%buf)

contains
	subroutine read_gomplt(fname,dj)
		character(*) :: fname
		type(gom_t) :: dj

		integer(4) :: i, j, k, l

		open(42,file=trim(fname),access='direct',recl=4)

		read(42,rec=1) dj%rank
		read(42,rec=2) dj%surface
		read(42,rec=3) dj%p3
		read(42,rec=4) dj%p2
		read(42,rec=5) dj%p1
		read(42,rec=6) dj%qmin(3)
		read(42,rec=7) dj%qmax(3)
		read(42,rec=8) dj%qmin(2)
		read(42,rec=9) dj%qmax(2)
		read(42,rec=10) dj%qmin(1)
		read(42,rec=11) dj%qmax(1)

		allocate(dj%buf(dj%p1,dj%p2,dj%p3))

		l=12
		do k=1,dj%p3
			do j=1,dj%p2
				do i=1,dj%p1
					read(42,rec=l) dj%buf(i,j,k)
					l=l+1
				end do
			end do
		end do

		close(42)
	end subroutine

	subroutine write_gomplt(fname,dj)
		character(*) :: fname
		type(gom_t) :: dj

		integer(4) :: i, j, k, l

		open(42,file=trim(fname),access='direct',recl=4)

		write(42,rec=1) dj%rank
		write(42,rec=2) dj%surface
		write(42,rec=3) dj%p3
		write(42,rec=4) dj%p2
		write(42,rec=5) dj%p1
		write(42,rec=6) dj%qmin(3)
		write(42,rec=7) dj%qmax(3)
		write(42,rec=8) dj%qmin(2)
		write(42,rec=9) dj%qmax(2)
		write(42,rec=10) dj%qmin(1)
		write(42,rec=11) dj%qmax(1)

		allocate(dj%buf(dj%p1,dj%p2,dj%p3))

		l=12
		do k=1,dj%p3
			do j=1,dj%p2
				do i=1,dj%p1
					write(42,rec=l) real(dj%buf(i,j,k))
					l=l+1
				end do
			end do
		end do

		close(42)
	end subroutine

	subroutine diffgom(gt1, gt2)
		type(gom_t) :: gt1, gt2
		
		logical :: err=.false.
		integer(4) :: i, j, k

		if (gt1%p1 /= gt2%p1) err=.true.
		if (gt1%p2 /= gt2%p2) err=.true.
		if (gt1%p3 /= gt2%p3) err=.true.

		do k=1,gt1%p3
			do j=1,gt1%p2
				do i=1,gt1%p1
					gt2%buf(i,j,k)=gt1%buf(i,j,k)-gt2%buf(i,j,k)
				end do
			end do
		end do
	end subroutine

end program
