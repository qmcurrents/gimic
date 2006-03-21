!
! $Id$
!

module integral_m
	use globals_m
	use grid_m
	use jfield_m
	use jtensor_m
	use gaussint_m
	use lip_m
	use teletype_m
	implicit none

	type integral_t
		type(jtensor_t), pointer :: jt
		type(jfield_t), pointer :: jf
		type(grid_t), pointer :: grid
		real(DP), dimension(3,3) :: tint
	end type

	public init_integral, del_integral, integral_t
	public integrate, int_t_direct, int_s_direct, int_t_2d, lipton
	public integral_master, integral_slave, write_integral
	
	private 
	
	integer(I4) :: nlip
contains
	subroutine init_integral(it, jt, jf, grid)
		type(integral_t) :: it
		type(jtensor_t), target :: jt
		type(jfield_t),  target :: jf
		type(grid_t), target :: grid
		
		it%jf=>jf
		it%jt=>jt
		it%grid=>grid
		it%tint=D0
	end subroutine

	subroutine del_integral(it)
		type(integral_t) :: it

		nullify(it%jf)
		nullify(it%jt)
		nullify(it%grid)
	end subroutine

	subroutine integrate(jf, grid)
		type(jfield_t), intent(inout) :: jf
		type(grid_t), intent(in) :: grid

		integer(I4) :: p1, p2, p3
		call get_grid_size(grid, p1, p2, p3)
		if (p3 == 1) then
			call integrate_s(jf,grid)
		else
			call integrate_t(jf,grid)
		end if
	end subroutine

	subroutine int_t_direct(it)
		type(integral_t) :: it

		integer(I4) :: i, j, k, p1, p2, p3
		real(DP), dimension(3) :: rr
		real(DP), dimension(3,3) :: xsum
		type(tensor_t), dimension(:), allocatable  :: jt1, jt2, jt3
		
        call eta(it%jt,it%grid)
		call get_grid_size(it%grid, p1, p2, p3)

		allocate(jt1(p1))
		allocate(jt2(p2))
		allocate(jt3(p3))

		do k=1,p3
			do j=1,p2
				do i=1,p1
					rr=gridpoint(it%grid, i, j, k)
					call jtensor(it%jt, rr, jt1(i))
				end do
				jt2(j)%t=int_t_1d(jt1,it%grid,1)
			end do
			jt3(k)%t=int_t_1d(jt2,it%grid,2)
		end do
		xsum=int_t_1d(jt3,it%grid,3)

!        call print_t_int(xsum)
		deallocate(jt1, jt2, jt3)
	end subroutine

	subroutine int_s_direct(it)
		type(integral_t), intent(inout) :: it

		integer(I4) :: i, j, k, p1, p2, p3
		real(DP), dimension(3) :: normal, rr, center, bb
		real(DP) :: psum, nsum, w, jp, r, bound
		real(DP) :: psum2, nsum2
		real(DP) :: psum3, nsum3
		real(DP) :: xsum, xsum2, xsum3
		type(vector_t) :: jvec
		type(tensor_t) :: jt
		
		call get_grid_size(it%grid, p1, p2, p3)
		call getkw('cdens.magnet', bb)

		normal=get_grid_normal(it%grid)

		bound=1.d+10
		call getkw('integral.radius', bound)
		call grid_center(it%grid,center)

		xsum3=0.d0
		psum3=0.d0
		nsum3=0.d0
		do k=1,p3
			xsum2=0.d0
			psum2=0.d0
			nsum2=0.d0
			do j=1,p2
				xsum=0.d0
				psum=0.d0
				nsum=0.d0
				do i=1,p1
					rr=gridpoint(it%grid, i, j, k)
					r=sqrt(sum((rr-center)**2))
					call jtensor(it%jt, rr, jt)
					jvec%v=matmul(jt%t,bb)
					if ( r > bound ) then
						w=0.d0
					else
						w=get_weight(it%grid, i, 1) 
						jp=dot_product(normal,jvec%v)
					end if
					xsum=xsum+jp*w
					if (jp > 0.d0) then
						psum=psum+jp*w
					else
						nsum=nsum+jp*w
					end if
				end do
				w=get_weight(it%grid,j,2)
				xsum2=xsum2+xsum*w
				psum2=psum2+psum*w
				nsum2=nsum2+nsum*w
			end do
			w=get_weight(it%grid,k,3)
			xsum3=xsum3+xsum2*w
			psum3=psum3+psum2*w
			nsum3=nsum3+nsum2*w
		end do

		call nl
		call msg_out(repeat('*', 60))
		write(str_g, '(a,f13.6)') '   Induced current (au)    :', xsum3
		call msg_out(str_g)
		write(str_g, '(a,f13.6,a,f11.6,a)') &
			'      Positive contribution:', psum3, '  (',au2si(psum3),' )'
		call msg_out(str_g)
		write(str_g, '(a,f13.6,a,f11.6,a)') &
			'      Negative contribution:', nsum3, '  (',au2si(nsum3),' )'
		call msg_out(str_g)
		call nl
		write(str_g, '(a,f13.6)') '   Induced current (nA/T)  :', au2si(xsum3)
		call msg_out(str_g)
		write(str_g, '(a,f13.6)') '      (conversion factor)  :', au2si(1.d0)
		call msg_out(str_g)
		call msg_out(repeat('*', 60))
	end subroutine

	subroutine int_t_2d_sigma(it, k)
		type(integral_t) :: it
		integer(I4), intent(in) :: k
		
		integer(I4) :: i, j, p1, p2
		real(DP), dimension(3) :: rr
		type(tensor_t), dimension(:), allocatable  :: jt1, jt2
		
		real(DP), dimension(3,3) :: s1
		real(DP) :: r3
		
		call get_grid_size(it%grid, p1, p2)

		allocate(jt1(p1))
		allocate(jt2(p2))

		do j=1,p2
			do i=1,p1
				rr=gridpoint(it%grid, i, j, k)
				call jtensor(it%jt, rr, jt1(i))
				r3=sqrt(sum(rr**2))**3

				if (r3 < 1.d-15) then
					jt1(i)%t=D0
					cycle
				end if
				s1=jt1(i)%t/r3
				
				s1(1,1)=rr(2)*s1(3,1)-rr(3)*s1(2,1)
				s1(1,2)=rr(2)*s1(3,2)-rr(3)*s1(2,2)
				s1(1,3)=rr(2)*s1(3,3)-rr(3)*s1(2,3)

				s1(2,1)=rr(3)*s1(1,1)-rr(1)*s1(3,1)
				s1(2,2)=rr(3)*s1(1,2)-rr(1)*s1(3,2)
				s1(2,3)=rr(3)*s1(1,3)-rr(1)*s1(3,3)

				s1(3,1)=rr(1)*s1(2,1)-rr(2)*s1(1,1)
				s1(3,2)=rr(1)*s1(2,2)-rr(2)*s1(1,2)
				s1(3,3)=rr(1)*s1(2,3)-rr(2)*s1(1,3)
				jt1(i)%t=s1
			end do
			jt2(j)%t=int_t_1d(jt1,it%grid,1)
		end do
		it%tint=int_t_1d(jt2,it%grid,2)

		deallocate(jt1, jt2)
	end subroutine

	subroutine int_t_2d(it, k)
		type(integral_t) :: it
		integer(I4), intent(in) :: k
		
		integer(I4) :: i, j, p1, p2
		real(DP), dimension(3) :: rr
		type(tensor_t), dimension(:), allocatable  :: jt1, jt2
		
		call get_grid_size(it%grid, p1, p2)

		allocate(jt1(p1))
		allocate(jt2(p2))

		do j=1,p2
			do i=1,p1
				rr=gridpoint(it%grid, i, j, k)
				call jtensor(it%jt, rr, jt1(i))
			end do
			jt2(j)%t=int_t_1d(jt1,it%grid,1)
		end do
		it%tint=int_t_1d(jt2,it%grid,2)

		deallocate(jt1, jt2)
	end subroutine

	function int_t_1d(jt, grid, axis) result(xsum)
		type(tensor_t), dimension(:), intent(in) :: jt
		type(grid_t), intent(in) :: grid
		integer(I4), intent(in) :: axis
		real(DP), dimension(3,3) :: xsum

		integer(I4) :: i, pts
		real(DP) :: w

		pts=size(jt)
		
		xsum=0.d0
		do i=1,pts
			w=get_weight(grid, i, axis) 
			xsum=xsum+jt(i)%t*w
		end do
	end function

	subroutine integrate_v(jf, grid)
		type(jfield_t), intent(inout) :: jf
		type(grid_t), intent(in) :: grid

		integer(I4) :: i, j, k, p1, p2, p3
		real(DP), dimension(3) :: rr
		real(DP) :: w
		real(DP), dimension(3) :: xsum, xsum2, xsum3
		type(vector_t), dimension(:,:), pointer :: jvec
		
		call get_grid_size(grid, p1, p2, p3)

		xsum3=0.d0
		do k=1,p3
			call get_jvectors(jf, k, jvec)
			xsum2=0.d0
			do j=1,p2
				xsum=0.d0
				do i=1,p1
					w=get_weight(grid, i, 1) 
					xsum=xsum+jvec(i,j)%v*w
				end do
				w=get_weight(grid,j,2)
				xsum2=xsum2+xsum*w
			end do
			w=get_weight(grid,k,3)
			xsum3=xsum3+xsum2*w
		end do

		call nl
		call msg_out(repeat('*', 60))
		write(str_g, '(a,3f13.6)') '   Induced current (au)    :', xsum3
		call msg_out(str_g)
		call nl
		write(str_g, '(a,3f13.6)') '   Induced current (nA/T)  :', &
			(au2si(xsum3(i)), i=1,3)
		call msg_out(str_g)
		call nl
		write(str_g, '(a,f13.6)') '        Conversion Factor   :', au2si(1.d0)
		call msg_out(str_g)
		call msg_out(repeat('*', 60))
	end subroutine

	subroutine integrate_t(jf, grid)
		type(jfield_t), intent(inout) :: jf
		type(grid_t), intent(in) :: grid

		integer(I4) :: i, j, k, p1, p2, p3
		real(DP), dimension(3) :: rr
		real(DP) :: w
		real(DP), dimension(3,3) :: xsum, xsum2, xsum3
		type(tensor_t), dimension(:,:), pointer :: jt
		
		call get_grid_size(grid, p1, p2, p3)

		xsum3=0.d0
		do k=1,p3
			call get_jtensors(jf,k,jt)
			xsum2=0.d0
			do j=1,p2
				xsum=0.d0
				do i=1,p1
					w=get_weight(grid, i, 1) 
					xsum=xsum+jt(i,j)%t*w
				end do
				w=get_weight(grid,j,2)
				xsum2=xsum2+xsum*w
			end do
			w=get_weight(grid,k,3)
			xsum3=xsum3+xsum2*w
		end do

!        call print_t_int(xsum3)
	end subroutine

	subroutine integrate_s(jf, grid)
		type(jfield_t), intent(inout) :: jf
		type(grid_t), intent(in) :: grid

		integer(I4) :: i, j, k, p1, p2, p3
		real(DP), dimension(3) :: normal, rr, center 
		real(DP) :: psum, nsum, w, jp, r, bound
		real(DP) :: psum2, nsum2
		real(DP) :: psum3, nsum3
		real(DP) :: xsum, xsum2, xsum3
		type(vector_t), dimension(:,:), pointer :: jvec
		
		call get_grid_size(grid, p1, p2, p3)

		normal=get_grid_normal(grid)

		bound=1.d+10
		call getkw('integral.radius', bound)
		call grid_center(grid,center)

		xsum3=0.d0
		psum3=0.d0
		nsum3=0.d0
		do k=1,p3
			call get_jvectors(jf, k, jvec)
			xsum2=0.d0
			psum2=0.d0
			nsum2=0.d0
			do j=1,p2
				xsum=0.d0
				psum=0.d0
				nsum=0.d0
				do i=1,p1
					rr=gridpoint(grid, i, j, k)
					r=sqrt(sum((rr-center)**2))
					if ( r > bound ) then
						w=0.d0
					else
						w=get_weight(grid, i, 1) 
						jp=dot_product(normal,jvec(i,j)%v)
					end if
					xsum=xsum+jp*w
					if (jp > 0.d0) then
						psum=psum+jp*w
					else
						nsum=nsum+jp*w
					end if
				end do
				w=get_weight(grid,j,2)
				xsum2=xsum2+xsum*w
				psum2=psum2+psum*w
				nsum2=nsum2+nsum*w
			end do
			w=get_weight(grid,k,3)
			xsum3=xsum3+xsum2*w
			psum3=psum3+psum2*w
			nsum3=nsum3+nsum2*w
		end do

		call nl
		call msg_out(repeat('*', 60))
		write(str_g, '(a,f13.6)') '   Induced current (au)    :', xsum3
		call msg_out(str_g)
		write(str_g, '(a,f13.6,a,f11.6,a)') &
			'      Positive contribution:', psum3, '  (',au2si(psum3),' )'
		call msg_out(str_g)
		write(str_g, '(a,f13.6,a,f11.6,a)') &
			'      Negative contribution:', nsum3, '  (',au2si(nsum3),' )'
		call msg_out(str_g)
		call nl
		write(str_g, '(a,f13.6)') '   Induced current (nA/T)  :', au2si(xsum3)
		call msg_out(str_g)
		write(str_g, '(a,f13.6)') '      (conversion factor)  :', au2si(1.d0)
		call msg_out(str_g)
		call msg_out(repeat('*', 60))
	end subroutine

	subroutine integrate_old(jf, grid)
		type(jfield_t), intent(inout) :: jf
		type(grid_t), intent(in) :: grid

		integer(I4) :: i, j, k, p1, p2, p3
		real(DP), dimension(3) :: normal, rr, center 
		real(DP), dimension(:), allocatable :: tmp, ptmp, ntmp
		real(DP) :: xsum, psum, nsum, w, jp, r, bound
		type(vector_t), dimension(:,:), pointer :: jvec
		
		call get_grid_size(grid, p1, p2, p3)
		allocate(tmp(p1))
		allocate(ptmp(p1))
		allocate(ntmp(p1))

		normal=get_grid_normal(grid)

		bound=1.d+10
		call getkw('integral.radius', bound)
		call grid_center(grid,center)

		call get_jvectors(jf, 1, jvec)
		do j=1,p2
			xsum=0.d0
			psum=0.d0
			nsum=0.d0
			do i=1,p1
				rr=gridpoint(grid, i, j, 1)
				r=sqrt(sum((rr-center)**2))
				if ( r > bound ) then
					w=0.d0
				else
					w=get_weight(grid, i, 1) 
                	jp=dot_product(normal,jvec(i,j)%v)
				end if
                xsum=xsum+jp*w
				if (jp > 0.d0) then
					psum=psum+jp*w
				else
					nsum=nsum+jp*w
				end if
			end do
			tmp(j)=xsum
			ptmp(j)=psum
			ntmp(j)=nsum
		end do
		
		xsum=0.d0
		psum=0.d0
		nsum=0.d0
		do j=1, p2
			w=get_weight(grid,j,2)
			xsum=xsum+tmp(j)*w
			psum=psum+ptmp(j)*w
			nsum=nsum+ntmp(j)*w
		end do

		deallocate(tmp, ptmp, ntmp)
		call nl
		call msg_out(repeat('*', 60))
		write(str_g, '(a,f13.6)') '   Induced current (au)    :', xsum
		call msg_out(str_g)
		write(str_g, '(a,f13.6,a,f11.6,a)') &
			'      Positive contribution:', psum, '  (',au2si(psum),' )'
		call msg_out(str_g)
		write(str_g, '(a,f13.6,a,f11.6,a)') &
			'      Negative contribution:', nsum, '  (',au2si(nsum),' )'
		call msg_out(str_g)
		call nl
		write(str_g, '(a,f13.6)') '   Induced current (nA/T)  :', au2si(xsum)
		call msg_out(str_g)
		write(str_g, '(a,f13.6)') '      (conversion factor)  :', au2si(1.d0)
		call msg_out(str_g)
		call msg_out(repeat('*', 60))
	end subroutine

	subroutine lipton(jf, g)
		type(jfield_t), intent(inout) :: jf
		type(grid_t), intent(inout) :: g

		integer(I4) :: i, j, p1, p2
		real(DP), dimension(3) :: normal, ljt
		integer(I4), dimension(3) :: ngp
		real(DP), dimension(:), allocatable :: tmp, gpts
		real(DP) :: xsum, w, jp, x, q
		type(vector_t), dimension(:,:), pointer :: jvec
		type(vector_t), dimension(:), allocatable :: z

		type(grid_t) :: lgrid
		

		!call getkw('integral.grid_points', g%npts)
		call getkw('integral.gauss_points', ngp)
		call getkw('integral.lip_order', nlip)
		call msg_info('lip_order kwyword ignored so far... using 3.')

		call copy_grid(g,lgrid)
		do i=1,3
        	call setup_lobby(0.d0, lgrid%l(i), ngp(i), lgrid%gdata(i))
		end do

		call get_grid_size(g, p1, p2)
        allocate(tmp(p2))
        allocate(gpts(p1))
        allocate(z(p1))

		normal=get_grid_normal(g)
		call get_jvectors(jf, 1, jvec)

		! dirt.
		gpts=g%gdata(1)%pts

		do j=1,p2
			do i=1,p1
				z(i)%v=jvec(i,j)%v !dirt. fix later.
			end do
			w=get_weight(lgrid, 1, 1)
			jp=dot_product(normal, z(1)%v)
			xsum=jp*w
			w=get_weight(lgrid, p1, 1)
			jp=dot_product(normal, z(p1)%v)
			xsum=xsum+jp*w
			do i=2,p1-1
				x=lgrid%gdata(1)%pts(i) !dirt.
				ljt=lip3n(x, gpts(i-1:i+1), z(i-1:i+1))
				jp=dot_product(normal, ljt)
				w=get_weight(lgrid, i, 1)
				xsum=xsum+jp*w
			end do
			tmp(j)=xsum
		end do
		
		! dirt.
		gpts=g%gdata(2)%pts

		w=get_weight(lgrid, 1, 2)
		xsum=w*tmp(1)
		w=get_weight(lgrid, p2, 2)
		xsum=xsum+w*tmp(p2)
		do i=2, p2-1
			x=lgrid%gdata(2)%pts(i) !dirt.
			q=lipn(x, gpts(i-1:i+1), tmp(i-1:i+1))
			w=get_weight(lgrid,i,2)
			xsum=xsum+tmp(i)*w
		end do

        deallocate(tmp, gpts, z)
        write(str_g, '(a,f)') 'Integrated current (LIP):', xsum
		call msg_out(str_g)
		call del_grid(lgrid)
	end subroutine

	function au2si(au) result(si)
		real(DP), intent(in) :: au
		real(DP) :: si
		
		real(DP) :: aulength, auspeedoflight, speedoflight, aucharge, hbar
		real(DP) :: autime, autesla, audjdb

		aulength=0.52917726D-10
		auspeedoflight=137.03599D0
		speedoflight=299792458.D0
		aucharge=1.60217733D-19
		hbar=1.05457267D-34

		autime=aulength*auspeedoflight/speedoflight
		autesla=hbar/aucharge/aulength/aulength
		audjdb=aucharge/autime/autesla

		si=au*audjdb*1.d+09 ! nA/T

!        write(6,*) 'The obtained conversion factors'
!        write(6,*) autime,' au time in seconds'
!        write(6,*) autesla,' au magnetic field in tesla'
!        write(6,*) audjdb*1.D+09,' au induced current in nanoampere/tesla'
	end function 

	subroutine print_t_int(xsum)
		real(DP), dimension(3,3) :: xsum
		
		integer(I4) :: i
		print *, xsum(1,:)
		print *, xsum(2,:)
		print *, xsum(3,:)
		
		call nl
		call msg_out(repeat('*', 70))
		write(str_g, '(a,3f13.6)') '   Induced current (au)    :', xsum(1,:)
		call msg_out(str_g)
		write(str_g, '(a,3f13.6)') '                            ', xsum(2,:)
		call msg_out(str_g)
		write(str_g, '(a,3f13.6)') '                            ', xsum(3,:)
		call msg_out(str_g)
		call nl
		write(str_g, '(a,3f13.6)') '   Induced current (nA/T)  :', &
			(au2si(xsum(1,i)), i=1,3)
		call msg_out(str_g)
		write(str_g, '(a,3f13.6)') '                            ', &
			(au2si(xsum(2,i)), i=1,3)
		call msg_out(str_g)
		write(str_g, '(a,3f13.6)') '                            ', &
			(au2si(xsum(3,i)), i=1,3)
		call msg_out(str_g)
		call nl
		write(str_g, '(a,f13.6)') '        Conversion Factor   :', au2si(1.d0)
		call msg_out(str_g)
		call msg_out(repeat('*', 70))

	end subroutine


	subroutine twinings
		return
	end subroutine

	subroutine integral_master(it, who)
		use mpi_m
		type(integral_t) :: it
		integer(I4), intent(in) :: who

#ifdef HAVE_MPI
		integer(I4), parameter :: BUFSZ=9
		integer(I4) :: jnum, ierr
		integer(I4), dimension(MPI_STATUS_SIZE) :: stat
		real(DP), dimension(3,3) :: buf
		real(DP) :: w

		call mpi_recv(jnum,1,MPI_INTEGER,who,INTGRL_TAG, &
				MPI_COMM_WORLD, stat, ierr)
		call mpi_recv(buf,BUFSZ,MPI_DOUBLE_PRECISION,who, &
				INTGRL_TAG, MPI_COMM_WORLD, stat, ierr)
			
		w=get_weight(it%grid, jnum, 3) 
		it%tint=it%tint+buf*w

#endif
	end subroutine

	subroutine integral_slave(it, jnum)
		use mpi_m
		type(integral_t) :: it
		integer(I4), intent(in) :: jnum

#ifdef HAVE_MPI
		integer(I4) :: ierr
		integer(I4), parameter :: INTBUFSZ=9
		integer(I4), dimension(MPI_STATUS_SIZE) :: stat
		real(DP), dimension(3,3) :: buf

		buf=D0

		call int_t_2d(it, jnum)
		call mpi_send(INTGRL_TAG, 1, MPI_INTEGER, 0, JOB_TAG,&
			MPI_COMM_WORLD, ierr)
		call mpi_send(jnum, 1, MPI_INTEGER, 0, INTGRL_TAG,&
			MPI_COMM_WORLD, ierr)
		call mpi_send(it%tint, INTBUFSZ, MPI_DOUBLE_PRECISION,0, INTGRL_TAG, &
			MPI_COMM_WORLD, ierr)
#endif
	end subroutine
	
	subroutine write_integral(it)
		type(integral_t) :: it
		integer(I4) :: tun

       	call print_t_int(it%tint)
		open(JINTFD, file='JINTEGRAL')
		call get_teletype_unit(tun)
		call set_teletype_unit(JINTFD)
       	call print_t_int(it%tint)
		call set_teletype_unit(tun)
		close(JINTFD)
	end subroutine
end module
