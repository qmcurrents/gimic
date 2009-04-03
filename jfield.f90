! This module calulates the current field by looping over gridpoints and 
! calling the appropriate routines in jtensor_class
!
! Ok, this is really a mess...
! This modules needs a seriuos clean up. I have been adding, and adding, and
! adding stuff and it's not pretty.
!

module jfield_class
	use globals_m
	use jtensor_class
	use grid_class
	use basis_class
	use teletype_m
	use parallel_m
	use magnet_m
	implicit none

	type jfield_t
		type(jtensor_t), pointer :: jt
		real(DP), dimension(3) :: b
		type(tensor_t), dimension(:,:), pointer :: jj, jja, jjb, jsd
		type(vector_t), dimension(:,:), pointer :: vv, vva, vvb, vsd
		integer(I4) :: zpos_j, zpos_v
		type(grid_t), pointer :: grid
	end type

	public init_jfield, del_jfield, jfield, jvectors, jfield_direct
	public jfield_t, get_jvectors, get_jtensors, jvector_plot
	public set_jtensors
	
	private

	real(DP) :: vec_scale=D1
	integer(I4) :: spin
	character(BUFLEN) :: jtensor_file, jvec_file
	character(BUFLEN) :: jmod_plt, jvec_plt, njvec_plt, jprj_plt
contains
	subroutine init_jfield(f, jt, g)
		type(jfield_t) :: f
		type(jtensor_t), target :: jt
		type(grid_t), target :: g

		integer(I4) :: i, j, k

		f%jt=>jt
		f%grid=>g

		call push_section(input, 'cdens')
		if (master_p) then
			call setup_files(f)
		end if

		call get_grid_size(g,i,j,k)
		allocate(f%jj(i,j))
		allocate(f%vv(i,j))
		if (uhf_p) then
			allocate(f%jja(i,j))
			allocate(f%jjb(i,j))
			allocate(f%jsd(i,j))
			allocate(f%vva(i,j))
			allocate(f%vvb(i,j))
			allocate(f%vsd(i,j))
		end if

		call getkw(input, 'scale_vectors', vec_scale)
		if ( vec_scale /= D1 ) then
			write(str_g, '(a,f7.3)') 'Vector scaling =', vec_scale
			call msg_note(str_g)
		end if
		call nl
		call pop_section(input)

		call push_section(input, 'cdens')
        call get_magnet(g, f%b)
		call pop_section(input)
	end subroutine

	subroutine setup_files(jf)
		type(jfield_t) :: jf
	
		logical :: plot_p
		integer(I4) :: jtrecl, jvrecl

		jtrecl=9*DP
		jvrecl=3*DP
		
		jtensor_file=''
		call getkw(input, 'jtensor', jtensor_file)
        if (jtensor_file == '') then
            call msg_error('setup_files(): No tensor file specified!')
            stop
        end if
		call getkw(input, 'jvector', jvec_file)
        if (jvec_file == '') then
            call msg_error('setup_files(): No vector file specified!')
            stop
        end if
        open(JTFD, file=jtensor_file, access='direct', recl=jtrecl)
		open(JVECFD, file=jvec_file, access='direct', recl=jvrecl)
		if (uhf_p) then
			open(JTFD+1, file=trim(jtensor_file) // '_A', &
			access='direct', recl=jtrecl)
			open(JTFD+2, file=trim(jtensor_file) // '_B', &
			access='direct', recl=jtrecl)
			open(JTFD+3, file=trim(jtensor_file) // '_SD', &
			access='direct', recl=jtrecl)
			open(JVECFD+1, file=trim(jvec_file) // '_A', &
			access='direct', recl=jvrecl)
			open(JVECFD+2, file=trim(jvec_file) // '_B', &
			access='direct', recl=jvrecl)
			open(JVECFD+3, file=trim(jvec_file) // '_SD', &
			access='direct', recl=jvrecl)
		end if

	end subroutine

	subroutine del_jfield(f)
		type(jfield_t), intent(inout) :: f
		
		f%b=0.d0
		deallocate(f%jj)
		deallocate(f%vv)
		nullify(f%jj)
		nullify(f%vv)
		if (uhf_p) then
			deallocate(f%jja)
			deallocate(f%jjb)
			deallocate(f%jsd)
			deallocate(f%vva)
			deallocate(f%vvb)
			deallocate(f%vsd)
		end if
		
		if (master_p) then
			close(JTFD)
			close(JVECFD)
			if (uhf_p) then
			close(JTFD+1)
			close(JTFD+2)
			close(JTFD+3)
			close(JVECFD+1)
			close(JVECFD+2)
			close(JVECFD+3)
			end if
		end if
	end subroutine
	
	subroutine jvectors(jf)
		type(jfield_t), intent(inout) :: jf

		integer(I4) :: i, j, k, p1, p2, p3

		call get_grid_size(jf%grid, p1, p2, p3)

		call msg_note( 'Contracting j-tensors with magnetic field.')
		call msg_note( 'Writing j-vectors in: ' // trim(jvec_file))
		call nl

		do k=1,p3
			call jtens_io(jf, k, 'r')
			jf%zpos_j=k
			do j=1,p2
				do i=1,p1
					jf%vv(i,j)%v=matmul(jf%jj(i,j)%t, jf%b)
					if (uhf_p) then
						jf%vva(i,j)%v=matmul(jf%jja(i,j)%t, jf%b)
						jf%vvb(i,j)%v=matmul(jf%jjb(i,j)%t, jf%b)
						jf%vsd(i,j)%v=matmul(jf%jsd(i,j)%t, jf%b)
					end if
				end do
			end do
			call jvec_io(jf, k, 'w')
			jf%zpos_v=k
		end do

	end subroutine

	subroutine jtens_io(jf, p3, op)
		type(jfield_t), intent(inout) :: jf
		integer(I4), intent(in) :: p3
		character(1), intent(in) :: op

		integer(I4) :: i, j, k, p1, p2 

		call get_grid_size(jf%grid, p1, p2)

		k=p1*p2*(p3-1)+1
		do j=1, p2
			do i=1, p1
				select case(op)
				case('r')
					read(JTFD, rec=k) jf%jj(i,j)%t
					if (uhf_p) then
						read(JTFD+1, rec=k) jf%jja(i,j)%t
						read(JTFD+2, rec=k) jf%jjb(i,j)%t
						read(JTFD+3, rec=k) jf%jsd(i,j)%t
					end if
				case('w')
					write(JTFD, rec=k) jf%jj(i,j)%t
					if (uhf_p) then
						write(JTFD+1, rec=k) jf%jja(i,j)%t
						write(JTFD+2, rec=k) jf%jjb(i,j)%t
						write(JTFD+3, rec=k) jf%jsd(i,j)%t
					end if
				case default
					call msg_error('jtens_io(): Invalid operation: ' // op)
					call exit(1)
				end select
				k=k+1
			end do
		end do
	end subroutine

	subroutine jvec_io(jf, p3, op)
		type(jfield_t), intent(inout) :: jf
		integer(I4), intent(in) :: p3
		character(1), intent(in) :: op

		integer(I4) :: i, j, k, p1, p2 

		call get_grid_size(jf%grid, p1, p2)

		k=p1*p2*(p3-1)+1
		do j=1,p2
			do i=1,p1
				select case(op)
				case('r')
					read(JVECFD, rec=k) jf%vv(i,j)%v
					if (uhf_p) then
						read(JVECFD+1, rec=k) jf%vva(i,j)%v
						read(JVECFD+2, rec=k) jf%vvb(i,j)%v
						read(JVECFD+3, rec=k) jf%vsd(i,j)%v
					end if
				case('w')
					write(JVECFD, rec=k) jf%vv(i,j)%v
					if (uhf_p) then
						write(JVECFD+1, rec=k) jf%vva(i,j)%v
						write(JVECFD+2, rec=k) jf%vvb(i,j)%v
						write(JVECFD+3, rec=k) jf%vsd(i,j)%v
					end if
				case default
					call msg_error('jvec_io(): Invalid operation: ' // op)
					call exit(1)
				end select
				k=k+1
			end do
		end do
	end subroutine

	subroutine wrt_jvec(rr,v,fd)
		real(DP), dimension(:), intent(in) :: rr, v
		integer(I4), intent(in) :: fd

		if (fd == 0) return

		write(fd, '(6f11.7)')  rr, v
	end subroutine	

	subroutine wrt_njvec(rr,v,fd)
		real(DP), dimension(:), intent(in) :: rr, v
		integer(I4), intent(in) :: fd

		real(DP) :: nfac

		if (fd == 0) return

		nfac=sqrt(sum(v(:)**2))
		if (nfac < 1.d-15) then
			write(fd, '(6f11.7)')  rr, 0.d0, 0.d0, 0.d0
		else 
			write(fd, '(6f11.7)')  rr, v/nfac
		end if
	end subroutine	

	subroutine wrt_jproj(rr,v,grid,fd)
		real(DP), dimension(:), intent(in) :: rr, v
		type(grid_t) :: grid
		integer(I4), intent(in) :: fd

		real(DP) :: jprj
		real(DP), dimension(3) :: norm

		if (fd == 0) return
		norm=get_grid_normal(grid)

		jprj=dot_product(norm,v)
		write(fd, '(3e19.12)') rr, jprj
	end subroutine	

	subroutine wrt_jmod(rr,v,fd)
		real(DP), dimension(:), intent(in) :: rr, v
		integer(I4), intent(in) :: fd

		real(DP) :: jmod

		if (fd == 0) return

		jmod=sqrt(sum(v**2))
		write(fd, '(6f11.7)')  rr, jmod
	end subroutine	

	function getvecs(jf, spin) result(jv)
		type(jfield_t), intent(inout) :: jf
		integer(I4) :: spin

		type(vector_t), dimension(:,:), pointer :: jv
		select case(spin)
			case(1)
				jv=>jf%vv
			case(2)
				jv=>jf%vva
			case(3)
				jv=>jf%vvb
			case(4)
				jv=>jf%vsd
		end select
	end function

	function open_plot(kname, spin) result(fd)
		character(*), intent(in) :: kname
		integer(I4), intent(in) :: spin

		integer(I4) :: fd
		character(132) :: fname

		fd=0
		if (.not.master_p) return
		call getkw(input, kname, fname)
		if (trim(fname) == '') return

		call getfd(fd)
		if (fd == 0) then
			stop 1
		end if

		select case(spin)
			case(1)
				fname=trim(fname)//'.txt'
			case(2)
				fname=trim(fname)//'_a'//'.txt'
			case(3)
				fname=trim(fname)//'_b'//'.txt'
			case(4)
				fname=trim(fname)//'_sd'//'.txt'
		end select

		write(str_g, *) 'Writing ', trim(kname), ' in ', trim(fname)
		call msg_note(str_g)
		open(fd,file=trim(fname),status='unknown')

		return 

	end function

	subroutine jvector_plot(jf)
		type(jfield_t), intent(inout) :: jf

		integer(I4) :: i, j, k, p1, p2, p3, spin, ispin
		integer(I4) :: fd1, fd2, fd3, fd4
		real(DP), dimension(3) :: v, rr
		type(vector_t), dimension(:,:), pointer :: jv

		spin=1
		if (uhf_p) then
			spin=4
		end if

		call push_section(input, 'cdens.plot')
		call get_grid_size(jf%grid, p1, p2, p3)

		do ispin=1,spin
			fd1= open_plot('vector',ispin)
			fd2= open_plot('modulus',ispin)
			fd3= open_plot('projection',ispin)
			fd4= open_plot('nvector',ispin)

			call jvec_io(jf, 1, 'r')
			jv=>getvecs(jf,ispin)
			jf%zpos_v=1
			do j=1,p2
				do i=1,p1
					rr=gridpoint(jf%grid, i,j,1)*AU2A
					v=jv(i,j)%v*AU2A
					call wrt_jvec(rr,v,fd1)
					call wrt_jmod(rr,v,fd2)
					call wrt_jproj(rr,v,jf%grid,fd3)
					call wrt_njvec(rr,v,fd4)
				end do
				if (fd1 /= 0) write(fd1, *)
				if (fd2 /= 0) write(fd2, *)
				if (fd3 /= 0) write(fd3, *)
				if (fd4 /= 0) write(fd4, *)
			end do

			call closefd(fd1)
			call closefd(fd2)
			call closefd(fd3)
			call closefd(fd4)
			call jmod_cubeplot(jf,ispin)
		end do
		call jmod_gopenmol(jf)
		call pop_section(input)
	end subroutine

	subroutine jfield(jf)
		type(jfield_t) :: jf

		integer(I4) :: i, j, k, p1, p2, p3
		real(DP), dimension(3) :: rr
		character(10) :: op
		integer(I4) :: lo, hi, npts

		call eta(jf%jt,jf%grid)
		call get_grid_size(jf%grid, p1, p2, p3)

		call schedule(p2, lo, hi)
	
		do k=1,p3
			do j=lo,hi
				do i=1,p1
					rr=gridpoint(jf%grid, i, j, k)
					call ctensor(jf%jt, rr, jf%jj(i,j), 'total')
					if (uhf_p) then
						call ctensor(jf%jt, rr, jf%jja(i,j), 'alpha')
						call ctensor(jf%jt, rr, jf%jjb(i,j), 'beta')
						call ctensor(jf%jt, rr, jf%jsd(i,j), 'spindens')
					end if
				end do
			end do
			call gather_data(jf%jj,jf%jj(:,lo:hi))
			if (uhf_p) then
				call gather_data(jf%jja,jf%jja(:,lo:hi))
				call gather_data(jf%jjb,jf%jjb(:,lo:hi))
				call gather_data(jf%jsd,jf%jsd(:,lo:hi))
			end if
			if (master_p) call jtens_io(jf, k, 'w')
			jf%zpos_j=k
		end do
	end subroutine 

	subroutine jfield_direct(jf, k)
		type(jfield_t) :: jf
		integer(I4), intent(in) :: k

		integer(I4) :: i, j, p1, p2
		integer(I4) :: lo, hi, npts
		real(DP), dimension(3) :: rr
		type(tensor_t), dimension(:,:), pointer :: jj
		
		call get_grid_size(jf%grid, p1, p2)
		jj=>jf%jj

		call schedule(p2, lo, hi)

		do j=lo,hi
			do i=1,p1
				rr=gridpoint(jf%grid, i, j, k)
				call ctensor(jf%jt, rr, jf%jj(i,j), 'total')
				if (uhf_p) then
					call ctensor(jf%jt, rr, jf%jja(i,j), 'alpha')
					call ctensor(jf%jt, rr, jf%jjb(i,j), 'beta')
					call ctensor(jf%jt, rr, jf%jsd(i,j), 'spindens')
				end if
			end do
		end do
		call gather_data(jf%jj,jf%jj(:,lo:hi))
		if (uhf_p) then
			call gather_data(jf%jja,jf%jja(:,lo:hi))
			call gather_data(jf%jjb,jf%jjb(:,lo:hi))
			call gather_data(jf%jsd,jf%jsd(:,lo:hi))
		end if
	end subroutine 

	subroutine get_jtensors(jf,z,jt)
		type(jfield_t), intent(inout) :: jf
		integer(I4), intent(in) :: z
		type(tensor_t), dimension(:,:), pointer :: jt

		if (jf%zpos_j /= z) then
			call jtens_io(jf, z, 'r')
		end if

		jt=>jf%jj
	end subroutine 

	subroutine set_jtensors(jf,z,jt)
		type(jfield_t), intent(inout) :: jf
		integer(I4), intent(in) :: z
		type(tensor_t), dimension(:,:) :: jt

		jf%jj=jt
		call jtens_io(jf, z, 'w')
		jf%zpos_j=z
	end subroutine 

	subroutine get_jvectors(jf,z,jvec)
		type(jfield_t), intent(inout) :: jf
		integer(I4), intent(in) :: z
		type(vector_t), dimension(:,:), pointer :: jvec

		if (jf%zpos_v /= z) then
			call jvec_io(jf, z, 'r')
		end if

		jvec=>jf%vv
	end subroutine 

	subroutine print_jt(rr, jt)
		real(DP), dimension(3), intent(in) :: rr
		real(DP), dimension(3,3), intent(in) :: jt

		integer(I4) :: l
		print *
		print '(a,3f12.8)', 'Current tensor at ', rr
		print *, '======================================='
		print '(3e19.12)', (jt(l,:), l=1,3)
		print *
		print '(a,e19.12)', ' Trace:', jt(1,1)+jt(2,2)+jt(3,3)
	end subroutine

	subroutine jmod_cubeplot(jf, spin)
		type(jfield_t) :: jf
		integer(I4), intent(in) :: spin

		integer(I4) :: p1, p2, p3, fd1, fd2, fd3
		integer(I4) :: i, j, k, l
		real(DP), dimension(3) :: qmin, qmax
		real(DP), dimension(3) :: norm, step
		real(DP) :: maxi, mini, val, sgn 
		integer(I4), dimension(3) :: npts
		type(vector_t), dimension(:,:), pointer :: buf
		character(BUFLEN) :: fname

		buf=>jf%vv
		fname=''

		call get_grid_size(jf%grid, p1, p2, p3)
		npts=(/p1,p2,p3/)
		norm=get_grid_normal(jf%grid)
		qmin=gridpoint(jf%grid,1,1,1)*AU2A
		qmax=gridpoint(jf%grid,p1,p2,p3)*AU2A

		step=(qmax-qmin)/(npts-1)

		call getkw(input, 'cube_mod', fname)
		fd1=opencube(trim(fname), spin, qmin, step, npts)
!        call getkw(input, 'cube_mordi', fname)
!        fd2=opencube(trim(fname)//'+', spin, qmin, step, npts)
!        fd3=opencube(trim(fname)//'-', spin, qmin, step, npts)


		maxi=0.d0
		mini=0.d0
		l=0
		do k=1,p3
			call jvec_io(jf, k, 'r')
			do j=1,p2
				do i=1,p1
					sgn=dot_product(norm,buf(i,j)%v)
					val=(sqrt(sum(buf(i,j)%v**2)))
					if (val > maxi) maxi=val
					if (val < mini) mini=val

					if (fd1 /= 0) then 
						write(fd1,'(f12.6)',advance='no') val
						if (mod(l,6) == 5) write(fd1,*)
					end if
!                    if (fd2 /= 0) then 
!                        if (sgn > 0 ) then
!                            write(fd2,'(f12.6)',advance='no') val
!                        else
!                            write(fd2,'(f12.6)',advance='no') 0.d0
!                        end if
!                        if (mod(l,6) == 5) write(fd2,*)
!                    end if
!                    if (fd3 /= 0) then 
!                        if (sgn < 0 ) then
!                            write(fd3,'(f12.6)',advance='no') -1.d0*val
!                        else
!                            write(fd3,'(f12.6)',advance='no') 0.d0
!                        end if
!                        if (mod(l,6) == 5) write(fd3,*)
!                    end if
					l=l+1
				end do
			end do
		end do
		print *, 'maximini:', maxi, mini
	end subroutine

	function opencube(kname, spin, origin, step, npts) result(fd)
		character(*), intent(in) :: kname
		integer(I4), intent(in) :: spin
		real(DP), dimension(3), intent(in) :: origin, step
		integer(I4), dimension(3), intent(in) :: npts

		integer(I4) :: fd
		character(132) :: fname

		fd=0
		fname=kname
		if (.not.master_p) return
		if (trim(fname) == '') return

		call getfd(fd)
		if (fd == 0) then
			stop 1
		end if

		select case(spin)
			case(1)
				fname=trim(fname)//'.cube'
			case(2)
				fname=trim(fname)//'a'//'.cube'
			case(3)
				fname=trim(fname)//'b'//'.cube'
			case(4)
				fname=trim(fname)//'sd'//'.cube'
		end select

		write(str_g, *) 'Writing ', trim(kname), ' in ', trim(fname)
		call msg_note(str_g)

		open(fd,file=trim(fname),form='formatted',status='unknown')
		write(fd,*) 'Gaussian cube data, generated by genpot'
		write(fd,*) 
		write(fd, '(i5,3f12.6)') 0, origin
		write(fd, '(i5,3f12.6)') npts(1), step(1), 0.d0, 0.d0
		write(fd, '(i5,3f12.6)') npts(2), 0.d0, step(2), 0.d0
		write(fd, '(i5,3f12.6)') npts(3), 0.d0, 0.d0, step(3)

		return 

	end function

	subroutine jmod_gopenmol(jf)
		type(jfield_t) :: jf

		integer(I4) :: surface, rank, p1, p2, p3
		integer(I4) :: i, j, k, l
		real(SP), dimension(3) :: qmin, qmax
		real(DP), dimension(3) :: norm
		real(SP) :: maxi, mini, val
		type(vector_t), dimension(:,:), pointer :: buf
		character(BUFLEN) :: gopen_file

		buf=>jf%vv
		gopen_file=''
		call getkw(input, 'cdens.plot.gopenmol', gopen_file)
		if (trim(gopen_file) == '') return
		open(GOPFD,file=trim(gopen_file),access='direct',recl=1)
		open(GOPFD2,file=trim('prj_'//gopen_file),access='direct',recl=1)

		surface=200
		rank=3

		call get_grid_size(jf%grid, p1, p2, p3)
		norm=get_grid_normal(jf%grid)
		qmin=real(gridpoint(jf%grid,1,1,1)*AU2A)
		qmax=real(gridpoint(jf%grid,p1,p2,p3)*AU2A)

		write(GOPFD,rec=1) rank;      write(GOPFD2,rec=1) rank
		write(GOPFD,rec=2) surface;   write(GOPFD2,rec=2) surface
		write(GOPFD,rec=3) p3;        write(GOPFD2,rec=3) p3
		write(GOPFD,rec=4) p2;        write(GOPFD2,rec=4) p2
		write(GOPFD,rec=5) p1;        write(GOPFD2,rec=5) p1
		write(GOPFD,rec=6) qmin(3);   write(GOPFD2,rec=6) qmin(3)
		write(GOPFD,rec=7) qmax(3);   write(GOPFD2,rec=7) qmax(3)
		write(GOPFD,rec=8) qmin(2);   write(GOPFD2,rec=8) qmin(2)
		write(GOPFD,rec=9) qmax(2);   write(GOPFD2,rec=9) qmax(2)
		write(GOPFD,rec=10) qmin(1);  write(GOPFD2,rec=10) qmin(1)
		write(GOPFD,rec=11) qmax(1);  write(GOPFD2,rec=11) qmax(1)

		print *, rank;     
		print *, surface;  
		print *, p3;       
		print *, p2;       
		print *, p1;       
		print *, qmin(3);  
		print *, qmax(3);  
		print *, qmin(2);  
		print *, qmax(2);  
		print *, qmin(1); 
		print *, qmax(1); 

		maxi=0.d0
		mini=0.d0
		l=12
		do k=1,p3
			call jvec_io(jf, k, 'r')
			do j=1,p2
				do i=1,p1
					val=real(sqrt(sum(buf(i,j)%v**2)))
					write(GOPFD,rec=l) val
					if (val > maxi) maxi=val
					if (val < mini) mini=val
					write(GOPFD2,rec=l) real(dot_product(norm,buf(i,j)%v))
					l=l+1
				end do
			end do
		end do
		print *, 'maximini:', maxi, mini

		close(GOPFD)
		close(GOPFD2)
	end subroutine

!    subroutine read_jdata(f, g) 
!        type(jfield_t), intent(inout) :: f
!        type(grid_t), intent(inout) :: g
!
!        integer(I4) :: p1, p2, p3, i, j, k, foo_p 
!        
!        open(JTFD, file=jtensor_file, status='old')
!        print *, '*** Reading grid data from: ', xtrim(jtensor_file)
!        print *
!        call read_grid(g,JTFD)
!        call get_grid_size(g,p1,p2,p3)
!
!        allocate(f%j3(p1,p2,p3))
!        allocate(f%v3(p1,p2,p3))
!        
!        do k=1,p3
!            do j=1,p2
!                do i=1,p1
!                    read(JTFD, *) f%j3(i,j,k)%r
!                    read(JTFD, *) f%j3(i,j,k)%j
!                end do
!            end do
!        end do
!        close(JTFD)
!    end subroutine
end module
