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

		plot_p=.false.
		jvec_file =''
		jmod_plt=''
		jvec_plt=''
		njvec_plt=''
		jprj_plt=''
		call getkw(input, 'plot', plot_p)
		if (plot_p) then
			call getkw(input, 'plot.vector', jvec_plt)
			call getkw(input, 'plot.nvector', njvec_plt)
			call getkw(input, 'plot.modulus', jmod_plt)
			call getkw(input, 'plot.projection', jprj_plt)
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

	subroutine jvector_plot(jf)
		type(jfield_t), intent(inout) :: jf

		integer(I4) :: i, j, k, p1, p2, p3, spin, ispin
		integer(I4), dimension(:), pointer :: z
		real(DP), dimension(3) :: foo, norm, rr
		real(DP), dimension(2) :: mr
		real(DP) :: jmod, jprj, nfac
		character(BUFLEN+20+5) :: tfil
		character(2) :: scase

		call get_grid_size(jf%grid, p1, p2, p3)
		norm=get_grid_normal(jf%grid)

		spin=1
		if (uhf_p) then
			spin=4
		end if

		do ispin=1,spin
			select case(ispin)
				case(1)
					scase=''
				case(2)
					scase='a'
				case(3)
					scase='b'
				case(4)
					scase='sd'
			end select
			if (jvec_plt /= '') then
				tfil=trim(jvec_plt)//trim(scase)//'.txt' 
				call msg_note('Writing j in        : '//trim(tfil))
				open(JVPFD, file=trim(tfil))
				write(JVPFD, *) '#********************************************#'
				write(JVPFD, *) '#                                            #'
				write(JVPFD, *) '#    GIMIC - THE HELLO WORLD PROGRAM         #'
				write(JVPFD, *) '#                                            #'
				write(JVPFD, *) '#********************************************#'
			end if
			if (jmod_plt /= '') then
				tfil=trim(jmod_plt)//trim(scase)//'.txt'
				call msg_note('Writing |j| in      : '//trim(tfil))
				open(MODFD, file=trim(tfil))
			end if
			if (njvec_plt /= '') then
				tfil=trim(njvec_plt)//trim(scase)//'.txt'
				call msg_note('Writing j/|j| in    : ' //trim(tfil))
				open(NJVFD, file=trim(tfil))
			end if
			if (jprj_plt /= '') then
				tfil=trim(jprj_plt)//trim(scase)//'.txt'
				call msg_note('Writing (j;n) in    : '//trim(tfil))
				open(JPRJFD, file=trim(tfil))
			end if

			call jvec_io(jf, 1, 'r')
			jf%zpos_v=1
			do j=1,p2
				do i=1,p1
					rr=gridpoint(jf%grid, i,j,1)*AU2A
					select case(ispin)
						case(1)
							foo=jf%vv(i,j)%v
						case(2)
							foo=jf%vva(i,j)%v
						case(3)
							foo=jf%vvb(i,j)%v
						case(4)
							foo=jf%vsd(i,j)%v
					end select
					foo=foo*AU2A
					nfac=sqrt(sum(foo(:)**2))
					jmod=sqrt(sum(foo**2))
					jprj=dot_product(norm,foo)
					if (jmod_plt /= '') then
						write(MODFD, '(4e18.10)') rr, jmod
					end if
					if (jvec_plt /= '') then
						write(JVPFD, '(6f11.7)')  rr, foo
					end if
					if (njvec_plt /= '') then
						if (nfac < 1.d-15) then
							write(NJVFD, '(6f11.7)')  rr, 0.d0, 0.d0, 0.d0
						else 
							write(NJVFD, '(6f11.7)')  rr, foo/nfac
						end if
					end if
					if (jprj_plt /= '') then
						write(JPRJFD, '(3e19.12)') rr, jprj
					end if
				end do
				if (jvec_plt /= '') then
					write(JVPFD,*)
				end if
				if (jmod_plt /= '') then
					write(MODFD,*)
				end if
				if (jprj_plt /= '') then
					write(JPRJFD,*)
				end if
			end do

			if (jmod_plt /= '') close(MODFD)
			if (jvec_plt /= '') close(JVPFD)
			if (njvec_plt /= '') close(NJVFD)
			if (jprj_plt /= '')  close(JPRJFD)
			call nl
		end do
		call jmod_gopenmol(jf)
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
