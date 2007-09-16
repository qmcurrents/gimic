! 
! GIMIC - a pretty advanced 'Hello World!' program.
!

! keeping the global stuff to a bare minimum... :)
program gimic
	use globals_m  
	use teletype_m
	use basis_m
	use timer_m
    implicit none 
	character(BUFLEN) :: buf
	type(molecule_t) :: mol
	
	call set_debug_level(3)

	call new_getkw(input)

	call getenv('MPIRUN', buf)
	if ( trim(buf) == '1' ) mpirun_p=.true.
	call getkw(input,'mpirun', mpirun_p)
	if (.not.mpirun_p) call program_header

	call initialize()
	call cdens()
	call finalize()

	call stockas_klocka()
	call msg_out('Hello World! (tm)')
!    call program_footer()
	
	if (bert_is_evil) then
		call msg_error('Bert is evil, and your results are wicked.')
		call nl
	end if

contains

	subroutine read_inpbuf(inpfile, inpbuf)
		character(*), intent(in) :: inpfile
		character, dimension(:), pointer :: inpbuf

		external fgetc
		integer(I4) :: i, n, err
		integer(I4) :: fgetc 

		write(str_g, 100) "Reading input from '", trim(inpfile),"'"
		call msg_info(str_g)
		open(99, file=trim(inpfile), status='old', action='read', iostat=err)

		if (err /= 0 ) then
			write(str_g, 100) "read_inpbuf(): open failed for '",&
				trim(inpfile),"'"
			call msg_error(str_g)
			stop
		end if

		n=getfsize(99)
		allocate(inpbuf(n))
		do i=1,n
			err=fgetc(99, inpbuf(i))
		end do
		close(99)
100 format(a,a,a)
	end subroutine

	subroutine read_molbuf(molbuf)
!
! read in the MOL file into a buffer which can be broadcast 
!
		character(MAX_LINE_LEN), dimension(:), pointer :: molbuf
		integer(I4) :: i, n, ios
		character(BUFLEN) :: basfile
		call getkw(input, 'basis', basfile)
		call msg_note('Basis sets read from file: ' // trim(basfile))
		call nl

		open(BASFD, file=basfile, status='old', action='read', &
			 form='formatted', iostat=ios)

		if (ios /= 0 ) then
			call msg_error('read_intgrl(): open failed.')
			stop
		end if

		n=getnlines(BASFD)
		allocate(molbuf(n))
		do i=1,n
			read(BASFD, '(200a)') molbuf(i) ! MAX_LINE_LEN
		end do
		close(BASFD)

	end subroutine

	subroutine initialize()
		use parallel_m

		external fdate, hostnm
		integer(I4) :: i, hostnm, rank, ierr
		character(BUFLEN) :: title, fdate, sys
		real(DP), dimension(3) :: center
!        character, dimension(:), pointer :: inpbuf
!        character(MAX_LINE_LEN), dimension(:), pointer :: molbuf

		if (mpirun_p) then
			rank=start_mpi()
		else
			master_p=.true.
			rank=0
		end if

!        nullify(inpbuf)
!        nullify(molbuf)

		ierr=hostnm(sys)
!        call read_inpbuf('standard input', inpbuf)
		if (master_p) then
			write(str_g, '(a,i3,a,a)') 'I''m master', rank,' on ', trim(sys)
			call msg_debug(str_g,2)
		else
            call set_teletype_unit(DEVNULL)
			write(str_g, '(a,i3,a,a)') 'I''m puppet', rank,' on ', trim(sys)
			call msg_debug(str_g,2)
		end if

!        if (mpirun_p) call bcast_inpbuf(inpbuf)

!        call init_parse_buf(inpbuf, ierr)
		if (ierr == 2) then
			call del_getkw(input)
			call msg_error('cdens(): parse error')
			call exit(1)
		end if
		
		i=0
		call getkw(input, 'debug', i)
		call set_debug_level(i)

		call msg_out(fdate())
		
		call getkw(input, 'title', title)
		call msg_out(' TITLE: '// trim(title))
		call nl

		spherical=.false.
		call getkw(input, 'spherical', spherical)
		
!        if (master_p) then
!            call read_molbuf(molbuf)
!        end if
!        if (mpirun_p) call bcast_molbuf(molbuf)

!        call init_basis(mol, molbuf)
        call init_basis(mol)

!        deallocate(inpbuf)
!        deallocate(molbuf)

	end subroutine

	subroutine finalize()
		use parallel_m
		call del_basis(mol)
		call stop_mpi()
		call del_getkw(input)
	end subroutine

	subroutine setup_grids(calc, grid, igrid, dgrid, egrid)
		use grid_m
		integer(I4), dimension(:) :: calc
		type(grid_t) :: grid, igrid, dgrid, egrid
		logical :: p, divj_p, int_p, cdens_p, edens_p
		integer(I4) :: i,k
		real(DP), dimension(3) :: center
		
		p=.false.; divj_p=.false.; int_p=.false.
		cdens_p=.false.; edens_p=.false.

		call plot_grid_xyz('koord.xyz', grid, mol, 0)
		do k=1,size(calc)
			select case(calc(k))
			case(CDENS_TAG)
				i=0
				call init_grid(grid, mol)
				call grid_center(grid,center)
				call getkw(input, 'grid.gridplot', i)
				if (master_p) then
					if (i > 0) then
						call plot_grid_xyz('grid.xyz', grid, mol, i)
					end if
				end if
!                call proper_coordsys(grid)
			case(DIVJ_TAG)
				i=0
				if (.not.section_is_set(input, 'divj')) then
					call msg_error('Divergence calculation requested, '//& 
					& 'but ''divj'' is undefined in input!')
					cycle
				end if
			
				p=section_is_set(input, 'divj.grid')
				if (p) then
					call push_section(input, 'divj')
					call init_grid(dgrid, mol)
					call grid_center(dgrid,center)
					call getkw(input, 'grid.gridplot', i)
					call pop_section(input)
				else
					call msg_info('Using default grid for divergence calucaltion')
					call init_grid(dgrid, mol)
				end if
				if (master_p) then
					if (i > 0) then
						call plot_grid_xyz('divj.xyz', dgrid, mol, i)
					end if
				end if
!                call proper_coordsys(grid)
			case(INTGRL_TAG)
				i=0
				p=keyword_is_set(input, 'integral.grid')
				if (p) then
					call push_section(input, 'integral')
					call init_grid(igrid, mol)
					call grid_center(igrid,center)
					call getkw(input, 'grid.gridplot', i)
					call pop_section(input)
				else
					call init_grid(igrid, mol)
				end if
				if (master_p) then
					if (i > 0) then
						call plot_grid_xyz('integral.xyz', igrid, mol, i)
					end if
				end if
!                call proper_coordsys(grid)
			case(EDENS_TAG)
				i=0
				p=keyword_is_set(input, 'edens.grid')
				if (p) then
					call push_section(input, 'edens')
					call init_grid(egrid, mol)
					call grid_center(egrid,center)
					call getkw(input, 'grid.gridplot', i)
					call pop_section(input)
				else
					call init_grid(egrid, mol)
				end if
				if (master_p) then
					if (i > 0) then
						call plot_grid_xyz('edens.xyz', egrid, mol, i)
					end if
				end if
			case(0)
				continue
			end select
		end do
		
	end subroutine

	subroutine cdens
		use globals_m
		use basis_m
		use cao2sao_m
		use dens_m
		use jtensor_m
		use jfield_m
		use caos_m
		use grid_m
		use gaussint_m
		use integral_m
		use divj_m
		use edens_m
		use parallel_m

		type(jfield_t) :: jf
		type(grid_t) :: grid, igrid, dgrid, egrid
		type(jtensor_t) :: jt
		type(dens_t) :: xdens, modens
		type(divj_t) :: dj
		type(edens_t) :: ed
		type(integral_t) :: it
		type(parallel_t) :: pt
		type(cao2sao_t) :: c2s

		integer(I4) :: i,j, ncalc
		integer(I4), dimension(4) :: calc
		logical :: divj_p, int_p, cdens_p, edens_p, rerun_p
		logical :: nike_p, xdens_p, modens_p
		character(LINELEN), dimension(:), pointer :: cstr

		divj_p=.false.; int_p=.false.
		cdens_p=.false.; edens_p=.false.
		rerun_p=.false.; nike_p=.true.
		xdens_p=.false.; modens_p=.false.

		if (spherical) then
			call init_c2sop(c2s,mol)
			call set_c2sop(mol, c2s)
		end if
		call getkw(input, 'rerun', rerun_p)
		
		! figure out work order
		nullify(cstr)
		call getkw(input, 'calc', cstr)
		ncalc=size(cstr)
		do i=1,ncalc
!            j=size(cstr(i))+1
!            write(str_g, *) cstr(i)
!            if (str_g(2:j) == 'cdens') then 
			if (cstr(i)(1:5) == 'cdens') then 
				calc(i)=CDENS_TAG
				cdens_p=.true.
				xdens_p=.true.
!            else if (str_g(2:j) == 'integral') then 
			else if (cstr(i)(1:8) == 'integral') then 
				calc(i)=INTGRL_TAG
				int_p=.true.
				xdens_p=.true.
!            else if (str_g(2:j) == 'divj') then 
			else if (cstr(i)(1:4) == 'divj') then 
				calc(i)=DIVJ_TAG
				divj_p=.true.
				xdens_p=.true.
!            else if (str_g(2:j) == 'edens') then 
			else if (cstr(i)(1:5) == 'edens') then 
				calc(i)=EDENS_TAG
				edens_p=.true.
				modens_p=.true.
			end if
		end do
		if (mpirun_p .or. rerun_p) nike_p=.false.

		if (xdens_p) call init_dens(xdens, mol)
		if (modens_p) call init_dens(modens, mol, modens_p)
		if (master_p) then
			if (xdens_p) call read_dens(xdens)
			if (modens_p) call read_modens(modens)
		end if
		if (mpirun_p) then
			if (xdens_p) call bcast_dens(xdens)
			if (modens_p) call bcast_dens(modens)
		end if

		call setup_grids(calc(1:ncalc), grid, igrid, dgrid, egrid)

		call init_jtensor(jt,mol,xdens)

		if (divj_p) call init_divj(dj, dgrid, jt)
		if (cdens_p) call init_jfield(jf, jt, grid)
		if (int_p) call init_integral(it, jt, jf, igrid)
		if (edens_p) call init_edens(ed, mol, modens, egrid)
		if (.not.rerun_p .and. mpirun_p) then
			call init_parallel(pt, calc, jf, it, dj, ed)
		end if

		if (master_p) then
			if (mpirun_p .and. .not.rerun_p) then
				call scheduler(pt)
			end if
			do i=1,ncalc
				select case(calc(i))
				case(CDENS_TAG)
					if (nike_p) call jfield(jf)
					! Contract the tensors with B
					call jvectors(jf)
					call jvector_plot(jf)
				case(INTGRL_TAG)
					if (nike_p) then
!                        call int_t_direct(it)
						call int_s_direct(it)
						call nl
						call msg_info('Integrating |J|')
						call int_mod_direct(it)
					end if
!                    call write_integral(it)
				case(DIVJ_TAG)
					if (nike_p) call divj(dj)
					call divj_plot(dj)
				case(EDENS_TAG)
					if (nike_p) call edens(ed)
					call edens_plot(ed)
				case default
					call msg_error('gimic(): Unknown operation!')
				end select
			end do
		else
			call ceo(pt)
		end if

		if (mpirun_p .and. .not.rerun_p) call del_parallel(pt)
		if (int_p) then
			call del_integral(it)
			call del_grid(igrid)
		end if
		if (divj_p) then
			call del_divj(dj)
			call del_grid(dgrid)
		end if
		if (cdens_p) then
			call del_jfield(jf)
			call del_grid(grid)
		end if
		if (edens_p) then
			call del_edens(ed)
			call del_grid(egrid)
		end if
		if (xdens_p) call del_dens(xdens)
		if (modens_p) call del_dens(modens)
		if (spherical) call del_c2sop(c2s)
	end subroutine

	subroutine program_header
		integer(I4) :: i,j,sz
		integer(I4), dimension(3) :: iti
		integer(I4), dimension(:), allocatable :: seed

		call random_seed(size=sz)
		allocate(seed(sz))
		call random_seed(get=seed)
		call itime(iti)
		j=sum(iti)
		do i=1,sz
			seed(i)=seed(i)*j
		end do
		call random_seed(put=seed)
		deallocate(seed)

call nl
call msg_out('****************************************************************')
call msg_out('***                                                          ***')
call msg_out('***                    GIMIC                                 ***')
call msg_out('***           Coded by Jonas Juselius                        ***')
call msg_out('***                                                          ***')
call msg_out('***  This software is copyright (c) 2003 by Jonas Juselius,  ***')
call msg_out('***  University of Helsinki.                                 ***')
call msg_out('***                                                          ***')
call msg_out('***  You are free to distribute this software under the      ***')
call msg_out('***  terms of the GNU General Public License                 ***')
call msg_out('***                                                          ***')
call msg_out('****************************************************************')
call nl
call msg_out(repeat('=',75))
call msg_out('Said about GIMIC in the press:')
call msg_out('  - A Pretty Advanced ''Hello World!'' Program')
call msg_out(repeat('=',75))
call nl
	end subroutine

	subroutine program_footer
		real(DP) :: rnd
		character(*), dimension(6), parameter :: raboof=(/ &
			'GIMIC - Grossly Irrelevant Magnetically Incuced Currents', &
			'GIMIC - Gone Interrailing, My Inspiration Croaked       ', &
			'GIMIC - Galenskap I Miniatyr, Ingen Censur              ', &
			'GIMIC - Gone Insane, My Indifferent Cosmos              ', &
			'GIMIC - Give Idiots More Ice-Coffee                     ', &
			'GIMIC - Gleaming Indubitably Mediocre Inapplicable Crap '/)

		call random_number(rnd)
		call nl
		call msg_out(raboof(nint(rnd*6.d0)))
		call nl
	end subroutine

	! this routine is deprecated. everything is handled by getkw now
	subroutine cmdline(inpfile)
		use parallel_m
		character(*), intent(out) :: inpfile
		external getarg, iargc, getenv

		integer(I4) :: iargc, fgetc 
		character(BUFLEN) :: fdate
		integer(I4) :: argc
		character(BUFLEN) :: mpienv

	! figure out command line...
		argc=iargc()
		select case(argc)
			case(0) ! open default file
				inpfile=DEFAULT_INPUT
			case(1)
				call getarg(1, inpfile)
				if (inpfile(1:1) == '-') then
					if (trim(inpfile) /= '--mpi') call usage
					mpirun_p=.true.
				else
					call usage
				end if
				inpfile=DEFAULT_INPUT
			case (2) 
				call getarg(1, inpfile)
				if (trim(inpfile) /= '--mpi') call usage
				mpirun_p=.true.
				call getarg(2,inpfile)
			case default 
				call usage
		end select
		call getenv('MPIRUN', mpienv)
		if ( trim(mpienv) == '1' ) mpirun_p=.true.
	end subroutine

	subroutine usage()
		call msg_out('usage: gimic [--mpi] [file]')
		call exit(1)
	end subroutine

end program 

