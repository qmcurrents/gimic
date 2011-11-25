!! @brief  
! GIMIC - a pretty advanced 'Hello World!' program.
!

#include "config.h"

program gimic
	use globals_m  
	use teletype_m
	use basis_class
	use timer_m
	use magnet_m
	use parallel_m
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
	subroutine initialize()
		integer(I4) :: i, hostnm, rank, ierr
		integer(I4) :: chdir, system
		character(BUFLEN) :: title, fdate, sys, molfil, denfil
		character(4) :: rankdir
		real(DP), dimension(3) :: center

		spherical=.false.
		call getkw(input, 'spherical', spherical)
		call getkw(input, 'basis', molfil)
		call getkw(input, 'density', denfil)

		if (mpirun_p) then
			rank=start_mpi()
            call rankname(rankdir)
		else
			master_p=.true.
			rank=0
		end if

		ierr=hostnm(sys)
		if (master_p) then
			write(str_g, '(a,i3,a,a)') 'MPI master', rank,' on ', trim(sys)
			call msg_debug(str_g,2)
		else
            call set_teletype_unit(DEVNULL)
			write(str_g, '(a,i3,a,a)') 'MPI slave', rank,' on ', trim(sys)
			call msg_debug(str_g,2)
		end if

		call getkw(input, 'debug', debug_level)
		call set_debug_level(debug_level)

		call msg_out(fdate())
		
		call getkw(input, 'title', title)
		call msg_out(' TITLE: '// trim(title))
		call nl
		
        call init_basis(mol,molfil)

	end subroutine

	subroutine finalize()
		call del_basis(mol)
		call stop_mpi()
		call del_getkw(input)
	end subroutine

	subroutine cdens
		use globals_m
		use basis_class
		use cao2sao_class
		use dens_class
		use jtensor_class
		use jfield_class
		use caos_m
		use grid_class
		use gaussint_m
		use integral_class
		use divj_class
		use edens_class

		type(jfield_t) :: jf
		type(grid_t) :: cgrid, igrid, dgrid, egrid
		type(jtensor_t) :: jt
		type(dens_t) :: xdens, modens
		type(divj_t) :: dj
		type(edens_t) :: ed
		type(integral_t) :: it
		type(cao2sao_t) :: c2s

		integer(I4) :: i,j, ncalc
		integer(I4), dimension(4) :: calc
		logical :: divj_p, int_p, cdens_p, edens_p
		logical :: xdens_p, modens_p, dryrun_p, imod_p
		character(LINELEN), dimension(:), pointer :: cstr


		divj_p=.false.; int_p=.false.
		cdens_p=.false.; edens_p=.false.; dryrun_p=.false.
		xdens_p=.false.; modens_p=.false.; imod_p=.false.

		call getkw(input,'dryrun', dryrun_p)

		if (spherical) then
			call init_c2sop(c2s,mol)
			call set_c2sop(mol, c2s)
		end if

		call getkw(input, 'openshell', uhf_p)
		if (uhf_p) then
			call msg_info('Open-shell calculation')
		else
			call msg_info('Closed-shell calculation')
		end if
		call nl

		
		! figure out work order
		nullify(cstr)
		call getkw(input, 'calc', cstr)
		ncalc=size(cstr)
		do i=1,ncalc
			if (cstr(i)(1:5) == 'cdens') then 
				calc(i)=CDENS_TAG
				cdens_p=.true.
				xdens_p=.true.
			else if (cstr(i)(1:8) == 'integral') then 
				calc(i)=INTGRL_TAG
				int_p=.true.
				xdens_p=.true.
			else if (cstr(i)(1:4) == 'divj') then 
				calc(i)=DIVJ_TAG
				divj_p=.true.
				xdens_p=.true.
			else if (cstr(i)(1:5) == 'edens') then 
				calc(i)=EDENS_TAG
				edens_p=.true.
				modens_p=.true.
			end if
		end do

		if (xdens_p) then
			call init_dens(xdens, mol)
			call read_dens(xdens)
		end if
		if (modens_p) then
			call init_dens(modens, mol, modens_p)
			call read_modens(modens)
		end if


		call init_jtensor(jt,mol,xdens)

		if (dryrun_p) then
			call msg_note('Dry run, not calculating...')
			call nl
		end if

		do i=1,ncalc
			select case(calc(i))
			case(CDENS_TAG)
				call msg_out('Calculating current density')
				call msg_out('*****************************************')
				call setup_grid(calc(i), cgrid)
				call init_jfield(jf, jt, cgrid)
				if (dryrun_p) cycle
				call jfield(jf)
					! Contract the tensors with B
				if (master_p) then
					call jvectors(jf)
					call jvector_plot(jf)
				end if
			case(INTGRL_TAG)
				call msg_out('Integrating current density')
				call msg_out('*****************************************')
				call setup_grid(calc(i), igrid)
				call init_integral(it, jt, jf, igrid)
				if (dryrun_p) cycle
				call int_s_direct(it)
				call nl
				call getkw(input, 'integral.modulus', imod_p)
				if (imod_p) then
					call msg_note('Integrating |J|')
					call int_mod_direct(it)
				end if
				call getkw(input, 'integral.tensor', imod_p)
				if (imod_p) then
					call msg_note('Integrating current tensor')
					call int_t_direct(it)  ! tensor integral
				end if
!                    call write_integral(it)
			case(DIVJ_TAG)
				call msg_out('Calculating divergence')
				call msg_out('*****************************************')
				call setup_grid(calc(i), dgrid)
				call init_divj(dj, dgrid, jt)
				if (dryrun_p) cycle
				call divj(dj)
				if (master_p) call divj_plot(dj)
			case(EDENS_TAG)
				call msg_out('Calculating charge density')
				call msg_out('*****************************************')
				call setup_grid(calc(i), egrid)
				call init_edens(ed, mol, modens, egrid)
				if (dryrun_p) cycle
				call edens(ed)
				if (master_p) call edens_plot(ed)
			case default
				call msg_error('gimic(): Unknown operation!')
			end select
		end do

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
			call del_grid(cgrid)
		end if
		if (edens_p) then
			call del_edens(ed)
			call del_grid(egrid)
		end if
		if (xdens_p) call del_dens(xdens)
		if (modens_p) call del_dens(modens)
		if (spherical) call del_c2sop(c2s)
		call del_jtensor(jt)
	end subroutine

	subroutine setup_grid(calc, grid)
		use grid_class
		integer(I4) :: calc
		type(grid_t) :: grid
		logical :: p, divj_p, int_p, cdens_p, edens_p
		integer(I4) :: i
		real(DP), dimension(3) :: center
		
		p=.false.; divj_p=.false.; int_p=.false.
		cdens_p=.false.; edens_p=.false.

!        call plot_grid_xyz(grid, 'koord.xyz', mol)
		select case(calc)
			case(CDENS_TAG)
				call push_section(input, 'cdens')
				call init_grid(grid, mol)
				call grid_center(grid,center)
				call pop_section(input)
				if (master_p) then
					call plot_grid_xyz(grid, 'grid.xyz',  mol)
				end if
			case(DIVJ_TAG)
				if (.not.section_is_set(input, 'divj')) then
					call msg_error('Divergence calculation requested, '//& 
					& 'but ''divj'' is undefined in input!')
				else
					p=section_is_set(input, 'divj.grid')
					if (p) then
						call push_section(input, 'divj')
						call init_grid(grid, mol)
						call grid_center(grid,center)
						call pop_section(input)
					else
						call msg_critical('No grid definition for section divj')
					end if
					if (master_p) then
						call plot_grid_xyz(grid, 'divj.xyz', mol)
					end if
				end if
			case(INTGRL_TAG)
				p=keyword_is_set(input, 'integral.grid')
				if (p) then
					call push_section(input, 'integral')
					call init_grid(grid, mol)
					call grid_center(grid,center)
					call pop_section(input)
				else
					call msg_critical('No grid definition for section integral')
				end if
				if (master_p) then
					call plot_grid_xyz(grid, 'integral.xyz', mol)
				end if
			case(EDENS_TAG)
				i=0
				p=keyword_is_set(input, 'edens.grid')
				if (p) then
					call push_section(input, 'edens')
					call init_grid(grid, mol)
					call grid_center(grid,center)
					call pop_section(input)
				else
					call msg_critical('No grid definition for section edens')
				end if
				if (master_p) then
					call plot_grid_xyz(grid, 'edens.xyz', mol)
				end if
			case(0)
				continue
		end select
		
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
call msg_out('***           GIMIC '// GIMIC_VERSION // &
                                      '                                    ***')
call msg_out('***              Written by Jonas Juselius                   ***')
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
		character(*), dimension(5), parameter :: raboof=(/ &
			'GIMIC - Grossly Irrelevant Magnetically Incuced Currents', &
			'GIMIC - Gone Interrailing, My Inspiration Croaked       ', &
			'GIMIC - Galenskap I Miniatyr, Ingen Censur              ', &
			'GIMIC - Gone Insane, My Indifferent Cosmos              ', &
			'GIMIC - Give Idiots More Ice-Coffee                     '/)

		call random_number(rnd)
		call nl
		call msg_out(raboof(nint(rnd*5.d0)))
		call nl
	end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Deprecated routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! this routine is deprecated. everything is handled by getkw now
	subroutine cmdline(inpfile)
		character(*), intent(out) :: inpfile

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

!
! read in the input file into a buffer which can be broadcast 
!

	subroutine read_inpbuf(inpfile, inpbuf)
		character(*), intent(in) :: inpfile
		character, dimension(:), pointer :: inpbuf

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

		n=getfsize(inpfile)
		allocate(inpbuf(n))
		do i=1,n
			err=fgetc(99, inpbuf(i))
		end do
		close(99)
100 format(a,a,a)
	end subroutine

!
! read in the MOL file into a buffer which can be broadcast 
!
	subroutine read_molbuf(molbuf)
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

end program 

