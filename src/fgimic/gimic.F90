!
! GIMIC - a pretty advanced 'Hello World!' program.
!

#include "config.h"

program gimic
    use globals_module
    use getkw_class
    use settings_module
    use teletype_module
    use basis_class
    use timer_module
    use magnet_module
    use parallel_module
    use cao2sao_class
    use grid_class
    use basis_class
    use dens_class
    use jtensor_class
    use jfield_class
    use caos_module
    use gaussint_module
    ! acid and jav stuff
    use acid_module
    implicit none

    type(molecule_t) :: mol
    type(getkw_t) :: input
    type(grid_t) :: grid
    type(jfield_t) :: jf
    type(jtensor_t) :: jt
    type(dens_t) :: xdens
    real(DP), dimension(3) :: magnet

    character(BUFLEN) :: buf

    call new_getkw(input)
    call set_debug_level(3)

    mpi_rank = init_mpi(settings%is_mpirun)

    call program_header
    call initialize()
    call driver()
    call finalize()

    call stockas_klocka()
    call msg_out('Hello World! (tm)')
    call program_footer()

    if (your_results_are_questionable) then
        call msg_error('Your results are questionable.')
        call nl
    end if

contains
    subroutine initialize()
        integer(I4) :: i, hostnm, ierr
        integer(I4) :: chdir, system
        character(BUFLEN) :: title, fdate, sys
        real(DP), dimension(3) :: magnet

        logical :: screen

        call getkw(input, 'debug', debug_level)
        call set_debug_level(debug_level)

        call getkw(input, 'title', settings%title)
        call getkw(input, 'basis', settings%basis)
        call getkw(input, 'xdens', settings%xdens)
        call getkw(input, 'density', settings%density)
        call getkw(input, 'magnet_axis', settings%magnet_axis)
        call getkw(input, 'magnet', settings%magnet)
        call getkw(input, 'openshell', settings%is_uhf)
        call getkw(input, 'dryrun', settings%dryrun)
        call getkw(input, 'calc', settings%calc)
        call getkw(input, 'xdens', settings%xdens)
        call getkw(input, 'density', settings%density)
        call getkw(input, 'mofile', settings%mofile)
        call getkw(input, 'mos', settings%morange)

        settings%use_spherical=.false.
        settings%screen_thrs = SCREEN_THRS
        call getkw(input, 'Advanced.spherical', settings%use_spherical)
        call getkw(input, 'Advanced.GIAO', settings%use_giao)
        call getkw(input, 'Advanced.diamag', settings%use_diamag)
        call getkw(input, 'Advanced.paramag',settings%use_paramag)
        call getkw(input, 'Advanced.screening', settings%use_screening)
        call getkw(input, 'Advanced.screening_thrs', settings%screen_thrs)

        call getkw(input, 'Essential.acid', settings%acid)
        call getkw(input, 'Essential.jmod', settings%jmod)
        call getkw(input, 'Essential.prop', settings%prop)

        ierr=hostnm(sys)
        if (mpi_rank == 0) then
            write(str_g, '(a,i3,a,a)') 'MPI master', mpi_rank,' on ', trim(sys)
            call msg_debug(str_g,2)
        else
            call set_teletype_unit(DEVNULL)
            write(str_g, '(a,i3,a,a)') 'MPI slave', mpi_rank,' on ', trim(sys)
            call msg_debug(str_g,2)
        end if

        call msg_out(fdate())

        call msg_out(' TITLE: '// trim(settings%title))
        call nl

        if (.not.settings%use_giao) then
            call msg_info('GIAOs not used!')
            call nl
        end if

        if (.not.settings%use_diamag) then
            call msg_info( 'Diamagnetic contributions not calculated!')
            call nl
        end if
        if (.not.settings%use_paramag) then
            call msg_info( 'Paramagnetic contributions not calculated!')
            call nl
        end if
        if ((.not.settings%use_diamag).and.(.not.settings%use_paramag)) then
            call msg_out( '    ...this does not make sense...')
            call nl
            call msg_critical( '    PLEASE SEEK PROFESSIONAL HELP, ASAP!  ')
            call nl
            stop
        end if
    end subroutine

    subroutine finalize()
        call del_basis(mol)
        call stop_mpi()
        call del_getkw(input)
    end subroutine

    subroutine driver
        type(cao2sao_t) :: c2s

        if (settings%use_screening) then
            call new_basis(mol, settings%basis, settings%screen_thrs)
        else
            call new_basis(mol, settings%basis, -1.d0)
        end if

        if (.not. settings%dryrun) then
            if (settings%use_spherical) then
                call new_c2sop(c2s,mol)
                call set_c2sop(mol, c2s)
            end if

            call new_dens(xdens, mol)
            call read_dens(xdens, settings%xdens)
        end if

        call new_grid(grid, input, mol)
        if (mpi_rank == 0) then
            call plot_grid_xyz(grid, 'grid.xyz', mol)
        end if

        call get_magnet(grid, magnet)

        if (settings%is_uhf) then
            call msg_info('Open-shell calculation')
        else
            call msg_info('Closed-shell calculation')
        end if
        call nl

        if (settings%dryrun) then
            call msg_note('Dry run, not calculating ...')
            call nl
        end if

        if (settings%calc(1:5) == 'cdens') then
            call run_cdens(jf,mol,xdens)
        else if (settings%calc(1:8) == 'integral') then
            call run_integral()
        else
            call msg_error('gimic(): Unknown operation!')
        end if

        if (.not. settings%dryrun) then
            if (settings%use_spherical) then
                call del_c2sop(c2s)
            end if
            call del_dens(xdens)
        end if
        call del_grid(grid)
    end subroutine

    subroutine run_cdens(jf,mol,xdens)
        type(jfield_t) :: jf
        type(dens_t) :: xdens
        type(molecule_t) :: mol
        call msg_out('Calculating current density')
        call msg_out('*****************************************')
        call new_jfield(jf, grid, magnet)

        if (settings%dryrun) return

        call calc_jvectors(jf, mol, xdens)
        call jvector_plots(jf)

        if (settings%is_uhf) then
            call calc_jvectors(jf, mol, xdens, 'alpha')
            call jvector_plots(jf, 'alpha')

            call calc_jvectors(jf, mol, xdens, 'beta')
            call jvector_plots(jf, 'beta')

            call calc_jvectors(jf, mol, xdens, 'spindens')
            call jvector_plots(jf, 'spindens')
        endif
        call del_jfield(jf)
    end subroutine

    subroutine run_integral()
        use integral_class
        type(integral_t) :: it

        call msg_out('Integrating current density')
        call msg_out('*****************************************')
        call new_integral(it, grid)

        if (settings%dryrun) return

        if (settings%jmod) then
          call msg_note('Integrating |J|')
          call integrate_modulus(it, mol, xdens)
          if (settings%is_uhf) then
              call integrate_modulus(it, mol, xdens, 'alpha')
              call integrate_modulus(it, mol, xdens, 'beta')
              call integrate_modulus(it, mol, xdens, 'spindens')
          end if
          call nl
        else
          write(*,*) "Jmod integration skipped."
        end if

        call msg_note('Integrating current')
        call integrate_current(it, mol, xdens)
        if (settings%is_uhf) then
            call integrate_current(it, mol, xdens, 'alpha')
            call integrate_current(it, mol, xdens, 'beta')
            call integrate_current(it, mol, xdens, 'spindens')
        end if
        call nl
        if (settings%acid) then
        ! add here ACID integration !
           call msg_note('Integrating ACID density')
           call integrate_acid(it, mol, xdens)
           call nl
        end if

        call del_integral(it)
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
call msg_out('***           GIMIC '// PROJECT_VERSION // &
                       ' (' // GIT_HASH // ')                              ***')
call msg_out('***              Written by Jonas Juselius                   ***')
call msg_out('***                                                          ***')
call msg_out('***  This software is copyright (c) 2003-2011 by             ***')
call msg_out('***  Jonas Juselius,  University of Tromso.                  ***')
call msg_out('***                                                          ***')
call msg_out('***  You are free to distribute this software under the      ***')
call msg_out('***  terms of the GNU General Public License                 ***')
call msg_out('***                                                          ***')
call msg_out('***  A Pretty Advanced ''Hello World!'' Program                ***')
call msg_out('****************************************************************')
call nl
    end subroutine

    subroutine program_footer
        real(DP) :: rnd
        character(*), dimension(5), parameter :: raboof=(/ &
            'GIMIC - Grossly Irrelevant Magnetically Induced Currents', &
            'GIMIC - Gone Interrailing, My Inspiration Croaked       ', &
            'GIMIC - Galenskap I Miniatyr, Ingen Censur              ', &
            'GIMIC - Gone Insane, My Indifferent Cosmos              ', &
            'GIMIC - Give Idiots More Ice-Coffee                     '/)

        call random_number(rnd)
        call nl
        !call msg_out(raboof(nint(rnd*5.d0)))
        call msg_out('done.')
        call nl
    end subroutine

end program

! vim:et:sw=4:ts=4
