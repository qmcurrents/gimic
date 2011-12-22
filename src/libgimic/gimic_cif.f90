!
! GIMIC Fortran to C interface
!
! The routinese are to be used from C/C++ to access GIMIC functionality
!
module gimic_cif
    use globals_m  
    use settings_m  
    use teletype_m
    use basis_class
    use timer_m
    use cao2sao_class
    use basis_class
    use dens_class
    use jtensor_class
    use caos_m
    use gaussint_m
    implicit none

    type(molecule_t), save :: mol
    type(dens_t), save :: xdens
    type(cao2sao_t) :: c2s
    real(DP), dimension(3) :: magnet
end module

subroutine gimic_init(molfile, densfile)
    use gimic_cif
    character(*), intent(in) :: molfile
    character(*), intent(in) :: densfile

    call set_debug_level(0)
    settings%title="libgimic"
    settings%basis=molfile
    settings%xdens=densfile
    settings%density=""
    settings%magnet_axis="X"
    settings%magnet=(/0.0,0.0,0.0/)
    settings%is_uhf=.false.
    settings%mofile=""
    settings%morange=0
    settings%use_spherical=.false.
    settings%screen_thrs = SCREEN_THRS
    settings%use_giao = .true.
    settings%use_diamag = .true.
    settings%use_paramag = .true.
    settings%use_screening = .true.
    settings%lip_order = 5
    call new_basis(mol, settings%basis, settings%screen_thrs)

    if (settings%use_screening) then
        call new_basis(mol, settings%basis, settings%screen_thrs)
    else
        call new_basis(mol, settings%basis, -1.d0)
    end if

    if (settings%use_spherical) then
        call new_c2sop(c2s,mol)
        call set_c2sop(mol, c2s)
    end if

    call new_dens(xdens, mol)
    call read_dens(xdens, settings%xdens)
end subroutine

subroutine gimic_finalize()
    use gimic_cif
    if (settings%use_spherical) then
        call del_c2sop(c2s)
    end if
    call del_dens(xdens)
    call del_basis(mol)
end subroutine
