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
    type(jtensor_t), save :: jt
    type(dens_t), save :: xdens
    real(DP), dimension(3) :: magnet

end module

subroutine gimic_init(archive)
    use gimic_cif
    character(*), intent(in) :: archive
    call new_gimic(run, archive)
end subroutine

subroutine gimic_finalize()
    use gimic_cif
    call del_gimic(run)
end subroutine

