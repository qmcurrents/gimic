!
! Calulate the electronic density in a point.
! Coded by Jonas Juselius <jonas@iki.fi> 2003 & 2011
! 

module edens_class
    use globals_m
    use settings_m
    use basis_class
    use bfeval_class
    use dens_class
    use teletype_m
    implicit none
    
    public new_edens, del_edens, edens_t
    public edens
    
    type edens_t
        real(DP), dimension(:), pointer :: tmp
        real(DP), dimension(:,:), pointer :: aodens
        type(bfeval_t) :: bf
    end type

    private

contains
    ! set up memory (once) for the different components
    subroutine new_edens(this, mol, dens)
        type(edens_t) :: this
        type(molecule_t) :: mol
        type(dens_t), target :: dens

        integer(I4) :: n

        call get_dens(dens, this%aodens)
        n=get_ncgto(mol)
        allocate(this%tmp(n))
        call new_bfeval(this%bf, mol)
    end subroutine

    subroutine del_edens(this)
        type(edens_t) :: this

        if (associated(this%tmp)) deallocate(this%tmp)
        call del_bfeval(this%bf)
    end subroutine

    function edens(this, rr) result(ed)
        type(edens_t) :: this
        real(DP), dimension(3), intent(in) :: rr
        real(DP) :: ed

        real(DP), dimension(:,:), pointer :: buf
        real(DP), dimension(:), pointer :: bfvec

        call bfeval(this%bf,rr, bfvec)
        this%tmp=matmul(this%aodens, bfvec)
        ed = dot_product(this%tmp, bfvec)
    end function
end module

! vim:et:sw=4:ts=4
