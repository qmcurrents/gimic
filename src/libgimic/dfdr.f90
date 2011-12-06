
module dfdr_class
    use globals_m
    use settings_m
    use basis_class
    use caos_m
    use gtodefs_m
    use cao2sao_class
    implicit none

    public new_dfdr, del_dfdr, dfdr, dfdr_t

    type dfdr_t
        type(molecule_t), pointer :: mol
        real(DP), dimension(:,:), pointer :: dr
        real(DP), dimension(:,:), pointer :: sdr
    end type

    private
    
    integer(I4) :: i, j, k, natoms, nctr, nccomp, axis
    type(atom_t), pointer :: atom
    type(basis_t), pointer :: basis
    type(contraction_t), pointer :: ctr
    real(DP), dimension(3) :: rr, coord
    integer(I4) :: idx
        
contains

    subroutine new_dfdr(this, mol)
        type(dfdr_t) :: this
        type(molecule_t), target :: mol

        allocate(this%dr(get_nccgto(mol),3))
        if (use_spherical) allocate(this%sdr(get_ncgto(mol),3))
        this%mol=>mol
!        call dfdr(this, (/0.d0, 0.d0, 0.d0/), this%dr)
    end subroutine

    subroutine del_dfdr(this)
        type(dfdr_t) :: this

        deallocate(this%dr)
        if (use_spherical) deallocate(this%sdr)
    end subroutine

    subroutine dfdr(this, r, drv)
        type(dfdr_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(:,:), pointer :: drv

        integer(I4), dimension(99) :: posvec
        integer(I4) :: idx1, idx2
        
        natoms=get_natoms(this%mol)
        
        idx2=0
        this%dr=0.d0
        do i=1,natoms
            call get_atom(this%mol,i,atom)
            call get_coord(atom, coord)
            rr=r-coord
            call get_basis(atom, basis)
            call filter_screened(basis, rr, posvec, nctr)
            do k=1,nctr 
                j=posvec(k)
                call get_contraction(atom, j, ctr)
                idx=idx2+get_ctridx(basis, j)
                do axis=1,3
                    call dcgto(rr, ctr, axis, this%dr(idx:, axis))
                end do
            end do
            idx2=idx2+get_ncgto(basis)
        end do
        if (use_spherical) then
            do axis=1,3
                call cao2sao(this%mol%c2s, this%dr(:, axis), this%sdr(:,axis))
            end do
            drv=>this%sdr
        else
            drv=>this%dr
        end if
    end subroutine 
end module

! vim:et:sw=4:ts=4
