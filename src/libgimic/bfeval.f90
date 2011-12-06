!
! Basis function evaluator
!

module bfeval_class
    use globals_m
    use settings_m
    use basis_class
    use gtodefs_m
    use caos_m
    use cao2sao_class
    implicit none

    public bfeval, new_bfeval, del_bfeval, bfeval_t

    type bfeval_t
        type(molecule_t), pointer :: mol
        real(DP), dimension(:), pointer :: bf
        real(DP), dimension(:), pointer :: sbf
    end type

    private
    
    integer(I4) :: i, j, k, natoms, nctr, ncomp, nccomp, idx
    type(atom_t), pointer :: atom
    type(basis_t), pointer :: basis
    type(contraction_t), pointer :: ctr
    real(DP), dimension(3) :: rr, coord
        
contains

    subroutine new_bfeval(this, mol)
        type(bfeval_t) :: this
        type(molecule_t), target :: mol
        
        allocate(this%bf(get_ncgto(mol)))
        if (main%use_spherical) allocate(this%sbf(get_nccgto(mol)))
        this%mol=>mol
    end subroutine

    subroutine del_bfeval(this)
        type(bfeval_t) :: this

        deallocate(this%bf)
        if (main%use_spherical) deallocate(this%sbf)
        nullify(this%sbf)
        nullify(this%mol)
    end subroutine
    
    subroutine bfeval(this, r, ans)
        type(bfeval_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(:), pointer :: ans

        integer(I4), dimension(99) :: posvec
        integer(I4) :: idx1, idx2
        
        natoms=get_natoms(this%mol)
        
        idx=1
        idx2=0
        posvec=0
        this%bf=0.d0
        if (main%use_spherical) this%sbf=0.d0
        do i=1,natoms
            call get_atom(this%mol,i,atom)
            call get_coord(atom, coord)
            call get_basis(atom, basis)
            rr=r-coord
            call filter_screened(basis, rr, posvec, nctr)
            do k=1,nctr 
                j=posvec(k)
                call get_contraction(atom, j, ctr)
                idx=idx2+get_ctridx(basis, j)
                call cgto(rr, ctr, this%bf(idx:))
            end do
            idx2=idx2+get_ncgto(basis)
        end do
        if (main%use_spherical) then
            call cao2sao(this%mol%c2s, this%bf, this%sbf)
            ans=>this%sbf
        else
            ans=>this%bf
        end if
    end subroutine 

    subroutine bfeval_ns(this, r, ans)
        type(bfeval_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(:), pointer :: ans

        integer(I4), dimension(99) :: posvec
        integer(I4) :: idx1, idx2

        natoms=get_natoms(this%mol)
        
        idx=1
        idx2=0
        this%bf=0.d0
        do i=1,natoms
            call get_atom(this%mol,i,atom)
            call get_coord(atom, coord)
            call get_basis(atom, basis)
            rr=r-coord
            nctr=get_nctr(basis)
            write(89,*) 'atm', i, rr
            do j=1,nctr 
                call get_contraction(atom, j, ctr)
                idx=idx2+get_ctridx(basis, j)
                call cgto(rr, ctr, this%bf(idx:))
                write(89,*) j, idx
                write(89,*) this%bf(idx:idx+ctr%nccomp-1)
                write(89,*)
            end do
            idx2=idx2+get_ncgto(basis)
        end do
        if (main%use_spherical) then
            call cao2sao(this%mol%c2s, this%bf, this%sbf)
            ans=>this%sbf
        else
            ans=>this%bf
        end if
    end subroutine 

end module

! vim:et:sw=4:ts=4
