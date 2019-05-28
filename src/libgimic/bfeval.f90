!
! Basis function evaluator
!

module bfeval_class
    use globals_module
    use settings_module
    use basis_class
    use gtodefs_module
    use caos_module
    use cao2sao_class
    implicit none

    public new_bfeval, del_bfeval, bfeval_t
    public bfeval, dfdr, calc_basis

    type bfeval_t
        private
        type(molecule_t), pointer :: mol
        real(DP), dimension(:), pointer :: bf
        real(DP), dimension(:), pointer :: sbf
        real(DP), dimension(:,:), pointer :: dbop
        real(DP), dimension(:,:), pointer :: d2
        real(DP), dimension(:,:), pointer :: db
        real(DP), dimension(:,:), pointer :: dr
        real(DP), dimension(:,:), pointer :: sdr
        integer(I4) :: ptn
    end type

    private
!    real(DP), dimension(:,:), allocatable, target, save :: xvec
contains

    subroutine new_bfeval(this, mol)
        type(bfeval_t) :: this
        type(molecule_t), target :: mol

        allocate(this%bf(get_ncgto(mol)))
        allocate(this%dbop(3,get_natoms(mol)))
        allocate(this%d2(get_ncgto(mol),9))
        allocate(this%db(get_ncgto(mol),3))
        allocate(this%dr(get_nccgto(mol),3))
        if (settings%use_spherical) allocate(this%sdr(get_ncgto(mol),3))
        if (settings%use_spherical) allocate(this%sbf(get_nccgto(mol)))
        this%mol=>mol

    end subroutine

    subroutine del_bfeval(this)
        type(bfeval_t) :: this

        deallocate(this%bf)
        deallocate(this%db)
        deallocate(this%dbop)
        deallocate(this%d2)
        deallocate(this%dr)
        if (settings%use_spherical) deallocate(this%sdr)
        if (settings%use_spherical) deallocate(this%sbf)
    end subroutine

    subroutine calc_basis(this, r, bf, dr, db, d2)
        type(bfeval_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(:), pointer :: bf
        real(DP), dimension(:,:), pointer, optional :: dr
        real(DP), dimension(:,:), pointer, optional :: db, d2

        call bfeval(this, r, bf)
        if (present(dr)) then
            call dfdr(this, r, dr)
        endif
        if (present(db)) then
            call mkdbop(this, r)
            call dfdb(this, bf, db)
        endif
        if (present(d2)) then
            call d2fdrdb(this, bf, dr, d2)
        endif
    end subroutine

    subroutine bfeval(this, r, ans)
        type(bfeval_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(:), pointer :: ans

        integer(I4), dimension(99) :: posvec
        integer(I4) :: idx1, idx2

        integer(I4) :: i, j, k, natoms, nctr, ncomp, nccomp, idx
        type(atom_t), pointer :: atom
        type(basis_t), pointer :: basis
        type(contraction_t), pointer :: ctr
        real(DP), dimension(3) :: rr, coord

        natoms=get_natoms(this%mol)

        idx=1
        idx2=0
        posvec=0
        this%bf=0.d0
        if (settings%use_spherical) this%sbf=0.d0
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
        if (settings%use_spherical) then
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

        integer(I4) :: i, j, k, natoms, nctr, ncomp, nccomp, idx
        type(atom_t), pointer :: atom
        type(basis_t), pointer :: basis
        type(contraction_t), pointer :: ctr
        real(DP), dimension(3) :: rr, coord

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
        if (settings%use_spherical) then
            call cao2sao(this%mol%c2s, this%bf, this%sbf)
            ans=>this%sbf
        else
            ans=>this%bf
        end if
    end subroutine

    subroutine mkdbop(this, r)
        type(bfeval_t) :: this
        real(DP), dimension(3), intent(in) :: r

        type(atom_t), pointer :: atom
        real(DP) :: Rx, Ry, Rz
        real(DP), dimension(3) :: coord
        integer(I4) :: i, natoms

        natoms=get_natoms(this%mol)

        do i=1,natoms
            call get_atom(this%mol, i, atom)
            call get_coord(atom, coord)
            Rx=coord(1)
            Ry=coord(2)
            Rz=coord(3)
            this%dbop(1,i)=(r(2)*Rz-r(3)*Ry)
            this%dbop(2,i)=(r(3)*Rx-r(1)*Rz)
            this%dbop(3,i)=(r(1)*Ry-r(2)*Rx)
        end do
    end subroutine

    subroutine d2fdrdb(this, bfvec, drvec, thisv)
        type(bfeval_t) :: this
        real(DP), dimension(:), intent(in) :: bfvec
        real(DP), dimension(:,:), intent(in) :: drvec
        real(DP), dimension(:,:), pointer :: thisv

        type(atom_t), dimension(:), pointer :: atoms
        integer(I4) :: i, j, k, natoms, nctr, ncomp, idx
        real(DP), dimension(3) :: coord
        type(atom_t), pointer :: atom
        type(basis_t), pointer :: basis
        type(contraction_t), pointer :: ctr

        integer(I4), dimension(99) :: posvec
        integer(I4) :: l, idx1
        real(DP), dimension(3) :: dbov
        real(DP) :: ror1, ror2, ror3

        natoms=get_natoms(this%mol)

        idx=1
        this%d2=0.d0
        do i=1,natoms
            call get_atom(this%mol,i,atom)
            dbov=this%dbop(:,i)
            call get_coord(atom, coord)
            ror1=coord(1)
            ror2=coord(2)
            ror3=coord(3)
            call get_basis(atom, basis)
            nctr=get_nctr(basis)
            do j=1,nctr
                call get_contraction(atom, j, ctr)
                ncomp=get_ncomp(ctr)
                do k=1,ncomp

                    ! dBx
                    this%d2(idx,1)=drvec(idx,1)*dbov(1)                  ! dx
                    this%d2(idx,2)=drvec(idx,2)*dbov(1)+ror3*bfvec(idx)  ! dy
                    this%d2(idx,3)=drvec(idx,3)*dbov(1)-ror2*bfvec(idx)  ! dz

                    ! dBy
                    this%d2(idx,4)=drvec(idx,1)*dbov(2)-ror3*bfvec(idx)
                    this%d2(idx,5)=drvec(idx,2)*dbov(2)
                    this%d2(idx,6)=drvec(idx,3)*dbov(2)+ror1*bfvec(idx)

                    ! dBz
                    this%d2(idx,7)=drvec(idx,1)*dbov(3)+ror2*bfvec(idx)
                    this%d2(idx,8)=drvec(idx,2)*dbov(3)-ror1*bfvec(idx)
                    this%d2(idx,9)=drvec(idx,3)*dbov(3)
                    idx=idx+1
                end do
            end do
        end do
        thisv=>this%d2
    end subroutine

    subroutine print_thisvec(foo)
        real(DP), dimension(:,:), intent(in) :: foo

        integer(I4) :: i
        character(3) :: funk

        select case(size(foo,2))
            case(1)
                funk='s ='
            case(3)
                funk='p ='
            case(5:6)
                funk='d ='
            case(7:10)
                funk='f ='
        end select

        do i=1,size(foo,2)
            print '(a, 9f15.11)', funk, foo(:,i)
            !print '(a, 9e)', funk, foo(:,i)
        end do
    end subroutine

    subroutine dfdb(this, bfvec, dbv)
        type(bfeval_t) :: this
        real(DP), dimension(:), intent(in) :: bfvec
        real(DP), dimension(:,:), pointer :: dbv

        type(atom_t), pointer :: atom
        type(basis_t), pointer :: basis
        integer(I4) :: i, j, k, natoms, ncgto

        natoms=get_natoms(this%mol)

        j=1
        do k=1,natoms
            call get_atom(this%mol, k, atom)
            call get_basis(atom, basis)
            ncgto=get_ncgto(basis)
            do i=1,ncgto
                this%db(j,:)=this%dbop(:,k)*bfvec(j)
                j=j+1
            end do
        end do
        dbv=>this%db
    end subroutine

    subroutine dfdr(this, r, drv)
        type(bfeval_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(:,:), pointer :: drv

        integer(I4), dimension(99) :: posvec
        integer(I4) :: idx1, idx2

        integer(I4) :: i, j, k, natoms, nctr, nccomp, axis
        type(atom_t), pointer :: atom
        type(basis_t), pointer :: basis
        type(contraction_t), pointer :: ctr
        real(DP), dimension(3) :: rr, coord
        integer(I4) :: idx

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
        if (settings%use_spherical) then
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
