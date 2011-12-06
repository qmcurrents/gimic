
module dbop_class
    use globals_m
    use settings_m
    use basis_class
    use teletype_m
    implicit none
    
    public new_dbop, del_dbop, mkdbop, dbop_t

    type dbop_t 
        type(molecule_t), pointer :: mol
        real(DP), dimension(:,:), pointer :: db
    end type
    private

contains
    subroutine new_dbop(this, mol)
        type(dbop_t) :: this
        type(molecule_t), target :: mol

        nullify(this%db)

        if (associated(this%db)) then
            call msg_warn('new_dbop(): already allocated!')
        else
            allocate(this%db(3,get_natoms(mol)))
        end if
        this%mol=>mol
    end subroutine

    subroutine del_dbop(this)
        type(dbop_t) :: this
        if (associated(this%db)) then
            deallocate(this%db)
            nullify(this%db)
            nullify(this%mol)
        else
            call msg_warn('del_dbop(): not allocated!')
        end if
    end subroutine
    
    subroutine mkdbop(this, r, dbop)
        type(dbop_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(:,:), pointer :: dbop
        
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
            this%db(1,i)=(r(2)*Rz-r(3)*Ry)
            this%db(2,i)=(r(3)*Rx-r(1)*Rz)
            this%db(3,i)=(r(1)*Ry-r(2)*Rx)
        end do
        dbop=>this%db
    end subroutine

end module

! vim:et:sw=4:ts=4
