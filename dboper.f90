!
! $Id$
!

module dboper_m
	use globals_m
	use basis_m
	use teletype_m
	implicit none
	
	public init_dbop, del_dbop, mkdbop, dbop_t

	type dbop_t 
		type(molecule_t), pointer :: mol
		real(DP), dimension(:,:), pointer :: db
		real(DP), dimension(3) :: r
	end type
	private

contains
	subroutine init_dbop(dop, mol)
		type(dbop_t) :: dop
		type(molecule_t), target :: mol
		if (associated(dop%db)) then
			call msg_warn('init_dbop(): already allocated!')
		else
			allocate(dop%db(3,get_natoms(mol)))
		end if
		dop%mol=>mol
		dop%r=INITRV
	end subroutine

	subroutine del_dbop(dop)
		type(dbop_t) :: dop
		if (associated(dop%db)) then
			deallocate(dop%db)
			nullify(dop%db)
			nullify(dop%mol)
		else
			call msg_warn('del_dbop(): not allocated!')
		end if
	end subroutine
	
	subroutine mkdbop(dop, r, dbop)
		type(dbop_t) :: dop
		real(DP), dimension(3), intent(in) :: r
		real(DP), dimension(:,:), pointer :: dbop
		
		type(atom_t), pointer :: atom
		real(DP) :: Rx, Ry, Rz
		real(DP), dimension(3) :: coord
		integer(I4) :: i, natoms
		
		if (r(1)==dop%r(1) .and. r(2)==dop%r(2) .and. r(3)==dop%r(3)) then
			dbop=>dop%db
			return
		end if
		
		dop%r=r
		natoms=get_natoms(dop%mol)

		do i=1,natoms
			call get_atom(dop%mol, i, atom)
			call get_coord(atom, coord)
			Rx=coord(1)
			Ry=coord(2)
			Rz=coord(3)
			dop%db(1,i)=(r(2)*Rz-r(3)*Ry)
			dop%db(2,i)=(r(3)*Rx-r(1)*Rz)
			dop%db(3,i)=(r(1)*Ry-r(2)*Rx)
		end do
		dbop=>dop%db
	end subroutine

end module
