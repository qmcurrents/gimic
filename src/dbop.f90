
module dbop_class
	use globals_m
	use basis_class
	use teletype_m
	implicit none
	
	public init_dbop, del_dbop, mkdbop, dbop_t

	type dbop_t 
		type(molecule_t), pointer :: mol
		real(DP), dimension(:,:), pointer :: db
	end type
	private

contains
	subroutine init_dbop(self, mol)
		type(dbop_t) :: self
		type(molecule_t), target :: mol

		nullify(self%db)

		if (associated(self%db)) then
			call msg_warn('init_dbop(): already allocated!')
		else
			allocate(self%db(3,get_natoms(mol)))
		end if
		self%mol=>mol
	end subroutine

	subroutine del_dbop(self)
		type(dbop_t) :: self
		if (associated(self%db)) then
			deallocate(self%db)
			nullify(self%db)
			nullify(self%mol)
		else
			call msg_warn('del_dbop(): not allocated!')
		end if
	end subroutine
	
	subroutine mkdbop(self, r, dbop)
		type(dbop_t) :: self
		real(DP), dimension(3), intent(in) :: r
		real(DP), dimension(:,:), pointer :: dbop
		
		type(atom_t), pointer :: atom
		real(DP) :: Rx, Ry, Rz
		real(DP), dimension(3) :: coord
		integer(I4) :: i, natoms
		
		natoms=get_natoms(self%mol)

		do i=1,natoms
			call get_atom(self%mol, i, atom)
			call get_coord(atom, coord)
			Rx=coord(1)
			Ry=coord(2)
			Rz=coord(3)
			self%db(1,i)=(r(2)*Rz-r(3)*Ry)
			self%db(2,i)=(r(3)*Rx-r(1)*Rz)
			self%db(3,i)=(r(1)*Ry-r(2)*Rx)
		end do
		dbop=>self%db
	end subroutine

end module
