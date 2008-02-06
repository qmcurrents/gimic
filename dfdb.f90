
module dfdb_class
	use globals_m
	use basis_class
	use dbop_class
	use bfeval_class
	use teletype_m
	implicit none
	
	public init_dfdb, del_dfdb, dfdb
	public dfdb_t
	type dfdb_t
		type(molecule_t), pointer :: mol
		type(dbop_t), pointer :: dop
		type(bfeval_t), pointer :: bfv
		real(DP), dimension(:,:), pointer :: db
		real(DP), dimension(3) :: r
	end type
	
	private
	real(DP), dimension(:,:), allocatable, target, save :: xvec

contains
	subroutine init_dfdb(self, mol, dop, bfv)
		type(dfdb_t) :: self
		type(molecule_t), target :: mol
		type(dbop_t), target :: dop
		type(bfeval_t), target :: bfv

		nullify(self%db)

		if (associated(self%db)) then
			call msg_warn('init_dfdb(): already allocated!')
		else
			allocate(self%db(get_ncgto(mol),3))
		end if
		self%mol=>mol
		self%dop=>dop
		self%bfv=>bfv
		self%r=INITRV
	end subroutine

	subroutine del_dfdb(self)
		type(dfdb_t) :: self
		if (associated(self%db)) then
			deallocate(self%db)
			nullify(self%db)
		else
			call msg_warn('del_dfdb(): not allocated!')
		end if
	end subroutine
	
	subroutine dfdb(self, r, dbv)
		type(dfdb_t) :: self
		real(DP), dimension(3), intent(in) :: r
		real(DP), dimension(:,:), pointer :: dbv
		
		type(atom_t), pointer :: atom
		type(basis_t), pointer :: basis
		real(DP), dimension(:,:), pointer :: dbop
		real(DP), dimension(:), pointer :: bfvec
		integer(I4) :: i, j, k, natoms, ncgto

		! Check if we already have the result
		if (r(1)==self%r(1) .and. r(2)==self%r(2) .and. r(3)==self%r(3)) then 
			dbv=>self%db
			return
		end if
		
		self%r=r
		natoms=get_natoms(self%mol)
		call mkdbop(self%dop,r, dbop)
		call bfeval(self%bfv, r, bfvec)
		
		j=1
		do k=1,natoms
			call get_atom(self%mol, k, atom)
			call get_basis(atom, basis)
			ncgto=get_ncgto(basis)
			do i=1,ncgto
				self%db(j,:)=dbop(:,k)*bfvec(j)
				j=j+1
			end do
		end do
		dbv=>self%db
	end subroutine 

end module
