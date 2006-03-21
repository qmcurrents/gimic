!
! $Id$
!

module dfdb_m
	use globals_m
	use basis_m
	use dboper_m
	use bfeval_m
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
subroutine init_dfdb(dbt, mol, dop, bfv)
		type(dfdb_t) :: dbt
		type(molecule_t), target :: mol
		type(dbop_t), target :: dop
		type(bfeval_t), target :: bfv
		if (associated(dbt%db)) then
			call msg_warn('init_dfdb(): already allocated!')
		else
			allocate(dbt%db(get_ncgto(mol),3))
		end if
		dbt%mol=>mol
		dbt%dop=>dop
		dbt%bfv=>bfv
		dbt%r=INITRV
	end subroutine

	subroutine del_dfdb(dbt)
		type(dfdb_t) :: dbt
		if (associated(dbt%db)) then
			deallocate(dbt%db)
			nullify(dbt%db)
		else
			call msg_warn('del_dfdb(): not allocated!')
		end if
	end subroutine
	
	subroutine dfdb(dbt, r, dbv)
		type(dfdb_t) :: dbt
		real(DP), dimension(3), intent(in) :: r
		real(DP), dimension(:,:), pointer :: dbv
		
		type(atom_t), pointer :: atom
		type(basis_t), pointer :: basis
		real(DP), dimension(:,:), pointer :: dbop
		real(DP), dimension(:), pointer :: bfvec
		integer(I4) :: i, j, k, natoms, ncgto

		! Check if we already have the result
		if (r(1)==dbt%r(1) .and. r(2)==dbt%r(2) .and. r(3)==dbt%r(3)) then 
			dbv=>dbt%db
			return
		end if
		
		dbt%r=r
		natoms=get_natoms(dbt%mol)
		call mkdbop(dbt%dop,r, dbop)
		call bfeval(dbt%bfv, r, bfvec)
		
		j=1
		do k=1,natoms
			call get_atom(dbt%mol, k, atom)
			call get_basis(atom, basis)
			ncgto=get_ncgto(basis)
			do i=1,ncgto
				dbt%db(j,:)=dbop(:,k)*bfvec(j)
				j=j+1
			end do
		end do
		dbv=>dbt%db
	end subroutine 

end module
