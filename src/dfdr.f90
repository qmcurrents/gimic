
module dfdr_class
	use globals_m
	use basis_class
	use caos_m
	use gtodefs_m
	use cao2sao_class
	implicit none

	public init_dfdr, del_dfdr, dfdr, dfdr_t

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

	subroutine init_dfdr(self, mol)
		type(dfdr_t) :: self
		type(molecule_t), target :: mol

		allocate(self%dr(get_nccgto(mol),3))
		if (spherical) allocate(self%sdr(get_ncgto(mol),3))
		self%mol=>mol
!        call dfdr(self, (/0.d0, 0.d0, 0.d0/), self%dr)
	end subroutine

	subroutine del_dfdr(self)
		type(dfdr_t) :: self

		deallocate(self%dr)
		if (spherical) deallocate(self%sdr)
	end subroutine

	subroutine dfdr(self, r, drv)
		type(dfdr_t) :: self
		real(DP), dimension(3), intent(in) :: r
		real(DP), dimension(:,:), pointer :: drv

		integer(I4), dimension(99) :: posvec
		integer(I4) :: idx1, idx2
		
		natoms=get_natoms(self%mol)
		
		idx2=0
		self%dr=0.d0
		do i=1,natoms
			call get_atom(self%mol,i,atom)
			call get_coord(atom, coord)
			rr=r-coord
			call get_basis(atom, basis)
			call filter_screened(basis, rr, posvec, nctr)
			do k=1,nctr 
				j=posvec(k)
				call get_contraction(atom, j, ctr)
				idx=idx2+get_ctridx(basis, j)
				do axis=1,3
					call dcgto(rr, ctr, axis, self%dr(idx:, axis))
				end do
			end do
			idx2=idx2+get_ncgto(basis)
		end do
		if (spherical) then
			do axis=1,3
				call cao2sao(self%mol%c2s, self%dr(:, axis), self%sdr(:,axis))
			end do
			drv=>self%sdr
		else
			drv=>self%dr
		end if
	end subroutine 
end module
