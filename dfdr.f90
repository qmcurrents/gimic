
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
		real(DP), dimension(3) :: r
	end type

	private
	
	integer(I4) :: i, j, k, natoms, nctr, nccomp, axis
	type(atom_t), pointer :: atom
	type(basis_t), pointer :: basis
	type(contraction_t), pointer :: ctr
	real(DP), dimension(3) :: rr, coord
	integer(I4) :: idx
		
contains

	subroutine init_dfdr(dfr, mol)
		type(dfdr_t) :: dfr
		type(molecule_t), target :: mol

		nullify(dfr%dr)

		if (associated(dfr%dr)) then
			call msg_warn('init_dfdr(): already allocated!')
		else
			allocate(dfr%dr(get_nccgto(mol),3))
			if (spherical) allocate(dfr%sdr(get_ncgto(mol),3))
		end if
		dfr%mol=>mol
		dfr%r=INITRV
	end subroutine

	subroutine del_dfdr(dfr)
		type(dfdr_t) :: dfr
		if (associated(dfr%dr)) then
			deallocate(dfr%dr)
			if (spherical) deallocate(dfr%sdr)
			nullify(dfr%dr)
			nullify(dfr%sdr)
		else
			call msg_warn('del_dfdr(): not allocated!')
		end if
	end subroutine

	subroutine dfdr(dfr, r, drv)
		type(dfdr_t) :: dfr
		real(DP), dimension(3), intent(in) :: r
		real(DP), dimension(:,:), pointer :: drv

		integer(I4), dimension(99) :: posvec
		integer(I4) :: idx1, idx2

		! Check if we already have the result
		if (r(1)==dfr%r(1) .and. r(2)==dfr%r(2) .and. r(3)==dfr%r(3)) then 
			drv=>dfr%dr
			return
		end if
		
		dfr%r=r
		natoms=get_natoms(dfr%mol)
		
		idx=1
		idx2=0
		dfr%dr=0.d0
		do i=1,natoms
			call get_atom(dfr%mol,i,atom)
			call get_coord(atom, coord)
			rr=r-coord
			call get_basis(atom, basis)
			call filter_screened(basis, rr, posvec, nctr)
			do k=1,nctr 
				j=posvec(k)
				call get_contraction(atom, j, ctr)
				idx=idx2+get_ctridx(basis, j)
				do axis=1,3
					call dcgto(rr, axis, ctr, &
					dfr%dr(idx:, axis))
				end do
			end do
			idx2=idx2+get_ncgto(basis)
		end do
		if (spherical) then
			do axis=1,3
				call cao2sao(dfr%mol%c2s, dfr%dr(:, axis), dfr%sdr(:,axis))
			end do
			drv=>dfr%sdr
		else
			drv=>dfr%dr
		end if
	end subroutine 
end module
