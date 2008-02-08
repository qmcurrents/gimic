!
! Basis function evaluator
!

module bfeval_class
	use globals_m
	use basis_class
	use gtodefs_m
	use caos_m
	use cao2sao_class
	implicit none

	public bfeval, init_bfeval, del_bfeval, bfeval_t

	type bfeval_t
		type(molecule_t), pointer :: mol
		real(DP), dimension(:), pointer :: bf
		real(DP), dimension(:), pointer :: sbf
		real(DP), dimension(3) :: r
	end type

	private
	
	integer(I4) :: i, j, k, natoms, nctr, ncomp, nccomp, idx
	type(atom_t), pointer :: atom
	type(basis_t), pointer :: basis
	type(contraction_t), pointer :: ctr
	real(DP), dimension(3) :: rr, coord
		
contains

	subroutine init_bfeval(self, mol)
		type(bfeval_t) :: self
		type(molecule_t), target :: mol
		
		allocate(self%bf(get_ncgto(mol)))
		if (spherical) allocate(self%sbf(get_nccgto(mol)))
		self%mol=>mol
		self%r=INITRV
	end subroutine

	subroutine del_bfeval(self)
		type(bfeval_t) :: self

		deallocate(self%bf)
		if (spherical) deallocate(self%sbf)
		nullify(self%sbf)
		nullify(self%mol)
	end subroutine
	
	subroutine bfeval(self, r, ans)
		type(bfeval_t) :: self
		real(DP), dimension(3), intent(in) :: r
		real(DP), dimension(:), pointer :: ans

		integer(I4), dimension(99) :: posvec
		integer(I4) :: idx1, idx2

		! Check if we already have the result
		if (r(1)==self%r(1) .and. r(2)==self%r(2) .and. r(3)==self%r(3)) then 
			ans=>self%bf
			return
		end if
		
		self%r=r
		natoms=get_natoms(self%mol)
		
		idx=1
		idx2=0
		self%bf=0.d0
		do i=1,natoms
			call get_atom(self%mol,i,atom)
			call get_coord(atom, coord)
			call get_basis(atom, basis)
			rr=r-coord
			call filter_screened(basis, rr, posvec, nctr)
!            nctr=get_nctr(basis)
!            do j=1,nctr 
			do k=1,nctr 
				j=posvec(k)
				call get_contraction(atom, j, ctr)
				idx=idx2+get_ctridx(basis, j)
				call cgto(rr, ctr, self%bf(idx:))
			end do
			idx2=idx2+get_ncgto(basis)
		end do
		if (spherical) then
			call cao2sao(self%mol%c2s, self%bf, self%sbf)
			ans=>self%sbf
		else
			ans=>self%bf
		end if
	end subroutine 

end module
