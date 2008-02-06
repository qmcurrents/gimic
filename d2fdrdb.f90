
module d2fdrdb_class
	use globals_m
	use basis_class
	use dbop_class
	use bfeval_class
	use dfdr_class
	use teletype_m
	implicit none

	public init_d2fdrdb, del_d2fdrdb, d2fdrdb, d2fdrdb_t

	type d2fdrdb_t
		type(molecule_t), pointer :: mol
		type(dbop_t), pointer :: dop
		type(dfdr_t), pointer :: dfr
		type(bfeval_t), pointer :: bfv
		real(DP), dimension(:,:), pointer :: d2
		real(DP), dimension(3) :: r
	end type

	private
	real(DP), dimension(:,:), pointer :: drvec, dbop
	real(DP), dimension(:), pointer :: bfvec
	real(DP), dimension(3) :: dbov
	real(DP) :: ror1, ror2, ror3
	
contains

	subroutine init_d2fdrdb(self, mol, dop, dfr, bfv)
		type(d2fdrdb_t) :: self
		type(molecule_t), target :: mol
		type(dbop_t), target :: dop
		type(dfdr_t), target :: dfr
		type(bfeval_t), target :: bfv

		nullify(self%d2)

		if (associated(self%d2)) then
			call msg_warn('init_selfdrdb(): already allocated!')
		else
			allocate(self%d2(get_ncgto(mol),9))
		end if
		self%mol=>mol
		self%dop=>dop
		self%dfr=>dfr
		self%bfv=>bfv
		self%r=INITRV
	end subroutine

	subroutine del_d2fdrdb(self)
		type(d2fdrdb_t) :: self

		if (associated(self%d2)) then
			deallocate(self%d2)
			nullify(self%d2)
			nullify(self%mol)
		else
			call msg_warn('del_selfdrdb(): not allocated!')
		end if
	end subroutine

	subroutine d2fdrdb(self, r, selfv)
		type(d2fdrdb_t) :: self
		real(DP), dimension(:), intent(in) :: r
		real(DP), dimension(:,:), pointer :: selfv
	
		type(atom_t), dimension(:), pointer :: atoms
		integer(I4) :: i, j, k, natoms, nctr, ncomp, idx
		real(DP), dimension(3) :: coord
		type(atom_t), pointer :: atom
		type(basis_t), pointer :: basis
		type(contraction_t), pointer :: ctr

		integer(I4), dimension(99) :: posvec
		integer(I4) :: l, idx1, idx2
		
		! Check if we already have the result
		if (r(1)==self%r(1) .and. r(2)==self%r(2) .and. r(3)==self%r(3)) then 
			selfv=>self%d2
			return
		end if
		
		self%r=r
		natoms=get_natoms(self%mol)
		call mkdbop(self%dop, r, dbop)
		call bfeval(self%bfv, r, bfvec)
		call dfdr(self%dfr, r, drvec)
		
		idx=1
		idx2=0
		self%d2=0.d0
		do i=1,natoms
			call get_atom(self%mol,i,atom)
			dbov=dbop(:,i)
			call get_coord(atom, coord)
			ror1=coord(1)
			ror2=coord(2)
			ror3=coord(3)
			call get_basis(atom, basis)
			call filter_screened(basis, coord, posvec, nctr)
			do l=1,nctr
				j=posvec(l)
				call get_contraction(atom, j, ctr)
				idx=idx2+get_ctridx(basis, j)
				ncomp=get_ncomp(ctr)
				do k=1,ncomp
					call self_comp(self%d2, idx)
					idx=idx+1
				end do
			end do
			idx2=idx2+get_ncgto(basis)
		end do
		selfv=>self%d2
	end subroutine 

	subroutine self_comp(d2,idx)
		real(DP), dimension(:,:) :: d2
		integer(I4), intent(in) :: idx

		real(DP) :: bfv
		real(DP) :: drvec1, drvec2, drvec3
		real(DP) :: dbov1, dbov2, dbov3

		bfv=bfvec(idx)
		drvec1=drvec(idx,1)
		drvec2=drvec(idx,2)
		drvec3=drvec(idx,3)
		dbov1=dbov(1)
		dbov2=dbov(2)
		dbov3=dbov(3)
		
		! dBx
		d2(idx,1)=drvec1*dbov1           ! dx
		d2(idx,2)=drvec2*dbov1+ror3*bfv  ! dy
		d2(idx,3)=drvec3*dbov1-ror2*bfv  ! dz
		
		! dBy
		d2(idx,4)=drvec1*dbov2-ror3*bfv
		d2(idx,5)=drvec2*dbov2
		d2(idx,6)=drvec3*dbov2+ror1*bfv
		
		! dBz
		d2(idx,7)=drvec1*dbov3+ror2*bfv
		d2(idx,8)=drvec2*dbov3-ror1*bfv
		d2(idx,9)=drvec3*dbov3
	end subroutine

	subroutine print_selfvec(foo)
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
	
end module
