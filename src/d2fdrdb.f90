
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
		real(DP), dimension(:,:), pointer :: d2
		integer(I4) :: ptn
	end type

	private
	real(DP), dimension(3) :: dbov
	real(DP) :: ror1, ror2, ror3
	
contains

	subroutine init_d2fdrdb(self, mol)
		type(d2fdrdb_t) :: self
		type(molecule_t), target :: mol

		nullify(self%d2)

		if (associated(self%d2)) then
			call msg_warn('init_selfdrdb(): already allocated!')
		else
			allocate(self%d2(get_ncgto(mol),9))
		end if
		self%mol=>mol
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

	subroutine d2fdrdb(self, r,  bfvec, drvec, dbop, selfv)
		type(d2fdrdb_t) :: self
		real(DP), dimension(:), intent(in) :: r
		real(DP), dimension(:,:), intent(in) :: drvec, dbop
		real(DP), dimension(:), intent(in) :: bfvec
		real(DP), dimension(:,:), pointer :: selfv
	
		type(atom_t), dimension(:), pointer :: atoms
		integer(I4) :: i, j, k, natoms, nctr, ncomp, idx
		real(DP), dimension(3) :: coord
		type(atom_t), pointer :: atom
		type(basis_t), pointer :: basis
		type(contraction_t), pointer :: ctr

		integer(I4), dimension(99) :: posvec
		integer(I4) :: l, idx1
		
		natoms=get_natoms(self%mol)
		
		idx=1
		self%d2=0.d0
		do i=1,natoms
			call get_atom(self%mol,i,atom)
			dbov=dbop(:,i)
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
					self%d2(idx,1)=drvec(idx,1)*dbov(1)                  ! dx
					self%d2(idx,2)=drvec(idx,2)*dbov(1)+ror3*bfvec(idx)  ! dy
					self%d2(idx,3)=drvec(idx,3)*dbov(1)-ror2*bfvec(idx)  ! dz

					! dBy
					self%d2(idx,4)=drvec(idx,1)*dbov(2)-ror3*bfvec(idx)
					self%d2(idx,5)=drvec(idx,2)*dbov(2)
					self%d2(idx,6)=drvec(idx,3)*dbov(2)+ror1*bfvec(idx)

					! dBz
					self%d2(idx,7)=drvec(idx,1)*dbov(3)+ror2*bfvec(idx)
					self%d2(idx,8)=drvec(idx,2)*dbov(3)-ror1*bfvec(idx)
					self%d2(idx,9)=drvec(idx,3)*dbov(3)
					idx=idx+1
				end do
			end do
		end do
		selfv=>self%d2
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
