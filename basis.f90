!
! Handles both segmented and general contractions
!
! Contractions are stored in individual arrays, which are "stored" in a
! vector containing all basis functions/contractions for a certain atom.
! Is would not be difficult to change the storage model to store everyting
! in one matrix, if need be.
!

module basis_class
	use globals_m
	use intgrl_m
	use gtodefs_m
	implicit none 

	public del_basis, init_basis, print_seg
	public get_natoms, get_ngto, get_ncgto, set_c2sop, get_c2sop
	public get_coord, get_symbol, get_basis, get_atom, get_nctr
	public get_contraction, get_ncomp, get_nccomp, get_nccgto
	public filter_screened, get_ctridx
	
	interface get_ncgto
		module procedure get_ncgto_mol, get_ncgto_bas
	end interface 

	private
contains
	subroutine init_basis(bt, molfil)
		type(molecule_t) :: bt
		character(*) :: molfil

		integer(I4) :: i, natoms
        type(atom_t), dimension(:), pointer :: atoms

		call read_intgrl(molfil,bt%atoms, natoms)
		call setup_gtos()
		bt%natoms=natoms
		call calc_basdim(bt)
		
		call msg_out('Normalizing basis')
		call msg_out(repeat('=', 45))
	    atoms=>bt%atoms
		do i=1,natoms
			call nl
			write(str_g, '(a,i4,a4)') 'Normalizing basis for atom: ', &
				i, atoms(i)%symbol
			call msg_out(str_g)
			call normalize(atoms(i)%basis)
		end do
		call nl
		call msg_out('Calculating screening thresholds')
		do i=1,natoms
			call setup_screening(atoms(i)%basis)
		end do

		do i=1,natoms
			call print_atom_data(atoms(i),i)
		end do
		write(str_g, 77) '  Total number of primitive  GTO''s ', bt%ngto
		call msg_out(str_g)
		write(str_g, 77) '  Total number of contracted GTO''s ', bt%ncgto
		call msg_out(str_g)
		call nl

77		format(a,i4)
	end subroutine
!
! Loop over ctr and set up screening thresholds
! 
	subroutine setup_screening(basis)
		type(basis_t), intent(inout) :: basis

		integer(I4) :: i,j,l
		real(DP) :: xp, min_xp, x, dist
		logical :: screen=.true.

		call getkw(input, 'screening', screen)
		if (.not.screen) then
			basis%thrs=1.d10
			return
		end if

		do i=1,basis%nctr
			min_xp=1.d+15
			x=1.d+15
			l=basis%ctr(i)%l
			do j=1,basis%ctr(i)%npf
				xp=basis%ctr(i)%xp(j)
				if ( xp < min_xp) min_xp=xp
			end do
			dist=0.d0
			do while (x > SCREEN_THRS)
				dist=dist+0.25
				x=dist**l*exp(-min_xp*dist**2)
			end do
			basis%thrs(i)=dist
		end do
	end subroutine 

!
! return an index vector with congtributing contractions after screening has
! been applied
!
	subroutine filter_screened(basis, r, idxv, n) 
		type(basis_t), intent(in) :: basis
		real(DP), dimension(3), intent(in) :: r
		integer(I4), dimension(:), intent(out) :: idxv
		integer(I4), intent(out) :: n

		real(DP) :: r2
		integer(I4) :: i,j

		r2=sqrt(sum(r**2))
		j=1
		do i=1,basis%nctr
			if (r2 <= basis%thrs(i) ) then
				idxv(j)=i
				j=j+1
			end if
		end do
		n=j-1
	end subroutine

!
! Loop over ctr and normalize
! 
	subroutine normalize(basis)
		type(basis_t), intent(inout) :: basis

		integer(I4) :: i

		do i=1,basis%nctr
			call norm_ctr(basis%ctr(i))
		end do
		write(str_g,'(a,i5,a)') '         Normalized', basis%nctr, &
		' contractions'
		call msg_out(str_g)
		call nl
	end subroutine 
	

! Loop over contractions and normalize
! 
! This is how things supposedly work:
! The cc:s are given for _normalized_ gaussians, so we first need to 
! normalize the gaussians, and then normalize the contraction. The 
! normalized contraction is then stored as new cc:s for the primitives.
!

	subroutine norm_ctr(ctr)
		type(contraction_t), intent(inout) :: ctr 
		
		integer(I4) :: l, m
		real(DP) :: n, t, e1, e2, c1, c2
		real(DP) :: j
		
		j=D1*(ctr%l+1)
		n=D0
		do l=1,ctr%npf  ! loop over number of primitives
			c1=ctr%cc(l)  
			e1=ctr%xp(l) 
			do m=1,l
				e2=ctr%xp(m)
				c2=ctr%cc(m)
				t=D2*sqrt(e1*e2)/(e1+e2)
				t=c1*c2*(t**(j+DP50))
				n=n+t
				if ( l /= m ) n=n+t    !  ab=ba <=> 2ab
			end do
		end do
		n=D1/sqrt(n)
		do l=1, ctr%npf
			c1=ctr%cc(l)
			e1=ctr%xp(l)
			ctr%ncc(l)=c1*n*(D4*e1)**(DP50*j+DP25)*(DP50/PII)**DP75
		end do
	end subroutine 

!
! debug function
!
	subroutine print_seg(ctr)
		type(contraction_t), intent(in) :: ctr

		integer(I4) :: i, tun

		call get_teletype_unit(tun)
		if (tun == DEVNULL) return
		do i=1,Ctr%npf
			write(tun, 20, advance='no') ctr%xp(i)
			write(tun, 21, advance='no') ctr%cc(i),' (',ctr%ncc(i),') '
			call nl
		end do
20		format(f20.10)
21		format(f13.10,a,f13.10,a)
	end subroutine 

!
! set the global variables ngto and ncgto so that we know
! how to dimension the basis function value vectors later. 
!
	subroutine calc_basdim(bt)
		type(molecule_t), intent(inout) :: bt

		integer(I4) :: k, l, n
		type(basis_t), pointer :: basis
		type(contraction_t), pointer :: ctr

		bt%ngto=0
		bt%ncgto=0
		bt%nccgto=0

		do n=1,bt%natoms
			basis=>bt%atoms(n)%basis

			basis%ngto=0
			basis%ncgto=0
			do k=1,basis%nctr
				ctr=>basis%ctr(k)

				basis%ngto=basis%ngto+ctr%npf*ctr%ncomp
				basis%ncgto=basis%ncgto+ctr%ncomp

				bt%ngto=bt%ngto+ctr%npf*ctr%ncomp
				bt%ncgto=bt%ncgto+ctr%ncomp
				bt%nccgto=bt%nccgto+ctr%nccomp
			end do
			basis%pos=1
			do k=2,basis%nctr
				basis%pos(k)=basis%pos(k-1)+basis%ctr(k-1)%nccomp
			end do
!            write(89,*) basis%pos
!            write(89,*)
		end do
	end subroutine

!
! Print atom and basis info
!
	subroutine print_atom_data(Atm, num)
		type(atom_t), intent(in) :: Atm
		integer(I4), intent(in) :: num

		type(basis_t), pointer :: abas
		integer(I4) :: i, j, k, npf, ncf

		abas=>atm%basis
		
		call nl
		write(str_g, 88) '      Atom number', num, ' (', atm%symbol, ')'
		call msg_out(str_g)
		write(str_g, 99) '     ----------------------'
		call msg_out(str_g)
		call nl
		write(str_g, 77) '      Nuclear charge:          ', atm%charge
		call msg_out(str_g)
		write(str_g, 66) '      Highest angular momentum:', abas%lmax,' (', &
			& shell_names(abas%lmax),')'
		call msg_out(str_g)
		call nl
		
		k=1
		do i=1,abas%nshells
			npf=0
			ncf=0
			do j=1,abas%nctrps(i)
				npf=npf+abas%ctr(k)%npf
				ncf=ncf+1
				!call print_seg(abas%ctr(k))
				!print *
				k=k+1
			end do
			write(str_g, 55) npf,' GTO''s and ', ncf, ' CGTO''s of ''',&
				& shell_names(i-1),''' type.'
			call msg_out(str_g)
		end do
		call nl

55      format(i3,a,i3,a,a1,a)
66		format(a,i3,a2,2a1)
77		format(a,f8.4)
88		format(a,i4,a,a2,a)
99		format(a)
	end subroutine 

!
! basis destructor
!
	subroutine del_ctr(abas)
		type(basis_t), intent(inout) :: abas

		type(contraction_t), pointer :: ctr
		integer(I4) :: j
		
		if ( associated(abas%ctr) ) then 
			do j=1,abas%nctr
                ctr=>abas%ctr(j)
                if (associated(ctr%xp)) then
					deallocate(ctr%xp)
					nullify(ctr%xp)
				end if
                if (associated(ctr%cc)) then
					deallocate(ctr%cc)
					nullify(ctr%cc)
				end if
                if (associated(ctr%ncc)) then
					deallocate(ctr%ncc)
					nullify(ctr%ncc)
				end if
			end do
			deallocate(abas%ctr)
			deallocate(abas%thrs)
			deallocate(abas%pos)
			nullify(abas%ctr)
		end if
	end subroutine 
	
!
! basis and atom list destructor 
! 
	subroutine del_basis(bt)
		type(molecule_t), intent(inout) :: bt

		integer(I4) :: i
		type(atom_t), pointer :: atm

		do i=1,bt%natoms
			atm=>bt%atoms(i)
			call del_ctr(atm%basis)
			deallocate(atm%basis)
			nullify(atm%basis)
		end do
		deallocate(bt%atoms)
		nullify(bt%atoms)
		call msg_info('Deallocated basis set and atom data')
		call nl
	end subroutine 

! 
! Access functions to get private data
!

	subroutine set_c2sop(mol, c2s)
		type(molecule_t) :: mol
		type(cao2sao_t), target :: c2s

		mol%c2s=>c2s
	end subroutine

	subroutine get_c2sop(mol, c2s)
		type(molecule_t) :: mol
		type(cao2sao_t), pointer :: c2s

		mol%c2s=>c2s
	end subroutine

	function get_natoms(mm) result(n)
		type(molecule_t) :: mm
		integer(I4) :: n

		n=mm%natoms
	end function 

	function get_ngto(mm) result(n)
		type(molecule_t) :: mm
		integer(I4) :: n

		n=mm%ncgto
	end function 
	
	function get_ncgto_mol(mm) result(n)
		type(molecule_t) :: mm
		integer(I4) :: n

		n=mm%ncgto
	end function 

	function get_nccgto(mm) result(n)
		type(molecule_t) :: mm
		integer(I4) :: n

		n=mm%nccgto
	end function 

	function get_ncgto_bas(mm) result(n)
		type(basis_t) :: mm
		integer(I4) :: n

		n=mm%ncgto
	end function 

	subroutine get_atom(mm, i, a)
		type(molecule_t), intent(in) :: mm
		integer(I4), intent(in) :: i
		type(atom_t), pointer :: a

		a=>mm%atoms(i)
	end subroutine

	subroutine get_basis(mm, bas)
		type(atom_t), intent(in) :: mm
		type(basis_t), pointer :: bas

		bas=>mm%basis
	end subroutine

	subroutine get_coord(mm, c)
		type(atom_t), intent(in) :: mm
		real(DP), dimension(3), intent(out) :: c

		c=mm%coord
	end subroutine

	subroutine get_symbol(mm, s)
		type(atom_t), intent(in) :: mm
		character(*), intent(out) :: s

		s=mm%symbol
	end subroutine

	subroutine get_contraction(mm, k, ctr)
		type(atom_t), intent(in) :: mm
		integer(I4), intent(in) :: k
		type(contraction_t), pointer :: ctr

		ctr=>mm%basis%ctr(k)
	end subroutine

	function get_ncomp(mm) result(n)
		type(contraction_t) :: mm
		integer(I4) :: n

		n=mm%ncomp
	end function 

	function get_nccomp(mm) result(n)
		type(contraction_t) :: mm
		integer(I4) :: n

		n=mm%nccomp
	end function 

	function get_nctr(mm) result(n)
		type(basis_t) :: mm
		integer(I4) :: n

		n=mm%nctr
	end function 

	function get_ctridx(bas, idx) result(ctridx)
		type(basis_t), intent(in) :: bas
		integer(I4), intent(in) :: idx
		integer(I4) :: ctridx

		ctridx=bas%pos(idx)
	end function

	! integrate s-funcs for debugging purpouses...
!	subroutine s_int(Ctr)
!		type(contraction_t), intent(in) :: Ctr 
!
!		real(DP) :: qq
!		integer(I4) :: i,j
!
!		qq=0.d0
!		do i=1,Ctr%npf
!			do j=1,Ctr%npf
!				qq=qq+Ctr%ncc(i)*Ctr%ncc(j)*(PI/(Ctr%xp(i)+Ctr%xp(j)))**1.5
!			end do
!		end do
!		print *, 's type integral:', qq
!		
!	end subroutine 
end module
