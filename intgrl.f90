!
! $Id$
!
! Now you are surprised! This file has nothing to do with integral
! processing. Here are the routines to parse the 'INTGRL' file containing
! atom and basis set info, used by ACES2, Dalton...
!

module intgrl_m
	use globals_m
	use teletype_m
	implicit none
	
	public read_intgrl

	integer(I4) :: line
	character(MAX_LINE_LEN), dimension(:), pointer :: molbuf
	private
contains
!
! Read in all atom and basis set info
!
	subroutine read_intgrl(mbuf, atoms, natoms)
		character(MAX_LINE_LEN), dimension(:), target :: mbuf
		type(atom_t), dimension(:), pointer :: atoms
		integer(I4), intent(out) :: natoms
		
		character(6) :: intgrlkw
		character(9) :: turbokw
		integer(I4) :: ngto, i, ios, n, fgetc
		type(atom_t), pointer :: a     ! short hand...
		type(basis_t), pointer :: b  ! short hand...
		integer(I4), dimension(:), allocatable :: bdim

		molbuf=>mbuf
		line=1

		read(molbuf(line), '(a6)') intgrlkw
		line=line+1
		if ( intgrlkw /= 'INTGRL' ) then
			call nl
			call msg_error('  Hmm... this doesn''t look like a ''INTGRL''&
				& file. Check your input.')
			call nl
			stop
		end if
		read(molbuf(line), '(a9)') turbokw
		line=line+1
		if ( turbokw == 'TURBOMOLE' ) then
			call msg_info('Detected TURBOMOLE input')
			call nl
			turbomole_p=.true.
		end if
		
		line=line+1

		read(molbuf(line), *) natoms 
		line=line+1
		write(str_g, '(a,i4)') 'Number of atoms =', natoms
		call msg_out(str_g)
		call nl
		allocate(atoms(natoms))

		line=line+1

		! ok, this is a bit complicated, as we are transforming a (possibly)
		! generally contracted basis into a segmented basis. 
		! If the basis is generally contracted we need to sort densities later
		! on... sigh.
		allocate(bdim(natoms))
		! get actual number of contractions
		call get_basdim(bdim)

		do i=1,natoms
			a=>atoms(i)
			allocate(a%basis)
			b=>a%basis
			call read_atom(a)            ! read general atom/basis info
			allocate(b%ctr(bdim(i))) ! mem for all contrs. for atom i
			call read_segs(b)    	! atom i
			b%nctr=bdim(i)
		end do

		deallocate(bdim)
	end subroutine 

!
! Read atom data (charge, coords, etc.), and some of the basis set defs
! 
	subroutine read_atom(atm)
		type(atom_t), intent(inout) :: atm

		integer(I4) :: i, n, ios

		! i is (am) a dummy..
		read(molbuf(line), *) atm%charge, i, atm%basis%nshells, &
			atm%basis%nctrps(1:atm%basis%nshells)
		line=line+1
		atm%basis%lmax=atm%basis%nshells-1

		if (atm%basis%lmax > MAX_L) then
			write(str_g, '(a,i3)') ' Largest allowed l-quantum number in &
				&basis exceeded:', atm%basis%lmax 
			call msg_error(str_g)
			stop
		end if
		
		atm%basis%nctr=sum(atm%basis%nctrps(1:atm%basis%nshells))

		atm%symbol=molbuf(line)(1:2)
		atm%id=molbuf(line)(3:4)
		read(molbuf(line)(5:), *) atm%coord
		line=line+1
		
	end subroutine 

!
!  Loop over all segment of contractions and fetch primitives and c-coefs.
!
	subroutine read_segs(bas)
		type(basis_t), intent(inout) :: bas
		
		type(contraction_t), pointer :: ctr
		integer(I4) :: i, j, k, ncf, nctrps
		
		k=1
        do i=1,bas%nshells
			nctrps=0
            do j=1,bas%nctrps(i)
				ctr=>bas%ctr(k)
				! i => l-quantum number +1
				ctr%l=i-1 ! save l-qnum 
				! store number of components 
				if (spherical) then
					ctr%ncomp=2*ctr%l+1
				else
					ctr%ncomp=(ctr%l+1)*(ctr%l+2)/2
				end if
				ctr%nccomp=(ctr%l+1)*(ctr%l+2)/2
!                call read_contraction(ctr)
				call read_contraction2(bas, k, ncf)
				nctrps=nctrps+ncf
                k=k+ncf
            end do
			bas%nctrps(i)=nctrps
		end do
	end subroutine 
!
! Read one contraction segment
!
	subroutine read_contraction(ctr)
		type(contraction_t), intent(inout) :: ctr

		integer(I4) :: i, j, ios
		
		read(molbuf(line), *) ctr%npf, i
		line=line+1

		if (i > 1) then
			call msg_error('General contractions not supported!!!')
			stop
		end if

		allocate(ctr%xp(ctr%npf)) 
		allocate(ctr%cc(ctr%npf)) 
		allocate(ctr%ncc(ctr%npf))
		
		do i=1,ctr%npf
			read(molbuf(line), *) ctr%xp(i), ctr%cc(i)
		line=line+1
		end do
	end subroutine 

	subroutine read_contraction2(b, idx, ncf)
		type(basis_t), intent(inout) :: b
		integer(I4), intent(in) :: idx
		integer(I4), intent(out) :: ncf

		integer(I4) :: i, j, l, ncomp, nccomp, npf, nlines, nncf
		real(DP) :: xp
		real(DP), dimension(:), allocatable :: cc
		type(contraction_t), pointer :: ctr
		
		read(molbuf(line), *) npf, ncf
		line=line+1
		allocate(cc(ncf))

		do i=1,ncf
			ctr=>b%ctr(idx+i-1)
			allocate(ctr%xp(npf)) 
			allocate(ctr%cc(npf)) 
			allocate(ctr%ncc(npf))
		end do

		ncomp=b%ctr(idx)%ncomp
		nccomp=b%ctr(idx)%ncomp
		l=b%ctr(idx)%l

		do i=1,ncf
			ctr=>b%ctr(idx+i-1)
			ctr%ncf=1
			ctr%npf=npf
			ctr%ncomp=ncomp
			ctr%nccomp=nccomp
			ctr%l=l
		end do

		b%ctr(idx)%ncf=ncf ! needed for degeneralization...
		
		! this is hell... 
		nlines=get_nctr_lines(ncf)-1
		do i=1,npf
			if (ncf > 3) then
				read(molbuf(line), *) xp, cc(1:3)
			else
				read(molbuf(line), *) xp, cc(1:ncf)
			end if
			line=line+1
			l=0
			nncf=ncf-3
			do j=1,nlines
				l=l+4
				if (nncf > 4) then
					read(molbuf(line), *) cc(l:l+3) 
				else
					read(molbuf(line), *) cc(l:l+nncf-1) 
				end if
				line=line+1
				nncf=nncf-4
			end do
			do j=1,ncf
				ctr=>b%ctr(idx+j-1)
				ctr%xp(i)=xp
				ctr%cc(i)=cc(j)
			end do
		end do
		deallocate(cc)
	end subroutine 


	! figure out how many lines the cc coefs span
	function get_nctr_lines(ncf) result(nlines)
		integer(I4), intent(in) :: ncf
		integer(I4) :: nlines

		integer(I4) :: n, modn

		nlines=1
		if (ncf > 3) then
			n=ncf-3
			modn=mod(n,4)
			nlines=nlines+(n-modn)/4
			if (modn > 0) nlines=nlines+1
		end if
	end function

	subroutine get_basdim(bdim)
		integer(I4), dimension(:), intent(out) :: bdim
		
		integer(I4) :: i, j, k, n, nlines
		integer(I4) :: ncf, npf, nsh
		integer(I4), dimension(MAX_L+1) :: nctrps
		real(DP) :: rx, ix

		bdim=0
		do n=1,size(bdim) !natoms
			read(molbuf(line), *) rx, ix, nsh, nctrps(1:nsh)
			line=line+1
			line=line+1
			do i=1,nsh
				do j=1,nctrps(i)
					read(molbuf(line), *) npf, ncf 
					line=line+1
					nlines=get_nctr_lines(ncf)
					bdim(n)=bdim(n)+ncf
					do k=1,npf
						line=line+nlines
					end do
				end do
			end do
		end do
		
		line=6       ! get back to were we were called
	end subroutine

end module
