!
! Now you are surprised! This file has nothing to do with integral
! processing. Here are the routines to parse the 'INTGRL' file containing
! atom and basis set info, used by ACES2, Dalton...
!

module intgrl_module
    use globals_module
    use settings_module
    use teletype_module
    implicit none
    
    public read_intgrl

    private
contains
!
! Read in all atom and basis set info
!
    subroutine read_intgrl(fname, atoms, natoms)
        character(*) :: fname
        type(atom_t), dimension(:), pointer :: atoms
        integer(I4), intent(out) :: natoms
        
        character(6) :: intgrlkw
        character(9) :: turbokw
        integer(I4) :: i, ios
        type(atom_t), pointer :: a     ! short hand...
        type(basis_t), pointer :: b  ! short hand...
        integer(I4), dimension(:), allocatable :: bdim

        open(BASFD, file=trim(fname), status='old', action='read', &
             form='formatted', iostat=ios)

        if (ios /= 0 ) then
            call msg_error('read_intgrl(): open failed.')
            stop
        end if

        read(BASFD, '(a6)') intgrlkw
        if ( intgrlkw /= 'INTGRL' ) then
            call nl
            call msg_error('  Hmm... this doesn''t look like a ''INTGRL''&
                & file. Check your input.')
            call nl
            stop
        end if
        read(BASFD, '(a9)') turbokw
        if ( turbokw == 'TURBOMOLE' ) then
            call msg_info('Detected TURBOMOLE input')
            call nl
            is_turbomole=.true.
        end if
        
        read(BASFD, *) 

        read(BASFD, *) natoms 
        write(str_g, '(a,i4)') 'Number of atoms =', natoms
        call msg_out(str_g)
        call nl
        allocate(atoms(natoms))

        read(BASFD, *) 

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
            allocate(b%thrs(bdim(i))) ! mem for screening thrs. for atom i
            allocate(b%pos(bdim(i))) ! mem for pos idx thrs. for ctrs
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

        integer(I4) :: i
        character(80) :: tmp

        ! i is (am) a dummy..
        read(BASFD, *) atm%charge, i, atm%basis%nshells, &
            atm%basis%nctrps(1:atm%basis%nshells)
        atm%basis%lmax=atm%basis%nshells-1

        if (atm%basis%lmax > MAX_L) then
            write(str_g, '(a,i3)') ' Largest allowed l-quantum number in &
                &basis exceeded:', atm%basis%lmax 
            call msg_error(str_g)
            stop
        end if
        
        atm%basis%nctr=sum(atm%basis%nctrps(1:atm%basis%nshells))

        read(BASFD, '(a)') tmp
        atm%symbol=tmp(1:2)
        atm%id=tmp(3:4)
        read(tmp(5:), *) atm%coord 
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
                if (settings%use_spherical) then
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
        
        read(BASFD, *) ctr%npf, i

        if (i > 1) then
            call msg_error('General contractions not supported!!!')
            stop
        end if

        allocate(ctr%xp(ctr%npf)) 
        allocate(ctr%cc(ctr%npf)) 
        allocate(ctr%ncc(ctr%npf))
        
        do i=1,ctr%npf
            read(BASFD, *) ctr%xp(i), ctr%cc(i)
        end do
    end subroutine 

    subroutine read_contraction2(b, idx, ncf)
        type(basis_t), intent(inout) :: b
        integer(I4), intent(in) :: idx
        integer(I4), intent(out) :: ncf

        integer(I4) :: i, j, l, ncomp, nccomp, npf, nncf
        real(DP) :: xp
        real(DP), dimension(:), allocatable :: cc
        type(contraction_t), pointer :: ctr
        
        read(BASFD, *) npf, ncf
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
        
        do i=1,npf
            read(BASFD, *) xp, cc(1:ncf)
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
        
        integer(I4) :: i, j, k, n
        integer(I4) :: ncf, npf, nsh
        integer(I4), dimension(MAX_L+1) :: nctrps
        real(DP) :: rx
        real(DP), dimension(100) :: slask
        integer(I4) :: ix

        bdim=0
        do n=1,size(bdim) !natoms
            read(BASFD, *) rx, ix, nsh, nctrps(1:nsh)
            read(BASFD, *)
            do i=1,nsh
                do j=1,nctrps(i)
                    read(BASFD, *) npf, ncf 
                    bdim(n)=bdim(n)+ncf
                    do k=1,npf
                        read(BASFD,*) rx, slask(1:ncf)
                    end do
                end do
            end do
        end do
        
        rewind(BASFD)
        do n=1,5
            read(BASFD,*)
        end do
    end subroutine

end module

! vim:et:sw=4:ts=4
