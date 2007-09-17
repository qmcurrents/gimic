
module dens_m
	use globals_m
	use basis_m
	use reorder_m
	implicit none

	public init_dens, read_dens, set_dens, set_pdens
	public get_dens, get_pdens, del_dens, dens_t
	public read_modens
	
	type dens_t
		type(molecule_t), pointer :: mol
        real(DP), dimension(:,:,:), pointer :: da, db
		logical :: pdens_p=.true.
		integer(I4) :: spin
	end type
	
	private
	
	real(DP), dimension(:,:,:), pointer :: dens
contains
	subroutine init_dens(dd, mol, modens_p)
		type(dens_t) :: dd
		type(molecule_t), target :: mol
		logical, optional :: modens_p

		integer(I4) :: ncgto

		dd%mol=>mol
		ncgto=get_ncgto(mol)

		if (present(modens_p)) then
			if (modens_p) then
				dd%pdens_p=.false.
			end if
		end if

		if (uhf_p) then
			dd%spin=2
		else
			dd%spin=1
		end if

		if (dd%pdens_p) then
			allocate(dd%da(ncgto,ncgto,0:3))
			if (uhf_p) allocate(dd%db(ncgto,ncgto,0:3))
		else
			allocate(dd%da(ncgto,ncgto,0:0))
			if (uhf_p) allocate(dd%db(ncgto,ncgto,0:0))
		end if

	end subroutine

	subroutine read_dens(dd)
		type(dens_t), intent(inout) :: dd
		character(80) :: xdens_file

		integer(I4) :: b, mo, ispin
		type(reorder_t) :: bofh

		if ( .not.associated(dd%da) ) then
			call msg_error('read_dens(): dens not allocated!')
			call exit(1)
		end if
		if ( uhf_p.and..not.associated(dd%db) ) then
			call msg_error('read_dens(): beta dens not allocated!')
			call exit(1)
		end if

		call getkw(input, 'density', xdens_file)

		open(XDFD, file=xdens_file, status='old', err=42)

		dens=>dd%da
		do ispin=1,dd%spin
			if (ispin == 2) dens=>dd%db
			if (dd%pdens_p) then
				do b=0,3
					read(XDFD,*) dens(:,:,b)
				end do
			else
				read(XDFD,*) dens(:,:,0)
			end if
		end do
		if (uhf_p) then
			call  msg_info('scaling perturbed densities by 0.d5')
			dd%da(:,:,1:3)= dd%da(:,:,1:3)/2.d0
			dd%db(:,:,1:3)= dd%db(:,:,1:3)/2.d0
		end if

		if (turbomole_p) then
			call init_reorder(bofh, dd%mol)
			call msg_info('Reordering densities [TURBOMOLE]')
			call turbo_reorder(bofh)
			call reorder_dens(bofh, dd)
			call del_reorder(bofh)
		end if

        close(XDFD)
        return

42      bert_is_evil=.true. 
        call msg_error('Density file not found, all densities set to 1')
		call nl
		dens=>dd%da
		do ispin=1,dd%spin
			if (ispin == 2) dens=>dd%db
			do b=0,3
				dens(:,:,b)=1.d0
			end do
		end do
	end subroutine

	subroutine del_dens(dd)
		type(dens_t) :: dd

		if (associated(dd%da)) then
			deallocate(dd%da)
		else
			call msg_warn('del_dens(): not allocated!')
		end if
		
		if (uhf_p.and.associated(dd%db)) then
			deallocate(dd%db)
		else
			call msg_warn('del_dens(): not allocated!')
		end if
		
	end subroutine

	subroutine get_dens(dd,a, spin)
		type(dens_t) :: dd
		real(DP), dimension(:,:), pointer :: a
		integer(I4), optional :: spin

		dens=>dd%da
		if (present(spin)) then
			if (spin == spin_b) dens=>dd%db
		end if
		
		a=>dens(:,:,0)
	end subroutine  

	subroutine set_dens(dd,a,spin)
		type(dens_t) :: dd
		real(DP), dimension(:,:) :: a
		integer(I4), optional :: spin
		
		dens=>dd%da
		if (present(spin)) then
			if (spin == spin_b) dens=>dd%db
		end if
		dens(:,:,0)=a
	end subroutine  

	subroutine  get_pdens(dd,b,a,spin)
		type(dens_t) :: dd
		integer(I4), intent(in) :: b
		real(DP), dimension(:,:), pointer :: a
		integer(I4), optional :: spin

		dens=>dd%da
		if (present(spin)) then
			if (spin == spin_b) dens=>dd%db
		end if
		a=>dens(:,:,b)
	end subroutine  

	subroutine  set_pdens(dd,b,a,spin)
		type(dens_t) :: dd
		integer(I4), intent(in) :: b
		real(DP), dimension(:,:) :: a
		integer(I4), optional :: spin

		dens=>dd%da
		if (present(spin)) then
			if (spin == spin_b) dens=>dd%db
		end if
		dens(:,:,b)=a
	end subroutine  

	subroutine reorder_dens(bofh, dd)
		type(reorder_t) :: bofh
		type(dens_t) :: dd

		integer(I4) :: i,b,p,ncgto,  ispin

		ncgto=get_ncgto(dd%mol)
		p=0
		if (dd%pdens_p) then
			p=3
		end if

		dens=>dd%da
		do ispin=1,dd%spin
			if (ispin == 2) dens=>dd%db
			do b=0,p
				call reorder_cols(bofh, dens(:,:,b))
			end do
			do i=1,ncgto
				do b=0,p
					call reorder_vec(bofh, dens(:,i,b))
				end do
			end do
		end do
	end subroutine

	subroutine moco(dd, mos, spin)
		type(dens_t) :: dd
		real(DP), dimension(:,:) :: mos
		integer(I4), optional :: spin

		integer(I4), dimension(2) :: moran
		integer(I4) :: ncgto, i, a, b

		dens=>dd%da
		if (present(spin)) then
			if (spin == spin_b) then
				if (.not.uhf_p) then
					call msg_error('gimic: moco(): &
					&invalid spin for closed shell')
					stop
				end if
				dens=>dd%db
			end if
		end if
		ncgto=get_ncgto(dd%mol)
		moran=0
		call getkw(input, 'edens.mos', moran)

		if (moran(1) < 1 .or. moran(2) < 1) then
			call msg_error('Invalid MO range!')
			call exit(1) 
		end if
		if (moran(1) > moran(2) ) then
			call msg_error('Invalid MO range!')
			call exit(1) 
		end if

		call nl
		call msg_note('Constructing closed shell AO density matrix')
		write(str_g, '(a,i4,a,i4)') 'Occupancy is', moran(1), ' -', moran(2)
		call msg_note(str_g)
		call nl
		do a=1,ncgto
			do b=1,ncgto
				do i=moran(1),moran(2)
					dens(a,b,0)=dens(a,b,0)+mos(a,i)*mos(b,i)
				end do
			end do
		end do
		dens(:,:,0)=dens(:,:,0)*2.d0
	end subroutine

	subroutine read_modens(dd)
		type(dens_t) :: dd
		
		integer(4) :: n, i,j
		real(DP), dimension(:,:), allocatable :: mos
		type(reorder_t) :: bofh
		character(BUFLEN) :: mofile

		if (.not.turbomole_p) then
			call read_dens(dd)
			return
		end if

		dens=>dd%da
		dens(:,:,0)=D0
		call getkw(input, 'edens.mofile', mofile)
		if (trim(mofile) == '') then
			call read_dens(dd)
		end if
		open(XDFD, file=trim(mofile), status='old', err=42)
		read(XDFD, *) 
		read(XDFD, *) 
		read(XDFD, *) 

		n=get_ncgto(dd%mol)
		allocate(mos(n,n))
		
		do i=1,n
			read(XDFD, *)
			read(XDFD, '(4f20.2)') mos(:,i)
		end do
		close(XDFD)

		call moco(dd, mos)

		call init_reorder(bofh, dd%mol)
		call msg_info('Reordering densities [TURBOMOLE]')
		call turbo_reorder(bofh)
		call reorder_dens(bofh, dd)
		call del_reorder(bofh)
		
		deallocate(mos)
		return

42      bert_is_evil=.true. 
		call msg_error('MO file not found, densitiy set to 1.d0')
		dens(:,:,0)=1.d0
	end subroutine
end module 
