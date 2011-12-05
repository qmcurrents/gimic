
module dens_class
    use globals_m
    use basis_class
    use reorder_class
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
    subroutine init_dens(self, mol, modens_p)
        type(dens_t) :: self
        type(molecule_t), target :: mol
        logical, optional :: modens_p

        integer(I4) :: ncgto

        self%mol=>mol
        ncgto=get_ncgto(mol)

        if (present(modens_p)) then
            if (modens_p) then
                self%pdens_p=.false.
            end if
        end if

        if (uhf_p) then
            self%spin=2
        else
            self%spin=1
        end if

        if (self%pdens_p) then
            allocate(self%da(ncgto,ncgto,0:3))
            if (uhf_p) allocate(self%db(ncgto,ncgto,0:3))
        else
            allocate(self%da(ncgto,ncgto,0:0))
            if (uhf_p) allocate(self%db(ncgto,ncgto,0:0))
        end if

    end subroutine

    subroutine read_dens(self)
        type(dens_t), intent(inout) :: self
        character(80) :: xdens_file

        integer(I4) :: b, mo, ispin, s
        type(reorder_t) :: bofh
        real(DP), dimension(:,:), allocatable :: kusse

        if ( .not.associated(self%da) ) then
            call msg_error('read_dens(): dens not allocated!')
            call exit(1)
        end if
        if ( uhf_p.and..not.associated(self%db) ) then
            call msg_error('read_dens(): beta dens not allocated!')
            call exit(1)
        end if

        call getkw(input, 'density', xdens_file)

        open(XDFD, file=xdens_file, status='old', err=42)

        dens=>self%da
        s=size(self%da(:,1,1))
        allocate(kusse(s,s))

        do ispin=1,self%spin
            if (ispin == 2) dens=>self%db
            if (self%pdens_p) then
                do b=0,3
                    read(XDFD,*) kusse
                    dens(:,:,b)=kusse
                end do
            else
                read(XDFD,*) kusse
                dens(:,:,0)=kusse
            end if
        end do

        deallocate(kusse)

        if (uhf_p) then
            call  msg_info('scaling perturbed densities by 0.d5')
            self%da(:,:,1:3)= self%da(:,:,1:3)/2.d0
            self%db(:,:,1:3)= self%db(:,:,1:3)/2.d0
        end if

        if (turbomole_p) then
            call init_reorder(bofh, self%mol)
            call msg_info('Reordering densities [TURBOMOLE]')
            call turbo_reorder(bofh)
            call reorder_dens(bofh, self)
            call del_reorder(bofh)
        end if

        close(XDFD)
        return

42      bert_is_evil=.true. 
        call get_debug_level(b) 
        if (debug_level < 10) then
            call msg_critical('Density file not found!')
            stop
        end if
            
        call msg_error('Density file not found, all densities set to 1')
        call nl
        dens=>self%da
        do ispin=1,self%spin
            if (ispin == 2) dens=>self%db
            do b=0,3
                dens(:,:,b)=1.d0
            end do
        end do
    end subroutine

    subroutine del_dens(self)
        type(dens_t) :: self

        if (associated(self%da)) then
            deallocate(self%da)
        else
            call msg_warn('del_dens(): not allocated!')
        end if
        
        if (uhf_p.and.associated(self%db)) then
            deallocate(self%db)
        end if
        
    end subroutine

    subroutine get_dens(self,a, spin)
        type(dens_t) :: self
        real(DP), dimension(:,:), pointer :: a
        integer(I4), optional :: spin

        dens=>self%da
        if (present(spin)) then
            if (spin == spin_b) dens=>self%db
        end if
        
        a=>dens(:,:,0)
    end subroutine  

    subroutine set_dens(self,a,spin)
        type(dens_t) :: self
        real(DP), dimension(:,:) :: a
        integer(I4), optional :: spin
        
        dens=>self%da
        if (present(spin)) then
            if (spin == spin_b) dens=>self%db
        end if
        dens(:,:,0)=a
    end subroutine  

    subroutine  get_pdens(self,b,a,spin)
        type(dens_t) :: self
        integer(I4), intent(in) :: b
        real(DP), dimension(:,:), pointer :: a
        integer(I4), optional :: spin

        dens=>self%da
        if (present(spin)) then
            if (spin == spin_b) dens=>self%db
        end if
        a=>dens(:,:,b)
    end subroutine  

    subroutine  set_pdens(self,b,a,spin)
        type(dens_t) :: self
        integer(I4), intent(in) :: b
        real(DP), dimension(:,:) :: a
        integer(I4), optional :: spin

        dens=>self%da
        if (present(spin)) then
            if (spin == spin_b) dens=>self%db
        end if
        dens(:,:,b)=a
    end subroutine  

    subroutine reorder_dens(bofh, self)
        type(reorder_t) :: bofh
        type(dens_t) :: self

        integer(I4) :: i,b,p,ncgto,  ispin

        ncgto=get_ncgto(self%mol)
        p=0
        if (self%pdens_p) then
            p=3
        end if

        dens=>self%da
        do ispin=1,self%spin
            if (ispin == 2) dens=>self%db
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

    subroutine moco(self, mos, spin)
        type(dens_t) :: self
        real(DP), dimension(:,:) :: mos
        integer(I4), optional :: spin

        integer(I4), dimension(2) :: moran
        integer(I4) :: ncgto, i, a, b

        dens=>self%da
        if (present(spin)) then
            if (spin == spin_b) then
                if (.not.uhf_p) then
                    call msg_error('gimic: moco(): &
                    &invalid spin for closed shell')
                    stop
                end if
                dens=>self%db
            end if
        end if
        ncgto=get_ncgto(self%mol)
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

    subroutine read_modens(self)
        type(dens_t) :: self
        
        integer(4) :: n, i,j
        real(DP), dimension(:,:), allocatable :: mos
        type(reorder_t) :: bofh
        character(BUFLEN) :: mofile

        if (.not.turbomole_p) then
            call read_dens(self)
            return
        end if

        dens=>self%da
        dens(:,:,0)=D0
        call getkw(input, 'edens.mofile', mofile)
        if (trim(mofile) == '') then
            call read_dens(self)
        end if
        open(XDFD, file=trim(mofile), status='old', err=42)
        read(XDFD, *) 
        read(XDFD, *) 
        read(XDFD, *) 

        n=get_ncgto(self%mol)
        allocate(mos(n,n))
        
        do i=1,n
            read(XDFD, *)
            read(XDFD, '(4f20.2)') mos(:,i)
        end do
        close(XDFD)

        call moco(self, mos)

        call init_reorder(bofh, self%mol)
        call msg_info('Reordering densities [TURBOMOLE]')
        call turbo_reorder(bofh)
        call reorder_dens(bofh, self)
        call del_reorder(bofh)
        
        deallocate(mos)
        return

42      bert_is_evil=.true. 
        if (debug_level < 10) then
            call msg_critical('MO file not found!')
            stop
        end if
        call msg_error('MO file not found, densitiy set to 1.d0')
        dens(:,:,0)=1.d0
    end subroutine
end module 

! vim:et:sw=4:ts=4
