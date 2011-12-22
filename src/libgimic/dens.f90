
module dens_class
    use globals_m
    use settings_m
    use basis_class
    use reorder_class
    implicit none

    public new_dens, read_dens, set_dens, set_pdens
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
    subroutine new_dens(this, mol, modens_p)
        type(dens_t) :: this
        type(molecule_t), target :: mol
        logical, optional :: modens_p

        integer(I4) :: ncgto

        this%mol=>mol
        ncgto=get_ncgto(mol)

        if (present(modens_p)) then
            if (modens_p) then
                this%pdens_p=.false.
            end if
        end if

        if (settings%is_uhf) then
            this%spin=2
        else
            this%spin=1
        end if

        if (this%pdens_p) then
            allocate(this%da(ncgto,ncgto,0:3))
            if (settings%is_uhf) allocate(this%db(ncgto,ncgto,0:3))
        else
            allocate(this%da(ncgto,ncgto,0:0))
            if (settings%is_uhf) allocate(this%db(ncgto,ncgto,0:0))
        end if

    end subroutine

    subroutine read_dens(this, xdens_file)
        type(dens_t), intent(inout) :: this
        character(80), intent(in) :: xdens_file

        integer(I4) :: b, mo, ispin, s
        type(reorder_t) :: bofh
        real(DP), dimension(:,:), allocatable :: kusse

        if ( .not.associated(this%da) ) then
            call msg_error('read_dens(): dens not allocated!')
            call exit(1)
        end if
        if ( settings%is_uhf.and..not.associated(this%db) ) then
            call msg_error('read_dens(): beta dens not allocated!')
            call exit(1)
        end if

        open(XDFD, file=xdens_file, status='old', err=42)

        dens=>this%da
        s=size(this%da(:,1,1))
        allocate(kusse(s,s))

        do ispin=1,this%spin
            if (ispin == 2) dens=>this%db
            if (this%pdens_p) then
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

        if (settings%is_uhf) then
            call  msg_info('scaling perturbed densities by 0.d5')
            this%da(:,:,1:3)= this%da(:,:,1:3)/2.d0
            this%db(:,:,1:3)= this%db(:,:,1:3)/2.d0
        end if

        if (is_turbomole) then
            call new_reorder(bofh, this%mol)
            call msg_info('Reordering densities [TURBOMOLE]')
            call turbo_reorder(bofh)
            call reorder_dens(bofh, this)
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
        dens=>this%da
        do ispin=1,this%spin
            if (ispin == 2) dens=>this%db
            do b=0,3
                dens(:,:,b)=1.d0
            end do
        end do
    end subroutine

    subroutine del_dens(this)
        type(dens_t) :: this

        if (associated(this%da)) then
            deallocate(this%da)
        else
            call msg_warn('del_dens(): not allocated!')
        end if
        
        if (settings%is_uhf.and.associated(this%db)) then
            deallocate(this%db)
        end if
        
    end subroutine

    subroutine get_dens(this,a, spin)
        type(dens_t) :: this
        real(DP), dimension(:,:), pointer :: a
        integer(I4), optional :: spin

        if (present(spin)) then
            if (spin == spin_b) then
                a=>this%db(:,:,0)
            else 
                a=>this%da(:,:,0)
            end if
        else
            a=>this%da(:,:,0)
        end if
    end subroutine  

    subroutine set_dens(this,a,spin)
        type(dens_t) :: this
        real(DP), dimension(:,:) :: a
        integer(I4), optional :: spin
        
        dens=>this%da
        if (present(spin)) then
            if (spin == spin_b) dens=>this%db
        end if
        dens(:,:,0)=a
    end subroutine  

    subroutine  get_pdens(this,b,a,spin)
        type(dens_t) :: this
        integer(I4), intent(in) :: b
        real(DP), dimension(:,:), pointer :: a
        integer(I4), optional :: spin

        if (present(spin)) then
            if (spin == spin_b) then
                a=>this%db(:,:,b)
            else 
                a=>this%da(:,:,b)
            end if
        else
            a=>this%da(:,:,b)
        end if
    end subroutine  

    subroutine  set_pdens(this,b,a,spin)
        type(dens_t) :: this
        integer(I4), intent(in) :: b
        real(DP), dimension(:,:) :: a
        integer(I4), optional :: spin

        dens=>this%da
        if (present(spin)) then
            if (spin == spin_b) dens=>this%db
        end if
        dens(:,:,b)=a
    end subroutine  

    subroutine reorder_dens(bofh, this)
        type(reorder_t) :: bofh
        type(dens_t) :: this

        integer(I4) :: i,b,p,ncgto,  ispin

        ncgto=get_ncgto(this%mol)
        p=0
        if (this%pdens_p) then
            p=3
        end if

        dens=>this%da
        do ispin=1,this%spin
            if (ispin == 2) dens=>this%db
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

    subroutine moco(this, mos, morange, spin)
        type(dens_t) :: this
        real(DP), dimension(:,:) :: mos
        integer(I4), dimension(2), optional :: morange
        integer(I4), optional :: spin

        integer(I4), dimension(2) :: moran
        integer(I4) :: ncgto, i, a, b

        dens=>this%da
        if (present(spin)) then
            if (spin == spin_b) then
                if (.not.settings%is_uhf) then
                    call msg_error('gimic: moco(): &
                    &invalid spin for closed shell')
                    stop
                end if
                dens=>this%db
            end if
        end if
        ncgto=get_ncgto(this%mol)
        moran=0
        if (present(morange)) then
            moran = morange
        end if

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

    subroutine read_modens(this, dens_file, mofile, morange)
        type(dens_t) :: this
        character(*), intent(in) :: dens_file
        character(*), intent(in) :: mofile
        integer(I4), dimension(2), optional :: morange
        
        integer(4) :: n, i,j
        real(DP), dimension(:,:), allocatable :: mos
        type(reorder_t) :: bofh

        if (.not.is_turbomole) then
            call read_dens(this, dens_file)
            return
        end if

        dens=>this%da
        dens(:,:,0)=D0
        if (trim(mofile) == '') then
            call read_dens(this, dens_file)
        end if
        open(XDFD, file=trim(mofile), status='old', err=42)
        read(XDFD, *) 
        read(XDFD, *) 
        read(XDFD, *) 

        n=get_ncgto(this%mol)
        allocate(mos(n,n))
        
        do i=1,n
            read(XDFD, *)
            read(XDFD, '(4f20.2)') mos(:,i)
        end do
        close(XDFD)

        if (present(morange)) then
            call moco(this, mos, morange)
        else
            call moco(this, mos)
        end if

        call new_reorder(bofh, this%mol)
        call msg_info('Reordering densities [TURBOMOLE]')
        call turbo_reorder(bofh)
        call reorder_dens(bofh, this)
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
