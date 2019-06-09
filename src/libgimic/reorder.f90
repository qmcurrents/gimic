module reorder_class
    use globals_module
    use settings_module
    use basis_class
    implicit none

    type reorder_t
        type(molecule_t), pointer :: mol
        integer(I4), dimension(:), pointer :: sv
        real(DP), dimension(:,:), pointer :: tmp
    end type

contains
    subroutine reorder_vec(this, vec)
        type(reorder_t), intent(inout) :: this
        real(DP), dimension(:), intent(inout) :: vec

        integer(I4) :: i, ncgto

        ncgto=get_ncgto(this%mol)
        do i=1,ncgto
            this%tmp(this%sv(i),1)=vec(i)
        end do
        vec=this%tmp(:,1)
    end subroutine

    subroutine reorder_cols(this, vec)
        type(reorder_t), intent(inout) :: this
        real(DP), dimension(:,:), intent(inout) :: vec

        integer(I4) :: i, ncgto

        ncgto=get_ncgto(this%mol)
        do i=1,ncgto
            this%tmp(:,this%sv(i))=vec(:,i)
        end do
        vec=this%tmp
    end subroutine

    subroutine new_reorder(this, mol)
        type(reorder_t) :: this
        type(molecule_t), target :: mol

        integer(I4) :: ncgto

        this%mol=>mol

        ncgto=get_ncgto(this%mol)

        allocate(this%sv(ncgto))
        allocate(this%tmp(ncgto,ncgto))
    end subroutine

    subroutine turbo_reorder(this)
        type(reorder_t) :: this

        integer(I4) :: i, j, k, l, m, c, q, natoms, ncgto
        type(atom_t),  pointer :: atom
        type(basis_t), pointer :: aob
        type(contraction_t), pointer :: ctr
        integer(I4), dimension(:), allocatable :: lvec

        natoms=get_natoms(this%mol)
        ncgto=get_ncgto(this%mol)

        allocate(lvec(ncgto))

        q=1
        do i=1,natoms
            call get_atom(this%mol, i, atom)
            call get_basis(atom,aob)
            k=1
            do l=1,aob%nshells
                do j=1,aob%nctrps(l)
                    ctr=>aob%ctr(k)
                    k=k+1
                    do c=1,ctr%ncomp
                        lvec(q)=ctr%l
                        q=q+1
                    end do
                end do
            end do
        end do
        l=0
        q=1
        do while (q <= ncgto)
            do i=1,ncgto
                if (lvec(i) == l) then
                    this%sv(q)=i
                    q=q+1
                end if
            end do
            l=l+1
        end do
        deallocate(lvec)
    end subroutine

    subroutine degeneralize(this)
        type(reorder_t) :: this

        integer(I4) :: i, j, k, l, m, base, q, natoms, ncgto
        type(atom_t),  pointer :: atom
        type(basis_t), pointer :: aob
        type(contraction_t), pointer :: ctr

        natoms=get_natoms(this%mol)
        ncgto=get_ncgto(this%mol)

        this%sv=0
        q=1
        do i=1,natoms
            call get_atom(this%mol, i, atom)
            call get_basis(atom,aob)
            j=1
            do while (j <= aob%nctr)
                base=q
                ctr=>aob%ctr(j)
                if (ctr%l == 0) then
                    this%sv(q)=q
!                    print *, q, base, this%sv(q)
                    q=q+1
                    j=j+1
                    cycle
                end if
                do k=1,ctr%ncf
                    do m=0,ctr%ncomp-1
                        this%sv(q)=base+ctr%ncf*m
!                        print *, q, base, this%sv(q)
                        q=q+1
                    end do
                    base=base+1
                end do
                j=j+ctr%ncf
            end do
        end do
    end subroutine

    subroutine del_reorder(this)
        type(reorder_t) :: this

        deallocate(this%sv)
        deallocate(this%tmp)
    end subroutine

!    subroutine ugly_reorder_hack(this, vec)
!        type(reorder_t), intent(inout) :: this
!        real(DP), dimension(:), intent(inout) :: vec
!
!        integer(I4) :: curl, ncgto, pos, i
!
!        ncgto=get_ncgto()
!        this%tmp=0.d0
!        curl=0
!        pos=1
!        do while (pos <= ncgto)
!            do i=1,ncgto
!                if (this%vec(i) == curl) then
!                    this%tmp(pos)=vec(i)
!                    this%tmp(i)=vec(pos)
!                    pos=pos+1
!                end if
!            end do
!            curl=curl+1
!        end do
!        vec=this%tmp
!    end subroutine


end module


! vim:et:sw=4:ts=4
