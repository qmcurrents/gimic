module reorder_class
    use globals_m
    use settings_m
    use basis_class
    implicit none

    type reorder_t
        type(molecule_t), pointer :: mol
        integer(I4), dimension(:), pointer :: sv
        real(DP), dimension(:,:), pointer :: tmp
    end type

contains
    subroutine reorder_vec(bofh, vec)
        type(reorder_t), intent(inout) :: bofh
        real(DP), dimension(:), intent(inout) :: vec

        integer(I4) :: i, ncgto
        
        ncgto=get_ncgto(bofh%mol)
        do i=1,ncgto
            bofh%tmp(bofh%sv(i),1)=vec(i)
        end do
        vec=bofh%tmp(:,1)
    end subroutine

    subroutine reorder_cols(bofh, vec)
        type(reorder_t), intent(inout) :: bofh
        real(DP), dimension(:,:), intent(inout) :: vec

        integer(I4) :: i, ncgto
        
        ncgto=get_ncgto(bofh%mol)
        do i=1,ncgto
            bofh%tmp(:,bofh%sv(i))=vec(:,i)
        end do
        vec=bofh%tmp
    end subroutine

    subroutine new_reorder(bofh, mol)
        type(reorder_t) :: bofh
        type(molecule_t), target :: mol

        integer(I4) :: ncgto

        bofh%mol=>mol

        ncgto=get_ncgto(bofh%mol)

        allocate(bofh%sv(ncgto))
        allocate(bofh%tmp(ncgto,ncgto))
    end subroutine

    subroutine turbo_reorder(bofh)
        type(reorder_t) :: bofh

        integer(I4) :: i, j, k, l, m, c, q, natoms, ncgto
        type(atom_t),  pointer :: atom
        type(basis_t), pointer :: aob
        type(contraction_t), pointer :: ctr
        integer(I4), dimension(:), allocatable :: lvec
        
        natoms=get_natoms(bofh%mol)
        ncgto=get_ncgto(bofh%mol)

        allocate(lvec(ncgto))
        
        q=1
        do i=1,natoms
            call get_atom(bofh%mol, i, atom)
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
                    bofh%sv(q)=i
                    q=q+1
                end if
            end do
            l=l+1
        end do
        deallocate(lvec)
    end subroutine

    subroutine degeneralize(bofh)
        type(reorder_t) :: bofh

        integer(I4) :: i, j, k, l, m, base, q, natoms, ncgto
        type(atom_t),  pointer :: atom
        type(basis_t), pointer :: aob
        type(contraction_t), pointer :: ctr
        
        natoms=get_natoms(bofh%mol)
        ncgto=get_ncgto(bofh%mol)
        
        bofh%sv=0
        q=1
        do i=1,natoms
            call get_atom(bofh%mol, i, atom)
            call get_basis(atom,aob)
            j=1
            do while (j <= aob%nctr)
                base=q
                ctr=>aob%ctr(j)
                if (ctr%l == 0) then
                    bofh%sv(q)=q
!                    print *, q, base, bofh%sv(q)
                    q=q+1
                    j=j+1
                    cycle
                end if
                do k=1,ctr%ncf
                    do m=0,ctr%ncomp-1
                        bofh%sv(q)=base+ctr%ncf*m
!                        print *, q, base, bofh%sv(q)
                        q=q+1
                    end do
                    base=base+1
                end do
                j=j+ctr%ncf
            end do
        end do
    end subroutine

    subroutine del_reorder(bofh)
        type(reorder_t) :: bofh

        deallocate(bofh%sv)
        deallocate(bofh%tmp)
    end subroutine 

!    subroutine ugly_reorder_hack(bofh, vec)
!        type(reorder_t), intent(inout) :: bofh
!        real(DP), dimension(:), intent(inout) :: vec
!
!        integer(I4) :: curl, ncgto, pos, i
!        
!        ncgto=get_ncgto()
!        bofh%tmp=0.d0
!        curl=0 
!        pos=1
!        do while (pos <= ncgto)
!            do i=1,ncgto
!                if (bofh%vec(i) == curl) then
!                    bofh%tmp(pos)=vec(i)
!                    bofh%tmp(i)=vec(pos)
!                    pos=pos+1
!                end if
!            end do
!            curl=curl+1
!        end do
!        vec=bofh%tmp
!    end subroutine


end module


! vim:et:sw=4:ts=4
