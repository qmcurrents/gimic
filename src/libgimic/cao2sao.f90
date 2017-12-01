!
! Set up the projection operators for transforming from cartesian to
! spherical angular basis as described in "Molecular Electronic-Structure
! Theory" by Helgaker et al. [Section 9.1.2].
!
! Well... This scheme does not work well with aces2, so the projection
! operators are set up statically by copy 'n paste from the aces2 sources.
! Sigh.
!

module cao2sao_class
    use globals_module
    use settings_module
    use gtodefs_module
    use factorial_module
    use basis_class
    implicit none

    public new_c2sop, del_c2sop, cao2sao
    private

    integer(I4), dimension(0:MAX_L, -MAX_L:MAX_L) :: sphmap
    type(cao2sao_t), dimension(0:MAX_L) :: c2s_oper

    integer(I4) :: am
    real(DP) :: vm

    real(DP), dimension(1,1), target :: sproj
    real(DP), dimension(3,3), target :: pproj
    real(DP), dimension(5,6), target :: dproj
    real(DP), dimension(7,10), target :: fproj
    real(DP), dimension(9,15), target :: gproj
    real(DP), dimension(11,21), target :: hproj

contains
    subroutine setup_c2soper()
        integer(I4) :: l, i,j

        ! Ordering of spherical harmonics in aces2
        ! (the ordering is... random!!?)
        if (is_turbomole) then ! these are NOT right!!!
            sphmap(0,0)=1
            sphmap(1,-1:1)= (/2,3,1/)
            sphmap(2,-2:2)= (/2,5,1,3,4/)
            sphmap(3,-3:3)= (/6,5,2,3,1,7,4/)
            sphmap(4,-4:4)= (/7,5,2,9,1,3,6,8,4/)
            sphmap(5,-5:5)= (/1,1,1,1,1,1,1,1,1,1,1/) !fake
        else
            sphmap(0,0)=1
            sphmap(1,-1:1)= (/2,3,1/)
            sphmap(2,-2:2)= (/2,5,1,3,4/)
            sphmap(3,-3:3)= (/6,5,2,3,1,7,4/)
            sphmap(4,-4:4)= (/7,5,2,9,1,3,6,8,4/)
            sphmap(5,-5:5)= (/1,1,1,1,1,1,1,1,1,1,1/) !fake
        end if

        c2s_oper(0)%po=>sproj
        c2s_oper(1)%po=>pproj
        c2s_oper(2)%po=>dproj
        c2s_oper(3)%po=>fproj
        c2s_oper(4)%po=>gproj
        c2s_oper(5)%po=>hproj
        do l=0,MAX_L
            c2s_oper(l)%po=D0
            call mkc2sop(l, C2S_OPER(l)%po)
        end do
    end subroutine

    subroutine new_c2sop(c2s, mol)
        type(cao2sao_t), target :: c2s
        type(molecule_t) :: mol

        integer(I4) :: nc, ncc, i,j
        integer(I4) :: natoms, nctr, ncomp, nccomp, idx, cidx
        type(atom_t), pointer :: atom
        type(basis_t), pointer :: basis
        type(contraction_t), pointer :: ctr

        call msg_note('Setting up CAO to SAO projection operators')
        call setup_c2soper

        natoms=get_natoms(mol)
        nc=get_ncgto(mol)
        ncc=get_nccgto(mol)
        allocate(c2s%po(nc,ncc))
        c2s%po=D0

        idx=1
        cidx=1
        do i=1,natoms
            call get_atom(mol,i,atom)
            call get_basis(atom, basis)
            nctr=get_nctr(basis)
            do j=1,nctr
                call get_contraction(atom, j, ctr)
                ncomp=get_ncomp(ctr)
                nccomp=get_nccomp(ctr)
                c2s%po(idx:idx+ncomp-1,cidx:cidx+nccomp-1)=c2s_oper(ctr%l)%po
                cidx=cidx+nccomp
                idx=idx+ncomp
            end do
        end do
!        do i=1,nc
!            do j=1,ncc
!                write(*, '(i4)', advance='no') c2s%po(i,j)
!            end do
!            print *
!        end do
    end subroutine

    subroutine del_c2sop(c2s)
        type(cao2sao_t) :: c2s
        deallocate(c2s%po)
    end subroutine

    subroutine cao2sao(c2s, bf, sbf)
        type(cao2sao_t) :: c2s
        real(DP), dimension(:), intent(in) :: bf
        real(DP), dimension(:), intent(out) :: sbf

        sbf=matmul(c2s%po,bf)
    end subroutine

    ! Calculate the normalization factor for SGTO_lm
    function nslm(l,m) result(q)
        integer(I4), intent(in) :: l, m

        real(DP) :: q, fac

        q=1.d0/(2.d0**am)
        q=q/dble(fact(l))
        fac=2.d0*dble(fact(l+am)*fact(l-am))
        if (am == 0) then
            fac=fac/2.d0
        end if
        q=q*sqrt(fac)
    end function

    ! Calculate the normalization factor for SGTO_lm
    function nslm2(l,m) result(q)
        integer(I4), intent(in) :: l, m

        real(DP) :: q, fac

        q=1.d0/(2.d0**am)
        q=q/dble(fact(l))
        fac=2.d0*dble(fact(l+am)*fact(l-am))
        if (am == 0) then
            fac=fac/2.d0
        end if
        q=q*sqrt(fac)
    end function

    ! Calculate the expansion coefficient for SGTO_lm
    function clmtuv(l,m,t,u,v) result(q)
        integer(I4), intent(in) :: l, m, t, u
        real(DP) :: v
        integer(I4) :: iq

        real(DP) :: q

        q=4.d0**t
        q=1.d0/q
        if (mod(dble(t)+v-vm, 2.d0) > 0.d0) then
            q=-q
        end if
        iq=int(binom(l,t)*binom(l-t,am+t)*binom(t,u)*binom(am,int(2.d0*v)))
        q=q*dble(iq)
    end function

    ! Construct the projection operator for basis functions with
    ! quantum number l.
    ! This routine has a bug...
    subroutine mkc2sop(l, xop)
        integer(I4), intent(in) :: l
        real(DP), dimension(:,:), intent(out) :: xop

        integer(I4) :: t, u, v
        integer(I4) :: i, j, k, m, n, idx
        real(DP) :: vr, xnslm

        xop=0.0
        idx=0
        do m=-l,l
!            idx=sphmap(l,m)
            idx=idx+1
            am=abs(m)
            vm=0.0
            if ( m < 0 ) vm=0.5d0
            xnslm=nslm(l,m)
            do t=0,(l-am)/2
                do u=0,t
                    do v=0,int(am*0.5d0-vm)
                        vr=dble(v)+vm
                        i=2*t+am-int(2.d0*(dble(u)+vr))
                        j=int(2.d0*(dble(u)+vr))
                        k=l-2*t-am
                        n=gtomap(i,j,k)
                        xop(idx,n)=clmtuv(l,m,t,u,vr)*xnslm
                    end do
                end do
            end do
        end do
        call renorm(xop)

    end subroutine

    ! Renormalize expansion coefficients to integers.
    ! First find the smallest expansion coefficient, and divide.
    ! The trick is then to find the smallest integer multiplicator that
    ! makes all expansion coefficients integer.
    subroutine renorm(xop)
        real(DP), dimension(:,:), intent(inout) :: xop

        integer(I4) :: i, j, k
        real(DP) :: amin, ava, avb, sim

        do i=1,size(xop,1)
            amin=1.d+10
            do j=1,size(xop,2)
                ava=abs(xop(i,j))
                if (ava > 0.d0) then
                    if (ava < amin) amin=ava
                end if
            end do
            xop(i,:)=xop(i,:)/amin
            sim=D1
            do j=1,size(xop,2)
                ava=abs(xop(i,j))
                if (ava > D0) then
                    avb=mod(ava,D1)
                    if (avb > 1.d-10) then
                        ava=D1/avb
                        do while(mod(ava,D1) > 1.d-10)
                            avb=mod(ava, D1)
                            if (avb > 1.d-10) ava=ava/avb
                        end do
                        if (ava > sim) sim=ava
                    end if
                end if
            end do
            xop(i,:)=xop(i,:)*sim
        end do

!        call print_xop(xop)
    end subroutine

    subroutine print_xop(xop)
        real(DP), dimension(:,:), intent(inout) :: xop

        integer(I4) :: i, j

        do i=1,size(xop,1)
            do j=1,size(xop,2)
                write(*, '(i4)', advance='no') xop(i,j)
            end do
            print *
        end do
        print *
    end subroutine
end module

! vim:et:sw=4:ts=4
