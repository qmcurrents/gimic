!
! This is the actual work horse, calculates the current tensor
! for a particular spin case. The cost of calling jtensor twice is
! very small, since no exponentials have to be calculated.
!
module jtensor_class
    use globals_module
    use settings_module
    use basis_class
    use bfeval_class
    use dens_class
    use teletype_module
    implicit none

!    intrinsic dtime

    type jtensor_t
        private
        type(molecule_t), pointer :: mol
        type(dens_t), pointer :: xdens
        type(bfeval_t) :: basv
        real(DP), dimension(:), pointer :: pdbf, denbf, dendb
        real(DP), dimension(3) :: rho
        ! scratch mem
        real(DP), dimension(:), pointer :: bfvec
        real(DP), dimension(:,:), pointer :: dbvec, drvec, d2fvec, dbop
        real(DP), dimension(:,:), pointer :: aodens, pdens
    end type

    public new_jtensor, del_jtensor, jtensor, jtensor2, get_jvector
    public ctensor, ctensor2, jvector
    public jtensor_t

    private
    integer(I4), parameter :: NOTIFICATION=1000

contains
    ! set up memory (once) for the different components
    subroutine new_jtensor(this, mol, xdens)
        type(jtensor_t) :: this
        type(molecule_t), target :: mol
        type(dens_t), target :: xdens
        integer(I4) ::  ncgto

        ncgto=get_ncgto(mol)

!        call new_dens(ddens)
        this%mol=>mol
        this%xdens=>xdens
        call new_bfeval(this%basv, this%mol)

        ! intermediates
        allocate(this%denbf(ncgto))
        allocate(this%dendb(ncgto))
        allocate(this%pdbf(ncgto))
    end subroutine

    subroutine del_jtensor(this)
        type(jtensor_t) :: this

        call del_bfeval(this%basv)
        deallocate(this%denbf)
        deallocate(this%pdbf)
        deallocate(this%dendb)
    end subroutine

    subroutine ctensor(this, r, j, op)
        type(jtensor_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(9), intent(inout) :: j
        character(*) :: op

        real(DP), dimension(9) :: jt1, jt2

        select case (op)
            case ('alpha')
                call jtensor(this, r, j, spin_a)
            case ('beta')
                if (settings%is_uhf) then
                    call jtensor(this, r, j, spin_b)
                else
                    call msg_error('ctensor(): &
                    &beta current requested, but not open-shell system!')
                    stop
                end if
            case ('total')
                if (settings%is_uhf) then
                    call jtensor(this, r, jt1, spin_a)
                    call jtensor(this, r, jt2, spin_b)
                    j=jt1+jt2
                else
                    call jtensor(this, r, j, spin_a)
                end if
            case ('spindens')
                if (.not.settings%is_uhf) then
                    call msg_error('ctensor(): &
                    &spindens requested, but not open-shell system!')
                    stop
                end if
                call jtensor(this, r, jt1, spin_a)
                call jtensor(this, r, jt2, spin_b)
                j=jt1-jt2
        end select
    end subroutine

    subroutine ctensor2(this, r, pj, dj, op)
        type(jtensor_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(9), intent(inout) :: pj, dj
        character(*) :: op

        real(DP), dimension(9) :: pj1, pj2, dj1, dj2

        select case (op)
            case ('alpha')
                call jtensor2(this, r, pj, dj, spin_a)
            case ('beta')
                if (settings%is_uhf) then
                    call jtensor2(this, r, pj, dj, spin_b)
                else
                    call msg_error('ctensor(): &
                    &beta current requested, but not open-shell system!')
                    stop
                end if
            case ('total')
                if (settings%is_uhf) then
                    call jtensor2(this, r, pj1, dj1, spin_a)
                    call jtensor2(this, r, pj2, dj2, spin_b)
                    pj=pj1+pj2
                    dj=dj1+dj2
                else
                    call jtensor2(this, r, pj, dj, spin_a)
                end if
            case ('spindens')
                if (.not.settings%is_uhf) then
                    call msg_error('ctensor(): &
                    &spindens requested, but not open-shell system!')
                    stop
                end if
                call jtensor2(this, r, pj1, dj1, spin_a)
                call jtensor2(this, r, pj2, dj2, spin_b)
                pj=pj1-dj1
                dj=pj2-dj2
        end select
    end subroutine

    subroutine jtensor(this, r, j, spin)
        type(jtensor_t) :: this
        real(DP), dimension(:), intent(in) :: r
        real(DP), dimension(9), intent(inout) :: j
        integer(I4) :: spin

        integer(I4) :: i, b
        integer(I4), save :: notify=1

        this%rho=DP50*r ! needed for diamag. contr.
        if (settings%use_giao) then
            call calc_basis(this%basv, r, this%bfvec, this%drvec, &
            this%dbvec, this%d2fvec)
        else
            call calc_basis(this%basv, r, this%bfvec, this%drvec)
        end if

        call contract(this, j, spin)
    end subroutine

    subroutine jtensor2(this, r, pj, dj, spin)
        type(jtensor_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(9), intent(inout) :: pj, dj
        integer(I4) :: spin

!        real(DP), dimension(:,:), pointer :: dbop
        integer(I4) :: i, b
        integer(I4), save :: notify=1

        this%rho=DP50*r ! needed for diamag. contr.
        if (settings%use_giao) then
            call calc_basis(this%basv, r, this%bfvec, this%drvec, &
            this%dbvec, this%d2fvec)
        else
            call calc_basis(this%basv, r, this%bfvec, this%drvec)
        end if

        call contract2(this, pj, dj, spin)

    end subroutine

    subroutine jvector(this, r, bb, jv, op)
        type(jtensor_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(:), intent(in) :: bb
        real(DP), dimension(:), intent(out) :: jv
        character(*) :: op

        real(DP), dimension(9) :: j
        call ctensor(this, r, j, op)
        jv=matmul(reshape(j,(/3,3/)), bb)
    end subroutine

    subroutine get_jvector(pj, dj, bb, jv)
        real(DP), dimension(:), intent(in) :: pj, dj
        real(DP), dimension(:), intent(in) :: bb
        real(DP), dimension(:), intent(out) :: jv

        jv=matmul(reshape(pj+dj,(/3,3/)), bb)
    end subroutine

    subroutine contract2(this, ctp, ctd, spin)
        type(jtensor_t) :: this
        real(DP), dimension(3,3), intent(out) :: ctp, ctd
        integer(I4), intent(in) :: spin

        integer(I4) :: i, j, k, ii,jj
        real(DP) :: prsp1, prsp2       ! paramagnetic wavefunction response
        real(DP) :: ppd                ! paramagnetic probability density
        real(DP), dimension(3) :: dpd  ! diamagnetic probability density
        real(DP) :: diapam

        call get_dens(this%xdens, this%aodens, spin)
        this%denbf=matmul(this%bfvec, this%aodens)

        k=1
        diapam=dot_product(this%denbf, this%bfvec)
        do i=1,3! dB <x,y,z>
            ! get perturbed densities: x,y,z
            call get_pdens(this%xdens, i, this%pdens,spin)
            this%pdbf=matmul(this%bfvec, this%pdens)
            this%dendb=matmul(this%dbvec(:,i), this%aodens)
            dpd(i)=diapam*this%rho(i) ! diamag. contr. to J
            do j=1,3 !dm <x,y,z>
                prsp1=-dot_product(this%dendb, this%drvec(:,j))    ! (-i)**2=-1
                prsp2=dot_product(this%denbf, this%d2fvec(:,k))
                ppd=dot_product(this%pdbf, this%drvec(:,j))
                ctp(j,i)=ZETA*ppd
                if (settings%use_giao) ctp(j,i)=ctp(j,i)+ZETA*(prsp1+prsp2)
                k=k+1
            end do
        end do

! Calculate the total J tensor. The diamagnetic pobability density only
! contributes to the diagonal.

        ctd(1,1)=D0; ctd(2,2)=D0; ctd(3,3)=D0;
        ctd(1,2)=dpd(3)
        ctd(1,3)=-1.d0*dpd(2)
        ctd(2,1)=-1.d0*dpd(3)
        ctd(2,3)=dpd(1)
        ctd(3,1)=dpd(2)
        ctd(3,2)=-1.d0*dpd(1)

    end subroutine

!
! Contract all contributions to J. This is where all the actual work is done.
!
    subroutine contract(this, ct, spin)
        type(jtensor_t) :: this
        real(DP), dimension(3,3), intent(out) :: ct
        integer(I4), intent(in) :: spin

        integer(I4) :: b, m, k, ii,jj
        real(DP) :: prsp1, prsp2       ! paramagnetic wavefunction response
        real(DP) :: ppd                ! paramagnetic probability density
        real(DP), dimension(3) :: dpd  ! diamagnetic probability density
        real(DP) :: diapam

        call get_dens(this%xdens, this%aodens, spin)
        this%denbf=matmul(this%bfvec, this%aodens)

        k=1
        diapam=dot_product(this%denbf, this%bfvec)
        do b=1,3! dB <x,y,z>
            ! get perturbed densities: x,y,z
            call get_pdens(this%xdens, b, this%pdens, spin)
            this%pdbf=matmul(this%bfvec, this%pdens)
            if (settings%use_giao) then
              this%dendb=matmul(this%dbvec(:,b), this%aodens)
            end if
            dpd(b)=diapam*this%rho(b) ! diamag. contr. to J
            do m=1,3 !dm <x,y,z>
              if (settings%use_giao) then
                prsp1=-dot_product(this%dendb, this%drvec(:,m)) ! (-i)**2=-1
                prsp2=dot_product(this%denbf, this%d2fvec(:,k))
              end if
              ppd=dot_product(this%pdbf, this%drvec(:,m))
              ct(m,b)=ZETA*ppd
              if (settings%use_giao) ct(m,b)=ct(m,b)+ZETA*(prsp1+prsp2)
!              print *, m,b
!              print *, ppd, prsp1, prsp2
              k=k+1
            end do
        end do


! Calculate the total J tensor. The diamagnetic pobability density only
! contributes to the diagonal.

        ! annihilate paramagnetic contribution
        if (.not.settings%use_paramag) then
            ct=D0
            bert_is_evil=.true.
        end if
        ! annihilate diamagnetic  contribution
        if (.not.settings%use_diamag) then
            dpd=D0
            bert_is_evil=.true.
        end if

        ct(1,2)=ct(1,2)+dpd(3)
        ct(1,3)=ct(1,3)-dpd(2)
        ct(2,1)=ct(2,1)-dpd(3)
        ct(2,3)=ct(2,3)+dpd(1)
        ct(3,1)=ct(3,1)+dpd(2)
        ct(3,2)=ct(3,2)-dpd(1)

    end subroutine

!    subroutine jdebug(this, r)
!        type(jtensor_t) :: this
!        real(DP), dimension(3), intent(in) :: r

!        integer(I4) :: i, b
!        integer(I4), save :: notify=1

!        call bfeval(this%bfv, r, this%bfvec)
!        call dfdr(this%dfr, r, this%drvec)
!        call mkdbop(this%dop, r, this%dbop)
!        call dfdb(this%dbt, r, this%bfvec, this%dbop, this%dbvec)
!        call d2fdrdb(this%d2f, r, this%bfvec, this%drvec, this%dbop, this%d2fvec)
!
!        print *, 'bfvec'
!        print *, repeat('-', 70)
!        print *, this%bfvec
!        print *

!        print *, 'drvec'
!        print *, repeat('-', 70)
!        print 45
!        print 41, this%drvec(1,:)
!        print 42, this%drvec(2,:)
!        print 43, this%drvec(3,:)
!        print 44, this%drvec(4,:)
!        print *

!        print *, 'dbvec'
!        print *, repeat('-', 70)
!        print 45
!        print 41, this%dbvec(1,:)
!        print 42, this%dbvec(2,:)
!        print 43, this%dbvec(3,:)
!        print 44, this%dbvec(4,:)
!        print *

!        print *, 'd2fvec'
!        print *, repeat('-', 70)
!        print *, 'X'
!        print 45
!        print 41, this%d2fvec(1,1:3)
!        print 42, this%d2fvec(2,1:3)
!        print 43, this%d2fvec(3,1:3)
!        print 44, this%d2fvec(4,1:3)
!        print *
!        print *, 'Y'
!        print 45
!        print 41, this%d2fvec(1,4:6)
!        print 42, this%d2fvec(2,4:6)
!        print 43, this%d2fvec(3,4:6)
!        print 44, this%d2fvec(4,4:6)
!        print *
!        print *, 'Z'
!        print 45
!        print 41, this%d2fvec(1,7:9)
!        print 42, this%d2fvec(2,7:9)
!        print 43, this%d2fvec(3,7:9)
!        print 44, this%d2fvec(4,7:9)
!41 format('s ',3f15.10)
!42 format('px',3f15.10)
!43 format('py',3f15.10)
!44 format('pz',3f15.10)
!45 format('===        x             y               z')
!    end subroutine

end module

! vim:et:sw=4:ts=4
