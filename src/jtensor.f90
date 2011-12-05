!
! This is the actual work horse, calculates the current tensor
! for a particular spin case. The cost of calling jtensor twice is 
! very small, since no exponentials have to be calculated.
! 
module jtensor_class
    use globals_m
    use basis_class
    use bfeval_class
    use dens_class
    use dfdb_class
    use dfdr_class
    use d2fdrdb_class
    use dbop_class
    use grid_class
    use teletype_m
    implicit none

!    intrinsic dtime
    
    type jtensor_t
        type(molecule_t), pointer :: mol
        type(dens_t), pointer :: xdens
        type(dbop_t) :: dop
        type(dfdb_t) :: dbt
        type(d2fdrdb_t) :: 	d2f
        type(dfdr_t) :: dfr
        type(bfeval_t) :: bfv
        real(DP), dimension(:), pointer :: pdbf, denbf, dendb
    end type

    public init_jtensor, del_jtensor, jtensor, jtensor2, jvector
    public ctensor, ctensor2
    public qtensor, qtensor2
    public jtensor_t, jdebug, eta
    
    private
    
    real(DP), dimension(3) :: rho

    real(DP), dimension(:), pointer :: bfvec
    real(DP), dimension(:,:), pointer :: dbvec, drvec, d2fvec, dbop
    real(DP), dimension(:,:), pointer :: aodens, pdens
    integer(I4), parameter :: NOTIFICATION=1000

contains
    ! set up memory (once) for the different components
    subroutine init_jtensor(self, mol, xdens)
        type(jtensor_t) :: self
        type(molecule_t), target :: mol
        type(dens_t), target :: xdens
        integer(I4) ::  ncgto
        
        ncgto=get_ncgto(mol)
        
!        call init_dens(ddens)
        self%mol=>mol
        self%xdens=>xdens
        call init_bfeval(self%bfv, self%mol)
        call init_dbop(self%dop, self%mol)
        call init_dfdb(self%dbt, self%mol)
        call init_dfdr(self%dfr, self%mol)
        call init_d2fdrdb(self%d2f, self%mol)

        ! intermediates
        allocate(self%denbf(ncgto))
        allocate(self%dendb(ncgto))
        allocate(self%pdbf(ncgto))
    end subroutine

    subroutine del_jtensor(self)
        type(jtensor_t) :: self

        call del_dbop(self%dop)
        call del_dfdb(self%dbt)
        call del_bfeval(self%bfv)
        call del_dfdr(self%dfr)
        call del_d2fdrdb(self%d2f)
        deallocate(self%denbf)
        deallocate(self%pdbf)
        deallocate(self%dendb)
    end subroutine


    subroutine qtensor(self, r)
        type(jtensor_t) :: self
        real(DP), dimension(3), intent(in) :: r

        type(tensor_t) :: jt1, jt2

        if (uhf_p) then
            call jtensor(self, r, jt1, spin_a)
            call jtensor(self, r, jt2, spin_b)
            write(81, *) jt1%t
            write(82, *) jt2%t
            write(83, *) jt1%t+jt2%t
            write(84, *) jt1%t-jt2%t
        else
            call jtensor(self, r, jt1, spin_a)
            write(85, *) jt1%t
        end if 
    end subroutine

    subroutine qtensor2(self, r)
        type(jtensor_t) :: self
        real(DP), dimension(3), intent(in) :: r

        type(tensor_t) :: pjt1, pjt2
        type(tensor_t) :: djt1, djt2

        if (uhf_p) then
            call jtensor2(self, r, pjt1, djt1, spin_a)
            call jtensor2(self, r, pjt2, djt2, spin_b)
            write(90, *) pjt1%t
            write(91, *) pjt2%t
            write(92, *) djt1%t
            write(93, *) djt2%t
            write(94, *) pjt1%t+pjt2%t
            write(95, *) djt1%t+djt2%t
            write(96, *) pjt1%t-djt1%t
            write(97, *) pjt2%t-djt2%t
        else
            call jtensor2(self, r, pjt1, djt1, spin_a)
            write(98, *) pjt1%t
            write(99, *) djt1%t
        end if 
    end subroutine

    subroutine ctensor(self, r, j, op)
        type(jtensor_t) :: self
        real(DP), dimension(3), intent(in) :: r
        type(tensor_t), intent(inout) :: j
        character(*) :: op

        type(tensor_t) :: jt1, jt2

        select case (op)
            case ('alpha')
                call jtensor(self, r, j, spin_a)
            case ('beta')
                if (uhf_p) then
                    call jtensor(self, r, j, spin_b)
                else
                    call msg_error('ctensor(): &
                    &beta current requested, but not open-shell system!')
                    stop
                end if
            case ('total')
                if (uhf_p) then
                    call jtensor(self, r, jt1, spin_a)
                    call jtensor(self, r, jt2, spin_b)
                    j%t=jt1%t+jt2%t
                else
                    call jtensor(self, r, j, spin_a)
                end if
            case ('spindens')
                if (.not.uhf_p) then
                    call msg_error('ctensor(): &
                    &spindens requested, but not open-shell system!')
                    stop
                end if
                call jtensor(self, r, jt1, spin_a)
                call jtensor(self, r, jt2, spin_b)
                j%t=jt1%t-jt2%t
        end select
    end subroutine

    subroutine ctensor2(self, r, pj, dj, op)
        type(jtensor_t) :: self
        real(DP), dimension(3), intent(in) :: r
        type(tensor_t), intent(inout) :: pj, dj
        character(*) :: op

        type(tensor_t) :: pj1, pj2, dj1, dj2

        select case (op)
            case ('alpha')
                call jtensor2(self, r, pj, dj, spin_a)
            case ('beta')
                if (uhf_p) then
                    call jtensor2(self, r, pj, dj, spin_b)
                else
                    call msg_error('ctensor(): &
                    &beta current requested, but not open-shell system!')
                    stop
                end if
            case ('total')
                if (uhf_p) then
                    call jtensor2(self, r, pj1, dj1, spin_a)
                    call jtensor2(self, r, pj2, dj2, spin_b)
                    pj%t=pj1%t+pj2%t
                    dj%t=dj1%t+dj2%t
                else
                    call jtensor2(self, r, pj, dj, spin_a)
                end if
            case ('spindens')
                if (.not.uhf_p) then
                    call msg_error('ctensor(): &
                    &spindens requested, but not open-shell system!')
                    stop
                end if
                call jtensor2(self, r, pj1, dj1, spin_a)
                call jtensor2(self, r, pj2, dj2, spin_b)
                pj%t=pj1%t-dj1%t
                dj%t=pj2%t-dj2%t
        end select
    end subroutine

    subroutine jtensor(self, r, j, spin)
        type(jtensor_t) :: self
        real(DP), dimension(:), intent(in) :: r
        type(tensor_t), intent(inout) :: j
        integer(I4) :: spin

        integer(I4) :: i, b
        integer(I4), save :: notify=1

        rho=DP50*r ! needed for diamag. contr.

        call bfeval(self%bfv, r, bfvec)
        call mkdbop(self%dop, r, dbop)
        call dfdr(self%dfr, r, drvec)
        if (giao_p) then
            call dfdb(self%dbt, r, bfvec, dbop, dbvec)
            call d2fdrdb(self%d2f, r, bfvec, drvec, dbop, d2fvec)
        end if

        call contract(self, j%t, spin)

!        write(75,*) jt1(:,:,1)
!        write(76,*) jt1(:,:,2)
!        write(77,*) j%t
!        write(78,*) j%t/2.d0


!		if (mod(notify,NOTIFICATION) == 0) then
!			print '(a, i6)', '* points done:', notify
!			call flush(6)
!		end i
!		notify=notify+1
        
    end subroutine

    subroutine jtensor2(self, r, pj, dj, spin)
        type(jtensor_t) :: self
        real(DP), dimension(3), intent(in) :: r
        type(tensor_t), intent(inout) :: pj, dj
        integer(I4) :: spin

        real(DP), dimension(:,:), pointer :: dbop
        integer(I4) :: i, b
        integer(I4), save :: notify=1

        rho=DP50*r ! needed for diamag. contr.

        call bfeval(self%bfv, r, bfvec)
        call dfdr(self%dfr, r, drvec)
        call mkdbop(self%dop, r, dbop)
        if (giao_p) then
            call dfdb(self%dbt, r, bfvec, dbop, dbvec)
            call d2fdrdb(self%d2f, r, bfvec, drvec, dbop, d2fvec)
        end if

        call contract2(self, pj%t, dj%t, spin)

    end subroutine

    subroutine jvector(pj, dj, bb, jv)
        real(DP), dimension(:,:), intent(in) :: pj, dj
        real(DP), dimension(:), intent(in) :: bb
        real(DP), dimension(:), intent(out) :: jv
        
        jv=matmul(pj+dj, bb)
    end subroutine

    subroutine contract2(self, ctp, ctd, spin)
        type(jtensor_t) :: self
        real(DP), dimension(3,3), intent(out) :: ctp, ctd
        integer(I4), intent(in) :: spin

        integer(I4) :: i, j, k, ii,jj
        real(DP) :: prsp1, prsp2       ! paramagnetic wavefunction response
        real(DP) :: ppd                ! paramagnetic probability density
        real(DP), dimension(3) :: dpd  ! diamagnetic probability density
        real(DP) :: diapam

        call get_dens(self%xdens, aodens, spin)  
        self%denbf=matmul(bfvec, aodens)

        k=1
        diapam=dot_product(self%denbf, bfvec)
        do i=1,3! dB <x,y,z>
            ! get perturbed densities: x,y,z
            call get_pdens(self%xdens, i, pdens,spin) 
            self%pdbf=matmul(bfvec, pdens)
            self%dendb=matmul(dbvec(:,i), aodens)
            dpd(i)=diapam*rho(i) ! diamag. contr. to J
            do j=1,3 !dm <x,y,z>
                prsp1=-dot_product(self%dendb, drvec(:,j))    ! (-i)**2=-1 
                prsp2=dot_product(self%denbf, d2fvec(:,k))
                ppd=dot_product(self%pdbf, drvec(:,j))
                ctp(j,i)=ZETA*ppd
                if (giao_p) ctp(j,i)=ctp(j,i)+ZETA*(prsp1+prsp2)
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
    subroutine contract(self, ct, spin)
        type(jtensor_t) :: self
        real(DP), dimension(3,3), intent(out) :: ct
        integer(I4), intent(in) :: spin

        integer(I4) :: b, m, k, ii,jj
        real(DP) :: prsp1, prsp2       ! paramagnetic wavefunction response
        real(DP) :: ppd                ! paramagnetic probability density
        real(DP), dimension(3) :: dpd  ! diamagnetic probability density
        real(DP) :: diapam

        call get_dens(self%xdens, aodens, spin)  
        
        self%denbf=matmul(bfvec, aodens)

        k=1
        diapam=dot_product(self%denbf, bfvec)
        do b=1,3! dB <x,y,z>
            ! get perturbed densities: x,y,z
            call get_pdens(self%xdens, b, pdens, spin) 
            self%pdbf=matmul(bfvec, pdens)
            self%dendb=matmul(dbvec(:,b), aodens)
            dpd(b)=diapam*rho(b) ! diamag. contr. to J
            do m=1,3 !dm <x,y,z>
                prsp1=-dot_product(self%dendb, drvec(:,m))    ! (-i)**2=-1 
                prsp2=dot_product(self%denbf, d2fvec(:,k))
                ppd=dot_product(self%pdbf, drvec(:,m))
                ct(m,b)=ZETA*ppd
                if (giao_p) ct(m,b)=ct(m,b)+ZETA*(prsp1+prsp2)
!                print *, m,b
!                print *, ppd, prsp1, prsp2
                k=k+1
            end do
        end do


! Calculate the total J tensor. The diamagnetic pobability density only
! contributes to the diagonal.

        ! annihilate paramagnetic contribution
        if (.not.paramag_p)	then
            ct=D0
            bert_is_evil=.true.
        end if
        ! annihilate diamagnetic  contribution
        if (.not.diamag_p) then
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

    subroutine eta(self, grid, fac)
        type(jtensor_t) :: self
        type(grid_t) :: grid
        real(DP), intent(in), optional :: fac

        integer(I4) :: i, j, p1, p2, p3
        type(tensor_t) :: foo
        real(DP) :: delta_t
        real(4) :: tim1, tim2
        real(4), dimension(2) :: times
        real(DP), dimension(3) :: bar=(/D1,D1,D1/)
        real(DP), dimension(3) :: foobar
        real(DP), parameter :: SC=0.25d0
        
        call get_grid_size(grid, p1, p2, p3)
        
        call etime(times, tim1)
        tim1=times(1)
        do i=1,100
            call jtensor(self, (/i*SC, i*SC, i*SC/), foo, spin_a)
            foobar=matmul(bar,foo%t)
        end do
        call etime(times, tim2)
        tim2=times(1)
        
        delta_t=tim2-tim1
        if ( present(fac) ) delta_t=delta_t*fac
        write(str_g, '(a,f9.2,a)') 'Estimated CPU time for &
            &calculation: ', delta_t*real(p1*p2*p3)/100.d0, ' sec'
        call msg_info(str_g)
        call nl
    end subroutine	

! Debug function...
!    function contract_killer(nuc, mag, k) result(jj)
!        integer(I4), intent(in) :: nuc, mag, k
!
!        real(DP) :: jj
!        integer(I4) :: mu, nu, n
!        real(DP) :: q
!
!        n=get_ncgto()
!        q=0.0
!        do mu=1,n
!            do nu=1,n 
!                q=q+aodens(nu,mu)*( dbvec(mu,mag)*drvec(nu,nuc)- &
!                  bfvec(mu)*(d2fvec(nu,k)))- &
!                  pdens(nu,mu)*bfvec(mu)*drvec(nu,nuc)
!            end do
!        end do
!        jj=q
!    end function
    
    subroutine jdebug(self)
        type(jtensor_t) :: self

        integer(I4) :: i, b
        integer(I4), save :: notify=1
        real(DP), dimension(3) :: r

        call getkw(input, 'grid.origin', r)

        call bfeval(self%bfv, r, bfvec)
        call dfdr(self%dfr, r, drvec)
        call mkdbop(self%dop, r, dbop)
        call dfdb(self%dbt, r, bfvec, dbop, dbvec)
        call d2fdrdb(self%d2f, r, bfvec, drvec, dbop, d2fvec)
        
        print *, 'bfvec'
        print *, repeat('-', 70)
        print *, bfvec
        print *

        print *, 'drvec'
        print *, repeat('-', 70)
        print 45
        print 41, drvec(1,:)
        print 42, drvec(2,:)
        print 43, drvec(3,:)
        print 44, drvec(4,:)
        print *

        print *, 'dbvec'
        print *, repeat('-', 70)
        print 45
        print 41, dbvec(1,:)
        print 42, dbvec(2,:)
        print 43, dbvec(3,:)
        print 44, dbvec(4,:)
        print *

        print *, 'd2fvec'
        print *, repeat('-', 70)
        print *, 'X'
        print 45
        print 41, d2fvec(1,1:3)
        print 42, d2fvec(2,1:3)
        print 43, d2fvec(3,1:3)
        print 44, d2fvec(4,1:3)
        print *
        print *, 'Y'
        print 45
        print 41, d2fvec(1,4:6)
        print 42, d2fvec(2,4:6)
        print 43, d2fvec(3,4:6)
        print 44, d2fvec(4,4:6)
        print *
        print *, 'Z'
        print 45
        print 41, d2fvec(1,7:9)
        print 42, d2fvec(2,7:9)
        print 43, d2fvec(3,7:9)
        print 44, d2fvec(4,7:9)
41 format('s ',3f15.10)
42 format('px',3f15.10)
43 format('py',3f15.10)
44 format('pz',3f15.10)
45 format('===        x             y               z')
    end subroutine
    
end module
