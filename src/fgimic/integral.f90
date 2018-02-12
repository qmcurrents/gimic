! Integrate the current density using quadrature.
!
! Written by Jonas Juselius, University of Helsinki, 2003 (I think).
! Improvements by jj, University of Tromsø, 2011 (I'm sure).
!

module integral_class
    use globals_module
    use settings_module
    use grid_class
    use jfield_class
    use dens_class
    use jtensor_class
    use gaussint_module
    use lip_module
    use teletype_module
    use parallel_module
    use magnet_module
    ! ACID stuff
    use acid_module
    use tensor_module
    implicit none

    type integral_t
        type(grid_t), pointer :: grid
    end type

    public new_integral, del_integral, integral_t
    public integrate_tensor_field, integrate_current
    public integrate_modulus
    public integrate_acid

    private

    character(8) :: spin='total'
    integer(I4) :: nlip
contains
    subroutine new_integral(this, grid)
        type(integral_t) :: this
        type(grid_t), target :: grid

        this%grid=>grid
    end subroutine

    subroutine del_integral(this)
        type(integral_t) :: this
    end subroutine

    subroutine integrate_current(this, mol, xdens, spinn)
        type(integral_t), intent(inout) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens
        character(*), optional :: spinn

        integer(I4) :: i, j, k, p1, p2, p3, lo, hi
        real(DP), dimension(3) :: normal, rr, center, bb, bb_prime, bb_second
        real(DP) :: psum, nsum, w, jp, r, bound
        real(DP) :: psum2, nsum2
        real(DP) :: psum3, nsum3
        real(DP) :: xsum, xsum2, xsum3

        real(DP) :: jp_prime, jp_second
        real(DP) :: p_prime_sum, n_prime_sum
        real(DP) :: p_prime_sum2, n_prime_sum2
        real(DP) :: p_prime_sum3, n_prime_sum3
        real(DP) :: x_prime_sum, x_prime_sum2, x_prime_sum3
        real(DP) :: p_second_sum, n_second_sum
        real(DP) :: p_second_sum2, n_second_sum2
        real(DP) :: p_second_sum3, n_second_sum3
        real(DP) :: x_second_sum, x_second_sum2, x_second_sum3

        real(DP), dimension(3) :: jvec, jvec_prime, jvec_second
        real(DP), dimension(9) :: tt
        type(jtensor_t) :: jt

        if (present(spinn)) spin = spinn

        if (settings%is_uhf) then
            select case(spin)
                case('total')
                    call msg_note("Integrating total density")
                case('alpha')
                    call msg_note("Integrating alpha density")
                case('beta')
                    call msg_note("Integrating beta density")
                case('spindens')
                    call msg_note("Integrating spin density")
                case default
                    call msg_error("Invalid spin: " // spin)
                    stop
            end select
        end if

        call get_grid_size(this%grid, p1, p2, p3)
        call get_magnet(this%grid, bb)
        !call jfield_eta(this%jf)

        normal=get_grid_normal(this%grid)

        bb_prime=normal
        bb_second=cross_product(bb,bb_prime);

        bound=1.d+10
        bound=this%grid%radius
        if (bound < 1.d+10) then
            write(str_g, *) 'Integration bound set to radius ', bound
            call msg_out(str_g)
        end if

        call grid_center(this%grid,center)
        call schedule(p2, lo, hi)
        if (mpi_rank == 0) then
            lo = 1
        end if

        xsum3=0.d0
        psum3=0.d0
        nsum3=0.d0

        x_prime_sum3=0.d0
        p_prime_sum3=0.d0
        n_prime_sum3=0.d0
        x_second_sum3=0.d0
        p_second_sum3=0.d0
        n_second_sum3=0.d0
!!$OMP PARALLEL PRIVATE(i,j,k,r,rr,xsum,psum,nsum) &
!!$OMP PRIVATE(jt,w,jp,tt,jvec) &
!!$OMP SHARED(xsum2,psum2,nsum2) &
!!$OMP SHARED(p1,p2,p3,this,center,spin,bb,normal,mol,xdens,lo,hi) &
!!$OMP REDUCTION(+:xsum3,psum3,nsum3)
        call new_jtensor(jt, mol, xdens)
        do k=1,p3
            xsum2=0.d0
            psum2=0.d0
            nsum2=0.d0

            x_prime_sum2=0.d0
            p_prime_sum2=0.d0
            n_prime_sum2=0.d0
            x_second_sum2=0.d0
            p_second_sum2=0.d0
            n_second_sum2=0.d0

            !!$OMP DO SCHEDULE(STATIC) REDUCTION(+:xsum2,psum2,nsum2)
            do j=lo,hi
                xsum=0.d0
                psum=0.d0
                nsum=0.d0

                x_prime_sum=0.d0
                p_prime_sum=0.d0
                n_prime_sum=0.d0
                x_second_sum=0.d0
                p_second_sum=0.d0
                n_second_sum=0.d0

                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    r=sqrt(sum((rr-center)**2))
                    call ctensor(jt, rr, tt, spin)
                    jvec=matmul(reshape(tt,(/3,3/)),bb)
                    jvec_prime=matmul(reshape(tt,(/3,3/)),bb_prime)
                    jvec_second=matmul(reshape(tt,(/3,3/)),bb_second)
                    if ( r > bound ) then
                        w=0.d0
                        jp=0.d0
                        jp_prime=0.d0
                        jp_second=0.d0
                    else
                        w=get_weight(this%grid, i, 1)
                        jp=dot_product(normal,jvec)*w
                        jp_prime=dot_product(normal,jvec_prime)*w
                        jp_second=dot_product(normal,jvec_second)*w
                    end if
                    ! total
                    xsum=xsum+jp
                    x_prime_sum = x_prime_sum + jp_prime
                    x_second_sum = x_second_sum + jp_second
                    if (jp > 0.d0) then
                        ! get positive contribution
                        psum=psum+jp
                        p_prime_sum = p_prime_sum + jp_prime
                        p_second_sum = p_second_sum + jp_second

                    else
                        ! get negative contribution
                        nsum=nsum+jp
                        n_prime_sum = n_prime_sum + jp_prime
                        n_second_sum = n_second_sum + jp_second
                    end if
                    end do
                w=get_weight(this%grid,j,2)
                xsum2=xsum2+xsum*w
                psum2=psum2+psum*w
                nsum2=nsum2+nsum*w
                
                x_prime_sum2=x_prime_sum2+x_prime_sum*w
                p_prime_sum2=p_prime_sum2+p_prime_sum*w
                n_prime_sum2=n_prime_sum2+n_prime_sum*w

                x_second_sum2=x_second_sum2+x_second_sum*w
                p_second_sum2=p_second_sum2+p_second_sum*w
                n_second_sum2=n_second_sum2+n_second_sum*w
                
            end do
            !!$OMP END DO
            !!$OMP MASTER
            call collect_sum(xsum2, xsum)
            call collect_sum(psum2, psum)
            call collect_sum(nsum2, nsum)
            w=get_weight(this%grid,k,3)
            xsum3=xsum3+xsum*w
            psum3=psum3+psum*w
            nsum3=nsum3+nsum*w

            x_prime_sum3=x_prime_sum3+x_prime_sum*w
            p_prime_sum3=p_prime_sum3+p_prime_sum*w
            n_prime_sum3=n_prime_sum3+n_prime_sum*w

            x_second_sum3=x_second_sum3+x_second_sum*w
            p_second_sum3=p_second_sum3+p_second_sum*w
            n_second_sum3=n_second_sum3+n_second_sum*w
            !!$OMP END MASTER
        end do
        call del_jtensor(jt)
!!$OMP END PARALLEL

write (*,*) 'DEBUGGING'

        call nl
        call msg_out(repeat('*', 60))
!        write(str_g, '(a,3f10.5)') '   Magnetic field <x,y,z> =', bb
!        call msg_out(str_g)
        write(str_g, '(a,f15.8)') '   Induced current (au)    :', xsum3
        call msg_out(str_g)
        write(str_g, '(a,f15.8,a,f13.8,a)') &
            '      Positive contribution:', psum3, '  (',au2si(psum3),' )'
        call msg_out(str_g)
        write(str_g, '(a,f15.8,a,f13.8,a)') &
            '      Negative contribution:', nsum3, '  (',au2si(nsum3),' )'
        call msg_out(str_g)
        call nl
        write(str_g, '(a,f15.8)') '   Induced current (nA/T)  :', au2si(xsum3)
        call msg_out(str_g)
        write(str_g, '(a,f15.8)') '      (conversion factor)  :', au2si(1.d0)
        call msg_out(str_g)
        call nl

        call nl
        call msg_out(repeat('*', 60))
        call nl
        write(str_g, '(a,3f10.5)') '   Magnetic field <x,y,z> =', bb_prime
        call msg_out(str_g)
        write(str_g, '(a,f15.8)') '   Induced j_prime current (au)    :', x_prime_sum3
        call msg_out(str_g)
        write(str_g, '(a,f15.8,a,f13.8,a)') &
            '      Positive j_prime contribution:', p_prime_sum3, '  (',au2si(p_prime_sum3),' )'
        call msg_out(str_g)
        write(str_g, '(a,f15.8,a,f13.8,a)') &
            '      Negative j_prime contribution:', n_prime_sum3, '  (',au2si(n_prime_sum3),' )'
        call msg_out(str_g)
        call nl
        write(str_g, '(a,f15.8)') '   Induced j_prime current (nA/T)  :', au2si(x_prime_sum3)
        call msg_out(str_g)
        call nl
        call msg_out(repeat('*', 60))
        call nl
        write(str_g, '(a,3f10.5)') '   Magnetic field <x,y,z> =', bb_second
        call msg_out(str_g)
        write(str_g, '(a,f15.8)') '   Induced j_second current (au)    :', x_second_sum3
        call msg_out(str_g)
        write(str_g, '(a,f15.8,a,f13.8,a)') &
            '      Positive j_second contribution:', p_second_sum3, '  (',au2si(p_second_sum3),' )'
        call msg_out(str_g)
        write(str_g, '(a,f15.8,a,f13.8,a)') &
            '      Negative j_second contribution:', n_second_sum3, '  (',au2si(n_second_sum3),' )'
        call msg_out(str_g)
        call nl
        write(str_g, '(a,f15.8)') '   Induced j_second current (nA/T)  :', au2si(x_second_sum3)
        call msg_out(str_g)
        call nl
        call msg_out(repeat('*', 60))

        spin = 'total'
    end subroutine

    ! integrate the modulus of the current, retaining the sign
    ! test version
    subroutine integrate_modulus(this, mol, xdens, spinn)
        type(integral_t), intent(inout) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens
        character(*), optional :: spinn

        integer(I4) :: i, j, k, p1, p2, p3, lo, hi
        real(DP), dimension(3) :: normal, rr, center, bb
        real(DP) :: psum, nsum, w, jp, r, bound, sgn
        real(DP) :: psum2, nsum2
        real(DP) :: psum3, nsum3
        real(DP) :: xsum, xsum2, xsum3
        real(DP), dimension(3) :: jvec
        real(DP), dimension(9) :: tt
        type(jtensor_t) :: jt

        if (present(spinn)) spin = spinn

        if (settings%is_uhf) then
            select case(spin)
                case('total')
                    call msg_note("Integrating total density")
                case('alpha')
                    call msg_note("Integrating alpha density")
                case('beta')
                    call msg_note("Integrating beta density")
                case('spindens')
                    call msg_note("Integrating spin density")
                case default
                    call msg_error("Invalid spin: " // spin)
                    stop
            end select
        end if

        call get_grid_size(this%grid, p1, p2, p3)
        call get_magnet(this%grid, bb)

        normal=get_grid_normal(this%grid)

        bound=this%grid%radius
        if (bound < 1.d+10) then
            write(str_g, *) 'Integration bound set to radius ', bound
            call msg_out(str_g)
        end if
        call grid_center(this%grid,center)
        call schedule(p2, lo, hi)
        if (mpi_rank == 0) then
            lo = 1
        end if

        sgn=1.d0
        xsum3=0.d0
        psum3=0.d0
        nsum3=0.d0
!$OMP PARALLEL PRIVATE(i,j,k,r,rr,sgn,xsum,psum,nsum) &
!$OMP PRIVATE(jt,w,jp,tt,jvec) &
!$OMP SHARED(xsum2,psum2,nsum2) &
!$OMP SHARED(p1,p2,p3,this,center,spin,bb,normal,lo,hi) &
!$OMP REDUCTION(+:xsum3,psum3,nsum3)
        call new_jtensor(jt, mol, xdens)
        do k=1,p3
            xsum2=0.d0
            psum2=0.d0
            nsum2=0.d0
            !$OMP DO SCHEDULE(STATIC) REDUCTION(+:xsum2,psum2,nsum2)
            do j=lo,hi
                xsum=0.d0
                psum=0.d0
                nsum=0.d0
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    r=sqrt(sum((rr-center)**2))
                    call ctensor(jt, rr, tt, spin)
                    jvec=matmul(reshape(tt,(/3,3/)),bb)
                    if ( r > bound ) then
                        w=0.d0
                    else
                        w=get_weight(this%grid, i, 1)
                        jp=dot_product(normal,jvec)
                        if (abs(jp) < 1.d-12) then ! prob. parallel component
                            sgn=0.d0
                        else if (jp > 0) then
                            sgn=1.d0
                        else
                            sgn=-1.d0
                        end if
                    end if
                    jp=sgn*sqrt(sum(jvec**2))
                    xsum=xsum+jp*w
                    if (jp > 0.d0) then
                        psum=psum+jp*w
                    else
                        nsum=nsum+jp*w
                    end if
                end do
                w=get_weight(this%grid,j,2)
                xsum2=xsum2+xsum*w
                psum2=psum2+psum*w
                nsum2=nsum2+nsum*w
            end do
            !$OMP END DO
            !$OMP MASTER
            call collect_sum(xsum2, xsum)
            call collect_sum(psum2, psum)
            call collect_sum(nsum2, nsum)
            w=get_weight(this%grid,k,3)
            xsum3=xsum3+xsum*w
            psum3=psum3+psum*w
            nsum3=nsum3+nsum*w
            !$OMP END MASTER
        end do
        call del_jtensor(jt)
!$OMP END PARALLEL

        call nl
        call msg_out(repeat('*', 60))
        write(str_g, '(a,f13.6)') 'Induced mod current (au)   :', xsum3
        call msg_out(str_g)
        write(str_g, '(a,f13.6,a,f11.6,a)') &
            '      Positive contribution:', psum3, '  (',au2si(psum3),' )'
        call msg_out(str_g)
        write(str_g, '(a,f13.6,a,f11.6,a)') &
            '      Negative contribution:', nsum3, '  (',au2si(nsum3),' )'
        call msg_out(str_g)
        call nl
        write(str_g, '(a,f13.6)') 'Induced mod current (nA/T) :', au2si(xsum3)
        call msg_out(str_g)
        write(str_g, '(a,f13.6)') '      (conversion factor)  :', au2si(1.d0)
        call msg_out(str_g)
        call msg_out(repeat('*', 60))
        call nl
        spin = 'total'
    end subroutine

    subroutine integrate_tensor_field(this, mol, xdens)
        type(integral_t) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens

        integer(I4) :: i, j, k, p1, p2, p3
        real(DP), dimension(3) :: rr
        real(DP), dimension(9) :: xsum
        real(DP), dimension(:,:), allocatable  :: jt1, jt2, jt3
        type(jtensor_t) :: jt

        !call jfield_eta(this%jf)
        call get_grid_size(this%grid, p1, p2, p3)

        call new_jtensor(jt, mol, xdens)
        allocate(jt1(9,p1))
        allocate(jt2(9,p2))
        allocate(jt3(9,p3))

        do k=1,p3
            do j=1,p2
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    call ctensor(jt, rr, jt1(:,i), 'total')
                end do
                jt2(:,j)=int_t_1d(jt1,this%grid,1)
            end do
            jt3(:,k)=int_t_1d(jt2,this%grid,2)
        end do
        xsum=int_t_1d(jt3,this%grid,3)

        call print_tensor_int(xsum)
        call del_jtensor(jt)
        deallocate(jt1, jt2, jt3)
    end subroutine

    function int_t_1d(jt, grid, axis) result(xsum)
        real(DP), dimension(:,:), intent(in) :: jt
        type(grid_t), intent(in) :: grid
        integer(I4), intent(in) :: axis
        real(DP), dimension(9) :: xsum

        integer(I4) :: i, pts
        real(DP) :: w

        pts=size(jt)

        xsum=0.d0
        do i=1,pts
            w=get_weight(grid, i, axis)
            xsum=xsum+jt(:,i)*w
        end do
    end function

    subroutine print_tensor_int(xsum)
        real(DP), dimension(3,3) :: xsum

        integer(I4) :: i
        print *, xsum(1,:)
        print *, xsum(2,:)
        print *, xsum(3,:)

        call nl
        call msg_out(repeat('*', 70))
        write(str_g, '(a,3f13.6)') '   Induced current (au)    :', xsum(1,:)
        call msg_out(str_g)
        write(str_g, '(a,3f13.6)') '                            ', xsum(2,:)
        call msg_out(str_g)
        write(str_g, '(a,3f13.6)') '                            ', xsum(3,:)
        call msg_out(str_g)
        call nl
        write(str_g, '(a,3f13.6)') '   Induced current (nA/T)  :', &
            (au2si(xsum(1,i)), i=1,3)
        call msg_out(str_g)
        write(str_g, '(a,3f13.6)') '                            ', &
            (au2si(xsum(2,i)), i=1,3)
        call msg_out(str_g)
        write(str_g, '(a,3f13.6)') '                            ', &
            (au2si(xsum(3,i)), i=1,3)
        call msg_out(str_g)
        call nl
        write(str_g, '(a,f13.6)') '        Conversion Factor   :', au2si(1.d0)
        call msg_out(str_g)
        call msg_out(repeat('*', 70))
        call nl
    end subroutine

    ! do here ACID integration here !
    subroutine integrate_acid(this, mol, xdens)
    use acid_module
    ! this subroutine is based on integrate_current !
        type(integral_t), intent(inout) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens

        integer(I4) :: i, j, k, p1, p2, p3, lo, hi
        real(DP), dimension(3) :: rr, center
        real(DP) :: w, r, bound
        real(DP) :: xsum, xsum2, xsum3
        real(DP) :: val, acid_val
        real(DP), dimension(9) :: tt
        type(jtensor_t) :: jt

        if (settings%is_uhf) then
            select case(spin)
                case('total')
                    call msg_note("Integrating total density")
                case('alpha')
                    call msg_note("Integrating alpha density")
                case('beta')
                    call msg_note("Integrating beta density")
                case('spindens')
                    call msg_note("Integrating spin density")
                case default
                    call msg_error("Invalid spin: " // spin)
                    stop
            end select
        end if

        call get_grid_size(this%grid, p1, p2, p3)

        bound=1.d+10
        bound=this%grid%radius
        if (bound < 1.d+10) then
            write(str_g, *) 'Integration bound set to radius ', bound
            call msg_out(str_g)
        end if

        call grid_center(this%grid,center)
        call schedule(p2, lo, hi)
        if (mpi_rank == 0) then
            lo = 1
        end if

        xsum3=0.d0
! NB! Untested parallelization, check with valgrind!
!$OMP PARALLEL PRIVATE(i,j,k,r,rr,xsum) &
!$OMP PRIVATE(jt,w,tt) &
!$OMP SHARED(xsum2) &
!$OMP SHARED(p1,p2,p3,this,center,spin,mol,xdens,lo,hi) &
!$OMP REDUCTION(+:xsum3)
        call new_jtensor(jt, mol, xdens)
        do k=1,p3
            xsum2=0.d0
            !$OMP DO SCHEDULE(STATIC) REDUCTION(+:xsum2)
            do j=lo,hi
                xsum=0.d0
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    r=sqrt(sum((rr-center)**2))
                    call ctensor(jt, rr, tt, spin)
                    ! attention: output of get_acid is in au !
                    val = get_acid(rr,tt)
                    if ( r > bound ) then
                        w=0.d0
                    else
                        w=get_weight(this%grid, i, 1)
                    end if
                    xsum=xsum+val*w
                end do
                w=get_weight(this%grid,j,2)
                xsum2=xsum2+xsum*w
            end do
            !$OMP END DO
            !$OMP MASTER
            call collect_sum(xsum2, xsum)
            w=get_weight(this%grid,k,3)
            xsum3=xsum3+xsum*w
            !$OMP END MASTER
        end do
        call del_jtensor(jt)
        acid_val = dsqrt(xsum3)
!$OMP END PARALLEL

        call nl
        call msg_out(repeat('*', 60))
        write(str_g, '(a,f13.6)') '   ACID (au) sqrt(delta J^2):', acid_val
        call msg_out(str_g)
        write(str_g, '(a,f13.6)') '   ACID (nA/T)              :', au2si(acid_val)
        call msg_out(str_g)
        call nl
        call msg_out(repeat('*', 60))
        call nl
    end subroutine

end module

! vim:et:sw=4:ts=4
