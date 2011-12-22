! Integrate the current density using quadrature.
!
! Written by Jonas Juselius, University of Helsinki, 2003 (I think).
!

module integral_class
    use globals_m
    use settings_m
    use grid_class
    use jfield_class
    use jtensor_class
    use gaussint_m
    use lip_m
    use teletype_m
    use parallel_m
    use magnet_m
    implicit none

    type integral_t
        type(jtensor_t), pointer :: jt
        type(grid_t), pointer :: grid
        real(DP), dimension(3,3) :: tint
    end type

    public new_integral, del_integral, integral_t
    public integrate_tensor_field, integrate_current
    public integrate_modulus
    
    private 
    
    character(8) :: spin='total'
    integer(I4) :: nlip
contains
    subroutine new_integral(this, jt, grid)
        type(integral_t) :: this
        type(jtensor_t), target :: jt
        type(grid_t), target :: grid
        
        this%jt=>jt
        this%grid=>grid
        this%tint=D0
    end subroutine

    subroutine del_integral(this)
        type(integral_t) :: this
    end subroutine

    subroutine integrate_current(this)
        type(integral_t), intent(inout) :: this

        integer(I4) :: i, j, k, p1, p2, p3
        real(DP), dimension(3) :: normal, rr, center, bb
        real(DP) :: psum, nsum, w, jp, r, bound
        real(DP) :: psum2, nsum2
        real(DP) :: psum3, nsum3
        real(DP) :: xsum, xsum2, xsum3
        type(vector_t) :: jvec
        type(tensor_t) :: jt

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

        bound=1.d+10
        bound=this%grid%radius
        if (bound < 1.d+10) then
            write(str_g, *) 'Integration bound set to radius ', bound
            call msg_out(str_g)
        end if
        
        call grid_center(this%grid,center)

        xsum3=0.d0
        psum3=0.d0
        nsum3=0.d0
        do k=1,p3
            xsum2=0.d0
            psum2=0.d0
            nsum2=0.d0
            do j=1,p2
                xsum=0.d0
                psum=0.d0
                nsum=0.d0
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    r=sqrt(sum((rr-center)**2))
                    call ctensor(this%jt, rr, jt, spin)
                    jvec%v=matmul(jt%t,bb)
                    if ( r > bound ) then
                        w=0.d0
                        jp=0.d0
                    else
                        w=get_weight(this%grid, i, 1) 
                        jp=dot_product(normal,jvec%v)*w
                    end if
                    xsum=xsum+jp
                    if (jp > 0.d0) then
                        psum=psum+jp
                    else
                        nsum=nsum+jp
                    end if
                end do
                w=get_weight(this%grid,j,2)
                xsum2=xsum2+xsum*w
                psum2=psum2+psum*w
                nsum2=nsum2+nsum*w
            end do
            w=get_weight(this%grid,k,3)
            xsum3=xsum3+xsum*w
            psum3=psum3+psum*w
            nsum3=nsum3+nsum*w
        end do

        call nl
        call msg_out(repeat('*', 60))
        write(str_g, '(a,f13.6)') '   Induced current (au)    :', xsum3
        call msg_out(str_g)
        write(str_g, '(a,f13.6,a,f11.6,a)') &
            '      Positive contribution:', psum3, '  (',au2si(psum3),' )'
        call msg_out(str_g)
        write(str_g, '(a,f13.6,a,f11.6,a)') &
            '      Negative contribution:', nsum3, '  (',au2si(nsum3),' )'
        call msg_out(str_g)
        call nl
        write(str_g, '(a,f13.6)') '   Induced current (nA/T)  :', au2si(xsum3)
        call msg_out(str_g)
        write(str_g, '(a,f13.6)') '      (conversion factor)  :', au2si(1.d0)
        call msg_out(str_g)
        call msg_out(repeat('*', 60))
    end subroutine

    ! integrate the modulus of the current, retaining the sign
    ! test version
    subroutine integrate_modulus(this)
        type(integral_t), intent(inout) :: this

        integer(I4) :: i, j, k, p1, p2, p3
        real(DP), dimension(3) :: normal, rr, center, bb
        real(DP) :: psum, nsum, w, jp, r, bound, sgn
        real(DP) :: psum2, nsum2
        real(DP) :: psum3, nsum3
        real(DP) :: xsum, xsum2, xsum3
        type(vector_t) :: jvec
        type(tensor_t) :: jt

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

        sgn=1.d0
        xsum3=0.d0
        psum3=0.d0
        nsum3=0.d0
        do k=1,p3
            xsum2=0.d0
            psum2=0.d0
            nsum2=0.d0
            do j=1,p2
                xsum=0.d0
                psum=0.d0
                nsum=0.d0
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    r=sqrt(sum((rr-center)**2))
                    call ctensor(this%jt, rr, jt, spin)
                    jvec%v=matmul(jt%t,bb)
                    if ( r > bound ) then
                        w=0.d0
                    else
                        w=get_weight(this%grid, i, 1) 
                        jp=dot_product(normal,jvec%v)
                        if (abs(jp) < 1.d-12) then ! prob. parallel component
                            sgn=0.d0
                        else if (jp > 0) then
                            sgn=1.d0
                        else
                            sgn=-1.d0
                        end if
                    end if
                    jp=sgn*sqrt(sum(jvec%v**2))
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
            w=get_weight(this%grid,k,3)
            xsum3=xsum3+xsum*w
            psum3=psum3+psum*w
            nsum3=nsum3+nsum*w
        end do

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
    end subroutine

    subroutine integrate_tensor_field(this)
        type(integral_t) :: this

        integer(I4) :: i, j, k, p1, p2, p3
        real(DP), dimension(3) :: rr
        real(DP), dimension(3,3) :: xsum
        type(tensor_t), dimension(:), allocatable  :: jt1, jt2, jt3
        
        !call jfield_eta(this%jf)
        call get_grid_size(this%grid, p1, p2, p3)

        allocate(jt1(p1))
        allocate(jt2(p2))
        allocate(jt3(p3))

        do k=1,p3
            do j=1,p2
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    call ctensor(this%jt, rr, jt1(i), 'total')
                end do
                jt2(j)%t=int_t_1d(jt1,this%grid,1)
            end do
            jt3(k)%t=int_t_1d(jt2,this%grid,2)
        end do
        xsum=int_t_1d(jt3,this%grid,3)

        call print_tensor_int(xsum)
        deallocate(jt1, jt2, jt3)
    end subroutine

    function int_t_1d(jt, grid, axis) result(xsum)
        type(tensor_t), dimension(:), intent(in) :: jt
        type(grid_t), intent(in) :: grid
        integer(I4), intent(in) :: axis
        real(DP), dimension(3,3) :: xsum

        integer(I4) :: i, pts
        real(DP) :: w

        pts=size(jt)

        xsum=0.d0
        do i=1,pts
            w=get_weight(grid, i, axis) 
            xsum=xsum+jt(i)%t*w
        end do
    end function

    function au2si(au) result(si)
        real(DP), intent(in) :: au
        real(DP) :: si
        
        real(DP) :: aulength, auspeedoflight, speedoflight, aucharge, hbar
        real(DP) :: autime, autesla, audjdb

        aulength=0.52917726D-10
        auspeedoflight=137.03599D0
        speedoflight=299792458.D0
        aucharge=1.60217733D-19
        hbar=1.05457267D-34

        autime=aulength*auspeedoflight/speedoflight
        autesla=hbar/aucharge/aulength/aulength
        audjdb=aucharge/autime/autesla

        si=au*audjdb*1.d+09 ! nA/T

!        write(6,*) 'The obtained conversion factors'
!        write(6,*) autime,' au time in seconds'
!        write(6,*) autesla,' au magnetic field in tesla'
!        write(6,*) audjdb*1.D+09,' au induced current in nanoampere/tesla'
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
    end subroutine
end module

! vim:et:sw=4:ts=4
