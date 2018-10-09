module gaussint_module
    use globals_module
    use settings_module
    use teletype_module
    implicit none

    type gdata_t
        real(DP), dimension(:), pointer :: pts, wgt
    end type

    public gdata_t, gaussl, lobattomy, legendrep, gaussgrid, gauscale
    public setup_gauss_data

    private
    real(DP), parameter :: EPS=3.0d-12
    integer(I4), parameter :: NEWTON_MAX_ITER=10
    integer(I4), parameter :: MAX_LOBATTO_POINTS=13 ! > 13 is unstable...

contains

    subroutine gaussl(a,b,pts,weight)
        real(DP), intent(in) :: a,b
        real(DP), dimension(:), intent(out) :: pts, weight

        integer(I4) :: n, m, i, iter
        real(DP) :: z, z1, xm, xl, lp, dlp

        n=size(pts)
        if (mod(n,2) == 0) then
            m=n/2
        else
            m=(n+1)/2
        end if

        xm=(b+a)*0.5
        xl=(b-a)*0.5

        do i=1,m
            z=cos(PII*(i-0.25)/(n+0.5))
            do iter=1,NEWTON_MAX_ITER
                call legendrep1(z,n,lp,dlp)
                z1=z
                z=z1-lp/dlp

                if (abs(z-z1) <= EPS) then
                    exit
                end if
            end do
            if ( iter >= NEWTON_MAX_ITER ) then
                stop '*** integration did not converge!'
            end if

            pts(i)=xm-xl*z
            pts(n+1-i)=xm+xl*z

            weight(i)=2.d0*xl/((1.d0-z**2)*dlp**2)
            weight(n+1-i)=weight(i)
        end do
    end subroutine gaussl

    subroutine lobattomy(a,b,pts,weight)
        real(DP), intent(in) :: a,b
        real(DP), dimension(:), intent(out) :: pts, weight

        integer(I4) :: n, i, iter
        real(DP) :: z, z1, xm, xl, lp, dlp, d2lp, damp

        n=size(pts)
        if ( n > MAX_LOBATTO_POINTS ) then
            write(str_g,'(a,i4)') &
                'lobattomy(): Warning, too many Gauss points:', n
            call msg_warn(str_g)
            call msg_warn('lobattomy(): &
                &using this many points is not stable.')
            call nl
        end if

        xm=(b+a)*0.5
        xl=(b-a)*0.5

        ! end points...
        pts(1)=xm-xl
        pts(n)=xm+xl
        weight(1)=2.d0*xl/(n**2-n)
        weight(n)=weight(1)

        do i=2,n-1
            ! Better inital guess needed, this one is for P not P'
            z=cos(PII*(i-0.25)/(n+0.5))
            do iter=1,NEWTON_MAX_ITER
                call legendrep2(z,n-1,lp,dlp,d2lp)
                z1=z
                z=z1-dlp/d2lp
                ! damping...
                damp=0.5d0
                do while (abs(z) > 1.d0)
                    z=z1-dlp/d2lp*damp
                    damp=damp**2
                end do

                if (abs(z-z1) <= EPS) then
                    exit
                end if
            end do
            if ( iter >= NEWTON_MAX_ITER ) then
                stop '*** lobattomy(): Newton step did not converge!'
            end if

            pts(i)=xm-xl*z
            weight(i)=2.d0*xl/((n**2-n)*lp**2)
        end do
    end subroutine lobattomy

    subroutine legendrep(x, n, y)
        real(DP), intent(in) :: x
        integer(I4), intent(in) :: n
        real(DP), intent(out) :: y

        real(DP) :: c1, c2, c4, ym, yp, dym, dyp, d2ym, d2yp
        integer(I4) :: i

        if ( n < 0 .or. abs(x) > 1.0 ) then
            write(str_g,'(a,i3,d19.12)') 'legendrep(): bad argument:', n, x
            call msg_error(str_g)
        end if

        y = 1.d0
        if (n == 0) return

        y = x
        if (n == 1) return

        yp = 1.d0
        do i=2,n
            c1 = real(i)
            c2 = c1 * 2.0 - 1.0
            c4 = c1 - 1.0
            ym = y
            y = (c2 * x * y - c4 * yp) / c1
            yp = ym
        end do
    end subroutine

    subroutine legendrep1(x, n, y, dy)
        real(DP), intent(in) :: x
        integer(I4), intent(in) :: n
        real(DP), intent(out) :: y
        real(DP), intent(out) :: dy

        real(DP) :: c1, c2, c4, ym, yp, dym, dyp, d2ym, d2yp
        integer(I4) :: i

        if ( n < 0 .or. abs(x) > 1.0 ) then
            write(str_g,'(a,i3,d19.12)') 'legendrep(): bad argument:', n, x
            call msg_error(str_g)
        end if

        y = 1.d0
        dy = 0.d0
        if (n == 0) return

        y = x
        dy = 1.d0
        if (n == 1) return

        yp = 1.d0
        dyp = 0.d0
        do i=2,n
            c1 = real(i)
            c2 = c1 * 2.0 - 1.0
            c4 = c1 - 1.0
            ym = y
            y = (c2 * x * y - c4 * yp) / c1
            yp = ym
            dym = dy
            dy = (c2 * x * dy - c4 * dyp + c2 * yp) / c1
            dyp = dym
        end do
    end subroutine

    subroutine legendrep2(x, n, y, dy, d2y)
        real(DP), intent(in) :: x
        integer(I4), intent(in) :: n
        real(DP), intent(out) :: y
        real(DP), intent(out) :: dy, d2y

        real(DP) :: c1, c2, c4, ym, yp, dym, dyp, d2ym, d2yp
        integer(I4) :: i

        if ( n < 0 .or. abs(x) > 1.0 ) then
            write(str_g,'(a,i3,d19.12)') 'legendrep(): bad argument:', n, x
            call msg_error(str_g)
        end if

        y = 1.d0
        dy = 0.d0
        d2y = 0.d0
        if (n == 0) return

        y = x
        dy = 1.d0
        d2y = 0.d0
        if (n == 1) return

        yp = 1.d0
        dyp = 0.d0
        d2yp = 0.d0
        do i=2,n
            c1 = real(i)
            c2 = c1 * 2.0 - 1.0
            c4 = c1 - 1.0
            ym = y
            y = (c2 * x * y - c4 * yp) / c1
            yp = ym
            dym = dy
            dy = (c2 * x * dy - c4 * dyp + c2 * yp) / c1
            dyp = dym
            d2ym = d2y
            d2y = (c2 * x * d2y - c4 * d2yp + c2 * 2.0 * dyp) / c1
            d2yp = d2ym
        end do
    end subroutine

    subroutine gaussgrid(a, b, ngp, xvec, wvec)
        real(DP), intent(in) :: a, b
        integer(I4), intent(in) :: ngp
        real(DP), dimension(:), intent(out) :: xvec, wvec

        integer(I4) :: npts, nblock, i, j, k
        real(DP) :: step, pos, xl
        real(DP), dimension(:), allocatable :: pts, wgt

        npts=size(xvec)
        allocate(pts(ngp))
        allocate(wgt(ngp))

        if (mod(npts,ngp) /= 0) then
            call msg_critical('gaussgrid(): npts is not dividable by ngp!')
            stop
        end if

        ! split range [a,b] into blocks for piecewise integration
        nblock=(npts)/ngp
        ! lenght of one block
        step=(b-a)/real(nblock)

        ! set up unscaled Lobatto points and weights ( interval [-1,1])
!        call lobattomy(-1.d0, 1.d0, pts, wgt)
        call gaussl(-1.d0, 1.d0, pts, wgt)

        k=0
        pos=a
        xl=step*0.5
        ! scale and translate Lobatto points and weights
        do i=1,nblock
            do j=1,ngp
                k=k+1
                xvec(k)=pts(j)*xl+pos+xl
                wvec(k)=wgt(j)*xl
            end do
            pos=pos+step
        end do

        deallocate(pts, wgt)
    end subroutine

    subroutine setup_gauss_data(a, b, ngp, gdata, quadr)
        real(DP), intent(in) :: a, b
        integer(I4), intent(in) :: ngp
        type(gdata_t), intent(inout) :: gdata
        character(*), intent(in) :: quadr

        integer(I4) :: nblock, i, ngau, npts, foo
        real(DP), dimension(:), allocatable :: tpts, twgt
        real(DP) :: xl, step
        logical :: lobatto

        npts=size(gdata%pts)

        if (npts == 1) then
            gdata%pts(1)=0.d0
            gdata%wgt(1)=1.d0
            return
        end if

        if (mod(npts,ngp) /= 0) then
            call msg_critical('gaussgrid(): npts is not dividable by ngp!')
            stop
        end if

        ! split range [a,b] into blocks for piecewise integration
        nblock=npts/ngp
        ! lenght of one block
        step=(b-a)/real(nblock)
        xl=step*0.5

        ! set up unscaled Lobatto points and weights ( interval [-1,1])
        allocate(tpts(ngp))
        allocate(twgt(ngp))

        select case (quadr)
            case ('gauss')
                call gaussl(-1.d0, 1.d0, tpts, twgt)
            case ('lobatto')
                call lobattomy(-1.d0, 1.d0, tpts, twgt)
            case default
                call msg_error('Invalid quadrature type: ' // quadr)
                stop
        end select
        foo=1
        tpts = tpts + 1.0
        do i=1,nblock
            gdata%pts(foo:foo+ngp-1)=tpts*xl+real(i-1)*step
            gdata%wgt(foo:foo+ngp-1)=twgt*xl
            foo=foo+ngp
        end do

        deallocate(tpts, twgt)
    end subroutine


    ! calculate scaling and translation factors for Gausspoints defined
    ! on the interval [-1,1]
    !   a, b: the interval to transform to
    !   s   : scaling factor
    !   t   : translation factor
    subroutine gauscale(a, b, s, t)
        real(DP), intent(in) :: a, b
        real(DP), intent(out) :: s, t

        s=(b-a)*0.5d0
        t=(b+a)*0.5d0
    end subroutine

!    subroutine itester
!        integer(I4) :: np, ngrid
!        integer(I4) :: i, j
!        real(DP), dimension(:) :: pts, weight
!        real(DP), dimension(:) :: tmp
!        real(DP) :: a, b, x
!        real(DP) :: ss, s, t
!        real(DP) :: p, dp, d2p
!        real(DP), parameter :: aexp=1.0d0
!        type(lobo_t) :: lobby
!        real(DP), dimension(2) :: qq
!
!        do np=23,23
!            allocate(pts(np))
!            allocate(weight(np))
!            a=1.5d0
!            call gaussl(-a, a, pts, weight)
!
!            ss=0.d0
!            do i=1,np
!                ss=ss+f(pts(i))*weight(i)
!            end do
!
!            print '(i4,3f22.15)', np, ss, sqrt(PI/aexp)
!            deallocate(pts, weight)
!        end do
!        print *
!        print *, repeat('*',72)
!        print *
!
!        do np=10,10
!            allocate(pts(np))
!            allocate(weight(np))
!            a=1.0d0
!            call lobattomy(-a, a, pts, weight)
!
!            ss=0.d0
!            b=-10.d0
!            a=abs(b+1.d-3)
!            do while (b < a)
!                call gauscale(b,b+2.d0,s,t)
!                do i=1,np
!                    x=pts(i)*s+t
!                    ss=ss+f(x)*weight(i)
!                end do
!                b=b+2.d0
!            end do
!
!            print '(i4,3f22.15)', np, ss, sqrt(PI/aexp)
!            deallocate(pts, weight)
!        end do
!        print *
!        print *, repeat('*',72)
!        print *
!
!        do np=10,10
!            ! adapt gridsize to nearest divisior
!            ngrid=200+np-mod(200,np)
!            allocate(pts(ngrid))
!            allocate(weight(ngrid))
!
!            call lobogrid(-10.d0, 10.d0, np, pts, weight)
!            lobby=new_lobby(-10.d0, 10.d0, np, ngrid)
!
!            ss=0.d0
!            do i=1,ngrid
!                ss=ss+f(pts(i))*weight(i)
!                print *, pts(i), weight(i)
!            end do
!
!            print '(i4,3f22.15)', np, ss, sqrt(PI/aexp)
!            deallocate(pts, weight)
!        end do
!        print *
!        print *, repeat('*',72)
!        print *
!
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        ! TEST LOBBY
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        np=12
!!        ngrid=200+np-mod(200,np)
!!        allocate(pts(ngrid))
!!        allocate(weight(ngrid))
!!
!!        call lobogrid(-5.d0, 5.d0, np, pts, weight)
!!        lobby=new_lobby(-5.d0, 5.d0, np, ngrid)
!!
!!        do i=1,ngrid
!!            qq=lobopoint(i,lobby)
!!            print '(4f20.13)', pts(i), qq(1), weight(i), qq(2)
!!            if (abs(pts(i)-qq(1)) > 1.d-10) print *, 'foo!'
!!            if (abs(weight(i)-qq(2)) > 1.d-10) print *, 'bar!'
!!        end do
!!
!!        deallocate(pts, weight)
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        np=10
!        ngrid=200+np-mod(200,np)
!        allocate(pts(ngrid))
!        allocate(weight(ngrid))
!        allocate(tmp(ngrid))
!        call lobogrid(-10.d0, 10.d0, np, pts, weight)
!        ! first integration
!        do np=1,ngrid
!            ss=0.d0
!            do i=1,ngrid
!                ss=ss+h(pts(i), pts(np))*weight(i)
!            end do
!            tmp(np)=ss
!        end do
!        ! second integration
!        ss=0.d0
!        do i=1,ngrid
!            ss=ss+tmp(i)*weight(i)
!        end do
!        print *
!        print *, '2D',ss
!        deallocate(pts, weight)
!        print *
!
!    contains
!        function f(x) result(y)
!            real(DP), intent(in) :: x
!            real(DP) :: y
!
!            y=exp(-aexp*x**2)!-x**2+2.d0
!        end function
!
!        function h(x,y) result(z)
!            real(DP), intent(in) :: x, y
!            real(DP) :: z
!
!            z=exp(-aexp*(x**2+y**2))
!        end function
!
!    end subroutine
end module


! vim:et:sw=4:ts=4
