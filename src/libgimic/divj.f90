!
!  Calculate the divergence of a vector field.
!  Useful for checking the convergce of the actual gauge invariance,
!  in contrast to the basis set convergence and gauge independence.
!

module divj_m
    use globals_m
    use settings_m
    use jtensor_class
    use teletype_m
    implicit none

    public divj
    private
    
    real(DP), dimension(5), parameter :: wgt=(/2.d0,-16.d0,0.d0,16.d0,-2.d0/)
    real(DP), parameter :: step=1.d-3
    real(DP), parameter :: hx=1.d0/(24.d0*step)
    real(DP), parameter :: hy=1.d0/(24.d0*step)
    real(DP), parameter :: hz=1.d0/(24.d0*step)

contains
    function divj(rr, magnet) result(div)
        real(DP), dimension(3), intent(in) :: rr, magnet
        real(DP) :: div
         
        div = divergence(rr, magnet)
    end function

    function divergence(rr, bb) result(div)
        real(DP), dimension(3), intent(in) :: rr, bb
        real(DP) :: div

        type(jtensor_t) :: jt
        integer(I4) :: q
        real(DP), dimension(9) :: jtx, jty, jtz
        real(DP) :: djx,djy,djz
        real(DP), dimension(5) :: jx,jy,jz
        real(DP), dimension(3) :: tvec

        jx=D0; jy=D0; jz=D0
        do q=-2,2
            if ( q == 0 ) cycle
            call ctensor(jt, rr+(/step*real(q),D0,D0/), jtx, 'total')
            call ctensor(jt, rr+(/D0,step*real(q),D0/), jty, 'total')
            call ctensor(jt, rr+(/D0,D0,step*real(q)/), jtz, 'total')
            tvec=matmul(reshape(jtx,(/3,3/)), bb)
            jx(q+3)=tvec(1)
            tvec=matmul(reshape(jty,(/3,3/)), bb)
            jy(q+3)=tvec(2)
            tvec=matmul(reshape(jtz,(/3,3/)), bb)
            jz(q+3)=tvec(3)
        end do
        djx=hx*sum(wgt*jx)
        djy=hy*sum(wgt*jy)
        djz=hz*sum(wgt*jz)
        div=djx+djy+djz
    end function

    function divergence2(rr, bb) result(div)
        real(DP), dimension(3), intent(in) :: rr, bb
        real(DP) :: div

        type(jtensor_t) :: jt
        integer(I4) :: q
        real(DP), dimension(9) :: jtxp, jtyp, jtzp
        real(DP), dimension(9) :: jtxd, jtyd, jtzd
        real(DP) :: djx,djy,djz
        real(DP), dimension(5) :: jx,jy,jz
        real(DP), dimension(3) :: tvec

        jx=D0; jy=D0; jz=D0
        do q=-2,2
            call ctensor2(jt, rr+(/step*real(q),D0,D0/), jtxp,jtxd, 'total')
            call ctensor2(jt, rr+(/D0,step*real(q),D0/), jtyp,jtyd, 'total')
            call ctensor2(jt, rr+(/D0,D0,step*real(q)/), jtzp,jtzd, 'total')
            call get_jvector(jtxp, jtxd, bb, tvec)
            jx(q+3)=tvec(1)
            call get_jvector(jtyp, jtyd, bb, tvec)
            jy(q+3)=tvec(2)
            call get_jvector(jtzp, jtzd, bb, tvec)
            jz(q+3)=tvec(3)
        end do
        djx=hx*sum(wgt*jx)
        djy=hy*sum(wgt*jy)
        djz=hz*sum(wgt*jz)
        div=djx+djy+djz
    end function
end module

! vim:et:sw=4:ts=4
