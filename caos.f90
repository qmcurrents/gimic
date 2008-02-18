!
! Basis function evaluator
!

module caos_m
	use globals_m
	use gtodefs_m
	use basis_class
	use cao2sao_class
	implicit none
	
	public  cgto, dcgto

	private

	real(DP) :: rr2
	type(contraction_t), pointer :: cc

contains
	
	! Evaluate one contracted CAO 
	function cao() result(ff)
		real(DP) :: ff

		integer(I4) :: i
		real(DP), dimension(:), pointer :: ncc, xp
		ncc=>cc%ncc
		xp=>cc%xp
		
		ff=D0
		do i=1,cc%npf
			ff=ff+ncc(i)*exp(-xp(i)*rr2)
		end do
	end function
	
	! Evaluate one differentiated CAO
	function dcao() result(ff)
		real(DP) :: ff

		integer(I4) :: i
		real(DP), dimension(:), pointer :: ncc, xp
		ncc=>cc%ncc
		xp=>cc%xp
		
		ff=D0
		do i=1,cc%npf
			ff=ff-xp(i)*ncc(i)*exp(-xp(i)*rr2)
		end do
	end function

	subroutine cao2(vcao, vdcao)
		real(DP), intent(out) :: vcao, vdcao
		
		integer(I4) :: i
		real(DP) :: q
		real(DP), dimension(:), pointer :: ncc, xp
		ncc=>cc%ncc
		xp=>cc%xp
		
		vcao=D0
		vdcao=D0
		do i=1,cc%npf
			q=ncc(i)*exp(-xp(i)*rr2)
!            if (abs(q) < 1.d-20) cycle
			vcao=vcao+q
			vdcao=vdcao+xp(i)*q
		end do
	end subroutine 

	subroutine cgto(r, ctr, val)
		real(DP), dimension(3), intent(in) :: r
		type(contraction_t), intent(in), target :: ctr
		real(DP), dimension(:), intent(out) :: val

		real(DP) :: p
		real(DP) :: q
		integer(I4) :: i
		real(DP), dimension(:,:), pointer :: f

		rr2=sum(r**2)

		call get_gto_nlm(ctr%l, f)
		cc=>ctr
		q=cao()
		do i=1,Ctr%nccomp 
			p=product(r**f(:,i))*q
!            if (abs(p) > 1.d-20) val(i)=p
			val(i)=p
		end do
	end subroutine
	
	subroutine dcgto(r, ax, ctr, val)
		real(DP), dimension(3), intent(in) :: r
		integer(I4), intent(in) :: ax
		type(contraction_t), intent(in), target :: Ctr
		real(DP), dimension(:), intent(out) :: val

		real(DP) :: bfval, dbfval
		real(DP), dimension(:,:), pointer :: f
		real(DP), dimension(3) :: df
		real(DP) :: up, down
		integer(I4) :: i, j

		rr2=sum(r**2)
		
		call get_gto_nlm(ctr%l,f)
		cc=>ctr
		call cao2(bfval, dbfval)
		do i=1,ctr%nccomp 
			df=f(:,i)
			df(ax)=df(ax)-1.d0
			if (df(ax) < D0) df(ax)=D0
			down=f(ax,i)*product(r**df)*bfval
			up=2.d0*r(ax)*product(r**f(:,i))*dbfval
			val(i)=down-up
		end do
	end subroutine

end module
