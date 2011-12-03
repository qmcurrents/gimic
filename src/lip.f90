
module lip_m
	use globals_m
	implicit none
	public lipn, lip3n, vector_t, lip3n_t

	integer(I2), parameter :: MAX_LIP_N=11

	! max lip order = 11, sanity check needed...
	type lip3n_t
		type(vector_t), dimension(MAX_LIP_N) :: grid
		type(vector_t), dimension(MAX_LIP_N) :: value
	end type
	
	private
contains
	
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
	! General N-point Lagrange Inperpolation Polynomial generator
	!
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function lipn(x, g, z) result(f)
		real(DP), intent(in) :: x
		real(DP), dimension(:), intent(in) :: g
		real(DP), dimension(:), intent(in) :: z

		real(DP) :: f, q
		integer(I4) :: j, k, npoints

!       g is the grid point
!       z is the func val

		npoints=size(g)

		f=0.d0
		do j=1,npoints
			q=1.d0
			do k=1,npoints
				if (k /= j) q=q*(x-g(k))/(g(j)-g(k))
			end do
			f=f+q*z(j)
		end do
	end function

	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
	! General N-point Lagrange Inperpolation Polynomial generator
	! for vector valued functions
	! 
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function lip3n(x, g, z) result(f)
		real(DP), intent(in) :: x
		real(DP), dimension(:), intent(in) :: g
		type(vector_t), dimension(:), intent(in) :: z
		real(DP), dimension(3) :: f

		real(DP) :: q
		integer(I4) :: j, k, npoints

		npoints=size(g)

		f=0.d0
		do j=1,npoints
			q=1.d0
			do k=1,npoints
				if (k /= j) q=q*(x-g(k))/(g(j)-g(k))
			end do
			f=f+q*z(j)%v
		end do
	end function

end module
