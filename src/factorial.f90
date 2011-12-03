
module factorial_m
	use kinds_m
	implicit none
	public fact, fact2, binom
	contains
	
	function fact(n) result(m)
		integer(I4), intent(in) :: n
		integer(I4) :: i, m

		if (n > 0 ) then
			m=n
			do i=m-1, 2, -1
				m=m*i
			end do
		else
			m=1
		endif
	end function
	
	subroutine fact2(n, m)
		integer(I4), intent(inout) :: n
		integer(I4), intent(in) :: m
		integer(I4) :: i
		
		if (n > 0 ) then
			do i=n-1, m, -1
				n=n*i
			end do
		else
			n=1
		endif
	end subroutine

	function binom(a, b) result(c)
		integer(I4), intent(in) :: a, b
		integer(I8) :: c
		
		c=fact(a)/(fact(b)*fact(a-b))
	end function

end module

