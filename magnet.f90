module magnet_m
	use kinds_m
	use globals_m
	use grid_class
	implicit none
	
	public get_magnet
	private
	real(DP), dimension(3), save :: magnet=0.d0
contains
	subroutine get_magnet(g,mag)
		type(grid_t) :: g
		real(DP), dimension(3), intent(out) :: mag

		character(2) :: axis
		real(DP), dimension(3) :: i, j, k
		real(DP) :: dir=1.d0

		if (keyword_is_set(input, 'cdens.magnet_axis')) then
			call getkw(input, 'cdens.magnet_axis', axis)
			if (axis(1:1) == '-') then 
				dir=-1.d0
				axis(1:1)=axis(2:2)
			end if

			call get_basv(g, i, j, k)
			select case(axis(1:1))
				case('i') 
					mag=i*dir
				case('j')
					mag=j*dir
				case('k')
					mag=k*dir
				case('x')
					mag=(/D1, D0, D0/)*dir
				case('y')
					mag=(/D0, D1, D0/)*dir
				case('z')
					mag=(/D0, D0, D1/)*dir
				case default
					call msg_error('Invalid axis specifier: ' // axis)
					stop
			end select
		else
			call getkw(input, 'cdens.magnet', mag)
		end if

		call check_field(g,mag)
		if (.not.vcmp(mag, magnet)) then
			magnet=mag
			write(str_g, '(a,3f10.5)') '   Magnetic field <x,y,z> =', mag
			call msg_out(str_g)
		end if
	end subroutine

	subroutine check_field(g,mag) 
		type(grid_t) :: g
		real(DP), dimension(3), intent(inout) :: mag

		real(DP) :: x

		x=dot_product(g%basv(:,3), mag)
		if (x < 0.d0) then
			mag=-mag
			call msg_info('Left handed coordinate system,&
		    	& reversing magnetic field')
		end if
		if (abs(x) /= 1.d0 .and. abs(x) > 1.d-10) then
			call msg_warn('Magnetic field not orthogonal to grid')
		end if
	end subroutine

	function is_diamagnetic(g)  result(dia)
		type(grid_t) :: g
		logical :: dia

		real(DP), dimension(3) :: magnet
		real(DP) :: x

		call push_section(input, '.')
		call get_magnet(g, magnet) 
		call pop_section(input)

		dia=.true.
		x=dot_product(g%basv(:,3), magnet)
		if (x < 0.d0) then
			dia=.false.
		end if
		
	end function

end module
