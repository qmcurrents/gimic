module magnet_m
	use kinds_m
	use globals_m
	use grid_class
	implicit none
	
	public get_magnet
	private
contains
	subroutine get_magnet(g,mag)
		type(grid_t) :: g
		real(DP), dimension(3), intent(out) :: mag

		character(2) :: axis
		real(DP), dimension(3) :: i, j, k, v
		real(DP) :: dir
		logical :: ortho

		ortho=.false.
		dir=1.d0
		mag=0.d0
		if (keyword_is_set(input, 'magnet_axis')) then
			call getkw(input, 'magnet_axis', axis)
			if (axis(1:1) == '-') then 
				dir=-1.d0
				axis(1:1)=axis(2:2)
			end if

			call get_basvec(g, i, j, k)
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
				case('T')
					ortho=.true.
					call get_ortho(g,mag)
					mag=mag*dir
				case default
					call msg_error('Invalid axis specifier: ' // axis)
					stop
			end select
		else if (keyword_is_set(input, 'magnet')) then
			call getkw(input, 'magnet', mag)
		end if

		if (vcmp(mag, NILL_VECTOR)) then
			call msg_critical('Magnetic field is zero, not wasting more CPU.')
			stop
		end if

		if (ortho) then
			call get_basvec(g, 2, v)
		else
			call get_basvec(g, 3, v)
		end if

		call check_field(v,mag)
		write(str_g, '(a,3f10.5)') '   Magnetic field <x,y,z> =', mag
		call msg_out(str_g)
		call nl
	end subroutine

	subroutine check_field(dir,mag) 
		real(DP), dimension(3), intent(in) :: dir
		real(DP), dimension(3), intent(inout) :: mag

		real(DP) :: x

		x=dot_product(dir, mag)
		if (x > 0.d0) then
			mag=-mag
			call msg_info('Left handed coordinate system,&
				& reversing magnetic field')
		end if
		if (abs(x)-1.d0 > 1.d-12 .and. abs(x) > 1.d-12) then
			call msg_warn('Magnetic field not orthogonal to grid')
		end if
	end subroutine

	function is_diamagnetic(g)  result(dia)
		type(grid_t) :: g
		logical :: dia

		real(DP), dimension(3) :: magnet,v
		real(DP) :: x

		call push_section(input, '.')
		call get_magnet(g, magnet) 
		call pop_section(input)
		call get_basvec(g, 3, v)

		stop 'is_diamagnetic(): not implemented yet'
		dia=.true.
		x=dot_product(v, magnet)
		if (x < 0.d0) then
			dia=.false.
		end if
		
	end function

end module
