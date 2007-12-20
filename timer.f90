module timer_m
	use globals_m
	use teletype_m
	implicit none

	character(BUFLEN) :: date
	
	real(4) :: delta_t
	real(4), dimension(2) :: times

	public stockas_klocka, swatch
	private
contains
	subroutine stockas_klocka()
		call etime(times, delta_t)
		call msg_out(repeat('-', 70))
		write(str_g, '(a,f9.2,a)') '   wall time:', delta_t, 'sec'
		call msg_out(str_g)
		write(str_g, '(a,f9.2,a)') '        user:', times(1), 'sec'
		call msg_out(str_g)
		write(str_g, '(a,f9.2,a)') '         sys:', times(2), 'sec'
		call msg_out(str_g)
		call msg_out(repeat('-', 70))
		call fdate(date)
		call msg_out(trim(date))
	end subroutine 

	subroutine swatch(foo)
		integer(I4), intent(in) :: foo
		
		call dtime(times, delta_t)
		if (foo > 0) then
			call msg_out(repeat('-', 70))
			write(str_g, '(a,f9.2,a)') '   wall time:', delta_t, 'sec'
			call msg_out(str_g)
			write(str_g, '(a,f9.2,a)') '        user:', times(1), 'sec'
			call msg_out(str_g)
			write(str_g, '(a,f9.2,a)') '         sys:', times(2), 'sec'
			call msg_out(str_g)
			call msg_out(repeat('-', 70))
		end if
	end subroutine 
end module

