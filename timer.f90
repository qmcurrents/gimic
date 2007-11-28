module timer_m
	use globals_m
	use teletype_m
	implicit none
	external etime, dtime, fdate

	real(4) :: etime, dtime
	character(BUFLEN) :: fdate
	
	real(4) :: delta_t
	real(4), dimension(2) :: times

	public stockas_klocka, swatch
	private
contains
	subroutine stockas_klocka()
		delta_t=etime(times)
		call msg_out(repeat('-', 70))
		write(str_g, '(a,f9.2,a)') '   wall time:', delta_t, 'sec'
		call msg_out(str_g)
		write(str_g, '(a,f9.2,a)') '        user:', times(1), 'sec'
		call msg_out(str_g)
		write(str_g, '(a,f9.2,a)') '         sys:', times(2), 'sec'
		call msg_out(str_g)
		call msg_out(repeat('-', 70))
		call msg_out(fdate())
	end subroutine 

	subroutine swatch(foo)
		integer(I4), intent(in) :: foo
		
		delta_t=dtime(times)
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

