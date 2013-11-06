module timer_module
    use globals_module
    use settings_module
    use teletype_module
    implicit none

    real(4) :: delta_t
    real(4), dimension(2) :: times

    public stockas_klocka, swatch
    private
contains
    subroutine stockas_klocka()
        real(4) :: etime
        character(BUFLEN) :: fdate

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

        real(4) :: dtime
        
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


! vim:et:sw=4:ts=4
