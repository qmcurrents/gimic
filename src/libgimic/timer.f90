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

        real(4) :: hours_per_sec = 1.0/3600.0
        real(4) :: days_per_hour = 1.0/24.0
        real(4) :: hours_per_day = 24.0

        delta_t=etime(times)
        call msg_out(repeat('-', 70))
        write(str_g, '(a,f9.2,a,a,f6.1,a)') '   wall time:', delta_t, 'sec', &
            ' (', delta_t*hours_per_sec, ' h )'
        if ( delta_t*hours_per_sec .gt. 2*hours_per_day ) then
            write(str_g, '(a,f4.1,a)') '(', delta_t*hours_per_sec*days_per_hour, ' days )' 
        end if
        call msg_out(str_g)
        write(str_g, '(a,f9.2,a,a,f6.1,a)') '        user:', times(1), 'sec', &
            ' (', times(1)*hours_per_sec, ' h )'
        if ( times(1)*hours_per_sec .gt. 2*hours_per_day ) then
            write(str_g, '(a,f4.1,a)') '(', times(1)*hours_per_sec*days_per_hour, ' days )' 
        end if
        call msg_out(str_g)
        write(str_g, '(a,f9.2,a,a,f6.1,a)') '         sys:', times(2), 'sec', &
            ' (', times(2)*hours_per_sec, ' h )'
        if ( times(2)*hours_per_sec .gt. 2*hours_per_day ) then
            write(str_g, '(a,f4.1,a)') '(', times(2)*hours_per_sec*days_per_hour, ' days )' 
        end if
        call msg_out(str_g)
        call msg_out(repeat('-', 70))
        call msg_out(fdate())
    end subroutine

    subroutine swatch(foo)
        integer(I4), intent(in) :: foo

        real(4) :: dtime
        real(4) :: hours_per_sec = 1.0/3600.0
        real(4) :: days_per_hour = 1.0/24.0
        real(4) :: hours_per_day = 24.0

        delta_t=dtime(times)
        if (foo > 0) then
            call msg_out(repeat('-', 70))
        write(str_g, '(a,f9.2,a,a,f6.1,a)') '   wall time:', delta_t, 'sec', &
            ' (', delta_t*hours_per_sec, ' h )'
        if ( delta_t*hours_per_sec .gt. 2*hours_per_day ) then
            write(str_g, '(a,f4.1,a)') '(', delta_t*hours_per_sec*days_per_hour, ' days )' 
        end if
        call msg_out(str_g)
        write(str_g, '(a,f9.2,a,a,f6.1,a)') '        user:', times(1), 'sec', &
            ' (', times(1)*hours_per_sec, ' h )'
        if ( times(1)*hours_per_sec .gt. 2*hours_per_day ) then
            write(str_g, '(a,f4.1,a)') '(', times(1)*hours_per_sec*days_per_hour, ' days )' 
        end if
        call msg_out(str_g)
        write(str_g, '(a,f9.2,a,a,f6.1,a)') '         sys:', times(2), 'sec', &
            ' (', times(2)*hours_per_sec, ' h )'
        if ( times(2)*hours_per_sec .gt. 2*hours_per_day ) then
            write(str_g, '(a,f4.1,a)') '(', times(2)*hours_per_sec*days_per_hour, ' days )' 
        end if
        call msg_out(str_g)
            call msg_out(repeat('-', 70))
        end if
    end subroutine
end module


! vim:et:sw=4:ts=4
