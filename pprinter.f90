!
! Pretty printer module
!

module pprinter_m
	use kinds_m
	implicit none

	integer(SP), parameter :: STDOUT=6
	integer(SP), parameter :: STDERR=0
	integer(SP), parameter :: DEVNULL=-1
	integer(SP), parameter :: NONSTDOUT=7

	integer(SP), parameter :: PPBUFSZ=8096
	
	character(PPBUFSZ) :: ppbuf 

	integer(SP) :: pplun=6
	integer(SP) :: ppwrn=6
	integer(SP) :: pperr=6
	integer(SP) :: ppdbg=6
	integer(SP), private :: level=0

	interface xstr
		module procedure xstr_isp
		module procedure xstr_rdp
		module procedure xstr_isp_vec
		module procedure xstr_rdp_vec
		module procedure xstr_frm_isp
		module procedure xstr_frm_rdp
	end interface
	
	interface pprint
		module procedure pprint_str
		module procedure nl
	end interface
contains
	subroutine pprint_str(str)
		character(*), intent(in) :: str

		if (pplun == DEVNULL) return
		
		write(pplun, 100) trim(str)
100 format(x,a)	
	end subroutine

	subroutine pnote(str)
		character(*), intent(in) :: str

		if (pplun == DEVNULL) return

		write(pplun, 100) ' *** ', trim(str)
100 format(a,a)	
	end subroutine
	
	subroutine pinfo(str)
		character(*), intent(in) :: str

		if (pplun == DEVNULL) return

!        write(pplun, 100) ' +++ INFO: ', trim(str)
		write(pplun, 100) ' +++ ', trim(str)
100 format(a,a)	
	end subroutine

	subroutine pwarn(str)
		character(*), intent(in) :: str

		if (ppwrn == DEVNULL) return

!        write(pplun, 100) ' ### WARNING: ', trim(str), ' ###'
		write(ppwrn, 100) ' >>> ', trim(str) 
100 format(a,a)	
	end subroutine

	subroutine pdebug(str, l)
		character(*), intent(in) :: str
		integer(SP) :: l

		if (level == 0) return
		if (l > level) return
		if (ppdbg == DEVNULL) return

!        write(pplun, 100) ' DEBUG: ', trim(str), ' @@@'
		write(ppdbg, 100) ' @@@ ', trim(str), ' @@@'
100 format(a,a,a)	
	end subroutine

	subroutine perror(str)
		character(*), intent(in) :: str

		write(STDERR, 102) ' <<<  ERROR: ', trim(str), ' >>>'

100 format(a)	
101 format(x,a)	
102 format(a,a,a)	
	end subroutine

	subroutine pcritical(str)
		character(*), intent(in) :: str

		write(STDERR, 101) repeat('>', 70)
		write(STDERR, 100) ' >'
		write(STDERR, 102) ' > ', trim(str)
		write(STDERR, 100) ' >'
		write(STDERR, 101) repeat('>', 70)

100 format(a)	
101 format(x,a)	
102 format(a,a)	
	end subroutine

	subroutine nl
		if (pplun == DEVNULL) return

		write(pplun, *) 
	end subroutine

	subroutine set_pplun(u, chan)
		integer(SP), intent(in) :: u
		character(*), optional :: chan

		logical :: o_p
		integer(SP) :: tun
		character(32) :: ch
		
		tun=u
		ch='other'
		if (tun /= DEVNULL) then
			inquire(unit=tun,opened=o_p)
			if (.not.o_p) then
				write(ppbuf, '(a,i3,a)') 'Teletype unit ', tun,' not open'
				call pwarn(ppbuf)
				tun=STDOUT
			end if
		end if

		if (present(chan)) ch=chan

		select case(ch)
			case('warn')
				ppwrn=tun
			case('debug')
				ppdbg=tun
			case default
				pplun=tun
		end select
	end subroutine

	subroutine xpplun(u)
		integer(SP), intent(out) :: u
		u=pplun
	end subroutine 

	subroutine set_debug_level(l)
		integer(SP), intent(in) :: l
		level=l
	end subroutine 

	subroutine get_debug_level(l)
		integer(SP), intent(out) :: l
		l=level
	end subroutine 

	subroutine disable_stdout()
		logical :: o_p

		inquire(STDOUT,opened=o_p)
		if (o_p) then
			close(STDOUT)
		end if
		open(STDOUT, file='/dev/null')

		inquire(NONSTDOUT,opened=o_p)
		if (.not. o_p) then
			open(NONSTDOUT, file='/dev/stdout')
		end if

		if ( pplun == STDOUT ) then
			pplun=NONSTDOUT
		end if
	end subroutine
	
	subroutine enable_stdout()
		logical :: o_p

		inquire(NONSTDOUT,opened=o_p)
		if (o_p) then
			close(NONSTDOUT)
		end if

		inquire(STDOUT,opened=o_p)
		if (.not. o_p) then
			open(STDOUT, file='/dev/stdout')
		end if

		if ( pplun == NONSTDOUT ) then
			pplun=STDOUT
		end if
	end subroutine

	function xstr_isp(arg) result(s)
		integer(SP), intent(in) :: arg
		character(PPBUFSZ) :: s

		write(s, *) arg
	end function

	function xstr_rdp(arg) result(s)
		real(DP), intent(in) :: arg
		character(PPBUFSZ) :: s

		write(s, *) arg
	end function

	function xstr_isp_vec(arg) result(s)
		integer(SP), dimension(:), intent(in) :: arg
		character(PPBUFSZ) :: s

		write(s, *) arg
	end function

	function xstr_rdp_vec(arg) result(s)
		real(DP), dimension(:), intent(in) :: arg
		character(PPBUFSZ) :: s

		write(s, *) arg
	end function

	function xstr_frm_isp(frm,arg) result(s)
		character(*) :: frm
		integer(SP), intent(in) :: arg
		character(PPBUFSZ) :: s

		write(s, frm) arg
	end function

	function xstr_frm_rdp(frm,arg) result(s)
		character(*) :: frm
		real(DP), intent(in) :: arg
		character(PPBUFSZ) ::  s

		write(s, frm) arg
	end function
end module
