!
! $Id$
!
!include 'string_m.h'
module globals_m
	use precision_m
	use teletype_m
	use getkw_class
	implicit none 
	!include 'getkwf.h'

	! file descriptors
	integer, parameter :: BASFD=42
	integer, parameter :: XDFD=43
	integer, parameter :: MODFD=45
	integer, parameter :: JVPFD=46
	integer, parameter :: NJVFD=47
	integer, parameter :: JPRJFD=49
	integer, parameter :: DIVJFD=50
	integer, parameter :: DJPFD=51
	integer, parameter :: JINTFD=52
	integer, parameter :: EDFD=53
	integer, parameter :: EDPFD=54
	integer, parameter :: GOPFD=55
	integer, parameter :: GOPFD2=56
	integer, parameter :: GRIDFD=57
	integer, parameter :: JTFD=70
	integer, parameter :: JVECFD=80

	
	integer, parameter :: MAX_L=5 ! Max angular momentum in basis (h)
	integer, parameter :: BUFLEN=80
	integer, parameter :: MAX_LINE_LEN=200

	integer(I4), parameter :: CDENS_TAG=10
	integer(I4), parameter :: INTGRL_TAG=11
	integer(I4), parameter :: DIVJ_TAG=12
	integer(I4), parameter :: EDENS_TAG=13

	integer(I4), parameter :: SPIN_A=1, SPIN_B=2

	real(DP), parameter :: PI=3.141592653589793D0
	!real(DP), parameter :: PI=atan(1.d0)*4.d0
	real(DP), parameter :: PI_FAC=6.349363593424098D-002 ! (0.5/PI)**0.75
	real(DP), parameter :: C_AU=137.035987 !137.0359895
	real(DP), parameter :: ALPHA=1.d0/C_AU !0.007297351
	real(DP), parameter :: ALPHA2=ALPHA**2
	real(DP), parameter :: ZETA=0.5d0

	real(DP), parameter :: AU2M = 0.52917726e-10
	real(DP), parameter :: AU2NM = 0.52917726e-1
	real(DP), parameter :: AU2A = 0.52917726
	real(DP), parameter :: AU2PM = 52.917726
	real(DP), parameter :: A2AU = 1.889725949
	real(DP), parameter :: NM2AU = 1.889725949e+1
	real(DP), parameter :: PM2AU = 1.889725949e+3

	real(DP), dimension(3), parameter :: INITRV = 123456789.d5
	
	real(DP), parameter :: D0=0.0D0,D1=1.0D0,D2=2.0D0,D4=4.0D0,D5=5.0D0
	real(DP), parameter :: DP25=0.25D0, DP50=0.50D0, DP75=0.75D0
	
	! default filenames
	character(*), parameter :: DEFAULT_INPUT='gimic.inp'
	character(*), parameter :: DEFAULT_DENSFILE='XDENS'
	character(*), parameter :: DEFAULT_MOFILE='mos'

	logical :: turbomole_p=.false.
	logical :: spherical=.false.
    logical :: bert_is_evil=.false.
	logical :: mpirun_p=.false.
	logical :: master_p=.false.
	logical :: uhf_p=.false.

	type(getkw_t) :: input

	character, dimension(0:7), parameter :: shell_names = &
    	& (/'s','p','d','f','g','h','i','j'/)

! The atom_t type holds _individual_ atom data. atom%basis is therefore a 
! pointer to a basis_t type, so that we don't have to duplicate the basis 
! for every atom, if it is the same.
	type atom_t
		character(2) :: symbol
		character(2) :: id     
		real(DP), dimension(3) :: coord
		real(DP) :: charge
		type(basis_t), pointer :: basis  
	end type

	type basis_t
		integer(I4) :: lmax, nctr, nshells, ncgto, ngto
		integer, dimension(MAX_L+1) :: nctrps ! contractions per shell
		type(contraction_t), dimension(:), pointer :: ctr ! s_n, p_n...
	end type

	type contraction_t
		integer(I4) :: npf,ncf,l,ncomp,nccomp ! n prims, n contrs, l-qnumber
		real(DP), dimension(:), pointer :: xp
		real(DP), dimension(:), pointer :: cc ! cc:s
		real(DP), dimension(:), pointer :: ncc ! norm cc:s
	end type

	type cao2sao_t
		real(DP), dimension(:,:), pointer :: po
	end type

	type tensor_t
		real(DP), dimension(3,3) :: t
	end type

	type vector_t
		real(DP), dimension(3) :: v
	end type

	interface xchar
		module procedure xchar_i, xchar_d
	end interface

contains
    function xtrim(str) result(r)
        character(BUFLEN), intent(in) :: str
        character(BUFLEN) :: r
        
        integer(I4) :: i
        character(1), dimension(BUFLEN) :: tmp

        tmp=transfer(str,tmp)
        
        do i=1,BUFLEN
            if (iachar(tmp(i)) == 0) then
                tmp(i) = ' '
                exit
            end if
        end do
        r=transfer(tmp,r)
        
    end function

	function getnlines(fd) result(n)
		integer(I4), intent(in) :: fd
		integer(I4) :: n
		
		n=0
		rewind(fd)
		do 
			read(fd,*,end=100)
			n=n+1
		end do
100		rewind(fd)
	end function

	function getfsize(fd) result(fs)
		integer(I4) :: fd, fs

		integer(I4) :: i, fseek, ftell

		i=fseek(fd,0,2)
		fs=ftell(fd)
		rewind(fd)
	end function

	function xchar_i(i) result(s)
		integer(I4) :: i
		character(20) :: s

		write(s, '(i)') i
		s=adjustl(s)
	end function

	function xchar_d(d) result(s)
		real(DP) :: d
		character(20) :: s

		write(s, '(f)') d
		s=adjustl(s)
	end function

	function hostname() result(hn)
		character(BUFLEN) :: hn
		external hostnm
		integer(I4) :: hostnm, ierr

		ierr=hostnm(hn)
	end function
end module
