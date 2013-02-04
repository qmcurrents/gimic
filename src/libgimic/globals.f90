!
! Global parameters, types and routines
!
module globals_module
    use kinds_module
    use settings_module
    use teletype_module
    use getkw_class
    implicit none 

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
    integer, parameter :: JVECFD=30

    
    integer, parameter :: MAX_L=5 ! Max angular momentum in basis (h)
    integer, parameter :: BUFLEN=80
    integer, parameter :: MAX_LINE_LEN=200

    integer(I4), parameter :: CDENS_TAG=10
    integer(I4), parameter :: INTGRL_TAG=11
    integer(I4), parameter :: DIVJ_TAG=12
    integer(I4), parameter :: EDENS_TAG=13

    integer(I4), parameter :: SPIN_A=1, SPIN_B=2

    real(DP), parameter :: PII=3.141592653589793D0
    !real(DP), parameter :: PI=atan(1.d0)*4.d0
    real(DP), parameter :: PII_FAC=6.349363593424098D-002 ! (0.5/PI)**0.75
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
    real(DP), parameter :: SCREEN_THRS = 1.d-6

    real(DP), dimension(3), parameter :: INITRV = 123456789.d5
    
    real(DP), parameter :: D0=0.0D0,D1=1.0D0,D2=2.0D0,D4=4.0D0,D5=5.0D0
    real(DP), parameter :: DP25=0.25D0, DP50=0.50D0, DP75=0.75D0
    real(DP), parameter :: DP33=0.3333333D0 
    real(DP), dimension(3), parameter :: NILL_VECTOR=(/0.d0, 0.d0, 0.d0/)

    
    ! default filenames
    character(*), parameter :: DEFAULT_INPUT='gimic.inp'
    character(*), parameter :: DEFAULT_DENSFILE='XDENS'
    character(*), parameter :: DEFAULT_MOFILE='mos'

    character, dimension(0:7), parameter :: shell_names = &
        & (/'s','p','d','f','g','h','i','j'/)

    integer(4), dimension(100), private, save :: afd=0
    logical :: is_turbomole = .false.
    logical :: bert_is_evil = .false.
    integer :: mpi_world_size = 1
    integer :: mpi_rank = -1

    public au2si

    type cao2sao_t
        real(DP), dimension(:,:), pointer :: po
    end type

    type tensor_t
        sequence
        real(DP), dimension(3,3) :: t
    end type

    type vector_t
        real(DP), dimension(3) :: v
    end type

    type contraction_t
        integer(I4) :: npf,ncf,l,ncomp,nccomp ! n prims, n contrs, l-qnumber
        real(DP), dimension(:), pointer :: xp
        real(DP), dimension(:), pointer :: cc ! cc:s
        real(DP), dimension(:), pointer :: ncc ! norm cc:s
    end type

    type basis_t
        integer(I4) :: lmax, nctr, nshells, ncgto, ngto
        integer, dimension(MAX_L+1) :: nctrps ! contractions per shell
        type(contraction_t), dimension(:), pointer :: ctr ! s_n, p_n...
        real(DP), dimension(:), pointer :: thrs ! screening
        integer(I4), dimension(:), pointer :: pos ! starting idx for n:th ctr
    end type

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

    type molecule_t
        integer(I4) :: natoms
        type(atom_t), dimension(:), pointer :: atoms
        integer(I4) :: ngto, ncgto, nccgto
        type(cao2sao_t), pointer :: c2s
    end type

    interface copy_tensor
        module procedure copy_tensor_s
        module procedure copy_tensor_1d
        module procedure copy_tensor_2d
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

    function getfsize(fil) result(fs)
        character(*), intent(in) :: fil
        integer(I4) :: fs
        integer(I4), dimension(13) :: buf

        integer(I4) :: st

        call stat(fil, buf,st)
        fs=buf(8)
    end function

    function enumfile(string, rank) result(fname)
        character(*), intent(in) :: string
        integer(4), intent(in) :: rank
        character(100) :: fname

        write(fname,'(a,i0)') trim(string) // '.', rank
    end function

    function vcmp(v1, v2) result(p)
        real(DP), dimension(:), intent(in) :: v1, v2
        logical :: p

        real(DP), parameter :: thrs=1.d-12
        integer(I4) :: i, s1, s2

        p=.false.
        s1=size(v1)
        s2=size(v2)
        if (s1 /= s2) then
            call msg_error('Trying to compare vectors of different size!')
            return
        end if

        do i=1,s1
            if (abs(v1(i))-abs(v2(i)) > thrs) return
        end do
        p=.true.
    end function

    subroutine copy_tensor_s(src,dest)
        type(tensor_t), intent(in) :: src
        type(tensor_t), intent(out) :: dest

        dest%t=src%t
    end subroutine

    subroutine copy_tensor_1d(src,dest)
        type(tensor_t), dimension(:), intent(in) :: src
        type(tensor_t), dimension(:), intent(out) :: dest

        integer(I4) :: i
        integer(I4) :: ni, nj

        ni=size(src)
        nj=size(dest)

        if (ni /= nj) then
            call msg_error('copy_tensor_1d(): dimesnion mismatch!')
            stop
        end if

        do i=1,ni
            dest(i)%t=src(i)%t
        end do

    end subroutine

    subroutine copy_tensor_2d(src,dest)
        type(tensor_t), dimension(:,:), intent(in) :: src
        type(tensor_t), dimension(:,:), intent(out) :: dest

        integer(I4) :: i, j
        integer(I4), dimension(2) :: ni, nj
        ni=shape(src)
        nj=shape(dest)

        if (ni(1) /= nj(1) .or. ni(2) /= nj(2)) then
            call msg_error('copy_tensor_2d(): dimesnion mismatch!')
            stop
        end if

        do i=1,ni(1)
            do j=1,nj(2)
                dest(i,j)%t=src(i,j)%t
            end do
        end do

    end subroutine

    subroutine getfd(fd) 
        integer(4), intent(out) :: fd

        integer(4) :: i

        do i=81,100
            if (afd(i) == 0) then
                afd(i) = i
                fd=i
                return 
            end if
        end do
        write(str_g, '(a)') 'getfd(): File descriptor table full!'
        call msg_error(str_g)
        stop 1
    end subroutine

    subroutine freefd(fd)
        integer(4), intent(inout) :: fd

        if (fd == 0) return

        if (fd < 80 .or. fd > 100) then
            write(str_g, '(a,i4)') &
              'freefd(): File descriptor out of bounds: ', fd
            call msg_warn(str_g)
        else 
            if (afd(fd) /= 0) then
                afd(fd) = 0
                fd=0
            else
                write(str_g, '(a,i4)') &
                  'freefd(): Invalid file descriptor: ', fd
                call msg_warn(str_g)
            end if
        end if
    end subroutine

    subroutine closefd(fd)
        integer(4), intent(inout) :: fd

        integer(4) :: ffd

        if (fd == 0) return
        ffd=fd
        call freefd(fd)
        if (fd == 0) then
            close(ffd)
        end if

    end subroutine

    function au2si(au) result(si)
        real(DP), intent(in) :: au
        real(DP) :: si
        
        real(DP) :: aulength, auspeedoflight, speedoflight, aucharge, hbar
        real(DP) :: autime, autesla, audjdb

        aulength=0.52917726D-10
        auspeedoflight=137.03599D0
        speedoflight=299792458.D0
        aucharge=1.60217733D-19
        hbar=1.05457267D-34

        autime=aulength*auspeedoflight/speedoflight
        autesla=hbar/aucharge/aulength/aulength
        audjdb=aucharge/autime/autesla

        si=au*audjdb*1.d+09 ! nA/T

!        write(6,*) 'The obtained conversion factors'
!        write(6,*) autime,' au time in seconds'
!        write(6,*) autesla,' au magnetic field in tesla'
!        write(6,*) audjdb*1.D+09,' au induced current in nanoampere/tesla'
    end function 


end module


! vim:et:sw=4:ts=4
