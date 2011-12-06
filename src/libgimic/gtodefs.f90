!
! Definitions of the angular part for cartesian GTOs. 
! GTOMAP maps the set of {i,j,k} to the correspnding index in the GTx arrays.
!

module gtodefs_m
    use globals_m
    use settings_m
    implicit none

    type gto_jar
        integer(I4) :: l, ncomp
        real(DP), dimension(:,:), pointer :: nlm
    end type
    
    integer(I4), dimension(0:MAX_L,0:MAX_L,0:MAX_L) :: GTO_MAP

    integer, parameter, public :: X_  = 1, Y_  = 2, Z_  = 3
    integer, parameter, public :: XX_ = 1, XY_ = 2, XZ_ = 3
    integer, parameter, public :: YY_ = 4, YZ_ = 5, ZZ_ = 6
    integer, parameter, public :: XXX_=1, XXY_=2, XXZ_=3, XYY_=4, XYZ_=5
    integer, parameter, public :: YYY_=6, YYZ_=7, ZZX_=8, ZZY_=9, ZZZ_=10 
    
    real(DP), dimension(3), parameter :: S_0 = (/0.0,0.0,0.0/)
    
    real(DP), dimension(3), parameter :: P_X = (/1.0,0.0,0.0/)
    real(DP), dimension(3), parameter :: P_Y = (/0.0,1.0,0.0/)
    real(DP), dimension(3), parameter :: P_Z = (/0.0,0.0,1.0/)
    
    real(DP), dimension(3), parameter :: D_XX = (/2.0,0.0,0.0/)
    real(DP), dimension(3), parameter :: D_XY = (/1.0,1.0,0.0/)
    real(DP), dimension(3), parameter :: D_XZ = (/1.0,0.0,1.0/)
    real(DP), dimension(3), parameter :: D_YY = (/0.0,2.0,0.0/)
    real(DP), dimension(3), parameter :: D_YZ = (/0.0,1.0,1.0/)
    real(DP), dimension(3), parameter :: D_ZZ = (/0.0,0.0,2.0/)
    
    real(DP), dimension(3), parameter :: F_XXX = (/3.0,0.0,0.0/)
    real(DP), dimension(3), parameter :: F_XXY = (/2.0,1.0,0.0/)
    real(DP), dimension(3), parameter :: F_XXZ = (/2.0,0.0,1.0/)
    real(DP), dimension(3), parameter :: F_XYY = (/1.0,2.0,0.0/)
    real(DP), dimension(3), parameter :: F_XYZ = (/1.0,1.0,1.0/)
    real(DP), dimension(3), parameter :: F_XZZ = (/1.0,0.0,2.0/)
    real(DP), dimension(3), parameter :: F_YYY = (/0.0,3.0,0.0/)
    real(DP), dimension(3), parameter :: F_YYZ = (/0.0,2.0,1.0/)
    real(DP), dimension(3), parameter :: F_YZZ = (/0.0,1.0,2.0/)
    real(DP), dimension(3), parameter :: F_ZZZ = (/0.0,0.0,3.0/)

    real(DP), dimension(3), parameter :: G_400 = (/4.0, 0.0, 0.0/)
    real(DP), dimension(3), parameter :: G_310 = (/3.0, 1.0, 0.0/)
    real(DP), dimension(3), parameter :: G_301 = (/3.0, 0.0, 1.0/)
    real(DP), dimension(3), parameter :: G_220 = (/2.0, 2.0, 0.0/)
    real(DP), dimension(3), parameter :: G_211 = (/2.0, 1.0, 1.0/)
    real(DP), dimension(3), parameter :: G_202 = (/2.0, 0.0, 2.0/)
    real(DP), dimension(3), parameter :: G_130 = (/1.0, 3.0, 0.0/)
    real(DP), dimension(3), parameter :: G_121 = (/1.0, 2.0, 1.0/)
    real(DP), dimension(3), parameter :: G_112 = (/1.0, 1.0, 2.0/)
    real(DP), dimension(3), parameter :: G_103 = (/1.0, 0.0, 3.0/)
    real(DP), dimension(3), parameter :: G_040 = (/0.0, 4.0, 0.0/)
    real(DP), dimension(3), parameter :: G_031 = (/0.0, 3.0, 1.0/)
    real(DP), dimension(3), parameter :: G_022 = (/0.0, 2.0, 2.0/)
    real(DP), dimension(3), parameter :: G_013 = (/0.0, 1.0, 3.0/)
    real(DP), dimension(3), parameter :: G_004 = (/0.0, 0.0, 4.0/)

    real(DP), dimension(3), parameter :: H_500 = (/5.0, 0.0, 0.0/)
    real(DP), dimension(3), parameter :: H_410 = (/4.0, 1.0, 0.0/)
    real(DP), dimension(3), parameter :: H_401 = (/4.0, 0.0, 1.0/)
    real(DP), dimension(3), parameter :: H_320 = (/3.0, 2.0, 0.0/)
    real(DP), dimension(3), parameter :: H_311 = (/3.0, 1.0, 1.0/)
    real(DP), dimension(3), parameter :: H_302 = (/3.0, 0.0, 2.0/)
    real(DP), dimension(3), parameter :: H_230 = (/2.0, 3.0, 0.0/)
    real(DP), dimension(3), parameter :: H_221 = (/2.0, 2.0, 1.0/)
    real(DP), dimension(3), parameter :: H_212 = (/2.0, 1.0, 2.0/)
    real(DP), dimension(3), parameter :: H_203 = (/2.0, 0.0, 3.0/)
    real(DP), dimension(3), parameter :: H_140 = (/1.0, 4.0, 0.0/)
    real(DP), dimension(3), parameter :: H_131 = (/1.0, 3.0, 1.0/)
    real(DP), dimension(3), parameter :: H_122 = (/1.0, 2.0, 2.0/)
    real(DP), dimension(3), parameter :: H_113 = (/1.0, 1.0, 3.0/)
    real(DP), dimension(3), parameter :: H_104 = (/1.0, 0.0, 4.0/)
    real(DP), dimension(3), parameter :: H_050 = (/0.0, 5.0, 0.0/)
    real(DP), dimension(3), parameter :: H_041 = (/0.0, 4.0, 1.0/)
    real(DP), dimension(3), parameter :: H_032 = (/0.0, 3.0, 2.0/)
    real(DP), dimension(3), parameter :: H_023 = (/0.0, 2.0, 3.0/)
    real(DP), dimension(3), parameter :: H_014 = (/0.0, 1.0, 4.0/)
    real(DP), dimension(3), parameter :: H_005 = (/0.0, 0.0, 5.0/)
     
    real(DP), dimension(3,1), target :: GTS = reshape & 
        ((/S_0/),(/3,1/))
    
    real(DP), dimension(3,3), target :: GTP = reshape &
        ((/P_X, P_Y, P_Z/), (/3,3/))
        
    real(DP), dimension(3,6), target :: GTD = reshape &
        ((/D_XX, D_XY, D_XZ, D_YY, D_YZ, D_ZZ /), (/3,6/))
        
    real(DP), dimension(3,10), target :: GTF = reshape &
        ((/F_XXX, F_XXY, F_XXZ, F_XYY, F_XYZ, &
        F_XZZ, F_YYY, F_YYZ, F_YZZ, F_ZZZ/), (/3,10/))

    real(DP), dimension(3,15), target :: GTG = reshape &
        ((/G_400, G_310, G_301, G_220, G_211, G_202, G_130, G_121, &
        G_112, G_103, G_040, G_031, G_022, G_013, G_004/), (/3,15/))

    real(DP), dimension(3,21), target :: GTH = reshape &
        ((/H_500, H_410, H_401, H_320, H_311, H_302, H_230, H_221, H_212, &
           H_203, H_140, H_131, H_122, H_113, H_104, H_050, H_041, H_032, &
           H_023, H_014, H_005 /),(/3,21/))
    
    ! TURBOMOLE ordering... (in symlib/pdfgtr.f)
    real(DP), dimension(3,6), target :: GTD_TM = reshape &
        ((/D_XX, D_YY, D_ZZ, D_XY, D_XZ, D_YZ /), (/3,6/))
        
    real(DP), dimension(3,10), target :: GTF_TM = reshape &
        ((/F_XXX, F_YYY, F_ZZZ, F_XXY, F_XXZ, &
        F_XYY, F_YYZ, F_XZZ, F_YZZ, F_XYZ/), (/3,10/))

    real(DP), dimension(3,15), target :: GTG_TM = reshape &
        ((/G_400, G_040, G_004, G_310, G_301, G_130, G_031, G_103, &
        G_013, G_220, G_202, G_022, G_211, G_121, G_112/), (/3,15/))

    real(DP), dimension(3,21), target :: GTH_TM = reshape &
        ((/H_500, H_050, H_005, H_410, H_401, H_140, H_041, H_104, H_014, &
           H_320, H_302, H_230, H_032, H_203, H_023, H_311, H_131, H_113, &
           H_221, H_212, H_122 /),(/3,21/))

    type(gto_jar), dimension(0:MAX_L):: GTO_DEFS
    
    public setup_gtos, get_gto_nlm, gtomap
    private
contains
    subroutine setup_gtos()
        integer(I4) :: a, i, k, l, m
        integer(I4) :: turbo
        
        GTO_DEFS(0)%l=0; GTO_DEFS(0)%ncomp=1   
        GTO_DEFS(1)%l=1; GTO_DEFS(1)%ncomp=3   
        GTO_DEFS(2)%l=2; GTO_DEFS(2)%ncomp=6   
        GTO_DEFS(3)%l=3; GTO_DEFS(3)%ncomp=10  
        GTO_DEFS(4)%l=4; GTO_DEFS(4)%ncomp=15  
        GTO_DEFS(5)%l=5; GTO_DEFS(5)%ncomp=21  

        GTO_DEFS(0)%nlm=>GTS
        GTO_DEFS(1)%nlm=>GTP
        
        if (is_turbomole) then
            GTO_DEFS(2)%nlm=>GTD_TM
            GTO_DEFS(3)%nlm=>GTF_TM
            GTO_DEFS(4)%nlm=>GTG_TM
            GTO_DEFS(5)%nlm=>GTH_TM
        else
            GTO_DEFS(2)%nlm=>GTD
            GTO_DEFS(3)%nlm=>GTF
            GTO_DEFS(4)%nlm=>GTG
            GTO_DEFS(5)%nlm=>GTH
        end if
        
        do a=0,MAX_L
            do i=1,GTO_DEFS(a)%ncomp
                k=int(GTO_DEFS(a)%nlm(1,i))
                l=int(GTO_DEFS(a)%nlm(2,i))
                m=int(GTO_DEFS(a)%nlm(3,i))
                GTO_MAP(k,l,m)=i
            end do
        end do
    end subroutine

    subroutine get_gto_nlm(l,g) 
        integer(I4), intent(in) :: l
        real(DP), dimension(:,:), pointer :: g

        g=>GTO_DEFS(l)%nlm
    end subroutine

    function gtomap(i,j,k) result(n)
        integer(I4), intent(in) :: i,j,k
        integer(I4) :: n

        n=GTO_MAP(i,j,k)
    end function
end module	

! vim:et:sw=4:ts=4
