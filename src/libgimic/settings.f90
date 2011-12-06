!
!  Tunable parameters
!
module settings_m
    use kinds_m
    implicit none 
    
    integer(I4) :: debug_level=0

    type main_settings_t
        character(256) :: title
        character(128) :: debug_file
        logical :: dryrun
        logical :: show_axis
        character(16) :: calc
        character(128) :: basis
        character(128) :: density
        real, dimension(3) :: magnet
        character(2) :: magnet_axis
        logical :: is_uhf=.false.
        logical :: is_mpirun=.false.
        logical :: use_spherical=.false.
        logical :: use_giao=.true.
        logical :: use_diamag=.true.
        logical :: use_paramag=.true.
        real :: screen_thrs
        logical :: use_screening
    end type

    type integral_settings_t
    end type
    
    type(main_settings_t), save :: main
    type(integral_settings_t), save :: integral_conf
end module

! vim:et:sw=4:ts=4
