!
!  Tunable parameters
!
module settings_module
    use kinds_module
    implicit none

    integer(I4) :: debug_level=0

    type settings_t
        character(256) :: title
        character(128) :: debug_file
        logical :: dryrun
        logical :: show_axis
        character(16) :: calc
        character(128) :: basis
        character(128) :: xdens
        character(128) :: mofile
        character(128) :: density
        real(DP), dimension(3) :: magnet
        character(2) :: magnet_axis
        logical :: is_uhf=.false.
        logical :: is_mpirun=.false.
        logical :: use_spherical=.false.
        logical :: use_giao=.true.
        logical :: use_diamag=.true.
        logical :: use_paramag=.true.
        logical :: acid=.true.
        logical :: jmod=.true.
        real(DP) :: screen_thrs
        logical :: use_screening
        integer(I4), dimension(2) :: morange
        integer(I4) :: lip_order
    end type

    type(settings_t), save :: settings
end module

! vim:et:sw=4:ts=4
