!
!  Tunable parameters
!
module settings_m
    use kinds_m
    implicit none 
    
    integer(I4) :: debug_level=0

    logical :: is_turbomole=.false.
    logical :: is_mpirun=.false.
    logical :: is_uhf=.false.

    logical :: use_spherical=.false.
    logical :: use_giao=.true.
    logical :: use_diamag=.true.
    logical :: use_paramag=.true.

    integer :: mpi_world_size = 1
    integer :: mpi_rank = -1
end module

! vim:et:sw=4:ts=4
