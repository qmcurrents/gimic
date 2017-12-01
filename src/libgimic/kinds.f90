module kinds_module
    implicit none
    integer, parameter :: MP = kind(0)
    integer, parameter :: SP = kind(0.0e0)
    integer, parameter :: DP = kind(0.0d0)

    integer, parameter :: I2 = 2
    integer, parameter :: I4 = 4
    integer, parameter :: I8 = 8
end module

! vim:et:sw=4:ts=4
