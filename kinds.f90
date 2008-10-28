module kinds_m
	implicit none
	integer, parameter :: MP = kind(0)
	integer, parameter :: SP = kind(0.0e0)
	integer, parameter :: DP = kind(0.0d0)
#ifdef HAVE_QP
 	integer, parameter :: QP = kind(0.0q0)
	integer, parameter :: QP = 16
#endif
	integer, parameter :: I2 = 2
	integer, parameter :: I4 = 4
	integer, parameter :: I8 = 8
end module 
