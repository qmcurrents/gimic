module acid_module
     ! H. Fliegl, University of Oslo (2013)
     use globals_module
     use teletype_module
     use tensor_module
     !use grid_class
     !use magnet_module
     implicit none

contains
     function get_acid(r, tens) result(dT2)
     ! purpose: do acid summation of T tensor
     ! see: R. Herges and D. Gruenich,
     !      J. Chem. Phys. A, 105, 3214 (2001)
     !
     !     t_xx t_xy t_xz   1 4 7
     ! T = t_yx t_yy t_yz   2 5 8
     !     t_zx t_zy t_zz   3 6 9
     !
     ! Delta T^2 =
     !   1/3 [ (t_xx - tyy)^2 + (t_yy - t_zz)^2 + (t_zz - t_xx)^2 ]
     ! + 1/2 [ (t_xy + t_yx)^2 + (t_xz + t_zx)^2 + (t_yz + t_zy)^2 ]
     !
     ! assume: T tensor is written on tens(9:number of points)
     !         first column of T eq. first entries of tens(1:3,npts)
     !
     real(DP), dimension(3), intent(in) :: r
     real(DP), dimension(9), intent(in) :: tens
     real(DP) :: dT2, tmp
     real(DP) :: xxmyy, yymzz, zzmxx, xypyx, xzpzx, yzpzy

     xxmyy = (tens(1) - tens(5))*(tens(1) - tens(5))
     yymzz = (tens(5) - tens(9))*(tens(5) - tens(9))
     zzmxx = (tens(9) - tens(1))*(tens(9) - tens(1))

     xypyx = (tens(4) + tens(2))*(tens(4) + tens(2))
     xzpzx = (tens(7) + tens(3))*(tens(7) + tens(3))
     yzpzy = (tens(8) + tens(6))*(tens(8) + tens(6))

     tmp = DP33*(xxmyy + yymzz + zzmxx) + DP50*(xypyx + xzpzx + yzpzy)
     ! conversion into SI units nA/T --> done later
     ! dT2 = au2si(tmp)
     dT2 = tmp
     end function

end module


! vim:et:sw=4:ts=4
