
module tensor_m
    use kinds_m
    implicit none
    
    interface operator(.tp.)
        module procedure d_tprod 
    end interface
    
    interface operator(.tp.)
        module procedure i_tprod 
    end interface

    interface cross_product
        module procedure d_cross_product, i_cross_product
    end interface
contains
    function d_tprod(A,B) result(C)
        real(DP), dimension(:), intent(in) :: A, B
        real(DP), dimension(size(A),size(B)) :: C
        
        integer(I4) :: i, sz
        sz=size(A)
        do i=1,sz
            C(:,i)=A(i)*B
        end do
    end function

    function i_tprod(A,B) result(C)
        integer(I4), dimension(:), intent(in) :: A, B
        integer(I4), dimension(size(A),size(B)) :: C
        
        integer(I4) :: i, sz
        sz=size(A)
        do i=1,sz
            C(:,i)=A(i)*B
        end do
    end function

    function d_cross_product(a,b) result(c)
        real(DP), dimension(:), intent(in) :: a,b
        real(DP), dimension(3) :: c

        c(1)=a(2)*b(3)-a(3)*b(2)
        c(2)=a(3)*b(1)-a(1)*b(3)
        c(3)=a(1)*b(2)-a(2)*b(1)
    end function

    function i_cross_product(a,b) result(c)
        integer(8), dimension(:), intent(in) :: a,b
        integer(8), dimension(3) :: c

        c(1)=a(2)*b(3)-a(3)*b(2)
        c(2)=a(1)*b(3)-a(3)*b(1)
        c(3)=a(1)*b(2)-a(2)*b(1)
    end function

end module
