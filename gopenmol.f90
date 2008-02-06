module
	use globals_m
	use grid_class
	implicit none

contains	 
	subroutine gom_plot_vec(vv, grid, fname)
		type(vector_t), dimension(:,:) :: vv
		type(grid_t) :: grid
		character(*) :: fname

		integer(I4) :: surface, rank, p1, p2, p3
		integer(I4) :: i, j, k, l
		real(SP), dimension(3) :: qmin, qmax
		real(DP), dimension(3) :: norm

		vv
		open(GOPFD,file=trim(fname),access='direct',recl=4)

		surface=200
		rank=3

		call get_grid_size(grid, p1, p2, p3)
		norm=get_grid_normal(grid)
		qmin=gridpoint(grid,1,1,1)*AU2A
		qmax=gridpoint(grid,p1,p2,p3)*AU2A

		write(GOPFD,rec=1) rank;      write(GOPFD2,rec=1) rank
		write(GOPFD,rec=2) surface;   write(GOPFD2,rec=2) surface
		write(GOPFD,rec=3) p3;        write(GOPFD2,rec=3) p3
		write(GOPFD,rec=4) p2;        write(GOPFD2,rec=4) p2
		write(GOPFD,rec=5) p1;        write(GOPFD2,rec=5) p1
		write(GOPFD,rec=6) qmin(3);   write(GOPFD2,rec=6) qmin(3)
		write(GOPFD,rec=7) qmax(3);   write(GOPFD2,rec=7) qmax(3)
		write(GOPFD,rec=8) qmin(2);   write(GOPFD2,rec=8) qmin(2)
		write(GOPFD,rec=9) qmax(2);   write(GOPFD2,rec=9) qmax(2)
		write(GOPFD,rec=10) qmin(1);  write(GOPFD2,rec=10) qmin(1)
		write(GOPFD,rec=11) qmax(1);  write(GOPFD2,rec=11) qmax(1)

		l=12
		do k=1,p3
!            call jvec_io(jf, k, 'r')
			do j=1,p2
				do i=1,p1
!                    write(GOPFD,rec=l) real(sqrt(sum(vv(i,j)%v**2)))
					l=l+1
				end do
			end do
		end do

		close(GOPFD)
		close(GOPFD2)
	end subroutine
end module
