!
! Write VTK ImageData XML files
!
module vtkplot_module
    use globals_module
    use settings_module
    use grid_class
    implicit none

contains
    subroutine write_vtk_imagedata(fname, grid, pdata)
        character(*), intent(in) :: fname
        type(grid_t) :: grid
        real(DP), dimension(:,:,:) :: pdata

        integer(I4) :: skp=1
        integer(I4) :: i,j,k,l, fd
        integer(I4), dimension(3) :: npts
        real(DP), dimension(2) :: qrange
        real(DP), dimension(3) :: qmin, qmax, step

        call getfd(fd)
        open(fd, file=trim(fname), form='formatted', status='unknown')

        call get_grid_size(grid, npts(1), npts(2), npts(3))
        qmin=gridpoint(grid,1,1,1)
        qmax=gridpoint(grid,npts(1),npts(2),npts(3))
        step=(qmax-qmin)/(npts-1)

        write(fd, '(a)') '<?xml version="1.0"?>'
        write(fd, *) '<VTKFile type="ImageData" ', &
            'version="0.1" byte_order="LittleEndian">'
        write(fd, *) '  <ImageData WholeExtent="', &
            0, npts(1)-1, &
            0, npts(2)-1, &
            0, npts(3)-1, &
            '" Origin="', qmin(1), qmin(2),qmin(3), '" Spacing="', &
            step(1), step(2), step(3), '">'
        write(fd, *) '  <Piece Extent="', &
            0, npts(1)-1, &
            0, npts(2)-1, &
            0, npts(3)-1, &
            '">'
        write(fd, *) '  <PointData Scalars="scalars">'
        write(fd, *) '  <DataArray Name="scalars" type="Float32" ', &
            'NumberOfComponents="1" Format="ascii">'

        l=0
        do i=1,npts(1)
            do j=1,npts(2)
                do k=1,npts(3)
                    write(fd,'(e14.6)',advance='no') pdata(i,j,k)
                    if (mod(l,4) == 0) write(fd,*)
                    l=l+1
                end do
            end do
        end do
        write(fd, *)
        write(fd, *) '   </DataArray>'
        write(fd, *) '   </PointData>'

!        write(fd, *) '   <CellData Scalaras="foo">'
!        write(fd, *) '   </CellData>'

        write(fd, *) '   </Piece>'
        write(fd, *) '   </ImageData>'
        write(fd, *) '</VTKFile>'

        call closefd(fd)

    end subroutine
end module

! vim:et:sw=4:ts=4
