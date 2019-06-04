!
! Write VTK ImageData XML files
!
module vtkplot_module
    use globals_module
    use settings_module
    use grid_class
    use ISO_FORTRAN_ENV

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
!        step=(qmax-qmin)/(npts-1) ! In a 2D plot qmin(3)-qmax(3) = 0 -> dividing it by npts(3) gives a NaN
        do i=1,3
           step(i) = qmax(i) - qmin(i)
           if ( step(i) > 1E-8 ) then    ! floating point numbers cannot be compared as ( step(i) == 0.0 )
               step(i) = step(i) / ( npts(i) - 1)  ! if ( step(i) != 0 ) divide, else leave it 0
           end if
        end do

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
        write(fd, *) '  <DataArray Name="scalars" type="Float64" ', &
            'NumberOfComponents="1" Format="ascii">'

        l=0
        !do i=1,npts(1)
        !    do j=1,npts(2)
        !        do k=1,npts(3)
        ! otherwise discrepancy with geometry when visualizing
        ! both - density and molecular geometry
        do k=1,npts(3)
            do j=1,npts(2)
                do i=1,npts(1)
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

    subroutine write_vtk_vector_imagedata(fname, grid, pdata)
        character(*), intent(in) :: fname
        type(grid_t) :: grid
        real(DP), dimension(:,:,:,:) :: pdata

        integer(I4) :: skp=1
        integer(I4) :: i,j,k,l, fd, idx
        integer(I4), dimension(3) :: npts
        real(DP), dimension(2) :: qrange
        real(DP), dimension(3) :: qmin, qmax, step
        real(DP) :: avg_x, avg_y, avg_z
        real(DP) :: norm

        call getfd(fd)
        open(fd, file=trim(fname), form='formatted', status='unknown')

        call get_grid_size(grid, npts(1), npts(2), npts(3))
        qmin=gridpoint(grid,1,1,1)
        qmax=gridpoint(grid,npts(1),npts(2),npts(3))
!        step=(qmax-qmin)/(npts-1) ! In a 2D plot qmin(3)-qmax(3) = 0 -> dividing it by npts(3) gives a NaN
        do i=1,3
           step(i) = qmax(i) - qmin(i)
           if ( step(i) > 1E-8 ) then    ! floating point numbers cannot be compared as ( step(i) == 0.0 )
               step(i) = step(i) / ( npts(i) - 1)  ! if ( step(i) != 0 ) divide, else leave it 0
           end if
        end do

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
            ! structure of vtk image data file
            ! <VTKFile type=" ImageData" ...>
            !<ImageData WholeExtent=" x1 x2 y1 y2 z1 z2"
            !Origin=" x0 y0 z0" Spacing=" dx dy dz">
            !<Piece Extent=" x1 x2 y1 y2 z1 z2">
            !<PointData>...</ PointData>
            !<CellData>...</ CellData>
            !</ Piece>
            !</ ImageData>
            !</ VTKFile>


        write(fd, *) '  <PointData Scalars="scalars">'
        ! data name data type nocomp format
        write(fd, *) '  <DataArray Name="vectors" type="Float64" ', &
            'NumberOfComponents="3" Format="ascii">'

        l=0
        ! note there is another loop structure for jvec!!!!
        !do i=1,npts(1)
        !    do j=1,npts(2)
        !        do k=1,npts(3)
        do k=1,npts(3)
            do j=1,npts(2)
                do i=1,npts(1)
                    !write(fd,'(e14.6)',advance='no') pdata(i,j,k,1:3)
                    !do idx = 1, 3
                      !write(fd,'(e14.6)', advance='no') pdata(i,j,k,idx)
                      !write(fd,'(e14.6)',advance='no') pdata(i,j,k,idx)
                      !if (mod(l,4) == 0) write(fd,*)
                      !l=l+1
                    !end do
                    !l=l+1
                    write(fd,'(3e14.6)') pdata(i,j,k,1), pdata(i,j,k,2), pdata(i,j,k,3)
                end do
            end do
        end do
        write(fd, *)
        write(fd, *) '   </DataArray>'
        write(fd, *) '   </PointData>'

        write(fd, *) '   <CellData Scalars="foo">'

        ! Assign the norm of the current vector to each cell.
        ! A cell is a volume which is---for ImageData---implicitly defined as little cubes and bounded by eight points,
        ! so we average the current vectors of the eight adjacent points, and take the norm of the result.
        ! In the special case of a cube of zero-thickness, cells are not volumes but rather areas, so we average over four vectors.
        ! https://lorensen.github.io/VTKExamples/site/VTKFileFormats/#imagedata
        if(npts(1) > 1 .and. npts(2) > 1 .and. npts(3) > 1) then
            do k=1,npts(3)-1
                do j=1,npts(2)-1
                    do i=1,npts(1)-1
                        avg_x = (pdata(i,j,k,1) + pdata(i+1,j,k,1) + pdata(i,j+1,k,1) + pdata(i,j,k+1,1) &
                              + pdata(i+1,j+1,k,1) + pdata(i+1,j,k+1,1) + pdata(i,j+1,k+1,1) + pdata(i+1,j+1,k+1,1))/8.0
                        avg_y = (pdata(i,j,k,2) + pdata(i+1,j,k,2) + pdata(i,j+1,k,2) + pdata(i,j,k+1,2) &
                              + pdata(i+1,j+1,k,2) + pdata(i+1,j,k+1,2) + pdata(i,j+1,k+1,2) + pdata(i+1,j+1,k+1,2))/8.0
                        avg_z = (pdata(i,j,k,3) + pdata(i+1,j,k,3) + pdata(i,j+1,k,3) + pdata(i,j,k+1,3) &
                              + pdata(i+1,j+1,k,3) + pdata(i+1,j,k+1,3) + pdata(i,j+1,k+1,3) + pdata(i+1,j+1,k+1,3))/8.0
                        norm = sqrt( avg_x**2 + avg_y**2 + avg_z**2 )
                        write(fd,'(e14.6)') norm
                    end do
                end do
            end do
        else if(npts(1) > 1 .and. npts(2) > 1 .and. npts(3) == 1) then
            k=1
            do j=1,npts(2)-1
                do i=1,npts(1)-1
                    avg_x = (pdata(i,j,k,1) + pdata(i+1,j,k,1) + pdata(i,j+1,k,1) + pdata(i+1,j+1,k,1))/4.0
                    avg_y = (pdata(i,j,k,2) + pdata(i+1,j,k,2) + pdata(i,j+1,k,2) + pdata(i+1,j+1,k,2))/4.0
                    avg_z = (pdata(i,j,k,3) + pdata(i+1,j,k,3) + pdata(i,j+1,k,3) + pdata(i+1,j+1,k,3))/4.0
                    norm = sqrt( avg_x**2 + avg_y**2 + avg_z**2 )
                    write(fd,'(e14.6)') norm
                end do
            end do
        else if(npts(1) > 1 .and. npts(2) == 1 .and. npts(3) > 1) then
            j=1
            do k=1,npts(3)-1
                do i=1,npts(1)-1
                    avg_x = (pdata(i,j,k,1) + pdata(i+1,j,k,1) + pdata(i,j,k+1,1) + pdata(i+1,j,k+1,1))/4.0
                    avg_y = (pdata(i,j,k,2) + pdata(i+1,j,k,2) + pdata(i,j,k+1,2) + pdata(i+1,j,k+1,2))/4.0
                    avg_z = (pdata(i,j,k,3) + pdata(i+1,j,k,3) + pdata(i,j,k+1,3) + pdata(i+1,j,k+1,3))/4.0
                    norm = sqrt( avg_x**2 + avg_y**2 + avg_z**2 )
                    write(fd,'(e14.6)') norm
                end do
            end do
        else if(npts(1) == 1 .and. npts(2) > 1 .and. npts(3) > 1) then
            i=1
            do k=1,npts(3)-1
                do j=1,npts(2)-1
                    avg_x = (pdata(i,j,k,1) + pdata(i,j+1,k,1) + pdata(i,j,k+1,1) + pdata(i,j+1,k+1,1))/4.0
                    avg_y = (pdata(i,j,k,2) + pdata(i,j+1,k,2) + pdata(i,j,k+1,2) + pdata(i,j+1,k+1,2))/4.0
                    avg_z = (pdata(i,j,k,3) + pdata(i,j+1,k,3) + pdata(i,j,k+1,3) + pdata(i,j+1,k+1,3))/4.0
                    norm = sqrt( avg_x**2 + avg_y**2 + avg_z**2 )
                    write(fd,'(e14.6)') norm
                end do
            end do
        end if

        write(fd, *) '   </CellData>'

        write(fd, *) '   </Piece>'
        write(fd, *) '   </ImageData>'
        write(fd, *) '</VTKFile>'

        call closefd(fd)
    end subroutine


    ! fname: output filename
    ! grid: grid point coordinates
    ! pdata: vectors at grid points
    ! cell definitions
    subroutine write_vtk_vector_unstructuredgrid(fname, grid, pdata, cells)
    ! for writing vtu files:
    ! https://lorensen.github.io/VTKExamples/site/VTKFileFormats/#unstructuredgrid
        character(*), intent(in)   :: fname
        real(real64), intent(in)   :: grid(:, :)  ! which is (3, npoints)
        real(real64), intent(in)   :: pdata(:, :) ! which is (npoints, 3)
        integer(int32), intent(in) :: cells(:, :) ! which is (4, ncells)

        integer(int32) :: fd
        integer(int32) :: p, c
        integer(int32) :: npoints, ncells

        npoints = size(grid, 2)
        ncells = size(cells, 2)
        write(*,*) npoints, ncells
      
        call getfd(fd)
        open(fd, file=trim(fname), form='formatted', status='unknown')
        write(*,*) fname, " opened"

        write(fd, '(a)') '<?xml version="1.0"?>'
        write(fd, '(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
        write(fd, '(a)') '  <UnstructuredGrid>'
        write(fd, '(a, i5, a, i5, a)') '    <Piece NumberOfPoints="', npoints, '" NumberOfCells="', ncells, '">'
        write(fd, '(a)') '      <Points>'
        write(fd, '(a)') '        <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
        do p=1, npoints ! loop over points, write coords
            write(fd, '(a, 3e20.10)') '        ', grid(1,p), grid(2,p), grid(3,p)
        end do
        write(fd, '(a)') '        </DataArray>'
        write(fd, '(a)') '      </Points>'
        write(fd, '(a)') '      <PointData Scalars="scalars">'
        write(fd, '(a)') '        <DataArray Name="vectors" type="Float64" NumberOfComponents="3" Format="ascii">'
        do p=1, npoints ! loop over points, write vectors
            write(fd, '(a, 3e20.10)') '        ', pdata(p,1), pdata(p,2), pdata(p,3)
        end do
        write(fd, '(a)') '        </DataArray>'
        write(fd, '(a)') '      </PointData>'
        write(fd, '(a)') '      <Cells>'
        write(fd, '(a)') '        <DataArray type="Int32" Name="connectivity" Format="ascii">'
        do c=1, ncells ! loop over cells, print node indices (base 0)
            write(fd, '(a, 4i10)') '        ', cells(1,c)-1, cells(2,c)-1, cells(3,c)-1, cells(4,c)-1
        end do
        write(fd, '(a)') '        </DataArray>'
        write(fd, '(a)') '        <DataArray type="Int32" Name="offsets" Format="ascii">'
        do c=1, ncells ! loop over cells, print offsets, which are effectively the index of the last element in each tetrahedron (base 1)
            write(fd, '(a,i4)') '        ', 4*c
        end do
        write(fd, '(a)') '        </DataArray>'
        write(fd, '(a)') '        <DataArray type="Int32" Name="types" Format="ascii">'
        do c=1, ncells ! loop over cells, print cell type.  tetrahedra have type '10'
            write(fd, '(a,i3)') '        ', 10
        end do
        write(fd, '(a)') '        </DataArray>'
        write(fd, '(a)') '      </Cells>'
        write(fd, '(a)') '      <CellData Scalars="foo">'
        do c=1, ncells ! loop over cells, print a value ... what about 0.0?
            write(fd, '(a,f4.1)') '        ', 0.0
        end do
        write(fd, '(a)') '      </CellData>'
        write(fd, '(a)') '    </Piece>'
        write(fd, '(a)') '  </UnstructuredGrid>'
        write(fd, '(a)') '</VTKFile>'

        call closefd(fd)
    end subroutine

end module

