!
! Write VTK ImageData XML files
!
module vtkplot_module
    use globals_module
    use settings_module
    use grid_class

    ! Maria debugging
    use teletype_module
    ! end debugging

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
        integer(I4) :: npoints ! 
        real(DP), dimension(2) :: qrange
        real(DP), dimension(3) :: qmin, qmax, step
        real(DP) :: norm 

        ! Maria debugging
        ! character(120) :: debug_str
        ! end debugging

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

! Maria debugging
!        write(debug_str, '(3(a,i0,x))') "### Debugging: npts 1", npts(1), "npts 2", npts(2), "npts 3", npts(3);
!        call msg_out(debug_str)
!        write(debug_str, '(3(a,f8.4,x))') "### Debugging: qmin 1", qmin(1), "qmin 2", qmin(2), "qmin 3", qmin(3);
!        call msg_out(debug_str)
!        write(debug_str, '(3(a,f8.4,x))') "### Debugging: qmax 1", qmax(1), "qmax 2", qmax(2), "qmax 3", qmax(3);
!        call msg_out(debug_str)
!        write(debug_str, '(i0)') npoints;
!        call msg_out(debug_str)
!        write(debug_str, '(3(a,f8.4,x))') "### Debugging: step 1", step(1), "step 2", step(2), "step 3", step(3);
!        call msg_out(debug_str)
! end debugging

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
                    write(fd,'(3e14.6)') pdata(i,j,k,1), &
                    pdata(i,j,k,2), pdata(i,j,k,3)
                end do
            end do
        end do
        write(fd, *) 
        write(fd, *) '   </DataArray>'
        write(fd, *) '   </PointData>'

        write(fd, *) '   <CellData Scalars="foo">'

        do k=1,npts(3)
            do j=1,npts(2)
                do i=1,npts(1)
                    norm = sqrt(pdata(i,j,k,1)**2       & 
                    + pdata(i,j,k,2)**2                 &
                    + pdata(i,j,k,3)**2)
                    write(fd,'(e14.6)') norm 
                end do
            end do
        end do

        write(fd, *) '   </CellData>'

        write(fd, *) '   </Piece>'
        write(fd, *) '   </ImageData>'
        write(fd, *) '</VTKFile>'

        call closefd(fd)

    end subroutine

!    subroutine write_cml(this, fname, mol)
!        type(grid_t), intent(inout) :: this
!        character(*), intent(in) :: fname
!        type(molecule_t) :: mol 
!
!        integer(I4) :: natoms, i
!        integer(I4) :: p1, p2, p3
!        real(DP), dimension(3) :: r, coord
!        character(2) :: symbol
!        type(atom_t), pointer :: atom
!        logical :: show_axis = .true.
!
!        natoms=get_natoms(mol)
!
!        call get_grid_size(this,p1,p2,p3)
!        !  call getkw(input, 'show_axis', show_axis)
!        i=0
!        if (show_axis) i=1 
!
!
!        write(str_g, '(2a)') 'Grid plot in ', trim(fname)
!        call msg_note(str_g)
!        open(77,file=trim(fname))
!
!        if (p3 > 1) then
!            write(77,*) natoms+8+i
!            write(77,*)
!        else if (p3 == 1) then
!            write(77,*) natoms+4+i
!            write(77,*)
!        else
!            write(77,*) natoms
!            write(77,*)
!        end if
!
!        do i=1,natoms
!            call get_atom(mol, i, atom)
!            call get_symbol(atom, symbol)
!            call get_coord(atom, coord)
!            write(77,'(a, 3f16.10)') symbol, coord*au2a
!        end do
!
!        if (p3 > 1) then
!            r=gridpoint(this,1,1,1)
!            write(77,'(a,3f16.10)') 'X ', r*au2a
!            r=gridpoint(this,p1,1,1)
!            write(77,'(a,3f16.10)') 'X ', r*au2a
!            r=gridpoint(this,1,p2,1)
!            write(77,'(a,3f16.10)') 'X ', r*au2a
!            r=gridpoint(this,1,1,p3)
!            write(77,'(a,3f16.10)') 'X ', r*au2a
!            r=gridpoint(this,p1,p2,1)
!            write(77,'(a,3f16.10)') 'X ', r*au2a
!            r=gridpoint(this,p1,1,p3)
!            write(77,'(a,3f16.10)') 'X ', r*au2a
!            r=gridpoint(this,1,p2,p3)
!            write(77,'(a,3f16.10)') 'X ', r*au2a
!            r=gridpoint(this,p1,p2,p3)
!            write(77,'(a,3f16.10)') 'X ', r*au2a
!        else if (p3 == 1) then
!            r=gridpoint(this,1,1,1)
!            write(77,'(a,3f16.10)') 'X ', r*au2a
!            r=gridpoint(this,p1,1,1)
!            write(77,'(a,3f16.10)') 'X ', r*au2a
!            r=gridpoint(this,1,p2,1)
!            write(77,'(a,3f16.10)') 'X ', r*au2a
!            r=gridpoint(this,p1,p2,1)
!            write(77,'(a,3f16.10)') 'X ', r*au2a
!        end if
!
!        if (show_axis) then
!            select case (trim(this%mode))
!                case ('std','base')
!                    write(77,'(a,3f16.10)') 'Be ', &
!                        (this%origin+this%basv(:,3)*2.0)*au2a
!                case ('bond')
!                    write(77,'(a,3f16.10)') 'Be ', &
!                        (this%origin+this%ortho*2.0)*au2a
!            end select
!        end if
!
!        close(77)
!    end subroutine

end module

! vim:et:sw=4:ts=4
