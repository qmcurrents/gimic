! This module calulates the current field by looping over gridpoints and 
! calling the appropriate routines in jtensor_class
!
! Ok, this is really a mess...
! This modules needs a seriuos clean up. I have been adding, and adding, and
! adding stuff and it's not pretty.
!

module jfield_class
    use globals_m
    use jtensor_class
    use grid_class
    use basis_class
    use teletype_m
    use parallel_m
    use tensor_m
    implicit none

    type jfield_t
        type(jtensor_t), pointer :: jt
        real(DP), dimension(3) :: b
        type(tensor_t), dimension(:,:), pointer :: jj, jja, jjb, jsd
        type(vector_t), dimension(:,:), pointer :: vv, vva, vvb, vsd
        integer(I4) :: zpos_j, zpos_v
        type(grid_t), pointer :: grid
    end type

    public new_jfield, del_jfield, jfield, jvectors, jfield_direct
    public jfield_t, get_jvectors, get_jtensors, jvector_plot, jfield_eta
    public set_jtensors
    
    private

    real(DP) :: vec_scale=D1
    integer(I4) :: spin
    character(BUFLEN) :: jtensor_file, jvec_file
    character(BUFLEN) :: jmod_plt, jvec_plt, njvec_plt, jprj_plt
contains
    subroutine new_jfield(this, jt, g, magnet)
        type(jfield_t) :: this
        type(jtensor_t), target :: jt
        type(grid_t), target :: g
        real(DP), dimension(3), intent(in) :: magnet

        integer(I4) :: i, j, k

        this%b=magnet
        this%jt=>jt
        this%grid=>g

        if (mpi_rank == 0) then
            call setup_files(this)
        end if

        call get_grid_size(g,i,j,k)
        allocate(this%jj(i,j))
        allocate(this%vv(i,j))
        if (uhf_p) then
            allocate(this%jja(i,j))
            allocate(this%jjb(i,j))
            allocate(this%jsd(i,j))
            allocate(this%vva(i,j))
            allocate(this%vvb(i,j))
            allocate(this%vsd(i,j))
        end if

!        call getkw(input, 'scale_vectors', vec_scale)
!        if ( vec_scale /= D1 ) then
!            write(str_g, '(a,f7.3)') 'Vector scaling =', vec_scale
!            call msg_note(str_g)
!        end if
!        call nl

    end subroutine

    subroutine setup_files(this)
        type(jfield_t) :: this
    
        logical :: plot_p
        integer(I4) :: jtrecl, jvrecl

        character(*), parameter :: jtensor_file = 'jtensor'
        character(*), parameter :: jvector_file = 'jvector'

        jtrecl=9*DP
        jvrecl=3*DP
        
        open(JTFD, file=trim(jtensor_file) // '.bin', &
            access='direct', recl=jtrecl)
        open(JVECFD, file=trim(jvector_file) // '.bin', &
            access='direct', recl=jvrecl)
        if (uhf_p) then
            open(JTFD+1, file=trim(jtensor_file) // '_A.bin', &
            access='direct', recl=jtrecl)
            open(JTFD+2, file=trim(jtensor_file) // '_B.bin', &
            access='direct', recl=jtrecl)
            open(JTFD+3, file=trim(jtensor_file) // '_SD.bin', &
            access='direct', recl=jtrecl)
            open(JVECFD+1, file=trim(jvec_file) // '_A.bin', &
            access='direct', recl=jvrecl)
            open(JVECFD+2, file=trim(jvec_file) // '_B.bin', &
            access='direct', recl=jvrecl)
            open(JVECFD+3, file=trim(jvec_file) // '_SD.bin', &
            access='direct', recl=jvrecl)
        end if

    end subroutine

    subroutine del_jfield(this)
        type(jfield_t), intent(inout) :: this
        
        this%b=0.d0
        deallocate(this%jj)
        deallocate(this%vv)
        nullify(this%jj)
        nullify(this%vv)
        if (uhf_p) then
            deallocate(this%jja)
            deallocate(this%jjb)
            deallocate(this%jsd)
            deallocate(this%vva)
            deallocate(this%vvb)
            deallocate(this%vsd)
        end if
        
        if (mpi_rank == 0) then
            close(JTFD)
            close(JVECFD)
            if (uhf_p) then
            close(JTFD+1)
            close(JTFD+2)
            close(JTFD+3)
            close(JVECFD+1)
            close(JVECFD+2)
            close(JVECFD+3)
            end if
        end if
    end subroutine
    
    subroutine jvectors(this)
        type(jfield_t), intent(inout) :: this

        integer(I4) :: i, j, k, p1, p2, p3

        call get_grid_size(this%grid, p1, p2, p3)

        call msg_note( 'Contracting j-tensors with magnetic field.')
        call msg_note( 'Writing j-vectors in: ' // trim(jvec_file))
        call nl

        do k=1,p3
            call jtens_io(this, k, 'r')
            this%zpos_j=k
            do j=1,p2
                do i=1,p1
                    this%vv(i,j)%v=matmul(this%jj(i,j)%t, this%b)
                    if (uhf_p) then
                        this%vva(i,j)%v=matmul(this%jja(i,j)%t, this%b)
                        this%vvb(i,j)%v=matmul(this%jjb(i,j)%t, this%b)
                        this%vsd(i,j)%v=matmul(this%jsd(i,j)%t, this%b)
                    end if
                end do
            end do
            call jvec_io(this, k, 'w')
            this%zpos_v=k
        end do

    end subroutine

    subroutine jtens_io(this, p3, op)
        type(jfield_t), intent(inout) :: this
        integer(I4), intent(in) :: p3
        character(1), intent(in) :: op

        integer(I4) :: i, j, k, p1, p2 

        call get_grid_size(this%grid, p1, p2)

        k=p1*p2*(p3-1)+1
        do j=1, p2
            do i=1, p1
                select case(op)
                case('r')
                    read(JTFD, rec=k) this%jj(i,j)%t
                    if (uhf_p) then
                        read(JTFD+1, rec=k) this%jja(i,j)%t
                        read(JTFD+2, rec=k) this%jjb(i,j)%t
                        read(JTFD+3, rec=k) this%jsd(i,j)%t
                    end if
                case('w')
                    write(JTFD, rec=k) this%jj(i,j)%t
                    if (uhf_p) then
                        write(JTFD+1, rec=k) this%jja(i,j)%t
                        write(JTFD+2, rec=k) this%jjb(i,j)%t
                        write(JTFD+3, rec=k) this%jsd(i,j)%t
                    end if
                case default
                    call msg_error('jtens_io(): Invalid operation: ' // op)
                    call exit(1)
                end select
                k=k+1
            end do
        end do
    end subroutine

    subroutine jfield_eta(this, fac)
        type(jfield_t) :: this
        real(DP), intent(in), optional :: fac

        integer(I4) :: i, j, p1, p2, p3
        type(tensor_t) :: foo
        real(DP) :: delta_t
        real(4) :: tim1, tim2
        real(4), dimension(2) :: times
        real(DP), dimension(3) :: bar=(/D1,D1,D1/)
        real(DP), dimension(3) :: foobar
        real(DP), parameter :: SC=0.25d0
        
        call get_grid_size(this%grid, p1, p2, p3)
        
        call etime(times, tim1)
        tim1=times(1)
        do i=1,100
            call jtensor(this%jt, (/i*SC, i*SC, i*SC/), foo, spin_a)
            foobar=matmul(bar,foo%t)
        end do
        call etime(times, tim2)
        tim2=times(1)
        
        delta_t=tim2-tim1
        if ( present(fac) ) delta_t=delta_t*fac
        write(str_g, '(a,f9.2,a)') 'Estimated CPU time for &
            &calculation: ', delta_t*real(p1*p2*p3)/100.d0, ' sec'
        call msg_info(str_g)
        call nl
    end subroutine	

    subroutine jvec_io(this, p3, op)
        type(jfield_t), intent(inout) :: this
        integer(I4), intent(in) :: p3
        character(1), intent(in) :: op

        integer(I4) :: i, j, k, p1, p2 

        call get_grid_size(this%grid, p1, p2)

        k=p1*p2*(p3-1)+1
        do j=1,p2
            do i=1,p1
                select case(op)
                case('r')
                    read(JVECFD, rec=k) this%vv(i,j)%v
                    if (uhf_p) then
                        read(JVECFD+1, rec=k) this%vva(i,j)%v
                        read(JVECFD+2, rec=k) this%vvb(i,j)%v
                        read(JVECFD+3, rec=k) this%vsd(i,j)%v
                    end if
                case('w')
                    write(JVECFD, rec=k) this%vv(i,j)%v
                    if (uhf_p) then
                        write(JVECFD+1, rec=k) this%vva(i,j)%v
                        write(JVECFD+2, rec=k) this%vvb(i,j)%v
                        write(JVECFD+3, rec=k) this%vsd(i,j)%v
                    end if
                case default
                    call msg_error('jvec_io(): Invalid operation: ' // op)
                    call exit(1)
                end select
                k=k+1
            end do
        end do
    end subroutine

    subroutine read_jvecs(buf, spin)
        type(vector_t), dimension(:,:,:), intent(out) :: buf
        integer(I4), intent(in) :: spin

        integer(I4) :: i, j, k, l, p1, p2, p3 

        p1=size(buf(:,1,1))
        p2=size(buf(1,:,1))
        p3=size(buf(1,1,:))

        if (spin > 1 .and. .not.uhf_p) then
            call msg_error("read_jvecs(): requested spin component for  & 
            &closed-shell!")
            stop 1
        end if

        l=1
        do k=1,p3
            do j=1,p2
                do i=1,p1
                    read(JVECFD+spin-1, rec=l) buf(i,j,k)%v
                    l=l+1
                end do
            end do
        end do
    end subroutine

    subroutine wrt_jvec(rr,v,fd)
        real(DP), dimension(:), intent(in) :: rr, v
        integer(I4), intent(in) :: fd

        if (fd == 0) return

        write(fd, '(6f11.7)')  rr, v
    end subroutine	

    subroutine wrt_njvec(rr,v,fd)
        real(DP), dimension(:), intent(in) :: rr, v
        integer(I4), intent(in) :: fd

        real(DP) :: nfac

        if (fd == 0) return

        nfac=sqrt(sum(v(:)**2))
        if (nfac < 1.d-15) then
            write(fd, '(6f11.7)')  rr, 0.d0, 0.d0, 0.d0
        else 
            write(fd, '(6f11.7)')  rr, v/nfac
        end if
    end subroutine	

    subroutine wrt_jproj(rr,v,grid,fd)
        real(DP), dimension(:), intent(in) :: rr, v
        type(grid_t) :: grid
        integer(I4), intent(in) :: fd

        real(DP) :: jprj
        real(DP), dimension(3) :: norm

        if (fd == 0) return
        norm=get_grid_normal(grid)

        jprj=dot_product(norm,v)
        write(fd, '(3e19.12)') rr, jprj
    end subroutine	

    subroutine wrt_jmod(rr,v,fd)
        real(DP), dimension(:), intent(in) :: rr, v
        integer(I4), intent(in) :: fd

        real(DP) :: jmod

        if (fd == 0) return

        jmod=sqrt(sum(v**2))
        write(fd, '(6f11.7)')  rr, jmod
    end subroutine	

    function getvecs(this, spin) result(jv)
        type(jfield_t), intent(inout) :: this
        integer(I4) :: spin

        type(vector_t), dimension(:,:), pointer :: jv
        select case(spin)
            case(1)
                jv=>this%vv
            case(2)
                jv=>this%vva
            case(3)
                jv=>this%vvb
            case(4)
                jv=>this%vsd
        end select
    end function

    function open_plot(basename, spin) result(fd)
        character(*), intent(in) :: basename
        integer(I4), intent(in) :: spin

        integer(I4) :: fd
        character(132) :: fname

        fd=0
        if (mpi_rank > 0) return
        if (trim(basename) == '') return

        call getfd(fd)
        if (fd == 0) then
            stop 1
        end if

        select case(spin)
            case(1)
                fname=trim(basename)//'.txt'
            case(2)
                fname=trim(basename)//'_a.txt'
            case(3)
                fname=trim(basename)//'_b.txt'
            case(4)
                fname=trim(basename)//'_sd.txt'
        end select

        write(str_g, *) 'Writing ', trim(basename), ' in ', trim(fname)
        call msg_note(str_g)
        open(fd,file=trim(fname),status='unknown')

        return 

    end function

    subroutine jvector_plot(this)
        type(jfield_t), intent(inout) :: this

        integer(I4) :: i, j, k, p1, p2, p3, spin, ispin
        integer(I4) :: fd1, fd2, fd3, fd4
        real(DP), dimension(3) :: v, rr
        type(vector_t), dimension(:,:), pointer :: jv

        spin=1
        if (uhf_p) then
            spin=4
        end if

        call get_grid_size(this%grid, p1, p2, p3)

        do ispin=1,spin
            fd1= open_plot('jvec',ispin)
            fd2= open_plot('jmod',ispin)
            !fd3= open_plot('jvec_projected',ispin)
            !fd4= open_plot('jvec_normalzied',ispin)

            do k=1,p3
                call jvec_io(this, k, 'r')
                jv=>getvecs(this,ispin)
                this%zpos_v=k
                do j=1,p2
                    do i=1,p1
                        rr=gridpoint(this%grid, i,j,k)*AU2A
                        v=jv(i,j)%v*AU2A
                        call wrt_jvec(rr,v,fd1)
                        call wrt_jmod(rr,v,fd2)
                        !call wrt_jproj(rr,v,this%grid,fd3)
                        !call wrt_njvec(rr,v,fd4)
                    end do
                    if (fd1 /= 0) write(fd1, *)
                    if (fd2 /= 0) write(fd2, *)
                    !if (fd3 /= 0) write(fd3, *)
                    !if (fd4 /= 0) write(fd4, *)
                end do
            end do

            call closefd(fd1)
            call closefd(fd2)
            !call closefd(fd3)
            !call closefd(fd4)
            if (p3 > 1) call jmod_cubeplot(this, 'jmod', ispin)
        end do
    end subroutine

    subroutine jfield(this)
        type(jfield_t) :: this

        integer(I4) :: i, j, k, p1, p2, p3
        real(DP), dimension(3) :: rr
        character(10) :: op
        integer(I4) :: lo, hi, npts

        call jfield_eta(this)
        call get_grid_size(this%grid, p1, p2, p3)

        call schedule(p2, lo, hi)
    
        do k=1,p3
            do j=lo,hi
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    call ctensor(this%jt, rr, this%jj(i,j), 'total')
                    if (uhf_p) then
                        call ctensor(this%jt, rr, this%jja(i,j), 'alpha')
                        call ctensor(this%jt, rr, this%jjb(i,j), 'beta')
                        call ctensor(this%jt, rr, this%jsd(i,j), 'spindens')
                    end if
                end do
            end do
            call gather_data(this%jj,this%jj(:,lo:hi))
            if (uhf_p) then
                call gather_data(this%jja,this%jja(:,lo:hi))
                call gather_data(this%jjb,this%jjb(:,lo:hi))
                call gather_data(this%jsd,this%jsd(:,lo:hi))
            end if
            if (mpi_rank == 0) call jtens_io(this, k, 'w')
            this%zpos_j=k
        end do
    end subroutine 

    subroutine jfield_direct(this, k)
        type(jfield_t) :: this
        integer(I4), intent(in) :: k

        integer(I4) :: i, j, p1, p2
        integer(I4) :: lo, hi, npts
        real(DP), dimension(3) :: rr
        type(tensor_t), dimension(:,:), pointer :: jj
        
        call get_grid_size(this%grid, p1, p2)
        jj=>this%jj

        call schedule(p2, lo, hi)

        do j=lo,hi
            do i=1,p1
                rr=gridpoint(this%grid, i, j, k)
                call ctensor(this%jt, rr, this%jj(i,j), 'total')
                if (uhf_p) then
                    call ctensor(this%jt, rr, this%jja(i,j), 'alpha')
                    call ctensor(this%jt, rr, this%jjb(i,j), 'beta')
                    call ctensor(this%jt, rr, this%jsd(i,j), 'spindens')
                end if
            end do
        end do
        call gather_data(this%jj,this%jj(:,lo:hi))
        if (uhf_p) then
            call gather_data(this%jja,this%jja(:,lo:hi))
            call gather_data(this%jjb,this%jjb(:,lo:hi))
            call gather_data(this%jsd,this%jsd(:,lo:hi))
        end if
    end subroutine 

    subroutine get_jtensors(this,z,jt)
        type(jfield_t), intent(inout) :: this
        integer(I4), intent(in) :: z
        type(tensor_t), dimension(:,:), pointer :: jt

        if (this%zpos_j /= z) then
            call jtens_io(this, z, 'r')
        end if

        jt=>this%jj
    end subroutine 

    subroutine set_jtensors(this,z,jt)
        type(jfield_t), intent(inout) :: this
        integer(I4), intent(in) :: z
        type(tensor_t), dimension(:,:) :: jt

        this%jj=jt
        call jtens_io(this, z, 'w')
        this%zpos_j=z
    end subroutine 

    subroutine get_jvectors(this,z,jvec)
        type(jfield_t), intent(inout) :: this
        integer(I4), intent(in) :: z
        type(vector_t), dimension(:,:), pointer :: jvec

        if (this%zpos_v /= z) then
            call jvec_io(this, z, 'r')
        end if

        jvec=>this%vv
    end subroutine 

    subroutine print_jt(rr, jt)
        real(DP), dimension(3), intent(in) :: rr
        real(DP), dimension(3,3), intent(in) :: jt

        integer(I4) :: l
        print *
        print '(a,3f12.8)', 'Current tensor at ', rr
        print *, '======================================='
        print '(3e19.12)', (jt(l,:), l=1,3)
        print *
        print '(a,e19.12)', ' Trace:', jt(1,1)+jt(2,2)+jt(3,3)
    end subroutine

    subroutine jmod_cubeplot(this, fname, spin)
        type(jfield_t) :: this
        character(*), intent(in) :: fname
        integer(I4), intent(in) :: spin

        integer(I4) :: p1, p2, p3, fd1, fd2
        integer(I4) :: i, j, k, l
        real(DP), dimension(3) :: qmin, qmax
        real(DP), dimension(3) :: norm, step, mag, v, rr
        real(DP) :: maxi, mini, val, sgn 
        integer(I4), dimension(3) :: npts
        type(vector_t), dimension(:,:,:), allocatable :: buf

        call get_grid_size(this%grid, p1, p2, p3)
        npts=(/p1,p2,p3/)
        norm=get_grid_normal(this%grid)
        qmin=gridpoint(this%grid,1,1,1)
        qmax=gridpoint(this%grid,p1,p2,p3)

        step=(qmax-qmin)/(npts-1)

        fd1=opencube(trim(fname), spin, qmin, step, npts)
        fd2=opencube(trim(fname)//'_quasi', spin, qmin, step, npts)

        mag = this%b

        allocate(buf(p1,p2,p3))
        call read_jvecs(buf,spin)
        maxi=0.d0
        mini=0.d0
        l=0
        do i=1,p1
            do j=1,p2
                do k=1,p3
                    v=buf(i,j,k)%v
                    rr=gridpoint(this%grid,i,j,k)
                    val=(sqrt(sum(v**2)))
                    rr=rr-dot_product(mag,rr)*mag
                    norm=cross_product(mag,rr)
                    sgn=dot_product(norm,v)
                    if (val > maxi) maxi=val
                    if (val < mini) mini=val
                    if (fd1 /= 0) then 
                        write(fd1,'(f12.6)',advance='no') val
                        if (mod(l,6) == 5) write(fd1,*)
                    end if
                    if (fd2 /= 0) then 
                        if (sgn >= 0.d0 ) then
                            write(fd2,'(f12.6)',advance='no') val
                        else
                            write(fd2,'(f12.6)',advance='no') -1.d0*val
                        end if
                        if (mod(l,6) == 5) write(fd2,*)
                    end if
                    l=l+1
                end do
            end do
        end do
        deallocate(buf)
        print *, 'maximini:', maxi, mini
    end subroutine

    function opencube(kname, spin, origin, step, npts) result(fd)
        character(*), intent(in) :: kname
        integer(I4), intent(in) :: spin
        real(DP), dimension(3), intent(in) :: origin, step
        integer(I4), dimension(3), intent(in) :: npts

        integer(I4) :: fd
        character(132) :: fname

        fd=0
        fname=kname
        if (mpi_rank > 0) return
        if (trim(fname) == '') return

        call getfd(fd)
        if (fd == 0) then
            stop 1
        end if

        select case(spin)
            case(1)
                fname=trim(fname)//'.cube'
            case(2)
                fname=trim(fname)//'a'//'.cube'
            case(3)
                fname=trim(fname)//'b'//'.cube'
            case(4)
                fname=trim(fname)//'sd'//'.cube'
        end select

        write(str_g, *) 'Writing ', trim(kname), ' in ', trim(fname)
        call msg_note(str_g)

        open(fd,file=trim(fname),form='formatted',status='unknown')
        write(fd,*) 'Gaussian cube data, generated by genpot'
        write(fd,*) 
        write(fd, '(i5,3f12.6)') 0, origin
        write(fd, '(i5,3f12.6)') npts(1), step(1), 0.d0, 0.d0
        write(fd, '(i5,3f12.6)') npts(2), 0.d0, step(2), 0.d0
        write(fd, '(i5,3f12.6)') npts(3), 0.d0, 0.d0, step(3)

        return 

    end function

!    subroutine read_jdata(f, g) 
!        type(jfield_t), intent(inout) :: f
!        type(grid_t), intent(inout) :: g
!
!        integer(I4) :: p1, p2, p3, i, j, k, foo_p 
!        
!        open(JTFD, file=jtensor_file, status='old')
!        print *, '*** Reading grid data from: ', xtrim(jtensor_file)
!        print *
!        call read_grid(g,JTFD)
!        call get_grid_size(g,p1,p2,p3)
!
!        allocate(f%j3(p1,p2,p3))
!        allocate(f%v3(p1,p2,p3))
!        
!        do k=1,p3
!            do j=1,p2
!                do i=1,p1
!                    read(JTFD, *) f%j3(i,j,k)%r
!                    read(JTFD, *) f%j3(i,j,k)%j
!                end do
!            end do
!        end do
!        close(JTFD)
!    end subroutine
end module

! vim:et:sw=4:ts=4
