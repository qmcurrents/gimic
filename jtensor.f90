!
! This is the actual work horse, calculates the current tensor
! for a particular spin case. The cost of calling jtensor twice is 
! very small, since no exponentials have to be calculated.
! 
module jtensor_class
	use globals_m
	use basis_class
	use bfeval_class
	use dens_class
	use dfdb_class
	use dfdr_class
	use d2fdrdb_class
	use dbop_class
	use grid_class
	use teletype_m
	implicit none
	
	type jtensor_t
		type(molecule_t), pointer :: mol
		type(dens_t), pointer :: xdens
		type(dbop_t) :: dop
		type(dfdb_t) :: dbt
		type(d2fdrdb_t) :: 	d2f
		type(dfdr_t) :: dfr
		type(bfeval_t) :: bfv
		real(DP), dimension(:), pointer :: pdbf, denbf, dendb
	end type

	public init_jtensor, del_jtensor, jtensor, jtensor2, jvector
	public ctensor, ctensor2
	public qtensor, qtensor2
	public jtensor_t, jdebug, eta
	
	private
	
	real(DP), dimension(3) :: rho

	real(DP), dimension(:), pointer :: bfvec
	real(DP), dimension(:,:), pointer :: dbvec, drvec, d2fvec
	real(DP), dimension(:,:), pointer :: aodens, pdens
	logical :: diamag_p, paramag_p, spin_density
	integer(I4), parameter :: NOTIFICATION=1000

contains
	! set up memory (once) for the different components
	subroutine init_jtensor(jt, mol, xdens)
		type(jtensor_t) :: jt
		type(molecule_t), target :: mol
		type(dens_t), target :: xdens
		integer(I4) ::  ncgto
		
		ncgto=get_ncgto(mol)
		
!        call init_dens(ddens)
        jt%mol=>mol
		jt%xdens=>xdens
		call init_bfeval(jt%bfv, jt%mol)
		call init_dbop(jt%dop, jt%mol)
		call init_dfdb(jt%dbt, jt%mol, jt%dop, jt%bfv)
		call init_dfdr(jt%dfr, jt%mol)
		call init_d2fdrdb(jt%d2f, jt%mol, jt%dop, jt%dfr, jt%bfv)

		diamag_p=.true.
		paramag_p=.true.
		spin_density=.false.
		call getkw(input, 'cdens.diamag', diamag_p)
		call getkw(input, 'cdens.paramag', paramag_p)
		if (uhf_p) call getkw(input, 'cdens.spin_density', spin_density)

		if (diamag_p == 0) then
			call msg_info( 'Diamagnetic contributions not calculated!')
			call nl
		end if
		if (paramag_p == 0) then
			call msg_info( 'Paramagnetic contributions not calculated!')
			call nl
		end if
		if (diamag_p == 0 .and. paramag_p == 0) then
			call msg_out( '    ...this does not make sense...')
			call nl
			call msg_critical( '    PLEASE SEEK PROFESSIONAL HELP, ASAP!  ')
			call nl
			stop
		end if
		! intermediates
		allocate(jt%denbf(ncgto))
		allocate(jt%dendb(ncgto))
		allocate(jt%pdbf(ncgto))
	end subroutine

	subroutine qtensor(jt, r)
		type(jtensor_t) :: jt
		real(DP), dimension(3), intent(in) :: r

		type(tensor_t) :: jt1, jt2

		if (uhf_p) then
			call jtensor(jt, r, jt1, spin_a)
			call jtensor(jt, r, jt2, spin_b)
			write(81, *) jt1%t
			write(82, *) jt2%t
			write(83, *) jt1%t+jt2%t
			write(84, *) jt1%t-jt2%t
		else
			call jtensor(jt, r, jt1, spin_a)
			write(85, *) jt1%t
		end if 
	end subroutine

	subroutine qtensor2(jt, r)
		type(jtensor_t) :: jt
		real(DP), dimension(3), intent(in) :: r

		type(tensor_t) :: pjt1, pjt2
		type(tensor_t) :: djt1, djt2

		if (uhf_p) then
			call jtensor2(jt, r, pjt1, djt1, spin_a)
			call jtensor2(jt, r, pjt2, djt2, spin_b)
			write(90, *) pjt1%t
			write(91, *) pjt2%t
			write(92, *) djt1%t
			write(93, *) djt2%t
			write(94, *) pjt1%t+pjt2%t
			write(95, *) djt1%t+djt2%t
			write(96, *) pjt1%t-djt1%t
			write(97, *) pjt2%t-djt2%t
		else
			call jtensor2(jt, r, pjt1, djt1, spin_a)
			write(98, *) pjt1%t
			write(99, *) djt1%t
		end if 
	end subroutine

	subroutine ctensor(jt, r, j, op)
		type(jtensor_t) :: jt
		real(DP), dimension(3), intent(in) :: r
		type(tensor_t), intent(inout) :: j
		character(*) :: op

		type(tensor_t) :: jt1, jt2

		select case (op)
			case ('alpha')
				call jtensor(jt, r, j, spin_a)
			case ('beta')
				if (uhf_p) then
					call jtensor(jt, r, j, spin_b)
				else
					call msg_error('ctensor(): &
					&beta current requested, but not open-shell system!')
					stop
				end if
			case ('total')
				if (uhf_p) then
					call jtensor(jt, r, jt1, spin_a)
					call jtensor(jt, r, jt2, spin_b)
					j%t=jt1%t+jt2%t
				else
					call jtensor(jt, r, j, spin_a)
				end if
			case ('spindens')
				if (.not.uhf_p) then
					call msg_error('ctensor(): &
					&spindens requested, but not open-shell system!')
					stop
				end if
				call jtensor(jt, r, jt1, spin_a)
				call jtensor(jt, r, jt2, spin_b)
				j%t=jt1%t-jt2%t
		end select
	end subroutine

	subroutine ctensor2(jt, r, pj, dj, op)
		type(jtensor_t) :: jt
		real(DP), dimension(3), intent(in) :: r
		type(tensor_t), intent(inout) :: pj, dj
		character(*) :: op

		type(tensor_t) :: pj1, pj2, dj1, dj2

		select case (op)
			case ('alpha')
				call jtensor2(jt, r, pj, dj, spin_a)
			case ('beta')
				if (uhf_p) then
					call jtensor2(jt, r, pj, dj, spin_b)
				else
					call msg_error('ctensor(): &
					&beta current requested, but not open-shell system!')
					stop
				end if
			case ('total')
				if (uhf_p) then
					call jtensor2(jt, r, pj1, dj1, spin_a)
					call jtensor2(jt, r, pj2, dj2, spin_b)
					pj%t=pj1%t+pj2%t
					dj%t=dj1%t+dj2%t
				else
					call jtensor2(jt, r, pj, dj, spin_a)
				end if
			case ('spindens')
				if (.not.uhf_p) then
					call msg_error('ctensor(): &
					&spindens requested, but not open-shell system!')
					stop
				end if
				call jtensor2(jt, r, pj1, dj1, spin_a)
				call jtensor2(jt, r, pj2, dj2, spin_b)
				pj%t=pj1%t-dj1%t
				dj%t=pj2%t-dj2%t
		end select
	end subroutine

	subroutine jtensor(jt, r, j, spin)
		type(jtensor_t) :: jt
		real(DP), dimension(3), intent(in) :: r
		type(tensor_t), intent(inout) :: j
		integer(I4) :: spin

		integer(I4) :: i, b
		integer(I4), save :: notify=1

		rho=DP50*r ! needed for diamag. contr.

		call bfeval(jt%bfv, r, bfvec)
		call dfdb(jt%dbt, r, dbvec)
		call dfdr(jt%dfr, r, drvec)
		call d2fdrdb(jt%d2f, r, d2fvec)

		call contract(jt, j%t, spin)

!        write(75,*) jt1(:,:,1)
!        write(76,*) jt1(:,:,2)
!        write(77,*) j%t
!        write(78,*) j%t/2.d0


!		if (mod(notify,NOTIFICATION) == 0) then
!			print '(a, i6)', '* points done:', notify
!			call flush(6)
!		end i
!		notify=notify+1
		
	end subroutine

	subroutine jtensor2(jt, r, pj, dj, spin)
		type(jtensor_t) :: jt
		real(DP), dimension(3), intent(in) :: r
		type(tensor_t), intent(inout) :: pj, dj
		integer(I4) :: spin

		real(DP), dimension(:,:), pointer :: dbop
		integer(I4) :: i, b
		integer(I4), save :: notify=1

		rho=DP50*r ! needed for diamag. contr.

		call bfeval(jt%bfv, r, bfvec)
		call dfdb(jt%dbt, r, dbvec)
		call dfdr(jt%dfr, r, drvec)
		call d2fdrdb(jt%d2f, r, d2fvec)

		call contract2(jt, pj%t, dj%t, spin)

	end subroutine

	subroutine jvector(pj, dj, bb, jv)
		real(DP), dimension(:,:), intent(in) :: pj, dj
		real(DP), dimension(:), intent(in) :: bb
		real(DP), dimension(:), intent(out) :: jv
		
		jv=matmul(pj+dj, bb)
	end subroutine

	subroutine contract2(jt, ctp, ctd, spin)
		type(jtensor_t) :: jt
		real(DP), dimension(3,3), intent(out) :: ctp, ctd
		integer(I4), intent(in) :: spin

		integer(I4) :: i, j, k, ii,jj
		real(DP) :: prsp1, prsp2       ! paramagnetic wavefunction response
		real(DP) :: ppd                ! paramagnetic probability density
		real(DP), dimension(3) :: dpd  ! diamagnetic probability density
		real(DP) :: diapam

		call get_dens(jt%xdens, aodens, spin)  
		jt%denbf=matmul(bfvec, aodens)

		k=1
		diapam=dot_product(jt%denbf, bfvec)
		do i=1,3! dB <x,y,z>
			! get perturbed densities: x,y,z
			call get_pdens(jt%xdens, i, pdens,spin) 
			jt%pdbf=matmul(bfvec, pdens)
			jt%dendb=matmul(dbvec(:,i), aodens)
			dpd(i)=diapam*rho(i) ! diamag. contr. to J
			do j=1,3 !dm <x,y,z>
				prsp1=-dot_product(jt%dendb, drvec(:,j))    ! (-i)**2=-1 
				prsp2=dot_product(jt%denbf, d2fvec(:,k))
				ppd=dot_product(jt%pdbf, drvec(:,j))
				ctp(j,i)=ZETA*(ppd+prsp1+prsp2)
!                ctp(j,i)=ZETA*ppd
				k=k+1
			end do
		end do

! Calculate the total J tensor. The diamagnetic pobability density only
! contributes to the diagonal.

		ctd(1,1)=D0; ctd(2,2)=D0; ctd(3,3)=D0;
		ctd(1,2)=dpd(3)
		ctd(1,3)=-1.d0*dpd(2)
		ctd(2,1)=-1.d0*dpd(3)
		ctd(2,3)=dpd(1)
		ctd(3,1)=dpd(2)
		ctd(3,2)=-1.d0*dpd(1)

	end subroutine

!
! Contract all contributions to J. This is where all the actual work is done.
!
	subroutine contract(jt, ct, spin)
		type(jtensor_t) :: jt
		real(DP), dimension(3,3), intent(out) :: ct
		integer(I4), intent(in) :: spin

		integer(I4) :: b, m, k, ii,jj
		real(DP) :: prsp1, prsp2       ! paramagnetic wavefunction response
		real(DP) :: ppd                ! paramagnetic probability density
		real(DP), dimension(3) :: dpd  ! diamagnetic probability density
		real(DP) :: diapam

		call get_dens(jt%xdens, aodens, spin)  
		
		jt%denbf=matmul(bfvec, aodens)

		k=1
		diapam=dot_product(jt%denbf, bfvec)
		do b=1,3! dB <x,y,z>
			! get perturbed densities: x,y,z
			call get_pdens(jt%xdens, b, pdens, spin) 
			jt%pdbf=matmul(bfvec, pdens)
			jt%dendb=matmul(dbvec(:,b), aodens)
			dpd(b)=diapam*rho(b) ! diamag. contr. to J
			do m=1,3 !dm <x,y,z>
				prsp1=-dot_product(jt%dendb, drvec(:,m))    ! (-i)**2=-1 
				prsp2=dot_product(jt%denbf, d2fvec(:,k))
				ppd=dot_product(jt%pdbf, drvec(:,m))
!                ct(m,b)=ZETA*ppd
				ct(m,b)=ZETA*(ppd+prsp1+prsp2)
				k=k+1
			end do
		end do

! Calculate the total J tensor. The diamagnetic pobability density only
! contributes to the diagonal.

		! annihilate paramagnetic contribution
		if (.not.paramag_p)	then
			ct=D0
			bert_is_evil=1
		end if
		! annihilate diamagnetic  contribution
		if (.not.diamag_p) then
			dpd=D0    
			bert_is_evil=1
		end if

		ct(1,2)=ct(1,2)+dpd(3)
		ct(1,3)=ct(1,3)-dpd(2)
		ct(2,1)=ct(2,1)-dpd(3)
		ct(2,3)=ct(2,3)+dpd(1)
		ct(3,1)=ct(3,1)+dpd(2)
		ct(3,2)=ct(3,2)-dpd(1)

	end subroutine

	subroutine del_jtensor(jt)
		type(jtensor_t) :: jt

		call del_dbop(jt%dop)
		call del_dfdb(jt%dbt)
		call del_bfeval(jt%bfv)
		call del_dfdr(jt%dfr)
		call del_d2fdrdb(jt%d2f)
		deallocate(jt%denbf)
		deallocate(jt%pdbf)
		deallocate(jt%dendb)
	end subroutine

	subroutine eta(jj, grid, fac)
		type(jtensor_t) :: jj
		type(grid_t) :: grid
		real(DP), intent(in), optional :: fac

		external dtime
		real(4) :: dtime
		integer(I4) :: i, j, p1, p2, p3
		type(tensor_t) :: foo
		real(4) :: delta_t
		real(4), dimension(2) :: times
		real(DP), dimension(3) :: bar=(/D1,D1,D1/)
		real(DP), dimension(3) :: foobar
		real(DP), parameter :: SC=1.d-2
		
!        return
		call get_grid_size(grid, p1, p2, p3)
		
		delta_t=dtime(times)
		do i=1,100
			call jtensor(jj, (/i*SC, i*SC, i*SC/), foo, spin_a)
			foobar=matmul(bar,foo%t)
		end do
		delta_t=dtime(times)
		if ( present(fac) ) times=times*fac
		write(str_g, '(a,f9.2,a)') '[Under]Estimated wall time for &
			&calculation: ', times(1)*real(p1*p2*p3)/100.d0, ' sec'
		call msg_info(str_g)
		call nl
	end subroutine	

! Debug function...
!    function contract_killer(nuc, mag, k) result(jj)
!        integer(I4), intent(in) :: nuc, mag, k
!
!        real(DP) :: jj
!        integer(I4) :: mu, nu, n
!        real(DP) :: q
!
!        n=get_ncgto()
!        q=0.0
!        do mu=1,n
!            do nu=1,n 
!                q=q+aodens(nu,mu)*( dbvec(mu,mag)*drvec(nu,nuc)- &
!                  bfvec(mu)*(d2fvec(nu,k)))- &
!                  pdens(nu,mu)*bfvec(mu)*drvec(nu,nuc)
!            end do
!        end do
!        jj=q
!    end function
	
	subroutine jdebug(jt)
		type(jtensor_t) :: jt

		integer(I4) :: i, b
		integer(I4), save :: notify=1
		real(DP), dimension(3) :: r

		call getkw(input, 'grid.origin', r)

		call bfeval(jt%bfv, r, bfvec)
		call dfdb(jt%dbt, r, dbvec)
		call dfdr(jt%dfr, r, drvec)
		call d2fdrdb(jt%d2f, r, d2fvec)
		
		print *, 'bfvec'
		print *, repeat('-', 70)
		print *, bfvec
		print *

		print *, 'drvec'
		print *, repeat('-', 70)
		print 45
		print 41, drvec(1,:)
		print 42, drvec(2,:)
		print 43, drvec(3,:)
		print 44, drvec(4,:)
		print *

		print *, 'dbvec'
		print *, repeat('-', 70)
		print 45
		print 41, dbvec(1,:)
		print 42, dbvec(2,:)
		print 43, dbvec(3,:)
		print 44, dbvec(4,:)
		print *

		print *, 'd2fvec'
		print *, repeat('-', 70)
		print *, 'X'
		print 45
		print 41, d2fvec(1,1:3)
		print 42, d2fvec(2,1:3)
		print 43, d2fvec(3,1:3)
		print 44, d2fvec(4,1:3)
		print *
		print *, 'Y'
		print 45
		print 41, d2fvec(1,4:6)
		print 42, d2fvec(2,4:6)
		print 43, d2fvec(3,4:6)
		print 44, d2fvec(4,4:6)
		print *
		print *, 'Z'
		print 45
		print 41, d2fvec(1,7:9)
		print 42, d2fvec(2,7:9)
		print 43, d2fvec(3,7:9)
		print 44, d2fvec(4,7:9)
41 format('s ',3f15.10)
42 format('px',3f15.10)
43 format('py',3f15.10)
44 format('pz',3f15.10)
45 format('===        x             y               z')
	end subroutine
	
end module
