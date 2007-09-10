!********************************************************************** 
!
! Common blocks for interfacing with ACES2
! 
!********************************************************************** 
module cpdcom
    use precision_m
	implicit none
	
	integer(I4), dimension(8,8) :: dirprd
	integer(I4), dimension(8,22) :: irpdpd
	integer(I4), dimension(100) :: iflags

	integer(I4) :: iintfp,nirrep,ntpert
	integer(I4) :: noca,nocb,nvrta,nvrtb
	integer(I4) :: iyzpert,ixzpert,ixypert

	! some global vars...
	integer(I4):: nbas, nbas2
	logical :: correlated, canonical

	! Unused...
	integer(I4), dimension(255,2) :: irreps
	integer(I4), dimension(2,500) :: isytyp
	integer(I4), dimension(18) :: id
	integer(I4), dimension(8) :: npert, ipert

	integer(I4) :: iintln,ifltln,ialone,ibitwd
	integer(I4) :: nstart
	integer(I4) :: ixpert,iypert,izpert
	integer(I4) :: itransx,itransy,itransz
	integer(I4) :: nucind									  

	common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
	common/syminf/nstart,nirrep,irreps,dirprd
	common/sympop/irpdpd,isytyp,id 
	common/flags/iflags
	common/info/noca,nocb,nvrta,nvrtb
	common/pert/ntpert,npert,ipert,ixpert,iypert,izpert, &
		     iyzpert,ixzpert,ixypert,itransx,itransy,itransz, &
		     nucind 
contains
	subroutine setup_commons
		if (iflags(2) == 0) then
			correlated=.false.
		else
			correlated=.true.
		end if
		
		if (iflags(64) == 0) then
			canonical=.false.
		else
			canonical=.true.
		end if
	end subroutine
end module


!********************************************************************** 
!
! This module handles unpacking of (symmetry) packed matrices
!
!********************************************************************** 
module sympack
    use precision_m
	use cpdcom
	implicit none

	public symc, matexp
	
	type symmat
		integer(I4) :: nbas
		real(DP), dimension(:,:), pointer :: mat
	end type
	private

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! This routine expands a compressed matrix A(p,q)
! p >= q to an array A(pq) with p,q. Note this routine 
! expects that the array A is symmetry packed.
!
!  INPUT : irrep  ...  the irrep of the corresponding part of a
!	       num ......  population vector for i and j
!	       a     ....  the matrix a
!	       ianti ..... 0 for symmetric and 1 for antisymmetric
!			    matrices
!
!          output : b .......  the expanded matrix a
!
! coded jg jan/91
! F90 verision jj/02
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine matexp(irrep,num,a,b,ianti)
		integer(I4), intent(in) :: irrep
		integer(I4), dimension(:), intent(in) :: num
		real(DP), dimension(:), intent(in) :: a
		real(DP), dimension(:,:), intent(inout) :: b
		integer(I4), intent(in) :: ianti
		
		integer(I4), dimension(8) :: ipold, ipnew
		
		integer(I4) :: i, j, irrepj, irrepi, numi, numj, ipo, ipn, ipnn
		integer(I4) :: ipoo, ij, ipi, ipj, idx
		integer(I4) :: ind1, ind2, ind3
		real(DP) :: anti

		if (ianti == 0) then
			anti=1.d0
		else
			anti=-1.d0
		end if

		if ( nirrep > 1) then
			print *, 'WARNING!!! matexp not tested properly with nirrep > 1!'
		end if
		!
		! take here care, if we are handling irrep=1 (total symmetric)
		! or irrep=1 (otherwise)
		!
		if (irrep == 1) then
			ipn=0
			ipo=0
			do irrepj=1,nirrep
				numj=num(irrepj)
		!
		! loop over orbitals, but also in backward order
		!
				do j=1,numj
					do i=1,j
						idx=(j*(j-1))/2+i+ipo
						b(ipn+i,ipn+j)=a(idx)
						b(ipn+j,ipn+i)=anti*a(idx)
					end do
				end do
				ipn=ipn+numj
				ipo=ipo+(numj**2+1)/2
			end do
		else
			ipi=0
			ipj=0
			ipo=0
			do irrepj=1,nirrep
				irrepi=dirprd(irrep,irrepj)
				numj=num(irrepj)
				numi=num(irrepi)
				if (irrepi < irrepj) then
					do j=1,numj
						do i=1,numi
							ipoo=(j-1)*numi+i+ipo
							b(ipi+i,ipj+j)=a(ipoo) !???
						end do
					end do
					ipo=ipo+numj*numi
				else
					do j=1,numj
						do i=1,numi
							ind1=(i-1)*numj+j+ipo
							b(ipi+i,ipj+j)=anti*a(ind1) !???
						end do
					end do
					ipo=ipo ! heh ;)
				end if
				ipi=ipi+numi
				ipj=ipj+numj
			end do
		end if
	end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  This subroutine returns the symmetry packed list of the input
!  matrix cfull. It is used in the integral derivative code to
!  symmetyry pack the eigen vector matrix.
!
!
!  cfull ..... input matrix ( size nbast*nbast)
!  csym ...... symmetry packed output matrix (square matrices within
!							   each irrep)
!  ibas ..... number of basis function per irrep( array of dim 8)
!
!   coded oct/90   jg
!   f90 version    jj (2002)
!
!   TODO: use types so that each irrep is a square matrix instead of
!         of an array.
!
!  type(irrmo)
!     integer(I4) :: nbas (basis funcs for this irrep)
!     real(DP), dimension(:,:), pointer :: imo
!  end type
!
!  now we can allocate an array of irrmo:s of size nirreps (actual). 
!  no more messing with hairy offsets :-)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine symc(cfull,csym,ibas)
		real(DP), dimension(:,:), intent(in) :: cfull
		real(DP), dimension(:), intent(out) :: csym
		integer(I4), dimension(8), intent(in) :: ibas

		integer(I4) ::  irrep, mo, l, ofs, off

		ofs=1
		off=1
		do irrep=1,nirrep
			l=ibas(irrep)
			csym(ofs:ofs+l*l-1)=reshape(cfull(off:off+l,off:off+l),(/l*l/))
			off=off+l
			ofs=ofs+l*l
		end do
	end subroutine
end module

!********************************************************************** 
!
! Set up the molecular orbitals 
! 
!********************************************************************** 
module mocoef
    use precision_m
	use cpdcom
	use sympack
	implicit none
	
	public init_mos, del_mos, getmos
	
	private
	real(DP), dimension(:,:), pointer :: mos
	logical :: debug=.false.
contains

	subroutine init_mos
		real(DP), dimension(:,:), allocatable :: scfevca0
		real(DP), dimension(:), allocatable :: foo

		allocate(mos(nbas,nbas))

		call getrec(20,'JOBARC','SCFEVCA0',nbas2*iintfp,mos) 
		! reformat (symmetry pack) mo coefficients (not needed right now...)
		! call symc(scfevca0,foo,ibas)

		if (debug) then
			print *, ' MO coefficients C(mu,p)' 
			call output(mos,1,nbas,1,nbas,nbas,nbas,1)
			print *
		end if 

	end subroutine
	
	subroutine del_mos
		deallocate(mos)
	end subroutine

	subroutine getmos(c)
		real(DP), dimension(:,:), pointer :: c
		c=>mos
	end subroutine

end module

!********************************************************************** 
!
! Construct the HF/correlated AO density matrix
!
!********************************************************************** 
module densmat
    use precision_m
	use cpdcom
	use mocoef
	implicit none
	
	public getdens
	
	private
	integer(I4), dimension(8) :: ip
	integer(I4), dimension(3) :: irrepb
	logical :: debug=.false.

contains
	subroutine getdens(dens)
		real(DP), dimension(:,:), intent(out) :: dens

		if (correlated) then
			call cdens(dens)  ! correlated cases
		else
			call hfdens(dens)    ! hf-scf case
		end if
	end subroutine

	subroutine hfdens(dens)
		real(DP), dimension(:,:), intent(out) :: dens

		call getrec(20,'JOBARC','SCFDENSA',nbas2*iintfp,dens)

		if (debug) then
			print *, ' HF-SCF density matrix d(mu,nu)'
			call output(dens,1,nbas,1,nbas,nbas,nbas,1)
			print *
		end if
	end subroutine
	
	subroutine cdens(dens)
		real(DP), dimension(:,:), intent(out) :: dens

		real(DP), dimension(:,:), allocatable :: tmp
		real(DP), dimension(:,:), pointer :: mos
		
		!  read total mo-density from jobarc
		call getrec(20,'JOBARC','RELDENSA',nbas2*iintfp,dens)

		if (debug) then
	  		print *, ' Correlated density matrix d(p,q)'
	  		call output(dens,1,nbas,1,nbas,nbas,nbas,1)
			print *
	 	end if
		
		allocate(tmp(nbas,nbas))

		! transform correlated (total) mo density to ao basis
		! MO--> AO transformation
		call getmos(mos)
		tmp=matmul(mos, dens)
		dens=matmul(tmp,transpose(mos))

		deallocate(tmp)

		if (debug) then
			print *, ' Correlated density matrix D(mu,nu)'
			call output(dens,1,nbas,1,nbas,nbas,nbas,1)
			print *
		end if
	end subroutine

end module


!********************************************************************** 
!
! This module handles the correlation contributions to the perturbed 
! densities.
!
!********************************************************************** 
module ccpdensmod
    use precision_m
	use cpdcom
	implicit none
	
	public readccpd
	private
	logical :: debug=.false.
contains

	subroutine readccpd(b, pdens)
		real(DP), dimension(:,:), intent(out) :: pdens
		integer(I4), intent(in) :: b ! magnetic perturbation
		
		real(DP), dimension(:,:), allocatable :: block
		real(DP), dimension(:), allocatable :: scr
		
		allocate(scr(nbas2))
		allocate(block(nbas,nbas))

		! add contribution from ddxdf .... ddoo, ddvv, ddvo
		call getpdc1(block,scr,b)
		pdens=block

		! add contribution from ddxdf1 ...  ddoo 
		! (perturbed canonical orbitals only) 
		if (canonical) then
			call getpdc2(block,scr,b)
			pdens=pdens+block
		end if
		
		! add contribution from ddxdf2 ... ddvo 
		call getpdc3(block,scr,b)
		pdens=pdens+block

		deallocate(scr)
		deallocate(block)

	end subroutine 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Unpack the PDENSCnm vector into a full matrix
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine unpackpd(vec, mat)
		real(DP), dimension(:), intent(in) :: vec
		real(DP), dimension(:,:), intent(inout) :: mat

		integer(I4) :: o, v, v1, v2, ov, oo, vv

		o=noca; v=nvrta
		v1=o+1
		v2=o+v
		oo=o*o; vv=v*v; ov=o*v

		if (debug) then
			call prvecr(vec,oo+vv+ov)
		end if
		
		! deal with occupied-occupied mat 
		mat(1:o,1:o)=reshape(vec(1:oo),(/noca,noca/))

		! deal with virtuial-virtual mat 
		mat(v1:v2,v1:v2)=reshape(vec(oo+1:oo+vv),(/nvrta,nvrta/))
		
		! deal with virtuial-occupied mats
		mat(v1:v2,1:o)=reshape(vec(oo+vv+1:oo+vv+ov),(/nvrta,noca/))
		mat(1:o,v1:v2)=-transpose(mat(v1:v1,1:o))

		if (debug) then
			call prvecr(mat,nbas2)
		end if

	end subroutine
	
	subroutine getpdc1(pdens, scr, b)
		real(DP), dimension(:,:), intent(inout) :: pdens
		real(DP), dimension(:), intent(out) :: scr
		integer(I4), intent(in) :: b

		integer(I4) :: nlength
		
		nlength=(noca**2+nvrta**2+nvrta*noca)*iintfp

		select case(b)
			case(1)
				call getrec(20,'JOBARC','PDENSC1X',nlength,scr)
			case(2)
				call getrec(20,'JOBARC','PDENSC1Y',nlength,scr)
			case(3)
				call getrec(20,'JOBARC','PDENSC1Z',nlength,scr)
		end select

		call unpackpd(scr, pdens)
	end subroutine

	subroutine getpdc2(pdens, scr, b)
		real(DP), dimension(:,:), intent(inout) :: pdens
		real(DP), dimension(:), intent(out) :: scr
		integer(I4), intent(in) :: b

		integer(I4) :: nlength

		nlength=(noca**2+nvrta**2+nvrta*noca)*iintfp

		select case(b)
			case(1)
				call getrec(20,'JOBARC','PDENSC2X',nlength,scr)
			case(2)
				call getrec(20,'JOBARC','PDENSC2Y',nlength,scr)
			case(3)
				call getrec(20,'JOBARC','PDENSC2Z',nlength,scr)
		end select

		call unpackpd(scr, pdens)
	end subroutine

	subroutine getpdc3(pdens, scr, b)
		real(DP), dimension(:,:), intent(inout) :: pdens
		real(DP), dimension(:), intent(out) :: scr
		integer(I4), intent(in) :: b

		integer(I4) :: nlength

		nlength=(noca**2+nvrta**2+nvrta*noca)*iintfp

		select case(b)
			case(1)
				call getrec(20,'JOBARC','PDENSC3X',nlength,scr)
			case(2)
				call getrec(20,'JOBARC','PDENSC3Y',nlength,scr)
			case(3)
				call getrec(20,'JOBARC','PDENSC3Z',nlength,scr)
		end select

		call unpackpd(scr, pdens)
	end subroutine

end module

!********************************************************************** 
!
! Set up the perturbed density marices
! 
!********************************************************************** 
module pdensmat
    use precision_m
	use cpdcom
	use mocoef
	use sympack
	use ccpdensmod
	implicit none
	intrinsic sqrt
	
	public init_pdens, del_pdens, getpdens, getpmos, psigma

	private
	real(DP), parameter :: FACT=26.62568098450085d+00
	real(DP), parameter :: HALFM=-.5d0
	integer(I4), dimension(3) :: irrepb
	real(DP), dimension(:,:), allocatable :: uai, uij
	real(DP), dimension(:,:), pointer :: dmo
	logical :: debug=.false.

contains
	
	subroutine init_pdens
		allocate(uai(nvrta,noca))
		allocate(uij(noca,noca))
		allocate(dmo(nbas,noca))

		! get magnetic field information			  
		call getrec(20,'JOBARC','NSYMPERT',27,ntpert)				     

		! irreps of b-field components
		irrepb(1)=iyzpert
		irrepb(2)=ixzpert
		irrepb(3)=ixypert

		! init counter for b-field perturbation in each irrep
	end subroutine

	subroutine del_pdens
		deallocate(uai)
		deallocate(uij)
		deallocate(dmo)
	end subroutine
	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! perturbed AO density for magnetic perturbation b
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine getpdens(b, pdens)
		integer(I4), intent(in) :: b
		real(DP), dimension(:,:), intent(out) :: pdens

		integer(I4), dimension(8), save :: ip=0
		integer(I4) :: irrepx
		real(DP) :: tmp

		! irrep of b-field component
		irrepx=irrepb(b)
		
		! increment perturbation counter
		ip(irrepx)=ip(irrepx)+1
		
		! read in U(a*,i)^X from moints file
		call getuai(ip, irrepx, uai)

		! read in D(i*,j)^X from moints file, and convert to U(i,j)^X
		call getuij(ip, irrepx, uij)

		! Construct perturbed occupied MO:s
		call mkpmos(uij, uai, dmo)

		! construct perturbed density matrix
		call mkpdens(dmo, pdens)

		! add electron-correlated contribution
		if (correlated) then
			call ccpdens(b, pdens)
		end if

	end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! perturbed MOS for magnetic perturbation b
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine getpmos(b, pmos)
		integer(I4), intent(in) :: b
		real(DP), dimension(:,:), intent(out) :: pmos

		integer(I4), dimension(8), save :: ip=0
		integer(I4) :: irrepx
		real(DP) :: tmp

		! irrep of b-field component
		irrepx=irrepb(b)
		
		! increment perturbation counter
		ip(irrepx)=ip(irrepx)+1
		
		! read in U(a*,i)^X from moints file
		call getuai(ip, irrepx, uai)

		! read in D(i*,j)^X from moints file, and convert to U(i,j)^X
		call getuij(ip, irrepx, uij)

		! Construct perturbed MO:s
		call mkpmos(uij, uai, pmos, nbas)
	end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! electron-correlation contribution to perturbed densities
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine ccpdens(b, pdens)
		integer(I4), intent(in) :: b
		real(DP), dimension(:,:), intent(inout) :: pdens

		real(DP), dimension(:,:), pointer :: mos
		real(DP), dimension(:,:), allocatable :: scr, tdens

		allocate(tdens(nbas,nbas))
		allocate(scr(nbas,nbas))

		call readccpd(b, tdens)
		call getmos(mos)

		! transform from mo to ao basis: C D' C^T
		scr=matmul(mos,tdens)
		tdens=matmul(scr,transpose(mos))

		! add correlation contribution to HF part
		pdens=pdens+tdens

		! antisymmetrize perturbed density matrix: 
		! d(mu,nu) = 1/2 (d(mu,nu) - d(nu,mu))
		pdens=0.5d0*(pdens-transpose(pdens))

		! print out total perturbed ao density matrix
		if (debug) then
			print *, ' Perturbed correlated density matrix', &
					'  dD(mu,nu)/dx'
			call output(pdens,1,nbas,1,nbas,nbas,nbas,1)
			print *
		end if

		deallocate(tdens, scr)

	end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Calculate paramagnetic contributions to sigma (for debugging).
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine psigma(pdens, b)
		real(DP), dimension(:,:), intent(in) :: pdens
		integer(I4), intent(in) :: b
		
		integer(I4) :: natoms, nuc, mu
		real(DP) :: ps
		real(DP), dimension(:), allocatable :: scr, scr2
		real(DP), dimension(:,:), allocatable :: dmdens
		integer(I4), dimension(8) :: irrbas

		allocate(dmdens(nbas,nbas))
		allocate(scr(nbas2))
		allocate(scr2(nbas2))

		irrbas=0
		irrbas(1)=nbas
		call getrec(20,'JOBARC','NREALATM',1,natoms)

		do nuc=1,3*natoms

			! get DH(MU,NU)/DM and expand to full square matrix
			call getlst(scr,nuc+3,1,1,1,102)
			call matexp(1,irrbas,scr,dmdens,1)

			if (debug) then
				print *, 'dh(mu,nu)/dm'
				call output(dmdens,1,nbas,1,nbas,nbas,nbas,1)
			end if
 
			! contraction   DH(MU,NU)/DM * Dd(MU,NU)/Db
			ps=0.d0
			do mu=1,nbas
				ps=ps+dot_product(dmdens(:,mu),pdens(:,mu))
			end do
			ps=FACT*ps
			print '(a,f12.6,a,i3,a,i3)', ' Paramagnetic sigma', &
				ps,' for B component', b, ' m component', nuc
		end do

		deallocate(dmdens, scr, scr2)

	end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! make perturbed occupied MOs:
! dC(mu,i)/dX =  sum_p C(mu,p) U(p,i)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine mkpmos(uij, uai, dmo, imos)
		real(DP), dimension(:,:), intent(in) :: uij, uai
		real(DP), dimension(:,:), intent(out) :: dmo
		integer(I4), optional :: imos
		
		real(DP) :: tmp
		integer(I4) :: i, mu, p, nmos
		real(DP), dimension(:,:), pointer :: mos

		call getmos(mos)
		
		if (.not.present(imos)) then
			nmos=noca
		else
			nmos=imos
		end if

		do i=1,nmos
			do mu=1,nbas
				tmp=0.d0
				do p=1,noca
					tmp=tmp+mos(mu,p)*uij(p,i)
				end do
				do p=1,nvrta
					tmp=tmp+mos(mu,p+noca)*uai(p,i)
				end do
				dmo(mu,i)=tmp
			end do
		end do

		if (debug) then 
			print *, ' Perturbed occupied MO coefficients dC(mu,i)/dx'
			call output(dmo,1,nbas,1,noca,nbas,noca,1)
			print *
		end if
	end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! construct perturbed density matrix:
! dD(nu,mu)/dx = sum_i [dC(nu,i)*/dx C(nu,i) + C(mu,i) dC(nu,i)/dx]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine mkpdens(dmo, pdens)
		real(DP), dimension(:,:), intent(in) :: dmo
		real(DP), dimension(:,:), intent(out) :: pdens
		
		real(DP) :: tmp
		integer(I4) :: i, mu, nu
		real(DP), dimension(:,:), pointer :: mos

		call getmos(mos)

		do mu=1,nbas
			do nu=1,nbas
				tmp=0.d0
				do i=1,noca
					tmp=tmp+2.d0*(mos(nu,i)*dmo(mu,i)-dmo(nu,i)*mos(mu,i))
				end do
				pdens(nu,mu)=tmp
			end do
		end do

		if (debug) then
			print *, 'Perturbed HF-SCF density matrix dD(mu,nu)/dx'
			call output(pdens,1,nbas,1,nbas,nbas,nbas,1)
		end if
	end subroutine
	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Get CPHF coefficients
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine getuai(ip, irrepx, uai)
		integer(I4), dimension(8), intent(in) :: ip
		integer(I4), intent(in) :: irrepx
		real(DP), dimension(:,:), intent(out) :: uai

		! read in U(a*,i)^X from moints file
		call getlst(uai,ip(irrepx),1,1,irrepx,182)

		if (debug) then
			print *, ' CPHF coefficients U(a,i)'
			call output(uai,1,nvrta,1,noca,nvrta,noca,1)
			print *
		end if
	end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Get CPHF coefficients
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine getuij(ip, irrepx, uij)
		integer(I4), dimension(8), intent(in) :: ip
		integer(I4), intent(in) :: irrepx
		real(DP), dimension(:,:), intent(out) :: uij

		! read in S(i*,j)^X from moints file
		call getlst(uij,ip(irrepx),1,1,irrepx,170)

		! convert S(i,j)^X to U(i,j)^X (factor -1/2) (not for perturbed 
		! canonical orbitals, iflags(64).ne.0)
		if (.not.canonical) then
			uij=HALFM*uij
		end if

		if (debug) then
			print *, ' CPHF coefficients U(j,i)'
			call output(uij,1,noca,1,noca,noca,noca,1)
			print *
		end if
	end subroutine
end module

!********************************************************************** 
!
! Quick and dirty program to get density matrix and perturbed
! density matriceS
!
! Original version: HF-SCF densities
! Modified version: correlated densities
! F90 version     : J. Juselius, July 2002
!
! TODO: So far no symmetry usage has been implemented
!       SPHERICAL -> CARTESIAN transformation
!
! J. Gauss University of Mainz, 2002
! 
!********************************************************************** 

program cpdens
    use precision_m
	use cpdcom
	use mocoef
	use densmat
	use pdensmat
	implicit none

	logical :: debug=.false.
	real(DP), dimension(1) :: scr

	integer(I4) :: b

	real(DP), dimension(:,:), allocatable :: dens, pdens, dmo
	real(DP), dimension(:,:), pointer :: mos
	real(DP), dimension(:), pointer :: eig

	! initialize aces2 module
	call crapsi(scr,b,0)

	! nbas number of basis functions
	nbas=noca+nvrta
	nbas2=nbas*nbas

	! output files
	open(42,file='XDENS')
	allocate(dens(nbas,nbas))
	
	if (.not.correlated) then
		open(43,file='XMOS')
		allocate(dmo(nbas,nbas))
	end if

	call setup_commons
	call init_mos()
	call init_pdens()
	
	! get unperturbed (electron-correlated) density matrix and MOs
	call getdens(dens)
	call getmos(mos)

	! write unperturbed density and MOs to disk

	write(42,'(d)') dens
	write(42,*) 
	if (.not.correlated) then
		write(43,*) mos
		write(43,*)
	end if
	
	allocate(pdens(nbas,nbas))

	if (debug) then
		print *, 'LOOP OVER MAGNETIC FIELD COMPONENTS'
		print *, '------------------------------------'
		print *
	end if
	
	! loop over magnetic field components ...
	do b=1,3

		if (debug) then
			print '(a,i3)', ' B component:', b
		end if

		! get perturbed (electron-correlated) density matrix for B comp b
		call getpdens(b, pdens)

		! write perturbed density and MOs to disk

!        write(42,*) pdens
		write(42,'(d)') pdens
		write(42,*) 

		! get perturbed MOs for B 
		if (.not.correlated) then
			!call getpmos(b, dmo)
			!write(43,*) dmo
			!write(43,*)
		end if

		! check perturbed density matrices by computing the 
		! paramagnetic shieldings
!        call psigma(pdens, b) 
		print *
	end do

	deallocate(dens, pdens)
	if (allocated(dmo)) deallocate(dmo)

	call del_mos
	call del_pdens

	close(42)
	if (.not.correlated) close(43)

	print *
	print *, 'Density and perturbed densities written in file ''XDENS'''
	if (.not.correlated) then
		print *, 'MOs and perturbed MOs written in file ''XMOS'''
	end if

	! all done, call crapso
	call crapso
end program

