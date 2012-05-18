module aces2_module
    integer, parameter :: ABEL=8

    public

contains

!
! THIS ROUTINE EXPANDS THE A COMPRESSED MATRIX A(P,Q)
! P >= Q TO AN ARRAY A(PQ) WITH P,Q. NOTE THIS ROUTINE 
! EXPECTS THAT THE ARRAY A IS SYMMETRY PACKED.
!
!  INPUT : IRREP  ...  THE IRREP OF THE CORRESPONDING PART OF A
!          NUM ......  POPULATION VECTOR FOR I AND J
!          DISSIZE ..  DISTRIBUTION SIZE OF A
!          A     ....  THE MATRIX A
!          IANTI ..... 0 FOR SYMMETRIC AND 1 FOR ANTISYMMETRIC
!                      MATRICES
!
!  OUTPUT : B .......  THE EXPANDED MATRIX A
!
!END
!
! CODED JG JAN/91
!
    subroutine matexp(irrep,num,a,b,ianti)
        implicit double precision (a-h,o-z)
        real(8), dimension(:) :: a, b
        real(8), dimension(ABEL) :: num
        dimension ipold(ABEL),ipnew(ABEL)
        common /syminf/nstart,nirrep,irrepa(255),irrepb(255), dirprd(ABEL,ABEL)

        integer :: dirprd

        ind(j,i)=(j*(j-1))/2+i

        data zero /0.0d0/
    !
    ! TAKE HERE CARE, IF WE ARE HANDLING IRREP=1 (TOTAL SYMMETRIC)
    ! OR IRREP=1 (OTHERWISE)
    !
        if(irrep.eq.1) then
    !
    ! GET FIRST POINTERS FOR OLD AND NEW INDICES
    !
            ipold(1)=0
            ipnew(1)=0
            do irrepj=1,(nirrep-1)
                ipold(irrepj+1)=ipold(irrepj)+(num(irrepj)*(num(irrepj)+1))/2
                ipnew(irrepj+1)=ipnew(irrepj)+num(irrepj)**2
            end do
    !
    ! NOW LOOP OVER ALL IRREPS
    !
            do  irrepj=1,nirrep
                numj=num(irrepj)
                ipo=ipold(irrepj)
                ipn=ipnew(irrepj)
    !
    ! LOOP OVER ORBITALS, BUT ALSO IN BACKWARD ORDER
    !
                if(ianti.eq.0) then
                    do  j=1,numj
                        do i=1,j
                            ind1=ind(j,i)+ipo
                            ind2=(j-1)*numj+i+ipn
                            ind3=(i-1)*numj+j+ipn
                            b(ind2)=a(ind1)
                            b(ind3)=a(ind1)
                        end do
                    end do
                else
                    do j=1,numj
                        do i=1,j
                            ind1=ind(j,i)+ipo
                            ind2=(j-1)*numj+i+ipn
                            ind3=(i-1)*numj+j+ipn 
                            b(ind2)=a(ind1)
                            b(ind3)=-a(ind1)
                        end do
                    end do
                endif
            end do
        else
    !
    ! FILL THE POINTERS OF THE OLD AND NEW ARRAY
    !
            ipold(1)=0
            ipnew(1)=0
            do irrepj=1,nirrep-1
                irrepi=dirprd(irrep,irrepj)
                numj=num(irrepj)
                numi=num(irrepi)
                ipnew(irrepj+1)=ipnew(irrepj)+numj*numi
                if(irrepi.lt.irrepj) then
                    ipold(irrepj+1)=ipold(irrepj)+numj*numi
                else
                    ipold(irrepj+1)=ipold(irrepj)
                endif
            end do
    !
    ! NOW COPY OLD ARRAYS TO NEW LOCATION
    !
            do irrepj=1,nirrep
                irrepi=dirprd(irrep,irrepj)
                numj=num(irrepj)
                numi=num(irrepi)
                if(irrepj.gt.irrepi) then
                    ipn=ipnew(irrepj)
                    ipo=ipold(irrepj)
                    do ij=1,numj*numi
                        ipnn=ipn+ij
                        ipoo=ipo+ij
                        b(ipnn)=a(ipoo)
                    end do
                else
                    if(ianti.eq.0) then
                        ipn=ipnew(irrepj)
                        ipo=ipold(irrepi)
                        do j=1,numj
                            do i=1,numi
                                ind1=(i-1)*numj+j+ipo
                                ind2=(j-1)*numi+i+ipn
                                b(ind2)=a(ind1)
                            end do
                        end do
                    else
                        ipn=ipnew(irrepj)
                        ipo=ipold(irrepi)
                        do j=1,numj
                            do i=1,numi
                                ind1=(i-1)*numj+j+ipo
                                ind2=(j-1)*numi+i+ipn
                                b(ind2)=-a(ind1)
                            end do
                        end do
                    endif
                endif
            end do
        endif
    end subroutine

!
! read the correlation contribution to the perturbed density from JOBARC
!
    subroutine readpd(ddens,scr,noca,nvrta,nbas,inumber,ib,irrepx, ispin)
        implicit double precision (a-h,o-z)
        integer pop,vrt,dirprd
        dimension ddens(*),scr(*),ioffo(8),ioffv(8)

        common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
        common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)
        common/syminf/nstart,nirrep,irreps(255,2),dirprd(8,8)
        common/sym/pop(8,2),vrt(8,2),nt(2),nf1(2),nf2(2)

        nlength=irpdpd(irrepx,21)+irpdpd(irrepx,19)+irpdpd(irrepx,9)

        ioff=0
        do irrep=1,nirrep
            ioffo(irrep)=ioff
            ioff=ioff+pop(irrep,ispin)
        end do
        do irrep=1,nirrep
            ioffv(irrep)=ioff
            ioff=ioff+vrt(irrep,ispin)
        end do

        if(inumber.eq.1) then
            if(ib.eq.1) call getrec(20,'JOBARC','PDENSC1X',nlength*iintfp,scr)
            if(ib.eq.2) call getrec(20,'JOBARC','PDENSC1Y',nlength*iintfp,scr)
            if(ib.eq.3) call getrec(20,'JOBARC','PDENSC1Z',nlength*iintfp,scr)
        else if(inumber.eq.2) then
            if(ib.eq.1) call getrec(20,'JOBARC','PDENSC2X',nlength*iintfp,scr)
            if(ib.eq.2) call getrec(20,'JOBARC','PDENSC2Y',nlength*iintfp,scr)
            if(ib.eq.3) call getrec(20,'JOBARC','PDENSC2Z',nlength*iintfp,scr)
        else if(inumber.eq.3) then
            if(ib.eq.1) call getrec(20,'JOBARC','PDENSC3X',nlength*iintfp,scr)
            if(ib.eq.2) call getrec(20,'JOBARC','PDENSC3Y',nlength*iintfp,scr)
            if(ib.eq.3) call getrec(20,'JOBARC','PDENSC3Z',nlength*iintfp,scr)
        endif
        !
        ! DEAL WITH OCCUPIED-OCCUPIED BLOCK 
        !
        ioff=0
        do irrepr=1,nirrep
            irrepl=dirprd(irrepx,irrepr)
            nocr=pop(irrepr,ispin)
            nocl=pop(irrepl,ispin)
            do i=1,nocr
                do j=1,nocl
                    ioff=ioff+1
                    index2=j+ioffo(irrepl)+(i+ioffo(irrepr)-1)*nbas
                    ddens(index2)=ddens(index2)+scr(ioff)
                end do
            end do
        end do
        !
        ! DEAL WITH VIRTUAL-VIRTUAL BLOCK 
        !
        do irrepr=1,nirrep
            irrepl=dirprd(irrepx,irrepr)
            nvrtr=vrt(irrepr,ispin)
            nvrtl=vrt(irrepl,ispin)
            do i=1,nvrtr
                do j=1,nvrtl
                    ioff=ioff+1
                    index2=j+ioffv(irrepl)+(i+ioffv(irrepr)-1)*nbas
                    ddens(index2)=ddens(index2)+scr(ioff)
                end do
            end do
        end do
        !
        ! DEAL WITH VIRTUAL-OCCUPIED BLOCK 
        !
        do irrepr=1,nirrep
            irrepl=dirprd(irrepx,irrepr)
            nocr=pop(irrepr,ispin)
            nvrtl=vrt(irrepl,ispin)
            do i=1,nocr
                do j=1,nvrtl
                    ioff=ioff+1
                    index2=j+ioffv(irrepl)+(i+ioffo(irrepr)-1)*nbas
                    index3=i+ioffo(irrepr)+(j+ioffv(irrepl)-1)*nbas
                    ddens(index2)=ddens(index2)+scr(ioff)
                    ddens(index3)=ddens(index3)-scr(ioff)
                end do
            end do
        end do
    end subroutine

!
!   THIS SUBROUTINE REORDERS A GIVEN MATRIX IN SUCH A WAY THAT
!   ABACUS (CALCULATOR AND IN THE FUTURE CRAY) CAN HANDLE THE
!   BACK TRANSFORMATION OF THE ONE-ELECTRON DENSITY MATRIX
!
!   A ... .... INPUT ARRAY
!   SCR ...... OUTPUT ARRAY
!   NBAST .... TOTAL NUMBER OF BASIS FUNCTIONS
!   POP ...... NUMBER OF OCCUPIED ORBITALS PER IRREP
!   VRT ...... NUMBER OF VIRTUAL ORBITALS PER IRREP
!   NBAS ..... NUMBER OF BASIS FUNCTIONS PER IRREP
!
!END
!
! CODED OCT/90 JG
!
    subroutine reordc(a,scr,nbast,pop,vrt,nbas)
        implicit double precision(a-h,o-z)
        integer dirprd,pop,vrt
        dimension a(nbast*nbast),scr(nbast*nbast)
        dimension pop(8),vrt(8),nbas(8)

        common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
        common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
        !
        ! COPY FIRT CFULL TO SCRATCH
        !
        call myicopy(a,scr,nbast*nbast*iintfp)
        !
        !  OFFSET FOR OCCUPIED AND VIRTUAL BLOCK
        !
        ioffo=1
        ioffv=1
        do irrep=1,nirrep
            ioffv=ioffv+pop(irrep)*nbast
        end do
        !
        !  OFFSET FOR TARGET ARRAY
        !
        iofft=1
        !
        !  LOOP OVER ALL IRREPS
        !
        do  irrep=1,nirrep
            !
            ! FILL FIRST WITH OCCUPIED ORBITALS OF THIS BLOCK
            !
            call myicopy(scr(ioffo),a(iofft),iintfp*pop(irrep)*nbast)

            ioffo=ioffo+nbast*pop(irrep)
            iofft=iofft+nbast*pop(irrep)
            !
            ! FILL NOW WITH VIRTUAL ORBITALS OF THIS BLOCK
            !
            print *, iofft,iintfp*vrt(irrep)*nbast,&
            iofft+iintfp*vrt(irrep)*nbast 
            call myicopy(scr(ioffv),a(iofft),iintfp*vrt(irrep)*nbast)

            ioffv=ioffv+nbast*vrt(irrep)
            iofft=iofft+nbast*vrt(irrep)
        end do
    end subroutine

!
!  THIS SUBROUTINE RETURNS THE SYMMETRY PACKED LIST OF THE INPUT
!  MATRIX CFULL. IT IS USED IN THE INTEGRAL DERIVATIVE CODE TO
!  SYMMETYRY PACK THE EIGEN VECTOR MATRIX.
!
!
!  CFULL ..... INPUT MATRIX ( SIZE NBAST*NBAST)
!  CSYM ...... SYMMETRY PACKED OUTPUT MATRIX (SQUARE MATRICES WITHIN
!                                             EACH IRREP)
!  NBAST ..... TOTAL NUMBER OF BASIS FUNCTION
!  NBAS ...... NUMBER OF BASIS FUNCTION PER IRREP( ARRAY OF DIM 8)
!  SCF ....... FLAG WHICH TELLS IF THE EIGENVECTORS HAVE BEEN REORDERED
!              OR NOT.
!              SCF = .TRUE.   NO REORDERING REQUIRED HERE
!                  = .FALSE.  REORDERING IS REQUIRED
!  ISPIN ..... SPIN CASE (AS USUAL)
!
!END
!
!   CODED OCT/90   JG
!
    subroutine symc(cfull,csym,nbast,nbas,scf,ispin)
        implicit double precision(a-h,o-z)
        integer dirprd,pop,popfull,vrt
        logical scf
        dimension cfull(nbast,nbast),csym(*),nbas(8),popfull(8)

        common/flags/iflags(100)
        common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
        common/sym/pop(8,2),vrt(8,2),nt(2),nf1(2),nf2(2)
        !
        !  FOR CORRELATION METHODS REORDER C SO THAT WE HAVE ALL EIGENVECTORS
        !  OF THE SAME IRREP TOGETHER
        !
        !print *, '!!!!check line 327'
        call reordc(cfull,csym,nbast,pop(1,ispin),vrt(1,ispin),nbas)

        ioffr=1
        ioffc=0
        ioffp=1

        do irrep=1,nirrep
            do imo=1,nbas(irrep)
                ioffc=ioffc+1
                length=nbas(irrep)
                call scopy(length,cfull(ioffr,ioffc),1,csym(ioffp),1)
                ioffp=ioffp+length
            end do
            ioffr=ioffr+nbas(irrep)
        end do
    end subroutine

!
!  THIS SUBROUTINE RETURNS THE SYMMETRY PACKED LIST OF THE INPUT
!  MATRIX CFULL. IT IS USED IN THE INTEGRAL DERIVATIVE CODE TO
!  SYMMETYRY PACK THE EIGEN VECTOR MATRIX.
!
!
!  CFULL ..... OUTPUT MATRIX ( SIZE NBAST*NBAST)
!  CSYM ...... SYMMETRY PACKED INPUT MATRIX (SQUARE MATRICES WITHIN
!                                             EACH IRREP)
!  NBAST ..... TOTAL NUMBER OF BASIS FUNCTION
!  NBAS ...... NUMBER OF BASIS FUNCTION PER IRREP( ARRAY OF DIM 8)
!  SCF ....... FLAG WHICH TELLS IF THE EIGENVECTORS HAVE BEEN REORDERED
!              OR NOT.
!              SCF = .TRUE.   NO REORDERING REQUIRED HERE
!                  = .FALSE.  REORDERING IS REQUIRED
!  ISPIN ..... SPIN CASE (AS USUAL)
!  IRREPX .... IRREP OF ARRAY TO BE UNPACKED
!END
!
!   CODED OCT/90   JG
!
    subroutine symc2(cfull,csym,nbast,nbas,scf,ispin,irrepx)
        implicit double precision(a-h,o-z)
        integer dirprd,pop,popfull,vrt
        logical scf
        dimension cfull(nbast,nbast),csym(*),nbas(8),popfull(8)

        common/flags/iflags(100)
        common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
        common/sym/pop(8,2),vrt(8,2),nt(2),nf1(2),nf2(2)

        ioffr=1
        ioffc=0
        ioffp=1

        call zero(cfull,nbast*nbast)
        do irrep=1,nirrep
            do imo=1,nbas(irrep)
                ioffc=ioffc+1
                irrep2=dirprd(irrep,irrepx)
                length=nbas(irrep2)
                ioffr=1
                do irrep1=1,irrep2-1
                    ioffr=ioffr+nbas(irrep1)
                end do
                call scopy(length,csym(ioffp),1,cfull(ioffr,ioffc),1)
                ioffp=ioffp+length
            end do
        end do
    end subroutine
end module

program xcpdens_sym
!
! QUICK AND DIRTY PROGRAM TO GET DENSITY MATRIX AND PERTURBED
! DENSITY MATRICES
!
! ORIGINAL VERSION: HF-SCF DENSITIES
! MODIFIED VERSION: CORRELATED DENSITIES
!
!END
!
! J. GAUSS, UNIVERSITY OF MAINZ, 2002
!
! Open-shells and F90 version 
! J. Juselius, University of Tromsø, 2007
!
    use aces2_module
    implicit double precision(a-h,o-z)

#ifdef MPI
    include 'mpif.h'
#endif
    integer, parameter :: maxcor=2000000, alpha=1, beta=2
    logical :: debug, magn
    real(8), dimension(maxcor) :: scr
    integer, dimension(3) :: irrepb
    integer, dimension(8) :: ioffmo,ioffdmo, nbasi, ip
    integer :: dirprd,pop,vrt
    real(8), dimension(:,:), pointer :: dmat

    common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
    common/syminf/nstart,nirrep,irreps(255,2),dirprd(8,8)
    common/sym/pop(8,2),vrt(8,2),nt(2),nf1(2),nf2(2)
    common/sympop/irpdpd(8,22),isytyp(2,500),id(18)
    common/symloc/isymoff(8,8,25) 
    common/flags/iflags(100)
    common/info/noca,nocb,nvrta,nvrtb
    common/pert/ntpert,npert(8),ipert(8),ixpert,iypert,izpert, &
    iyzpert,ixzpert,ixypert,itransx,itransy,itransz,&
    nucind                                                        

    data fact/26.62568098450085d+00/
    data halfm/-.5d0/

    integer :: uhf, ispin

    ! INIT ACES2 MODULE
#ifdef MPI
    call initenv('MODULE')
#endif
    call crapsi(scr,iuhf,0)
    debug=.false.

    if (debug) print *, '/machsp/', iintfp, iintln

    select case(iflags(18))
        case(3,4,5,6,15,16,17)
            magn=.true.
        case default
            magn=.false.
    end select
    
    if (magn) then
        open(42, file='XDENS', status='unknown')
    else
        open(42, file='CAODENS', status='unknown')
    end if

    do ispin=1,iuhf+1
        call gendens(ispin)
    end do
    close(42)

    if (.not.magn) then
        nbas=noca+nvrta
        allocate(dmat(nbas,nbas))
        call wrtltdm(dmat)
        deallocate(dmat)
    end if


    ! ALL DONE, CALL CRAPSO
    call crapso()
    call finalize()
contains
!
!  Write the lower triangular CAO density matrix
!
    subroutine wrtltdm(dmat)
        real(8), dimension(:,:), intent(out) :: dmat

        integer :: i, j, nbas
        
        nbas=size(dmat(:,1))
        
        open(42, file='CAODENS', status='unknown')
        read(42,*) dmat
        close(42)

        open(42, file='CAODENS2', status='unknown')
        write(42,*) nbas
        do j=1,nbas
            do i=j,nbas
                write(42,*) dmat(i,j)
            end do
        end do
        close(42)
    end subroutine
    
    subroutine gendens(ispin)
        integer, intent(in) :: ispin
!
! NBAS NUMBER OF BASIS FUNCTIONS
!
        if (ispin == 1) then
            nbas=noca+nvrta
        else
            nbas=nocb+nvrtb
        end if
        nbas2=nbas*nbas
!
! GET NUMBER OF BASIS FUNCTIONS CONSIDERING A CARTESIAN BASIS
!
        ione=1
        CALL GETREC(20,'JOBARC','NAOBASFN',IONE,NBASC)
        if (debug) print *, 'NBAS[C]', NBAS, NBASC
!
! SYMMETRY STUFF, NUMBER OF BASIS FUNCTIONS PER IRREP 
!
        do irrep=1,nirrep
            nbasi(irrep)=pop(irrep,ispin)+vrt(irrep,ispin)
        end do
!
! UNPERTURBED HF-SCF DENSITY MATRIX ...
!******************************************************
        idens=1
        istart=idens+nbasc*nbasc
!
! READ FROM JOBARC 
!
        if(iflags(2).eq.0) then
!
! HF-SCF CASE
!
!  READ AO-DENSITY FROM JOBARC
!
            if (ispin == 1) then
                call getrec(20,'JOBARC','SCFDENSA',nbas2*iintfp,scr(idens))
            else
                call getrec(20,'JOBARC','SCFDENSB',nbas2*iintfp,scr(idens))
            end if
!
! SAO --> CAO TRANSFORMATION; TWO XGEMM CALLS
!
            istart2=istart
            istart3=istart2+nbas*nbasc
            call getrec(20,'JOBARC','CMP2ZMAT',nbas*nbasc*iintfp,&
            scr(istart2))
            call xgemm('N','N',nbasc,nbas,nbas,1.d0,scr(istart2), &
            nbasc,scr(idens),nbas,0.d0,scr(istart3),nbasc)
            call xgemm('N','T',nbasc,nbasc,nbas,1.d0,scr(istart3),nbasc,&
            scr(istart2),nbasc,0.d0,scr(idens),nbasc)

!
! PRINT OUT AO HF-SCF DENSITY MATRIX ...
! **********************************
            if(debug) then
                write(*,*) ' HF-SCF density matrix D(mu,nu)'
                call output(scr(idens),1,nbasc,1,nbasc,nbasc,nbasc,1)
            endif

            atol=1.d-32
            do ijj=idens,idens+nbasc*nbasc-1
                if (dabs(scr(ijj)) < atol) then
                    write(42,*) 0.d0
                else
                    write(42,*) scr(ijj)
                end if
            end do
            write(42,*) 


        else
!
! CORRELATED UNPERTURBED DENSITY MATRIX
! *************************************
!
!  READ TOTAL MO-DENSITY FROM JOBARC
!
            if (ispin == 1) then
                call getrec(20,'JOBARC','RELDENSA',nbas2*iintfp,scr(idens))
            else
                call getrec(20,'JOBARC','RELDENSB',nbas2*iintfp,scr(idens))
            endif
!
! PRINT OUT MO DENSITY MATRIX ...
!
            if(debug) then
                write(*,*) ' Correlated density matrix D(p,q)'
                call output(scr(idens),1,nbas,1,nbas,nbas,nbas,1)
            endif
!
! TRANSFORM CORRELATED (TOTAL) MO DENSITY TO AO BASIS
!
! MEMORY ALLOCATION FOR MOS (KEEP ALSO FOR LATER USE, ADDRES IMO0 !)
! ISTART2 USED FOR SCRATCH
!
            imo=istart
            imo0=imo
            istart=imo+nbas*nbas
            istart2=istart

            if (ispin == 1) then
                call getrec(20,'JOBARC','SCFEVCA0',nbas2*iintfp,scr(imo)) 
            else
                call getrec(20,'JOBARC','SCFEVCB0',nbas2*iintfp,scr(imo)) 
            endif

!
! REFORMAT MO COEFFICIENTS
!
            if(debug) then
                write(*,*) ' MO coefficients c(mu,p)' 
                call output(scr(imo),1,nbas,1,nbas,nbas,nbas,1)
            endif 
!
! MO--> (S)AO TRANSFORMATION: TWO XGEMM CALLS
!
            call xgemm('n','n',nbas,nbas,nbas,1.d0,scr(imo),nbas,scr(idens),&
            nbas,0.d0,scr(istart2),nbas)
            call xgemm('n','t',nbas,nbas,nbas,1.d0,scr(istart2),nbas,&
            scr(imo),nbas,0.d0,scr(idens),nbas)
!        write(*,*) ' Correlated density matrix D(mu,nu)'
!        call output(scr(idens),1,nbas,1,nbas,nbas,nbas,1)
!
! SAO --> CAO TRANSFORMATION; TWO XGEMM CALLS
!
            istart3=istart2+nbasc*nbas
            call getrec(20,'JOBARC','CMP2ZMAT',nbas*nbasc*iintfp,&
            &               scr(istart2))
            call xgemm('n','n',nbasc,nbas,nbas,1.d0,scr(istart2),&
            nbasc,scr(idens),nbas,0.d0,scr(istart3),nbasc)
            call xgemm('n','t',nbasc,nbasc,nbas,1.d0,scr(istart3),nbasc,&
            scr(istart2),nbasc,0.d0,scr(idens),nbasc)
!
! PRINT (TOTAL) CORRELATED AO DENSITY MATRIX
! ******************************************
!
            if(debug) then
                write(*,*) ' Correlated density matrix D(mu,nu)'
                call output(scr(idens),1,nbasc,1,nbasc,nbasc,nbasc,1)
            endif

            atol=1.d-32
            do ijj=idens,idens+nbasc*nbasc-1
                if (dabs(scr(ijj)) < atol) then
                    write(42,*) 0.d0
                else
                    write(42,*) scr(ijj)
                end if
            end do
            write(42,*) 

        endif

!! return here if we don't have perturbed densities
    if (.not.magn) return

!
! PERTURBED DENSITY MATRIX .... MORE TO DO
! ****************************************
!
!  GET FIRST UNPERTURBED MO COEFFICIENTS
!
        imo=istart
        istart=imo+nbas*nbas

        if (ispin == 1) then
            call getrec(20,'JOBARC','SCFEVCA0',nbas2*iintfp,scr(istart)) 
        else
            call getrec(20,'JOBARC','SCFEVCB0',nbas2*iintfp,scr(istart)) 
        endif

!
! REFORMAT MO COEFFICIENTS
!
        call symc(scr(istart),scr(imo),nbas,nbasi,.false.,ispin)
        if(debug) then
            write(*,*) ' MO coefficients c(mu,p)'
            !         CALL OUTPUT(SCR(IMO),1,NBAS,1,NBAS,NBAS,NBAS,1)
            call prvecr(scr(imo),nbas*nbas)
        endif
!
! GET MAGNETIC FIELD INFORMATION                    
!
!        print *, '/kusse/', ntpert
        call getrec(20,'JOBARC','NSYMPERT',27,ntpert)                             
!
! IRREPS OF B-FIELD COMPONENTS
!
        irrepb(1)=iyzpert
        irrepb(2)=ixzpert
        irrepb(3)=ixypert
!
! COUNTER FOR B-FIELD PERTURBATION IN EACH IRREP
!
        do irrep=1,nirrep
            ip(irrep)=0
        end do
!
!   LOOP OVER MAGNETIC FIELD COMPONENTS ...
!   ****************************************
!
        istart0=istart

        do ib=1,3
!
! IRREP OF B-FIELD COMPONENT
!
            irrepx=irrepb(ib)
!
! INCREMENT PERTURBATION COUNTER
!
            ip(irrepx)=ip(irrepx)+1
!
! MEMORY FOR CPHF COEFFICIENTS AND OVERLAP MATRIX DERIVATIVES
!
            iuai=istart0
            isij=iuai+irpdpd(irrepx,9)
            istart=isij+irpdpd(irrepx,21)
! 
! READ IN U(A*,I)^x FROM MOINTS FILE
!     
!            print *, '!!! check getlst() 182'
            if (ispin == 1) then
                call getlst(scr(iuai),ip(irrepx),1,1,irrepx,182)
            else
                call getlst(scr(iuai),ip(irrepx),1,1,irrepx,183)
            end if
!
! READ IN S(I*,J)^x FROM MOINTS FILE
!
!            print *, '!!! check getlst() 170'
            if (ispin == 1) then
                call getlst(scr(isij),ip(irrepx),1,1,irrepx,170)
            else
                call getlst(scr(isij),ip(irrepx),1,1,irrepx,171)
            end if
!
! CONVERT S(I,J)^x TO U(I,J)^x (FACTOR -1/2) (NOT FOR PERTURBED 
! CANONICAL ORBITALS, IFLAGS(64).NE.0)
!
            if(iflags(64).eq.0) then
                call sscal(irpdpd(irrepx,21),halfm,scr(isij),1)
            endif
!
! MAKE PERTURBED OCCUPIED MOs
!
!  dc(mu,i)/dx =  sum_p c(mu,p) U(p,i)
! 
            idmo=istart
            istart=idmo+nbas*noca

            ioff1=idmo
            ioff2=imo
            do irrep=1,nirrep
                ioffdmo(irrep)=ioff1
                ioffmo(irrep)=ioff2
                ioff1=ioff1+pop(irrep,ispin)*nbasi(dirprd(irrep,irrepx))
                ioff2=ioff2+nbasi(irrep)*nbasi(irrep)
            end do

            do irrep=1,nirrep

                irrepmu=irrep
                irrepp=irrep
                irrepi=dirprd(irrepx,irrep)

                noci=pop(irrepi,ispin)
                nb=nbasi(irrepmu)
                nocp=pop(irrepp,ispin)
                nvrtp=vrt(irrepp,ispin)

                ioffu=iuai-1+isymoff(irrepi,irrepx,9)
                ioffs=isij-1+isymoff(irrepi,irrepx,21)
                ioffdm=ioffdmo(irrepi)
                ioffm=ioffmo(irrepp)

                do  i=1,noci 
                    do  imu=1,nb
                        index=ioffdm-1+imu+(i-1)*nb
                        sum=0.d0
                        do  m=1,nocp
                            index1=ioffm-1+imu+(m-1)*nb
                            index2=ioffs-1+m+(i-1)*nocp
                            sum=sum+scr(index1)*scr(index2)
                        end do
                        do  ie=1,nvrtp
                            index1=ioffm-1+imu+(ie+nocp-1)*nb
                            index2=ioffu-1+ie+(i-1)*nvrtp
                            sum=sum+scr(index1)*scr(index2)
                        end do
                        scr(index)=sum
                    end do
                end do
            end do

            if(debug) then 
                write(*,*) ' Perturbed MO coefficients dC(mu,i)/dx'
                call output(scr(idmo),1,nbas,1,noca,nbas,noca,1)
            endif
!
! CONSTRUCT PERTURBED DENSITY MATRIX
!
!  dD(nu,mu)/dx = sum_i [dc(mu,i)*/dx c(nu,i) + c(mu,i) dc(nu,i)/dx]
!
            iddens=istart
            istart=iddens+nbasc*nbasc
            ioffdd=iddens
            do  irrep=1,nirrep
                irrepmu=irrep
                irrepnu=dirprd(irrep,irrepx)
                nbasmu=nbasi(irrepmu)
                nbasnu=nbasi(irrepnu)
                nocmu=pop(irrepmu,ispin)
                nocnu=pop(irrepnu,ispin)
                do imu=1,nbasmu
                    do inu=1,nbasnu
                        index=ioffdd+inu-1+(imu-1)*nbasnu
                        sum=0.d0
                        do i=1,nocnu
                            index2=ioffmo(irrepnu)+inu-1+(i-1)*nbasnu
                            indexd1=ioffdmo(irrepnu)+imu-1+(i-1)*nbasmu
                            sum=sum+2.d0*(scr(index2)*scr(indexd1))
                        end do
                        do i=1,nocmu
                            index1=ioffmo(irrepmu)+imu-1+(i-1)*nbasmu
                            indexd2=ioffdmo(irrepmu)+inu-1+(i-1)*nbasnu
                            sum=sum-2.d0*(scr(indexd2)*scr(index1))
                        end do
                        scr(index)=sum
                    end do
                end do
                ioffdd=ioffdd+nbasmu*nbasnu
            end do
!
! SAO --> CAO TRANSFORMATION; TWO XGEMM CALLS
!
            call symc2(scr(istart),scr(iddens),nbas,nbasi,.false.,ispin,irrepx)

            call scopy(nbas*nbas,scr(istart),1,scr(iddens),1)

            istart2=istart
            istart3=istart2+nbasc*nbas
            call getrec(20,'JOBARC','CMP2ZMAT',nbas*nbasc*iintfp,&
            scr(istart2))
            call xgemm('n','n',nbasc,nbas,nbas,1.d0,scr(istart2),&
            nbasc,scr(iddens),nbas,0.d0,scr(istart3),nbasc)
            call xgemm('n','t',nbasc,nbasc,nbas,1.d0,scr(istart3),nbasc,&
            scr(istart2),nbasc,0.d0,scr(iddens),nbasc)
!
! PERTURBED HF-SCF DENSITY MATRIX
! *******************************
!
            if(iflags(2).eq.0.or.debug) then
            write(*,*) ' final perturbed HF-SCF density matrix dD(mu,nu)/dx'
!                call output(scr(iddens),1,nbasc,1,nbasc,nbasc,nbasc,1)
            endif

            if (iflags(2).eq.0) then
                do ijj=iddens,iddens+nbasc*nbasc-1
                    if (dabs(scr(ijj)) < atol) then
                        write(42,*) 0.d0
                    else
                        write(42,*) scr(ijj)
                    end if
                end do
                write(42,*)
            endif

            if(iflags(2).ne.0) then
!
! ADD ELECTRON-CORRELATED CONTRIBUTION
!
! READ CONTRIBUTION FROM DDXDF .... DDOO, DDVV, DDVO
!
                iscr1=istart
                iscr2=iscr1+nbasc*nbasc
                call zero(scr(iscr1),nbas*nbas)

                call readpd(scr(iscr1),scr(iscr2),noca,nvrta,nbas,1,&
                ib,irrepx, ispin)
!
! READ CONTRIBUTION FROM DDXDF1 ... DDOO (PERTURBED CANONICAL ORBITALS ONLY) 
!
                if(iflags(64).ne.0) then
                    call readpd(scr(iscr1),scr(iscr2),noca,nvrta,nbas,2,ib,&
                    irrepx, ispin)
                endif
!
! READ CONTRIBUTION FROM DDXDF2 ... DDVO 
!
                call readpd(scr(iscr1),scr(iscr2),noca,nvrta,nbas,3,ib,&
                irrepx, ispin)
!
!  TRANSFORM FROM MO TO AO BASIS: TWO XGEMM CALLS
!
                call xgemm('n','n',nbas,nbas,nbas,1.d0,scr(imo0),nbas,&
                    scr(iscr1), nbas,0.d0,scr(iscr2),nbas)
                call xgemm('n','t',nbas,nbas,nbas,1.d0,scr(iscr2),nbas,&
                    scr(imo0),nbas,0.d0,scr(iscr1),nbas)
                istart2=iscr2
                istart3=istart2+nbas*nbasc
                call getrec(20,'JOBARC','CMP2ZMAT',nbas*nbasc*iintfp,&
                    scr(istart2))
                call xgemm('n','n',nbasc,nbas,nbas,1.d0,scr(istart2),&
                    nbasc,scr(iscr1),nbas,0.d0,scr(istart3),nbasc)
                call xgemm('n','t',nbasc,nbasc,nbas,1.d0,scr(istart3),nbasc,&
                    scr(istart2),nbasc,0.d0,scr(iscr1),nbasc)
!
! ADD CORRELATION CONTRIBUTION TO HF PART
!
                call saxpy(nbasc*nbasc,1.d0,scr(iscr1),1,scr(iddens),1)
!
! ANTISYMMETRIZE PERTURBED DENSITY MATRIX: D(mu,nu) = 1/2 (D(mu,nu) - D(nu,mu))
!
                do  i=1,nbasc
                    do j=1,nbasc
                        index1=iddens-1+i+(j-1)*nbasc
                        index2=iddens-1+j+(i-1)*nbasc
                        a=(scr(index1)-scr(index2))/2.d0
                        scr(index1)=a
                        scr(index2)=-a
                    end do
                end do
!
! PERTURBED CORRELATED DENSITY MATRIX
! ***********************************
!
                if(debug) then
                    write(*,*) ' perturbed correlated density matrix',&
                    '  dD(mu,nu)/dx'
                    call output(scr(iddens),1,nbasc,1,nbasc,nbasc,nbasc,1)
                endif

                do ijj=iddens,iddens+nbasc*nbasc-1
                    if (dabs(scr(ijj)) < atol) then
                        write(42,*) 0.d0
                    else
                        write(42,*) scr(ijj)
                    end if
                end do
                write(42,*) 
            endif
!
! CHECK PERTURBED DENSITY MATRICES BY COMPUTING PARAMAGNETIC SHIELDINGS
!
!test         CALL GETREC(20,'JOBARC','NREALATM',1,NATOMS)
!test         DO 300 INUC=1,3*NATOMS
!
! GET dh(mu,nu)/dm AND EXPAND TO FULL SQUARE MATRIX
!
!test            CALL GETLST(SCR(ISTART),INUC+3,1,1,1,102)
!test            CALL MATEXP(1,NBAS,SCR(ISTART),SCR(ISTART+NBAS2),1)
!test            IF(DEBUG) THEN
!test               call output(scr(istart+nbas2),1,nbas,1,nbas,nbas,nbas,1)
!test            ENDIF
! 
! CONTRACTION   dh(mu,nu)/dm * dD(mu,nu)/dB
!
!test            A=FACT*SDOT(NBAS2,SCR(ISTART+NBAS2),1,SCR(IDDENS),1)
!test            write(*,*) ' paramagnetic sigma',a,' for B component',ib,&
!test     &                 ' m component',inuc
!test300      CONTINUE
        end do
    end subroutine
end program

! vim:et:sw=4:ts=4
