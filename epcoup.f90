module epcoup
  use prec
  use constants
  use fileio
  use couplings

  implicit none

  type epCoupling
    integer :: nbands, nkpts, nmodes, nqpts, natepc, natmd
    integer, allocatable, dimension(:) :: atnum
    !! atom number of each md cell atom in phonon cell
    integer, allocatable, dimension(:,:) :: R
    !! cell numbers for each atoms of md cell
    integer, allocatable, dimension(:) :: qqmap
    !! map q points in phonon calculation to q points in e-p calcultion.
    integer, allocatable, dimension(:,:,:) :: kkqmap
    !! map ik & jk of electronic states to q or -q of e-p matrix.
    real(kind=DP), allocatable, dimension(:,:) :: cellepc
    real(kind=DP), allocatable, dimension(:,:) :: kptsepc
    real(kind=DP), allocatable, dimension(:,:) :: qptsepc
    real(kind=DP), allocatable, dimension(:,:) :: energy
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: epmat
    real(kind=DP), allocatable, dimension(:,:) :: freq
    real(kind=DP), allocatable, dimension(:,:) :: qptsph
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: phmodes
    real(kind=DP), allocatable, dimension(:,:) :: cellmd
    real(kind=DP), allocatable, dimension(:,:,:) :: displ
    real(kind=DP), allocatable, dimension(:,:,:) :: vel
    real(kind=DP), allocatable, dimension(:,:,:) :: phproj
    ! Phonon projection of displacement in MD.
  end type

  contains

  subroutine readEPC(inp, epc)
    ! Read informations about e-p couplings from .epc file.
    ! Including cell for epc, k, q points, eigen energies and e-p matrix.

    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: ierr, i, j
    integer :: nb, nk, nm, nq, na
    integer :: iband, jband, ikpt, imode, iqpt
    character(len=72) :: filepc
    
    filepc = 'graphene.epc'

    open(unit=90, file=filepc, status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "epc file does NOT exist!"
      stop
    end if

    read(unit=90, fmt=*)
    read(unit=90, fmt=*) nb, nk, nm, nq, na
    epc%nbands = nb
    epc%nkpts  = nk
    epc%nmodes = nm
    epc%nqpts  = nq
    epc%natepc = na
    allocate(epc%cellepc(na+3,3), &
             epc%kptsepc(nk,3), &
             epc%qptsepc(nq,3), &
             epc%energy(nk,nb), &
             epc%epmat(nb,nb,nk,nm,nq))

    read(unit=90, fmt=*)
    do i=1,3
      read(unit=90, fmt=*) (epc%cellepc(i,j), j=1,3)
    enddo
    read(unit=90, fmt=*)
    do i=1,na
      read(unit=90, fmt=*) (epc%cellepc(i+3,j), j=1,3)
    enddo

    read(unit=90, fmt=*)
    do i=1,nk
      read(unit=90, fmt=*) (epc%kptsepc(i,j), j=1,3)
    enddo

    read(unit=90, fmt=*)
    do i=1,nq
      read(unit=90, fmt=*) (epc%qptsepc(i,j), j=1,3)
    enddo

    read(unit=90, fmt=*)
    do i=1,nk
      read(unit=90, fmt=*) (epc%energy(i,j), j=1,nb)
    enddo

    read(unit=90, fmt=*)
    do iband=1,nb
      do jband=1,nb
        do ikpt=1,nk
          do imode=1,nm
            do iqpt=1,nq
              read(unit=90, fmt='(2f15.10)') epc%epmat(iband,jband,ikpt,imode,iqpt)
            end do
          end do
        end do
      end do
    end do

    close(90)

  end subroutine readEPC

  subroutine readEPMat(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: ierr, line
    integer :: nbands, nkpts, nmodes, nqs
    integer :: iband, jband, ikpt, imode, iq
    character(len=255) :: epmdat, hash

    epmdat = "epmatwp.dat"
    
    open(unit=908, file=epmdat, status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "epmatwp.dat does NOT exist!"
      stop
    end if

    do line=1,2
      read(unit=908, fmt=*)
    end do
    read(unit=908, fmt=*) hash, nbands, nkpts, nmodes, nqs

    allocate(epc%epmat(nbands, nbands, nkpts, nmodes, nqs))
    do iband=1,nbands
      do jband=1,nbands
        do ikpt=1,nkpts
          do imode=1,nmodes
            do iq=1,nqs
              read(unit=908, fmt=*) epc%epmat(iband, jband, ikpt, imode, iq)
            end do
          end do
        end do
      end do
    end do

  end subroutine readEPMat

  subroutine readPhmodes(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: ierr, i
    integer :: iq, imode, iatom, iaxis
    integer :: nqs, nmodes, nat, naxis
    character(len=255) :: filphmodes, fmtDISPL
    character(len=24) :: charac, bra, ket
    complex(kind=DP) :: temp
    real(kind=DP), allocatable, dimension(:) :: qpt
    real(kind=DP), allocatable, dimension(:,:) :: at

    filphmodes = 'graphene.matdyn.modes'

    open(unit=909, file=filphmodes, status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "phonon modes file does NOT exist!"
      stop
    end if

    nat = epc%natepc
    nmodes = epc%nmodes
    nqs = epc%nqpts
    naxis = 3 ! 3 dimension in xyz space.
    allocate(epc%phmodes(nqs, nmodes, nat, naxis))
    allocate(epc%qptsph(nqs, naxis))
    allocate(epc%freq(nqs, nmodes))
    allocate(qpt(naxis), at(naxis, naxis))

    at = epc%cellepc(1:naxis,1:naxis) / epc%cellepc(1,1)
    ! write(*,'(3f12.7)') ((at(i,iaxis),i=1,3), iaxis=1,3)

    do iq=1,nqs
      read(unit=909, fmt=*)
      read(unit=909, fmt=*)
      read(unit=909, fmt=9019) charac, (qpt(iaxis), iaxis=1,naxis)
      do iaxis=1,naxis
        epc%qptsph(iq,iaxis) = at(iaxis,1)*qpt(1) + at(iaxis,2)*qpt(2) &
                             + at(iaxis,3)*qpt(3)
      end do
      ! write(*, 9019) charac, epc%qptsph(iq,:)
      read(unit=909, fmt=*)
      do imode=1,nmodes
        read(unit=909, fmt=9011) charac, epc%freq(iq, imode)
        ! write(*,9011) charac, epc%freq(iq, imode)
        do iatom=1,nat
          read(unit=909, fmt=9021) bra, (epc%phmodes(iq, imode, iatom, iaxis), &
                                         iaxis=1,naxis), ket
          ! write(*, 9021) bra, (epc%phmodes(iq, imode, iatom, iaxis), &
          !                      iaxis=1,naxis), ket
        end do
      end do
      read(unit=909, fmt=*)
    end do

  9011 format ( 5x, a14, f15.6 )
  9019 format ( 1x, a4, 3f12.4 )
  9021 format ( 1x, a1, 3(f10.6,1x,f10.6,3x), a1 )

  end subroutine readPhmodes

  subroutine readDISPL(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    real(kind=q) :: scal, dr
    real(kind=q), allocatable, dimension(:,:) :: supercell
    character(24) :: filename
    integer :: ierr, i, iaxis, iatom, time
    integer :: nat, naxis, mdtime

    mdtime = 200
    filename = 'XDATCAR'

    open(unit=33, file=filename, status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "XDATCAR file does NOT exist!"
      stop
    end if

    naxis = 3
    allocate(supercell(naxis,naxis))

    read(unit=33, fmt=*)
    read(unit=33, fmt=*) scal
    do i=1,naxis
      read(unit=33, fmt=*) (supercell(i, iaxis), iaxis=1,naxis)
      supercell(i,:) = supercell(i,:) * scal
    end do
    read(unit=33, fmt=*)
    read(unit=33, fmt=*) nat
    epc%natmd = nat

    allocate(epc%displ(mdtime, nat, naxis), &
             epc%vel(mdtime, nat, naxis))
    do time=1,mdtime
      read(unit=33, fmt=*)
      do iatom=1,nat
        read(unit=33, fmt=*) (epc%displ(time, iatom, iaxis), iaxis=1, naxis)
        if (time==1) then
          epc%vel(time,iatom,:) = 0
        else
          ! Fix atom position if atom moves to other cell.
          do iaxis=1,naxis
            dr = epc%displ(time, iatom, iaxis) - epc%displ(time-1, iatom, iaxis)
            if (dr>0.9) then
              epc%displ(time, iatom, iaxis) = epc%displ(time, iatom, iaxis) - 1
            else if (dr<-0.9) then
              epc%displ(time, iatom, iaxis) = epc%displ(time, iatom, iaxis) + 1
            endif
          enddo
          epc%vel(time,iatom,:) &
          = ( epc%displ(time, iatom, :) - epc%displ(time-1, iatom, :) ) / inp%POTIM
        endif
      end do
    end do

    allocate(epc%cellmd(nat+naxis, naxis))
    epc%cellmd(:3,:) = supercell
    ! epc%cellmd(4:,:) = epc%displ(1,:,:)
    do iatom=1,nat
      do iaxis=1,naxis
        epc%cellmd(iatom+naxis, iaxis) = SUM(epc%displ(:, iatom, iaxis)) / mdtime
      end do
    end do

    do time=1,mdtime
      do iatom=1,nat
        epc%displ(time,iatom,:) &
        = epc%displ(time,iatom,:) - epc%cellmd(iatom+naxis, :)
      end do
    end do

    close(33)

  end subroutine readDISPL

  subroutine cellPROJ(epc)
    ! Calculate cell numbers for each atom of md cell
    ! Project from phonon cell to md cell.
    implicit none

    integer :: naxis
    integer :: iaxis, iatom, jatom, i
    integer, allocatable, dimension(:) :: N
    !! scale numbers of md cell compared with phonon cell
    real(kind=DP), allocatable, dimension(:) :: dr
    real(kind=DP), allocatable, dimension(:) :: temp1, temp2

    type(epCoupling), intent(inout) :: epc

    naxis = 3
    allocate(epc%R(epc%natmd,naxis), epc%atnum(epc%natmd), &
             N(naxis), dr(naxis), temp1(naxis), temp2(naxis))

    do iaxis=1,naxis
      N(iaxis) = NINT( SUM(epc%cellmd(iaxis,:)) / SUM(epc%cellepc(iaxis,:)) )
    enddo

    do iatom=1,epc%natmd
      temp1 = 9999.9
      do jatom=1,epc%natepc
        dr = epc%cellmd(iatom+naxis,:) * N - epc%cellepc(jatom+naxis,:)
        do iaxis=1,naxis
          temp2(iaxis) = ABS( dr(iaxis) - NINT(dr(iaxis)) )
        enddo
        if (SUM(temp2)<SUM(temp1)) then
          temp1 = temp2
          epc%atnum(iatom) = jatom
          epc%R(iatom,:) = (/(MOD(NINT(dr(i)),N(i)), i=1,naxis)/)
        endif
      enddo
    enddo

    deallocate(N, dr, temp1, temp2)

  end subroutine cellPROJ

  subroutine phDecomp(inp, epc)
    ! Project the atomic motion from an MD simulation to phonon modes.
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    ! The normal mode coordinate and its derivative.
    complex(kind=DP) :: Q, dQ
    ! The potential and kinetic energies of the normal mode.
    real(kind=DP) :: U, T
    real(kind=DP) :: theta
    real(kind=DP), allocatable, dimension(:) :: displ, vel
    integer :: iq, imode, iatom, iaxis, time
    integer :: nqs, nmodes, nat, naxis

    call readPhmodes(inp, epc)
    call readDISPL(inp, epc)
    call cellPROJ(epc)

    naxis = 3
    nmodes = epc%nmodes
    nqs = epc%nqpts
    allocate(epc%phproj(inp%NSW, nmodes, nqs), displ(naxis), vel(naxis))

    do time=1,inp%NSW
      do iq=1,epc%nqpts
        do imode=1,epc%nmodes
          Q = (0.0, 0.0)
          dQ = (0.0, 0.0)
          do iatom=1,epc%natmd
            displ = 0.0; vel = 0.0
            do iaxis=1,naxis
              displ = displ + epc%displ(time, iatom, iaxis) * epc%cellmd(iaxis, :) 
              vel = vel + epc%vel(time, iatom, iaxis) * epc%cellmd(iaxis, :) 
            end do
            theta = 2 * PI * DOT_PRODUCT( epc%qptsph(iq,:), epc%R(iatom,:) )
            Q  =  Q + EXP(imgUnit*theta) * DOT_PRODUCT( &
                      CONJG(epc%phmodes(iq, imode, epc%atnum(iatom), :)), displ )
            dQ = dQ + EXP(imgUnit*theta) * DOT_PRODUCT( &
                      CONJG(epc%phmodes(iq, imode, epc%atnum(iatom), :)), vel )
          end do
          ! mass of C = 12.011
          U = 0.5 * 12.011 * epc%freq(iq, imode)**2 * CONJG(Q) * Q / epc%natmd &
            * 1.66 * 6.2415 / 100000
          T = 0.5 * 12.011 * CONJG(dQ) * dQ / epc%natmd &
            * 1.66 * 6.2415 * 10
          epc%phproj(time, imode, iq) = (U + T) / ( hbar * epc%freq(iq, imode) / 1000 )
          !if (epc%phproj(time, imode, iq) > 5) then
          !write(*,'(3I4)') time,iq, imode
          !write(*,'(4f16.7)') Q,dQ
          !write(*,'(2f19.10)') U,T
          !write(*,'(2f15.10)') epc%phproj(time, imode, iq), epc%freq(iq, imode)
          !write(*,*)
          !end if
        end do
      end do
    end do

  end subroutine phDecomp

  subroutine kqMatch(epc)
    implicit none

    type(epCoupling), intent(inout) :: epc

    real(kind=DP) :: Norm
    real(kind=DP) :: dkq1(3), dkq2(3), dq(3)
    integer :: ik, jk, iq, jq, iaxis, naxis

    naxis = 3
    Norm = 0.001

    allocate(epc%kkqmap(epc%nkpts,epc%nkpts,2), epc%qqmap(epc%nqpts))

    epc%kkqmap = -1
    epc%qqmap = -1

    do ik=1,epc%nkpts
      do jk=1,epc%nkpts
        do iq=1,epc%nqpts
          dkq1 = epc%kptsepc(ik,:) - epc%kptsepc(jk,:) - epc%qptsepc(iq,:)
          dkq2 = epc%kptsepc(ik,:) - epc%kptsepc(jk,:) + epc%qptsepc(iq,:)
          do iaxis=1,naxis
            dkq1(iaxis) = ABS(dkq1(iaxis)-NINT(dkq1(iaxis)))
            dkq2(iaxis) = ABS(dkq2(iaxis)-NINT(dkq2(iaxis)))
          end do
          if (SUM(dkq1)<Norm) epc%kkqmap(ik,jk,1) = iq
          if (SUM(dkq2)<Norm) epc%kkqmap(ik,jk,2) = iq
          if (epc%kkqmap(ik,jk,1)>0 .and. epc%kkqmap(ik,jk,2)>0) exit
        end do
      end do
    end do

    do iq=1,epc%nqpts
      do jq=1,epc%nqpts
        dq = epc%qptsepc(iq,:) - epc%qptsph(jq,:)
        do iaxis=1,naxis
          dq(iaxis) = ABS(dq(iaxis)-NINT(dq(iaxis)))
        end do
        if (SUM(dq)<Norm) then
          epc%qqmap(iq) = jq
          exit
        end if
      end do
    end do

  end subroutine kqMatch

  subroutine TDepCoupIJ(olap, olap_sec, inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap
    type(overlap), intent(inout) :: olap_sec

    real(kind=DP) :: proj, norm, phn1, phn2
    real(kind=DP), allocatable, dimension(:) :: dkq, dq
    integer :: iq, jq1, jq2, imode, iatom, iaxis
    integer :: nqs, nmodes, nat, naxis
    integer :: time, ik, jk, iband, jband
    complex(kind=q) :: temp
    logical :: lcoup

    ! Initialization
    olap%NBANDS = inp%NBANDS * inp%NKPOINTS
    olap%TSTEPS = inp%NSW
    olap%dt = inp%POTIM
    allocate(olap%Dij(olap%NBANDS, olap%NBANDS, olap%TSTEPS-1))
    allocate(olap%Eig(olap%NBANDS, olap%TSTEPS-1))
    olap%Dij = cero
    olap%EIg = cero

    olap_sec%NBANDS = inp%NBASIS * inp%NKPOINTS
    olap_sec%TSTEPS = inp%NSW
    olap_sec%dt = inp%POTIM
    allocate(olap_sec%Dij(olap_sec%NBANDS, olap_sec%NBANDS, olap_sec%TSTEPS-1))
    allocate(olap_sec%Eig(olap_sec%NBANDS, olap_sec%TSTEPS-1))

    inquire(file='COUPCAR', exist=lcoup)
    if (lcoup) then
      ! file containing couplings exists, then read it
      if (inp%LCPTXT) then
        call readNaEig(olap_sec, inp)
      else
        call CoupFromFile(olap)
        call copyToSec(olap, olap_sec, inp)
        call writeNaEig(olap_sec, inp)
      end if
    else

      write(*,'(A)') "------------------------------------------------------------"
      write(*,*) "Calculating couplings from e-p matrix..."

      call readEPC(inp, epc)
      call phDecomp(inp, epc)
      call kqMatch(epc)
     
      naxis = 3
      norm = 0.01
      allocate(dkq(naxis), dq(naxis))
     
      do time=1,inp%NSW-1
        do iband=1,inp%NBANDS-0
          do jband=1,inp%NBANDS-0
            do ik=1,inp%NKPOINTS-0
              do jk=1,inp%NKPOINTS-0

                iq = epc%kkqmap(ik,jk,1)
                if (iq<0) exit

                temp = (0.0, 0.0)
                jq1 = epc%qqmap(epc%kkqmap(ik,jk,1)) ! ik-jk =  q
                jq2 = epc%qqmap(epc%kkqmap(ik,jk,2)) ! ik-jk = -q

                do imode=1,epc%nmodes
                  phn1 = ABS(epc%phproj(time, imode, jq1))
                  phn2 = ABS(epc%phproj(time, imode, jq2))
                  temp = temp + ( SQRT(phn1) + SQRT(phn2+1) ) &
                              * epc%epmat(iband, jband, jk, imode, iq)
                end do

                temp = temp / SQRT(1.0 * epc%natmd / epc%natepc) &
                       * 2 * inp%POTIM * imgUnit / hbar
                olap%Dij(inp%NBANDS*(ik-1)+iband, jband+inp%NBANDS*(jk-1), time) = temp

              end do
            end do
          end do
        end do
      end do
     
      do iband=1,inp%NBANDS
        do ik=1,inp%NKPOINTS
          olap%Eig(iband+inp%NBANDS*(ik-1),:) = epc%energy(ik,iband)
        enddo
      enddo

      write(*,*) "Done..."
      write(*,'(A)') "------------------------------------------------------------"
     
      call copyToSec(olap, olap_sec, inp)
      call CoupToFile(olap)
      call writeNaEig(olap_sec, inp)

    end if

  end subroutine TDepCoupIJ

  subroutine copyToSec(olap, olap_sec, inp)
    implicit none
    type(overlap), intent(inout) :: olap_sec
    type(overlap), intent(in) :: olap
    type(namdInfo), intent(in) :: inp

    integer :: ik, jk, NB, NBas

    NB = inp%NBANDS
    NBas = inp%NBASIS
    do ik=1,inp%NKPOINTS
      olap_sec%Eig((ik-1)*NBas+1:ik*NBas, :) = &
      olap%Eig((ik-1)*NB+inp%BMIN:(ik-1)*NB+inp%BMAX, :)
      do jk=1,inp%NKPOINTS
        olap_sec%Dij((ik-1)*NBas+1:ik*NBas,(jk-1)*NBas+1:jk*NBas, :) = &
        olap%Dij((ik-1)*NB+inp%BMIN:(ik-1)*NB+inp%BMAX, &
                 (jk-1)*NB+inp%BMIN:(jk-1)*NB+inp%BMAX, :)
      end do
    end do

  end subroutine copyToSec


end module epcoup
