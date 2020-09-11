module epcoup
  use prec
  use fileio
  use couplings

  implicit none

  type epCoupling
    integer :: nbands, nkpts, nmodes, nqpts, natepc, natmd
    integer, allocatable, dimension(:,:) :: R
    !! cell numbers for each atoms of md cell
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
    real(kind=DP), allocatable, dimension(:,:) :: atompos

    filphmodes = 'graphene.matdyn.modes'

    open(unit=909, file=filphmodes, status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "phonon modes file does NOT exist!"
      stop
    end if

    nqs = 241
    nmodes = 6
    nat = 2
    naxis = 3 ! 3 dimension in xyz space.
    allocate(epc%phmodes(nqs, nmodes, nat, naxis))
    allocate(epc%qptsph(nqs, naxis))
    allocate(epc%freq(nqs, nmodes))
    allocate(atompos(nat, naxis))

    atompos(1,:) = (/0.0, 0.0, 0.5/)
    atompos(2,:) = (/0.333333333, 0.666666667, 0.5/)

    do iq=1,nqs
      read(unit=909, fmt=*)
      read(unit=909, fmt=*)
      read(unit=909, fmt=9019) charac, (epc%qptsph(iq, iaxis), &
                                                iaxis=1,naxis)
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

    real(kind=q) :: scal
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

    allocate(epc%displ(mdtime, nat, naxis))
    do time=1,mdtime
      read(unit=33, fmt=*)
      do iatom=1,nat
        read(unit=33, fmt=*) (epc%displ(time, iatom, iaxis), iaxis=1, naxis)
      end do
    end do

    allocate(epc%cellmd(nat+naxis, naxis))
    epc%cellmd(1:3,:) = supercell
    do iatom=1,nat
      do iaxis=1,naxis
        epc%cellmd(iatom+3,iaxis) = SUM(epc%displ(:,iatom, iaxis)) / mdtime
      enddo
    enddo

    do i=1,nat+3
      write(*,'(3f12.7)') (epc%cellmd(i,iaxis), iaxis=1,3)
    enddo

    close(33)

  end subroutine readDISPL

  subroutine cellPROJ(epc)
    ! Calculate cell numbers for each atom of md cell
    ! Project from phonon cell to md cell.
    implicit none

    integer :: naxis
    integer :: iaxis, iatom, jatom
    integer, allocatable, dimension(:) :: N
    !! scale numbers of md cell compared with phonon cell
    real(kind=DP), allocatable, dimension(:) :: dr
    real(kind=DP), allocatable, dimension(:) :: temp1, temp2

    type(epCoupling), intent(inout) :: epc

    naxis = 3
    allocate(R(epc%natmd,naxis), N(naxis), &
             dr(naxis), temp1(naxis), temp2(naxis))

    do iaxis=1,naxis
      N(iaxis) = NINT( SUM(epc%cellmd(iaixs,:)) / SUM(epc%cellepc(iaixs,:)) )
    enddo

    do iatom=1,epc%natmd
      temp1 = 9999.9
      do jatom=1,epc%natepc
        dr = epc%cellmd(iatom,:) * N - epc%cellepc(jatom,:)
        do iaxis=1,naxis
          temp2(iaxis) = ABS( dr(iaxis) - NINT(dr(axis)) )
        enddo
        if (SUM(temp2)<SUM(temp1)) then
          temp1 = temp2
          R(iatom,:) = (NINT(dr(i)), i=1,naxis)
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
    integer :: iq, imode, iatom, iaxis, time
    integer :: nqs, nmodes, nat, naxis

    call readPhmodes(inp, epc)
    call readDISPL(inp, epc)

    allocate(epc%phproj(inp%NSW, nmodes, nqs))

    do time=1,inp%NSW
      do iq=1,nqs
        do imode=1,nmodes
          Q = (0.0, 0.0)
          dQ = (0.0, 0.0)
          do iatom=1,nat
          Q = Q + dot_product( CONJG(epc%phmodes(iq, imode, iatom, :)), &
                               epc%displ(time, iatom, :) )
          dQ = dQ + dot_product( CONJG(epc%phmodes(iq, imode, iatom, :)), &
                                 epc%vel(time, iatom, :) )
          end do
          U = 0.5 * epc%freq(iq, imode)**2 * CONJG(Q) * Q
          T = 0.5 * CONJG(dQ) * dQ
          epc%phproj(time, iq, imode) = (U + T) / epc%freq(iq, imode)
        end do
      end do
    end do

  end subroutine phDecomp

  subroutine TDepCoupIJ(olap, olap_sec, inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap
    type(overlap), intent(inout) :: olap_sec

    real(kind=DP) :: proj
    integer :: iq, imode, iatom, iaxis
    integer :: nqs, nmodes, nat, naxis
    integer :: time, ik, jk, iband, jband
    complex(kind=q) :: temp

    ! Initialization
    olap%NBANDS = inp%NBANDS
    olap%TSTEPS = inp%NSW
    olap%dt = inp%POTIM
    allocate(olap%Dij(olap%NBANDS, olap%NBANDS, olap%TSTEPS-1))
    allocate(olap%Eig(olap%NBANDS, olap%TSTEPS-1))

    call readEPMat(inp, epc)
    call phDecomp(inp, epc)

    do time=1,inp%NSW-1
      do ik=1,inp%NKPOINTS
        do jk=1,inp%NKPOINTS
          do iband=1,inp%NBANDS
            do jband=1,inp%NBANDS
              do iq=1, nqs
                if (iq==(jk-ik)) then
                  temp = (0.0, 0.0)
                  do imode=1,nmodes
                    temp = temp + epc%phproj(time, iq, imode) * &
                           epc%epmat(iband, jband, imode, ik, iq)
                  end do
                  olap%Dij(iband*(ik-1), jband*(jk-1), time) = temp
                end if
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine TDepCoupIJ

end module epcoup
