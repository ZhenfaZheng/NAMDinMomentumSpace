module epcoup
  use prec
  use fileio
  use couplings

  implicit none

  type epCoupling
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: epmat
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: phmodes
    real(kind=DP), allocatable, dimension(:,:,:) :: displacement
  end type

  contains

  subroutine ReadEPMat(inp, epc)
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

  end subroutine ReadEPMat

  subroutine ReadPhmodes(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: ierr, i
    integer :: iq, imode, iatom, iaxis
    integer :: nqs, nmodes, nat, naxis
    character(len=255) :: filphmodes, fmtDISPL, line
    character(len=1) :: bra, ket

    filphmodes = 'graphene.matdyn.modes'

    open(unit=909, file=filphmodes, status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "phonon modes file does NOT exist!"
      stop
    end if

    nqs = 10
    nmodes = 6
    nat = 2
    naxis = 3
    allocate(epc%phmodes(nqs, nmodes, nat, naxis))

    do iq=1,nqs
      do i=1,3
        read(unit=909, fmt=*) line
      end do
      do imode=1,nmodes
        read(unit=909, fmt=*) line
        do iatom=1,nat
          read(unit=909, fmt=9021) bra, (epc%phmodes(iq, imode, iatom, iaxis), &
                                         iaxis=1,naxis), ket
        end do
      end do
      read(unit=909, fmt=*) line
    end do

  9020 format ( 1x, '(', 3(f10.6,1x,f10.6,3x), ')' )
  9021 format ( 1x, a, 3(f10.6,1x,f10.6,3x), a )

  end subroutine ReadPhmodes

  subroutine ReadDISPL(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
  end subroutine ReadDISPL

  subroutine TDepCoupIJ(inp, epc, olap)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap
  end subroutine TDepCoupIJ

end module epcoup
