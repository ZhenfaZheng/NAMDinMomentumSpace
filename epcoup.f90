module epcoup
  use prec
  use fileio
  use couplings

  implicit none

  type epCoupling
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: epmat
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: phmodes
    real(kind=DP), allocatable, dimension(:,:,:,:) :: phmodesR
    real(kind=DP), allocatable, dimension(:,:,:) :: displ
    ! Phonon projection of displacement in MD.
    real(kind=DP), allocatable, dimension(:,:,:) :: phproj
  end type

  contains

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
    character(len=255) :: filphmodes, fmtDISPL, line
    character(len=1) :: bra, ket
    complex(kind=DP) :: temp

    filphmodes = 'graphene.matdyn.modes'

    open(unit=909, file=filphmodes, status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "phonon modes file does NOT exist!"
      stop
    end if

    nqs = 10
    nmodes = 6
    nat = 2
    naxis = 3 ! 3 dimension in xyz space.
    allocate(epc%phmodes(nqs, nmodes, nat, naxis))
    allocate(epc%phmodesR(nqs, nmodes, nat, naxis))

    do iq=1,nqs
      do i=1,3
        read(unit=909, fmt=*) line
      end do
      do imode=1,nmodes
        read(unit=909, fmt=*) line
        do iatom=1,nat
          read(unit=909, fmt=9021) bra, (epc%phmodes(iq, imode, iatom, iaxis), &
                                         iaxis=1,naxis), ket
          do iaxis=1,naxis
            temp = epc%phmodes(iq, imode, iatom, iaxis)
            if (abs(real(temp)) >= abs(aimag(temp))) then
              epc%phmodesR(iq, imode, iatom, iaxis) &
              = abs(temp) * sign(1.0d0, real(temp))
            else
              epc%phmodesR(iq, imode, iatom, iaxis) &
              = abs(temp) * sign(1.0d0, aimag(temp))
            end if
          end do
        end do
      end do
      read(unit=909, fmt=*) line
    end do
    ! write(*,*) epc%phmodesR(1,1,:,:)

  9020 format ( 1x, '(', 3(f10.6,1x,f10.6,3x), ')' )
  9021 format ( 1x, a, 3(f10.6,1x,f10.6,3x), a )

  end subroutine readPhmodes

  subroutine readDISPL(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
  end subroutine readDISPL

  subroutine TDepCoupIJ(inp, epc, olap)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap

    real(kind=DP) :: proj
    integer :: iq, imode, iatom, iaxis, time
    integer :: nqs, nmodes, nat, naxis

    call readEPMat(inp, epc)
    call readPhmodes(inp, epc)
    call readDISPL(inp, epc)

    nqs = 10
    nmodes = 6
    nat = 2
    naxis = 3

    allocate(epc%phproj(inp%NSW, nmodes, nqs))

    do time=1,inp%NSW
      do iq=1,nqs
        do imode=1,nmodes
          proj = 0.0d0
          do iatom=1,nat
          proj = proj + dot_product( epc%phmodesR(iq, imode, iatom, :), &
                                     epc%displ(time, iatom, :))
          end do
          epc%phproj(time, imode, iq) = proj
        end do
      end do
    end do

  end subroutine TDepCoupIJ

end module epcoup
