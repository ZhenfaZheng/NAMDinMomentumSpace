module epcoup
  use prec
  use fileio
  use couplings

  implicit none

  type epCoupling
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: epcq
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: epcr
  end type

  contains

  subroutine ReadEPC(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: ierr, line
    integer :: nbands, nkpts, nmodes, nqs
    integer :: iband, jband, ikpt, imode, iq
    character(len=255) :: epcdat, hash

    epcdat = "epmatwp.dat"
    
    open(unit=908, file=epcdat, status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "epmatwp.dat does NOT exist!"
      stop
    end if

    do line=1,2
      read(unit=908, fmt=*)
    end do
    read(unit=908, fmt=*) hash, nbands, nkpts, nmodes, nqs

    allocate(epc%epcq(nbands, nbands, nkpts, nmodes, nqs))
    do iband=1,nbands
      do jband=1,nbands
        do ikpt=1,nkpts
          do imode=1,nmodes
            do iq=1,nqs
              read(unit=908, fmt=*) epc%epcq(iband, jband, ikpt, imode, iq)
            end do
          end do
        end do
      end do
    end do

  end subroutine

  subroutine EPCq2r(epc)
    implicit none

    type(epCoupling), intent(inout) :: epc
  end subroutine

  subroutine TDepCoupIJ(inp, epc, olap)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap
  end subroutine

end module epcoup
