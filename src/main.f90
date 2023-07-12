Program main
  use prec
  use lattice
  use wavecar
  use couplings
  use epcoup
  use hamil
  use shop
  use shotf
  use fileio
  use TimeProp
  use mpi

  implicit none

  type(namdInfo) :: inp
  type(TDKS) :: ks
  type(overlap) :: olap, olap_sec
  type(epCoupling) :: epc

  real(kind=q) :: start, fin
  integer :: irank, ierr
  integer :: t1, t2, rate
  integer :: ns

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)

  if (irank==0) write(*,*)
  if (irank==0) write(*,*) "Hefei-NAMD (epc version 2.1.15, Jul 13, 2023)"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First, get user inputs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getUserInp(inp)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  ! call printUserInp(inp)
  call system_clock(count_rate=rate)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Secondly, get couplings
  ! In the very first run, the following subroutine will calculate the
  ! NA-couplings from WAVECARs and then write it to a binary file called
  ! COUPCAR.  From the second run on, the subroutine will just read the
  ! NA-couplings from the file. However, for a general NAMD run, the file is way
  ! too huge, the solution is to write only the information we need to another
  ! plain text file. If such files exist (set LCPTXT = .TRUE. in the inp), then
  ! we may skip the huge binary file and read the plain text file instead. This
  ! is done in the subroutine 'initTDKS'.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (inp%LEPC) then
    call TDepCoupIJ(olap, olap_sec, inp, epc)
  else
    call TDCoupIJ(trim(inp%rundir), olap, olap_sec, inp)
  end if
  ! write(*,*) "T_coup: ", fin - start
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! if (irank==0) then
  do ns=1, inp%NSAMPLE
    inp%NAMDTINI = inp%NAMDTINI_A(ns)
    inp%INIBAND  = inp%INIBAND_A(ns,:)
    inp%INIKPT  = inp%INIKPT_A(ns,:)
    if (irank==0) call printUserInp(inp)
    call system_clock(t1)
    if (inp%LCPROP .OR. ns==1) then
      ! initiate KS matrix
      call initTDKS(ks, inp, olap_sec, epc)
    ! Time propagation
    ! call Propagation(ks, inp, olap_sec, epc)
    end if
    ! Run surface hopping
    if (inp%LSHP) then
    ! call runSH(ks, inp, olap_sec, epc)
      call runSH_EPC(ks, inp, olap_sec, epc)
    ! call runSHotf(ks, inp, olap_sec, epc)
      ! if (irank==0) call printSH(ks, inp)
      ! if (inp%LEPC) call phQ2R(inp, olap_sec, epc, ks%ph_pops)
      ! if (inp%LEPC) call saveXDAT(inp, olap_sec, epc)
      ! if (inp%LEPC) call printPHPROP(ks, inp, olap_sec, ns)
    end if
    call system_clock(t2)
    if (irank==0) write(*,'(A, F10.2)') "CPU Time [s]:", (t2-t1)/real(rate)
  end do
  ! end if

  call MPI_FINALIZE(ierr)

end Program
