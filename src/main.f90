Program main
  use prec
  use lattice
  use wavecar
  use couplings
  use epcoup
  use hamil
  use shop
  use fileio
  use TimeProp

  implicit none

  type(namdInfo) :: inp
  type(TDKS) :: ks
  type(overlap) :: olap, olap_sec
  type(epCoupling) :: epc

  real(kind=q) :: start, fin
  integer :: t1, t2, rate
  integer :: ns

  write(*,*)
  write(*,*) "Hefei-NAMD (epc version 1.9.1, May 04, 2022)"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First, get user inputs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getUserInp(inp)
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ns=1, inp%NSAMPLE
    inp%NAMDTINI = inp%NAMDTINI_A(ns)
    inp%INIBAND  = inp%INIBAND_A(ns,:)
    inp%INIKPT  = inp%INIKPT_A(ns,:)
    call printUserInp(inp)
    ! initiate KS matrix
    call system_clock(t1)
    call initTDKS(ks, inp, olap_sec)
    ! Time propagation
    call Propagation(ks, inp, olap_sec)
    ! Run surface hopping
    if (inp%LSHP) then
      call runSH(ks, inp, olap_sec)
      call printSH(ks, inp)
      if (inp%LEPC) call printPHPROP(ks, inp, olap_sec, ns)
    end if
    call system_clock(t2)
    write(*,'(A, F10.2)') "CPU Time [s]:", (t2-t1)/real(rate)
  end do

end Program
