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

  call getUserInp(inp)
  !call TDepCoupIJ(olap, olap_sec, inp, epc)
  !call readEPC(inp, epc)
  !call phDecomp(inp, epc)
  !call kqMatch(epc)
  call readEPCf(inp, epc, olap)


end Program

