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

  call readEPC(inp,epc)
  call readDISPL(inp, epc)
  call cellPROJ(epc)


end Program

