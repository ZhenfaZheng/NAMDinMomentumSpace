Program main
  use prec
  use fileio
  use epcoup

  implicit none

  type(namdInfo) :: inp
  type(epCoupling) :: epc

  call readEPC(inp, epc)
  call readPhmodes(inp, epc)
  call readDISPL(inp, epc)
  call cellPROJ(epc)

end Program
