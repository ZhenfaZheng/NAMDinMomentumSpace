Program main
  use prec
  use fileio
  use epcoup

  implicit none

  type(namdInfo) :: inp
  type(epCoupling) :: epc

  call ReadPhmodes(inp, epc)

end Program
