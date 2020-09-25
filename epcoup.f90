module epcoup
  use prec
  use constants
  use fileio
  use couplings

  implicit none

  type epCoupling
    integer :: nbands, nkpts, nmodes, nqpts, natepc, natmd
    integer, allocatable, dimension(:) :: atnum
    !! atom number of each md cell atom in phonon cell
    integer, allocatable, dimension(:,:) :: R
    !! cell numbers for each atoms of md cell
    integer, allocatable, dimension(:) :: qqmap
    !! map q points in phonon calculation to q points in e-p calcultion.
    integer, allocatable, dimension(:,:,:) :: kkqmap
    !! map ik & jk of electronic states to q or -q of e-p matrix.
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
    integer :: ib, jb, ik, im, iq

    open(unit=90, file=inp%FILEPC, status='unknown', action='read', iostat=ierr)
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
    do ib=1,nb
      do jb=1,nb
        do ik=1,nk
          do im=1,nm
            do iq=1,nq
              read(unit=90, fmt='(2f15.10)') epc%epmat(ib,jb,ik,im,iq)
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
    integer :: ib, jb, ik, im, iq
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
    do ib=1,nbands
      do jb=1,nbands
        do ik=1,nkpts
          do im=1,nmodes
            do iq=1,nqs
              read(unit=908, fmt=*) epc%epmat(ib, jb, ik, im, iq)
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
    integer :: iq, im, ia, iax
    integer :: nqs, nmodes, nat
    character(len=24) :: charac, bra, ket
    complex(kind=DP) :: temp
    real(kind=DP), allocatable, dimension(:) :: qpt
    real(kind=DP), allocatable, dimension(:,:) :: at

    open(unit=909, file=inp%FILPH, status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "phonon modes file does NOT exist!"
      stop
    end if

    nat = epc%natepc
    nmodes = epc%nmodes
    nqs = epc%nqpts
    allocate(epc%phmodes(nqs, nmodes, nat, 3))
    allocate(epc%qptsph(nqs, 3))
    allocate(epc%freq(nqs, nmodes))
    allocate(qpt(3), at(3, 3))

    at = epc%cellepc(1:3,1:3) / epc%cellepc(1,1)
    ! write(*,'(3f12.7)') ((at(i,iax),i=1,3), iax=1,3)

    do iq=1,nqs
      read(unit=909, fmt=*)
      read(unit=909, fmt=*)
      read(unit=909, fmt=9019) charac, (qpt(iax), iax=1,3)
      do iax=1,3
        epc%qptsph(iq,iax) = at(iax,1)*qpt(1) + at(iax,2)*qpt(2) &
                             + at(iax,3)*qpt(3)
      end do
      ! write(*, 9019) charac, epc%qptsph(iq,:)
      read(unit=909, fmt=*)
      do im=1,nmodes
        read(unit=909, fmt=9011) charac, epc%freq(iq, im)
        ! write(*,9011) charac, epc%freq(iq, im)
        do ia=1,nat
          read(unit=909, fmt=9021) bra, (epc%phmodes(iq, im, ia, iax), &
                                         iax=1,3), ket
          ! write(*, 9021) bra, (epc%phmodes(iq, im, ia, iax), &
          !                      iax=1,3), ket
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

    real(kind=q) :: scal, dr
    real(kind=q), allocatable, dimension(:,:) :: supercell
    integer :: ierr, i, iax, ia, t
    integer :: nat, mdtime

    mdtime = inp%NSW

    open(unit=33, file=inp%FILMD, status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "XDATCAR file does NOT exist!"
      stop
    end if

    allocate(supercell(3,3))

    read(unit=33, fmt=*)
    read(unit=33, fmt=*) scal
    do i=1,3
      read(unit=33, fmt=*) (supercell(i, iax), iax=1,3)
      supercell(i,:) = supercell(i,:) * scal
    end do
    read(unit=33, fmt=*)
    read(unit=33, fmt=*) nat
    epc%natmd = nat

    allocate(epc%displ(mdtime, nat, 3), &
             epc%vel(mdtime, nat, 3))
    do t=1,mdtime
      read(unit=33, fmt=*)
      do ia=1,nat
        read(unit=33, fmt=*) (epc%displ(t, ia, iax), iax=1, 3)
        if (t==1) then
          epc%vel(t,ia,:) = 0
        else
          ! Fix atom position if atom moves to other cell.
          do iax=1,3
            dr = epc%displ(t, ia, iax) - epc%displ(t-1, ia, iax)
            if (dr>0.9) then
              epc%displ(t, ia, iax) = epc%displ(t, ia, iax) - 1
            else if (dr<-0.9) then
              epc%displ(t, ia, iax) = epc%displ(t, ia, iax) + 1
            endif
          enddo
          epc%vel(t,ia,:) &
          = ( epc%displ(t, ia, :) - epc%displ(t-1, ia, :) ) / inp%POTIM
        endif
      end do
    end do

    allocate(epc%cellmd(nat+3, 3))
    epc%cellmd(:3,:) = supercell
    ! epc%cellmd(4:,:) = epc%displ(1,:,:)
    do ia=1,nat
      do iax=1,3
        epc%cellmd(ia+3, iax) = SUM(epc%displ(:, ia, iax)) / mdtime
      end do
    end do

    do t=1,mdtime
      do ia=1,nat
        epc%displ(t,ia,:) &
        = epc%displ(t,ia,:) - epc%cellmd(ia+3, :)
      end do
    end do

    close(33)

  end subroutine readDISPL


  subroutine cellPROJ(epc)
    ! Calculate cell numbers for each atom of md cell
    ! Project from phonon cell to md cell.
    implicit none

    integer :: iax, ia, ja, i
    integer, allocatable, dimension(:) :: N
    !! scale numbers of md cell compared with phonon cell
    real(kind=DP), allocatable, dimension(:) :: dr
    real(kind=DP), allocatable, dimension(:) :: temp1, temp2

    type(epCoupling), intent(inout) :: epc

    allocate(epc%R(epc%natmd,3), epc%atnum(epc%natmd), &
             N(3), dr(3), temp1(3), temp2(3))

    do iax=1,3
      N(iax) = NINT( SUM(epc%cellmd(iax,:)) / SUM(epc%cellepc(iax,:)) )
    enddo

    do ia=1,epc%natmd
      temp1 = 9999.9
      do ja=1,epc%natepc
        dr = epc%cellmd(ia+3,:) * N - epc%cellepc(ja+3,:)
        do iax=1,3
          temp2(iax) = ABS( dr(iax) - NINT(dr(iax)) )
        enddo
        if (SUM(temp2)<SUM(temp1)) then
          temp1 = temp2
          epc%atnum(ia) = ja
          epc%R(ia,:) = (/(MOD(NINT(dr(i)),N(i)), i=1,3)/)
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
    real(kind=DP) :: Eu, Ek
    real(kind=DP) :: theta
    real(kind=DP), allocatable, dimension(:) :: displ, vel
    integer :: iq, im, ia, iax, t
    integer :: nqs, nmodes, nat

    call readPhmodes(inp, epc)
    call readDISPL(inp, epc)
    call cellPROJ(epc)

    nmodes = epc%nmodes
    nqs = epc%nqpts
    allocate(epc%phproj(inp%NSW, nmodes, nqs), displ(3), vel(3))

    do t=1,inp%NSW
      do iq=1,epc%nqpts
        do im=1,epc%nmodes
          Q = (0.0, 0.0)
          dQ = (0.0, 0.0)
          do ia=1,epc%natmd
            displ = 0.0; vel = 0.0
            do iax=1,3
              displ = displ + epc%displ(t, ia, iax) * epc%cellmd(iax, :) 
              vel = vel + epc%vel(t, ia, iax) * epc%cellmd(iax, :) 
            end do
            theta = 2 * PI * DOT_PRODUCT( epc%qptsph(iq,:), epc%R(ia,:) )
            Q  =  Q + EXP(imgUnit*theta) * DOT_PRODUCT( &
                      CONJG(epc%phmodes(iq, im, epc%atnum(ia), :)), displ )
            dQ = dQ + EXP(imgUnit*theta) * DOT_PRODUCT( &
                      CONJG(epc%phmodes(iq, im, epc%atnum(ia), :)), vel )
          end do
          ! mass of C = 12.011
          Eu = 0.5 * 12.011 * epc%freq(iq, im)**2 * CONJG(Q) * Q / epc%natmd &
            * 1.66 * 6.2415 / 100000
          Ek = 0.5 * 12.011 * CONJG(dQ) * dQ / epc%natmd &
            * 1.66 * 6.2415 * 10
          epc%phproj(t, im, iq) = (Eu + Ek) / ( hbar * epc%freq(iq, im) / 1000 )
        end do
      end do
    end do

  end subroutine phDecomp


  subroutine kqMatch(epc)
    implicit none

    type(epCoupling), intent(inout) :: epc

    real(kind=DP) :: norm
    real(kind=DP) :: dkq1(3), dkq2(3), dq(3)
    integer :: ik, jk, iq, jq, iax

    norm = 0.005
    ! If k1-k2 < norm, recognize k1 and k2 as same k point.
    ! So, number of kx, ky, kz or qx, qy, qz must not supass 1/norm = 200

    allocate(epc%kkqmap(epc%nkpts,epc%nkpts,2), epc%qqmap(epc%nqpts))

    epc%kkqmap = -1
    epc%qqmap = -1

    do ik=1,epc%nkpts
      do jk=1,epc%nkpts
        do iq=1,epc%nqpts
          dkq1 = epc%kptsepc(ik,:) - epc%kptsepc(jk,:) - epc%qptsepc(iq,:)
          dkq2 = epc%kptsepc(ik,:) - epc%kptsepc(jk,:) + epc%qptsepc(iq,:)
          do iax=1,3
            dkq1(iax) = ABS(dkq1(iax)-NINT(dkq1(iax)))
            dkq2(iax) = ABS(dkq2(iax)-NINT(dkq2(iax)))
          end do
          if (SUM(dkq1)<norm) epc%kkqmap(ik,jk,1) = iq
          if (SUM(dkq2)<norm) epc%kkqmap(ik,jk,2) = iq
          if (epc%kkqmap(ik,jk,1)>0 .and. epc%kkqmap(ik,jk,2)>0) exit
        end do
      end do
    end do

    do iq=1,epc%nqpts
      do jq=1,epc%nqpts
        dq = epc%qptsepc(iq,:) - epc%qptsph(jq,:)
        do iax=1,3
          dq(iax) = ABS(dq(iax)-NINT(dq(iax)))
        end do
        if (SUM(dq)<norm) then
          epc%qqmap(iq) = jq
          exit
        end if
      end do
    end do

  end subroutine kqMatch


  subroutine TDepCoupIJ(olap, olap_sec, inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap
    type(overlap), intent(inout) :: olap_sec

    real(kind=DP) :: proj, norm, phn
    real(kind=DP), allocatable, dimension(:) :: dkq, dq
    integer :: iq, jq1, jq2, im
    integer :: nqs, nmodes, nat
    integer :: t, ik, jk, ib, jb
    complex(kind=q) :: temp
    logical :: lcoup

    ! Initialization
    olap%NBANDS = inp%NBANDS * inp%NKPOINTS
    olap%TSTEPS = inp%NSW
    olap%dt = inp%POTIM
    allocate(olap%Dij(olap%NBANDS, olap%NBANDS, olap%TSTEPS-1))
    allocate(olap%Eig(olap%NBANDS, olap%TSTEPS-1))
    olap%Dij = cero
    olap%EIg = cero

    olap_sec%NBANDS = inp%NBASIS * inp%NKPOINTS
    olap_sec%TSTEPS = inp%NSW
    olap_sec%dt = inp%POTIM
    allocate(olap_sec%Dij(olap_sec%NBANDS, olap_sec%NBANDS, olap_sec%TSTEPS-1))
    allocate(olap_sec%Eig(olap_sec%NBANDS, olap_sec%TSTEPS-1))

    inquire(file='COUPCAR', exist=lcoup)
    if (lcoup) then
      ! file containing couplings exists, then read it
      if (inp%LCPTXT) then
        call readNaEig(olap_sec, inp)
      else
        call CoupFromFile(olap)
        call copyToSec(olap, olap_sec, inp)
        call writeNaEig(olap_sec, inp)
      end if
    else

      write(*,'(A)') "------------------------------------------------------------"
      write(*,*) "Calculating couplings from e-p matrix..."

      call readEPC(inp, epc)
      call phDecomp(inp, epc)
      call kqMatch(epc)
     
      allocate(dkq(3), dq(3))
     
      do t=1,inp%NSW-1
        do ib=1,inp%NBANDS
          do jb=1,inp%NBANDS
            do ik=1,inp%NKPOINTS
              do jk=1,inp%NKPOINTS

                iq = epc%kkqmap(ik,jk,1)

                if (iq>0) then

                  temp = (0.0, 0.0)
                  jq1 = epc%qqmap(epc%kkqmap(ik,jk,1)) ! ik-jk =  q
                  jq2 = epc%qqmap(epc%kkqmap(ik,jk,2)) ! ik-jk = -q

                  if (jq1>0) then
                    do im=1,epc%nmodes
                      phn = ABS(epc%phproj(t, im, jq1))
                      temp = temp + SQRT(phn) * epc%epmat(ib, jb, jk, im, iq)
                    end do
                  end if
                 
                  if (jq2>0) then
                    do im=1,epc%nmodes
                      phn = ABS(epc%phproj(t, im, jq2))
                      temp = temp + SQRT(phn+1) * epc%epmat(ib, jb, jk, im, iq)
                    end do
                  end if
                 
                  temp = temp / SQRT(1.0 * epc%natmd / epc%natepc) &
                         * 2 * inp%POTIM * imgUnit / hbar
                  olap%Dij(inp%NBANDS*(ik-1)+ib, jb+inp%NBANDS*(jk-1), t) = temp

                end if

              end do
            end do
          end do
        end do
      end do
     
      do ib=1,inp%NBANDS
        do ik=1,inp%NKPOINTS
          olap%Eig(ib+inp%NBANDS*(ik-1),:) = epc%energy(ik,ib)
        enddo
      enddo

      write(*,*) "Done..."
      write(*,'(A)') "------------------------------------------------------------"
     
      call copyToSec(olap, olap_sec, inp)
      call CoupToFile(olap)
      call writeNaEig(olap_sec, inp)

    end if

  end subroutine TDepCoupIJ


  subroutine copyToSec(olap, olap_sec, inp)
    implicit none
    type(overlap), intent(inout) :: olap_sec
    type(overlap), intent(in) :: olap
    type(namdInfo), intent(in) :: inp

    integer :: ik, jk, NB, NBas

    NB = inp%NBANDS
    NBas = inp%NBASIS
    do ik=1,inp%NKPOINTS
      olap_sec%Eig((ik-1)*NBas+1:ik*NBas, :) = &
      olap%Eig((ik-1)*NB+inp%BMIN:(ik-1)*NB+inp%BMAX, :)
      do jk=1,inp%NKPOINTS
        olap_sec%Dij((ik-1)*NBas+1:ik*NBas,(jk-1)*NBas+1:jk*NBas, :) = &
        olap%Dij((ik-1)*NB+inp%BMIN:(ik-1)*NB+inp%BMAX, &
                 (jk-1)*NB+inp%BMIN:(jk-1)*NB+inp%BMAX, :)
      end do
    end do

  end subroutine copyToSec


end module epcoup
