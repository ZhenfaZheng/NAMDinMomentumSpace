module epcoup
  use prec
  use constants
  use fileio
  use couplings
  use hdf5

  implicit none

  type epCoupling
    integer :: nbands, nkpts, nmodes, nqpts, natepc, natmd
    integer, allocatable, dimension(:) :: atmap
    !! mapping of each atom in MD cell to atom in unit cell for epc calculation.
    integer, allocatable, dimension(:,:) :: Rp
    !! number of lattice vector of each atom in MD cell, Rp(natmd,3)
    integer, allocatable, dimension(:) :: qqmap
    !! map q points in phonon calculation to q points in e-p calcultion.
    integer, allocatable, dimension(:,:) :: kkqmap
    !! map ik & jk of electronic states to q or -q of e-p matrix.
    real(kind=q), allocatable, dimension(:) :: mass
    real(kind=q), allocatable, dimension(:,:) :: energy
    real(kind=q), allocatable, dimension(:,:) :: cellep
    real(kind=q), allocatable, dimension(:,:) :: kptsep
    real(kind=q), allocatable, dimension(:,:) :: qptsep
    real(kind=q), allocatable, dimension(:,:) :: freqep
    real(kind=q), allocatable, dimension(:,:) :: freqph
    complex(kind=q), allocatable, dimension(:,:,:) :: normcoord
    complex(kind=q), allocatable, dimension(:,:,:,:,:) :: epmat
    complex(kind=q), allocatable, dimension(:,:,:,:) :: phmodes
    real(kind=q), allocatable, dimension(:,:) :: cellmd
    real(kind=q), allocatable, dimension(:,:,:) :: displ
  end type

  contains

  subroutine releaseEPC(epc)
    implicit none

    type(epCoupling), intent(inout) :: epc

    if ( allocated(epc%Rp) ) deallocate(epc%Rp)
    if ( allocated(epc%atmap) ) deallocate(epc%atmap)
    if ( allocated(epc%qqmap) ) deallocate(epc%qqmap)
    if ( allocated(epc%kkqmap) ) deallocate(epc%kkqmap)
    if ( allocated(epc%energy) ) deallocate(epc%energy)
    if ( allocated(epc%cellep) ) deallocate(epc%cellep)
    if ( allocated(epc%kptsep) ) deallocate(epc%kptsep)
    if ( allocated(epc%qptsep) ) deallocate(epc%qptsep)
    if ( allocated(epc%freqep) ) deallocate(epc%freqep)
    if ( allocated(epc%freqph) ) deallocate(epc%freqph)
    if ( allocated(epc%epmat) ) deallocate(epc%epmat)
    if ( allocated(epc%phmodes) ) deallocate(epc%phmodes)
    if ( allocated(epc%cellmd) ) deallocate(epc%cellmd)
    if ( allocated(epc%displ) ) deallocate(epc%displ)

  end subroutine releaseEPC


  subroutine readEPCpert(inp, epc, olap)
    ! Read information of e-ph coupling from perturbo.x output file
    ! (prefix_ephmat_p1.h5), which is in form of hdf5.
    ! Extract informations include: el bands, ph dispersion, k & q lists, ephmat

    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap

    integer :: hdferror, info(4)
    integer :: nk, nq, nb, nm, nat
    integer :: ib, jb, ik, jk, im, iq, it, ibas, jbas
    integer(hid_t) :: file_id, gr_id, dset_id
    integer(hsize_t) :: dim1(1), dim2(2), dim4(4)
    character(len=72) :: tagk, fname, grname, dsetname
    real(kind=q) :: kbT, phn, dE
    real(kind=q), allocatable, dimension(:,:) :: kqltemp, pos, lattvec
    real(kind=q), allocatable, dimension(:,:,:,:) :: eptemp_r, eptemp_i, phmtemp
    complex(kind=q), allocatable, dimension(:,:,:,:) :: eptemp
    complex :: eptemp_s, iomega

    fname = inp%FILEPM
    call h5open_f(hdferror)
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, hdferror)

    ! Read el bands, ph dispersion, k & q lists.
    grname = 'el_ph_band_info'
    call h5gopen_f(file_id, grname, gr_id, hdferror)

    dsetname = 'information'
    dim1 = shape(info, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, info, dim1, hdferror)
    call h5dclose_f(dset_id, hdferror)
    nk = info(1); epc%nkpts  = nk
    nq = info(2); epc%nqpts  = nq
    nb = info(3); epc%nbands = nb
    nm = info(4); epc%nmodes = nm
    nat = nm/3  ; epc%natepc = nat
    ! write(*,*) info

    dsetname = 'mass_a.u.'
    allocate(epc%mass(nat))
    dim1 = shape(epc%mass, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, epc%mass, dim1, hdferror)
    call h5dclose_f(dset_id, hdferror)

    dsetname = 'lattice_vec_angstrom'
    allocate(epc%cellep(nat+3,3), lattvec(3,3))
    dim2 = shape(lattvec, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, lattvec, dim2, hdferror)
    call h5dclose_f(dset_id, hdferror)
    epc%cellep(:3,:) = transpose(lattvec)

    dsetname = 'atom_pos'
    allocate(pos(3,nat))
    dim2 = shape(pos, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pos, dim2, hdferror)
    call h5dclose_f(dset_id, hdferror)
    epc%cellep(4:,:) = transpose(pos)

    dsetname = 'k_list'
    allocate(kqltemp(3,nk), epc%kptsep(nk,3))
    dim2 = shape(kqltemp, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, kqltemp, dim2, hdferror)
    call h5dclose_f(dset_id, hdferror)
    epc%kptsep = transpose(kqltemp)
    deallocate(kqltemp)

    dsetname = 'q_list'
    allocate(kqltemp(3,nq), epc%qptsep(nq,3))
    dim2 = shape(kqltemp, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, kqltemp, dim2, hdferror)
    call h5dclose_f(dset_id, hdferror)
    epc%qptsep = transpose(kqltemp)
    deallocate(kqltemp)

    dsetname = 'el_band_eV'
    allocate(epc%energy(nb,nk))
    dim2 = shape(epc%energy, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, epc%energy, dim2, hdferror)
    call h5dclose_f(dset_id, hdferror)

    do ik=1,nk
      do ib=1,nb
        olap%Eig(ib+nb*(ik-1),:) = epc%energy(ib,ik)
      end do
    end do

    dsetname = 'ph_disp_meV'
    allocate(epc%freqep(nm,nq))
    dim2 = shape(epc%freqep, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, epc%freqep, dim2, hdferror)
    call h5dclose_f(dset_id, hdferror)
    epc%freqep = epc%freqep / 1000.0_q ! transform unit to eV

    dsetname = 'phmod_ev_r'
    allocate(epc%phmodes(nq, nm, nat, 3))
    allocate(phmtemp(nq, nm, nat, 3))
    dim4 = shape(phmtemp, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, phmtemp, dim2, hdferror)
    call h5dclose_f(dset_id, hdferror)
    epc%phmodes = phmtemp

    dsetname = 'phmod_ev_i'
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, phmtemp, dim2, hdferror)
    call h5dclose_f(dset_id, hdferror)
    epc%phmodes = epc%phmodes + imgUnit * phmtemp

    call h5gclose_f(gr_id, hdferror)

    ! Read e-ph matrix for every k.
    grname = 'g_ephmat_total_meV'
    call h5gopen_f(file_id, grname, gr_id, hdferror)
    allocate(eptemp(nb, nb, nm, nq))
    allocate(eptemp_r(nb, nb, nm, nq))
    allocate(eptemp_i(nb, nb, nm, nq))
    allocate(epc%epmat(nb, nb, nk, nm, nq))

    call kqMatch(epc)
    kbT = inp%TEMP * BOLKEV

    do ik=1,nk

      write(tagk, '(I8)') ik

      dsetname = 'g_ik_r_' // trim( adjustl(tagk) )
      dim4 = shape(eptemp, kind=hsize_t)
      call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eptemp_r, dim4, hdferror)
      call h5dclose_f(dset_id, hdferror)
      ! if (ik==1) write(*,*) eptemp(1,1,:,2)

      dsetname = 'g_ik_i_' // trim( adjustl(tagk) )
      dim4 = shape(eptemp, kind=hsize_t)
      call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eptemp_i, dim4, hdferror)
      call h5dclose_f(dset_id, hdferror)

      eptemp = ( eptemp_r + imgUnit * eptemp_i ) / 1000.0_q
      epc%epmat(:,:,ik,:,:) = eptemp

      do jk=1,nk

        iq = epc%kkqmap(ik, jk)
        if (iq<0) cycle

        do im=1,nm

          if (epc%freqep(im,iq)<5.0E-3_q) cycle
          phn = 1.0 / (exp(epc%freqep(im,iq)/kbT)-1)
          iomega = imgUnit * epc%freqep(im,iq)/hbar

          do jb=1,nb
            do ib=1,nb

              ibas = nb*(ik-1)+ib
              jbas = nb*(jk-1)+jb
              eptemp_s = eptemp(ib, jb, im, iq)

              do it=1,inp%NSW-1
                olap%Dij(ibas, jbas, it) &
                = olap%Dij(ibas, jbas, it) + eptemp_s * &
                ( SQRT(phn)*exp(-iomega*it) + SQRT(phn+1)*exp(iomega*it) )
              end do ! it loop

            end do ! ib loop
          end do ! jb loop

        enddo ! im loop

      enddo ! jk loop

    enddo ! ik loop

    call h5gclose_f(gr_id, hdferror)

    call h5fclose_f(file_id, hdferror)

    olap%Dij = olap%Dij / SQRT(1.0*nq)
    do ib=1,nb*nk
      ! olap%Eig(ib,:) = olap%Eig(ib,:) + real(olap%Dij(ib,ib,:))
      ! olap%Dij(ib,ib,:) = cero
      olap%Dij(ib,ib,:) = real(olap%Dij(ib,ib,:))
      do jb=ib+1,nb*nk
        olap%Dij(jb,ib,:) = CONJG(olap%Dij(ib,jb,:))
      enddo
    enddo

  end subroutine readEPCpert


  subroutine readDISPL(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    real(kind=q) :: scal, dr, lattvec(3,3)
    integer :: ierr, i, iax, ia, t
    integer :: nat, mdtime

    mdtime = inp%NSW

    open(unit=33, file=inp%FILMD, status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "XDATCAR file does NOT exist!"
      stop
    end if

    read(unit=33, fmt=*)
    read(unit=33, fmt=*) scal
    do i=1,3
      read(unit=33, fmt=*) (lattvec(i, iax), iax=1,3)
      lattvec(i,:) = lattvec(i,:) * scal
    end do
    read(unit=33, fmt=*)
    read(unit=33, fmt=*) nat
    epc%natmd = nat

    allocate(epc%displ(mdtime, nat, 3))

    do t=1,mdtime
      read(unit=33, fmt=*)
      do ia=1,nat
        read(unit=33, fmt=*) (epc%displ(t, ia, iax), iax=1, 3)
        do iax=1,3
          if (t==1) cycle
          ! Modify atom position if atom moves to other cell.
          dr = epc%displ(t, ia, iax) - epc%displ(t-1, ia, iax)
          if (dr>0.9) then
            epc%displ(t, ia, iax) = epc%displ(t, ia, iax) - 1
          else if (dr<-0.9) then
            epc%displ(t, ia, iax) = epc%displ(t, ia, iax) + 1
          endif
        enddo
      end do
    end do

    allocate(epc%cellmd(nat+3, 3))
    epc%cellmd(:3,:) = lattvec
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

    integer :: iax, ia, ja, i, N(3)
    !! scale numbers of md cell compared with phonon cell
    real(kind=q) :: dr(3), temp1(3), temp2(3)

    type(epCoupling), intent(inout) :: epc

    allocate(epc%Rp(epc%natmd,3), epc%atmap(epc%natmd))

    do iax=1,3
      ! ONLY suit for sample shape cells !
      N(iax) = NINT( SUM(epc%cellmd(iax,:)) / SUM(epc%cellep(iax,:)) )
    enddo

    do ia=1,epc%natmd
      temp1 = 9999.9
      do ja=1,epc%natepc
        dr = epc%cellmd(ia+3,:) * N - epc%cellep(ja+3,:)
        do iax=1,3
          temp2(iax) = ABS( dr(iax) - NINT(dr(iax)) )
        enddo
        if (SUM(temp2)<SUM(temp1)) then
          temp1 = temp2
          epc%atmap(ia) = ja
          epc%Rp(ia,:) = (/(MOD(NINT(dr(i)),N(i)), i=1,3)/)
        endif
      enddo
    enddo

  end subroutine cellPROJ


  subroutine phDecomp(inp, epc)
    ! Project the atomic motion from an MD simulation to phonon modes.
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    ! Temporary normal mode coordinate
    complex(kind=q) :: tempQ
    real(kind=q) :: theta, displ(3), Np
    integer :: iq, im, ia, ja, iax, t
    integer :: nqs, nmodes, nat

    call readDISPL(inp, epc)
    call cellPROJ(epc)

    nqs = epc%nqpts
    nmodes = epc%nmodes

    allocate(epc%normcoord(inp%NSW, nmodes, nqs))

    do t=1,inp%NSW
      do iq=1,epc%nqpts
        do im=1,epc%nmodes

          tempQ = cero
          do ia=1,epc%natmd
            displ = 0.0
            ja = epc%atmap(ia)
            do iax=1,3
              displ = displ + epc%displ(t, ia, iax) * epc%cellmd(iax, :) 
            end do
            theta = 2 * PI * DOT_PRODUCT( epc%qptsep(iq,:), epc%Rp(ia,:) )
            tempQ  = tempQ + EXP(-imgUnit*theta) * SQRT(epc%mass(ja)) * &
                     DOT_PRODUCT( CONJG(epc%phmodes(iq, im, ja, :)), displ )
          end do
          epc%normcoord(t, im, iq) = tempQ

        end do
      end do
    end do

    Np = 1.0_q * epc%natmd / epc%natepc
    ! Np is the number of unit cells in MD cell.
    epc%normcoord = epc%normcoord / SQRT(Np)
    close(unit=41)

  end subroutine phDecomp


  subroutine kqMatch(epc)
    implicit none

    type(epCoupling), intent(inout) :: epc

    real(kind=q) :: norm
    real(kind=q) :: dkq(3), dq(3)
    integer :: ik, jk, iq, jq, iax

    norm = 0.001
    ! If k1-k2 < norm, recognize k1 and k2 as same k point.
    ! So, number of kx, ky, kz or qx, qy, qz must not supass 1/norm = 1000

    allocate(epc%kkqmap(epc%nkpts,epc%nkpts), epc%qqmap(epc%nqpts))

    epc%kkqmap = -1
    epc%qqmap = -1

    do ik=1,epc%nkpts
      do jk=1,epc%nkpts
        do iq=1,epc%nqpts
          dkq = epc%kptsep(ik,:) - epc%kptsep(jk,:) - epc%qptsep(iq,:)
          do iax=1,3
            dkq(iax) = ABS(dkq(iax)-NINT(dkq(iax)))
          end do
          if (SUM(dkq)<norm) then
              epc%kkqmap(ik,jk) = iq
              exit
          end if
        end do
      end do
    end do

  end subroutine kqMatch


  subroutine calcEPC(olap, inp, epc)
    implicit none

    type(overlap), intent(inout) :: olap
    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(in) :: epc

    integer :: nsw, nk, nb, nq, nm, Np
    integer :: it, ik, jk, ib, jb, iq, im, ibas, jbas
    real(kind=q) :: lqvtemp, lqv ! zero-point displacement amplitude
    real(kind=q) :: kbT, phn

    nsw = inp%NSW
    nk  = inp%NKPOINTS
    nb  = inp%NBANDS
    nm  = epc%nmodes

    lqvtemp = hbar * SQRT( EVTOJ / (2.0_q*AMTOKG) ) * 1.0E-5_q

    do ik=1,nk
      do jk=1,nk

        iq = epc%kkqmap(ik,jk)
        if (iq<0) cycle

        do ib=1,nb
          do jb=1,nb

            ibas = nb*(ik-1)+ib
            jbas = nb*(jk-1)+jb

            do im=1,nm
              lqv = lqvtemp / SQRT(epc%freqep(im, iq)) ! in unit of Angstrom
              do it=1,nsw-1
                olap%Dij(ibas, jbas, it) = olap%Dij(ibas, jbas, it) &
                + epc%epmat(ib,jb,ik,im,iq) * epc%normcoord(it,im,iq) / lqv
              end do
            end do

          end do
        end do

      end do
    end do

    do ib=1,nb
      do ik=1,nk
        ibas = nb*(ik-1)+ib
        olap%Eig(ibas,:) = epc%energy(ib,ik)
        ! olap%Dij(ibas,ibas,:) = ABS(olap%Dij(ib,ib,:))
        ! if (olap%Eig(ibas,1)<-2.8 .and. olap%Eig(ibas,1)>-3.2) &
        ! olap%Eig(ibas,:) = epc%energy(ib,ik) + (ik-25) * 0.02
        ! olap%Eig(ibas,:) = -4.4
        ! olap%Eig(ibas,:) = olap%Eig(ibas,:) + ABS(olap%Dij(ibas,ibas,200)) * 4
        ! olap%Dij(ibas,ibas,:) = cero
        do jbas=ibas+1,nb*nk
          olap%Dij(jbas,ibas,:) = CONJG(olap%Dij(ibas,jbas,:))
        enddo
      enddo
    enddo

    Np = epc%natmd / epc%natepc
    olap%Dij = olap%Dij / SQRT(1.0 * Np)
    
  end subroutine calcEPC


  subroutine TDepCoupIJ(olap, olap_sec, inp, epc)
    implicit none

    type(namdInfo), intent(inout) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap
    type(overlap), intent(inout) :: olap_sec

    real(kind=q) :: proj, norm, phn
    real(kind=q), allocatable, dimension(:) :: dkq, dq
    integer :: iq, jq1, jq2, im
    integer :: nqs, nmodes, nat
    integer :: t, ik, jk, ib, jb
    complex(kind=q) :: temp
    logical :: lcoup

    ! Initialization
    if ( (.not. inp%LCPTXT) .or. (.not. inp%LBASSEL) ) &
      call initOlap(olap, inp, inp%NBANDS * inp%NKPOINTS)

    inquire(file='COUPCAR', exist=lcoup)
    if (lcoup) then
      ! file containing couplings exists, then read it
      if (inp%LCPTXT .and. inp%LBASSEL) then
        call readBasis(inp)
        call initOlap(olap_sec, inp, inp%NBASIS)
        call readNaEig(olap_sec, inp)
      else
        call CoupFromFile(olap)
        call copyToSec(olap, olap_sec, inp)
        call writeNaEig(olap_sec, inp)
      end if

    else if (inp%EPCTYPE==1) then

      write(*,'(A)') "------------------------------------------------------------"
      write(*,*) "Calculating couplings from dense e-p matrix..."

      call readEPCpert(inp, epc, olap)
      call copyToSec(olap, olap_sec, inp)
      call CoupToFile(olap)
      call writeNaEig(olap_sec, inp)

      write(*,*) "Done!"
      write(*,'(A)') "------------------------------------------------------------"

    else

      write(*,'(A)') "------------------------------------------------------------"
      write(*,*) "Calculating couplings from e-p matrix..."

      call readEPCpert(inp, epc, olap)
      call phDecomp(inp, epc)
      ! call kqMatch(epc)
      call calcEPC(olap, inp, epc)

      write(*,*) "Done!"
      write(*,'(A)') "------------------------------------------------------------"

      call copyToSec(olap, olap_sec, inp)
      call CoupToFile(olap)
      call writeNaEig(olap_sec, inp)

    end if

    call releaseEPC(epc)
    if ( allocated(olap%Dij) )deallocate(olap%Dij, olap%Eig)

  end subroutine TDepCoupIJ


  subroutine readBasis(inp)
    implicit none

    type(namdInfo), intent(inout) :: inp
    integer :: ierr, ibas, ik, ib

    open(unit=38, file='BASSEL', status='unknown', action='read', iostat=ierr)

    if (ierr /= 0) then
      write(*,*) "XDATCAR file does NOT exist!"
      stop
    end if

    inp%BASSEL = 0
    read(unit=38, fmt=*) inp%NBASIS
    do ibas=1,inp%NBASIS
      read(unit=38, fmt=*) ik, ib
      inp%BASSEL(ik, ib) = ibas
    end do

    close(unit=38)

  end subroutine readBasis


  subroutine selBasis(inp, olap)
    implicit none

    type(overlap), intent(in) :: olap
    type(namdInfo), intent(inout) :: inp

    integer :: ik, ib, nb
    real :: emin, emax

    nb = inp%NBANDS
    emin = inp%EMIN
    emax = inp%EMAX

    inp%NBASIS = 0
    inp%BASSEL = 0

    do ik=inp%KMIN, inp%KMAX
      do ib=inp%BMIN, inp%BMAX
        if ( olap%Eig((ik-1)*nb+ib, 1)>emin .and. &
             olap%Eig((ik-1)*nb+ib, 1)<emax ) then
          inp%NBASIS = inp%NBASIS + 1
          inp%BASSEL(ik,ib) = inp%NBASIS
        end if
      end do
    end do

    open(unit=39, file='BASSEL', status='unknown', action='write')

    write(unit=39, fmt='(I18)') inp%NBASIS

    do ik=inp%KMIN, inp%KMAX
      do ib=inp%BMIN, inp%BMAX
        if (inp%BASSEL(ik,ib)>0) &
          write(unit=39, fmt='(2I12)') ik, ib
      end do
    end do

    close(unit=39)

    !! For array element inp%BASSEL(ik,ib),
    !! 0 represent ik,ib state isn't selected as basis;
    !! number >0 represent ik,ib state is selected as basis,
    !! and the number is the state's serial number among the basises.

  end subroutine selBasis


  subroutine initOlap(olap, inp, nb)
    implicit none
    type(overlap), intent(inout) :: olap
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: nb

    olap%NBANDS = nb
    olap%TSTEPS = inp%NSW
    olap%dt = inp%POTIM

    allocate(olap%Dij(nb, nb, inp%NSW-1))
    allocate(olap%Eig(nb, inp%NSW-1))

    olap%Dij = cero
    olap%Eig = 0.0_q

  end subroutine initOlap


  subroutine copyToSec(olap, olap_sec, inp)
    implicit none
    type(namdInfo), intent(inout) :: inp
    type(overlap), intent(inout) :: olap_sec
    type(overlap), intent(in) :: olap

    integer :: ik, jk, ib, jb, nb, iBas, jBas

    if (inp%LBASSEL) then
      call readBasis(inp)
    else
      call selBasis(inp, olap)
    end if

    nb = inp%NBANDS

    call initOlap(olap_sec, inp, inp%NBASIS)

    do ik=inp%KMIN, inp%KMAX
      do ib=inp%BMIN, inp%BMAX

        if (inp%BASSEL(ik,ib)>0) then

          iBas = inp%BASSEL(ik,ib)
          olap_sec%Eig(iBas, :) = olap%Eig((ik-1)*nb+ib, :)

          do jk=inp%KMIN, inp%KMAX
            do jb=inp%BMIN, inp%BMAX

              if (inp%BASSEL(jk,jb)>0) then

                jBas = inp%BASSEL(jk,jb)
                olap_sec%Dij(iBas, jBas, :) = &
                olap%Dij( (ik-1)*nb+ib, (jk-1)*nb+jb, : )

              end if

            end do
          end do

        end if

      end do
    end do

  end subroutine copyToSec


end module epcoup
