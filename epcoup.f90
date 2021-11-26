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
    if ( allocated(epc%phmodes) ) deallocate(epc%phmodes)
    if ( allocated(epc%cellmd) ) deallocate(epc%cellmd)
    if ( allocated(epc%displ) ) deallocate(epc%displ)

  end subroutine


  subroutine readEPCinfo(inp, epc)
    implicit none

    type(namdInfo), intent(inout) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: hdferror, info(4)
    integer(hsize_t) :: dim1(1)
    integer :: nk, nq, nb, nm, nat
    integer(hid_t) :: file_id, gr_id, dset_id
    character(len=256) :: fname, grname, dsetname

    fname = inp%FILEPM
    call h5open_f(hdferror)
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, hdferror)
    grname = 'el_ph_band_info'
    call h5gopen_f(file_id, grname, gr_id, hdferror)
    dsetname = 'information'
    dim1 = shape(info, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)

    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, info, dim1, hdferror)

    call h5dclose_f(dset_id, hdferror)
    call h5gclose_f(gr_id, hdferror)
    call h5fclose_f(file_id, hdferror)

    nk = info(1); epc%nkpts  = nk
    nq = info(2); epc%nqpts  = nq
    nb = info(3); epc%nbands = nb
    nm = info(4); epc%nmodes = nm
    nat = nm/3  ; epc%natepc = nat

    if (inp%NBANDS .NE. nb) then
      write(*,*) 'NBANDS seems to be wrong!'
      stop
    end if
    if (inp%NKPOINTS .NE. nk) then
      write(*,*) 'NKPOINTS seems to be wrong!'
      stop
    end if
    if ( (inp%NMODES .NE. 1) .AND. (inp%NMODES .NE. nm) ) then
      write(*,*) 'NMODES seems to be wrong!'
      stop
    end if
    if ( inp%Np==1 ) inp%Np = nq
    inp%NMODES = nm

  end subroutine


  subroutine readEPCpert(inp, epc, olap)
    ! Read information of e-ph coupling from perturbo.x output file
    ! (prefix_ephmat_p1.h5), which is in form of hdf5.
    ! Extract informations include: el bands, ph dispersion, k & q lists, ephmat

    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap

    integer :: hdferror
    integer :: nk, nq, nb, nm, nat
    integer :: ib, jb, ik, jk, im, iq, it, ibas, jbas
    integer(hid_t) :: file_id, gr_id, dset_id
    integer(hsize_t) :: dim1(1), dim2(2), dim4(4)
    character(len=72) :: tagk, fname, grname, dsetname
    real(kind=q), allocatable, dimension(:,:) :: kqltemp, pos, lattvec
    real(kind=q), allocatable, dimension(:,:,:,:) :: eptemp_r, eptemp_i, phmtemp
    complex(kind=q), allocatable, dimension(:,:,:,:) :: eptemp

    nk  = epc%nkpts
    nq  = epc%nqpts
    nb  = epc%nbands
    nm  = epc%nmodes
    nat = epc%natepc

    write(*,*) "Reading ephmat.h5 file."

    fname = inp%FILEPM

    call h5open_f(hdferror)
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, hdferror)

    ! Read el bands, ph dispersion, k & q lists.
    grname = 'el_ph_band_info'
    call h5gopen_f(file_id, grname, gr_id, hdferror)

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


    if (inp%EPCTYPE==2) then
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
    end if

    call h5gclose_f(gr_id, hdferror)


    ! Read e-ph matrix for every k.
    grname = 'g_ephmat_total_meV'
    call h5gopen_f(file_id, grname, gr_id, hdferror)
    allocate(eptemp(nb, nb, nm, nq))
    allocate(eptemp_r(nb, nb, nm, nq))
    allocate(eptemp_i(nb, nb, nm, nq))

    call kqMatch(epc)

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

      do jk=1,nk

        iq = epc%kkqmap(ik, jk)
        if (iq<0) cycle

        ! print '(4(3(F10.5,3x)))', &
        ! epc%kptsep(ik,:), epc%kptsep(jk,:), epc%qptsep(iq,:), &
        ! epc%kptsep(ik,:)-epc%kptsep(jk,:)-epc%qptsep(iq,:)

        do im=1,nm

          if (epc%freqep(im,iq)<5.0E-3_q) cycle

          do jb=1,nb
            do ib=1,nb

              ibas = nb*(ik-1)+ib
              jbas = nb*(jk-1)+jb
              olap%gij(ibas, jbas, im) = eptemp(ib, jb, im, iq)
              olap%Phfreq(ibas, jbas, im) = epc%freqep(im,iq)

            end do ! ib loop
          end do ! jb loop

        enddo ! im loop

      enddo ! jk loop

    enddo ! ik loop

    call h5gclose_f(gr_id, hdferror)

    call h5fclose_f(file_id, hdferror)

  end subroutine


  subroutine readDISPL(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    real(kind=q) :: scal, dr, lattvec(3,3), temp(3)
    integer :: ierr, i, iax, ia, it
    integer :: nat, nsw

    nsw = inp%NSW
    write(*,*) 'Reading MD traj.'

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

    allocate(epc%displ(nsw, nat, 3))

    do it=1,nsw
      read(unit=33, fmt=*)
      do ia=1,nat
        read(unit=33, fmt=*) (epc%displ(it, ia, iax), iax=1, 3)
        do iax=1,3
          if (it==1) cycle
          ! Modify atom position if atom moves to other cell.
          dr = epc%displ(it, ia, iax) - epc%displ(it-1, ia, iax)
          if (dr>0.9) then
            epc%displ(it, ia, iax) = epc%displ(it, ia, iax) - 1
          else if (dr<-0.9) then
            epc%displ(it, ia, iax) = epc%displ(it, ia, iax) + 1
          endif
        enddo
      end do
    end do

    allocate(epc%cellmd(nat+3, 3))
    epc%cellmd(:3,:) = lattvec
    ! epc%cellmd(4:,:) = epc%displ(1,:,:)
    do ia=1,nat
      do iax=1,3
        epc%cellmd(ia+3, iax) = SUM(epc%displ(:, ia, iax)) / nsw
      end do
    end do

    do it=1,nsw
      do ia=1,nat
        temp = 0.0_q
        do iax=1,3
          ! minus average coordinates.
          ! change from crastal coordinates to cartesian coordinates.
          temp = temp + ( epc%displ(it,ia,iax) - epc%cellmd(ia+3, iax) ) * &
                 lattvec(iax, :)
        end do
        epc%displ(it,ia,:) = temp
      end do
    end do

    close(33)

  end subroutine


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

  end subroutine


  subroutine phDecomp(inp, olap, epc)
    ! Project the atomic motion from an MD simulation to phonon modes.
    implicit none

    type(namdInfo), intent(in) :: inp
    type(overlap), intent(inout) :: olap
    type(epCoupling), intent(inout) :: epc

    ! Temporary normal mode coordinate
    complex(kind=q) :: temp
    complex(kind=q), allocatable, dimension(:,:) :: eiqR
    complex(kind=q), allocatable, dimension(:,:,:) :: tempQ
    real(kind=q) :: theta, displ(3), Np
    integer :: iq, im, ia, ja, iax, it
    integer :: nqs, nmodes, nat, nsw
    integer :: ik, jk, nks, ib, jb, ibas, jbas, nb

    call readDISPL(inp, epc)
    call cellPROJ(epc)

    write(*,*) "Decomposing phonon modes from MD traj."

    nks = epc%nkpts
    nb = inp%NBANDS
    nqs = epc%nqpts
    nmodes = epc%nmodes
    nsw = inp%NSW
    nat = epc%natmd

    allocate(eiqR(nqs,nat), tempQ(nqs, nmodes, nsw-1))

    do iq=1,nqs
      do ia=1,nat
        theta = 2 * PI * DOT_PRODUCT( epc%qptsep(iq,:), epc%Rp(ia,:) )
        eiqR(iq, ia) = EXP(-imgUnit*theta)
      end do
    end do

    do iq=1,nqs
      do im=1,nmodes
        do it=1,nsw-1
          temp = cero
          do ia=1,nat
            ja = epc%atmap(ia)
            temp = temp + eiqR(iq,ia) * SQRT(epc%mass(ja)) * &
                   DOT_PRODUCT( CONJG(epc%phmodes(iq, im, ja, :)), &
                   epc%displ(it, ia, :) )
          end do
          tempQ(iq,im,it) = temp
        end do
      end do
    end do

    do ik=1,nks
      do jk=1,nks
        iq = epc%kkqmap(ik, jk)
        do ib=1,nb
          ibas = nb*(ik-1)+ib
          do jb=1,nb
            jbas = nb*(jk-1)+jb
            olap%PhQ(ibas, jbas, :, :) = tempQ(iq,:,:)
          end do
        end do
      end do
    end do

    deallocate(tempQ)

    ! Np is the number of unit cells in MD cell.
    Np = epc%natmd / epc%natepc
    olap%PhQ = olap%PhQ / SQRT(Np)

  end subroutine


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

  end subroutine


  subroutine calcEPC(olap, inp)
    implicit none

    type(overlap), intent(inout) :: olap
    type(namdInfo), intent(inout) :: inp

    integer :: nb, nm, nsw
    integer :: ib, jb, im, it
    ! zero-point displacement amplitude
    real(kind=q) :: lqvtemp, lqv
    real(kind=q) :: kbT, phn
    real(kind=q) :: dE, dE1, dE2
    complex(kind=q) :: iw, iwtemp, idwt, eptemp

    nsw = inp%NSW
    nb  = olap%NBANDS
    nm  = olap%NMODES
    kbT = inp%TEMP * BOLKEV
    idwt = imgUnit / kbT
    ! iwt = i * (dE/hbar) * (hbar/kbT)

    write(*,*) "Calculating e-ph couplings."

    allocate(olap%Dij(nb, nb, nsw-1))
    olap%Dij = cero

    if (olap%COUPTYPE==1) then

     iwtemp = imgUnit * inp%POTIM / hbar
     allocate(olap%EPcoup(nb, nb, nm, 2, nsw-1))
     olap%EPcoup = cero

     do ib=1,nb
      do jb=ib,nb
        dE = olap%Eig(ib,1) - olap%Eig(jb,1)
        do im=1,nm

          if (olap%Phfreq(ib, jb, im)<5.0E-3_q) cycle
          phn = 1.0 / ( exp(olap%Phfreq(ib, jb, im)/kbT) - 1.0 )
          eptemp = olap%gij(ib, jb, im) * SQRT(2.0 * phn + 1.0) * 0.5
          if (jb==ib)  eptemp = ABS(eptemp)
          iw = iwtemp * olap%Phfreq(ib, jb, im)

          do it=1,nsw-1
            olap%EPcoup(ib, jb, im, 1, it) =  eptemp * exp(-iw*it)
            olap%EPcoup(ib, jb, im, 2, it) =  eptemp * exp(iw*it)
          end do ! it loop

          olap%Dij(ib,jb,:) = olap%Dij(ib,jb,:) + &
              SUM(olap%EPcoup(ib,jb,im,:,:), dim=1)

          dE1 = dE - olap%Phfreq(ib, jb, im) - 1.0E-8_q
          olap%EPcoup(ib,jb,im,1,:) = olap%EPcoup(ib,jb,im,1,:) * &
            ( exp(idwt * dE1) - 1.0 ) / ( idwt * dE1 )
          dE2 = dE + olap%Phfreq(ib, jb, im) + 1.0E-8_q
          olap%EPcoup(ib,jb,im,2,:) = olap%EPcoup(ib,jb,im,2,:) * &
            ( exp(idwt * dE2) - 1.0 ) / ( idwt * dE2 )

        end do ! im loop

        olap%Dij(jb,ib,:) = CONJG(olap%Dij(ib,jb,:))
        olap%EPcoup(jb,ib,:,:,:) = CONJG(olap%EPcoup(ib,jb,:,:,:))

      end do
     end do
     olap%EPcoup = olap%EPcoup / SQRT(olap%Np)
     olap%Dij = olap%Dij / SQRT(olap%Np)

    else if (olap%COUPTYPE==2) then

     lqvtemp = hbar * SQRT( EVTOJ / (2.0_q*AMTOKG) ) * 1.0E-5_q
     allocate(olap%EPcoup(nb, nb, nm, 1, nsw-1))
     olap%EPcoup = cero

     do ib=1,nb
      do jb=1,nb

        dE = olap%Eig(ib,1) - olap%Eig(jb,1)
        do im=1,nm
          ! unit of Angstrom
          lqv = lqvtemp / SQRT(olap%Phfreq(ib, jb, im))
          olap%EPcoup(ib,jb,im,1,:) &
            = olap%gij(ib, jb, im) * olap%PhQ(ib, jb, im, :) / lqv

          if (jb==ib) &
            olap%EPcoup(ib,ib,im,1,:) = ABS(olap%EPcoup(ib,ib,im,1,:))

          olap%Dij(ib,jb,:) = olap%Dij(ib,jb,:) + olap%EPcoup(ib,jb,im,1,:)

          dE1 = dE - olap%Phfreq(ib, jb, im) - 1.0E-8_q
          dE2 = dE + olap%Phfreq(ib, jb, im) + 1.0E-8_q
          olap%EPcoup(ib,jb,im,1,:) = olap%EPcoup(ib,jb,im,1,:) * &
            ( (exp(idwt * dE1) - 1.0) / (idwt * dE1) + &
              (exp(idwt * dE2) - 1.0) / (idwt * dE2) ) / 2.0
        end do

        olap%Dij(jb,ib,:) = CONJG(olap%Dij(ib,jb,:))
        olap%EPcoup(jb,ib,:,:,:) = CONJG(olap%EPcoup(ib,jb,:,:,:))

      end do
     end do
     olap%EPcoup = olap%EPcoup / SQRT(olap%Np)
     olap%Dij = olap%Dij / SQRT(olap%Np)

    end if

  end subroutine


  subroutine TDepCoupIJ(olap, olap_sec, inp, epc)
    implicit none

    type(namdInfo), intent(inout) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap
    type(overlap), intent(inout) :: olap_sec

    integer :: nmodes, nb, nsw
    logical :: lcoup

    write(*,*)
    write(*,'(A)') "------------------------------------------------------------"

    call readEPCinfo(inp, epc)

    ! Initialization
    if ( .not. ( inp%LCPTXT .and. inp%LBASSEL) ) &
      call initOlap(olap, inp, inp%NBANDS * inp%NKPOINTS)

    inquire(file='EPCAR', exist=lcoup)
    if (lcoup) then
      ! file containing couplings exists, then read it
      if (inp%LCPTXT .and. inp%LBASSEL) then
        call readBasis(inp)
        nb = inp%NBASIS
        call initOlap(olap_sec, inp, nb)

        nsw = inp%NSW
        nmodes = inp%NMODES
        allocate(olap_sec%Dij(nb, nb, nsw-1))
        allocate(olap_sec%EPcoup(nb, nb, nmodes, 1, nsw-1))
        olap%Dij = cero ; olap%EPcoup = cero

        call readNaEig(olap_sec, inp)
        call readEP(olap_sec)
      else
        call CoupFromFileEP(olap)
        call copyToSec(olap, olap_sec, inp)
        call calcEPC(olap_sec, inp)
        call writeNaEig(olap_sec, inp)
        call writeEP(olap_sec)
      end if

    else if (inp%EPCTYPE==1) then

      write(*,*) "TypeI e-ph coupling calculation."

      call readEPCpert(inp, epc, olap)
      call CoupToFileEP(olap)
      call copyToSec(olap, olap_sec, inp)
      call calcEPC(olap_sec, inp)
      call writeNaEig(olap_sec, inp)
      call writeEP(olap_sec)

    else

      write(*,*) "TypeII e-ph coupling calculation."

      call readEPCpert(inp, epc, olap)
      call phDecomp(inp, olap, epc)
      call CoupToFileEP(olap)
      call copyToSec(olap, olap_sec, inp)
      call calcEPC(olap_sec, inp)
      call writeNaEig(olap_sec, inp)
      call writeEP(olap_sec)

    end if

    call releaseEPC(epc)
    if ( allocated(olap%Dij) )deallocate(olap%Dij, olap%Eig)

    write(*,*) "Done!"
    write(*,'(A)') "------------------------------------------------------------"
    write(*,*)

  end subroutine


  subroutine readBasis(inp)
    implicit none

    type(namdInfo), intent(inout) :: inp
    integer :: ierr, ibas, ik, ib

    open(unit=38, file='BASSEL', status='unknown', action='read', iostat=ierr)

    if (ierr /= 0) then
      write(*,*) "BASSEL file does NOT exist!"
      stop
    end if

    inp%BASSEL = -1
    read(unit=38, fmt=*) inp%NBASIS
    do ibas=1,inp%NBASIS
      read(unit=38, fmt=*) ik, ib
      inp%BASSEL(ik, ib) = ibas
    end do

    close(unit=38)

  end subroutine


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
    inp%BASSEL = -1

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

  end subroutine


  subroutine sortBasis(inp, olap)
    implicit none

    type(overlap), intent(in) :: olap
    type(namdInfo), intent(inout) :: inp

    integer :: ibas, i, j, ik, nk, ib, nb, nbas
    real(kind=q), dimension(:,:), allocatable :: basis
    real(kind=q) :: temp(3)

    nb = inp%NBANDS
    nk = inp%NKPOINTS
    nbas = inp%NBASIS
    allocate(basis(nbas, 3))

    ibas = 0
    do ik=1,nk
      do ib=1,nb
        if (inp%BASSEL(ik,ib)>0) then
          ibas = ibas + 1
          basis(ibas,1) = ik
          basis(ibas,2) = ib
          basis(ibas,3) = olap%Eig((ik-1)*nb+ib, 1)
        end if
      end do
    end do

    do i=1,nbas-1
      do j=1,nbas-i
        if ( basis(j+1,3) < basis(j,3) ) then
          temp = basis(j,:)
          basis(j,:) = basis(j+1,:)
          basis(j+1,:) = temp
        end if
      end do
    end do

    open(unit=39, file='BASSEL', status='unknown', action='write')
    write(unit=39, fmt='(I18)') inp%NBASIS

    do ibas=1,nbas
      ik = NINT(basis(ibas,1))
      ib = NINT(basis(ibas,2))
      inp%BASSEL(ik, ib) = ibas
      write(unit=39, fmt='(2I12)') ik, ib
    end do
    close(unit=39)

  end subroutine


  subroutine initOlap(olap, inp, nb)
    implicit none
    type(overlap), intent(inout) :: olap
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: nb
    integer :: nmodes, nsw

    nmodes = inp%NMODES
    nsw = inp%NSW

    olap%NBANDS = nb
    olap%TSTEPS = inp%NSW
    olap%dt = inp%POTIM
    olap%NMODES = nmodes
    olap%COUPTYPE = inp%EPCTYPE
    olap%Np = inp%Np

    allocate(olap%Eig(nb, nsw-1))
    allocate(olap%gij(nb, nb, nmodes))
    allocate(olap%Phfreq(nb, nb, nmodes))

    if (olap%COUPTYPE==2) then
      allocate(olap%PhQ(nb, nb, nmodes, nsw-1))
      olap%PhQ = cero
    end if

    olap%gij = cero
    olap%Phfreq = 0.0_q
    olap%Eig = 0.0_q

  end subroutine


  subroutine copyToSec(olap, olap_sec, inp)
    implicit none
    type(namdInfo), intent(inout) :: inp
    type(overlap), intent(inout) :: olap_sec
    type(overlap), intent(in) :: olap

    integer :: ik, jk, ib, jb, nb, iBas, jBas
    integer :: kmin, kmax, bmin, bmax

    if (inp%LBASSEL) then
    ! If .TRUE., BMIN, BMAX, KMIN, KMAX, EMIN, EMAX are ignored!!!
      call readBasis(inp)
    else
    ! selected basises whose energies are between EMIN ~ EMAX,
    ! in the range BMIN ~ BMAX and KMIN ~ KMAX.
      call selBasis(inp, olap)
    end if

    if (inp%LSORT) call sortBasis(inp, olap)
    call initOlap(olap_sec, inp, inp%NBASIS)

    nb = inp%NBANDS
    kmin = inp%KMIN; kmax = inp%KMAX
    bmin = inp%BMIN; bmax = inp%BMAX

    do ik=kmin, kmax
      do ib=bmin, bmax

        if (inp%BASSEL(ik,ib)<0) cycle
        iBas = inp%BASSEL(ik,ib)
        olap_sec%Eig(iBas, :) = olap%Eig((ik-1)*nb+ib, :)

        do jk=kmin, kmax
          do jb=bmin, bmax

            if (inp%BASSEL(jk,jb)<0) cycle
            jBas = inp%BASSEL(jk,jb)

            ! olap_sec%Dij(iBas, jBas, :) = &
            ! olap%Dij( (ik-1)*nb+ib, (jk-1)*nb+jb, : )

            olap_sec%gij(iBas, jBas, :) = &
            olap%gij( (ik-1)*nb+ib, (jk-1)*nb+jb, : )

            olap_sec%Phfreq(iBas, jBas, :) = &
            olap%Phfreq( (ik-1)*nb+ib, (jk-1)*nb+jb, : )

            if (inp%EPCTYPE==2) &
            olap_sec%PhQ(iBas, jBas, :, :) = &
            olap%PhQ( (ik-1)*nb+ib, (jk-1)*nb+jb, :, : )

          end do
        end do

      end do
    end do

  end subroutine


  subroutine CoupToFileEP(olap)
    implicit none
    type(overlap), intent(in) :: olap

    ! Couplings are save to a binary file
    integer :: recordL, ierr, irec
    integer :: i, j, nb, im, nmodes, it
    ! to find out the record length
    complex(kind=q), allocatable, dimension(:)  :: values

    allocate(values(olap%NBANDS * olap%NBANDS))
    inquire (iolength=recordL) values
    deallocate(values)

    if (recordL<200) recordL=200 ! at least 25 real numbers (8 bytes)

    open(unit=30, file='EPCAR', access='direct', form='unformatted', &
         status='unknown', recl=recordL, iostat=ierr)
    if(ierr /= 0) then
        write(*,*) "File I/O error with EPCAR"
        stop
    end if

    irec = 1
    write(unit=30, rec=irec) &
        real(recordL,       kind=q), real(olap%NBANDS, kind=q), &
        real(olap%TSTEPS,   kind=q), real(olap%dt, kind=q), &
        real(olap%COUPTYPE, kind=q), real(olap%Np, kind=q), &
        real(olap%NMODES,   kind=q)

    nb = olap%NBANDS
    nmodes = olap%NMODES
    do im=1, nmodes
      irec = irec + 1
      ! Here we only save Eig(:,1)
      ! We suppose Eig(i,t) for any t does not change.
      write(unit=30, rec=irec) (olap%Eig(i,1), i=1,nb), &
          ((olap%Phfreq(i,j,im), i=1,nb), j=1,nb)
    end do
    do im=1, nmodes
      irec = irec + 1
      write(unit=30, rec=irec) ((olap%gij(i,j,im), i=1,nb), j=1,nb)
    end do
    do im=1, nmodes
      irec = irec + 1
      write(unit=30, rec=irec) ((olap%gij(i,j,im), i=1,nb), j=1,nb)
    end do
    if (olap%COUPTYPE==2) then
      do it=1,olap%TSTEPS-1
        do im=1, nmodes
          irec = irec + 1
          write(unit=30, rec=irec) ((olap%PhQ(i,j,im,it), i=1,nb), j=1,nb)
        end do
      end do
    end if

    close(unit=30)

  end subroutine


  subroutine CoupFromFileEP(olap)
    implicit none
    type(overlap), intent(inout) :: olap

    ! Couplings are save to a binary file
    integer :: irecordL, ierr, irec
    integer :: i, j, nb, im, nmodes, it
    real(kind=q) :: recordL, rnbands, rnsw, rdt, rctype, rnp, rnmodes

    open(unit=31, file='EPCAR', access='direct', form='unformatted', &
         status = 'unknown', recl=256, iostat=ierr)
    if(ierr /= 0) then
      write(*,*) "File I/O error with EPCAR"
      stop
    end if

    read(unit=31,rec=1) recordL, rnbands, rnsw, rdt, rctype, rnp, rnmodes
    ! write(*,*) recordL, rnbands, rnsw, rdt

    if (olap%NBANDS /= NINT(rnbands) .or. olap%TSTEPS /= NINT(rnsw) .or. &
       olap%COUPTYPE /= NINT(rctype) .or. olap%NMODES /= NINT(rnmodes)) then
      ! write(*,*) olap%NBANDS, NINT(rnbands), olap%TSTEPS, NINT(rnsw)
      write(*,*) "The EPCAR seems to be wrong..."
      stop
    end if

    close(31)

    irecordL = NINT(recordL)
    open(unit=31, file='EPCAR', access='direct', form='unformatted', &
         status = 'unknown', recl=irecordL, iostat=ierr)

    write(*,*) "Reading couplings from EPCAR..."

    irec = 1
    nb = olap%NBANDS
    nmodes = olap%NMODES
    olap%Np = rnp

    do im=1, nmodes
      irec = irec + 1
      ! We suppose Eig(i,t) for any t does not change.
      read(unit=31, rec=irec) (olap%Eig(i,1), i=1,nb), &
          ((olap%Phfreq(i,j,im), i=1,nb), j=1,nb)
      do i=1,nb
        olap%Eig(i,:) = olap%Eig(i,1)
      end do
    end do
    do im=1, nmodes
      irec = irec + 1
      read(unit=31, rec=irec) ((olap%gij(i,j,im), i=1,nb), j=1,nb)
    end do
    do im=1, nmodes
      irec = irec + 1
      read(unit=31, rec=irec) ((olap%gij(i,j,im), i=1,nb), j=1,nb)
    end do
    if (olap%COUPTYPE==2) then
      do it=1,olap%TSTEPS-1
        do im=1, nmodes
          irec = irec + 1
          read(unit=31, rec=irec) ((olap%PhQ(i,j,im,it), i=1,nb), j=1,nb)
        end do
      end do
    end if

    close(unit=31)

  end subroutine


  subroutine writeEP(olap)
    implicit none

    type(overlap), intent(inout) :: olap
    integer :: im, nmodes, ib, jb, nb, it, nsw, ierr

    nb = olap%NBANDS
    nsw = olap%TSTEPS
    nmodes = olap%NMODES

    open(unit=32, file='EPTXT', status='unknown', action='write')

    write(unit=32, fmt='(I6)') nmodes

    do im=1,nmodes
      do it=1,nsw-1
        write(unit=32, fmt='(*(f15.9))') &
          (( SUM(olap%EPcoup(ib,jb,im,:,it)), jb=1,nb ), ib=1,nb)
      end do
      write(unit=32,fmt=*)
    end do

    close(unit=32)

  end subroutine


  subroutine readEP(olap)
    implicit none

    type(overlap), intent(inout) :: olap
    integer :: im, nmodes, ib, jb, nb, it, nsw, ierr

    open(unit=33, file='EPTXT', status='unknown', &
         action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "EPTXT does NOT exist!"
      stop
    end if

    nb = olap%NBANDS
    nsw = olap%TSTEPS
    read(unit=33, fmt=*) nmodes

    do im=1,nmodes
      do it=1,nsw-1
        read(unit=33, fmt='(*(f15.9))') &
          (( olap%EPcoup(ib,jb,im,1,it), jb=1,nb ), ib=1,nb)
      end do
      read(unit=33,fmt=*)
    end do

    close(unit=33)

  end subroutine


end module epcoup
