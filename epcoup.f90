module epcoup
  use prec
  use constants
  use fileio
  use couplings
  use hdf5

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
    real(kind=q), allocatable, dimension(:,:) :: energy
    real(kind=q), allocatable, dimension(:,:) :: cellep
    real(kind=q), allocatable, dimension(:,:) :: kptsep
    real(kind=q), allocatable, dimension(:,:) :: qptsep
    real(kind=q), allocatable, dimension(:,:) :: freqep
    real(kind=q), allocatable, dimension(:,:) :: freqph
    real(kind=q), allocatable, dimension(:,:) :: qptsph
    complex(kind=q), allocatable, dimension(:,:,:,:,:) :: epmat
    complex(kind=q), allocatable, dimension(:,:,:,:) :: phmodes
    real(kind=q), allocatable, dimension(:,:) :: cellmd
    real(kind=q), allocatable, dimension(:,:,:) :: displ
    real(kind=q), allocatable, dimension(:,:,:) :: vel
    real(kind=q), allocatable, dimension(:,:,:) :: phproj
    ! Phonon projection of displacement in MD.
  end type

  contains

  subroutine releaseEPC(epc)
    implicit none

    type(epCoupling), intent(inout) :: epc

    if ( allocated(epc%R) ) deallocate(epc%R)
    if ( allocated(epc%atnum) ) deallocate(epc%atnum)
    if ( allocated(epc%qqmap) ) deallocate(epc%qqmap)
    if ( allocated(epc%kkqmap) ) deallocate(epc%kkqmap)
    if ( allocated(epc%energy) ) deallocate(epc%energy)
    if ( allocated(epc%cellep) ) deallocate(epc%cellep)
    if ( allocated(epc%kptsep) ) deallocate(epc%kptsep)
    if ( allocated(epc%qptsep) ) deallocate(epc%qptsep)
    if ( allocated(epc%freqep) ) deallocate(epc%freqep)
    if ( allocated(epc%freqph) ) deallocate(epc%freqph)
    if ( allocated(epc%qptsph) ) deallocate(epc%qptsph)
    if ( allocated(epc%epmat) ) deallocate(epc%epmat)
    if ( allocated(epc%phmodes) ) deallocate(epc%phmodes)
    if ( allocated(epc%cellmd) ) deallocate(epc%cellmd)
    if ( allocated(epc%displ) ) deallocate(epc%displ)
    if ( allocated(epc%vel) ) deallocate(epc%vel)
    if ( allocated(epc%phproj) ) deallocate(epc%phproj)

  end subroutine releaseEPC


  subroutine readEPC(inp, epc, olap)
    ! Read informations about e-p couplings from files in epc/.
    ! Folder epc/ include egnv (band energies), freq (phonon frequencies)
    ! and ephmat* (inteploted e-p matrix in dense k q mesh) files.

    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap

    integer :: ierr
    integer :: nm, nq, na
    integer :: bndmin, bndmax, nb
    integer :: nk1, nk2, nk3, nktot, nk
    integer :: i, j, ik, jk, ib, jb, iq, im
    integer :: ipool, npool, pool, lmpi
    integer :: ikf, nkf
    integer, allocatable :: nkq(:)
    real(kind=q) :: ef, kbT, phn, dE
    complex(kind=q) :: eptemp
    character(len=72) :: filinfo, filegnv, filfreq, filephmat, tag

    filinfo = trim(inp%FILEPC) // '/info'
    filegnv = trim(inp%FILEPC) // '/egnv'
    filfreq = trim(inp%FILEPC) // '/freq'

    open(unit=30, file=filinfo, action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "info file does NOT exist!"
      stop
    end if

    read(unit=30, fmt=*)
    read(unit=30, fmt=*) ef
    read(unit=30, fmt=*)
    read(unit=30, fmt=*) bndmin, bndmax, nk, nm, nq, na, npool, lmpi
    read(unit=30, fmt=*)
    allocate(nkq(npool))
    do ipool=1,npool
      read(unit=30, fmt=*) nkq(ipool)
    end do
    allocate(epc%cellep(na+3,3))
    read(unit=30, fmt=*)
    do i=1,3
      read(unit=30, fmt=*) (epc%cellep(i,j), j=1,3)
    enddo
    read(unit=30, fmt=*)
    do i=1,na
      read(unit=30, fmt=*) (epc%cellep(i+3,j), j=1,3)
    enddo

    nb = bndmax - bndmin + 1
    epc%nkpts  = nk
    epc%nbands = nb
    epc%nqpts  = nq
    epc%nmodes = nm
    epc%natepc = na

    close(30)
     
    open(unit=31, file=filegnv, action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "egnv file does NOT exist!"
      stop
    end if

    allocate(epc%kptsep(nk,3), epc%energy(nk,nb))
    read(unit=31, fmt=*)
    do ik=1,nk
      read(unit=31, fmt=*)
      read(unit=31, fmt=*) epc%kptsep(ik,:)
      do ib=1,nb
        read(unit=31, fmt=*) epc%energy(ik,ib)
        olap%Eig(ib+nb*(ik-1),:) = epc%energy(ik,ib)
      end do
    end do

    close(31)
     
    open(unit=32, file=filfreq, action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "freq file does NOT exist!"
      stop
    end if

    allocate(epc%qptsep(nq,3), epc%freqep(nq,nm))
    read(unit=32, fmt=*)
    do iq=1,epc%nqpts
      read(unit=32, fmt=*)
      read(unit=32, fmt=*) epc%qptsep(iq,:)
      do im=1,epc%nmodes
        read(unit=32, fmt=*) epc%freqep(iq,im)
      end do
    end do

    close(32)

    olap%Dij = cero
    allocate(epc%epmat(nb,nb,nk,nm,nq))
    do ipool=1,npool

      if (lmpi==0) then
        filephmat = trim(inp%FILEPC) // '/ephmat'
      else
        write(tag, *) ipool
        filephmat = trim(inp%FILEPC) // '/ephmat' // trim(adjustl(tag))
      end if

      open(unit=33, file=filephmat, action='read', iostat=ierr)
      if (ierr /= 0) then
        write(*,*) "ephmat file does NOT exist!"
        stop
      end if

      read(unit=33, fmt=*) nkf

      kbT = inp%TEMP * BOLKEV

      do ikf=1,nkq(ipool)

        read(unit=33, fmt=*) ik, jk, iq
        do im=1,epc%nmodes
          do jb=1,epc%nbands
            do ib=1,epc%nbands
              read(unit=33, fmt='(2ES20.10)') eptemp
              dE = epc%energy(ik,ib) - epc%energy(jk,jb) - epc%freqep(iq,im)/8065.541
              eptemp = eptemp * 0.015**2 / (dE**2 + 0.015**2) / PI
              epc%epmat(ib,jb,jk,im,iq) = eptemp
              phn = 1.0 / (exp(abs(epc%freqep(iq,im)/8065.541)/kbT)-1)
              eptemp = eptemp * (sqrt(phn) + sqrt(phn+1))
              olap%Dij(nb*(ik-1)+ib, jb+nb*(jk-1), :) &
              = olap%Dij(nb*(ik-1)+ib, jb+nb*(jk-1), :) + eptemp 
            end do
          end do
        end do
        
      end do

      close(33)

    end do

  end subroutine readEPC

  subroutine readEPCpert(inp, epc, olap)
    ! Read information of e-ph coupling from perturbo.x output file
    ! (prefix_ephmat_p1.h5), which is in form of hdf5.
    ! Extract informations include: el bands, ph dispersion, k & q lists, ephmat

    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap

    integer :: hdferror, info(4)
    integer :: nk, nq, nb, nm
    integer :: ib, jb, ik, jk, im, iq
    integer(hid_t) :: file_id, gr_id, dset_id
    integer(hsize_t) :: dim1(1), dim2(2), dim4(4)
    character(len=72) :: tagk, fname, grname, dsetname
    real(kind=q) :: kbT, phn, dE
    real(kind=q), allocatable, dimension(:,:) :: kqltemp
    real(kind=q), allocatable, dimension(:,:,:,:) :: eptemp_r, eptemp_i
    complex(kind=q), allocatable, dimension(:,:,:,:) :: eptemp
    complex :: eptemp_s

    call h5open_f(hdferror)
    fname = 'graphene_ephmat_p1.h5'
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
    ! write(*,*) info

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

      eptemp = ( eptemp_r + imgUnit * eptemp_i ) / 1000
      epc%epmat(:,:,ik,:,:) = eptemp

      do jk=1,nk
        iq = epc%kkqmap(ik, jk, 1)
        ! write(*,'(6F6.3)') epc%kptsep(ik,1:2), epc%kptsep(jk,1:2), epc%qptsep(iq,1:2)
        if (iq<0) exit
        do im=1,nm
          do jb=1,epc%nbands
            do ib=1,epc%nbands
              eptemp_s = eptemp(ib, jb, im, iq)
              dE = epc%energy(ib,ik) - epc%energy(jb,jk) - epc%freqep(im,iq)/1000
              ! write(*,'(8I4)') ik, ib, jk, jb, im, iq
              ! write(*,'(4F8.3)') dE, epc%energy(ib,ik), epc%energy(jb,jk), epc%freqep(im,iq)
              eptemp_s = eptemp_s * 0.015**2 / (dE**2 + 0.015**2) / PI
              if (epc%freqep(im,iq)>5) then
                phn = 1.0 / (exp(epc%freqep(im,iq)/1000/kbT)-1)
              else
                phn = 0.0
              end if
              ! write(*,'(2F8.3)') 0.015**2 / (dE**2 + 0.015**2) / PI, phn
              eptemp_s = eptemp_s * (sqrt(phn) + sqrt(phn+1))
              olap%Dij(nb*(ik-1)+ib, jb+nb*(jk-1), :) &
              = olap%Dij(nb*(ik-1)+ib, jb+nb*(jk-1), :) + eptemp_s
            end do
          end do
        enddo
      enddo

    enddo

    call h5gclose_f(gr_id, hdferror)

    call h5fclose_f(file_id, hdferror)

    do ib=1,nb*nk
      olap%Dij(ib,ib,:) = real(olap%Dij(ib,ib,:))
      do jb=ib+1,nb*nk
        olap%Dij(jb,ib,:) = CONJG(olap%Dij(ib,jb,:))
      enddo
    enddo

  end subroutine readEPCpert


  subroutine readEPCold(inp, epc)
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
    allocate(epc%cellep(na+3,3), &
             epc%kptsep(nk,3), &
             epc%qptsep(nq,3), &
             epc%energy(nk,nb), &
             epc%epmat(nb,nb,nk,nm,nq))

    read(unit=90, fmt=*)
    do i=1,3
      read(unit=90, fmt=*) (epc%cellep(i,j), j=1,3)
    enddo
    read(unit=90, fmt=*)
    do i=1,na
      read(unit=90, fmt=*) (epc%cellep(i+3,j), j=1,3)
    enddo

    read(unit=90, fmt=*)
    do i=1,nk
      read(unit=90, fmt=*) (epc%kptsep(i,j), j=1,3)
    enddo

    read(unit=90, fmt=*)
    do i=1,nq
      read(unit=90, fmt=*) (epc%qptsep(i,j), j=1,3)
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

  end subroutine readEPCold


  subroutine readPhmodes(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: ierr, i
    integer :: iq, im, ia, iax
    integer :: nqs, nmodes, nat
    character(len=24) :: charac, bra, ket
    complex(kind=q) :: temp
    real(kind=q), allocatable, dimension(:) :: qpt
    real(kind=q), allocatable, dimension(:,:) :: at

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
    allocate(epc%freqph(nqs, nmodes))
    allocate(qpt(3), at(3, 3))

    at = epc%cellep(1:3,1:3) / epc%cellep(1,1)
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
        read(unit=909, fmt=9011) charac, epc%freqph(iq, im)
        if (epc%freqph(iq,im)==0.0) epc%freqph(iq,im) = 0.00000001
        ! write(*,9011) charac, epc%freqph(iq, im)
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
    real(kind=q), allocatable, dimension(:) :: dr
    real(kind=q), allocatable, dimension(:) :: temp1, temp2

    type(epCoupling), intent(inout) :: epc

    allocate(epc%R(epc%natmd,3), epc%atnum(epc%natmd), &
             N(3), dr(3), temp1(3), temp2(3))

    do iax=1,3
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
    complex(kind=q) :: Qn, dQn
    ! The potential and kinetic energies of the normal mode.
    real(kind=q) :: Eu, Ek
    real(kind=q) :: theta
    real(kind=q), allocatable, dimension(:) :: displ, vel
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
          Qn = cero; dQn = cero
          do ia=1,epc%natmd
            displ = 0.0; vel = 0.0
            do iax=1,3
              displ = displ + epc%displ(t, ia, iax) * epc%cellmd(iax, :) 
              vel = vel + epc%vel(t, ia, iax) * epc%cellmd(iax, :) 
            end do
            theta = 2 * PI * DOT_PRODUCT( epc%qptsph(iq,:), epc%R(ia,:) )
            Qn  =  Qn + EXP(imgUnit*theta) * DOT_PRODUCT( &
                      CONJG(epc%phmodes(iq, im, epc%atnum(ia), :)), displ )
            dQn = dQn + EXP(imgUnit*theta) * DOT_PRODUCT( &
                      CONJG(epc%phmodes(iq, im, epc%atnum(ia), :)), vel )
          end do
          ! mass of C = 12.011
          Eu = 0.5 * 12.011 * epc%freqph(iq, im)**2 * CONJG(Qn) * Qn / epc%natmd &
            * 1.66 * 6.2415 / 100000
          Ek = 0.5 * 12.011 * CONJG(dQn) * dQn / epc%natmd &
            * 1.66 * 6.2415 * 10
          epc%phproj(t, im, iq) = (Eu + Ek) / ( hbar * epc%freqph(iq, im) / 1000 )
        end do
      end do
    end do

    call savePhp(inp, epc)

  end subroutine phDecomp


  subroutine savePhp(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(in) :: epc

    character(len=72) :: filename
    integer :: t, iq, im

    filename = 'phproj.dat'
    open(unit=38, file=filename, status='unknown', action='write')

    do t=1, inp%NSW
      write(unit=38, fmt='(*(f12.6))') &
          ((epc%phproj(t, im, iq), iq=1,epc%nqpts), im=1,epc%nmodes)
    end do

  end subroutine savePhp


  subroutine kqMatch(epc)
    implicit none

    type(epCoupling), intent(inout) :: epc

    real(kind=q) :: norm
    real(kind=q) :: dkq1(3), dkq2(3), dq(3)
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
          dkq1 = epc%kptsep(ik,:) - epc%kptsep(jk,:) - epc%qptsep(iq,:)
          dkq2 = epc%kptsep(ik,:) - epc%kptsep(jk,:) + epc%qptsep(iq,:)
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

    if (allocated(epc%qptsph)) then
      do iq=1,epc%nqpts
        do jq=1,epc%nqpts
          dq = epc%qptsep(iq,:) - epc%qptsph(jq,:)
          do iax=1,3
            dq(iax) = ABS(dq(iax)-NINT(dq(iax)))
          end do
          if (SUM(dq)<norm) then
            epc%qqmap(iq) = jq
            exit
          end if
        end do
      end do
    end if

  end subroutine kqMatch


  subroutine TDepCoupIJ(olap, olap_sec, inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
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
    olap%NBANDS = inp%NBANDS * inp%NKPOINTS
    olap%TSTEPS = inp%NSW
    olap%dt = inp%POTIM
    allocate(olap%Dij(olap%NBANDS, olap%NBANDS, olap%TSTEPS-1))
    allocate(olap%Eig(olap%NBANDS, olap%TSTEPS-1))
    olap%Dij = cero
    olap%EIg = cero

    olap_sec%NBANDS = inp%NBASIS * inp%NKSEL
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

    else if (inp%EPCTYPE==1) then

      write(*,'(A)') "------------------------------------------------------------"
      write(*,*) "Calculating couplings from dense e-p matrix..."

      call readEPCpert(inp, epc, olap)
      call copyToSec(olap, olap_sec, inp)
      call CoupToFile(olap)
      call writeNaEig(olap_sec, inp)

    else

      write(*,'(A)') "------------------------------------------------------------"
      write(*,*) "Calculating couplings from e-p matrix..."

      ! call readEPCold(inp, epc)
      call readEPC(inp, epc, olap)
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
                 
                  temp = temp / SQRT(1.0 * epc%natmd / epc%natepc)
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

    call releaseEPC(epc)
    deallocate(olap%Dij, olap%Eig)

  end subroutine TDepCoupIJ


  subroutine copyToSec(olap, olap_sec, inp)
    implicit none
    type(overlap), intent(inout) :: olap_sec
    type(overlap), intent(in) :: olap
    type(namdInfo), intent(in) :: inp

    integer :: ik, jk, NB, NBas

    NB = inp%NBANDS
    NBas = inp%NBASIS

    do ik=inp%KMIN,inp%KMAX
      do jk=inp%KMIN,inp%KMAX
        olap_sec%Dij( (ik-inp%KMIN)*NBas+1:(ik-inp%KMIN+1)*NBas, &
                      (jk-inp%KMIN)*NBas+1:(jk-inp%KMIN+1)*NBas, : ) = &
            olap%Dij( (ik-1)*NB+inp%BMIN:(ik-1)*NB+inp%BMAX, &
                      (jk-1)*NB+inp%BMIN:(jk-1)*NB+inp%BMAX, : )
      end do
      olap_sec%Eig( (ik-inp%KMIN)*NBas+1:(ik-inp%KMIN+1)*NBas, : ) = &
          olap%Eig( (ik-1)*NB+inp%BMIN:(ik-1)*NB+inp%BMAX, : )
    end do

  end subroutine copyToSec


end module epcoup
