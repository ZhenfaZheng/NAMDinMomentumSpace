module epcoup
  use prec
  use constants
  use fileio
  use couplings
  use hdf5
  use mpi

  implicit none

  type epCoupling
    integer :: nbands, nkpts, nmodes, nqpts, natepc, natmd, ipart
    integer, allocatable, dimension(:) :: nkpts_ps
    integer, allocatable, dimension(:) :: atmap
    !! mapping of each atom in MD cell to atom in unit cell for epc calculation.
    integer, allocatable, dimension(:,:) :: Rp
    !! number of lattice vector of each atom in MD cell, Rp(natmd,3)
    integer, allocatable, dimension(:,:) :: kkqmap
    !! map ik & jk of electronic states to q of phonon modes.
    real(kind=q), allocatable, dimension(:) :: mass, eig
    real(kind=q), allocatable, dimension(:,:,:) :: displ
    real(kind=q), allocatable, dimension(:,:) :: kpts, qpts
    real(kind=q), allocatable, dimension(:,:) :: elen, phfreq
    real(kind=q), allocatable, dimension(:,:) :: cellep, cellmd
    real(kind=q), allocatable, dimension(:,:,:) :: PhQ
    complex(kind=q), allocatable, dimension(:,:,:) :: eiwdt
    complex(kind=q), allocatable, dimension(:,:,:,:) :: phmodes
    real(kind=q), allocatable, dimension(:,:,:,:) :: epcec
    complex(kind=q), allocatable, dimension(:,:,:) :: gij
  end type

  contains


  subroutine TDepCoupIJ(olap, olap_sec, inp, epc)
    implicit none

    type(namdInfo), intent(inout) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap
    type(overlap), intent(inout) :: olap_sec

    integer :: nmodes, nb, nsw
    integer :: irank, nrank, ierr
    logical :: lcoup

    call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nrank, ierr)

    if (irank==0) write(*,*)
    if (irank==0) write(*,'(A)') &
      "------------------------------------------------------------"

    if (inp%EPCTYPE==1) then
      if (irank==0) write(*,*) "TypeI e-ph coupling calculation."
    else
      if (irank==0) write(*,*) "TypeII e-ph coupling calculation."
    end if

    if (irank==0) write(*,*) "Reading el & ph information."
    call readBasicInfo(inp, epc)
    call initEPC(inp, epc)
    call readEPHinfo(inp, epc)

    call selBasis(inp, epc)
    if (inp%LSORT) call sortBasis(inp, epc)
    if (irank==0) call outputBAS(inp, epc)

    call mpi_split_bas(inp)
    call initEPC2(inp, epc)

    call initOlap(olap_sec, inp, epc%nqpts, inp%NBASIS)
    call epcToOlapSec(inp, epc, olap_sec)

    if (irank==0) write(*,*) "Mapping k & k' points with q point."
    call kqMatchSec(inp, epc, olap_sec)
    if (irank==0) write(*,*) "Reading e-ph matrix elements."
    call readEPHmatSec(inp, epc, olap_sec)
    call HemitGij(inp, epc)

    if (inp%EPCTYPE==1) then
      if (irank==0) write(*,*) "Calculating TD phonon normal modes."
      call symPhfreq(inp, epc)
      call calcPhQ(olap_sec, inp, epc)
    else
      if (irank==0) write(*,*) "Decomposing phonon modes from MD traj."
      call phDecomp(inp, olap_sec, epc)
    end if

    ! call savePhQ(olap_sec)

    if (irank==0) write(*,*) "Calculating e-ph couplings."
  ! if (inp%LARGEBS) then
      call calcEPC(olap_sec, epc, inp)
      call writeEPELTXT(inp, epc)
      call writeEPPHTXT(inp, epc, olap_sec)
      ! call writeTXT_LBS(olap_sec)
      ! call writeTXT_LBS_Mode(olap_sec)
      ! call writeKKQMap(olap_sec)
  ! else
  !   call calcNAC(olap_sec, inp)
  !   call writeTXT(olap_sec)
  !   deallocate(olap_sec%Dij, olap_sec%EPcoup)
  !   call calcEPC(olap_sec, epc, inp)
  ! end if

  ! call releaseEPC(epc)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (irank==0) write(*,*) "Done!"
    if (irank==0) write(*,'(A)') &
      "------------------------------------------------------------"
    if (irank==0) write(*,*)

  end subroutine


  subroutine initEPC(inp, epc)
    implicit none

    type(namdInfo), intent(inout) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: nk, nq, nb, nm, nat

    nk  = SUM(epc%nkpts_ps) ; epc%nkpts = nk
    if (inp%NKPOINTS .NE. epc%nkpts) then
      write(*,*) 'NKPOINTS seems to be wrong!'
      stop
    end if

    nq  = epc%nqpts
    nb  = epc%nbands
    nm  = epc%nmodes
    nat = epc%natepc

    allocate(epc%mass(nat))
    allocate(epc%cellep(nat+3,3))
    allocate(epc%kpts(nk, 3))
    allocate(epc%elen(nk, nb))
    allocate(epc%qpts(nq, 3))
    allocate(epc%phfreq(nq, nm))
    allocate(epc%PhQ(nq, nm, 2))
    allocate(epc%eiwdt(nq, nm, 2))
    allocate(epc%phmodes(nq, nm, nat, 3))

  end subroutine


  subroutine initEPC2(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    integer :: ibas, nbas, nbas_p, nmodes, ik, ib

    nbas = inp%NBASIS
    nbas_p = inp%NBASIS_P
    nmodes = epc%nmodes

    allocate(epc%eig(nbas))
    allocate(epc%gij(nbas_p, nbas, nmodes))
    allocate(epc%epcec(nbas_p, nbas, nmodes, 2))

    epc%gij = cero
    epc%epcec = 0.0_q

    do ibas=1,nbas
      ik = inp%BASLIST(ibas,1)
      ib = inp%BASLIST(ibas,2)
      epc%eig(ibas) = epc%elen(ik, ib)
    end do

  end subroutine


  subroutine initOlap(olap, inp, nq, nb)
    implicit none
    type(overlap), intent(inout) :: olap
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: nq, nb
    integer :: nb_p, nmodes, nsw

    nb_p = inp%NBASIS_P
    nmodes = inp%NMODES
    nsw = inp%NSW

    olap%NBANDS = nb
    olap%TSTEPS = inp%NSW
    olap%dt = inp%POTIM
    olap%NMODES = nmodes
    olap%COUPTYPE = inp%EPCTYPE
    olap%Np = inp%Np
    olap%NQ = nq

    allocate(olap%Eig(nb, nsw-1))
    allocate(olap%gij(nb_p, nb, nmodes))
    ! allocate(olap%Phfreq(nq, nmodes))
    ! allocate(olap%PhQ(nq, nmodes, 2, nsw-1))
    allocate(olap%kkqmap(nb,nb))

    olap%gij = cero
    ! olap%PhQ = cero
    ! olap%Phfreq = 0.0_q
    olap%Eig = 0.0_q
    olap%kkqmap = -1

  end subroutine


  subroutine releaseEPC(epc)
    implicit none

    type(epCoupling), intent(inout) :: epc

    if ( allocated(epc%Rp) ) deallocate(epc%Rp)
    if ( allocated(epc%atmap) ) deallocate(epc%atmap)
    if ( allocated(epc%kpts) ) deallocate(epc%kpts)
    if ( allocated(epc%qpts) ) deallocate(epc%qpts)
    if ( allocated(epc%elen) ) deallocate(epc%elen)
    if ( allocated(epc%phfreq) ) deallocate(epc%phfreq)
    if ( allocated(epc%cellep) ) deallocate(epc%cellep)
    if ( allocated(epc%cellmd) ) deallocate(epc%cellmd)
    if ( allocated(epc%displ) ) deallocate(epc%displ)
    if ( allocated(epc%phmodes) ) deallocate(epc%phmodes)

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
    allocate(inp%BASLIST(inp%NBASIS, 3))
    do ibas=1,inp%NBASIS
      read(unit=38, fmt=*) ik, ib
      inp%BASSEL(ik, ib) = ibas
      inp%BASLIST(ibas, 1) = ik
      inp%BASLIST(ibas, 2) = ib
    end do

    close(unit=38)

  end subroutine


  subroutine selBasis(inp, epc)
    implicit none

    type(namdInfo), intent(inout) :: inp
    type(epCoupling), intent(in) :: epc

    integer :: ik, ib, nb, ibas
    real :: emin, emax

    if (inp%LBASSEL) then
      call readBasis(inp)
      return
    end if

    nb = inp%NBANDS
    emin = inp%EMIN
    emax = inp%EMAX

    inp%NBASIS = 0
    inp%BASSEL = -1

    do ik=inp%KMIN, inp%KMAX
      do ib=inp%BMIN, inp%BMAX
        if ( epc%elen(ik,ib)>emin .and. epc%elen(ik,ib)<emax ) then
          inp%NBASIS = inp%NBASIS + 1
          inp%BASSEL(ik,ib) = inp%NBASIS
        end if
      end do
    end do

    allocate(inp%BASLIST(inp%NBASIS, 3))

    ibas = 0
    do ik=inp%KMIN, inp%KMAX
      do ib=inp%BMIN, inp%BMAX
        if (inp%BASSEL(ik,ib)>0) then
          ibas = ibas + 1
          inp%BASLIST(ibas, 1) = ik
          inp%BASLIST(ibas, 2) = ib
        end if
      end do
    end do

    !! For array element inp%BASSEL(ik,ib),
    !! 0 represent ik,ib state isn't selected as basis;
    !! number >0 represent ik,ib state is selected as basis,
    !! and the number is the state's serial number among the basises.

  end subroutine


  subroutine sortBasis(inp, epc)
    implicit none

    type(namdInfo), intent(inout) :: inp
    type(epCoupling), intent(in) :: epc

    integer :: ibas, i, j, ik, ib, nb, nbas
    integer :: temp(3)

    nb = inp%NBANDS
    nbas = inp%NBASIS
    do ibas=1,nbas
      ik = inp%BASLIST(ibas,1)
      ib = inp%BASLIST(ibas,2)
      inp%BASLIST(ibas,3) = epc%elen(ik, ib) * 1000
    end do

    do i=1,nbas-1
      do j=1,nbas-i
        if ( inp%BASLIST(j+1,3) < inp%BASLIST(j,3) ) then
          temp = inp%BASLIST(j,:)
          inp%BASLIST(j,:) = inp%BASLIST(j+1,:)
          inp%BASLIST(j+1,:) = temp
        end if
      end do
    end do

    do ibas=1,nbas
      ik = inp%BASLIST(ibas,1)
      ib = inp%BASLIST(ibas,2)
      inp%BASSEL(ik, ib) = ibas
    end do

  end subroutine


  subroutine elen2eig(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    integer :: ibas, nbas, ik, ib

    nbas = inp%NBASIS
    allocate(epc%eig(nbas))

    do ibas=1,nbas
      ik = inp%BASLIST(ibas,1)
      ib = inp%BASLIST(ibas,2)
      epc%eig(ibas) = epc%elen(ik, ib)
    end do

  end subroutine


  subroutine outputBAS(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(in) :: epc
    integer :: ibas, nbas, ik, ib

    nbas = inp%NBASIS

    open(unit=39, file='BASSEL', status='unknown', action='write')
    write(unit=39, fmt='(I18)') nbas

    do ibas=1,nbas
      ik = inp%BASLIST(ibas,1)
      ib = inp%BASLIST(ibas,2)
      write(unit=39, fmt='(2I12, F20.8)') ik, ib, epc%elen(ik,ib)
    end do

    close(unit=39)

  end subroutine


  subroutine readBasicInfo(inp, epc)
    implicit none

    type(namdInfo), intent(inout) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: ip, nparts, nk

    nk = inp%NKPOINTS
    nparts = inp%NPARTS
    allocate(epc%nkpts_ps(nparts))

    do ip=1, nparts
      epc%ipart = ip
      call readBasicInfo_p(inp, epc)
    end do

  end subroutine


  subroutine readEPHinfo(inp, epc)
    implicit none

    type(namdInfo), intent(inout) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: ip, nparts, nk

    nparts = inp%NPARTS
    do ip=1, nparts
      epc%ipart = ip
      call readEPHinfo_p(inp, epc)
    end do

  end subroutine


  subroutine readEPHmatSec(inp, epc, olap_sec)
    implicit none

    type(namdInfo), intent(inout) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap_sec

    integer :: ip, nparts, nk

    nparts = inp%NPARTS
    do ip=1, nparts
      epc%ipart = ip
      call readEPHmatSec_p(inp, epc, olap_sec)
    end do

  end subroutine


  subroutine readBasicInfo_p(inp, epc)
    implicit none

    type(namdInfo), intent(inout) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: hdferror, info(4)
    integer(hsize_t) :: dim1(1)
    integer :: nq, nb, nm, nat
    integer(hid_t) :: file_id, gr_id, dset_id
    character(len=256) :: fname, grname, dsetname
    character(len=256) :: epmdir, prefix, iptag

    epmdir = inp%EPMDIR
    prefix = inp%EPMPREF
    write(iptag, '(I8)') epc%ipart
    fname = trim(epmdir) // trim(prefix) // '_ephmat_p' &
         // trim( adjustl(iptag) ) // '.h5'

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

    if(allocated(epc%nkpts_ps)) epc%nkpts_ps(epc%ipart) = info(1)

    nq = info(2); epc%nqpts  = nq
    nb = info(3); epc%nbands = nb
    nm = info(4); epc%nmodes = nm
    nat = nm/3  ; epc%natepc = nat

    if (inp%NBANDS .NE. nb) then
      write(*,*) 'NBANDS seems to be wrong!'
      stop
    end if
    if ( (inp%NMODES .NE. 1) .AND. (inp%NMODES .NE. nm) ) then
      write(*,*) 'NMODES seems to be wrong!'
      stop
    end if
    if ( inp%Np==1 ) inp%Np = nq
    inp%NQPOINTS = nq
    inp%NMODES = nm

  end subroutine


  subroutine readEPHinfo_p(inp, epc)
    ! Read information of e-ph coupling from perturbo.x output file
    ! (prefix_ephmat_p1.h5), which is in form of hdf5.
    ! Extract informations include: el bands, ph dispersion, k & q lists, ephmat

    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: hdferror
    integer :: nq, nb, nm, nat
    integer :: ipart, kst, kend, nk
    integer :: ib, jb, ik, jk, im, iq, it, ibas, jbas
    integer(hid_t) :: file_id, gr_id, dset_id
    integer(hsize_t) :: dim1(1), dim2(2), dim4(4)
    character(len=72) :: fname, grname, dsetname
    character(len=256) :: epmdir, prefix, iptag
    real(kind=q), allocatable, dimension(:,:) :: entemp, freqtemp
    real(kind=q), allocatable, dimension(:,:) :: kqltemp, pos, lattvec
    real(kind=q), allocatable, dimension(:,:,:,:) :: eptemp_r, eptemp_i, phmtemp
    complex(kind=q), allocatable, dimension(:,:,:,:) :: eptemp

    nq  = epc%nqpts
    nb  = epc%nbands
    nm  = epc%nmodes
    nat = epc%natepc

    ipart = epc%ipart
    if (ipart==1) then
      kst = 1
    else
      kst = SUM(epc%nkpts_ps(1:ipart-1)) + 1
    end if
    nk  = epc%nkpts_ps(ipart)
    kend = kst + nk - 1

    epmdir = inp%EPMDIR
    prefix = inp%EPMPREF
    write(iptag, '(I8)') ipart
    fname = trim(epmdir) // trim(prefix) // '_ephmat_p' &
         // trim( adjustl(iptag) ) // '.h5'

    call h5open_f(hdferror)
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, hdferror)

    ! Read el bands, ph dispersion, k & q lists.
    grname = 'el_ph_band_info'
    call h5gopen_f(file_id, grname, gr_id, hdferror)

    dsetname = 'k_list'
    allocate(kqltemp(3,nk))
    dim2 = shape(kqltemp, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, kqltemp, dim2, hdferror)
    call h5dclose_f(dset_id, hdferror)
    epc%kpts(kst:kend, :) = transpose(kqltemp)
    deallocate(kqltemp)

    dsetname = 'q_list'
    allocate(kqltemp(3,nq))
    dim2 = shape(kqltemp, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, kqltemp, dim2, hdferror)
    call h5dclose_f(dset_id, hdferror)
    epc%qpts = transpose(kqltemp)
    deallocate(kqltemp)

    dsetname = 'el_band_eV'
    allocate(entemp(nb,nk))
    dim2 = shape(entemp, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, entemp, dim2, hdferror)
    call h5dclose_f(dset_id, hdferror)
    epc%elen(kst:kend, :) = transpose(entemp)

    dsetname = 'ph_disp_meV'
    allocate(freqtemp(nm, nq))
    dim2 = shape(freqtemp, kind=hsize_t)
    call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, freqtemp, dim2, hdferror)
    call h5dclose_f(dset_id, hdferror)
    freqtemp = freqtemp / 1000.0_q ! transform unit to eV
    epc%phfreq = transpose(freqtemp)


  ! if (inp%EPCTYPE==2) then
      dsetname = 'mass_a.u.'
      dim1 = shape(epc%mass, kind=hsize_t)
      call h5dopen_f(gr_id, dsetname, dset_id, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, epc%mass, dim1, hdferror)
      call h5dclose_f(dset_id, hdferror)

      dsetname = 'lattice_vec_angstrom'
      allocate(lattvec(3,3))
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
  ! end if

    call h5gclose_f(gr_id, hdferror)

    call h5fclose_f(file_id, hdferror)

  end subroutine


  subroutine readEPHmat_p(inp, epc, olap)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap

    integer :: hdferror
    integer :: nq, nb, nm
    integer :: ipart, kst, nk, nk_tot
    integer :: ib, jb, ik, jk, im, iq, it, ibas, jbas
    integer(hid_t) :: file_id, gr_id, dset_id
    integer(hsize_t) :: dim1(1), dim2(2), dim4(4)
    character(len=72) :: tagk, fname, grname, dsetname
    character(len=256) :: epmdir, prefix, iptag
    real(kind=q), allocatable, dimension(:,:) :: entemp
    real(kind=q), allocatable, dimension(:,:) :: kqltemp, pos, lattvec
    real(kind=q), allocatable, dimension(:,:,:,:) :: eptemp_r, eptemp_i, phmtemp
    complex(kind=q), allocatable, dimension(:,:,:,:) :: eptemp

    nq  = epc%nqpts
    nb  = epc%nbands
    nm  = epc%nmodes

    ipart = epc%ipart
    if (ipart==1) then
      kst = 1
    else
      kst = SUM(epc%nkpts_ps(1:ipart-1)) + 1
    end if
    nk  = epc%nkpts_ps(ipart)
    nk_tot = epc%nkpts

    epmdir = inp%EPMDIR
    prefix = inp%EPMPREF
    write(iptag, '(I8)') epc%ipart
    fname = trim(epmdir) // trim(prefix) // '_ephmat_p' &
         // trim( adjustl(iptag) ) // '.h5'

    call h5open_f(hdferror)
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, hdferror)

    ! Read e-ph matrix for every k.
    grname = 'g_ephmat_total_meV'
    call h5gopen_f(file_id, grname, gr_id, hdferror)
    allocate(eptemp(nb, nb, nm, nq))
    allocate(eptemp_r(nb, nb, nm, nq))
    allocate(eptemp_i(nb, nb, nm, nq))

    do ik=1, nk

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

      do jk=1, nk_tot

        iq = epc%kkqmap(jk, ik+kst-1)
        ! kpts[jk] = kpts[ik] + qpts[iq]
        if (iq<0) cycle

        do im=1,nm

          if (epc%phfreq(iq,im)<inp%PHCUT) cycle

          do jb=1,nb
            do ib=1,nb

              ibas = nb*(ik+kst-2)+ib
              jbas = nb*(jk-1)+jb
              olap%kkqmap(jbas, ibas) = iq
              ! olap%gij(jbas, ibas, im) = eptemp(jb, ib, im, iq)
              epc%gij(jbas, ibas, im) = eptemp(jb, ib, im, iq)
              ! ibas ~ initial state; jbas ~ final state
              ! gij(jbas, ibas, im) = < psi_{jb,jk} | d_{iq,im} V | psi_{ib,ik} >
              ! kpts[jk] = kpts[ik] + qpts[iq]

            end do ! ib loop
          end do ! jb loop

        enddo ! im loop

      enddo ! jk loop

    enddo ! ik loop

    call h5gclose_f(gr_id, hdferror)

    call h5fclose_f(file_id, hdferror)

  end subroutine


  subroutine readEPHmatSec_p(inp, epc, olap_sec)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc
    type(overlap), intent(inout) :: olap_sec

    integer :: hdferror
    integer :: nq, nb, nm
    integer :: jrank, jst, ierr
    integer :: ipart, kst, kend, nk, nk_tot
    integer :: ib, jb, ik, jk, im, iq, it, ibas, jbas
    integer(hid_t) :: file_id, gr_id, dset_id
    integer(hsize_t) :: dim1(1), dim2(2), dim4(4)
    character(len=72) :: tagk, fname, grname, dsetname
    character(len=256) :: epmdir, prefix, iptag
    real(kind=q), allocatable, dimension(:,:) :: entemp
    real(kind=q), allocatable, dimension(:,:) :: kqltemp, pos, lattvec
    real(kind=q), allocatable, dimension(:,:,:,:) :: eptemp_r, eptemp_i, phmtemp
    complex(kind=q), allocatable, dimension(:,:,:,:) :: eptemp

    nq  = epc%nqpts
    nb  = epc%nbands
    nm  = epc%nmodes

    ipart = epc%ipart
    if (ipart==1) then
      kst = 1
    else
      kst = SUM(epc%nkpts_ps(1:ipart-1)) + 1
    end if
    nk  = epc%nkpts_ps(ipart)
    kend = kst + nk - 1
    nk_tot = epc%nkpts

    epmdir = inp%EPMDIR
    prefix = inp%EPMPREF
    write(iptag, '(I8)') epc%ipart
    fname = trim(epmdir) // trim(prefix) // '_ephmat_p' &
         // trim( adjustl(iptag) ) // '.h5'

    call h5open_f(hdferror)
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, hdferror)

    ! Read e-ph matrix for every k.
    grname = 'g_ephmat_total_meV'
    call h5gopen_f(file_id, grname, gr_id, hdferror)
    allocate(eptemp(nb, nb, nm, nq))
    allocate(eptemp_r(nb, nb, nm, nq))
    allocate(eptemp_i(nb, nb, nm, nq))

    call MPI_COMM_RANK(MPI_COMM_WORLD, jrank, ierr)
    jst = inp%ISTS(jrank+1)

    do ik=kst, kend

      if ( SUM(inp%BASSEL(ik,:)) == -nb ) cycle

      write(tagk, '(I8)') ik-kst+1

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

      do ib=1,nb
        ibas = inp%BASSEL(ik,ib)
        if (ibas==-1) cycle

        do jk=1, nk_tot
          do jb=1,nb

            jbas = inp%BASSEL_P(jk,jb)
            if (jbas==-1) cycle
            iq = epc%kkqmap(jbas, ibas)
            ! kpts[jk] = kpts[ik] + qpts[iq]
            if (iq<0) cycle

            do im=1,nm
              if (epc%phfreq(iq,im)<inp%PHCUT) cycle
              ! olap_sec%gij(jbas-jst+1, ibas, im) = eptemp(jb, ib, im, iq)
              epc%gij(jbas-jst+1, ibas, im) = eptemp(jb, ib, im, iq)
              ! ibas ~ initial state; jbas ~ final state
              ! gij(jbas, ibas, im) = < psi_{jb,jk} | d_{iq,im} V | psi_{ib,ik} >
              ! kpts[jk] = kpts[ik] + qpts[iq]
            end do ! im loop

          end do ! jb loop
        enddo ! jk loop

      enddo ! ib loop

    enddo ! ik loop

    call h5gclose_f(gr_id, hdferror)

    call h5fclose_f(file_id, hdferror)

  end subroutine


  subroutine HemitGij(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: crank, irank, jrank, xrank, yrank, nrank, ierr, tag
    integer :: ist, iend, jst, jend, xst, xend, yst, yend
    integer :: ibas, jbas, nbas, nbas_p, im, nmodes
    complex(kind=q), allocatable :: gtemp(:,:), gtemp2(:,:), gij(:), gji(:), gabs(:), theta(:)

    nbas = inp%NBASIS
    nbas_p = inp%NBASIS_P
    nmodes = epc%nmodes
    allocate(gij(nmodes), gji(nmodes))
    allocate(gtemp(nbas, nmodes))
    allocate(gtemp2(nbas, nmodes))

    call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nrank, ierr)

    do jrank = 0, nrank-1

      jst = inp%ISTS(jrank+1)
      jend = inp%IENDS(jrank+1)

      do jbas = jst, jend

        if (irank == jrank) then
          gtemp = epc%gij(jbas-jst+1, :, :)
        end if

        do im=1,nmodes
          call MPI_Bcast(gtemp(:,im), nbas, MPI_DOUBLE_COMPLEX, jrank, MPI_COMM_WORLD, ierr)
        end do
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        ist = inp%ISTS(irank+1)
        iend = inp%IENDS(irank+1)

        do ibas = ist, iend
          if (ibas<jbas) cycle
          gji = gtemp(ibas, :)
          gij = epc%gij(ibas-ist+1, jbas, :)
          ! here hermit gij and gji
          if (ibas == jbas) then
            gij = ABS(gij); gji = ABS(gij)
          else
            gji = CONJG(gij)
          end if
          epc%gij(ibas-ist+1, jbas, :) = gij
          gtemp(ibas,:) = gji
          ! print *, 'current at', irank, 'ij is', ibas, jbas, 'gij', gij(5)
        end do

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        do im=1,nmodes
          if ((jrank+1) > (nrank-1)) cycle
          do xrank = jrank+1, nrank-1
            if (irank == xrank) then
              tag = xrank
              call MPI_SEND(gtemp(ist:iend, im), nbas_p, MPI_DOUBLE_COMPLEX, &
                            jrank, tag, MPI_COMM_WORLD, ierr)
            else if (irank == jrank) then
              tag = xrank
              xst = inp%ISTS(xrank+1)
              xend = inp%IENDS(xrank+1)
              call MPI_RECV(gtemp(xst:xend, im), xend-xst+1, MPI_DOUBLE_COMPLEX, &
                            xrank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            end if
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          end do
        end do

        if (irank==jrank) then
          epc%gij(jbas-jst+1,jbas:nbas, :) = gtemp(jbas:nbas, :)
          ! print *, 'current at', irank, 'jbas is', jbas
          ! print *, 'gji', epc%gij(jbas-jst+1, :, 5)
        end if

      end do

    end do


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

    open(unit=35, file=inp%FILMD, status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "XDATCAR file does NOT exist!"
      stop
    end if

    read(unit=35, fmt=*)
    read(unit=35, fmt=*) scal
    do i=1,3
      read(unit=35, fmt=*) (lattvec(i, iax), iax=1,3)
      lattvec(i,:) = lattvec(i,:) * scal
    end do
    read(unit=35, fmt=*)
    read(unit=35, fmt=*) nat
    epc%natmd = nat

    allocate(epc%displ(nsw, nat, 3))

    do it=1,nsw
      read(unit=35, fmt=*)
      do ia=1,nat
        read(unit=35, fmt=*) (epc%displ(it, ia, iax), iax=1, 3)
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
          ! change from crystal coordinates to cartesian coordinates.
          temp = temp + ( epc%displ(it,ia,iax) - epc%cellmd(ia+3, iax) ) * &
                 lattvec(iax, :)
        end do
        epc%displ(it,ia,:) = temp
      end do
    end do

    close(unit=35)

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


  subroutine saveXDAT(inp, olap, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(overlap), intent(in) :: olap
    type(epCoupling), intent(in) :: epc

    integer :: ierr, i, j, it
    integer :: nat, nsw, nx, ny, nz
    real(kind=q) :: cell(3,3)
    character(40) :: FORM

    nx = 9; ny = 9; nz = 1
    nsw = inp%NAMDTIME / inp%POTIM
    nat = epc%natepc * nx * ny * nz

    cell = epc%cellep(1:3,:)
    cell(:,1) = cell(:,1) * nx
    cell(:,2) = cell(:,2) * ny
    cell(:,3) = cell(:,3) * nz

    open(unit=310, file='XDATCAR', status='unknown', &
         action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "XDATCAR file I/O error!"
      stop
    end if

    ! Write out initial header for XDATCAR
    write(unit=310, fmt='(A40)') inp%EPMPREF
    write(unit=310, fmt='(I12)') 1 ! scale
    write(unit=310, fmt='(1X,3F12.6)') ((cell(i,j), j=1,3), i=1,3)
    write(unit=310, fmt='(*(A8))') 'C'
    write(unit=310, fmt='(*(I8))') 162

    do it=1,nsw
      write(unit=310, fmt='(A,I6)') 'Direct configuration=', it
      write(unit=310, fmt='(1X,3F12.8)') ((olap%Rt(it,i,j),j=1,3),i=1,nat)
    end do

    close(unit=310)

  end subroutine


  subroutine phQ2R(inp, olap, epc, dQ)
    ! Calculate atomic motion R(t) from normal coordinates Q_qv(t)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(overlap), intent(inout) :: olap
    type(epCoupling), intent(in) :: epc
    real(kind=q), dimension(:,:,:), intent(in) :: dQ

    integer :: iq, im, ia, ja, iax, it
    integer :: nqs, nmodes, nat, nsw
    integer :: ix, iy, iz, nx, ny, nz, ntot ! numbers of cells in supercell
    real(kind=q) :: lqvtemp, lqv, theta, qR, eiqR
    real(kind=q) :: B(3,3)
    real(kind=q), allocatable, dimension(:,:,:) :: dR, dR_frac
    complex(kind=q), allocatable, dimension(:,:,:) :: Qt, dQcum

    nqs = epc%nqpts
    nmodes = epc%nmodes
    nat = epc%natepc
    nsw = inp%NSW
    nsw = inp%NAMDTIME / inp%POTIM

    nx = 9; ny = 9; nz = 1
    ntot = nx * ny * nz

    lqvtemp = hbar * SQRT( EVTOJ / (2.0_q*AMTOKG) ) * 1.0E-5_q

    allocate(dR(nsw, ntot*nat, 3))
    allocate(dR_frac(nsw, ntot*nat, 3))
    allocate(olap%Rt(nsw, ntot*nat, 3))
    allocate(Qt(nqs, nmodes, nsw))
    allocate(dQcum(nqs, nmodes, nsw))

    dQcum(:,:,1) = dQ(:,:,1)
    do it=2,nsw
      dQcum(:,:,it) = dQcum(:,:,it-1) + dQ(:,:,it)
    end do
    Qt = dQcum + olap%PhQ(:,:,1,1:nsw) + olap%PhQ(:,:,2,1:nsw)

    dR = 0.0; dR_frac = 0.0

    do ix=1,nx
    do iy=1,ny
    do iz=1,nz

      do iq=1,nqs

        qR = epc%qpts(iq,1) * (ix - 1) + epc%qpts(iq,2) * (iy - 1) &
           + epc%qpts(iq,3) * (iz - 1)
        eiqR = EXP( -imgUnit * 2 * PI * qR )

        do im=1,nmodes
          lqv = lqvtemp / SQRT(epc%phfreq(iq,im)) ! unit of Angstrom
          do it=1,nsw
            do ia=1,nat
              ja = ny*nz*nat*(ix-1) + nz*nat*(iy-1) + nat*(iz-1) + ia
              dR(it, ja, :) = dR(it, ja, :) &
                  + eiqR * lqv * Qt(iq, im, it) / SQRT(epc%mass(ia)) * &
                    epc%phmodes(iq, im, ia, :)
            end do
          end do
        end do
      end do

    end do
    end do
    end do

    dR = dR / SQRT(olap%Np)

    call calcRecVec(epc%cellep(:3,:), B)
    do it=1,nsw
      ! calc fractional coordinates from Cartisian coord.
      dR_frac(it,:,:) = matmul(dR(it,:,:), B)
    end do

    do ix=1,nx
    do iy=1,ny
    do iz=1,nz
      do ia=1,nat
        ja = ny*nz*nat*(ix-1) + nz*nat*(iy-1) + nat*(iz-1) + ia
        do it=1,nsw
          olap%Rt(it, ja, 1) = epc%cellep(ia+3, 1) + ix - 1 + dR_frac(it, ja, 1)
          olap%Rt(it, ja, 2) = epc%cellep(ia+3, 2) + iy - 1 + dR_frac(it, ja, 2)
          olap%Rt(it, ja, 3) = epc%cellep(ia+3, 3) + iz - 1 + dR_frac(it, ja, 3)
        end do
      end do
    end do
    end do
    end do

    olap%Rt(:,:,1) = olap%Rt(:,:,1) / nx
    olap%Rt(:,:,2) = olap%Rt(:,:,2) / ny
    olap%Rt(:,:,3) = olap%Rt(:,:,3) / nz

  end subroutine


  subroutine calcRecVec(A,B)
    ! Calculate reciprocal lattice vectors.
    implicit none

    real(kind=q), intent(in) :: A(3,3)
    real(kind=q), intent(inout) :: B(3,3)

    real(kind=q) :: a1(3), a2(3), a3(3)
    real(kind=q) :: a12(3), a23(3), a31(3) ! cross product
    real(kind=q) :: V

    a1 = A(1,:)
    a2 = A(2,:)
    a3 = A(3,:)

    call cross(a1, a2, a12)
    call cross(a2, a3, a23)
    call cross(a3, a1, a31)

    V = DOT_PRODUCT(a1, a23)

    B(1,:) = a23 / V
    B(2,:) = a31 / V
    B(3,:) = a12 / V

  end subroutine


  subroutine cross(a, b, c)
    ! Calculate cross product of 3D vector a and b.
    implicit none
    real(kind=q), intent(in) :: a(3)
    real(kind=q), intent(in) :: b(3)
    real(kind=q), intent(inout) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

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
    real(kind=q) :: lqvtemp, lqv, theta, displ(3), Np
    integer :: iq, im, ia, ja, iax, it
    integer :: nqs, nmodes, nat, nsw
    integer :: ik, jk, nks, ib, jb, ibas, jbas, nb

    call readDISPL(inp, epc)
    call cellPROJ(epc)

    nks = epc%nkpts
    nb = inp%NBANDS
    nqs = epc%nqpts
    nmodes = epc%nmodes
    nsw = inp%NSW
    nat = epc%natmd

    allocate(eiqR(nqs,nat))

    do iq=1,nqs
      do ia=1,nat
        theta = 2 * PI * DOT_PRODUCT( epc%qpts(iq,:), epc%Rp(ia,:) )
        eiqR(iq, ia) = EXP(-imgUnit*theta)
      end do
    end do

    lqvtemp = hbar * SQRT( EVTOJ / (2.0_q*AMTOKG) ) * 1.0E-5_q

    do iq=1,nqs
      do im=1,nmodes
        ! unit of Angstrom
        lqv = lqvtemp / SQRT(epc%phfreq(iq,im))
        do it=1,nsw-1
          temp = cero
          do ia=1,nat
            ja = epc%atmap(ia)
            temp = temp + eiqR(iq,ia) * SQRT(epc%mass(ja)) * &
                   DOT_PRODUCT( CONJG(epc%phmodes(iq, im, ja, :)), &
                   epc%displ(it, ia, :) )
          end do
          olap%PhQ(iq,im,:,it) = temp / lqv / 2.0_q
        end do
      end do
    end do

    ! Np is the number of unit cells in MD cell.
    Np = epc%natmd / epc%natepc
    olap%PhQ = olap%PhQ / SQRT(Np)

  end subroutine


  subroutine savePhQ(olap)
    implicit none

    type(overlap), intent(in) :: olap

    integer :: it, im, iq
    integer :: nsw, nmodes, nqs
    integer :: order, ntemp, steps, ierr
    character(len=48) :: mtag
    real(kind=q) :: potim

    nsw = olap%TSTEPS
    nmodes = olap%NMODES
    nqs = olap%NQ
    potim = olap%dt

    if (nsw<1000) then
      steps = 1
    else
      order = int( LOG10(REAL(nsw)) )
      ntemp = int( nsw / (10**order) )
      if (ntemp < 2) then
        steps = 1 * 10**(order-3)
      else if (ntemp < 5) then
        steps = 2 * 10**(order-3)
      else
        steps = 5 * 10**(order-3)
      end if
    end if

    do im=1, nmodes

      write(mtag, *) im

      open(unit=36, file='PHQ.' // trim(adjustl(mtag)), &
           status='unknown', action='write', iostat=ierr)
      if (ierr /= 0) then
        write(*,*) "PHQ file I/O error!"
        stop
      end if

      do it=1, nsw-1, steps
        write(unit=36, fmt='(*(G20.10))') &
          it*potim, (olap%PhQ(iq, im, 1, it), iq=1,nqs)
      end do

      close(36)

    end do

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

    epc%kkqmap = -1

    do ik=1,epc%nkpts
      do jk=1,epc%nkpts
        do iq=1,epc%nqpts
          dkq = epc%kpts(jk,:) - epc%kpts(ik,:) - epc%qpts(iq,:)
          do iax=1,3
            dkq(iax) = ABS(dkq(iax)-NINT(dkq(iax)))
          end do
          if (SUM(dkq)<norm) then
              epc%kkqmap(jk,ik) = iq
              ! kpts[jk] = kpts[ik] + qpts[iq]
              exit
          end if
        end do
      end do
    end do

  end subroutine



  subroutine kqMatchSec(inp, epc, olap_sec)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(overlap), intent(inout) :: olap_sec
    type(epCoupling), intent(inout) :: epc

    real(kind=q) :: norm
    real(kind=q) :: dkq(3), dq(3)
    integer :: ik, jk, iq, jq, iax
    integer :: ibas, jbas, nbas

    norm = 0.001
    ! If k1-k2 < norm, recognize k1 and k2 as same k point.
    ! So, number of kx, ky, kz or qx, qy, qz must not supass 1/norm = 1000

    nbas = inp%NBASIS
    olap_sec%kkqmap = -1

    do ibas=1,nbas
      ik = inp%BASLIST(ibas,1)
      ! if ( ibas>1 .and. ik==inp%BASLIST(ibas-1,1)) then
      !   olap_sec%kkqmap(ibas,:) = olap_sec%kkqmap(ibas-1,:)
      !   cycle
      ! end if
      do jbas=1,nbas
        jk = inp%BASLIST(jbas,1)
        ! if ( jbas>1 .and. jk==inp%BASLIST(jbas-1,1)) then
        !   olap_sec%kkqmap(ibas,jbas) = olap_sec%kkqmap(ibas,jbas-1)
        !   cycle
        ! end if
        do iq=1,epc%nqpts
          dkq = epc%kpts(jk,:) - epc%kpts(ik,:) - epc%qpts(iq,:)
          do iax=1,3
            dkq(iax) = ABS(dkq(iax)-NINT(dkq(iax)))
          end do
          if (SUM(dkq)<norm) then
            olap_sec%kkqmap(jbas,ibas) = iq
            ! kpts[jk] = kpts[ik] + qpts[iq]
            exit
          end if
        end do
      end do
    end do

    allocate(epc%kkqmap(nbas, nbas))
    epc%kkqmap = olap_sec%kkqmap

  end subroutine


  subroutine calcPhQ(olap, inp, epc)
    implicit none

    type(overlap), intent(inout) :: olap
    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: iq, nq, im, nmodes, it, nsw
    real(kind=q) :: kbT, phn, phQtemp
    complex(kind=q) :: iwdt, iwtemp

    nq = olap%NQ
    nmodes = olap%NMODES
    nsw = olap%TSTEPS

    kbT = inp%TEMP * BOLKEV
    iwtemp = imgUnit * inp%POTIM / hbar

    do iq=1, nq
      do im=1, nmodes
        if (epc%phfreq(iq,im)<inp%PHCUT) cycle
        phn = 1.0 / ( exp(epc%phfreq(iq,im)/kbT) - 1.0 )
        ! Here epc%PhQ = Q_qv / l_qv = sqrt(n_qv + 0.5)
        epc%PhQ(iq, im, :) = SQRT(phn+0.5)
        iwdt = iwtemp * epc%phfreq(iq,im)
        epc%eiwdt(iq, im, 1) = exp(-iwdt)
        epc%eiwdt(iq, im, 2) = exp( iwdt)
      end do
    end do

  end subroutine


  subroutine symPhfreq(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(inout) :: epc

    integer :: ibas, jbas, nbas, iq, jq
    real(kind=q), allocatable :: freqtemp(:)

    nbas = inp%NBASIS
    allocate(freqtemp(epc%nmodes))

    if (nbas<2) return

    do ibas = 1, nbas-1
      do jbas = ibas+1, nbas
        iq = epc%kkqmap(ibas, jbas)
        if (iq<0) cycle
        jq = epc%kkqmap(jbas, ibas)
        if (jq<0) cycle
        freqtemp = ( epc%phfreq(iq,:) + epc%phfreq(jq,:) ) / 2.0
        epc%phfreq(iq,:) = freqtemp
        epc%phfreq(jq,:) = freqtemp
      end do
    end do

  end subroutine


  subroutine calcEPC(olap, epc, inp)
  ! EPC calculations for large basis set.
    implicit none

    type(overlap), intent(inout) :: olap
    type(epCoupling), intent(inout) :: epc
    type(namdInfo), intent(in) :: inp

    integer :: irank, ierr
    integer :: ist, iend
    integer :: nb, nb_p, nm
    integer :: ib, jb, im, iq
    real(kind=q) :: dE, dE1, dE2
    real(kind=q) :: kbT, sigma, T0
    complex(kind=q) :: idwt

    nb_p = inp%NBASIS_P
    nb  = olap%NBANDS
    nm  = olap%NMODES
    ! kbT = inp%TEMP * BOLKEV
    ! idwt = imgUnit / kbT * TPI
    ! T0 = TPI / kbT / 2.0
    ! T = TPI * hbar / kbT in unit of fs.
    ! dw * T/2.0 = dE * T/2.0 / hbar
    ! T0 = T / hbar / 2.0
    if (inp%SIGMA==0) then
      sigma = inp%TEMP * BOLKEV
    else
      sigma = inp%SIGMA
    end if

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
    ist = inp%ISTS(irank+1)
    iend = inp%IENDS(irank+1)

    do ib=ist, iend
     do jb=1,nb

       dE = epc%eig(ib) - epc%eig(jb)
       iq = olap%kkqmap(ib,jb)
       if (iq<0) cycle

       do im=1,nm
         dE1 = dE - epc%phfreq(iq,im) - 1.0E-8_q
         epc%epcec(ib-ist+1,jb,im,1) &
           = ( ABS(epc%gij(ib-ist+1,jb,im)) * ABS(epc%PhQ(iq, im, 1)) ) ** 2 &
             * exp(-0.5 * (dE1/sigma)**2)
         dE2 = dE + epc%phfreq(iq,im) + 1.0E-8_q
         epc%epcec(ib-ist+1,jb,im,2) &
           = ( ABS(epc%gij(ib-ist+1,jb,im)) * ABS(epc%PhQ(iq, im, 1)) ) ** 2 &
             * exp(-0.5 * (dE2/sigma)**2)
       end do ! im loop

     end do
    end do
    epc%gij = epc%gij / SQRT(olap%Np)
    epc%epcec = epc%epcec/ olap%Np * SQRT(TPI) / sigma

  end subroutine


  subroutine epcToOlapSec(inp, epc, olap_sec)
    implicit none
    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(in) :: epc
    type(overlap), intent(inout) :: olap_sec

    integer :: ik, jk, ib, jb, nb, ibas, jbas, nbas

    nbas = inp%NBASIS

    nb = inp%NBANDS
    do ibas=1, nbas
      ik = inp%BASLIST(ibas,1)
      ib = inp%BASLIST(ibas,2)
      olap_sec%Eig(ibas, :) = epc%elen(ik, ib)
    end do

    olap_sec%Phfreq = epc%phfreq

  end subroutine


  subroutine copyToSec(olap, olap_sec, inp)
    implicit none
    type(namdInfo), intent(inout) :: inp
    type(overlap), intent(inout) :: olap_sec
    type(overlap), intent(in) :: olap

    integer :: ik, jk, ib, jb, nb, ibas, jbas, nbas
    integer :: kmin, kmax, bmin, bmax
    real :: t

    nb = inp%NBANDS
    nbas = inp%NBASIS

    do ibas=1,nbas
      ik = inp%BASLIST(ibas,1)
      ib = inp%BASLIST(ibas,2)
      olap_sec%Eig(ibas, :) = olap%Eig((ik-1)*nb+ib, :)
      do jbas=1,nbas
        jk = inp%BASLIST(jbas,1)
        jb = inp%BASLIST(jbas,2)

        olap_sec%kkqmap(ibas, jbas) = &
        olap%kkqmap( (ik-1)*nb+ib, (jk-1)*nb+jb )

        olap_sec%gij(ibas, jbas, :) = &
        olap%gij( (ik-1)*nb+ib, (jk-1)*nb+jb, : )

      end do
    end do

    olap_sec%Phfreq = olap%Phfreq
    olap_sec%PhQ = olap%PhQ

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
        write(*,*)
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
      ! write(unit=30, rec=irec) (olap%Eig(i,1), i=1,nb), &
      !     ((olap%Phfreq(i,j,im), i=1,nb), j=1,nb)
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
      write(*,*)
      stop
    end if

    read(unit=31,rec=1) recordL, rnbands, rnsw, rdt, rctype, rnp, rnmodes
    ! write(*,*) recordL, rnbands, rnsw, rdt

    if (olap%NBANDS /= NINT(rnbands) .or. olap%TSTEPS /= NINT(rnsw) .or. &
       olap%COUPTYPE /= NINT(rctype) .or. olap%NMODES /= NINT(rnmodes)) then
      ! write(*,*) olap%NBANDS, NINT(rnbands), olap%TSTEPS, NINT(rnsw)
      write(*,*) "The EPCAR seems to be wrong..."
      write(*,*)
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
      ! read(unit=31, rec=irec) (olap%Eig(i,1), i=1,nb), &
      !     ((olap%Phfreq(i,j,im), i=1,nb), j=1,nb)
      do i=1,nb
        olap%Eig(i,:) = olap%Eig(i,1)
      end do
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

  subroutine writeEPELTXT(inp, epc)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(in) :: epc

    integer :: irank, jrank, nrank, ierr
    integer :: ibas, jbas, nbas, nbas_p

    call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nrank, ierr)

    nbas   = inp%NBASIS
    nbas_p = inp%NBASIS_P

    ! if (irank==0) then
    !   open(unit=34, file='EPECTXT', status='unknown', action='write')
    !   close(unit=34)
    ! end if
    ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    do jrank = 0, nrank-1

      if (irank==jrank) then

        if (irank==0) then
          open(unit=34, file='EPELTXT', status='unknown', action='write')
        else
          open(unit=34, file='EPELTXT', status='old', position='append', action='write')
        end if

        do ibas = 1, nbas_p
          write(unit=34, fmt='(*(f15.9))') &
               ( SUM(epc%epcec(ibas,jbas,:,:)), jbas=1,nbas )
        end do

        close(unit=34)

      end if

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    end do

  end subroutine


  subroutine writeEPPHTXT(inp, epc, olap)
    implicit none

    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(in) :: epc
    type(overlap), intent(in) :: olap

    integer :: irank, ierr
    integer :: ibas, jbas, nbas
    integer :: iq, nqs, im, nmodes
    integer :: ist, iend
    real(kind=q), allocatable :: epcph(:,:), epcph_tot(:,:)

    call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)

    nbas   = inp%NBASIS
    nqs    = epc%nqpts
    nmodes = epc%nmodes

    ist = inp%ISTS(irank+1)
    iend = inp%IENDS(irank+1)

    allocate(epcph(nqs, nmodes))
    if (irank==0) allocate(epcph_tot(nqs, nmodes))
    if (irank==0) epcph_tot = 0.0_q
    epcph = 0.0_q

    do ibas=ist, iend
      do jbas=1, nbas
        iq = olap%kkqmap(ibas,jbas)
        if (iq<0) cycle
        epcph(iq,:) = epcph(iq,:) + SUM( epc%epcec(ibas-ist+1, jbas, :, :), DIM=2 )
      end do
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do im = 1, nmodes
      call MPI_Reduce(epcph(:,im), epcph_tot(:,im), nqs, &
           MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end do

    if (irank==0) then
      open(unit=35, file='EPPHTXT', status='unknown', action='write')
      do iq=1,nqs
        write(unit=35, fmt='(*(f15.9))') (epcph_tot(iq,im), im=1,nmodes)
      end do
      close(unit=35)
    end if

  end subroutine


  subroutine writeTXT_LBS(olap)
    implicit none

    type(overlap), intent(inout) :: olap
    integer :: iq, jq, nqs, im, nmodes, ib, jb, nb, ierr, it, nsw
    real(kind=q), allocatable :: epc(:,:), epcec(:,:), epcph(:,:)

    nb = olap%NBANDS
    nqs = olap%NQ
    nmodes = olap%NMODES
    nsw = olap%TSTEPS

    allocate(epc(nb,nb), epcec(nb,nb), epcph(nqs, nmodes))
    epc = 0.0_q; epcec = 0.0_q; epcph = 0.0_q

    do ib=1,nb
      do jb=ib,nb
        iq = olap%kkqmap(ib,jb)
        if (iq<0) cycle
        do it=1,nsw-1
          epc(ib,jb) = epc(ib,jb) + ABS( SUM(olap%gij(ib,jb,:) * &
            SUM(olap%PhQ(iq,:,:,it), dim=2)) )
          epcec(ib,jb) = epcec(ib,jb) + &
            ! ABS( SUM(olap%EPcoup(ib,jb,:,:,1) * olap%PhQ(iq,:,:,it)) )
            ABS( SUM(olap%EPcoup(ib,jb,:,:,1) * (olap%PhQ(iq,:,:,it)**2)) )
          epcph(iq,:) = epcph(iq,:) + &
            ABS( SUM(olap%EPcoup(ib,jb,:,:,1) * (olap%PhQ(iq,:,:,it)**2), DIM=2) )
        end do
        epc(jb,ib) = epc(ib,jb)
        epcec(jb,ib) = epcec(ib,jb)
        jq = olap%kkqmap(jb,ib)
        if (jq<0) cycle
        epcph(jq,:) = epcph(iq,:)
      end do
    end do
    epc = epc / (nsw-1)
    epcec = epcec / (nsw-1)
    epcph = epcph / (nsw-1)

    open(unit=32, file='EIGTXT', status='unknown', action='write')
    open(unit=33, file='EPTXT', status='unknown', action='write')
    open(unit=34, file='EPECTXT', status='unknown', action='write')
    open(unit=35, file='EPPHTXT', status='unknown', action='write')


    write(unit=32, fmt='(*(f12.6))') (olap%Eig(ib,1), ib=1,nb)
    write(unit=33, fmt='(*(f15.9))') (( epc(ib,jb), jb=1,nb ), ib=1,nb)
    write(unit=34, fmt='(*(f15.9))') (( epcec(ib,jb), jb=1,nb ), ib=1,nb)
    do iq=1,nqs
      write(unit=35, fmt='(*(f15.9))') (epcph(iq,im), im=1,nmodes)
    end do

    close(unit=32)
    close(unit=33)
    close(unit=34)
    close(unit=35)

  end subroutine


  subroutine writeTXT_LBS_Mode(olap)
    implicit none

    type(overlap), intent(inout) :: olap
    integer :: im, nmodes, ib, jb, nb, iq, ierr, it, nsw
    real(kind=q), allocatable :: epc(:,:), epcec(:,:)
    character(len=48) :: mtag, filename

    nb = olap%NBANDS
    nmodes = olap%NMODES
    nsw = olap%TSTEPS

    allocate(epc(nb,nb), epcec(nb,nb))
    epc = 0.0_q; epcec = 0.0_q

    open(unit=32, file='EIGTXT', status='unknown', action='write')
    write(unit=32, fmt='(*(f12.6))') (olap%Eig(ib,1), ib=1,nb)
    close(unit=32)

    do im=1,nmodes
      do ib=1,nb
        do jb=ib,nb
          iq = olap%kkqmap(ib,jb)
          do it=1,nsw-1
            epc(ib,jb) = epc(ib,jb) + ABS( olap%gij(ib,jb,im) * &
              SUM(olap%PhQ(iq,im,:,it)) )
            epcec(ib,jb) = epcec(ib,jb) + &
              ABS( SUM(olap%EPcoup(ib,jb,im,:,1) * olap%PhQ(iq,im,:,it)) )
          end do
          epc(jb,ib) = epc(ib,jb)
          epcec(jb,ib) = epcec(ib,jb)
        end do
      end do
      epc = epc / (nsw-1)
      epcec = epcec / (nsw-1)

      write(mtag, *) im

      filename = 'EPTXT.' // trim(adjustl(mtag))
      open(unit=33, file=filename, status='unknown', action='write')
      filename = 'EPECTXT.' // trim(adjustl(mtag))
      open(unit=34, file=filename, status='unknown', action='write')

      write(unit=33, fmt='(*(f15.9))') (( epc(ib,jb), jb=1,nb ), ib=1,nb)
      write(unit=34, fmt='(*(f15.9))') (( epcec(ib,jb), jb=1,nb ), ib=1,nb)

      close(unit=33)
      close(unit=34)

    end do

  end subroutine

  subroutine writeKKQMap(olap)
    implicit none

    type(overlap), intent(in) :: olap
    integer :: ib, jb, nbas

    nbas = olap%NBANDS

    open(unit=35, file='KKQMAP', status='unknown', action='write')

    do ib=1,nbas
      write(unit=35, fmt='(*(I8))') (olap%kkqmap(ib, jb), jb=1,nbas)
    end do

    close(unit=35)

  end subroutine


  subroutine writeTXT(olap)
    implicit none

    type(overlap), intent(in) :: olap
    integer :: im, nmodes, ib, jb, nb, it, nsw, ierr

    nb = olap%NBANDS
    nsw = olap%TSTEPS
    nmodes = olap%NMODES

    open(unit=32, file='EIGTXT',  status='unknown', action='write')
    open(unit=33, file='EPTXT',   status='unknown', action='write')
    open(unit=34, file='EPECTXT', status='unknown', action='write')

    do it=1,nsw-1
      write(unit=32, fmt='(*(f12.6))') (olap%Eig(ib,it), ib=1,nb)
      write(unit=33, fmt='(*(f15.9))') &
        ((olap%Dij(ib,jb,it), jb=1,nb), ib=1,nb)
      write(unit=34, fmt='(*(f15.9))') &
        (( SUM(olap%EPcoup(ib,jb,:,:,it)), jb=1,nb ), ib=1,nb)
    end do

    close(unit=32)
    close(unit=33)
    close(unit=34)

  end subroutine


  subroutine readEPTXTs(olap)
    implicit none

    type(overlap), intent(inout) :: olap
    integer :: im, nmodes, ib, jb, nb, it, nsw, ierr

    open(unit=32, file='EIGTXT', status='unknown', &
         action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "EIGTXT does NOT exist!"
      write(*,*)
      stop
    end if
    open(unit=33, file='EPTXT', status='unknown', &
         action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "EPTXT does NOT exist!"
      write(*,*)
      stop
    end if
    open(unit=34, file='EPECTXT', status='unknown', &
         action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "EPECTXT does NOT exist!"
      write(*,*)
      stop
    end if

    nb = olap%NBANDS
    nsw = olap%TSTEPS
    do it=1,nsw-1
      read(unit=32, fmt=*) (olap%Eig(ib,it), ib=1,nb)
      read(unit=33, fmt='(*(f15.9))') &
        ((olap%Dij(ib,jb,it), jb=1,nb), ib=1,nb)
      read(unit=34, fmt='(*(f15.9))') &
        (( olap%EPcoup(ib,jb,1,1,it), jb=1,nb ), ib=1,nb)
    end do

    close(unit=32)
    close(unit=33)
    close(unit=34)

  end subroutine


  subroutine mpi_split_bas(inp)

    implicit none
    type(namdInfo), intent(inout) :: inp

    integer :: irank, nrank, ierr
    integer :: ik, ib, ibas, nbas, nbas_p
    integer, allocatable :: ists(:), iends(:)
    integer :: ist, iend ! local start and end index

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nrank, ierr)

    nbas = inp%NBASIS
    allocate(ists(nrank), iends(nrank))
    allocate(inp%ISTS(nrank), inp%IENDS(nrank))
    CALL mpi_split_procs(nbas, nrank, ists, iends)
    inp%ISTS = ists; inp%IENDS = iends

    ist = ists(irank+1)
    iend = iends(irank+1)
    nbas_p = iend - ist + 1
    inp%NBASIS_P = nbas_p
    allocate(inp%BASLIST_P(nbas_p, 3))

    inp%BASSEL_P = -1
    inp%BASLIST_P = inp%BASLIST(ist:iend,:)
    do ibas = ist, iend
      ik = inp%BASLIST(ibas, 1)
      ib = inp%BASLIST(ibas, 2)
      inp%BASSEL_P(ik, ib) = ibas
    end do

  end subroutine


  subroutine mpi_split_procs(num, nrank, ists, iends)

    implicit none
    integer, intent(in) :: num, nrank
    integer, intent(out) :: ists(nrank), iends(nrank)

    integer :: base, rest
    integer :: i

    base = num / nrank
    rest = MOD(num, nrank)

    do i = 1, nrank
      if(i <= rest) then
         ists(i) = (i - 1) * (base + 1) + 1
         iends(i) = ists(i) + base
      else
         ists(i) = (i - 1) * base + merge(rest, 0, base>0) + 1
         iends(i) = ists(i) + base - 1
      end if
    enddo

  end subroutine mpi_split_procs


  subroutine testMat(mat, N)
    implicit none

    complex(kind=q), intent(in) :: mat(:,:)
    integer, intent(in) :: N

    integer :: ii, jj
    real(kind=q) :: temp

    temp = 0.0_q
    do ii = 1, N
      do jj = ii, N
        ! print '(*(I6))', ii, jj
        ! print '(*(F12.6))', mat(ii,jj), mat(jj,ii)
        ! print '(*(F12.6))', ABS(mat(ii,jj))-ABS(mat(jj,ii))
        temp = temp + ABS(mat(ii,jj))-ABS(mat(jj,ii))
      end do
    end do
    print '(F12.6)', temp
    ! print *
  end subroutine


end module epcoup
