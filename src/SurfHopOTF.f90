module shotf
  use prec
  use fileio
  use hamil
  use couplings
  use shop
  implicit none

  contains

  subroutine runSHotf(ks, inp, olap, epc)
  ! Surface Hopping Simulation without CPA (On-The-Fly type)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    type(overlap), intent(in) :: olap
    type(epCoupling), intent(in) :: epc

    integer :: cstat, nstat
    integer, allocatable :: cstat_all(:), occb(:)
    integer :: i, j, ibas, nbas, tion, Nt, iq
    complex(kind=q), allocatable :: dQ(:,:,:), Qtemp(:,:,:)
    logical :: lhop

    ks%sh_pops = 0
    ks%sh_prop = 0
    nbas = ks%ndim
    Nt = inp%NAMDTIME / inp%POTIM

    allocate(cstat_all(inp%NINIBS))
    allocate(occb(inp%NBASIS)) ! tag of basis occupied
    allocate(dQ(olap%NQ, olap%NMODES, 2))
    allocate(Qtemp(olap%NQ, olap%NMODES, 2))

    ! initialize the random seed for ramdom number production
    call init_random_seed()
    ! call calcBftot(ks, inp)
    ks%sh_Bfactor = 1.0
    lhop = .TRUE.

    do i=1, inp%NTRAJ

      ks%psi_c = cero
      ks%psi_p = cero
      ks%psi_n = cero

      occb = 0
      do j=1, inp%NINIBS
        ibas = inp%BASSEL(inp%INIKPT(j), inp%INIBAND(j))
        cstat_all(j) = ibas
        occb(ibas) = 1
        ks%psi_c(ibas) = uno
      end do

      ks%psi_c = SQRT(ks%psi_c / REAL(inp%NINIBS))

      dQ = cero

      do tion=1,Nt

        Qtemp = ks%PhQ(:,:,:,tion)
        ks%PhQ(:,:,:,tion) = ks%PhQ(:,:,:,tion) + dQ

        if (i==1) then
        ks%pop_a(:,tion) = CONJG(ks%psi_c) * ks%psi_c
        ks%norm(tion) = SUM(ks%pop_a(:,tion))
        ks%psi_a(:,tion) = ks%psi_c
        call CProp(tion, ks, inp, olap, epc)
        end if

        do j=1, inp%NINIBS

          cstat = cstat_all(j)
          call calcprop_EPC(tion, cstat, ks, inp, olap, epc)
        ! call calcprop_OTF(tion, cstat, ks, inp, olap)
          call whichToHop(cstat, nstat, ks)

          iq = olap%kkqmap(cstat, nstat)
          if (occb(nstat)>0 .OR. iq<0) then
            ks%sh_pops(cstat, tion) = ks%sh_pops(cstat, tion) + 1
          else
          ! call calcPHprop(tion, cstat, nstat, ks, olap)
            call calcDQ(tion, cstat, nstat, dQ, ks, olap, lhop)
            if (lhop) then
              occb(cstat) = 0; occb(nstat) = 1
              ks%sh_pops(nstat, tion) = ks%sh_pops(nstat, tion) + 1
              cstat_all(j) = nstat
              ks%ph_pops(iq, :, tion) &
                = ks%ph_pops(iq, :, tion) &
                + SUM(ks%ph_prop(cstat, nstat, :, :), dim=2)
            else
              ks%sh_pops(cstat, tion) = ks%sh_pops(cstat, tion) + 1
            end if
          end if

        end do

        ks%PhQ(:,:,:,tion) = Qtemp

      end do

    end do

    ks%sh_pops = ks%sh_pops / inp%NTRAJ

  end subroutine


  subroutine calcprop_OTF(tion, cstat, ks, inp, olap)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion
    integer, intent(in) :: cstat
    type(overlap), intent(in) :: olap

    integer :: iq, i
    real(kind=q) :: Akk
    complex(kind=q), allocatable :: epcoup(:)

    allocate(epcoup(ks%ndim))
    do i=1,ks%ndim
      iq = olap%kkqmap(cstat, i)
      epcoup(i) = SUM( olap%gij(cstat,i,:) * &
                       SUM(ks%PhQ(iq,:,:,tion), dim=2) )
    end do


    Akk = CONJG(ks%psi_a(cstat, tion)) * ks%psi_a(cstat, tion)
    ks%Bkm = -2. / hbar * AIMAG( CONJG(ks%psi_a(cstat, tion)) * &
             ks%psi_a(:, tion) * epcoup(:) )

    ks%sh_prop(cstat,:) = ks%Bkm / Akk * inp%POTIM
    forall (i=1:ks%ndim, ks%sh_prop(cstat,i) < 0) ks%sh_prop(cstat,i) = 0

  end subroutine


  subroutine calcPHprop(tion, cstat, nstat, ks, olap)
    implicit none
    integer, intent(in) :: tion, cstat, nstat
    type(TDKS), intent(inout)  :: ks
    type(overlap), intent(in) :: olap

    integer :: iq, im, nmodes
    real(kind=q), allocatable :: eptemp(:,:)
    real(kind=q) :: norm

    nmodes = olap%NMODES
    iq = olap%kkqmap(cstat, nstat)
    allocate(eptemp(nmodes, 2))

    do im=1,nmodes
      eptemp(im,:) = ABS(olap%gij(cstat,nstat,im) &
                         * ks%PhQ(iq,im,:,tion)) ** 2
      ! eptemp(im,:) = ABS(olap%EPcoup(cstat,nstat,im,:,1) &
      !                    * (ks%PhQ(iq,im,:,tion) ** 2))
    end do
    norm = SUM(eptemp)
    ks%ph_prop(cstat,nstat,:,1) = -eptemp(:,1) / norm
    ks%ph_prop(cstat,nstat,:,2) =  eptemp(:,2) / norm

  end subroutine


  subroutine calcDQ(tion, cstat, nstat, dQ, ks, olap, lhop)
    implicit none
    integer, intent(in) :: tion, cstat, nstat
    complex(kind=q), intent(inout) :: dQ(:,:,:)
    type(TDKS), intent(in)  :: ks
    type(overlap), intent(in) :: olap
    logical, intent(inout) :: lhop

    integer :: iq, im, nmodes, i
    real(kind=q), allocatable :: phn(:,:)

    lhop = .TRUE.
    nmodes = olap%NMODES
    iq = olap%kkqmap(cstat, nstat)

    allocate(phn(nmodes,2))
    phn = ABS(ks%PhQ(iq, :, :, tion)) ** 2

    mloop: do im=1,nmodes
      do i=1,2
        if (phn(im, i) + ks%ph_prop(cstat, nstat, im, i) < 0) then
          lhop = .FALSE.
          exit mloop
        end if
      end do
    end do mloop

    if (lhop) then
      do im=1,nmodes
        do i=1,2
          if (ks%ph_prop(cstat, nstat, im, i)>0) then
            dQ(iq, im, i) = dQ(iq, im, i) + ks%PhQ(iq, im, i, tion) &
                          * SQRT(ks%ph_prop(cstat, nstat, im, i)/phn(im, i))
          else
            dQ(iq, im, i) = dQ(iq, im, i) - ks%PhQ(iq, im, i, tion) &
                          * SQRT(-ks%ph_prop(cstat, nstat, im, i)/phn(im, i))
          end if
        end do
      end do
    end if

  end subroutine


  subroutine CProp(tion, ks, inp, olap, epc)
    implicit none
    integer, intent(in) :: tion
    type(TDKS), intent(inout)  :: ks
    type(namdInfo), intent(in) :: inp
    type(overlap), intent(in) :: olap
    type(epCoupling), intent(in) :: epc

    integer :: tele
    real(kind=q) :: edt

    edt = inp%POTIM / inp%NELM

    do tele = 1, inp%NELM-1
      ! construct hamiltonian matrix
      if (inp%LEPC) then
        call make_hamil_EPC(tion, tele, ks, inp, olap, epc)
      else
        call make_hamil(tion, tele, ks, inp)
      end if
      ! apply hamiltonian to state vector
      ! call hamil_act(ks)
      ks%hpsi = matmul(ks%ham_c, ks%psi_c)
      if (tion == 1 .AND. tele == 1) then
        ! write(*,*) ((ks%ham_c(i,j), j=1, ks%ndim), i=1, ks%ndim)
        ! This is the very first step of the time propagation
        ! use first order difference
        ! [c,n,p] meas current, next, previous respectively
        ks%psi_n = ks%psi_c - imgUnit * edt * ks%hpsi / hbar
      else
        ! use second order difference
        ks%psi_n = ks%psi_p - 2 * imgUnit * edt * ks%hpsi / hbar
      end if
      ks%psi_p = ks%psi_c
      ks%psi_c = ks%psi_n
    end do
 
  end subroutine


end module
