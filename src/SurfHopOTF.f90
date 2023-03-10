module shotf
  use prec
  use fileio
  use hamil
  use couplings
  use shop
  implicit none

  contains

  subroutine runSHotf(ks, inp, olap)
  ! Surface Hopping Simulation without CPA (On-The-Fly type)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    type(overlap), intent(in) :: olap

    integer :: i, j, ibas, nbas, tion, Nt, iq
    integer, allocatable :: cstat_all(:), occb(:)
    integer :: cstat, nstat

    ks%sh_pops = 0
    ks%sh_prop = 0
    nbas = ks%ndim
    Nt = inp%NAMDTIME / inp%POTIM

    allocate(cstat_all(inp%NINIBS))
    allocate(occb(inp%NBASIS))
    ! tag of basis occupied

    ! initialize the random seed for ramdom number production
    call init_random_seed()
    call calcBftot(ks, inp)

    do i=1, inp%NTRAJ

      occb = 0
      do j=1, inp%NINIBS
        ibas = inp%BASSEL(inp%INIKPT(j), inp%INIBAND(j))
        cstat_all(j) = ibas
        occb(ibas) = 1
      end do

      do tion=1,Nt

        if (i==1) then
        ks%pop_a(:,tion) = CONJG(ks%psi_c) * ks%psi_c
        ks%norm(tion) = SUM(ks%pop_a(:,tion))
        ks%psi_a(:,tion) = ks%psi_c
        call CProp(tion, ks, inp, olap)
        end if

        do j=1, inp%NINIBS

          cstat = cstat_all(j)
          call calcprop_EPC(tion, cstat, ks, inp, olap)
          call whichToHop(cstat, nstat, ks)

          if (nstat /= cstat .AND. occb(nstat)>0) cycle
          occb(cstat) = 0; occb(nstat) = 1
          ks%sh_pops(nstat, tion) = ks%sh_pops(nstat, tion) + 1
          cstat_all(j) = nstat

          if (nstat == cstat) cycle
          iq = olap%kkqmap(cstat, nstat)
          if (iq>0) then
            ks%ph_pops(iq, :, tion) &
              = ks%ph_pops(iq, :, tion) + ks%ph_prop(cstat, nstat, :)
          end if

        end do

      end do

    end do

    ks%sh_pops = ks%sh_pops / inp%NTRAJ

  end subroutine


  subroutine CProp(tion, ks, inp, olap)
    implicit none
    integer, intent(in) :: tion
    type(TDKS), intent(inout)  :: ks
    type(namdInfo), intent(in) :: inp
    type(overlap), intent(in) :: olap

    integer :: tele
    real(kind=q) :: edt

    edt = inp%POTIM / inp%NELM

    do tele = 1, inp%NELM-1
      ! construct hamiltonian matrix
      if (inp%LEPC) then
        call make_hamil_EPC(tion, tele, ks, inp, olap)
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
