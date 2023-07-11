module hamil
  use prec
  use fileio
  use couplings
  use constants
  use epcoup
  implicit none

  type TDKS
    integer :: ndim
    ! _[c,p,n] means current, previous, next
    complex(kind=q), allocatable, dimension(:) :: psi_c
    complex(kind=q), allocatable, dimension(:) :: psi_p
    complex(kind=q), allocatable, dimension(:) :: psi_n
    complex(kind=q), allocatable, dimension(:,:) :: psi_a
    ! the result of hamiltonian acting on a vector
    complex(kind=q), allocatable, dimension(:) :: hpsi
    ! population
    real(kind=q), allocatable, dimension(:,:) :: pop_a
    real(kind=q), allocatable, dimension(:) :: norm

    complex(kind=q), allocatable, dimension(:,:) :: ham_c
    complex(kind=q), allocatable, dimension(:,:) :: ham_p
    complex(kind=q), allocatable, dimension(:,:) :: ham_n

    ! KS eigenvalues
    real(kind=q), allocatable, dimension(:,:) :: eigKs
    ! Non-adiabatic couplings
    complex(kind=q), allocatable, dimension(:,:,:) :: NAcoup

    complex(kind=q), allocatable, dimension(:,:,:) :: PhQtemp
    complex(kind=q), allocatable, dimension(:,:,:,:) :: PhQ

    ! surface hopping related

    ! Bkm = REAL(CONJG(Akm) * Ckm)
    real(kind=q), allocatable, dimension(:) :: Bkm
    real(kind=q), allocatable, dimension(:,:) :: sh_pops
    real(kind=q), allocatable, dimension(:,:) :: sh_prop
    real(kind=q), allocatable, dimension(:,:,:) :: ph_pops
    real(kind=q), allocatable, dimension(:,:,:,:) :: ph_prop
    ! Blotzmann factor for SH probability scaling.
    real(kind=q), allocatable, dimension(:,:) :: sh_Bfactor

    ! whether the memory has been allocated
    logical :: LALLO = .FALSE.

  end type

  contains

  subroutine initTDKS(ks, inp, olap, epc)
    implicit none

    type(TDKS), intent(inout)  :: ks
    type(overlap), intent(in)  :: olap
    type(namdInfo), intent(in) :: inp
    type(epCoupling), intent(in)  :: epc

    real(kind=q) :: norm
    real(kind=q), allocatable, dimension(:,:) :: eptemp
    integer :: i, j, ib, iq, im, t, nsteps, N, Nt, nmodes, nqs
    integer :: initstep
    integer :: irank, ierr
    integer :: ist, iend, N_p

    ! memory allocation

    N = inp%NBASIS
    N_p = inp%NBASIS_P
    ks%ndim = inp%NBASIS
    Nt = inp%NAMDTIME / inp%POTIM

    if (.NOT. ks%LALLO) then
      allocate(ks%psi_c(N))
      allocate(ks%psi_p(N))
      allocate(ks%psi_n(N))
      allocate(ks%hpsi(N_p))
      ! allocate(ks%psi_a(N, Nt))
      ! allocate(ks%pop_a(N, Nt))
      allocate(ks%norm(Nt))

      allocate(ks%ham_c(N_p,N))
      allocate(ks%ham_p(N_p,N))
      allocate(ks%ham_n(N_p,N))

      allocate(ks%eigKs(N, Nt))


      allocate(ks%sh_pops(N, 1))
      allocate(ks%sh_prop(N,N))
      allocate(ks%sh_Bfactor(N,N))
      allocate(ks%Bkm(N))

      if (.NOT. inp%LEPC) then

        allocate(ks%NAcoup(N,N, Nt))

      else

        nqs = olap%NQ; nmodes = olap%NMODES
        allocate(ks%PhQtemp(nqs, nmodes, 2))
        ! allocate(ks%PhQ(nqs, nmodes, 2, Nt))
        allocate(ks%ph_pops(nqs, nmodes, Nt))
        allocate(ks%ph_prop(N, N, nmodes, 2))

        ks%ph_pops = 0.0; ks%ph_prop = 0.0
        allocate(eptemp(nmodes, 2))

        CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
        ist = inp%ISTS(irank+1)
        iend = inp%IENDS(irank+1)

        do i=ist, iend
          do j=1,N
            iq = olap%kkqmap(i,j)
            eptemp = epc%epcec(i-ist+1,j,:,:)
            norm = SUM(eptemp)
            if (norm>0) then
            ! ks%ph_prop(i,j,:) = (eptemp(:,2) - eptemp(:,1)) / norm
              ks%ph_prop(i,j,:,1) =  eptemp(:,1) / norm
              ks%ph_prop(i,j,:,2) = -eptemp(:,2) / norm
            end if
          end do
        end do

      end if

      ks%LALLO = .TRUE.

    end if

    ks%psi_c = cero
    ks%psi_p = cero
    ks%psi_n = cero

    if (inp%LCPROP) then
      do i=1, inp%NINIBS
        ib = inp%BASSEL(inp%INIKPT(i), inp%INIBAND(i))
        ks%psi_c(ib) = uno
      end do
      ks%psi_c = SQRT(ks%psi_c / REAL(inp%NINIBS))
    else
      do i=1, inp%NSAMPLE
        do j=1, inp%NINIBS
          ib = inp%BASSEL(inp%INIKPT_A(i,j), inp%INIBAND_A(i,j))
          ks%psi_c(ib) = ks%psi_c(ib) + uno
        end do
      end do
      ks%psi_c = SQRT(ks%psi_c / REAL(inp%NSAMPLE * inp%NINIBS))
    end if

    nsteps = inp%NSW - 1
    initstep = inp%NAMDTINI / inp%POTIM - 2

    if (inp%LEPC) then

      do t=1, Nt
        i = MOD(initstep+t, nsteps) + 1
        ks%eigKs(:,t) = olap%Eig(:,i)
      end do
      ks%PhQtemp = epc%PhQ * (epc%eiwdt ** (initstep + 1))

    else

      do t=1, Nt
        ! If time step > NSW-1, use Eig & couplings
        ! from initial time repeatedly.
        i = MOD(initstep+t, nsteps) + 1
        ! We don't need all the information, only a section of it
        ks%eigKs(:,t) = olap%Eig(:, i)
        ! Divide by 2 * POTIM here,
        ! because we didn't do this in the calculation of couplings
        ks%NAcoup(:,:,t) = olap%Dij(:,:, i) / (2*inp%POTIM)
      end do

    end if

  end subroutine


  ! constructing the hamiltonian
  subroutine make_hamil(TION, TELE, ks, inp)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: TION, TELE

    integer :: i

    ! the hamiltonian contains two parts, which are obtained by
    ! interpolation method between two ionic tims step

    ! The non-adiabatic coupling part
    ks%ham_c(:,:) = ks%NAcoup(:,:,TION) + &
      (ks%NAcoup(:,:,TION+1) - ks%NAcoup(:,:,TION)) * TELE / inp%NELM

    ! multiply by -i * hbar
    if (.not. inp%LEPC) ks%ham_c = -imgUnit * hbar * ks%ham_c

    ! the energy eigenvalue part
    do i=1, ks%ndim
      ks%ham_c(i,i) = ks%ham_c(i,i) + ks%eigKs(i,TION) + &
        (ks%eigKs(i,TION+1) - ks%eigKs(i,TION)) * TELE / inp%NELM
    end do
  end subroutine


  subroutine make_hamil_EPC(tion, tele, ks, inp, olap, epc)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion, tele
    type(overlap), intent(in) :: olap
    type(epCoupling), intent(in) :: epc

    integer :: ib, jb, iq
    complex(kind=q), allocatable :: tempQ(:,:,:)

    if (tion==1 .AND. tele==1) call calc_hamil_EPC(tion, ks, inp, olap, epc)
    if (tele==1) call calc_hamil_EPC(tion+1, ks, inp, olap, epc)

    ks%ham_c = ks%ham_p + (ks%ham_n - ks%ham_p) * TELE / inp%NELM

  end subroutine


  subroutine calc_hamil_EPC(tion, ks, inp, olap, epc)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: TION
    type(overlap), intent(in) :: olap
    type(epCoupling), intent(in) :: epc

    integer :: ib, jb, iq
    integer :: irank, ierr
    integer :: ist, iend

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
    ist = inp%ISTS(irank+1)
    iend = inp%IENDS(irank+1)

    ! ks%PhQtemp = ks%PhQtemp * epc%eiwdt

    ks%ham_p = ks%ham_n
    do ib=ist, iend
      ! do jb=ib,ks%ndim
      do jb=1,ks%ndim
        iq = olap%kkqmap(ib,jb)
        if (iq<0) cycle
        ks%ham_n(ib-ist+1,jb) = SUM( epc%gij(ib-ist+1,jb,:) * &
            SUM(ks%PhQtemp(iq,:,:), dim=2) )
        ! ks%ham_n(jb,ib) = CONJG(ks%ham_n(ib,jb))
      end do
      ks%ham_n(ib-ist+1,ib) = ks%ham_n(ib-ist+1,ib) + ks%eigKs(ib,1)
    end do

  end subroutine

  ! Acting the hamiltonian on the state vector
  subroutine hamil_act(ks)
    implicit none
    type(TDKS), intent(inout) :: ks
    integer :: i, j, N
    complex(kind=q) :: tmp

    N = ks%ndim
    do i=1, N
      tmp = cero
      do j=1, N
        tmp = tmp + ks%ham_c(i,j) * ks%psi_c(j)
      end do
      ks%hpsi(i) = tmp
    end do

  end subroutine

end module
