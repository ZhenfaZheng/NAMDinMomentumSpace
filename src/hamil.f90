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

    complex(kind=q), allocatable, dimension(:,:,:) :: EPcoup
    complex(kind=q), allocatable, dimension(:,:,:,:) :: PhQ

    ! surface hopping related

    ! Bkm = REAL(CONJG(Akm) * Ckm)
    real(kind=q), allocatable, dimension(:) :: Bkm
    real(kind=q), allocatable, dimension(:,:) :: sh_pops
    real(kind=q), allocatable, dimension(:,:) :: sh_prop
    real(kind=q), allocatable, dimension(:,:,:) :: ph_pops
    real(kind=q), allocatable, dimension(:,:,:) :: ph_prop
    ! Blotzmann factor for SH probability scaling.
    real(kind=q), allocatable, dimension(:,:) :: sh_Bfactor

    ! whether the memory has been allocated
    logical :: LALLO = .FALSE.

  end type

  contains

  subroutine initTDKS(ks, inp, olap)
    implicit none

    type(TDKS), intent(inout)  :: ks
    type(overlap), intent(in)  :: olap
    type(namdInfo), intent(in) :: inp

    real(kind=q) :: norm
    real(kind=q), allocatable, dimension(:,:) :: eptemp
    integer :: i, j, iq, im, t, nsteps, N, Nt
    integer :: initstep

    ! memory allocation

    N = inp%NBASIS
    ks%ndim = inp%NBASIS
    Nt = inp%NAMDTIME / inp%POTIM

    if (.NOT. ks%LALLO) then
      allocate(ks%psi_c(N))
      allocate(ks%psi_p(N))
      allocate(ks%psi_n(N))
      allocate(ks%hpsi(N))
      allocate(ks%psi_a(N, Nt))
      allocate(ks%pop_a(N, Nt))
      allocate(ks%norm(Nt))

      allocate(ks%ham_c(N,N))
      allocate(ks%ham_p(N,N))
      allocate(ks%ham_n(N,N))

      allocate(ks%eigKs(N, Nt))

      if (inp%LEPC) then
        allocate(ks%PhQ(olap%NQ, olap%NMODES, 2, Nt))
        allocate(ks%ph_pops(olap%NQ, olap%NMODES, Nt))
        allocate(ks%ph_prop(N,N, olap%NMODES))
      else
        allocate(ks%NAcoup(N,N, Nt))
      end if
      ! if (.NOT. inp%LARGEBS) then
      !   allocate(ks%NAcoup(N,N, Nt))
      !   if (inp%LEPC) allocate(ks%EPcoup(N,N, Nt))
      ! else
      !   allocate(ks%PhQ(olap%NQ, olap%NMODES, 2, Nt))
      ! end if

      allocate(ks%sh_pops(N, Nt))
      allocate(ks%sh_prop(N,N))
      allocate(ks%sh_Bfactor(N,N))
      allocate(ks%Bkm(N))

      ks%LALLO = .TRUE.
    end if

    ! cero = (0, 0)
    ! uno = (1, 0)
    ks%psi_c = cero
    ks%psi_p = cero
    ks%psi_n = cero
    ! ks%ham_c = cero
    ! ks%ham_p = cero
    ! ks%ham_n = cero
    ks%psi_c( inp%BASSEL(inp%INIKPT, inp%INIBAND) ) = uno
    initstep = inp%NAMDTINI / inp%POTIM - 2
    ! initstep = MOD(initstep-1, nsw-1) + 1
    nsteps = inp%NSW - 1

    if (inp%LEPC) then

      do t=1, Nt
        i = MOD(initstep+t, nsteps) + 1
        ks%eigKs(:,t) = olap%Eig(:,i)
        ks%PhQ(:,:,:,t) = olap%PhQ(:,:,:,i)
      end do
      allocate(eptemp(inp%NMODES, 2))
      do i=1,N
        do j=1,N
          iq = olap%kkqmap(i,j)
          eptemp = ABS(olap%EPcoup(i,j,:,:,1) * ks%PhQ(iq,:,:,1)) ** 2
          norm = SUM(eptemp)
          if (norm>0) ks%ph_prop(i,j,:) = (eptemp(:,2) - eptemp(:,1)) / norm
        end do
      end do

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


  subroutine make_hamil_LBS(tion, tele, ks, inp, olap)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion, tele
    type(overlap), intent(in) :: olap

    integer :: ib, jb, iq
    complex(kind=q), allocatable :: tempQ(:,:,:)

    if (tion==1 .AND. tele==1) call calc_hamil_LBS(tion, ks, inp, olap)
    if (tele==1) call calc_hamil_LBS(tion+1, ks, inp, olap)

    ks%ham_c = ks%ham_p + (ks%ham_n - ks%ham_p) * TELE / inp%NELM

  end subroutine


  subroutine calc_hamil_LBS(tion, ks, inp, olap)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: TION
    type(overlap), intent(in) :: olap

    integer :: ib, jb, iq

    ks%ham_p = ks%ham_n
    do ib=1,ks%ndim
      do jb=ib,ks%ndim
        iq = olap%kkqmap(ib,jb)
        ks%ham_n(ib,jb) = SUM( olap%gij(ib,jb,:) * &
            SUM(ks%PhQ(iq,:,:,tion), dim=2) )
        ks%ham_n(jb,ib) = CONJG(ks%ham_n(ib,jb))
      end do
      ks%ham_n(ib,ib) = ks%ham_n(ib,ib) + ks%eigKs(ib,1)
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
