module TimeProp
  use prec
  use couplings
  use hamil

  implicit none

  contains 

  subroutine Propagation(ks, inp, olap, epc)
    implicit none
    type(TDKS), intent(inout)  :: ks
    type(namdInfo), intent(in) :: inp
    type(overlap), intent(in) :: olap
    type(epCoupling), intent(in) :: epc

    integer :: tion, tele, Nt
    integer :: i, j
    real(kind=q) :: edt
    real(kind=q) :: start, fin
    integer :: irank, ierr
    integer :: ist, iend, nbas, nbas_p

    complex(kind=q), allocatable :: psi_n_local(:)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
    ist = inp%ISTS(irank+1)
    iend = inp%IENDS(irank+1)

    nbas   = inp%NBASIS
    nbas_p = inp%NBASIS_P
    allocate(psi_n_local(nbas_p))

    Nt = inp%NAMDTIME / inp%POTIM
    edt = inp%POTIM / inp%NELM
    ! write(*,*) inp%POTIM, inp%NELM, edt

    ! the OUTER loop
    do tion = 1, Nt - 1
      ks%pop_a(:,tion) = CONJG(ks%psi_c) * ks%psi_c
      ks%norm(tion) = SUM(ks%pop_a(:,tion))
      ks%psi_a(:,tion) = ks%psi_c

      ! check the norm of the state
      ! write(*,*) tion, ks%norm(tion), ks%psi_c(:)
      ! call cpu_time(start)
      call CProp_mpi(tion, ks, inp, olap, epc)
      ! the INNER loop
      ! do tele = 1, inp%NELM-1
      ! do tele = 1, inp%NELM-1
      !   ! construct hamiltonian matrix
      !   if (inp%LEPC) then
      !     call make_hamil_EPC(tion, tele, ks, inp, olap, epc)
      !   else
      !     call make_hamil(tion, tele, ks, inp)
      !   end if
      !   ! apply hamiltonian to state vector
      !   ! call hamil_act(ks)
      !   ks%hpsi = matmul(ks%ham_c, ks%psi_c)
      !   if (tion == 1 .AND. tele == 1) then
      !     ! write(*,*) ((ks%ham_c(i,j), j=1, ks%ndim), i=1, ks%ndim)
      !     ! This is the very first step of the time propagation
      !     ! use first order difference
      !     ! [c,n,p] meas current, next, previous respectively
      !     psi_n_local = ks%psi_c(ist:iend) - imgUnit * edt * ks%hpsi / hbar
      !   else
      !     ! use second order difference
      !     psi_n_local = ks%psi_p(ist:iend) - 2 * imgUnit * edt * ks%hpsi / hbar
      !   end if
      !   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !   ! CALL MPI_ALLgather(psi_n_local, nbas_p, MPI_DOUBLE_COMPLEX, &
      !   !                    ks%psi_n,   nbas_p, MPI_DOUBLE_COMPLEX, &
      !   !                    MPI_COMM_WORLD, ierr)
      !   CALL MPI_ALLgatherv(psi_n_local, nbas_p, MPI_DOUBLE_COMPLEX, &
      !                      ks%psi_n, inp%IENDS-inp%ISTS+1, inp%ISTS-inp%ISTS(1), MPI_DOUBLE_COMPLEX, &
      !                      MPI_COMM_WORLD, ierr)
      !   ks%psi_p = ks%psi_c
      !   ks%psi_c = ks%psi_n
      ! end do
      ! end of the INNER loop
      ! call cpu_time(fin)
      ! write(*,*) "T_ion ", tion, fin - start

    end do
    tion = Nt
    ks%pop_a(:,tion) = CONJG(ks%psi_c) * ks%psi_c
    ks%norm(tion) = SUM(ks%pop_a(:,tion))
    ks%psi_a(:,tion) = ks%psi_c
    ! end of the OUTER loop
  end subroutine

  subroutine CProp_mpi(tion, ks, inp, olap, epc)
    implicit none
    integer, intent(in) :: tion
    type(TDKS), intent(inout)  :: ks
    type(namdInfo), intent(in) :: inp
    type(overlap), intent(in) :: olap
    type(epCoupling), intent(in) :: epc

    integer :: tele
    integer :: irank, nrank, ierr
    integer :: ist, iend, nbas_p
    integer, allocatable :: sendcounts(:), displs(:)
    complex(kind=q), allocatable :: psi_n_local(:)
    real(kind=q) :: edt

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nrank, ierr)

    nbas_p = inp%NBASIS_P
    ist = inp%ISTS(irank+1)
    iend = inp%IENDS(irank+1)
    edt = inp%POTIM / inp%NELM

    allocate(psi_n_local(nbas_p))
    allocate(sendcounts(nrank), displs(nrank))
    sendcounts = inp%IENDS - inp%ISTS + 1
    displs = inp%ISTS - inp%ISTS(1)

    do tele = 1, inp%NELM-1
      ! construct hamiltonian matrix
      if (inp%LEPC) then
        call make_hamil_EPC(tion, tele, ks, inp, olap, epc)
      else
        call make_hamil(tion, tele, ks, inp)
      end if
      ! apply hamiltonian to state vector
      ks%hpsi = matmul(ks%ham_c, ks%psi_c)
      if (tion == 1 .AND. tele == 1) then
        ! This is the very first step of the time propagation
        ! use first order difference
        ! [c,n,p] meas current, next, previous respectively
        psi_n_local = ks%psi_c(ist:iend) - imgUnit * edt * ks%hpsi / hbar
      else
        ! use second order difference
        psi_n_local = ks%psi_p(ist:iend) - 2 * imgUnit * edt * ks%hpsi / hbar
      end if

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      CALL MPI_ALLgatherv(psi_n_local, nbas_p, MPI_DOUBLE_COMPLEX, &
                         ks%psi_n, sendcounts, displs, MPI_DOUBLE_COMPLEX, &
                         MPI_COMM_WORLD, ierr)

      ks%psi_p = ks%psi_c
      ks%psi_c = ks%psi_n
    end do

  end subroutine

end module
