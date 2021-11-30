module shop
  use prec
  use fileio
  use hamil
  use couplings
  implicit none

  contains

  ! initialize the random seed from the system clock
  ! code from: http://fortranwiki.org/fortran/show/random_seed
  subroutine init_random_seed()
    implicit none
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)

    deallocate(seed)
  end subroutine

  subroutine whichToHop(tion, ks, which)
    implicit none

    integer, intent(in) :: tion
    integer, intent(inout) :: which
    type(TDKS), intent(in) :: ks

    integer :: i
    real(kind=q) :: lower, upper, r

    which = 0
    call random_number(r)

    do i=1, ks%ndim
      if (i == 1) then
        lower = 0
        upper = ks%sh_prop(i,tion)
      else
        lower = upper
        upper = upper + ks%sh_prop(i,tion)
      end if
      if (lower <= r .AND. r < upper) then
        which = i
        exit
      end if
    end do

  end subroutine

  subroutine calcprop_LBS(tion, cstat, ks, inp, olap)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion
    integer, intent(in) :: cstat
    type(overlap), intent(in) :: olap

    real(kind=q) :: Akk
    complex(kind=q), allocatable :: epcoup(:)
    integer :: i, iq

    allocate(epcoup(ks%ndim))
    do i=1,ks%ndim
      iq = olap%kkqmap(cstat, i)
      epcoup(i) = SUM(olap%EPcoup(cstat,i,:,:,1) * ks%PhQ(iq,:,:,tion))
    end do

    Akk = CONJG(ks%psi_a(cstat, tion)) * ks%psi_a(cstat, tion)
    ks%Bkm = -2. / hbar * AIMAG( CONJG(ks%psi_a(cstat, tion)) * &
             ks%psi_a(:, tion) * epcoup(:) )

    ks%sh_prop(:,tion) = ks%Bkm / Akk * inp%POTIM * ks%sh_Bfactor(cstat,:)
    forall (i=1:ks%ndim, ks%sh_prop(i,tion) < 0) ks%sh_prop(i,tion) = 0

  end subroutine

  subroutine calcprop(tion, cstat, ks, inp)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion
    integer, intent(in) :: cstat

    integer :: i, j
    real(kind=q) :: Akk

    Akk = CONJG(ks%psi_a(cstat, tion)) * ks%psi_a(cstat, tion)

    if (inp%LEPC) then
      ks%Bkm = -2. / hbar * AIMAG( CONJG(ks%psi_a(cstat, tion)) * &
               ks%psi_a(:, tion) * ks%EPcoup(cstat, :, tion) )
    else
      ks%Bkm = 2. * REAL(CONJG(ks%psi_a(cstat, tion)) * ks%psi_a(:, tion) * &
                    ks%NAcoup(cstat, :, tion))
      call calcBfactor(ks, inp, cstat, tion)
    end if

    ks%sh_prop(:,tion) = ks%Bkm / Akk * inp%POTIM * ks%sh_Bfactor(cstat,:)
    forall (i=1:ks%ndim, ks%sh_prop(i,tion) < 0) ks%sh_prop(i,tion) = 0
    ! write(*,*) (ks%Bkm(i), i=1, ks%ndim)
    ! write(*,*) (ks%sh_prop(i, tion), i=1, ks%ndim)

  end subroutine

  ! calculate surface hopping probabilities
  subroutine runSH(ks, inp, olap)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    type(overlap), intent(in) :: olap
    integer :: i, j, tion, Nt
    integer :: istat, cstat, which

    ks%sh_pops = 0
    ks%sh_prop = 0
    istat = inp%BASSEL(inp%INIKPT, inp%INIBAND)
    Nt = inp%NAMDTIME / inp%POTIM

    ! initialize the random seed for ramdom number production
    call init_random_seed()

    if (inp%LEPC) call calcBftot(ks, inp)

    do i=1, inp%NTRAJ
      ! in the first step, current step always equal initial step
      cstat = istat
      do tion=1, Nt
        if (inp%LARGEBS) then
          call calcprop_LBS(tion, cstat, ks, inp, olap)
        else
          call calcprop(tion, cstat, ks, inp)
        end if
        call whichToHop(tion, ks, which)
        if (which > 0) cstat = which
        ks%sh_pops(cstat, tion) = ks%sh_pops(cstat, tion) + 1
      end do
    end do

    ks%sh_pops = ks%sh_pops / inp%NTRAJ
    ! ks%sh_prop = ks%sh_prop / inp%NTRAJ

    ! do tion=1, Nt
    !   write(*,*) (ks%sh_pops(i,tion), i=1, ks%ndim)
    ! end do
  end subroutine

  ! Calculate Boltzmann factors for SH probability correction.
  subroutine calcBfactor(ks, inp, cstat, tion)
    implicit none
    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion
    integer, intent(in) :: cstat

    integer :: jb, nb
    real(kind=q) :: dE, kbT

    nb = ks%ndim
    kbT = inp%TEMP * BOLKEV
    ks%sh_Bfactor(cstat,:) = 1.0_q

    if (inp%LHOLE) then
      do jb=1, nb
        dE = ks%eigKs(jb, tion) - ks%eigKs(cstat,tion)
        if (dE<0) then
          ks%sh_Bfactor(cstat,jb) = exp(dE / kbT)
        end if
      end do
    else
      do jb=1, nb
        dE = ks%eigKs(jb,tion) - ks%eigKs(cstat, tion)
        if (dE>0) then
          ks%sh_Bfactor(cstat,jb) = exp(-dE / kbT)
        end if
      end do
    end if

  end subroutine

  subroutine calcBftot(ks, inp)
    ! Calculate total Boltzmann factors only once
    ! if ks eigs are invariable.
    implicit none
    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp

    integer :: ib, nb

    nb = ks%ndim
    do ib=1, nb
      call calcBfactor(ks, inp, ib, 1)
    end do
  end subroutine

  ! need a subroutine here to write the results we need
  subroutine printSH(ks, inp)
    implicit none
    type(TDKS), intent(in) :: ks
    type(namdInfo), intent(in) :: inp

    integer :: i, j, tion, Nt, ierr, io
    character(len=48) :: buf

    write(buf, *) inp%NAMDTINI
    open(unit=24, file='SHPROP.' // trim(adjustl(buf)), &
         status='unknown', action='write', iostat=ierr)
    open(unit=25, file='PSICT.' // trim(adjustl(buf)), &
         status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "SHPROP file I/O error!"
      stop
    end if

    if (.NOT. inp%LEPC) then

      do io = 24, 25
        write(io,'(A,A12,A3,I5)') '#', 'BMIN',     ' = ', inp%BMIN
        write(io,'(A,A12,A3,I5)') '#', 'BMAX',     ' = ', inp%BMAX
        write(io,'(A,A12,A3,I5)') '#', 'INIBAND',  ' = ', inp%INIBAND
        write(io,'(A,A12,A3,I5)') '#', 'NBANDS',   ' = ', inp%NBANDS

        write(io,'(A,A12,A3,I5)')   '#', 'NSW',    ' = ', inp%NSW
        write(io,'(A,A12,A3,F5.1)') '#', 'POTIM',  ' = ', inp%POTIM
        write(io,'(A,A12,A3,F5.1)') '#', 'TEMP',   ' = ', inp%TEMP

        write(io,'(A,A12,A3,I5)') '#', 'NAMDTINI', ' = ', inp%NAMDTINI
        write(io,'(A,A12,A3,I5)') '#', 'NAMDTIME', ' = ', inp%NAMDTIME
        write(io,'(A,A12,A3,I5)') '#', 'NTRAJ',    ' = ', inp%NTRAJ
        write(io,'(A,A12,A3,I5)') '#', 'NELM',     ' = ', inp%NELM

        write(io,'(A,A12,A3,A)')  '#', 'RUNDIR',   ' = ', TRIM(ADJUSTL(inp%rundir))
        write(io,'(A,A12,A3,L5)') '#', 'LHOLE',    ' = ', inp%LHOLE
        write(io,'(A,A12,A3,L5)') '#', 'LSHP',     ' = ', inp%LSHP
        write(io,'(A,A12,A3,L5)') '#', 'LCPTXT',   ' = ', inp%LCPTXT
        write(io,'(A,A12,A3,L5)') '#', 'LGAMMA',   ' = ', inp%LGAMMA
      end do

    else

      do io = 24, 25
        write(io,'(A,A12,A3,I6)')     '#', 'BMIN',     ' = ', inp%BMIN
        write(io,'(A,A12,A3,I6)')     '#', 'BMAX',     ' = ', inp%BMAX
        write(io,'(A,A12,A3,I6)')     '#', 'KMIN',     ' = ', inp%KMIN
        write(io,'(A,A12,A3,I6)')     '#', 'KMAX',     ' = ', inp%KMAX
        if (inp%EMIN > -1.0E5_q) &
          write(io,'(A,A12,A3,F6.2)') '#', 'EMIN',     ' = ', inp%EMIN
        if (inp%EMAX <  1.0E5_q) &
          write(io,'(A,A12,A3,F6.2)') '#', 'EMAX',     ' = ', inp%EMAX

        write(io,'(A,A12,A3,I6)')     '#', 'NBASIS',   ' = ', inp%NBASIS
        write(io,'(A,A12,A3,I6)')     '#', 'NBANDS',   ' = ', inp%NBANDS
        write(io,'(A,A12,A3,I6)')     '#', 'NKPOINTS', ' = ', inp%NKPOINTS
        write(io,'(A,A12,A3,I6)')     '#', 'INIBAND',  ' = ', inp%INIBAND
        write(io,'(A,A12,A3,I6)')     '#', 'INIKPT',   ' = ', inp%INIKPT

        write(io,'(A,A12,A3,I6)')     '#', 'NSW',      ' = ', inp%NSW
        write(io,'(A,A12,A3,F6.1)')   '#', 'POTIM',    ' = ', inp%POTIM
        write(io,'(A,A12,A3,F6.1)')   '#', 'TEMP',     ' = ', inp%TEMP
        write(io,'(A,A12,A3,I6)')     '#', 'NAMDTINI', ' = ', inp%NAMDTINI
        write(io,'(A,A12,A3,I6)')     '#', 'NAMDTIME', ' = ', inp%NAMDTIME
        write(io,'(A,A12,A3,I6)')     '#', 'NTRAJ',    ' = ', inp%NTRAJ
        write(io,'(A,A12,A3,I6)')     '#', 'NELM',     ' = ', inp%NELM

        write(io,'(A,A12,A3,L6)')     '#', 'LEPC',     ' = ', inp%LEPC
        write(io,'(A,A12,A3,I6)')     '#', 'EPCTYPE',  ' = ', inp%EPCTYPE
        write(io,'(A,A12,A3,L6)')     '#', 'LBASSEL',  ' = ', inp%LBASSEL
        write(io,'(A,A12,A3,L6)')     '#', 'LSORT',    ' = ', inp%LSORT
        write(io,'(A,A12,A3,L6)')     '#', 'LCPTXT',   ' = ', inp%LCPTXT
        write(io,'(A,A12,A3,L6)')     '#', 'LHOLE',    ' = ', inp%LHOLE

        if (inp%EPCTYPE==2) &
          write(io,'(A,A12,A3,A)') '#', 'MDFIL', ' = ', TRIM(ADJUSTL(inp%FILMD))
        write(io,'(A,A12,A3,A)') '#', 'EPMFIL', ' = ', TRIM(ADJUSTL(inp%FILEPM))
      end do

    end if

    Nt = inp%NAMDTIME / inp%POTIM
    do tion=1, Nt
      write(unit=24, fmt='(*(G20.10))') &
            tion * inp%POTIM, SUM(ks%eigKs(:,tion) * ks%sh_pops(:,tion)), &
            (ks%sh_pops(i,tion), i=1, ks%ndim)
      ! write(unit=25, fmt="(2G20.10, *( ' ( ',G20.10,' , ',G20.10,' ) ' ) )") &
      write(unit=25, fmt='(*(G20.10))') &
            tion * inp%POTIM, SUM(ks%eigKs(:,tion) * ks%pop_a(:,tion)), &
          ! (ks%psi_a(i,tion), i=1, ks%ndim)
            (ks%pop_a(i,tion), i=1, ks%ndim)
    end do

    close(24)
    close(25)

  end subroutine

end module
