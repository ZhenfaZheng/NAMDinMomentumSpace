module fileio
  use prec

  implicit none

  type namdInfo
    integer :: BMIN
    integer :: BMAX
    integer :: NBASIS      ! No. of adiabatic states as basis
    integer :: NBANDS      ! No. of band of the system
    integer :: INIBAND     ! inititial adiabatic state of excited electron/hole
    integer :: NSW         ! No. of MD steps
    integer :: NAMDTINI    ! Initial time step of NAMD
    integer :: NAMDTIME    ! No. of steps of NAMD
    integer, allocatable, dimension(:) :: NAMDTINI_A    ! No. of steps of NAMD
    integer, allocatable, dimension(:) :: INIBAND_A     ! No. of steps of NAMD
    integer :: NTRAJ       ! No. of surface hopping trajectories
    integer :: NELM        ! No. of steps of electron wave propagation
    integer :: NSAMPLE     ! No. of steps of electron wave propagation

    real(kind=q) :: POTIM  ! Time step of MD run
    real(kind=q) :: TEMP   ! MD Temperature

    ! hole or electron surface hopping
    logical :: LHOLE
    ! whether to perform surface hopping, right now the value is .TRUE.
    logical :: LSHP
    logical :: LCPTXT

    ! running directories
    character(len=256) :: RUNDIR
    character(len=256) :: TBINIT

    ! whether the WAVECARs come from gamma version VASP or not.
    ! for other version WAVECARs, the NA couplings are complex number.
    ! for gamma version WAVECARs, the NA couplings are real number.
    ! NA coupling = <psi_i(t)| d/dt |(psi_j(t))>
    logical :: LGAMMA

    ! whether to use e-p matrix as NA couplings
    logical :: LEPC
    integer :: EPCTYPE
    integer :: KMIN
    integer :: KMAX
    integer :: NKPOINTS    ! No. of k-points of the system
    integer :: INIKPT      ! inititial k point of excited electron/hole
    integer, allocatable, dimension(:) :: INIKPT_A      ! No. of steps of NAMD
    integer, allocatable, dimension(:,:) :: BASSEL
    ! selected basises among the nk*nb eigen states
    real(kind=q) :: EMIN
    real(kind=q) :: EMAX

    character(len=256) :: FILEPC  ! epc file from EPW package, if LEPC=.TRUE.
    character(len=256) :: FILPH   ! phonon modes file from QE-PH, if LEPC=.TRUE.
    character(len=256) :: FILMD   ! MD trajetory (XDATCAR) from VASP, if LEPC=.TRUE.

  end type

  contains

    subroutine getUserInp(inp)
      implicit none

      type(namdInfo), intent(inout) :: inp

      ! local variables with the same name as those in "inp"
      integer :: bmin
      integer :: bmax
      integer :: nsw
      integer :: iniband
      integer :: nbands
      integer :: nkpoints
      integer :: namdtime
      ! integer :: namdtini    
      integer :: ntraj
      integer :: nelm
      integer :: nsample
      real(kind=q) :: potim
      real(kind=q) :: temp

      ! hole or electron surface hopping
      logical :: lhole
      ! surface hopping?
      logical :: lshp
      logical :: lcpext
      logical :: lgamma
      logical :: lepc
      ! running directories
      character(len=256) :: rundir
      character(len=256) :: tbinit
      ! files for e-p coupling calculation, if LEPC=.TRUE.
      character(len=256) :: filepc, filph, filmd
      integer :: epctype
      integer :: kmin
      integer :: kmax
      integer :: ik, ib, num
      real(kind=q) :: emin
      real(kind=q) :: emax


      namelist /NAMDPARA/ bmin, bmax, nsw,      &
                          nbands, nkpoints,     &
                          potim, ntraj, nelm,   &
                          temp, rundir,         &
                          lhole, lshp, lcpext,  &
                          lgamma, lepc, epctype,&
                          filepc, filph, filmd, &
                          kmin, kmax, emin, emax, &
                          namdtime,             &
                          nsample, tbinit

      integer :: ierr, i
      logical :: lext

      ! set default values for thos parameters
      rundir = 'run'
      tbinit = 'INICON'
      bmin = 0
      bmax = 0
      nbands = 0
      nkpoints = 1
      ! iniband = 0
      ntraj = 1000
      nelm = 1000
      lhole = .FALSE.
      lshp = .TRUE.
      ! namdtini = 1
      namdtime = 200
      potim = 1.0_q
      temp = 300.
      lcpext = .FALSE.
      lgamma = .TRUE.
      lepc = .FALSE.
      epctype = 1
      filepc = 'prefix.epc'
      filph = 'prefix.matdyn.modes'
      filmd = 'XDATCAR'
      kmin = 1
      kmax = 0
      emin = -99999.9
      emax =  99999.9

      open(file="inp", unit=8, status='unknown', action='read', iostat=ierr)
      if ( ierr /= 0 ) then
        write(*,*) "I/O error with input file: 'inp'"
      end if

      read(unit=8, nml=NAMDPARA)
      close(unit=8)

      ! if No. of k-points more than 1, lgamma is automatically set to False.
      if ( nkpoints > 1) lgamma = .FALSE.
      if ( kmax == 0 ) kmax = nkpoints

      allocate(inp%BASSEL(nkpoints,nbands))
      allocate(inp%INIBAND_A(nsample), inp%NAMDTINI_A(nsample))
      allocate(inp%INIKPT_A(nsample))

      inquire(file=tbinit, exist=lext)
      if (.NOT. lext) then
        write(*,*) "File containing initial conditions does NOT exist!"
      else
        open(unit=9, file=tbinit, action='read')
        if (lepc) then
          do i=1, nsample
            read(unit=9,fmt=*) inp%NAMDTINI_A(i), inp%INIBAND_A(i), inp%INIKPT_A(i)
          end do
        else
          inp%INIKPT_A = 1
          do i=1, nsample
            read(unit=9,fmt=*) inp%NAMDTINI_A(i), inp%INIBAND_A(i)
          end do
        end if
        close(9)
      end if

      ! do some checking...
      ! put the following checks in the future version
      ! if (bmin <= 0 .OR. bmax <= 0 .OR. bmin >= bmax) then
      !   write(*,*) "Please specify the correct BMIN/BMAX"
      !   stop
      ! end if

      ! if (iniband == 0 .OR. iniband < bmin .OR. iniband > bmax) then
      !   write(*,*) "Please specify the correct initial band!"
      !   stop
      ! end if

      ! if (nbands == 0) then
      !   write(*,*) "I need the No. of bands..."
      !   stop
      ! end if

      ! if (namdtini + namdtime - 1 > nsw) then
      !   write(*,*) "NAMDTIME too long..."
      !   stop
      ! end if
      ! here ends the currently simplified version of parameter checking

      ! assign the parameters
      inp%BMIN     = bmin
      inp%BMAX     = bmax
      inp%NBASIS   = bmax - bmin + 1
      inp%NSW      = nsw
      inp%NBANDS   = nbands
      inp%NKPOINTS = nkpoints
      inp%NAMDTIME = namdtime
      inp%NTRAJ    = ntraj
      inp%NELM     = nelm
      inp%LHOLE    = lhole
      inp%LSHP     = lshp
      inp%RUNDIR   = trim(rundir)
      inp%TBINIT   = trim(tbinit)
      inp%NSAMPLE  = nsample
      inp%POTIM    = potim
      inp%LCPTXT   = lcpext
      inp%LGAMMA   = lgamma
      inp%LEPC     = lepc
      inp%EPCTYPE  = epctype
      inp%TEMP     = temp
      inp%FILEPC   = trim(filepc)
      inp%FILPH    = trim(filph)
      inp%FILMD    = trim(filmd)
      inp%KMIN     = kmin
      inp%KMAX     = kmax
      inp%NBASIS   = inp%NBASIS * ( kmax - kmin + 1 )
      inp%EMIN     = emin
      inp%EMAX     = emax

      num = 0
      inp%BASSEL = -1 
      !! For array element inp%BASSEL(ik,ib),
      !! -1 represent ik,ib state isn't selected as basis;
      !! number >0 represent ik,ib state is selected as basis,
      !! and the number is the state's serial number among the basises.
      do ik = kmin, kmax
        do ib = bmin, bmax
          num = num + 1
          inp%BASSEL(ik,ib) = num
        end do
      end do

    end subroutine

    ! Need a subroutine to print out all the input parameters
    subroutine printUserInp(inp)
      implicit none
      type(namdInfo), intent(in) :: inp

      write(*,'(A)') "------------------------------------------------------------"
      write(*,'(A30,A3,I5)') 'BMIN',     ' = ', inp%BMIN
      write(*,'(A30,A3,I5)') 'BMAX',     ' = ', inp%BMAX
      write(*,'(A30,A3,I5)') 'INIBAND',  ' = ', inp%INIBAND
      write(*,'(A30,A3,I5)') 'INIKPT',   ' = ', inp%INIKPT
      write(*,'(A30,A3,I5)') 'NBANDS',   ' = ', inp%NBANDS
      write(*,'(A30,A3,I5)') 'NKPOINTS', ' = ', inp%NKPOINTS

      write(*,'(A30,A3,I5)')   'NSW',    ' = ', inp%NSW
      write(*,'(A30,A3,F5.1)') 'POTIM',  ' = ', inp%POTIM
      write(*,'(A30,A3,F5.1)') 'TEMP',   ' = ', inp%TEMP

      write(*,'(A30,A3,I5)') 'NAMDTINI', ' = ', inp%NAMDTINI
      write(*,'(A30,A3,I5)') 'NAMDTIME', ' = ', inp%NAMDTIME
      write(*,'(A30,A3,I5)') 'NTRAJ',    ' = ', inp%NTRAJ
      write(*,'(A30,A3,I5)') 'NELM',     ' = ', inp%NELM

      write(*,'(A30,A3,A)')  'RUNDIR',   ' = ', TRIM(ADJUSTL(inp%rundir))
      write(*,'(A30,A3,L5)') 'LHOLE',    ' = ', inp%LHOLE
      write(*,'(A30,A3,L5)') 'LSHP',     ' = ', inp%LSHP
      write(*,'(A30,A3,L5)') 'LCPTXT',   ' = ', inp%LCPTXT
      write(*,'(A30,A3,L5)') 'LGAMMA',   ' = ', inp%LGAMMA
      write(*,'(A30,A3,L5)') 'LEPC',     ' = ', inp%LEPC
      write(*,'(A30,A3,I5)') 'EPCTYPE',  ' = ', inp%EPCTYPE
      write(*,'(A30,A3,I5)') 'KMIN',     ' = ', inp%KMIN
      write(*,'(A30,A3,I5)') 'KMAX',     ' = ', inp%KMAX

      write(*,'(A)') "------------------------------------------------------------"
    end subroutine 

end module
