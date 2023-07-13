  subroutine GULP_initcomms(MPI_comm_in)
!
!  Initialises MPI if necessary, finds own taskid and number of tasks
!
!  Modified to handle possible precision issues - pass local scalars
!  to get procid/nprocs for benefit of Cray
!
!  12/03 Silent option added in which ioproc is set to false
!   3/07 Chemshell modifications added and renamed 
!   3/09 MPI communicator changed to MPI_comm_GULP based on passed argument
!   6/09 MPI barrier changed to MPI_comm_world instead of MPI_comm_GULP
!  11/16 Blacs initialisation added
!   4/17 ChemShell interaction modified
!   6/17 Old ChemShell restored as a compile option
!   9/20 Blacs context now set to be MPI_Comm_GULP for ChemShell
!
!  Julian Gale, CIC, Curtin University, September 2020
!

!
!  Modules
!
  use iochannels
  use gulpchemsh
  use parallel
  ! wkjee - in modules.F90, module parallel: MPI_comm_GULP

  implicit none
!
!  Passed variables
!
  integer*4,   intent(in)   :: MPI_comm_in
!
#ifdef MPI
  include 'mpif.h'
  integer ierr,lprocid,lnprocs
  logical lmpiinit

#ifdef OLDCS
  if (ichemsh_qm .lt. 0) then
#else
  ! wkjee - gulp standard default - ichemsh_link = 0 currently  (from module gulpchemsh - modules.F90)
  ! ichemsh_link = 0 - used for the use of KLMC
  if (ichemsh_link .eq. 0) then
#endif
!
!  Non-ChemShell case - initialise MPI
!
    ! wkjee - builtin MPI_Initilaized() checking if MPI_Init() done
    call MPI_Initialized(lmpiinit, ierr)
#ifdef KLMC_DEBUG
    write(*,'(A,L4)')   "in initcomms : MPI_Init() is already called / lmpiinit<logical> : ", lmpiinit
#endif
    ! wkjee - default mode, lmpiinit = .false.
    ! wkjee - if communicator on then lmpiint = .true. - i.e., if MPI_Init() done in the taskfarm part, it should be: lmpiinit = .true.
    ! wkjee - here, lmpiinit = .true. must be set if -DKLMC flag on , forcing this subroutine to take MPI_comm_GULP = MPI_comm_in
    !
    ! wkjee - klmc modification
#ifdef KLMC 
    if ( lmpiinit ) then
#ifdef KLMC_DEBUG
      write(*,'(A)') "in initcomms : MPI_comm_GULP set to MPI_comm_in (KLMC used)"
#endif
      MPI_comm_GULP = MPI_comm_in
    end if
#endif
    ! wkjee - end klmc modification

    if (.not. lmpiinit) then
       ! wkjee
#ifdef KLMC_DEBUG
       write(*,'(A)') "in initcomms : MPI_Init() is called (must not happen)"
#endif
       call MPI_init(ierr)
!
!  Set communicator for MPI based on MPI_comm_world
!
       ! wkjee - MPI_comm_GULP in module parallel (modules.F90)
#ifdef KLMC_DEBUG
       write(*,'(A,I)') "in initcomms : MPI_comm_in (must not happen)", MPI_comm_in
#endif
       MPI_comm_GULP = MPI_comm_world
    end if
  else
!
!  ChemShell case - MPI is assumed to be already running and communicator set using argument
!
    ! wkjee - MPI_comm_GULP in module parallel (modules.F90)
    ! with condition ( ichemsh_link .eq. 0 ) lines below does not occur.
#ifdef KLMC_DEBUG
    write(*,'(A,I)')       "in initcomms : (must not happen) MPI_comm_in", MPI_comm_in
    write(*,'(A,I16,I16)') "in initcomms : (must not happen) catching MPI_comm_in / MPI_comm_GULP : ", MPI_comm_in, MPI_comm_GULP
#endif
    MPI_comm_GULP = MPI_comm_in
  endif

  call MPI_comm_rank(MPI_comm_GULP,lprocid,ierr)
  call MPI_comm_size(MPI_comm_GULP,lnprocs,ierr)

#ifdef KLMC_DEBUG
  ! wkjee - invalid communicator debugging - solved
  if(lprocid.eq.0) then
    write(*,'(A)')   "in initcomms : before getting rank / size"
    write(*,'(A,I)') "in initcomms : MPI_comm_in   <important>", MPI_comm_in
    write(*,'(A,I)') "in initcomms : MPI_comm_GULP <important>", MPI_comm_GULP
  end if
#endif
  ! wkjee


  ! wkjee - procid / nproc (module parallel in modules.F90)
  ! procid : Number for local processor
  ! nrpoc  : Number of processors for parallel execution
#ifdef KLMC_DEBUG
  write(*,'(A,I4,I4)') "in initcomms : procid / nprocs", lprocid, lnprocs
#endif
  procid  = lprocid
  nprocs  = lnprocs
!
!  Initialise Blacs for use by pblas/scalapack
!
  iBlacsContext = MPI_comm_GULP
  call blacs_gridinit( iBlacsContext, 'C', 1, lnprocs)
#else
  ! ifndef MPI - wkjee
  procid  = 0
  nprocs  = 1
#endif
  ! wkjee - MPI default
  ! wkjee - lsilent<logical> (from iomod.F90) default: .false.
  if (lsilent) then
    ioproc = .false.
  else
    ! wkjee - setting io process
    ! ioproc<logical> : If true, this is the I/O processor
    ioproc = (procid.eq.0)
  endif

#ifdef MPI
  call MPI_barrier(MPI_comm_GULP,ierr)
#endif

  return
  end
