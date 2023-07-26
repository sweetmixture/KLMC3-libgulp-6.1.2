  subroutine gulpmain(iret, ichemsh_link_loc, MPI_comm_in )
!
!  General Utility Lattice Program including
!
!     (1) Electrostatic potential calculations
!     (2) Electronegativity equalisation method
!     (3) Interatomic potential fitting
!     (4) Structure optimisation
!     (5) Phonon calculations
!     (6) Defect calculations
!     (7) Molecular dynamics
!     (8) Free energy minimisation
!     (9) Calculations on 2-D and 3-D systems
!    (10) Monte Carlo calculations
!    (11) Continuum solvation models
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  11/96 Dynamic memory for large arrays added to f77 version
!  11/96 Scan moved into main routine for memory handling
!   1/97 Dynamic memory introduced by direct call of malloc
!        or via a C subroutine call
!  12/97 Reduction of memory added when only EEM/QEq uses derv2
!        and passing of derv2 to MD routines added so that 
!        variable charges can be used in MD.
!   4/98 Input is now passed through twice and stored temporarily
!        on channel 4.
!   5/98 Dynamic memory allocation for Linux phonons moved to
!        top level as this is the only way to make this OS work.
!   5/98 Structure of defect calls rearranged as this avoids 
!        problems under Linux - avoid any call to free as this
!        causes Linux to core dump
!   3/99 MPI parallelisation of MD introduced into standard version
!   6/00 Genetic Algorithm and Simulated Annealing Routine added (SMW)
!   7/00 Moved call to predict to within loop over configurations (SMW)
!  11/99 Modifications for neutron scattering added
!  10/00 Program converted to f90
!  12/00 2-D systems added
!   1/01 Monte Carlo calculations added
!   4/01 Method/derivative checking routine added
!   4/01 Option calls moved to separate subroutine
!   6/01 COSMO solvation model added
!   6/03 XML modifications added
!   9/03 Parallel I/O modified
!  11/04 Manu changes - see HISTORY
!   9/06 Setup & finish separated for QM/MM benefit
!   3/07 Modified to fit in with Chemshell
!   5/07 GULPfinish call removed and placed in top level routine
!   5/07 Option to explicitly open inputfile added (VM)
!   6/07 Handling of inputfile removed from this routine
!   7/07 Handling of inputfile added back
!   3/09 Communicator for MPI now passed in as argument
!   6/09 Call to options renamed to gulp_options
!   4/17 ChemShell interaction modified
!   6/17 Old ChemShell restored as a compile option
!
!  Julian Gale, CIC, Curtin University, June 2017
!
!  With contributions by :
!
!  Ivan Walton,         HPCI, Southampton,     March  1999
!  Scott Woodley,       Royal Institution,     June   2000
!
  use control
  use gulpchemsh
  use iochannels
  use parallel
#ifdef ACCELRYS
  use license
#endif
#ifdef TRACE
  use trace,     only : init_trace, trace_in, trace_out
#endif
  implicit none
  integer(i4), intent(out)  :: iret
  integer(i4), intent(in)   :: ichemsh_link_loc
  integer*4,   intent(in)   :: MPI_comm_in
!
!  Local variables
!
#ifdef ACCELRYS
  integer(i4)               :: ierror
#endif
!
!  ChemShell additions
!
  character(len=132)        :: infile
  character(len=132)        :: outfile
  logical                   :: linquire
  logical                   :: opc_root
  external opc_root

!*****************************
!  ChemShell initialisation  *
!*****************************
#ifdef OLDCS
  ichemsh_qm = ichemsh_link_loc
! wkjee - ichemsh_qm = -1
#else
  ichemsh_link = ichemsh_link_loc
  ! wkjee - = 0 : used arg
  ! wkjee - ichemsh <attribute::save> - in module gulpchemsh (in file modules.F90)
#endif

  loprt = .true.
  ! wkjee - loprt - in module gulpchemsh (in file modules.F90)

! =============================================================
! wkjee - block blow not working in standard gulp mode
! =============================================================
#ifdef OLDCS
  if (opc_root() .and. ichemsh_qm .ge. 0) then
#else
  if (opc_root() .and. ichemsh_link .eq. 1) then

    ! wkjee - block within this 'if' is not necessary
    ! write(*,'(A)') "in gulpmain, this is not even working"
#endif
    ioin = 15
    inquire(unit=ioin,opened=linquire)
    if (linquire) close(ioin)
    ioout = 16
    inquire(unit=ioout,opened=linquire)
    if (linquire) close(ioout)
     
    call GetGulpFileNames(infile, outfile, 132 )
     
    if (outfile .ne. "stdout")then
      open(unit=ioout,file=outfile,form='formatted')
    else
      ioout = 6
    endif
    open(unit=ioin,file=infile,form='formatted')
  endif
! =============================================================
! wkjee - block above not working in standard gulp mode
! =============================================================

! wkjee - anyway the initcomms called

#ifdef KLMC_DEBUG_GULPMAIN
  write(*,'(A)') "in gulpmain: before calling GULP_initcomms() before"
#endif
  call GULP_initcomms(MPI_comm_in)

#ifdef KLMC_DEBUG_GULPMAIN
  write(*,'(A)') "in gulpmain: before calling GULP_initcomms() after"
#endif

#ifdef OLDCS
  if (ichemsh_qm .lt. 0) then
! wkjee - ichemsh_qm = -1 .. here is the case
#else
  if (ichemsh_link .eq. 0) then
  ! wkjee - ichemsh_link = 0 .. or here is the case? - i.e., GULP standard mode is using
#ifdef KLMC_DEBUG_GULPMAIN
    if(ioproc) then
      write(*,'(A,I4)') "in gulpmain, ichemsh_link (klmc_link): ", ichemsh_link
    end if
#endif

#endif

    !
    ! wkjee - setupinputoutput routine should be modified - in order to take the input and output files as its argument
    !
#ifdef KLMC_DEBUG_GULPMAIN
    if(ioproc) then
      write(*,'(A)') "in gulpmain: before call setupinputoutput" 
    end if
#endif
    call setupinputoutput
#ifdef KLMC_DEBUG_GULPMAIN
    if(ioproc) then
      write(*,'(A)') "in gulpmain: after  call setupinputoutput" 
    end if
#endif
    ! wkjee end

  endif

#ifdef TRACE
! wkjee - TRACE is not defined in the default mode
!*************************
!  Trace initialisation  *
!*************************
  call init_trace
  call trace_in('gulpmain')
  ! wkjee - testing before setup
  ! write(*,'(A)') "in gulpmain, after calling 'init_trace' and 'trace_in'"
#endif
!***************
!  Setup GULP  *
!***************
  ! wkjee- check if 'gulp_setup' does inputfile read
#ifdef KLMC_DEBUG_GULPMAIN
  if(ioproc) write(*,'(A)') "in gulpmain: before call gulp_setup" 
#endif
  call gulp_setup
  ! wkjee - gulp_setup ... doing with reading input files
#ifdef KLMC_DEBUG_GULPMAIN
  if(ioproc) write(*,'(A)') "in gulpmain: after call gulp_setup"
#endif
!*****************
!  Call options  *
!*****************
  if (lfit) then
    call fit
#ifdef ACCELRYS
    call sendHeartbeat(ierror)
    if (ierror /= 0) call gulpfinish
#endif
  endif
      
#ifdef OLDCS
  if (ichemsh_qm.ne.99) call gulp_options
#else
  if (ichemsh_output.ne.1) call gulp_options
#endif

#ifdef ACCELRYS
!
!  Perform final tasks - cannot close scratch files under Linux as this causes a crash sometime!
!
  if (ioproc) then
    call license_checkin(ierror)
  endif
#endif

!
!  Final Chemshell tasks
!
  iret = 0
#ifdef OLDCS
  if (ichemsh_qm.ge.0) then
#else
  if (ichemsh_link.eq.1) then
#endif
    if (ioout.ne.6) close (ioout)
    close (ioin)
  endif
#ifdef TRACE
  call trace_out('gulpmain')
#endif

  !
  ! wkjee - may need to close ioin / ioout here
  !

  end subroutine gulpmain

!    
!  03/Jul/2007 aperlov:
!  This routine opens input/output files  if seedname is an argument of the command line
!  (output file is open only for the ioproc)
!  All this stuff is actually needed for parallel execution under vista where hp-mpi
!  is not capable of dealing with redirection properly, but is harmless and can be used for any OS.
!

!  wkjee - see iomod.F90 / module iochannels / module gulp_lengths
!  wkjee - see modules.F90 / module parallel
!  wkjee - seems doesn't do much on read input files // actual standard read inputfile is in 'call gulp_setup' - gulpsetup.F90

  subroutine setupinputoutput

  use parallel,   only : ioproc
  use iochannels, only : ioin,ioout,lioproconly
#ifdef KLMC
  use klmc
#endif

  implicit none
  integer            :: num_args
  character(len=132) :: seedname

#ifdef KLMC
  ! module klmc character :: gulp_klmc_iopath
  ! character(len=256) :: tmpstring
  character(len=256) :: klmc_gulp_input
  character(len=256) :: klmc_gulp_output
  logical lopen
#endif

  num_args = command_argument_count()
  seedname = ""
  if (num_args > 0) then
    call get_command_argument(1,seedname)
  endif

  ! wkjee - testing
  ! write(*,'(A,A)')  "in gulpmain: setupinputoutput (seedname): ", seedname
  ! write(*,'(A,I4)') "in gulpmain: setupinputoutput (num_args): ", num_args
  ! wkjee
!  
!  Open input file explicitly
! 

  ! wkjee - klmc 05 July 2023
#ifdef KLMC
  ! wkjee - using canonical form of inputfile name "gulp_klmc" always
  seedname = "gulp_klmc"
  if ( gulp_klmc_iopath .ne. "" ) then
    ! wkjee - setting gulpklmc inputoutput files
    klmc_gulp_input = trim(adjustl(gulp_klmc_iopath))//"/"//trim(seedname)//".gin"
    klmc_gulp_output = trim(adjustl(gulp_klmc_iopath))//"/"//trim(seedname)//".gout"
  endif
#endif
  ! wkjee - klmc io tested 05 July 2023
  ! 	seedname is set to 'abcd', if in the directory there is a file 'abcd.gin' then gulp will run normally
  !     no need to modify 'gulpsetup.F90' responsible for the default io
  ! end wkjee
  if (seedname .ne. "") then
!
!  Set a flag that indicates that all I/O in parallel should be via the I/O proc
!
    ! wkjee - check if seedname working
#ifdef KLMC_DEBUG_GULPMAIN
    if(ioproc) then
      write(*,'(A,A)') "in gulpmain: gulp_klmc_iopath  : ", gulp_klmc_iopath
      write(*,'(A,A)') "in gulpmain: klmc_gulp_input   : ", klmc_gulp_input
      write(*,'(A,A)') "in gulpmain: klmc_gulp_output: : ", klmc_gulp_output
    end if
#endif
    ! wkjee - lioproconly (module iochannels from iomod.F90) default = .false.
    lioproconly = .true.
    ! wkjee - set this .true. for cmmand line inputfile name getting
    ! wkjee - channels ioin / ioout responsible for read and writing !!! / iomod.F90 , module iochannels
    ! wkjee - if hijacked
    ! wkjee - 'ioin' and 'ioout' will be set and go through the subroutines
    !	channels, firstpass, secondpass
    ! end wkjee
    ! wkjee 
    ! write(*,'(A)') "in gulpmain: if ioin open works before"
#ifdef KLMC
    open(unit=ioin,file=klmc_gulp_input,form='formatted',status='unknown')
#else
    open(unit=ioin,file=trim(seedname)//".gin",form='formatted')
#endif
    ! wkjee - 06 July 2023 - instead of passing the parameters, but explicit path, it worked !! why??
    ! write(*,'(A)') "in gulpmain: if ioin open works after"
    ! wkjee
    if (ioproc) then
      ! wkjee
      ! inquire(unit=ioout,opened=lopen)
      ! write(*,'(A,L)') "in gulpmain: lopen<logical>: ", lopen
      ! wkjee - checked, lopen = .false.
      ! if( lopen ) then
      !   write(*,'(A)') "in gulpmain: opened"
      ! end if
      ! wkjee
      ! wkjee - logical 'ioproc'
      ! write(*,'(A)') "in gulpmain: if ioout open works before"
#ifdef KLMC
      ! open(unit=ioout,file=klmc_gulp_output,form='formatted',status='new')
      open(unit=ioout,file=klmc_gulp_output,form='formatted',status='unknown')
#else
      open(unit=ioout,file=trim(seedname)//".gout",form='formatted')
#endif
      ! open(unit=ioout,file=trim(seedname)//".gout",form='formatted',status='unknown')
      ! open(unit=ioout,file="/work/e05/e05/wkjee/Software/gulp-6.1.2/Src/Custom/path_test/gulp_klmc.gout",form='formatted',status='unknown')
      ! wkjee
      ! open(unit=ioout,file=trim(adjustl(seedname))//".gout",form='formatted')
      ! write(*,'(A)') "in gulpmain: if ioout open works after"
    endif
  endif

  end subroutine setupinputoutput
