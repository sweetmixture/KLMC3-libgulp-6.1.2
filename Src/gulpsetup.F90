  subroutine gulp_setup
!
!  Performs setup tasks for GULP
!
!   9/06 Created from gulp.F
!   3/07 initcomms renamed to GULP_initcomms
!   8/07 Call to GULP_initcomms removed to avoid MPI error
!   1/08 Declaration of ierror wrapped with ifdef
!   2/09 Modified to accommodate new version of FoX and gulp_cml
!   6/09 Call to banner changed to gulp_banner
!   6/09 Renamed to gulp_setup for consistency with other names
!   6/09 Output of extra species info added for PDF/CML
!   9/09 Number for buffer channel now accessed from module
!   7/11 Version incremented to 4.0
!   8/11 Output of hostname added to setup info
!  12/12 Duplicate use of iochannels removed
!   9/15 Opening of defect channels 41, 42, 48 removed since these
!        are no longer used
!   7/16 Modified to allow for one off allocation of pkim_model
!   4/17 ChemShell interaction modified
!   6/17 Old ChemShell restored as a compile option
!   1/18 Memory allocations checked after keywords have been set
!        through a second call to initmemory
!   1/18 Trace added
!   5/18 EEM setup added
!   8/18 Version number updated
!   8/18 KIM handling removed for version 2.
!   7/20 Version number updated
!   8/20 GFNFF added
!  10/21 lgfnff moved to control
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
!  Copyright Curtin University 2021
!
!  Julian Gale, CIC, Curtin University, October 2021
!
  use configurations
  use control
  use current
  use general
  use gulp_cml,        only : lcml, gulp_cml_init, gulp_cml_outkey
  use gulp_cml_phonon, only : gulp_cml_outspec
  use gulpchemsh
  use iochannels,      only : iotmp, ioout
  use m_pdfneutron,    only : lpdfout, outpdfin
  use parallel
#ifdef TRACE
  use trace,           only : trace_in, trace_out
#endif
#ifdef ACCELRYS
  use license
#endif
  implicit none
!
  character(len=40) :: hostname
  integer(i4)       :: iline
#ifdef ACCELRYS
  integer(i4)       :: ierror
#endif
  integer(i4)       :: hlength
  integer(i4)       :: status
  logical           :: lopened
  ! ChemShell additions
  logical           :: opc_root
  external opc_root
#ifndef OLDCS
  integer           :: ipunch
  integer           :: i
#endif
#ifdef MPI
  include 'mpif.h'
#endif
#ifdef TRACE
  call trace_in('gulp_setup')
#endif

!*************************
!  Nullify all pointers  *
!*************************
  ! wkjee - even if the pointer hasn't been previously assigned - from nullpointer.F90
  call nullpointer
!*************************
!  Initialise memory     *
!*************************
  call initmemory
!*************************
!  Initialise variables  *
!*************************
  iline = 0
  call initial
!******************
!  Output header  *
!******************
  version = 6.1_dp
#ifdef ACCELRYS
  call license_checkout(nprocs,ierror)
  if (ierror.ne.0) call gulpfinish
! Set traps for signals
  
  call setup_traps()
#endif
  call gulp_banner
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulp_setup: after banner"
#endif
!**************************
!  Get local information  *
!**************************
  call local
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulp_setup: after local"
#endif
!****************************
!  Set element information  *
!****************************
  call setele(lopened)
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A,L)') "in gulp_setup: after setele, lopened: ", lopened
#endif
  ! wkjee - lopened = .false.
!*******************************************************
!  Pre-processor passes before reading input properly  *
!*******************************************************

  ! wkjee - looks like channels (channels.F90) is responsible for reading *.gin)
  ! iotmp<integer(i4)>: temporal input buffer channel .. iomod.F90 (module iochannels)
  call channels(iotmp,.true.)
  ! wkjee - nothing really happened in the subroutine 'channels' - 05 July 2023
  ! wkjee - if commandline inputfile used in gulpmain 'setupinputoutput' this doesn't do anything

  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - before firstpass"
#endif
  ! wkjee
  call firstpass
  rewind(iotmp)
  ! wkjee - rewind does repositions the file connected to I/O unit=iotmp to the beginning of the file
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - before secondpass"
#endif
  ! wkjee
  call secondpass
  rewind(iotmp)
!*********************
!  Read in keywords  *
!*********************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - before getkeyward"
#endif
  ! wkjee
  call getkeyword(iline)
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - before setkeyward 1"
#endif
  ! wkjee
  call setkeyword
!***********************
!  Process main input  *
!***********************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - before inword"
#endif
  ! wkjee
  call inword(iline)
!*******************************************************
!  Set keywords again in case any were in the library  *
!*******************************************************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - before setkeyword 2 "
#endif
  ! wkjee
  call setkeyword
  close(iotmp,status='delete')
!****************************************
!  Set charge equilibration parameters  *
!****************************************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - before seteem "
#endif
  ! wkjee
  call seteem
!***************************
!  Output keyword details  *
!***************************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - before outkey "
#endif
  ! wkjee
  call outkey
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - after  outkey "
#endif
  ! wkjee
!*********************************************************
!  Re-initialise memory as keywords may change settings  *
!*********************************************************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - before initmemory "
#endif
  ! wkjee
  call initmemory
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - after initmemory "
#endif
  ! wkjee
!**************
!  Site name  *
!**************
  if (ioproc) then
    if (site.ne.' ') then
      write(ioout,'(''* '',a76,'' *'')') site(1:76)
      write(ioout,'(''********************************************************************************'')')
    endif
    call datetime(1_i4)
    write(ioout,'(''  Number of CPUs = '',i5,/)') nprocs
    hostname = ' '
    call get_environment_variable('HOSTNAME',hostname,hlength,status)
    if (status.eq.0) then
      write(ioout,'(''  Host name      = '',a40,/)') hostname
    elseif (status.eq.1) then
      call get_environment_variable('HOST',hostname,hlength,status)
      if (status.eq.0) then
        write(ioout,'(''  Host name      = '',a40,/)') hostname
      endif
    endif
  endif
!***********************
!  CML initialisation  *
!***********************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A,L)') "in gulpsetup - cml before: cml: ", lcml
#endif
  ! wkjee
  if (lcml) then
    call gulp_cml_init
    call gulp_cml_outkey
  endif
!*************************
!  GFNFF initialisation  *
!*************************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A,L)') "in gulpsetup - lgfnff before: lgfnff: ", lgfnff
#endif
  ! wkjee
  if (lgfnff) then
#ifndef NOGFNFF
    call gulp_gfnff_init
    call pgfnff_init
#endif
  endif
!**************************
!  One off initial set up *
!**************************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - setcfg before"
#endif
  ! wkjee
  call setcfg
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - setcfg after"
#endif
  ! wkjee
!*******************
!  Species output  *
!*******************
  if (ioproc) then
    call outspec
    if (lcml) call gulp_cml_outspec
    if (lpdfout) call outspec_pdf
  endif
!*****************************
!  Electronegativity output  *
!*****************************
  if (ioproc) then
    call outeem
  endif
!***********************************************************************
!  Check PDF related settings, convert frequencies and output PDF info *
!***********************************************************************
  if (lpdfout) call outpdfin
!************************
!  Output general info  *
!************************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - outget before"
#endif
  ! wkjee
  if (ioproc) call outgen
!*******************************
!  Output polarisability info  *
!*******************************
  if (ioproc) call outpolar
!*********************
!  Check potentials  *
!*********************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - checkpot before"
#endif
  ! wkjee
  call checkpot
!**********************
!  Output potentials  *
!**********************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - outpot before"
#endif
  ! wkjee
  call outpot
!*****************************
!  Check options for method  *
!*****************************
  ! wkjee
#ifdef KLMC_DEBUG_GULPSETUP
  write(*,'(A)') "in gulpsetup - method before"
#endif
  ! wkjee
  call methodok
!*****************************
!  ChemShell charge-only run *
!*****************************
#ifdef OLDCS
  if (ichemsh_qm .eq. 99 .and. opc_root()) then
    call ExportGulpCharges(numat, qlcfg)
  endif
#else
  if (ichemsh_output .eq. 1 .and. opc_root()) then

    ipunch = 7

    open(ipunch,file='gulp.charges',form='formatted')

    write(ipunch,*) "block=matrix records=0"
    write(ipunch,*) "block=matrix_title records=1"
    write(ipunch,*) "Gulp charges"

    write(ipunch,101) numat, 1, numat

    do i = 1,numat
        write(ipunch,100) qlcfg(i)
    enddo

    close(unit=ipunch)

  endif
100  format(f28.14)
101  format('block=dense_real_matrix records=',i6,' dimensions=',2i6)
#endif
#ifdef TRACE
  call trace_out('gulp_setup')
#endif
!
  return
  end
