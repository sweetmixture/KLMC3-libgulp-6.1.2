  subroutine stopnow(routine)
!
!  Stops program after flushing and closing output files
!
!   9/16 Length of string for routine increased in format statement
!
  use iochannels
  use parallel
  implicit none
!
  character(len=*) :: routine
#ifdef MPI
  include 'mpif.h'
  integer ierr
#endif

  if (ioproc) then
    write (ioout,'(/,a,i5,a,a20,/)') ' Program terminated by processor ',procid,' in ',routine
    call gflush(ioout)
#ifdef ACCELRYS
    call create_killfile()
#endif
    close(4,status='delete')
    ! wkjee
#ifdef KLMC
    ! wkjee - do not close processors
    if(ioproc) then
      write (ioout,'(A)') "KLMC> in stopnow(routine): something went wrong: checkout the .gout"
      call gflush(ioout)
    end if
    ! return
#endif
  endif
#ifdef MPI
!******************************
!  Close down all processors  *
!******************************
#ifndef KLMC
! wkjee - if KLMC then do not MPI_Finalize() and do not stop
  call MPI_finalize(ierr)
  stop
#endif

#ifdef KLMC
  ! wkjee - instead 'stop' simply return
  if(ioproc) then
    write(ioout,'(A)') "KLMC> in stopnow(return): something went wrong: checkout the .gout"
    call gflush(ioout)
  end if
  return
#endif

#else
  stop 'GULP terminated with an error'
#endif

  end

#ifdef ACCELRYS
  subroutine create_killfile ()
!=========================================================================C
! This creates an empty file called killfile in the current directory     C
!-------------------------------------------------------------------------C
! Arguments:                  none                                        C
!-------------------------------------------------------------------------C
!-------------------------------------------------------------------------C
! Written by Victor Milman, version 0.01, 01/05/02                        C
!=========================================================================C

  implicit none

  integer :: ios, lunit

  lunit = 60
  open(unit=lunit,form='FORMATTED',status='UNKNOWN', &
       access='SEQUENTIAL',file="killfile",iostat=ios)

  write (lunit,*) "Job terminated"
  close (lunit)

  return
  end subroutine create_killfile
#endif
