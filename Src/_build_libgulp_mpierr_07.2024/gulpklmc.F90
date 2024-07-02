!program gulpklmc_program
!
!end program gulpklmc_program ! this is not working ont Cray HPE

subroutine gulpklmc( MPI_comm_klmc, klmc_task_iopath, klmc_task_id, klmc_worker_id ) bind (C,name="gulpklmc")
  !
  ! wkjee - KLMC development 07/2023
  !
  use datatypes,  only : i4
  use parallel,   only : ioproc
  use klmc
  use iso_c_binding
  use mpi
  use iochannels
  ! able to use 'gulp_klmc_iopath' module from modules.F90 (added)
  ! include 'mpif.h'
  ! able to use 'gulp_klmc_iopath' module from modules.F90 (added)
  implicit none


  ! MPI_comm_klmc				: communicator passed by KLMC taskfarm
  ! klmc_task_iopath				: path to read/write gulpt input/output(standard *gout)
  ! klmc_task_id				: task_id
  ! klmc_worker_id				: taskfarm worker id

  integer*4,          intent(in) :: MPI_comm_klmc
  integer,            intent(in) :: klmc_task_id
  integer,            intent(in) :: klmc_worker_id
  ! character(kind=c_char,len=512), intent(in) :: klmc_task_iopath
  ! character(kind=c_char,len=512), intent(in) :: klmc_task_iopath
  ! character(kind=c_char), intent(in) :: klmc_task_iopath(512)
  ! character(kind=c_char), dimension(512), intent(in) :: klmc_task_iopath
  ! wkjee - generic string C->Fortran passer update
  character(kind=c_char), dimension(*), intent(in) :: klmc_task_iopath
  character(len=512) :: klmc_task_iopath_dummy
  
  ! iopath_length - string length of the file path
  integer iopath_length
  ! integer, value :: iopath_length
  integer ierr
  integer rank
  integer cpu_count

  integer(i4)        :: iret
  integer(i4)        :: klmc_link
  ! klmc_link .eq. ichemsh_link
  ! integer(i4)      :: ichemsh_link
  integer*4          :: MPI_comm_togulp
  
  ! timing
  character(len=40) :: start_timestamp, end_timestamp
  double precision  :: mpi_tstart, mpi_tend, mpi_elapsed_t

  ! wkjee - do not touch this block . related to gulp default ============
  iret = 0
  ! wkjee - standard  gulp mode ichemsh_link(klmc_link) = 0: using gulp standalone mode
  ! wkjee - chemshell      mode ichemsh_link(klmc_link) = 1: using chain gulp runs: without memory nullifying 
  klmc_link = 0

  !if ( rank .eq. 0 ) then
  !if (ioproc) then
  !  write(*,*) " * [1] passing klmc_comm"
  !endif
  call MPI_Barrier(MPI_comm_klmc,iret)

  MPI_comm_togulp = MPI_comm_klmc
  ! end wkjee ============================================================

  !if ( rank .eq. 0 ) then
  !if (ioproc) then
  !  write(*,*) " * [2] get rank"
  !endif
  call MPI_Barrier(MPI_comm_klmc,iret)

  call MPI_Comm_rank(MPI_comm_togulp,rank,ierr)
  !if ( rank .eq. 0 ) then
  !if (ioproc) then
  !  write(*,*) " * [3] get size"
  !endif
  call MPI_Barrier(MPI_comm_klmc,iret)

  call MPI_Comm_size(MPI_comm_togulp,cpu_count,ierr)

  ! wkjee - setup iopath (possibly a message sent from the taskfarm <master>)
  ! from module klmc in modules.F90 (GULP main source), attribute<save>

  ! passing the working directory - the directory must include 'gulp_klmc.gin' (standard gulp input file)
  ! klmc_task_iopath = "/work/e05/e05/wkjee/Software/gulp-6.1.2/Src/Custom/path_test"
  ! gulp_klmc_iopath = ""
  ! gulp_klmc_iopath = klmc_task_iopath(1:60)
  ! gulp_klmc_iopath = klmc_task_iopath
  ! gulp_klmc_iopath = ""

  ! wkjee - solid working without bulshit ending character
  ! <IMPORTANT> iopath_length - remove the terminating zero '0' at the end 

  ! wkjee - refactoring klmc input path
  klmc_task_iopath_dummy = ""
  iopath_length = 0
  do
     if(klmc_task_iopath(iopath_length+1) == C_NULL_CHAR) exit
     iopath_length = iopath_length + 1
     klmc_task_iopath_dummy(iopath_length:iopath_length) = klmc_task_iopath(iopath_length)
  end do
  ! iopath_length = len_trim(klmc_task_iopath)-1
  gulp_klmc_iopath = klmc_task_iopath_dummy(1:iopath_length)
  ! write(*,'(A,I4)') "in gulpklmc: klmc_task_iopath_length: ", len_trim(klmc_task_iopath)-1
  !write(*,'(A,A)') "in gulpklmc.F90 : KLMC klmc_task_iopath: ", klmc_task_iopath_dummy
  !write(*,'(A,A)') "in gulpklmc.F90 : KLMC gulp_klmc_iopath: ", gulp_klmc_iopath 

  ! passing taskfarm info
  gulp_klmc_task_id = klmc_task_id
  gulp_klmc_worker_id = klmc_worker_id

  !if ( rank .eq. 0 ) then
  !if (ioproc) then
  !  write(*,*) " * [4] get wtime"
  !endif
  call MPI_Barrier(MPI_comm_klmc,iret)

  ! wallclock time measure
  mpi_tstart = MPI_Wtime()
  start_timestamp = getCurrentDateTime()

!======================
! Launch GULP-KLMC 
!======================

  !if ( rank .eq. 0 ) then
  !if (ioproc) then
  !  write(ioout,*) " * [5] call initmax"
  !endif
  call MPI_Barrier(MPI_comm_klmc,iret)

  call gulpklmc_initmax

! write(*,'(A,L)') "in gulpklmc.F90: lklmcfresrun before gulpmain() : ", lklmcfreshrun

  !if ( rank .eq. 0 ) then
  !if (ioproc) then
  !  write(ioout,*) " * [6] call gulpmain"
  !endif
  call MPI_Barrier(MPI_comm_klmc,iret)

  call gulpmain(iret, klmc_link, MPI_comm_togulp)

! write(*,'(A,L)') "in gulpklmc.F90: lklmcfresrun after  gulpmain() : ", lklmcfreshrun

!======================
! Finish GULP-KLMC run 
!======================

! call reinitialise

  !if ( rank .eq. 0 ) then
  if (ioproc) then
    write(ioout,*) " * [1] call gulpfinish"
  endif
  call MPI_Barrier(MPI_comm_klmc,iret)

  call gulpfinish

  !if ( rank .eq. 0 ) then
  if (ioproc) then
    write(ioout,*) " * [2] call dealloc_all"
  endif
  call MPI_Barrier(MPI_comm_klmc,iret)

  call gulpklmc_deallocate_all
  
!=====================
! Print out log
!=====================
  mpi_tend = MPI_Wtime()
  end_timestamp = getCurrentDateTime()
  mpi_elapsed_t = mpi_tend - mpi_tstart

  !if ( rank .eq. 0 ) then
  if (ioproc) then
    write(ioout,'(A,A,A,A,A,F12.8,A,I4,A,I4,A,I4,A,A)') "KLMC_FINALIZE> start: ", trim(adjustl(start_timestamp)), &
    " end: ", trim(adjustl(end_timestamp)), " elapsed_time: ", mpi_elapsed_t, " cpu_count: ", cpu_count, &
    " task_id: ", klmc_task_id, " worker_id: ", klmc_worker_id, " io_path: ", gulp_klmc_iopath
  end if
  call MPI_Barrier(MPI_comm_klmc,iret)

  !if ( rank .eq. 0 ) then
  if (ioproc) then
    write(ioout,*) " * [3] JUST BEFORE RETURN !!!"
  endif
  call MPI_Barrier(MPI_comm_klmc,iret)

  return

! ---------------------------------------------------------------------------
  contains
! ---------------------------------------------------------------------------
    function getCurrentDateTime() result(dateTimeStr)
      character(len=40) :: dateTimeStr
      integer :: ierr
      character(len=19) :: timeStr
      character(len=21) :: dateStr
  
      call date_and_time(date=dateStr, time=timeStr)
      dateTimeStr = trim(dateStr) // ' ' // trim(timeStr)
  
      return
    end function getCurrentDateTime
  
    subroutine sflag( rank, flag, task_id )
      integer, intent(in) :: rank,flag,task_id
      if( rank == 0 ) then
        write(6,'(a,I4,I4)') "F> flag task_id: ", flag, task_id
      end if
    end subroutine

end subroutine gulpklmc
