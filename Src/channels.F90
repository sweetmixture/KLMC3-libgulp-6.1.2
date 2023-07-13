  subroutine channels(nunit,lformatted)
!
!  Opens scratch file channels with a unique name
!  to avoid compilcations when two jobs are running
!  in the same directory. Has 9999 tries at finding
!  a different name - crude but simple! Anyone who
!  has more than 10000 scratch files clearly has
!  too big a disk quota and has never listed the
!  contents of their directroy!
!
!  nunit      = unit number to be opened as a temporary file
!  lformatted = if .true. then file will be opened for 
!               formatted i/o
!
!  wkjee - formatted means the formatted input from ... < input.gin > gulp.out ?
!
!   3/97 Created
!  10/97 Modified to handle input general channel number
!   9/03 Node number included in file name
!   5/06 Handling of negative unit numbers removed
!  11/06 Error in error format statement fixed
!  11/06 Call to ftow replaced with itow
!  11/07 Unused variables cleaned up
!   6/12 File handling changed to allow for spaces in names
!   8/13 Prefix option added
!   1/19 maxwordlength added
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, January 2019
!
  use gulp_lengths
  use iochannels
  use parallel
  implicit none
!
!  Passed variables
!
! wkjee - nunit -> iomod.F90 (module) iotmp passed from channel.F90 (subroutine)
! wkjee - lformatted -> passed from gulpsetup.F90 (subroutine), default = .true.
  integer(i4),     intent(in)  :: nunit
  logical,         intent(in)  :: lformatted
!
!  Local variables
!
  integer(i4)                  :: i
  integer(i4)                  :: ii
  integer(i4)                  :: ind
  integer(i4)                  :: indp
  integer(i4)                  :: inq
  integer(i4)                  :: n
  integer(i4)                  :: nfn
  logical                      :: lexist
  character(len=maxwordlength) :: filename
  character(len=5)             :: cnum
  character(len=5)             :: cproc
  character(len=5)             :: cunit
!
#ifdef MPI
  include 'mpif.h'
#endif
!********************
!  Build file name  *
!********************
  filename = ' '
!
!  Check for prefix length
!
  ! wkjee - prefix : from (iomod.F90: module iochannels)
  ! write(*,'(A,A)') "in channels: prefix (1): ", prefix
  ! wkjee - prefix : empty in the line above
  ! wkjee - nfn: integer(i4)
  nfn = index(prefix,' ') - 1
  ! wkjee - check
  ! write(*,'(A,A,I4)') "in channels: prefix (2) / nfn: : ", prefix, nfn
  ! wkjee - nfn = 0
  if (nfn.gt.0) then
    filename(1:nfn) = prefix(1:nfn)
  endif
  filename(nfn+1:nfn+8) = 'gulptmp_'
  ! wkjee - check filename
  ! write(*,'(A,A)') "in channels: filename (1): ", filename
  ! wkjee - filename = 'gulptmp_'
!
!  For parallel runs include node number
!
  ! wkjee - procid / nproc (module parallel in modules.F90)
  ! procid : Number for local processor
  ! nrpoc  : Number of processors for parallel execution
  ! cproc<character(len=5)>

  ! wkjee - nproc > 1, i.e., parallel
  ! wkjee - blow, do writing filename ... gulptmp_x_x ...
  ! wkjee
  ! write(*,'(A,I4,I4)') "in channel procid / nprocs : ", procid, nprocs

  if (nprocs.gt.1) then
    if (procid.gt.999) then
      write(cproc,'(i4)') procid
      filename(nfn+9:nfn+12) = cproc(1:4)
      ! wkjee
      ! write(*,'(A,A,A)') "in channel procid / filename (1): ", cproc, filename
      indp = nfn + 12
    elseif (procid.gt.99) then
      write(cproc,'(i3)') procid
      filename(nfn+9:nfn+11) = cproc(1:3)
      indp = nfn + 11
    elseif (procid.gt.9) then
      write(cproc,'(i2)') procid
      filename(nfn+9:nfn+10) = cproc(1:2)
      indp = nfn + 10
    else
      write(cproc,'(i1)') procid
      filename(nfn+9:nfn+9) = cproc(1:1)
      indp = nfn + 9
    endif
    indp = indp + 1
    filename(indp:indp) = '_'
  else
    indp = nfn + 8
  endif
  ! wkjee
  ! write(*,'(A,A,A)') "in channel procid / filename (2): ", cproc, filename
  ! wkjee
!
  if (nunit.gt.9999) then
    write(cunit,'(i5)') nunit
  elseif (nunit.gt.999) then
    write(cunit,'(i4)') nunit
  elseif (nunit.gt.99) then
    write(cunit,'(i3)') nunit
  elseif (nunit.gt.9) then
    write(cunit,'(i2)') nunit
  else
    write(cunit,'(i1)') nunit
  endif
  filename(indp+1:indp+5) = cunit
  ! wkjee
  ! write(*,'(A,A,A)') "in channel procid / filename (3): ", cproc, filename
  ! wkjee
  call endstring(filename,len(filename),ind)
  filename(ind:ind) = '_'
  do n = 1,9999
    call itow(cnum,n,5_i4)
    ii = index(cnum,'.')
    if (ii.gt.0) then
      do i = ii,5
        cnum(i:i) = ' '
      enddo
    endif
    filename(ind+1:ind+5) = cnum
    inquire(file=filename,exist=lexist,iostat=inq)
    ! wkjee
    ! write(*,'(A,A,A20,L4,I4)') "in channel procid / filename / lexist / inq (4): ", cproc, filename, lexist, inq
    ! wkjee - lexist = .false. / inq = 0 -> therefore goto 20 
    if (inq.ne.0) goto 10
    if (.not. lexist) goto 20
  enddo
!
!  If execution reaches this point no suitable filename has been found
!
  ! wkjee
  ! when did srun gulp < input.gin > gulp.out ... this line has been reached
  ! write(*,'(A,I)') "in channel procid: ", procid
  ! write(*,'(A)')   "in channel could not find unique filename"
  ! wkjee

  write(ioout,'(a)') ' Could not find unique filename'
  goto 999

10 write(ioout,'(a)') ' Inquire failed'
  goto 999

20 if (lformatted) then
    ! wkjee - lformatted<logical> expected to be .true.
    ! write(*,'(A,A20,L4)') "in channel filename / lformatted (1): ", filename, lformatted
    ! wkjee
    open(nunit,file=filename,form='formatted',status='new',err=998)

    ! wkjee - Q! if nunit is the file?

  else
    ! wkjee - lformatted<logical> expected to be .true.
    ! write(*,'(A,A20,L4)') "in channel filename / lformatted (2): ", filename, lformatted
    ! wkjee
    open(nunit,file=filename,form='unformatted',status='new',err=998)
  endif
  rewind(nunit)
  ! wkjee
  ! write(*,'(A,I4,I4)') "in channel before return, procid / nunit", procid, nunit
  ! wkjee
  return

998 if (ioproc) write(ioout,'(a)') 'error opening file '//filename

999 if (ioproc) then
    write(ioout,'(a,i5)') ' Program terminated in subroutine channels attempting to open unit ',nunit
  endif
  call stopnow('channels')
  end
