  subroutine changemaxlib
!
!  Alters the size of the arrays associated with maxlib
!
!   9/10 Initialisations now performed in a subroutine
!   1/19 Maxwordlength changes
!   1/19 Use of general string reallocate added
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
  use library
  use reallocate
#ifdef KLMC
  ! 07/23 wkjee
  ! reset internal counter
  use klmc, only : lklmc_maxlib
#endif
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxlib = 0
#ifdef KLMC
  if(lklmc_maxlib) then
    ! set it as gulpdefault
    maxlib = 4
    oldmaxlib = 0
    lklmc_maxlib = .false.
  end if
#endif  
!
!  Library data
!
  call realloc_ch(maxwordlength,libname,maxlib,ierror)
  if (ierror.ne.0) call outofmemory('changemaxlib','libname')
!
!  Initialise defaults for new part of array
!
  if (maxlib.gt.oldmaxlib) then
    do i = oldmaxlib+1,maxlib
      call initmaxlibdefaults(i)
    enddo
  endif
!
!  Save current value of maxlib for next call
!
  oldmaxlib = maxlib
!
  return
  end
