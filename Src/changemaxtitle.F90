  subroutine changemaxtitle
!
!  Alters the size of the arrays associated with maxtitle
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
  use general
  use gulp_lengths
  use reallocate
#ifdef KLMC
  ! 07/23 wkjee
  ! reset internal counter
  use klmc, only : lklmc_maxtitle
#endif
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxtitle = 0
#ifdef KLMC
  if(lklmc_maxtitle) then
    ! set it as gulpdefault
    maxtitle = 20
    oldmaxtitle = 0
    lklmc_maxtitle = .false.
  end if
#endif  
!
!  Library data
!
  call realloc_ch(maxlinelength,titleword,maxtitle,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtitle','titleword')
!
!  Initialise defaults for new part of array
!
  if (maxtitle.gt.oldmaxtitle) then
    do i = oldmaxtitle+1,maxtitle
      call initmaxtitledefaults(i)
    enddo
  endif
!
!  Save current value of maxtitle for next call
!
  oldmaxtitle = maxtitle
!
  return
  end
