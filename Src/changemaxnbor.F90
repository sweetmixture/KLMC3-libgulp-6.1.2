  subroutine changemaxnboR
!
!  Alters the size of the arrays associated with maxnboR
! 
!   6/10 Option to for second atom specifier added to repulsive terms.
!   9/10 Initialisations now performed in a subroutine
!   1/14 lBOzrlR added
!   9/15 BOccoeffR made into a 2-D array
!  11/20 Tersoff reorganised
!   6/21 lBOgik arrays added
!   9/21 Symmetric flags added for bond order potentials
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
!  Julian Gale, CIC, Curtin University, September 2021
!
  use bondorderdata
  use reallocate
#ifdef KLMC
  ! 07/23 wkjee
  ! reset internal counter
  use klmc, only : lklmc_maxnbor
#endif
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxnboR = 0
#ifdef KLMC
  if(lklmc_maxnbor) then
    ! set it as gulpdefault
    maxnboR = 1
    oldmaxnboR = 0
    lklmc_maxnbor = .false.
  end if
#endif  
!
!  Bond order potential data
!
  call realloc(BOccoeffR,5_i4,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','BOccoeffR')
  call realloc(BOhcoeffR,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','BOhcoeffR')
  call realloc(BOlcoeffR,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','BOlcoeffR')
  call realloc(BOmcoeffR,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','BOmcoeffR')
  call realloc(BOocoeffR,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','BOocoeffR')
  call realloc(nBOspecR1,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','nBOspecR1')
  call realloc(nBOspecR2,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','nBOspecR2')
  call realloc(nBOtypR1,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','nBOtypR1')
  call realloc(nBOtypR2,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','nBOtypR2')
  call realloc(nBOtypeR,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','nBOtypeR')
  call realloc(lBOgikR,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','lBOgikR')
  call realloc(lBOsymR,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','lBOsymR')
  call realloc(lBOzrlR,maxnboR,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnboR','lBOzrlR')
!
!  Initialise defaults for new part of array
!
  if (maxnboR.gt.oldmaxnboR) then
    do i = oldmaxnboR+1,maxnboR
      call initmaxnbordefaults(i)
    enddo
  endif
!
!  Save current value of maxnboR for next call
!
  oldmaxnboR = maxnboR
!
  return
  end
