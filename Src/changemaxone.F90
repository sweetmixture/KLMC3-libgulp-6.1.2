  subroutine changemaxone
!
!  Alters the size of the arrays associated with maxone
!
!   1/10 Created from changemaxpot
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
  use reallocate
  use one
#ifdef KLMC
  ! 07/23 wkjee
  ! reset internal counter
  use klmc, only : lklmc_maxone
#endif
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxone = 0
#ifdef KLMC
  if(lklmc_maxone) then
    ! set it as gulpdefault
    maxone = 10
    oldmaxone = 0
    lklmc_maxone = .false.
  end if
#endif  
!
!  Potential data
!
  call realloc(onepot,maxone,ierror)
  if (ierror.ne.0) call outofmemory('changemaxone','onepot')
  call realloc(nspec11,maxone,ierror)
  if (ierror.ne.0) call outofmemory('changemaxone','nspec11')
  call realloc(nptyp11,maxone,ierror)
  if (ierror.ne.0) call outofmemory('changemaxone','nptyp11')
  call realloc_ch(5_i4,symbol1,maxone,ierror)
  if (ierror.ne.0) call outofmemory('changemaxone','symbol1')
!
!  Initialise defaults for new part of array
!
  if (maxone.gt.oldmaxone) then
    do i = oldmaxone+1,maxone
      call init1bodydefaults(i)
    enddo
  endif
!
!  Save current value of maxone for next call
!
  oldmaxone = maxone
!
  return
  end
