  subroutine changemaxplanepot
!
!  Alters the size of the arrays associated with maxplanepot
!
!   4/07 Created from changemaxUFFspec
!   9/10 Initialisations now performed in a subroutine
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
  use library
  use reallocate
  use plane
#ifdef KLMC
  ! 07/23 wkjee
  ! reset internal counter
  use klmc, only : lklmc_maxplanepot
#endif
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxplanepot = 0
#ifdef KLMC
  if(lklmc_maxplanepot) then
    ! set it as gulpdefault
    maxplanepot = 0
    oldmaxplanepot = 0
    lklmc_maxplanepot = .false.
  end if
#endif  
!
!  Species data
!
  call realloc_ch(5_i4,planepotsymbol,maxplanepot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxplanepot','planepotsymbol')
  call realloc(natplanepot,maxplanepot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxplanepot','natplanepot')
  call realloc(ntypplanepot,maxplanepot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxplanepot','ntypplanepot')
  call realloc(nplanepotpower,2_i4,maxplanepot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxplanepot','nplanepotpower')
  call realloc(nplanepottype,maxplanepot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxplanepot','nplanepottype')
  call realloc(planepot,3_i4,maxplanepot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxplanepot','planepot')
  call realloc(planepotrmin,maxplanepot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxplanepot','planepotrmin')
  call realloc(planepotrmax,maxplanepot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxplanepot','planepotrmax')
!
!  Initialise new parts of data arrays
!
  if (maxplanepot.gt.oldmaxplanepot) then
    do i = oldmaxplanepot+1,maxplanepot
      call initmaxplanepotdefaults(i)
    enddo
  endif
!
!  Save current value of maxplanepot for next call
!
  oldmaxplanepot = maxplanepot
!
  return
  end
