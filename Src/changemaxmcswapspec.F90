  subroutine changemaxmcswapspec
!
!  Alters the size of the arrays associated with maxmcswapspec
!
!   1/09 Created
!   9/10 Initialisations now performed in a subroutine
!   5/16 Multiple mcswaps added
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
!  Copyright Curtin University 2016
!
!  Julian Gale, CIC, Curtin University, May 2016
!
  use montecarlo
  use reallocate
#ifdef KLMC
  ! 07/23 wkjee
  ! reset internal counter
  use klmc, only : lklmc_maxmcswapspec
#endif
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxmcswapspec = 0
#ifdef KLMC
  if(lklmc_maxmcswapspec) then
    ! set it as gulpdefault
    maxmcswapspec = 1
    oldmaxmcswapspec = 0
    lklmc_maxmcswapspec = .false.
  end if
#endif  
!
  call realloc(nmcswapnat,maxmcswapspec,maxmcswaps,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmcswapspec','nmcswapnat')
  call realloc(nmcswaptype,maxmcswapspec,maxmcswaps,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmcswapspec','nmcswaptype')
!
!  Initialise new parts of data arrays
!
  if (maxmcswapspec.gt.oldmaxmcswapspec) then
    do i = oldmaxmcswapspec+1,maxmcswapspec
      call initmaxmcswapspecdefaults(i)
    enddo
  endif
!
!  Save current value of maxmcswapspec for next call
!
  oldmaxmcswapspec = maxmcswapspec
!
  return
  end
