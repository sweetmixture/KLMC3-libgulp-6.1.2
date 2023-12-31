  subroutine changemaxmctrans
!
!  Alters the size of the arrays associated with maxmctrans
!
!  10/21 Created from changemaxmctrans
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
  use montecarlo
  use reallocate
#ifdef KLMC
  ! 07/23 wkjee
  ! reset internal counter
  use klmc, only : lklmc_maxmctrans
#endif
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxmctrans = 0
#ifdef KLMC
  if(lklmc_maxmctrans) then
    ! set it as gulpdefault
    maxmctrans = 2
    oldmaxmctrans = 0
    lklmc_maxmctrans = .false.
  end if
#endif  
!
  call realloc(nmctrannat,2_i4,maxmctrans,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmctrans','nmctrannat')
  call realloc(nmctranspec,2_i4,maxmctrans,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmctrans','nmctranspec')
  call realloc(nmctrantype,2_i4,maxmctrans,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmctrans','nmctrantype')
  call realloc(ntranable,maxmctrans,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmctrans','ntranable')
  call realloc(ptran,maxmctrans,ierror)
  if (ierror.ne.0) call outofmemory('changemaxmctrans','pswap')
!
!  Initialise new parts of data arrays
!
  if (maxmctrans.gt.oldmaxmctrans) then
    do i = oldmaxmctrans+1,maxmctrans
      call initmaxmctransdefaults(i)
    enddo
  endif
!
!  Save current value of maxmctrans for next call
!
  oldmaxmctrans = maxmctrans
!
  return
  end
