  subroutine changemaxbondQ
!
!  Alters the size of the arrays associated with maxbondQ
!
!   5/07 Created from changemaxnboa.f90
!   9/10 Initialisations now performed in a subroutin
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
  use bondcharge
  use reallocate
#ifdef KLMC
  ! 07/23 wkjee
  ! reset internal counter
  use klmc, only : lklmc_maxbondq
#endif
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxbondQ = 0
#ifdef KLMC
  if(lklmc_maxbondq) then
  ! set it as gulpdefault
    maxbondQ = 1
    oldmaxbondQ = 0
    lklmc_maxbondq = .false.
  end if
#endif
!
!  Bond order potential data
!
  call realloc_ch(5_i4,symbolbondQ,2_i4,maxbondQ,ierror)
  if (ierror.ne.0) call outofmemory('changemaxbondQ','symbolbondQ')
  call realloc(nbondQspec1,maxbondQ,ierror)
  if (ierror.ne.0) call outofmemory('changemaxbondQ','nbondQspec1')
  call realloc(nbondQspec2,maxbondQ,ierror)
  if (ierror.ne.0) call outofmemory('changemaxbondQ','nbondQspec2')
  call realloc(nbondQtyp1,maxbondQ,ierror)
  if (ierror.ne.0) call outofmemory('changemaxbondQ','nbondQtyp1')
  call realloc(nbondQtyp2,maxbondQ,ierror)
  if (ierror.ne.0) call outofmemory('changemaxbondQ','nbondQtyp2')
  call realloc(bondQincrement,maxbondQ,ierror)
  if (ierror.ne.0) call outofmemory('changemaxbondQ','bondQincrement')
!
!  Initialise defaults for new part of array
!
  if (maxbondQ.gt.oldmaxbondQ) then
    do i = oldmaxbondQ+1,maxbondQ
      call initmaxbondqdefaults(i)
    enddo
  endif
!
!  Save current value of maxbondQ for next call
!
  oldmaxbondQ = maxbondQ
!
  return
  end
