  subroutine changemaxvalener
!
!  Alters the size of the arrays associated with maxvalener
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
!  Julian Gale, CIC, Curtin University, July 2021
!
  use bondvalence
  use reallocate
  implicit none
!
  integer(i4)       :: ierror
!
  call realloc(nVBspecE,maxvalener,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalener','nVBspecE')
  call realloc(nVBtypeE,maxvalener,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalener','nVBtypeE')
  call realloc(nvalenertype,maxvalener,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalener','nvalenertype')
  call realloc(VBparE,maxVBparE,maxvalener,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalener','VBparE')
!
  return
  end
