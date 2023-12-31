  subroutine changemaxspcellbo
!
!  Alters the size of the arrays associated with maxspcellbo
!
!   6/07 Structure of arrays changed to make atom and cell
!        arrays 1-D and a pointer to where the cell data
!        starts added
!   9/10 Initialisations now performed in a subroutine
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
!  Copyright Curtin University 2010
!
!  Julian Gale, CIC, Curtin University, September 2010
!
  use spatialbo
  use reallocate
#ifdef KLMC
  ! 07/23 wkjee
  ! reset internal counter
  use klmc, only : lklmc_maxspcellbo
#endif
!
  implicit none
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxspcellbo = 0
#ifdef KLMC
  if(lklmc_maxspcellbo) then
    ! set it as gulpdefault
    maxspcellbo = 0
    oldmaxspcellbo = 0
    lklmc_maxspcellbo = .false.
  end if
#endif  
!
!  Spatial decomposition data
!
  call realloc(nspcellatbo,maxspcellbo,ierror)
  if (ierror.ne.0) call outofmemory('changemaxspcellbo','nspcellatbo')
  call realloc(nspcellat1ptrbo,maxspcellbo,ierror)
  if (ierror.ne.0) call outofmemory('changemaxspcellbo','nspcellatptrbo')
!
!  Initialise defaults for new part of array
!
  if (maxspcellbo.gt.oldmaxspcellbo) then
    do i = oldmaxspcellbo+1,maxspcellbo
      call initmaxspcellbodefaults(i)
    enddo
  endif
!
!  Save current value of maxspcellbo for next call
!
  oldmaxspcellbo = maxspcellbo
!
  return
  end
