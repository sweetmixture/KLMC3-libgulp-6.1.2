  subroutine setpote
!
!  Setup tasks for twobody potentials
!
!   7/95 repcut added - limit range of exponential terms in
!        Buckingham potential to save expense
!   8/95 setup for 1/r**6 sum added
!   8/95 handling of repcut for lennard-jones added
!  12/00 apot now abs in recput calculation for exponential
!        forms to avoid negative log problems
!   4/05 Mods for cosh-spring potential added
!  11/05 Call to set taper parameters for Voter style added
!   8/07 Possibility of non-integer powers in L-J potentials added
!  11/07 Unused variables cleaned up
!   1/09 Arrays apot, bpot, cpot & dpot replaced by new 2-D array twopot
!   2/15 setlenn renamed to setcombine
!   6/19 nspecptr1/nspecptr2 added
!   9/19 Handling of unknown species added by adding species
!   6/21 Check on Baskes self potentials moved here from potword22
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
!  Julian Gale, CIC, Curtin University, June 2021
!
  use control
  use element
  use general
  use molecule
  use shells
  use species
  use two
  implicit none
!
!  Local variables
!
  integer(i4)        :: i
  integer(i4)        :: ni
  integer(i4)        :: nj
  integer(i4)        :: np
  integer(i4)        :: ns
  logical            :: lfound
  logical            :: lmatch
  real(dp)           :: ri
  real(dp)           :: rj
  real(dp)           :: rmpt
  real(dp)           :: rtol2
!
!  Set potential cutoffs for core-shell spring potential
!  This must be done after input is complete in case
!  cuts value is changed.
!  Also set maximum potential cutoff
!
  rpmax = 0.0_dp
  do i = 1,npote
    if (nptype(i).eq.8.or.nptype(i).eq.33) then
      rpot(i) = cuts
      rpot2(i) = 0.0_dp
    endif
    if (rpot(i).gt.rpmax) rpmax = rpot(i)
  enddo
  if (rpmax.gt.cutp) rpmax = cutp
!
!  Setup potentials with combination rules
!
  call setcombine
!
!  Set repulsive exponential cutoffs
!
  if (index(keyword,'norep').eq.0) then
    do i = 1,npote
      if (nptype(i).eq.1.or.nptype(i).eq.7) then
        if (abs(twopot(1,i)).lt.1.0d-6) then
          repcut(i) = rpot2(i)
        elseif (abs(twopot(2,i)).lt.1.0d-6) then
          repcut(i) = rpot(i)
        else
          repcut(i) = - twopot(2,i)*log((10.0**(-accuracy)/abs(twopot(1,i))))
          repcut(i) = min(repcut(i),rpot(i))
        endif
      elseif (nptype(i).eq.2) then
        if (abs(twopot(1,i)).lt.1.0d-6) then
          repcut(i) = rpot2(i)
        else
          rmpt = tpot(1,i)
          rmpt = 1.0_dp/rmpt
          repcut(i) = (abs(twopot(1,i)*10.0_dp**accuracy))**rmpt
          repcut(i) = min(repcut(i),rpot(i))
        endif
      else
        repcut(i) = rpot(i)
      endif
    enddo
  else
    do i = 1,npote
      repcut(i) = rpot(i)
    enddo
  endif
!
!  Setup C6 terms
!
  if (lc6) call setc6
!
!  Assign dummy cutoffs for two-body terms based on covalent
!  radii for bonded potentials. 
!
  rtol2 = 1.6_dp*rtol
  do i = 1,npote
    if (mmexc(i).eq.1) then
      ni = nspec1(i)
      nj = nspec2(i)
      if (ni.gt.maxele) ni = ni - maxele
      if (nj.gt.maxele) nj = nj - maxele
      ri = rcov(ni)
      rj = rcov(nj)
      rpot2(i) = 0.0_dp
      rpot(i) = rtol2*(ri + rj)
    endif
  enddo
!
!  Set up taper parameters for Voter style
!
  if (tapertype.eq.3) call settaper
!
!  Set up inverse rho values
!
  call rhoinv
!
!  Set pointers from each potential species to species list
!
  do i = 1,npote
!
!  Find species index for atom1
!
    lfound = .false.
    ns = 0
    do while (.not.lfound.and.ns.lt.nspec)
      ns = ns + 1
      lfound = lmatch(natspec(ns),ntypspec(ns),nspec1(i),nptyp1(i),.false.)   
    enddo
    if (lfound) then
      nspecptr1(i) = ns
    else
!
!  New species - add to list
!
      nspec = nspec + 1
      if (nspec.gt.maxspec) then
        maxspec = nspec + 20
        call changemaxspec
      endif
      natspec(nspec)  = nspec1(i)
      ntypspec(nspec) = nptyp1(i)
      qlspec(nspec) = 0.0_dp
      radspec(nspec) = 0.0_dp
      lqinspec(nspec) = .false.
    endif
!
!  Find species index for atom2
!
    lfound = .false.
    ns = 0
    do while (.not.lfound.and.ns.lt.nspec)
      ns = ns + 1
      lfound = lmatch(natspec(ns),ntypspec(ns),nspec2(i),nptyp2(i),.false.)
    enddo
    if (lfound) then
      nspecptr2(i) = ns
    else
!
!  New species - add to list
!
      nspec = nspec + 1
      if (nspec.gt.maxspec) then
        maxspec = nspec + 20
        call changemaxspec
      endif
      natspec(nspec)  = nspec2(i)
      ntypspec(nspec) = nptyp2(i)
      qlspec(nspec) = 0.0_dp
      radspec(nspec) = 0.0_dp
      lqinspec(nspec) = .false.
    endif
  enddo
!
!  Setup Baskes self term pointers if needed
!
  do i = 1,npote
    if (nptype(i).eq.45.or.nptype(i).eq.55) then
!
!  Search for self potentials
!
      lfound = .false.
      np = 0
      do while (.not.lfound.and.np.lt.npote)
        np = np + 1
        if (nptype(np).eq.45.or.nptype(np).eq.55) then
          lfound = (lmatch(nspec1(i),nptyp1(i),nspec1(np),nptyp1(np),.true.).and. &
                    lmatch(nspec1(i),nptyp1(i),nspec2(np),nptyp2(np),.true.))
        endif
      enddo
      if (lfound) then
        ipot(3,i) = np
      else
        call outerror('Baskes potential has no matching self term',0_i4)
        call stopnow('setpote')
      endif
!
      lfound = .false.
      np = 0
      do while (.not.lfound.and.np.lt.npote)
        np = np + 1
        if (nptype(np).eq.45.or.nptype(np).eq.55) then
          lfound = (lmatch(nspec2(i),nptyp2(i),nspec1(np),nptyp1(np),.true.).and. &
                    lmatch(nspec2(i),nptyp2(i),nspec2(np),nptyp2(np),.true.))
        endif
      enddo
      if (lfound) then
        ipot(4,i) = np
      else
        call outerror('Baskes potential has no matching self term',0_i4)
        call stopnow('setpote')
      endif
    endif
  enddo
!
  return
  end
