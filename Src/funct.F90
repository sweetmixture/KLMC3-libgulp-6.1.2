  subroutine funct(iflag,n,xc,fc,gc)
!
!  Supplies the function and derivatives
!
!   5/95 Modified to handle symmetrised second derivatives
!   6/95 Modified to allow for additive constraints
!   8/97 Modified to allow evaluation of cost function - SMW RIGB
!   9/97 Modified to allow evaluation 1st dervs of cost function
!   3/99 Option to write out derivatives to a file added
!  11/06 xc now passed getderv1
!   3/07 Intent added
!   5/08 lgeometryOK added as argument to xctox0 & iflag set to -2
!        on return if geometry is not OK.
!   5/09 Calls to x0tostr routines moved from energy to calling routine
!   7/09 Modified for empirical valence bond theory
!   1/11 Force minimisation added
!   5/11 Call to cutscheck added
!  10/13 Trap to prevent energy calculation with bad geometry added
!   9/16 Outderv, outdrv and outfrc now called for all processors
!   1/18 Trace added
!   5/20 strfin only called if lgrad1 is true
!  10/21 Output of CFG file added
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
  use control
  use gulp_files
  use general
  use genetic,     only : lgacost
  use parallel
  use symmetry
#ifdef KLMC
  use klmc,        only : lklmc_return_gulp
  use iochannels
#endif
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(inout) :: iflag
  integer(i4), intent(in)    :: n
  real(dp),    intent(out)   :: fc
  real(dp),    intent(in)    :: xc(*)
  real(dp),    intent(out)   :: gc(*)
!
!  Local variables
!
  integer(i4)                :: i
  logical                    :: lgeometryOK
  logical                    :: lgrad1
  logical                    :: lgrad2
  real(dp)                   :: g_cpu_time
  real(dp),             save :: tdmax = 0.0_dp
  real(dp)                   :: t1
  real(dp)                   :: t2
#ifdef TRACE
  call trace_in('funct')
#endif
!
  t1 = g_cpu_time()
  lgrad1 = (iflag.ge.1.or.lforcemin)
  lgrad2 = (iflag.ge.2)
!
  if (lgacost) then
!******************
!  Cost function  *
!******************
    call costfn(n,xc,fc,iflag)
    if (lgrad1) then
      call getderv1(n,xc,gc,lgrad2,.false.)
    endif
    return
  else
!***********
!  Energy  *
!***********
!
!  Convert optimisation variables to linear structure array
!
    call xctox0(n,xc,lgeometryOK)
    if (lgeometryOK) then
!
!  Convert linear structure array to main structure arrays
!
      if (lx0centroid) then
        call x0tostrcentroid
      else
        call x0tostr
      endif
!
!  Check core-shell distances
!
      call cutscheck
!
! KLMC TESTING BLOCK 10/24 wkjee
!#ifdef KLMC
!      ! 10/24 wkjee
!      ! return if core-shell check failed - TESTING BLOCK 08.10.2024
!      if (lklmc_return_gulp) then
!        if (ioproc) then
!          write(ioout,'(A)') " > KLMC_MESSAGE : abnormal stop requested : funct()"
!          call gflush(ioproc)
!        endif
!        return
!      endif
!#endif
!

!
      lfirst = .true.
!
!  Evaluate function and first derivatives
!
      call energy(fc,lgrad1,lgrad2)
!
!  For surface, get surface energy
!
      if (lseok) call surfaceenergy(fc)
!
!  Complete strain derivatives
!
      if (lgrad1.and.lstr) call strfin(lgrad2)
!
!  Output second derivatives 
!
      if (lgrad2) call outderv
!
!  Output energy and derivatives to a file
!
      if (ldrv) call outdrv(fc,lgrad1,lgrad2)
!
!  Output energy, forces, charges and/or charge differences to a file for Ovito/AtomEye
!
      if (lcfg) call outcfg(fc,lgrad1)
!
!  Option to write out a .frc file for QMPOT
!
      if (lfrc) call outfrc(fc,lgrad1,.false.)
!
!  First derivative handling
!
      if (lgrad1) then
        call getderv1(n,xc,gc,lgrad2,.false.)
      endif
!
!  Force minimisation option
!
      if (lforcemin) then
        fc = 0.0_dp
        if (n.gt.0) then
          do i = 1,n
            fc = fc + gc(i)**2
          enddo
          fc = sqrt(fc)/dble(n)
        endif
      endif
    else
!
!  Set dummy energy for a bad structure
!
      fc = 1.0d8
    endif
  endif
!
  t2 = g_cpu_time()
  if (t2-t1.gt.tdmax) tdmax = t2 - t1
  if (timmax.gt.0.0_dp) then
    if ((timmax-t2).lt.tdmax.and..not.lrelax) iflag = -1
  endif
  if (.not.lgeometryOK) iflag = -2
#ifdef TRACE
  call trace_out('funct')
#endif
!

#ifdef KLMC
  !
  ! 10/24 wkjee
  !
  if (lklmc_return_gulp) then
    if (ioproc) then
      write(ioout,'(A)') " > KLMC_MESSAGE : call stack : subroutine funct() > see funct.F90"
      call gflush(ioproc)
    endif
    return
  endif
!
! 10/24 wkjee
! KLMC developmenet purpose
! case: call stopnow() 
!
! funct() called at:
!
! functn.F90:    call funct(ifcall,n,xc,fc,gc)
! functn.F90:      call funct(ifcall,n,xc,fcf,gcf)
! functn.F90:      call funct(ifcall,n,xc,fcb,gcb)
! functn.F90:    call funct(ifcall,n,xc,fc,gc)
! gaexpd.F90:      call funct(iflag,nvar,xc,fc,gc)
! gaexpd.F90:    call funct(iflag,nvar,xc,fc,gc)
! gaexpd.F90:       call funct(iflag,nvar,xc,fc,gc)
! gaexpd.F90:       call funct(iflag,nvar,xc,fc,gc)
! gasubs.F90:    call funct(iflag,nvar,xc,fc,gcd)
! gasubs.F90:    call funct(iflag,nvar,xc,fc,gcd)
! gasubs.F90:    call funct(iflag,nvar,xc,fc,gcd)
! harmonicrelax.F90:    call funct(iflag,nvar,xc,fc,gc)
! harmonicrelax.F90:    call funct(iflag,nvar,xc,funct1,gc)
! lmbfgs.F90:    call funct(iflag,N,X,F,G)
! lmbfgs.F90:        call funct(iflag,N,X,F,G)
! minimise.F90:        call funct(iflag,nvar,xc,fc,gc,"min1")
! minimise.F90:              call funct(iflag,nvar,xc,funct1,gc,"min2")
! minimise.F90:        call funct(iflag,nvar,xc,fc,gc,"min3")
! numhess.F90:  call funct(iflag,nvar,xvar,funct2,gvar)
! olinmin.F90:    call funct(iflag,nvar,xparam,phi(2),grad,"olin1")
! olinmin.F90:    call funct(iflag,nvar,xparam,funct1,grad,"olin2")
! olinmin.F90:       call funct(iflag,nvar,xparam,funct1,grad,"olin3")
! optim.F90:          call funct(iflag,nvar,xc,fc,gc)
! optim.F90:            call funct(1_i4,nvar,xc,fc,gc)
! optim.F90:          call funct(iflag,nvar,xc,fc,gc)
! optim.F90:          call funct(iflag,nvar,xc,fc,gc)
! optim.F90:        call funct(iflag,nvar,xc,fc,gc)
! optim.F90:        call funct(iflag,nvar,xc,fc,gc)
! optim.F90:          call funct(iflag,nvar,xc,fc,gc)
! optim.F90:      call funct(iflag,nvar,xc,cost,gc)
! predict.F90:      call funct(iflag,nvar,xc,fc,gc)
! predict.F90:      call funct(iflag,nvar,xc,fc,gc)
!
! * call tree *
!
! optim.F90
!   > funct.F90
!   > functn.F90
!     > funct.F90
!   > minimise.F90
!     > funct.F90
!     > functn.F90
!       > funct.F90
!   > olinmin.F90
!     > funct.F90
!   > harmonicrealx.F90
!     > funct.F90
! ----
! funct.F90
!   # core-shell case only
!   > cutscheck.F90
!     > stopnow.F90
#endif
  return
  end
