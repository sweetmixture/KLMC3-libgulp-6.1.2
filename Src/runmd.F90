  subroutine runmd
!
!  Main subroutine for Molecular Dynamics run control
!
!  MD variables:
!
!  nmdeq      = no. of time steps for equilibration
!  nmdprod    = no. of time steps for production run
!  nmdsamp    = no. of time steps between sampling properties
!  nmdwrite   = no. of time steps between writing MD dumpfile
!  tstep      = time step
!  velx/y/z   = velocity components in cartesian space
!  x2-5       = higher order derivatives for predictor-corrector
!  xalat      = cartesian coordinates (including cell shifts)
!  xclat      = cartesian coordinates of image in main cell
!
!  Units of velocity internally are Angs/ps
!  Common block /internl/ which normally contains the fractional
!  coordinates is overwritten during MD runs.
!
!   2/97 Modifications from JRH added for GC / Shell model MD
!   3/97 Breathing shell MD with opt at each step added
!   3/97 NVT ensemble added 
!   8/97 Scaling of C-S temperature during production added
!   6/01 Reset of averages only performed at end of equilibration
!        if there is to be production - this enables restarting
!        of the equilibration phase only.
!  11/02 Centre of mass now constrained through out equilibration.
!   2/04 Variable temperature for annealing introduced
!   6/04 Removal of linear and angular momentum now included in production
!  12/04 Order of deallocations reversed
!   1/05 Logical to compute interatomic vector table added in call
!        to mdfunct
!   5/06 Modified to handle case where restart is called starting from
!        from first step of production - call to mdreset now occurs
!  12/07 Unused variables removed
!  12/07 Target temperature now added as a variable
!   3/08 New MD integrator added
!   4/08 Call to mdfunct added before start of run to initialise forces 
!        and stress
!   8/08 Metadynamics introduced directly into MD routine
!   2/09 New argument added to mdprop call
!   2/09 Modified to accommodate new version of FoX and gulp_cml
!   7/11 Production flag added to call of mdpredict
!   8/11 Call to optshe and mdfunct modified to pass step number
!   9/11 Metadynamics internal code replaced with Plumed
!  12/15 Structure returned to configuration arrays at the end of run
!   3/16 Modified to work with Plumed2
!   3/16 Plumed finalisation added
!   1/18 Trace added
!  11/21 Modifications for TI added
!  12/21 Parallel modifications for TI added
!  12/21 ntinow moved to module
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
!  Julian Gale, CIC, Curtin University, December 2021
!
  use configurations, only : rvcfg, xcfg, ycfg, zcfg, radcfg, ndimen
  use control
  use current
  use derivatives,    only : maxd1
  use distances,      only : ndistancereset
  use g_constants,    only : kjmtoev
  use general
  use gulp_cml,       only : lcml, gulp_cml_startmodule, gulp_cml_endmodule
  use iochannels
  use mdlogic,        only : ladiabatic, lmd_inc_rot
  use moldyn
#ifdef PLUMED
  use m_plumed,       only : lplumed
  use m_plumed,       only : init_plumed, finalise_plumed
#endif
  use m_pr,           only : taub, taut, taubcfg, tautcfg
  use m_ti
  use parallel
  use shells
#ifdef ACCELRYS
  use license
#endif
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                               :: i
  integer(i4)                               :: it
  integer(i4)                               :: iter
  integer(i4)                               :: n
  integer(i4)                               :: nmde
  integer(i4)                               :: nmdp
  integer(i4)                               :: nmds
  integer(i4)                               :: nmdw
  integer(i4)                               :: nsint
  integer(i4)                               :: nstart
  integer(i4)                               :: status
  logical                                   :: loptsh
  logical                                   :: lprod
  logical                                   :: lreset
  logical                                   :: lsetmol
  logical                                   :: lwritecycle
  real(dp), dimension(:), allocatable       :: bspring
  real(dp)                                  :: g_cpu_time
  real(dp)                                  :: deltaG
  real(dp)                                  :: dUdl
  real(dp)                                  :: fc
  real(dp), dimension(:), allocatable       :: spring
  real(dp)                                  :: t1
  real(dp)                                  :: t2
  real(dp)                                  :: targettemperature
  real(dp),                            save :: tdmax = 0.0_dp
  real(dp)                                  :: tscale
  real(dp)                                  :: tstp
#ifdef TRACE
  call trace_in('runmd')
#endif
!
  tstp = tstep(ncf)
!
!  Allocate local memory
!
  allocate(spring(max(1,nshell)),stat=status)
  if (status/=0) call outofmemory('runmd','spring')
  allocate(bspring(max(1,nbsmat)),stat=status)
  if (status/=0) call outofmemory('runmd','bspring')
!
!  Check size of first derivative arrays
!
  if (numat.gt.maxd1) then
    maxd1 = numat
    call changemaxd1
  endif
!*******************************************
!  Setup parameters and initialise arrays  *
!*******************************************
  if (timesofar.gt.0.and.ncf.eq.ncfmdrun) then
!
!  Restart previous run
!
    call mdrestart
    nstart = nint(timesofar/tstp) + 1
  else
!
!  New run
!
    call setmd(.true.)
    timesofar = 0.0_dp
    nstart = 1
  endif
#ifdef PLUMED
!**************************
!  PLUMED initialisation  *
!**************************
  if (lplumed) then
    call init_plumed
  endif
#endif
!
!  Local variables
!
  nmde = nmdeq(ncf)
  nmdp = nmdprod(ncf)
  nmds = nmdsamp(ncf)
  nmdw = nmdwrite(ncf)
  taub = taubcfg(ncf)
  taut = tautcfg(ncf)
  tscale = tmdscale(ncf)
  nsint = max(nint(tmdscint(ncf)/tstp),1)
  tmdscint(ncf) = nsint*tstp
  loptsh = ((nshell.ne.0.and.ladiabatic).or.nbsmat.ne.0)
  lsetmol = .false.
  if ((nshell+nbsmat).ne.0) then
    call setspring(spring,bspring)
!
!  If adiabatic calculation optimise shell positions/radii before 
!  MD is started to avoid large forces in the beginning of the MD 
!  run
!
    if (ladiabatic) call optshe(0_i4,fc,iter,spring,bspring)
  endif
!
!  Output MD parameters on input
!
  if (ioproc) call outmdin
!
!  First call of mdfunct to initialise forces and stress
!
   call mdfunct(0_i4,fc,lsetmol,lreset,.true.)
!*******************
!  Equilibration  *
!*******************
  if (nmde.ge.nstart) then
    lprod = .false.
    if (ioproc) then
      write(ioout,'(/,''  Molecular dynamics equilibration :'',/)')
      if (lcml) call gulp_cml_startmodule(title='Molecular dynamics equilibration', dictRef='gulp:MDinit')
    endif
    t1 = g_cpu_time()
    iter = 0
    do i = nstart,nmde
      if (ioproc) call gflush(ioout)
      lwritecycle = (mod(i,nmdw).eq.0)
      call mdtemperature(i,targettemperature)
      call mdpredict(.false.)
      if (loptsh) then
        call optshe(i,fc,it,spring,bspring)
        iter = iter + it
      else
        lreset = (mod(i,ndistancereset(ncf)).eq.0)
        call mdfunct(i,fc,lsetmol,lreset,.true.)
        lsetmol = .false.
      endif
      call mdcorrect(.true.,.false.)
      call mdprop(fc,i,nmds,targettemperature,.false.)
      if (nmdintegrator.ne.4) then
        if (timesofar.le.tscale.and.mod(i,nsint).eq.0) call mdscale(lprod)
      endif 
      call mdvelcor(ndim.eq.0.and..not.lmd_inc_rot,.true.,.false.)
      timesofar = timesofar + tstp
      if (lwritecycle) then
        if (ioproc) call mdwrite(fc,lprod,i)
        t2 = g_cpu_time()
        tdmax = max(1.5_dp*(t2-t1),tdmax)
        if (timmax.gt.0.0_dp) then
          if ((timmax-t2).lt.tdmax) call mdstop
        endif
        t1 = t2
      endif
#ifdef ACCELRYS
      call sendHeartbeat(status)
      if (status /= 0) call gulpfinish
#endif
    enddo
    if (lcml.and.ioproc) call gulp_cml_endmodule()
    if (loptsh.and.ioproc) write(ioout,'(/,''  Average number of iterations to optimise shell positions : '',f5.2)') &
      dble(iter)/dble(nmde-nstart+1)
!
!  Reset properties after equilibration
!
    if (nmdp.gt.0) call mdreset
  elseif (nmde.eq.0.or.nstart.eq.nmde+1) then
!
!  Initialise properties for production
!
    if (nmdp.gt.0) call mdreset
  endif
!***************
!  Production  *
!***************
!
!  If this is a TI then turn on TI
!
  if (ltirun) then
    lti = .true.
  endif
  if (nmdp.gt.0) then
    lprod = .true.
    if (ioproc) then
      write(ioout,'(/,''  Molecular dynamics production :'',/)')
      if (lcml) call gulp_cml_startmodule(title='Molecular dynamics production', dictRef='gulp:MDprod')
    endif
    nstart = max((nstart-nmde),1)
    t1 = g_cpu_time()
    iter = 0
    do i = nstart,nmdp
!
!  Thermodynamic integration
!
      if (lti) then
!
!  Check whether lambda needs to change at this step
!
        ntinow = ntinow + 1
        if (ntinow.gt.ntistep) then
          ntinow = 1
          if (nprocs.gt.1) then
            call sumone(dUdlambda(nlambdanow),dUdl,"runmd","dUdlambda")
            dUdlambda(nlambdanow) = dUdl
          else
            dUdl = dUdlambda(nlambdanow)
          endif
          if (ioproc) then
            write(ioout,'(''  TI: lambda dUdlambda : '',f12.9,1x,f21.12,'' eV'')') lambda,dUdl
          endif
          if (nlambdanow.lt.nlambda+1) then
            nlambdanow = nlambdanow + 1
            lambda = lambda + lambda_step
          elseif (nlambdanow.eq.nlambda+1) then
            if (ioproc) then
              deltaG = 0.0_dp
              do n = 2,nlambda
                deltaG = deltaG + dUdlambda(n)
              enddo
              deltaG = deltaG + 0.5_dp*dUdlambda(1)
              deltaG = deltaG + 0.5_dp*dUdlambda(nlambda+1)
              deltaG = deltaG*lambda_step
              write(ioout,'(''  TI: deltaG = '',f18.9,'' eV'')') deltaG
              write(ioout,'(''  TI: deltaG = '',f18.6,'' kJ/mol'')') deltaG/kjmtoev
            endif
!
!  Turn off TI now that we are finished
!
            if (ltirun) then
              lti = .false.
            endif
          endif
          dUdlambda(nlambdanow) = 0.0_dp
        endif
      endif
!
!  Temperature handling
!
      call mdtemperature(i+nmde,targettemperature)
      lwritecycle = (mod(i,nmdw).eq.0)
      call mdpredict(.true.)
!
      if (loptsh) then
        call optshe(i,fc,it,spring,bspring)
        iter = iter + it
      else
        lreset = (mod(i,ndistancereset(ncf)).eq.0)
        call mdfunct(i,fc,lsetmol,lreset,.true.)
        lsetmol = .false.
      endif
      call mdcorrect(.false.,.true.)
      call mdprop(fc,i,nmds,targettemperature,.true.)
      if (nmdintegrator.ne.4) then
        if (timesofar.le.tscale.and.mod(i,nsint).eq.0) then
          call mdscale(lprod)
        elseif (.not.loptsh.and.nshell.gt.0) then
          call mdshscale
        endif
      endif
      call mdvelcor(ndim.eq.0.and..not.lmd_inc_rot,.false.,.true.)
      timesofar = timesofar + tstp
      if (lwritecycle) then
        if (ioproc) call mdwrite(fc,lprod,i)
        t2 = g_cpu_time()
        tdmax = max(1.5_dp*(t2-t1),tdmax)
        if (timmax.gt.0.0_dp) then
          if ((timmax-t2).lt.tdmax) call mdstop
        endif
        t1 = t2
      endif
#ifdef ACCELRYS
      call sendHeartbeat(status)
      if (status /= 0) call gulpfinish
#endif
    enddo
    if (lcml.and.ioproc) call gulp_cml_endmodule()
    if (loptsh.and.ioproc) write(ioout,'(/,''  Average number of iterations to optimise shell positions : '',f5.2)') &
      dble(iter)/dble(nmdp-nstart+1)
  endif
!
!  Output final results
!
  if (ioproc) call outmd
!
!  Return structure to configuration arrays
!
  if (ndim.eq.3) then
    do i = 1,3
      rvcfg(1,i,ncf) = rv(1,i)
      rvcfg(2,i,ncf) = rv(2,i)
      rvcfg(3,i,ncf) = rv(3,i)
    enddo
  elseif (ndim.eq.2) then
    do i = 1,2
      rvcfg(1,i,ncf) = rv(1,i)
      rvcfg(2,i,ncf) = rv(2,i)
    enddo
  elseif (ndim.eq.1) then
    rvcfg(1,1,ncf) = rv(1,1)
  endif
  do i = 1,numat
    radcfg(nsft+i) = rada(i)
    call cart2frac(ndimen(ncf),xclat(i),yclat(i),zclat(i),rv,xcfg(nsft+i),ycfg(nsft+i),zcfg(nsft+i),icosx(i),icosy(i),icosz(i))
  enddo
#ifdef PLUMED
!************************
!  PLUMED finalisation  *
!************************
  if (lplumed) then
    call finalise_plumed
  endif
#endif
!
!  Free local memory
!
  deallocate(bspring,stat=status)
  if (status/=0) call deallocate_error('runmd','bspring')
  deallocate(spring,stat=status)
  if (status/=0) call deallocate_error('runmd','spring')
#ifdef TRACE
  call trace_out('runmd')
#endif
!
  return
  end
