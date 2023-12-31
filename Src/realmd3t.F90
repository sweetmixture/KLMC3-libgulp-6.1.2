  subroutine realmd3t(eatom,ereal,erecip,ec6,eqeq,esregion12,esregion2,eattach,lgrad1)
!
!  Subroutine for calculating real space energy and gradients using pretabulated 
!  interatomic vectors.
!
!  lfreeze = if .true. only do derivatives for atoms with a degree of freedom
!
!   1/05 Created from realmd3
!   4/05 Mods for cosh-spring added
!   7/05 Streitz and Mintmire modifications added
!   3/07 Printing of twobody energies added as an option
!   5/07 QM/MM scheme added
!   5/07 Argument list for twobody call modified
!  11/07 Unused variables cleaned up
!  11/07 Modifications for reaxFF Coulomb term added
!   1/08 lreaxFFqreal removed
!   3/08 ReaxFF charge scheme removed as this is handled elsewhere
!  12/08 scrho, and defect versions, converted to 2-D array for MEAM
!   1/09 Arrays apot, bpot, cpot & dpot replaced by new 2-D array twopot
!   2/09 Expanding and contracting of maxdis arrays removed
!   3/09 lorder12loc added to twobodymd argument list
!   4/09 MEAM density modifications removed since these are now handled in 
!        a separate routine.
!   6/09 Site energy and virials added
!  11/09 Region derivatives added
!  11/11 Region-region energy contributions stored
!   4/12 Explicit virial calculation removed as no longer needed
!   4/12 xvir, yvir and zvir removed
!   5/12 Atomic stresses added
!   7/13 Call to twobody modified to include core-shell occupancy factor
!   3/14 lorder12loc removed twobody argument list
!  12/14 rtrm1 changed from scalar to array
!   2/15 MM3buck added
!   2/18 Trace added
!   5/18 Modified for q ranges
!   8/18 Call to twostrterms introduced for setting of rpd
!   9/18 Strain module introduced
!   5/19 Skipping of ij loop now checks for polarisation
!   8/19 Short range damping of polarisation added
!   8/19 Trap for neemrptr being zero added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  11/19 Rigid molecule modifications added
!   3/20 Electrostatic cutoff only included where the either charge exceeds the threshold
!   4/20 derv3c and d2r2dsdc added for benefit of rigid molecules
!   4/20 derv3c changes reversed as they are no longer required
!   7/20 Modifications for gfortran v10
!  11/21 Species added to twobody calls
!  12/21 pGFNFF modifications added
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
  use configurations, only : lbsmat,lsliceatom,nregions,nregionno, nregiontype, QMMMmode
  use g_constants
  use control
  use current
  use datatypes
  use derivatives
  use distances,      only : ndistanceind, ndistanceij
  use eam,            only : lMEAMden
  use eemdata
  use element
  use energies,       only : siteenergy, eregion2region
  use general,        only : cutw
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_rad
  use iochannels,     only : ioout
  use kspace
  use m_strain,       only : twostrterms
  use mdlogic
  use molecule
  use optimisation
  use parallel
  use polarise
  use potentialxyz
  use realvectors
  use shells
  use sutton
  use symmetry
  use thresholds,     only : thresh_q
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  logical,     intent(in)                      :: lgrad1
  real(dp),    intent(inout)                   :: eatom
  real(dp),    intent(inout)                   :: eattach
  real(dp),    intent(inout)                   :: ec6
  real(dp),    intent(inout)                   :: eqeq
  real(dp),    intent(inout)                   :: ereal
  real(dp),    intent(inout)                   :: erecip
  real(dp),    intent(inout)                   :: esregion12
  real(dp),    intent(inout)                   :: esregion2
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ilast
  integer(i4)                                  :: ind
  integer(i4)                                  :: ix
  integer(i4)                                  :: j
  integer(i4)                                  :: jx
  integer(i4)                                  :: k
  integer(i4)                                  :: kl
  integer(i4)                                  :: ks
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: natjlast
  integer(i4)                                  :: nmolonly
  integer(i4)                                  :: nor
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: npots
  integer(i4)                                  :: nqri
  integer(i4)                                  :: nqrj
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: nspeci
  integer(i4)                                  :: nspecj
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: ntypjlast
  integer(i4)                                  :: status
  logical                                      :: lattach
  logical                                      :: lc6loc
  logical                                      :: lcspair
  logical                                      :: lewaldtype
  logical                                      :: lmatch
  logical                                      :: lmdl
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lorder12loc
  logical                                      :: lQMMMelectro
  logical                                      :: lreg2one
  logical                                      :: lreg2pair
  logical                                      :: lself
  logical                                      :: lsg1  
  logical                                      :: lslicei
  logical                                      :: lslicej
  real(dp)                                     :: c6tot
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2e
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2q
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2s
  real(dp)                                     :: cut2w
  real(dp)                                     :: d2self
  real(dp)                                     :: derive0self
  real(dp)                                     :: dqi(1)
  real(dp)                                     :: dqj(1)
  real(dp)                                     :: eatom2
  real(dp)                                     :: ec62 
  real(dp)                                     :: eqeq2  
  real(dp)                                     :: ereal2
  real(dp)                                     :: erecip2
  real(dp)                                     :: esum 
  real(dp)                                     :: factor
  real(dp)                                     :: fct
  real(dp)                                     :: half
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ospfct
  real(dp)                                     :: polfct
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: radi
  real(dp)                                     :: radj
  real(dp)                                     :: radsum
  real(dp)                                     :: rp
  real(dp)                                     :: rq
  real(dp)                                     :: rq2
  real(dp)                                     :: rstrdloc(6)
  real(dp)                                     :: sctrm1
  real(dp)                                     :: sctrm2
  real(dp),    dimension(:), allocatable       :: sum
  real(dp),    dimension(:), allocatable       :: sum2
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: tsum0
  real(dp)                                     :: tsuml
  real(dp)                                     :: zetaij
#ifdef TRACE
  call trace_in('realmd3t')
#endif
!
  time1 = g_cpu_time()
!
!  Local variables
!
  lmdl = lmd
  if (.not.lgrad1) lmdl = .false.
  lc6loc = (lc6.and.ndim.eq.3)
!
!  Set the Coulomb term type based on dimensionality :
!
!  1-D         => 1/r
!  2-D and 3-D => erfc(seta*r)/r
!
  lewaldtype = (ndim.ne.1.or.lwolf)
!
  tsuml = 0.0_dp
  lsg1 = (lgrad1.and.lstr)
!
!  Set up cutoffs
!
  cut2p = cutp*cutp
  cut2e = rmx2
  cut2s = cuts*cuts
  cut2w = cutw*cutw
  if (lwolf) then
    cut2q = cut2w
  else
    cut2q = cut2e
  endif
!
!  Set radius for Coulomb correction
!
  if (lgfnff) then
    rq = gfnff_eeq_rad
  else
    rq = rqeq
  endif
  rq2 = rq**2
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('realmd3t','npotl')
  allocate(sum(numat),stat=status)
  if (status/=0) call outofmemory('realmd3t','sum')
  allocate(sum2(numat),stat=status)
  if (status/=0) call outofmemory('realmd3t','sum2')
!
  if (.not.lnoreal) then
!
!  Opening banner for energy decomposition
!
    if (lPrintTwo) then
      call mpbarrier
      if (ioproc) then
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''  Two : Atom No. 1  Atom No. 2    Short-range energy (eV)   Coulomb energy (eV) '')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    endif
!***************************************************************
!  Atomistic and real space electrostatic component of energy  *
!***************************************************************
    ilast = - 1
    natjlast = - 1
    ntypjlast = - 1
!
!  Outer loop over sites
!
    do i = procid+1,numat,nprocs
!
!  Inner loop over second site
!
      nati = nat(i)
      ntypi = nftype(i)
      nspeci = nspecptr(i)
      qli = qf(i)
      oci = occuf(i)
      nregioni = nregionno(nsft+nrelf2a(i))
      nregiontypi = nregiontype(nregioni,ncf)
      lopi = (.not.lfreeze.or.lopf(nrelf2a(i)))
      lslicei = lsliceatom(nsft + nrelf2a(i))
      if (lbsmat(nrelf2a(i)+nsft)) then
        radi = radf(i)
      else
        radi = 0.0_dp
      endif
      if (leem.and.lmultiqrange.and.neemrptr(i).ne.0) then
        nqri = nqrnow(neemrptr(i))
      else
        nqri = 1
      endif
      ind = ndistanceind(i)
!
!  Start of loop over second atom
!
      do j = 1,i
!
!  Set self term factor
!
        if (i.eq.j) then
          half = 0.5_dp
        else
          half = 1.0_dp
        endif
!
!  Set overall number of distances for i - j
!
        nor = ndistanceij(j,i)
!
!  Freezing flag
!
        lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
        if (.not.lopi.and..not.lopj) goto 1110
!
!  Arrange nat1 to be lower than nat2. If they are equal then arrange that ntyp1 is lower than ntyp2
!
        natj = nat(j)
        ntypj = nftype(j)
        nspecj = nspecptr(j)
        if (nati.eq.natj) then
          nat1 = nati
          nat2 = natj
          if (ntypi.lt.ntypj) then
            lorder12loc = .true.
            ntyp1 = ntypi
            ntyp2 = ntypj
          else
            lorder12loc = .false.
            ntyp1 = ntypj
            ntyp2 = ntypi
          endif
        elseif (nati.lt.natj) then
          lorder12loc = .true.
          nat1 = nati
          nat2 = nat(j)
          ntyp1 = ntypi
          ntyp2 = nftype(j)
        else
          lorder12loc = .false.
          nat1 = nat(j)
          nat2 = nati
          ntyp1 = nftype(j)
          ntyp2 = ntypi
        endif
        nregionj = nregionno(nsft+nrelf2a(j))
        nregiontypj = nregiontype(nregionj,ncf)
!
!  Set region 2 pair flag
!
        lreg2one  = .false.
        lreg2pair = .false.
        if (lseok.and.nregions(ncf).ge.2) then
          lreg2pair = (nregioni.eq.2.and.nregionj.eq.2)
          if (.not.lreg2pair) lreg2one = (nregioni.eq.2.or.nregionj.eq.2)
        endif
        lslicej = lsliceatom(nsft + nrelf2a(j))
        lattach = (lslicei.and..not.lslicej.or.lslicej.and..not.lslicei)
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
        if (QMMMmode(ncf).gt.0) then
          if (nregiontypi.eq.1.and.nregiontypj.eq.1) goto 1110
        endif
!           
!  QM/MM : Set electrostatic embedding flag : If either i or j are QM atoms => exclude electrostatics
!       
        lQMMMelectro = (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1))
!
        qlj = qf(j)
        ocj = occuf(j)
        if (lbsmat(nsft+nrelf2a(j))) then
          radj = radf(j)
        else
          radj = 0.0_dp
        endif
        radsum = radi + radj
        ofct = half*oci*ocj
        fct = ofct*angstoev
        factor = qli*qlj*fct
        if (lpolar) then
          polfct = abs(qli*dpolar(j)) + abs(qlj*dpolar(i))
        else
          polfct = 0.0_dp
        endif
        if (leem.and.lmultiqrange.and.neemrptr(j).ne.0) then
          nqrj = nqrnow(neemrptr(j))
        else
          nqrj = 1
        endif
!
!  Possible core-shell flag
!
        lcspair = (abs(nat1-nat2).eq.maxele.or.(oci+ocj).lt.1.0001d0)
        if (lcspair) then
          ospfct = half*oci
        else
          ospfct = ofct
        endif
!
!  Locate potential number
!  Check whether potential requires specific types
!  Calculate sum of all dispersion terms for pair
!
        if (i.ne.ilast.or.natj.ne.natjlast.or.ntypj.ne.ntypjlast) then
          rp = 0.0_dp
          npots = 0
          c6tot = 0.0_dp
          do n = 1,npote
            if (lmatch(nat1,ntyp1,nspec1(n),nptyp1(n),.true.)) then
              if (lmatch(nat2,ntyp2,nspec2(n),nptyp2(n),.true.)) then
                npots = npots + 1
                npotl(npots) = n
                if (nptype(n).eq.8.or.nptype(n).eq.33) then
                  if (cuts.gt.rp) rp = cuts
                elseif (lc6loc) then
                  if (nptype(n).eq.1.or.nptype(n).eq.7) then
                    c6tot = c6tot + twopot(3,n)
                    if (repcut(n).gt.rp) rp = repcut(n)
                  elseif (nptype(n).eq.2.or.nptype(n).eq.21) then
                    c6tot = c6tot + twopot(2,n)
                    if (repcut(n).gt.rp) rp = repcut(n)
                  elseif (nptype(n).eq.57) then
                    c6tot = c6tot + twopot(3,n)*twopot(4,n)*twopot(5,n)**6
                    if (repcut(n).gt.rp) rp = repcut(n)
                  else
                    if (rpot(n).gt.rp) rp = rpot(n)
                  endif
                else
                  if (rpot(n).gt.rp) rp = rpot(n)
                endif
              endif
            endif
          enddo
          cut2r = rp*rp
          if (cut2r.gt.cut2p) cut2r = cut2p
          cut2 = cut2r
!
!  If both charges are less than threshold then exclude electrostatics from cutoff
!
          if (abs(qli*oci)+abs(qlj*ocj).gt.thresh_q) then
            if (cut2e.gt.cut2.and.lewald) cut2 = cut2e
            if (cut2w.gt.cut2.and.lwolf) cut2 = cut2w
          endif
          if (lqeq.or.lSandM.or.lgfnff) cut2 = max(cut2,rq2)
!
          ilast = i
          natjlast = natj
          ntypjlast = ntypj
        endif
!
!  If no valid potentials, charge product is zero, and no polarisation then skip loop
!
        if (npots.eq.0.and.abs(factor).lt.1.0d-8.and.polfct.lt.1.0d-8) goto 1110
!
!  Find data for twobody interactions from within total store
!
        call getdistances(nor,ind,i,cut2,lself,nmolonly)
!
        if (lself) then
          d2self = 0.0_dp
          ereal2 = 0.0_dp
          erecip2 = 0.0_dp
          ec62 = 0.0_dp
          derive0self = 0.0_dp
          call selfterm(ereal2,erecip2,ec62,derive0self,factor,fct,ofct,ospfct,1.0_dp,npotl,npots, &
                        c6tot,d2self,lgrad1,.false.,i,j,ix,jx,1.0_dp,0.5_dp,lewaldtype,qli,qlj)
          esum = ereal2 + erecip2 + ec62
          if (lreg2one) then
            esregion12 = esregion12 + esum
          elseif (lreg2pair) then
            esregion2 = esregion2 + esum
          else
            ereal = ereal + ereal2
            erecip = erecip + erecip2
            ec6 = ec6 + ec62
          endif
          if (lattach) eattach = eattach + esum
!
          eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + esum
!
          siteenergy(i) = siteenergy(i) + 0.5_dp*esum
          siteenergy(j) = siteenergy(j) + 0.5_dp*esum
        endif
!*******************************************
!  Charge first derivatives for self-term  *
!*******************************************
        if (lself.and.lgrad1.and.lDoQDeriv1) then
          dqi(1) = derive0self*qli
          dqj(1) = derive0self*qlj
          call d1charge(i,j,lopi,lopj,1_i4,dqj,dqi)
        endif
!
        if (nor.eq.0) goto 1110
!*************************************************
!  Calculate twobody contribution in real space  *
!*************************************************
        eatom2 = 0.0_dp
        ereal2 = 0.0_dp
        ec62 = 0.0_dp
        call twobodymd(eatom2,ereal2,ec62,nor,1_i4,npots,npotl,cut2r,cut2q,cut2s,nmolonly,factor, &
                       ofct,ospfct,radsum,sctrm1,sctrm2,qli,qlj,nspeci,nspecj,lcspair,lewaldtype, &
                       .false.,lQMMMelectro)
        esum = eatom2 + ereal2 + ec62
        if (lreg2one) then
          esregion12 = esregion12 + esum
        elseif (lreg2pair) then
          esregion2 = esregion2 + esum
        else
          eatom = eatom + eatom2
          ereal = ereal + ereal2
          ec6 = ec6 + ec62
        endif
        if (lattach) eattach = eattach + esum
!
        eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + esum
!
        siteenergy(i) = siteenergy(i) + 0.5_dp*esum
        siteenergy(j) = siteenergy(j) + 0.5_dp*esum
!
        if (lPrintTwo) then
          write(ioout,'(4x,i12,i12,f22.10,1x,f22.10)') i,j,eatom2+ec62,ereal2
        endif
!
        if ((lDoQDeriv1.or.lDoQDeriv2).and.lgrad1) then
          do k = 1,nor
            d0i(k) = d0i(k) + derive0(k)*qlj
            d0j(k) = d0j(k) + derive0(k)*qli
          enddo
        endif
        if (leem.or.lgfnff) then
          if (lqeq.or.lSandM.or.lgfnff) then
            eqeq2 = 0.0_dp
            if (lqeq) then
              call qeqbody(eqeq2,lgrad1,.false.,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
            elseif (lSandM) then
              call smbody(eqeq2,lgrad1,.false.,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
            elseif (lgfnff) then
              zetaij = 1.0_dp/sqrt(gfnff_eeq_alp(i)+gfnff_eeq_alp(j))
              call gfnffbody(eqeq2,zetaij,lgrad1,.false.,nor,1_i4,fct,qli,qlj)
            endif
            if (lreg2one) then
              esregion12 = esregion12 + eqeq2
            elseif (lreg2pair) then
              esregion2 = esregion2 + eqeq2
            else
              eqeq = eqeq + eqeq2
            endif
            if (lattach) eattach = eattach + eqeq2
!
            eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eqeq2
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*eqeq2
            siteenergy(j) = siteenergy(j) + 0.5_dp*eqeq2
          endif
        endif
        if (lsuttonc) then
          if (.not.lMEAMden) then
            if (lorder12loc) then
              scrho(1,i) = scrho(1,i) + sctrm1*ocj
              scrho(1,j) = scrho(1,j) + sctrm2*oci
              if (lattach) then
                scrho12(1,i) = scrho12(1,i) + sctrm1*ocj
                scrho12(1,j) = scrho12(1,j) + sctrm2*oci
              endif
            else
              scrho(1,i) = scrho(1,i) + sctrm2*ocj
              scrho(1,j) = scrho(1,j) + sctrm1*oci
              if (lattach) then
                scrho12(1,i) = scrho12(1,i) + sctrm2*ocj
                scrho12(1,j) = scrho12(1,j) + sctrm1*oci
              endif
            endif
          endif
        endif
!
!  Generate products for derivatives
!
        if (lmdl.or.lsg1) then
          call twostrterms(ndim,maxdis,nor,xtmp,ytmp,ztmp,0.0_dp,0.0_dp,0.0_dp,dr2ds,rpd,d2r2dsdx,d2r2ds2,.false.)
        endif
!*****************************
!  Charge first derivatives  *
!*****************************
        if (lgrad1.and.lDoQDeriv1) then
          call d1charge(i,j,lopi,lopj,nor,d0i,d0j)
        endif
!*************************************
!  Electrostatic potential on-sites  *
!*************************************
        if (lpolar) then
          do k = 1,nor
            vx(i) = vx(i) - qlj*derivqd(k)*xtmp(k)
            vy(i) = vy(i) - qlj*derivqd(k)*ytmp(k)
            vz(i) = vz(i) - qlj*derivqd(k)*ztmp(k)
            vx(j) = vx(j) + qli*derivqd(k)*xtmp(k)
            vy(j) = vy(j) + qli*derivqd(k)*ytmp(k)
            vz(j) = vz(j) + qli*derivqd(k)*ztmp(k)
          enddo
          if (lattach) then
            do k = 1,nor
              vx12(i) = vx12(i) - qlj*derivqd(k)*xtmp(k)
              vy12(i) = vy12(i) - qlj*derivqd(k)*ytmp(k)
              vz12(i) = vz12(i) - qlj*derivqd(k)*ztmp(k)
              vx12(j) = vx12(j) + qli*derivqd(k)*xtmp(k)
              vy12(j) = vy12(j) + qli*derivqd(k)*ytmp(k)
              vz12(j) = vz12(j) + qli*derivqd(k)*ztmp(k)
            enddo
          endif
        endif
!***********************
!  Radial derivatives  *
!***********************
        if (lgrad1) then
          if (radi.gt.0.0_dp.and.lopi) then
            do k = 1,nor
              raderv(i) = raderv(i) + rtrm1(k)
            enddo
          endif
          if (radj.gt.0.0_dp.and.lopj) then
            do k = 1,nor
              raderv(j) = raderv(j) + rtrm1(k)
            enddo
          endif
        endif
!************************
!  Internal Derivatives *
!************************
!
!  First derivatives
!
        if (lgrad1) then
          if (lopi) then
            do k = 1,nor
              xdrv(i) = xdrv(i) - deriv(k)*xtmp(k)
              ydrv(i) = ydrv(i) - deriv(k)*ytmp(k)
              zdrv(i) = zdrv(i) - deriv(k)*ztmp(k)
            enddo
          endif
          if (lopj) then
            do k = 1,nor
              xdrv(j) = xdrv(j) + deriv(k)*xtmp(k)
              ydrv(j) = ydrv(j) + deriv(k)*ytmp(k)
              zdrv(j) = zdrv(j) + deriv(k)*ztmp(k)
            enddo
          endif
          if (nregioni.ne.nregionj) then
            do k = 1,nor
              xregdrv(nregioni) = xregdrv(nregioni) - deriv(k)*xtmp(k)
              yregdrv(nregioni) = yregdrv(nregioni) - deriv(k)*ytmp(k)
              zregdrv(nregioni) = zregdrv(nregioni) - deriv(k)*ztmp(k)
              xregdrv(nregionj) = xregdrv(nregionj) + deriv(k)*xtmp(k)
              yregdrv(nregionj) = yregdrv(nregionj) + deriv(k)*ytmp(k)
              zregdrv(nregionj) = zregdrv(nregionj) + deriv(k)*ztmp(k)
            enddo
          endif
        endif
!***********************
!  Strain derivatives  *
!***********************
!
!  Strain only terms 
!
!  First derivatives 
!
        if (lsg1) then
          rstrdloc(1:nstrains) = 0.0_dp
          do kl = 1,nstrains
            ks = nstrptr(kl)
            do k = 1,nor
              rstrdloc(kl) = rstrdloc(kl) + deriv(k)*dr2ds(k,ks)
            enddo
            rstrd(kl) = rstrd(kl) + rstrdloc(kl)
          enddo
          if (latomicstress) then
            do kl = 1,nstrains
              atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*rstrdloc(kl)
              atomicstress(kl,j) = atomicstress(kl,j) + 0.5_dp*rstrdloc(kl)
            enddo
          endif
        endif
1110    continue
        ind = ind + ndistanceij(j,i)
!
!  End of loop over j
!
      enddo
!
!  End of loop over i
!
    enddo
!****************
!  Global sums  *
!****************
    tsum0 = g_cpu_time()
    if (lsuttonc.and.nprocs.gt.1) then
      if (.not.lMEAMden) then
        do i = 1,numat
          sum2(i) = scrho(1,i)
        enddo
        call sumall(sum2,sum,numat,"realmd3","scrho")
        do i = 1,numat
          scrho(1,i) = sum(i)
        enddo
        do i = 1,numat
          sum2(i) = scrho12(1,i)
        enddo
        call sumall(sum2,sum,numat,"realmd3","scrho12")
        do i = 1,numat
          scrho12(1,i) = sum(i)
        enddo
      endif
    endif
    tsuml = g_cpu_time() - tsum0
    tsum = tsum + tsuml
!
!  Closing banner for energy decomposition
!
    if (lPrintTwo) then
      call mpbarrier
      if (ioproc) then
        write(ioout,'(''--------------------------------------------------------------------------------'',/)')
      endif
    endif
!
!  Exit point
!
  endif
!
!  Free local memory
!
  deallocate(sum2,stat=status)
  if (status/=0) call deallocate_error('realmd3t','sum2')
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('realmd3t','sum')
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('realmd3t','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  tatom = tatom + time2 - time1 - tsuml
#ifdef TRACE
  call trace_out('realmd3t')
#endif
!
  return
  end
