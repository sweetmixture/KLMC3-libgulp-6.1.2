  subroutine realsd2(eatom,ereal,erecip,ec6,eqeq,lgrad1,lgrad2)
!
!  Radial second derivs need sorting!!
!
!  Subroutine for calculating real space energy, 1st and 2nd derivatives using symmetry.
!
!  Originally created 4/95
!
!   8/95 Ewald sum for dispersion terms added
!  11/96 Compression of second derivatives added when lfreeze=.true.
!        Because i(opt)-j(frozen) d2 blocks are stored in i(asym)-
!        i(full) block, it is necessary to exclude self-terms in the
!        the second derivatives.
!   2/97 BSM exponential potential added
!   4/97 Sutton-Chen modifications added
!  12/97 Modification of energies and derivatives made purely local
!   1/98 QEq modifications added
!   4/98 ESFF form of Lennard-Jones now allowed for
!   1/99 1-4 interaction scaling added
!   7/00 Searching for valid vectors and self terms moved to subroutines
!   2/01 Structure of rpd changed to suit 2-D case
!   2/01 Calculation of electric field for polarisation added
!   9/01 lmolq calculations accelerated using lneedmol 
!   2/02 lneedmol algorithm corrected
!   5/02 BSM made core/shell specific
!  10/02 ReaxFF modifications added
!  11/02 Wildcard atoms added
!   1/03 Wolf modifications made
!   1/03 Call to selfterm altered
!   3/03 lgrad1p introduced to correct vx/vy/vz
!   6/04 nreli.ne.j condition removed from lmolok
!   9/04 Call to selfterm modified due to addition of charge first derivatives
!   9/04 Charge first derivatives added
!   1/05 rp no longer sqrt'd and passed to rsearch routines
!   4/05 Mods for cosh-spring added
!   7/05 Streitz and Mintmire modifications added
!   3/07 Printing of twobody energies added as an option
!   5/07 Argument list for twobody call modified
!  11/07 Modifications for reaxFF Coulomb term added
!   1/08 lreaxFFqreal removed
!   3/08 ReaxFF charge scheme removed as this is handled elsewhere
!  12/08 scrho, and defect versions, converted to 2-D array for MEAM
!   1/09 Arrays apot, bpot, cpot & dpot replaced by new 2-D array twopot
!   2/09 Expanding and contracting of maxdis arrays removed
!   3/09 lorder12loc added to twobody argument list
!   4/09 MEAM density modifications removed since these are now handled in 
!        a separate routine.
!   5/12 Atomic stresses added
!   5/12 Atomic stresses removed for routines involving symmetry
!   7/13 Call to twobody modified to include core-shell occupancy factor
!   3/14 lorder12loc removed twobody argument list
!  12/14 rtrm1 changed from scalar to array
!   2/15 MM3buck added
!   2/18 Trace added
!   5/18 Modified for q ranges
!   7/18 Setting of nqri corrected for symmetry
!   8/18 Call to twostrterms introduced for setting of rpd
!   9/18 Strain module introduced
!   9/18 Second derivatives completed here rather than in strfin
!   5/19 Skipping of ij loop now checks for polarisation
!   8/19 Short range damping of polarisation added
!   8/19 Trap for neemrptr being zero added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  11/19 Rigid molecule modifications added
!   2/20 Breathing shell contribution now made conditional on lopi
!   3/20 Electrostatic cutoff only included where the either charge exceeds the threshold
!   4/20 derv3c and d2r2dsdc added for benefit of rigid molecules
!   4/20 derv3c changes reversed as they are no longer required
!  11/21 Species added to twobody calls
!  12/21 pGFNFF modifications added
!
!  NOTE : freezing cannot be used with EEM/QEq and symmetry - algorithm
!         would involve performing loops over all pairs of atoms and so
!         it looses most of the advantage of the algorithm. Also symmetry
!         is very problematic with charge derivatives and hence this
!         algorithm is not used for EEM/QEq
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
  use configurations, only : lbsmat
  use g_constants
  use control
  use current
  use datatypes
  use derivatives
  use eam,            only : lMEAMden
  use eemdata
  use element
  use general,        only : cutw
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_rad
  use iochannels,     only : ioout
  use kspace
  use m_strain,       only : twostrterms
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
  logical,     intent(in)                      :: lgrad2
  real(dp),    intent(inout)                   :: eatom
  real(dp),    intent(inout)                   :: ec6
  real(dp),    intent(inout)                   :: eqeq
  real(dp),    intent(inout)                   :: ereal
  real(dp),    intent(inout)                   :: erecip
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indri
  integer(i4)                                  :: indrif
  integer(i4)                                  :: indrj
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixf
  integer(i4)                                  :: iyf
  integer(i4)                                  :: izf
  integer(i4)                                  :: ixfo
  integer(i4)                                  :: iyfo
  integer(i4)                                  :: izfo
  integer(i4)                                  :: ixi
  integer(i4)                                  :: iyi
  integer(i4)                                  :: izi
  integer(i4)                                  :: ixj
  integer(i4)                                  :: iyj
  integer(i4)                                  :: izj
  integer(i4)                                  :: j
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: jxc
  integer(i4)                                  :: jyc
  integer(i4)                                  :: jzc
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: ks
  integer(i4)                                  :: kt
  integer(i4)                                  :: m
  integer(i4)                                  :: n
  integer(i4)                                  :: n3a
  integer(i4)                                  :: n3f
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: nfa
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmj
  integer(i4)                                  :: nmolonly
  integer(i4)                                  :: nor
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: npots
  integer(i4)                                  :: nqri
  integer(i4)                                  :: nqrj
  integer(i4)                                  :: nff
  integer(i4)                                  :: nfi
  integer(i4)                                  :: nfif
  integer(i4)                                  :: nfj
  integer(i4)                                  :: nreli
  integer(i4)                                  :: nspeci
  integer(i4)                                  :: nspecj
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: lc6loc
  logical                                      :: lcspair
  logical                                      :: lewaldtype
  logical                                      :: lgrad1p
  logical                                      :: lmatch
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lorder12loc
  logical                                      :: lself
  logical                                      :: lsg1
  logical                                      :: lsg2
  real(dp)                                     :: apt
  real(dp)                                     :: bpt
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
  real(dp)                                     :: dfct
  real(dp)                                     :: eatomr
  real(dp)                                     :: eatom_before
  real(dp)                                     :: ec6r
  real(dp)                                     :: ec6_before
  real(dp)                                     :: eqeqr
  real(dp)                                     :: erealr
  real(dp)                                     :: ereal_before
  real(dp)                                     :: etrm
  real(dp)                                     :: etrm1
  real(dp)                                     :: etrm2
  real(dp)                                     :: factor
  real(dp)                                     :: fct
  real(dp)                                     :: fcti
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ofct2
  real(dp)                                     :: ospfct
  real(dp)                                     :: ospfct2
  real(dp)                                     :: polfct
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: radi
  real(dp)                                     :: radj
  real(dp)                                     :: radsum
  real(dp)                                     :: rdiff
  real(dp)                                     :: rneqi
  real(dp)                                     :: rp
  real(dp)                                     :: rq
  real(dp)                                     :: rq2
  real(dp)                                     :: rt2
  real(dp)                                     :: sctrm1
  real(dp)                                     :: sctrm2
  real(dp)                                     :: sderv2r(6,6)
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: wrk(6)
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcom
  real(dp)                                     :: ycom
  real(dp)                                     :: zcom
  real(dp)                                     :: xcomi
  real(dp)                                     :: ycomi
  real(dp)                                     :: zcomi
  real(dp)                                     :: xcrd
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
  real(dp)                                     :: zetaij
#ifdef TRACE
  call trace_in('realsd2')
#endif
!
  time1 = g_cpu_time()
!
!  Local variables
!
  lc6loc = (lc6.and.ndim.eq.3)
  lgrad1p = (lgrad1.or.lpolar)
  if (lfreeze) then
    nfa = 0
    nff = 0
    do i = 1,nasym
      if (lopf(i)) then
        nfa = nfa + 1
        nff = nff + neqv(i)
      endif
    enddo
    n3a = 3*nfa
    n3f = 3*nff
  else
    n3a = 3*nasym
    n3f = 3*numat
  endif
!
  eatomr = 0.0_dp
  erealr = 0.0_dp
  ec6r = 0.0_dp
  eqeqr = 0.0_dp
  lsg1 = (lgrad1.and.lstr)
  lsg2 = (lgrad2.and.lstr)
  if (lsg1) then
    do i = 1,nstrains
      wrk(i) = 0.0_dp
    enddo
    if (lsg2) then
      do i = 1,nstrains
        do j = 1,nstrains
          sderv2r(j,i) = 0.0_dp
        enddo
      enddo
    endif
  endif
!
!  Set the Coulomb term type based on dimensionality :
!
!  1-D         => 1/r
!  2-D and 3-D => erfc(seta*r)/r
!
  lewaldtype = (ndim.ne.1.or.lwolf)
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
  if (lnoreal) then
#ifdef TRACE
    call trace_out('realsd2')
#endif
    return
  endif
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
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('realsd2','npotl')
!***************************************************************
!  Atomistic and real space electrostatic component of energy  *
!***************************************************************
!
!  Outer loop over sites
!
  ix = - 2
  iy = - 1
  iz =   0
  ixf = 1
  iyf = 2
  izf = 3
  ixfo = 1
  iyfo = 2
  izfo = 3
  nfi = 0
  nfif = 1
  do i = 1,nasym
    lopi = (.not.lfreeze.or.lopf(i))
    if (.not.lopi) goto 1100
    nfi = nfi + 1
!
!  Inner loop over second site
!
    xal = xalat(i)
    yal = yalat(i)
    zal = zalat(i)
    nati = iatn(i)
    ntypi = natype(i)
    qli = qa(i)
    oci = occua(i)
    if (lbsmat(nsft+i)) then
      radi = rada(i)
      indri = n3a + nfi
      indrif = n3f + nfif
    else
      radi = 0.0_dp
    endif
!
    if (leem.and.lmultiqrange) then
      if (leembond.or..not.lsymopt) then
        if (neemrptr(nreli).ne.0) then
          nqri = nqrnow(neemrptr(nreli))
        else
          nqri = 1
        endif
      else
        if (neemrptr(i).ne.0) then
          nqri = nqrnow(neemrptr(i))
        else
          nqri = 1
        endif
      endif
    else
      nqri = 1
    endif
!
    nfif = nfif + neqv(i)
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    nreli = nrela2f(i)
    nspeci = nspecptr(nreli)
    fcti = oci*neqv(i)
    if (abs(fcti).gt.1.0d-15) then
      rneqi = 1.0_dp/fcti
    else
      rneqi = 1.0_dp/dble(neqv(i))
    endif
    ixf = ixfo
    iyf = iyfo
    izf = izfo
    ixfo = ixfo + 3*neqv(i)
    iyfo = iyfo + 3*neqv(i)
    izfo = izfo + 3*neqv(i)
!
!  Molecule handling
!
    if (lmol) then
      nmi = natmol(nreli)
      indm = nmolind(nreli)
      call mindtoijk(indm,ixi,iyi,izi)
!
!  Set COM coordinates
!
      if (lrigid.and.nmi.gt.0) then
        xcomi = molxyz(1,natinmol(nreli),nmi)
        ycomi = molxyz(2,natinmol(nreli),nmi)
        zcomi = molxyz(3,natinmol(nreli),nmi)
      else
        xcomi = 0.0_dp
        ycomi = 0.0_dp
        zcomi = 0.0_dp
      endif
    else
      xcomi = 0.0_dp
      ycomi = 0.0_dp
      zcomi = 0.0_dp
    endif
!
!  Start of second atom loop
!
    jxc = -2
    jyc = -1
    jzc =  0
    nfj = 0
    do j = 1,numat
      natj = nat(j)
      ntypj = nftype(j)
      nspecj = nspecptr(j)
!
!  Arrange nat1 to be lower than nat2. If they are equal then arrange that ntyp1 is lower than ntyp2
!
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
        ntyp2 = ntypj
      else
        lorder12loc = .false.
        nat1 = nat(j)
        nat2 = nati
        ntyp1 = ntypj
        ntyp2 = ntypi
      endif
!
!  Freeze flag
!
      lopj = (.not.lfreeze.or.lopf(nrelf2a(j)))
      ofct = fcti
      dfct = 1.0_dp
      if (lopj) nfj = nfj + 1
      xcrd = xclat(j) - xal
      ycrd = yclat(j) - yal
      zcrd = zclat(j) - zal
      qlj = qf(j)
      ocj = occuf(j)
!
      if (leem.and.lmultiqrange) then
        if (leembond.or..not.lsymopt) then
          if (neemrptr(j).ne.0) then
            nqrj = nqrnow(neemrptr(j))
          else
            nqrj = 1
          endif
        else
          if (neemrptr(nrelf2a(j)).ne.0) then
            nqrj = nqrnow(neemrptr(nrelf2a(j)))
          else
            nqrj = 1
          endif
        endif
      else
        nqrj = 1
      endif
!
      if (lbsmat(nsft+nrelf2a(j)).and.lopj) then
        radj = radf(j)
        indrj = n3f + nfj
      else
        radj = 0.0_dp
      endif
      radsum = radi + radj
      if (.not.lfreeze.or.lopj) then
        jxc = jxc + 3
        jyc = jyc + 3
        jzc = jzc + 3
        jx = jxc
        jy = jyc
        jz = jzc
      else
        jx = ixf
        jy = iyf
        jz = izf
      endif
      ofct = ofct*ocj
      fct = ofct*angstoev
      factor = qli*qlj*fct
      if (lpolar) then
        polfct = abs(qli*dpolar(nrelf2a(j))) + abs(qlj*dpolar(i))
      else
        polfct = 0.0_dp
      endif
!
!  Possible core-shell pair flag
!
      lcspair = (abs(nat1-nat2).eq.maxele.or.(oci+ocj).lt.1.0001_dp)
      if (lcspair) then
        ospfct = fcti
      else
        ospfct = ofct
      endif
!
!  Molecule handling
!
      if (lmol) then
        nmj = natmol(j)
        indmj = nmolind(j)
        call mindtoijk(indmj,ixj,iyj,izj)
        ixj = ixj - ixi
        iyj = iyj - iyi
        izj = izj - izi
        lmolok = (nmi.eq.nmj.and.nmi.gt.0)
!
!  Set COM coordinates
!
        if (lrigid.and.nmj.gt.0) then
          xcom = molxyz(1,natinmol(j),nmj) - xcomi
          ycom = molxyz(2,natinmol(j),nmj) - ycomi
          zcom = molxyz(3,natinmol(j),nmj) - zcomi
        else
          xcom = - xcomi
          ycom = - ycomi
          zcom = - zcomi
        endif
      else
        lmolok = .false.
        xcom = - xcomi
        ycom = - ycomi
        zcom = - zcomi
      endif
!
!  Locate potential number
!
      rp = 0.0_dp
      npots = 0
      c6tot = 0.0_dp
      lneedmol = (lmol.and..not.lmolq)
      do n = 1,npote
        if (lmatch(nat1,ntyp1,nspec1(n),nptyp1(n),.true.)) then
          if (lmatch(nat2,ntyp2,nspec2(n),nptyp2(n),.true.)) then
            npots = npots + 1
            npotl(npots) = n
            if (mmexc(n).gt.0.or.(lintra(n).and..not.linter(n).or..not.lintra(n).and.linter(n))) lneedmol = .true.
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
!
!  If no valid potentials, charge product is zero, and no polarisation then skip loop
!
      if (npots.eq.0.and.abs(factor).lt.1.0d-8.and.polfct.lt.1.0d-8) goto 1000
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
!  If the Coulomb term and the potentials don't need molecule info for
!  for this pair of atoms, then turn off molecule checking to save time
!
      if (.not.lneedmol) lmolok = .false.
!
      d2self = 0.0_dp
!***********************
!  Find valid vectors  *
!***********************
      if (ndim.eq.3) then
        call rsearch3D(xcrd,ycrd,zcrd,lmolok,lcspair,nreli,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
      elseif (ndim.eq.2) then
        call rsearch2D(xcrd,ycrd,zcrd,lmolok,lcspair,nreli,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
      elseif (ndim.eq.1) then
        call rsearch1D(xcrd,ycrd,zcrd,lmolok,lcspair,nreli,j,nmi,ixj,iyj,izj,nor,nmolonly,lself,cut2)
      endif
!
      if (lself) then
        derive0self = 0.0_dp
        if (lfreeze) then
          call selfterm(ereal,erecip,ec6r,derive0self,factor,fct,ofct,ospfct,dfct,npotl,npots,c6tot,d2self, &
                        lgrad1,lgrad2,nreli,j,ix,jx,1.0_dp,1.0_dp,.true.,qli,qlj)
        else
          call selfterm(ereal,erecip,ec6r,derive0self,factor,fct,ofct,ospfct,dfct,npotl,npots,c6tot,d2self, &
                        lgrad1,lgrad2,nreli,j,ix,jx,0.5_dp,1.0_dp,.true.,qli,qlj)
        endif
      endif
!
      if (nor.eq.0) goto 1000
!**********************************
!  Evaluate twobody contribution  *
!**********************************
      if (lfreeze.and.lopj) then
        ofct2 = ofct*0.5_dp
        ospfct2 = ospfct*0.5_dp
        dfct = 2.0_dp
      else
        ofct2 = ofct
        ospfct2 = ospfct
        dfct = 1.0_dp
      endif
!
!  Sqrt distances
!
      do k = 1,nor
        dist(k) = sqrt(dist(k))
      enddo
!
      eatom_before = eatomr
      ereal_before = erealr
      ec6_before   = ec6r
!
      call twobody(eatomr,erealr,ec6r,lgrad1p,lgrad2,.false.,nor,1_i4,npots,npotl,cut2r,cut2q,cut2s, &
                   nmolonly,factor,ofct2,ospfct2,radsum,sctrm1,sctrm2,qli,qlj,nspeci,nspecj,lcspair, &
                   lewaldtype,.false.,.false.,.false.)
!
      if (lPrintTwo) then
#ifdef KLMC_DEBUG_A
        ! wkjee
        write(ioout,'(A)') "wkjee realsd flag "
#endif
        write(ioout,'(4x,i12,i12,f22.10,1x,f22.10)') nreli,j,eatomr+ec6r-eatom_before-ec6_before,erealr-ereal_before
      endif
!
      if (leem.or.lgfnff) then
        if (lqeq) then
          call qeqbody(eqeqr,lgrad1,lgrad2,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
        elseif (lSandM) then
          call smbody(eqeqr,lgrad1,lgrad2,nor,1_i4,fct,qli,qlj,nati,natj,nqri,nqrj)
        elseif (lgfnff) then
          zetaij = 1.0_dp/sqrt(gfnff_eeq_alp(nreli)+gfnff_eeq_alp(j))
          call gfnffbody(eqeqr,zetaij,lgrad1,lgrad2,nor,1_i4,fct,qli,qlj)
        endif
      endif
!
      if (lsuttonc) then
        if (.not.lMEAMden) then
          if (lorder12loc) then
            scrho(1,i) = scrho(1,i) + sctrm1*ocj
          else
            scrho(1,i) = scrho(1,i) + sctrm2*ocj
          endif
        endif
      endif
      if (lgrad1) then
        do k = 1,nor
          deriv(k) = dfct*deriv(k)
        enddo
        if (lgrad2) then
          do k = 1,nor
            deriv2(k) = dfct*deriv2(k)
            rderiv(k) = dfct*rderiv(k)
          enddo
        endif
      endif
      if (lsg1.or.lgrad2) then
!
!  Generate products for derivatives
!
        call twostrterms(ndim,maxdis,nor,xtmp,ytmp,ztmp,xcom,ycom,zcom,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
      endif
!*************************************
!  Electrostatic potential on-sites  *
!*************************************
      if (lpolar.and.j.ne.nreli) then
        do k = 1,nor
          vx(i) = vx(i) - qlj*derivqd(k)*xtmp(k)*rneqi
          vy(i) = vy(i) - qlj*derivqd(k)*ytmp(k)*rneqi
          vz(i) = vz(i) - qlj*derivqd(k)*ztmp(k)*rneqi
        enddo
      endif
!***********************
!  Radial derivatives  *
!***********************
      if (lgrad1) then
        if (lgrad2) then
          rt2 = 0.0_dp
          do k = 1,nor
            rt2 = rt2 + rtrm2(k)
          enddo
          rt2 = dfct*rt2
        endif
        if (lbsmat(nsft+i)) then
          do k = 1,nor
            raderv(i) = raderv(i) + dfct*rtrm1(k)
          enddo
          if (lgrad2) then
            derv2(indrif,indri) = derv2(indrif,indri) + rt2
          endif
        endif
        if (lgrad2.and.lopj) then
          if (lbsmat(nsft+i).and.radj.gt.0.0_dp) then
            derv2(indrj,indri) = derv2(indrj,indri) + rt2
          endif
        endif
      endif
!************************
!  Internal Derivatives *
!************************
!
!  First derivatives
!
      if (lgrad1) then
        do k = 1,nor
          xdrv(i) = xdrv(i) - deriv(k)*xtmp(k)
          ydrv(i) = ydrv(i) - deriv(k)*ytmp(k)
          zdrv(i) = zdrv(i) - deriv(k)*ztmp(k)
        enddo
        if (lgrad2.and.lbsmat(nsft+i)) then
          do k = 1,nor
            derv2(ixf,indri) = derv2(ixf,indri) - rderiv(k)*xtmp(k)
            derv2(iyf,indri) = derv2(iyf,indri) - rderiv(k)*ytmp(k)
            derv2(izf,indri) = derv2(izf,indri) - rderiv(k)*ztmp(k)
            derv2(indrif,ix) = derv2(indrif,ix) - rderiv(k)*xtmp(k)
            derv2(indrif,iy) = derv2(indrif,iy) - rderiv(k)*ytmp(k)
            derv2(indrif,iz) = derv2(indrif,iz) - rderiv(k)*ztmp(k)
          enddo
        endif
      endif
!
!  Second derivatives
!
      if (lgrad2.and.j.ne.nreli) then
!
!  Coordinate only
!
        do k = 1,nor
          derv2(jx,ix) = derv2(jx,ix) - deriv2(k)*d2r2dx2(k,1)
          derv2(jy,ix) = derv2(jy,ix) - deriv2(k)*d2r2dx2(k,6)
          derv2(jz,ix) = derv2(jz,ix) - deriv2(k)*d2r2dx2(k,5)
          derv2(jx,iy) = derv2(jx,iy) - deriv2(k)*d2r2dx2(k,6)
          derv2(jy,iy) = derv2(jy,iy) - deriv2(k)*d2r2dx2(k,2)
          derv2(jz,iy) = derv2(jz,iy) - deriv2(k)*d2r2dx2(k,4)
          derv2(jx,iz) = derv2(jx,iz) - deriv2(k)*d2r2dx2(k,5)
          derv2(jy,iz) = derv2(jy,iz) - deriv2(k)*d2r2dx2(k,4)
          derv2(jz,iz) = derv2(jz,iz) - deriv2(k)*d2r2dx2(k,3)
          derv2(jx,ix) = derv2(jx,ix) - deriv(k)
          derv2(jy,iy) = derv2(jy,iy) - deriv(k)
          derv2(jz,iz) = derv2(jz,iz) - deriv(k)
        enddo
!
!  Coordinate - radius mixed
!
        if (radj.gt.0.0_dp.and.lopj) then
          do k = 1,nor
            derv2(indrj,ix) = derv2(indrj,ix) - rderiv(k)*xtmp(k)
            derv2(indrj,iy) = derv2(indrj,iy) - rderiv(k)*ytmp(k)
            derv2(indrj,iz) = derv2(indrj,iz) - rderiv(k)*ztmp(k)
          enddo
        endif
        if (lbsmat(nsft+i).and.lopj) then
          do k = 1,nor
            derv2(jx,indri) = derv2(jx,indri) + rderiv(k)*xtmp(k)
            derv2(jy,indri) = derv2(jy,indri) + rderiv(k)*ytmp(k)
            derv2(jz,indri) = derv2(jz,indri) + rderiv(k)*ztmp(k)
          enddo
        endif
      endif
!***********************
!  Strain derivatives  *
!***********************
      if (lsg2) then
!
!  Mixed derivatives
!
        do kl = 1,nstrains
          ks = nstrptr(kl)
          do k = 1,nor
            derv3(ix,kl) = derv3(ix,kl) - deriv(k)*d2r2dsdx(k,ks,1) - xtmp(k)*deriv2(k)*dr2ds(k,ks)
            derv3(iy,kl) = derv3(iy,kl) - deriv(k)*d2r2dsdx(k,ks,2) - ytmp(k)*deriv2(k)*dr2ds(k,ks)
            derv3(iz,kl) = derv3(iz,kl) - deriv(k)*d2r2dsdx(k,ks,3) - ztmp(k)*deriv2(k)*dr2ds(k,ks)
          enddo
        enddo
        if (radi.gt.0.0_dp) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            do k = 1,nor
              derv3(indri,kl) = derv3(indri,kl) + rderiv(k)*dr2ds(k,ks)
            enddo
          enddo
        endif
      endif
!
!  First derivatives 
!
      if (lsg1) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          do k = 1,nor
            wrk(kl) = wrk(kl) + deriv(k)*dr2ds(k,ks)
          enddo
        enddo
!
!  Second derivatives
!
        if (lsg2) then
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              do k = 1,nor
                sderv2r(kl,kk) = sderv2r(kl,kk) + deriv2(k)*dr2ds(k,kt)*dr2ds(k,ks)
                sderv2r(kl,kk) = sderv2r(kl,kk) + deriv(k)*d2r2ds2(k,kt,ks)
              enddo
            enddo
          enddo
        endif
      endif
1000  continue
    enddo
!
!  Breathing shell self term
!
    if (lbsmat(nsft+i).and.lopi) then
      do m = 1,npote
        if (nptype(m).eq.14) then
          if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
            apt = twopot(1,m)*fcti
            rdiff = radi - twopot(2,m)
            if (lfreeze) then
              eatomr = eatomr + 0.5_dp*apt*rdiff*rdiff
            else
              eatomr = eatomr + apt*rdiff*rdiff
            endif
            if (lgrad1) then
              raderv(i) = raderv(i) + apt*rdiff
              if (lgrad2) then
                derv2(indrif,indri) = derv2(indrif,indri) + apt
              endif
            endif
          endif
        elseif (nptype(m).eq.17) then
          if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
            apt = twopot(1,m)*fcti
            bpt = twopot(2,m)
            rdiff = radi - twopot(3,m)
            etrm1 = exp(bpt*rdiff)
            etrm2 = 1.0_dp/etrm1
            etrm = apt*(etrm1 + etrm2)
            if (lfreeze) then
              eatomr = eatomr + etrm
            else
              eatomr = eatomr + 2.0_dp*etrm
            endif
            if (lgrad1) then
              raderv(i) = raderv(i) + apt*bpt*(etrm1-etrm2)
              if (lgrad2) then
                derv2(indrif,indri) = derv2(indrif,indri) + bpt*bpt*etrm
              endif
            endif
          endif
        elseif (nptype(m).eq.31) then
          if (lmatch(nati,ntypi,nspec1(m),nptyp1(m),.true.)) then
            apt = twopot(1,m)*fcti
            bpt = twopot(2,m)
            rdiff = radi - twopot(3,m)
            etrm1 = exp(bpt*rdiff)
            etrm = apt*etrm1
            if (lfreeze) then
              eatomr = eatomr + etrm
            else
              eatomr = eatomr + 2.0_dp*etrm
            endif
            if (lgrad1) then
              raderv(i) = raderv(i) + apt*bpt*etrm1
              if (lgrad2) then
                derv2(indrif,indri) = derv2(indrif,indri) + bpt*bpt*etrm
              endif
            endif
          endif
        endif
      enddo
    endif
!
!  Skip to here if i is frozen
!
1100 continue
  enddo
!
!  Double counting correction
!
  if (.not.lfreeze) then
    eatomr = 0.5_dp*eatomr
    erealr = 0.5_dp*erealr
    ec6r = 0.5_dp*ec6r
    eqeqr = 0.5_dp*eqeqr
    if (lsg1) then
      do i = 1,nstrains
        wrk(i) = 0.5_dp*wrk(i)
      enddo
    endif
    if (lsg2) then
      do i = 1,nstrains
        do j = 1,i
          sderv2r(i,j) = 0.5_dp*sderv2r(i,j)
        enddo
      enddo
    endif
  endif
!
!  Add local quantities to overall totals
!
  eatom = eatom + eatomr
  ereal = ereal + erealr
  ec6 = ec6 + ec6r
  eqeq = eqeq + eqeqr
  if (lsg1) then
    do i = 1,nstrains
      rstrd(i) = rstrd(i) + wrk(i)
    enddo
    if (lsg2) then
      do i = 1,nstrains
        do j = 1,nstrains
          sderv2(j,i) = sderv2(j,i) + sderv2r(j,i)
        enddo
      enddo
    endif
  endif
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
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('realsd2','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  trls = trls + time2 - time1
#ifdef TRACE
  call trace_out('realsd2')
#endif
!
  return
  end
