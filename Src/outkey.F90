  subroutine outkey
!
!  Output keywords
!
!   6/97 added keywords 'pred' and 'cost' (SMW)
!   5/98 modified keywords 'cost' etc. (SMW)
!   5/02 lfbfgs added
!  12/02 nomolecularinternalke added
!   1/03 spatial keyword added
!   2/03 Style updated
!   6/03 llbfgs added
!   4/05 pureq keyword added
!   7/05 Streitz and Mintmire modifications added
!   8/06 libff keyword added
!   9/06 libdump keyword added
!  11/06 NEB modifications added
!   1/07 Gasteiger added
!   1/07 lnoautobond added
!   3/07 lPrintEAM keyword added
!   3/07 lPrintTwo keyword added
!   3/07 dynamical keyword added
!   3/07 lpreserveQ added
!   5/07 lmeanke added
!   5/07 qbond added
!   7/07 lmeta added
!   4/08 Simplex keyword added
!  10/08 COSMO/COSMIC keywords added
!  11/08 lPrintFour keyword added
!  11/08 lPrintThree keyword added
!   1/09 nomcediff keyword added
!   5/09 numerical keyword added
!   6/09 Module name changed from three to m_three
!   6/09 PDF additions made
!   7/09 EVB added
!   7/09 Modified for lnoksym keyword
!  12/09 pregionforce keyword added
!   4/10 qtpie keyword added
!   8/10 lconvert and lphase keywords removed
!   8/10 lfix1atom added
!  12/10 Hiding of shells added
!   1/11 Force minimisation added
!   2/11 Output for norxQ keyword added
!   3/11 lstressout added
!   7/11 noaddshells added
!   8/11 nosderv keyword now tested for in keyword line
!   9/11 Metadynamics internal code replaced with Plumed
!   9/11 Madelung correction added
!  11/11 eregion added
!   3/12 Keyword added to suppress printing of fractions
!   5/12 Atomic stress keyword added
!   6/12 Thermal conductivity adde
!   7/12 Oldinten keyword added
!   9/12 Pacha added
!   9/12 Eckart transformation added
!  10/12 Option added to allow two atoms co-exist on a site with an
!        occupancy greater than one.
!  12/12 Trap for minimum image and spatial both being specified added
!  12/12 site_energy keyword added
!   7/13 Trap for minimum image and spatial corrected
!   7/13 Steepest descent added as an option on Hessian reset
!   8/13 linclude_imaginary added
!   8/13 lfindsym added
!   8/13 Raman keyword added
!  10/13 lsopt keyword added
!   2/14 abs keyword added
!   3/14 Harmonic relaxation keyword added
!   6/14 lnoexclude flag used
!   8/14 lgroupvelocity keyword added
!   1/15 nowrap added as a pseudonym for nomodcoord
!   3/15 kcal and kjmol keywords added
!   3/15 lnoqeem added
!   5/15 num3 added
!   7/15 shopt added
!   8/15 New lower options added
!   5/16 SPME added
!   1/17 lPrintSix added
!   3/17 loldvarorder added
!   5/17 lnewdefalg added
!   6/17 Module files renamed to gulp_files
!   8/17 Use of lnumerical added
!  10/17 Delta_dipole added
!  11/17 lmsd added
!  12/17 Fastfd keyword added
!   1/18 ghostcell keyword added
!   1/18 Grueneisen parameters added
!   2/18 lalamode added
!   3/18 Frame keyword added
!   4/18 dcharge keyword added
!   4/18 Trap for RFO and unit added
!   4/18 Setting of leigen based on lower or intensity removed
!   4/18 Split bond EEM added
!   6/18 lallbonds added
!   6/18 Strain cell added
!  11/18 fsprop added
!   2/19 Rigid molecules added
!   7/19 lceigen added
!   8/19 Output of c6 keyword added
!   9/19 lnotorXduplicates flag added
!   9/19 lnoshellzero added
!  10/19 Trap added for rigid with no molecule keyword
!  11/19 lxray added
!   3/20 lPrintVar added
!   3/20 lmolstdframe added
!   4/20 lmolcom_mass added
!   4/20 lnorotate added
!   5/20 lphasecom added
!   7/20 Use of lkcal and lkjmol added
!   8/20 lzdipole added
!   9/20 lpolarisation added
!   9/21 lPrintVB option added
!  11/21 Modifications for TI added
!  11/21 lmd_inc_com added
!  12/21 pGFNFF modifications added
!   2/22 lsync added
!   2/22 lequipartition added
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
!  Copyright Curtin University 2022
!
!  Julian Gale, CIC, Curtin University, February 2022
!  Scott Woodley, R.I.G.B., June 1997
!
  use bondvalence,    only : lPrintVB
  use configurations, only : ncfg
  use control
  use cosmic,         only : lcosmic
  use distances,      only : lStoreVectors
  use eam,            only : lPrintEAM
  use element
  use four,           only : lPrintFour
  use general
  use gulp_gfnff,     only : lgfnff_topowolf
  use iochannels
  use m_pdfneutron,   only : lmakeeigarray, lcoreinfo, lpdf
  use m_pdfneutron,   only : lnowidth, lpartial, lfreqcut, lkeepcut, lnoksym
  use m_precondition, only : luseprecon
  use m_three,        only : lPrintThree
  use m_ti,           only : ltirun
  use mdlogic
  use molecule
  use g_neb,          only : lnebclimbingimage, lnebdoublynudged
  use optimisation,   only : loptcellpar, lfix1atom, loldvarorder, lPrintVar
  use parallel
  use phonout,        only : linclude_imaginary
  use six,            only : lPrintSix
  use spme,           only : lspme
  use symmetry
  use two,            only : lPrintTwo
  use wolfcosmo,      only : lPureCoulomb0D
  implicit none
!
!  Local variables
!
  integer(i4) :: i
!
!  If this is not the I/O node then skip printing
!
  if (.not.ioproc) goto 10
!***************************
!  Output keyword details  *
!***************************
  if (lfit) then
    write(ioout,'(''*  fit          - perform fitting run                                          *'')')
    if (lopt) then
      write(ioout,'(''*  optimise     - perform optimisation run after fitting                       *'')')
    endif
    if (index(keyword,'simp').eq.1.or.index(keyword,' simp').ne.0) then
      write(ioout,'(''*  simplex      - use the simplex algorithm for fitting                        *'')')
    endif
  elseif (lpredict) then
    continue
  elseif (ltran) then
    write(ioout,'(''*  transition   - transition state search by rfo method                        *'')')
  elseif (lsync) then
    write(ioout,'(''*  synchronous  - transition state search using the synchronous transit method *'')')
  elseif (lnebclimbingimage) then
    write(ioout,'(''*  cineb        - TS search by climbing-image nudged elastic band method       *'')')
    if (.not.lnebdoublynudged) then
      write(ioout,'(''*  nodneb       - no double nudge of elastic band                              *'')')
    endif
  elseif (lsync) then
    write(ioout,'(''*  synchronous  - transition state search using the synchronous transit method *'')')
  elseif (lneb) then
    write(ioout,'(''*  neb          - transition state search by nudged elastic band method        *'')')
    if (.not.lnebdoublynudged) then
      write(ioout,'(''*  nodneb       - no double nudge of elastic band                              *'')')
    endif
  elseif (lopt.or.lrfo) then
    write(ioout,'(''*  optimise     - perform optimisation run                                     *'')')
  elseif (lgrad) then
    write(ioout,'(''*  gradient     - perform gradient run                                         *'')')
  elseif (lharmrelax) then
    write(ioout,'(''*  eharmonic    - compute harmonic relaxation energy                           *'')')
  elseif (lnoenergy) then
    write(ioout,'(''*  noenergy     - perform a single point run but do not calculate the energy   *'')')
  else
    write(ioout,'(''*  single       - perform a single point run                                   *'')')
  endif
  if (lforcemin) then
    write(ioout,'(''*  force_minim  - perform a force minimisation run                             *'')')
  endif
  if (lmc) then
    write(ioout,'(''*  montecarlo   - perform a Monte Carlo run                                    *'')')
    if (index(keyword,' nomc').ne.0.or.index(keyword,'nomc').eq.1) then
      write(ioout,'(''*  nomcediff    - use full energy calculation at each Monte Carlo step         *'')')
    endif
  endif
  if (lmd) then
    write(ioout,'(''*  md           - perform molecular dynamics run                               *'')')
    if (lmd_inc_com) then
      write(ioout,'(''*  inc_com      - include centre of mass motion in molecular dynamics run      *'')')
    endif
    if (lmd_inc_rot) then
      write(ioout,'(''*  inc_rot      - include rotation in molecular dynamics run                   *'')')
    endif
  endif
  if (ltirun) then
    write(ioout,'(''*  ti           - perform thermodynamic integration                            *'')')
  endif
  if (lsopt) then
    write(ioout,'(''*  sopt         - output shell iterations during molecular dynamics run        *'')')
  endif
  if (lminimage) then
    write(ioout,'(''*  minimum_image- perform calculation using the minimum image approach         *'')')
  endif
  if (lfree) then
    write(ioout,'(''*  free_energy  - use Gibbs free energy instead of internal energy             *'')')
  endif
  if (lequipartition) then
    write(ioout,'(''*  equipartition- use classical equipartition Gibbs free energy                *'')')
  endif
  if (lzsisa) then
    write(ioout,'(''*  zsisa        - only calculate free energy derivatives w.r.t. strain         *'')')
  endif
  if (lstaticfirst) then
    write(ioout,'(''*  static_first - optimise static energy before free energy                    *'')')
  endif
  if (lbulknoopt) then
    write(ioout,'(''*  bulk_noopt   - do not optimise during bulk calculation                      *'')')
  endif
  if (lconp) then
    write(ioout,'(''*  conp         - constant pressure calculation                                *'')')
  elseif (lconv) then
    write(ioout,'(''*  conv         - constant volume calculation                                  *'')')
  elseif (lcello) then
    write(ioout,'(''*  cell_only    - only optimise unit cell                                      *'')')
  elseif (lnoflags) then
    write(ioout,'(''*  noflags      - no flags to be read in, assumed to all be zero               *'')')
  endif
  if (lshello) then
    write(ioout,'(''*  shell        - only shells and radii to be optimised                        *'')')
  endif
  if (loptcellpar) then
    write(ioout,'(''*  ocell        - use cell parameters as the optimisation variables            *'')')
  endif
  if (lstraincell) then
    write(ioout,'(''*  strain       - use strains with respect to reference cell                   *'')')
  endif
  if (.not.lfix1atom) then
    write(ioout,'(''*  unfix        - do not fix one atom during optimisation                      *'')')
  endif
  if (index(keyword,' brea').ne.0.or.index(keyword,'brea').eq.1) then
    write(ioout,'(''*  breathe      - only radii to be optimised                                   *'')')
  endif
  if (index(keyword,'nobr').ne.0) then
    write(ioout,'(''*  nobreathe    - radii not to be optimised                                    *'')')
  endif
  if (index(keyword,'iso').ne.0) then
    write(ioout,'(''*  isotropic    - only isotropic cell expansion to be allowed                  *'')')
  elseif (index(keyword,'ort').ne.0) then
    write(ioout,'(''*  orthorhombic - only cell length expansion to be allowed                     *'')')
  endif
  if (index(keyword,'noel').ne.0) then
    write(ioout,'(''*  noelectro    - do not include electrostatic terms despite charges present   *'')')
  endif
  if (index(keyword,'qok').ne.0) then
    write(ioout,'(''*  qok          - running with non charge neutral cell is OK                   *'')')
    if (lmadelung) then
      write(ioout,'(''*  madelung     - add simple cubic Madelung correction to charged system       *'')')
    endif
  endif
  if (lc6) then
    write(ioout,'(''*  c6           - use real and reciprocal space sum for C6 terms               *'')')
  endif
  if (lgfnff) then
    write(ioout,'(''*  gfnff        - use the GFNFF of Spicher and Grimme                          *'')')
  endif
  if (lgfnff_topowolf) then
    write(ioout,'(''*  gwolf        - apply a Wolf sum to the topological charges in GFNFF         *'')')
  endif
  if (lspme) then
    write(ioout,'(''*  spme         - use smooth particle mesh Ewald where possible                *'')')
  endif
  if (lpreserveQ) then
    write(ioout,'(''*  preserve_Q   - preserve coordinate line charges over species input          *'')')
  endif
  if (lallowgt1) then
    write(ioout,'(''*  allow_gt_1   - allow sum of partial occupancies on a site to exceed one     *'')')
  endif
  if (lddipole) then
    write(ioout,'(''*  delta_dipole - use dipolar correction energy term based on change in dipole *'')')
  elseif (ldipole) then
    write(ioout,'(''*  dipole       - use dipolar correction energy term                           *'')')
  elseif (lzdipole) then
    write(ioout,'(''*  zdipole      - use zdipole correction energy term                           *'')')
  endif
  if (lframe) then
    write(ioout,'(''*  frame        - rotate molecules to the standard frame of reference          *'')')
  endif
  if (lcosmic) then
    write(ioout,'(''*  cosmic       - use COSMIC solvation energy (integral charges)               *'')')
  elseif (lcosmo) then
    write(ioout,'(''*  cosmo        - use COSMO solvation energy                                   *'')')
  endif
  if ((lcosmo.or.lcosmic).and.literativeQ) then
    write(ioout,'(''*  qiterative   - use iterative solution of SAS charges in COSMO/COSMIC        *'')')
  elseif (literativeQ) then
    write(ioout,'(''*  qiterative   - use iterative solution of charges in ReaxFF / EEM / QEq      *'')')
  endif
  if (luseprecon) then
    write(ioout,'(''*  precondition - use inverse matrix as EEM preconditioner                     *'')')
  endif
  if (.not.lreaxFFQ) then
    write(ioout,'(''*  norxQ        - do not compute charges in ReaxFF                             *'')')
  endif
  if (lPureCoulomb0D) then
    write(ioout,'(''*  pureQ        - use pure Coulomb potential for A matrix in COSMO/IC for 0D   *'')')
  endif
  if (lprop) then
    write(ioout,'(''*  property     - calculate properties for final geometry                      *'')')
  endif
  if (lpolarisation) then
    write(ioout,'(''*  polarisation - calculate bulk polarisation as part of properties            *'')')
  endif
  if (lstraincellprop) then
    write(ioout,'(''*  fsprop       - calculate properties using finite strain derivatives         *'')')
  endif
  if (lceigen) then
    write(ioout,'(''*  ceigen       - output eigen properties of elastic constant tensor           *'')')
  endif
  if (lnumerical) then
    write(ioout,'(''*  numerical    - calculate properties using numerical second derivatives      *'')')
  endif
  if (lfastfd) then
    write(ioout,'(''*  fastfd       - use fast algorithm for finite difference second derivatives  *'')')
  endif
  if (index(keyword,'num3').ne.0) then
    write(ioout,'(''*  num3         - calculate force constants using numerical third derivatives  *'')')
  endif
  if (index(keyword,'shop').ne.0) then
    write(ioout,'(''*  shopt        - optimise shells for numerical third derivatives              *'')')
  endif
  if (index(keyword,'oldu').ne.0) then
    write(ioout,'(''*  oldunits     - use 10^11 dyne/cm2 as default units for mechanical data      *'')')
  endif
  if (lphon) then
    write(ioout,'(''*  phonon       - calculate phonons for final geometry                         *'')')
  endif
  if (lghost) then
    write(ioout,'(''*  ghostcell    - use smaller ghost cell for phonon calculation                *'')')
  endif
  if (lalamode) then
    write(ioout,'(''*  alamode      - generate force constants for use with Alamode                *'')')
  endif
  if (leckart) then
    write(ioout,'(''*  eckart       - remove rotation/translation from vibrations of molecules     *'')')
  endif
  if (index(keyword,'nono').ne.0) then
    write(ioout,'(''*  nononanal    - exclude non-analytic correction to gamma point phonons       *'')')
  endif
  if (index(keyword,'dyna').ne.0) then
    write(ioout,'(''*  dynamical    - output dynamical matrix                                      *'')')
  endif
  if (lmeanke) then
    write(ioout,'(''*  meanke       - calculate mean kinetic energy per atom during phonon calc    *'')')
  endif
  if (lthermal) then
    write(ioout,'(''*  thermal      - calculate thermal conductivity for final geometry            *'')')
  endif
  if (lgroupvelocity) then
    write(ioout,'(''*  groupvelocity- group velocities to be computed at k points                  *'')')
  endif
  if (lgrueneisen) then
    write(ioout,'(''*  grueneisen   - Grueneisen parameters to be computed at k points             *'')')
  endif
  if (lkfull) then
    write(ioout,'(''*  kfull        - generate K points for the full centred cell                  *'')')
  endif
  if (ldefect) then
    write(ioout,'(''*  defect       - perform defect calculation after bulk run                    *'')')
  endif
  if (ldefect.and.lfreq) then
    write(ioout,'(''*  frequencies  - calculate defect frequencies                                 *'')')
  endif
  if (ldefect.and.(index(keyword,'regi').ne.0)) then
    write(ioout,'(''*  regi_before  - output region 1 at start of defect calculation               *'')')
  endif
  if (ldefect.and.(index(keyword,'r234').ne.0)) then
    write(ioout,'(''*  r234d        - use region 2a displacements in 3 and 4 body energy           *'')')
  endif
  if (ldefect.and.(index(keyword,'noan').ne.0)) then
    write(ioout,'(''*  noanisotropic- do not treat region 2b anisotropically                       *'')')
  endif
  if (ldefect.and.lnewdefalg) then
    write(ioout,'(''*  newda        - new modified algorithm to be used for defect calculations    *'')')
  endif
  if (leigen) then
    write(ioout,'(''*  eigenvectors - output phonon eigenvectors                                   *'')')
  endif
  if (lga) then
    write(ioout,'(''*  genetic      - use genetic algorithm to search conformation space           *'')')
  endif
  if (lanneal) then
    write(ioout,'(''*  annealling   - use simulated annealing to search conformation space         *'')')
  endif
  if (lpredict) then
    write(ioout,'(''*  predict      - use global minimiser                                         *'')')
    if (lopt) then
      write(ioout,'(''*  optimise     - will optimise best possible structures found w.r.t. energy   *'')')
    else
      write(ioout,'(''*               - will not optimise best possible structures found wrt energy  *'')')
    endif
    if (index(keyword,'cost').ne.0) then
      write(ioout,'(''*  cost         - will use Pannetier type cost function in global minimisers   *'')')
    endif
  else
    if (index(keyword,'cost').ne.0) then
      if (lopt) then      
        write(ioout,'(''*  cost         - will use Pannetier type cost function rather than energy     *'')')
      else
        write(ioout,'(''*  cost         - will calculate the cost function and energy for the lattice  *'')')
      endif
    endif      
  endif
  if (loldvarorder) then
    write(ioout,'(''*  oldvarorder  - use old variable order for optimisation                      *'')')
  endif
  if (lnotorXduplicates) then
    write(ioout,'(''*  no4duplicate - do not duplicate torsions with wildcard potentials           *'')')
  endif
  if (linclude_imaginary) then
    write(ioout,'(''*  include_imag - include imaginary frequencies in DOS and disperion output    *'')')
  endif
  if (linten) then
    write(ioout,'(''*  intensity    - calculate phonon eigenvectors and estimate IR intensities    *'')')
  endif
  if (loldinten) then
    write(ioout,'(''*  oldintensity - use old formula for IR intensities                           *'')')
  endif
  if (lmsd) then
    write(ioout,'(''*  msd          - calculate phonon mean squared displacements                  *'')')
  endif
  if (lraman) then
    write(ioout,'(''*  raman        - calculate Raman susceptibility                               *'')')
  endif
  if (llower) then
    write(ioout,'(''*  lower_symtry - reduce symmetry according to imaginary phonon modes          *'')')
  endif
  if (loptlower) then
    write(ioout,'(''*  optlower     - optimise after lowering symmetry due to imaginary modes      *'')')
  endif
  if (leregion) then
    write(ioout,'(''*  eregion      - output region - region interaction energies                  *'')')
  endif
  if (lsiteenergy) then
    write(ioout,'(''*  site_energy  - output energy decomposition per atomic site                  *'')')
  endif
  if (index(keyword,'pot').ne.0) then
    write(ioout,'(''*  potential    - calculate electrostatic site potentials                      *'')')
  endif
  if (index(keyword,'zer').eq.1.or.index(keyword,' zer').ne.0) then
    write(ioout,'(''*  zero_pot     - set the average site potential to zero                       *'')')
  endif
  if (index(keyword,'nodp').ne.0) then
    write(ioout,'(''*  nodpsym      - calculate electrostatic site potentials for all of region 1  *'')')
  endif
  if (index(keyword,'nod3').ne.0) then
    write(ioout,'(''*  nod3         - do not calculate third derivatives for ShengBTE              *'')')
  endif
  if (index(keyword,'efg').ne.0) then
    write(ioout,'(''*  efg          - calculate electrostatic electric field gradient              *'')')
  endif
  if (.not.lsym) then
    write(ioout,'(''*  nosymmetry   - turn off symmetry after initial structure generation         *'')')
    if (index(keyword,' full').ne.0.or.index(keyword,'full').eq.1) then
      write(ioout,'(''*  full         - generate full unit cell when symmetry is removed             *'')')
    endif
  elseif (index(keyword,' full').ne.0.or.index(keyword,'full').eq.1) then
    write(ioout,'(''*  full         - this directive has no effect unless nosymmetry is specified  *'')')
  endif
  if (lsymoff) then
    write(ioout,'(''*  symoff       - turn symmetry off after building structure                   *'')')
  endif
  if (lfindsym) then
    write(ioout,'(''*  find_symmetry- try to find the space group from the initial structure       *'')')
  endif
  if (lStoreVectors) then
    write(ioout,'(''*  storevector  - store a list of interatomic vectors                          *'')')
  endif
  if (.not.lmodco) then
    if (index(keyword,' nowr').ne.0.or.index(keyword,'nowr').eq.1) then
      write(ioout,'(''*  nowrap       - do not wrap coordinates                                      *'')')
    else
      write(ioout,'(''*  nomodcoord   - do not mod input coordinates                                 *'')')
    endif
  endif
  if (index(keyword,'deci').eq.1.or.index(keyword,' deci').ne.0) then
    write(ioout,'(''*  decimal_only - only output decimal coordinates and not fractions            *'')')
  endif
  if (index(keyword,'simu').ne.0) then
    write(ioout,'(''*  simultaneous - relax shell positions and radii during fitting               *'')')
  endif
  if (lnoshellzero) then
    write(ioout,'(''*  noshellzero  - displace shell coordinates away from zero during fitting     *'')')
  endif
  if (lrelax) then
    write(ioout,'(''*  relax        - relax structure during fitting                               *'')')
  endif
  if (index(keyword,'abs ').ne.0.or.index(keyword,'abs').eq.1) then
    write(ioout,'(''*  abs          - use absolute value for selected variables in fitting         *'')')
  endif
  if (index(keyword,' nosd').ne.0.or.index(keyword,'nosd').eq.1) then
    write(ioout,'(''*  nosderv      - symmetry not to be used in first derivative calculation      *'')')
  endif
  if (index(keyword,'outc').ne.0) then
    write(ioout,'(''*  outcon       - dump constraints to restart file                             *'')')
  endif
  if (index(keyword,'nod2').ne.0) then
    write(ioout,'(''*  nod2sym      - symmetry not to be used for second derivatives               *'')')
  endif
  if (index(keyword,'nods').ne.0) then
    write(ioout,'(''*  nodsym       - symmetry not to be used in defect calculation                *'')')
  endif
  if (lspatial) then
    write(ioout,'(''*  spatial      - use spatial decomposition algorithm                          *'')')
  endif
  if (lgasteiger) then
    write(ioout,'(''*  gasteiger    - compute Gasteiger charges                                    *'')')
  endif
  if (lqbond) then
    write(ioout,'(''*  qbond        - compute bond increment charges                               *'')')
  endif
  if (leem) then
    if (lqtpie) then
      write(ioout,'(''*  QTPie        - use QTPie formulation of charge transfer to damp long-range  *'')')
    endif
    if (lqeq) then
      write(ioout,'(''*  QEq          - electronegativity equalisation (Rappe and Goddard III method)*'')')
    elseif (index(keyword,'oldeem').ne.0) then
      write(ioout,'(''*  oldeem       - electronegativity equalisation method (original parameters)  *'')')
    elseif (lSandM.and..not.lSandMnoZZ) then
      write(ioout,'(''*  smzz         - electronegativity equalisation (Streitz and Mintmire) + Z-Z  *'')')
    elseif (lSandM) then
      write(ioout,'(''*  sm           - electronegativity equalisation method (Streitz and Mintmire) *'')')
    elseif (lpacha) then
      write(ioout,'(''*  pacha        - electronegativity equalisation method (Pacha of Marc Henry)  *'')')
    else
      write(ioout,'(''*  eem          - electronegativity equalisation method (new parameters)       *'')')
    endif
  endif
  if (leembond) then
    write(ioout,'(''*  eembond      - split bond electronegativity equalisation method to be used  *'')')
  endif
  if (lallbonds) then
    write(ioout,'(''*  allbonds     - use all bonds in split bond electronegativity method         *'')')
  endif
  if (index(keyword,'dcha').ne.0) then
    write(ioout,'(''*  dcharge      - output first derivatives of charge distribution              *'')')
  endif
  if (loldd2q) then
    write(ioout,'(''*  oldd2q       - use old algorithm for variable charge second derivatives     *'')')
  endif
  if (lnoqeem) then
    write(ioout,'(''*  noQeem       - exclude Coulomb contribution from EEM calculation            *'')')
  endif
  if (lbond) then
    write(ioout,'(''*  bond         - calculate bond lengths                                       *'')')
  endif
  if (laver) then
    write(ioout,'(''*  average      - calculate average bond lengths after bond analysis           *'')')
  endif
  if (ldist) then
    write(ioout,'(''*  distance     - calculate distances                                          *'')')
  endif
  if (langle) then
    write(ioout,'(''*  angle        - calculate angles for valid three body terms                  *'')')
  endif
  if (ltors) then
    write(ioout,'(''*  torsion      - calculate torsion angles for valid four body terms           *'')')
  endif
  if (lmol) then
    if (lmolq) then
      write(ioout,'(''*  molq         - molecule option activated, but with Coulomb terms retained   *'')')
    elseif (lmolmec) then
      write(ioout,'(''*  molmec       - molecule option activated, Coulomb subtract 1-2/1-3 terms    *'')')
    else
      write(ioout,'(''*  molecule     - molecule option activated, Coulomb subtract within molecule  *'')')
    endif
    if (lmolfix) then
      write(ioout,'(''*  fix_molecule - fix atoms of molecule as at the start of the run             *'')')
    endif
    if (lrigid) then
      write(ioout,'(''*  rigid        - molecules will be treated as rigid bodies                    *'')')
      if (.not.lphasecom) then
        write(ioout,'(''*  nocomphase   - do not phase rigid molecules using centre of mass            *'')')
      endif
      if (lnorotate) then
        write(ioout,'(''*  norotate     - rotation of rigid molecules to be excluded from optimisation *'')')
      endif
      if (.not.lmolcom_mass) then
        write(ioout,'(''*  com_nomass   - centre of molecule is not to be mass-weighted                *'')')
      endif
    endif
    if (lnoautobond) then
      write(ioout,'(''*  noautobond   - do not compute bonds based on covalent radii                 *'')')
    endif
  endif
  if (lcomp.and.lopt) then
    write(ioout,'(''*  compare      - compare initial and final structures                         *'')')
  endif
  if (lunit) then
    write(ioout,'(''*  unit         - approximate initial Hessian by unit diagonal matrix          *'')')
  endif
  if (lfbfgs) then
    write(ioout,'(''*  fbfgs        - use full numerical Hessian in BFGS for fitting               *'')')
  endif
  if (lconj) then
    write(ioout,'(''*  conjugate    - use conjugate gradients minimiser                            *'')')
  endif
  if (index(keyword,'posi').ne.0) then
    write(ioout,'(''*  positive     - ensure Hessian acts as positive-definite in Newton-Raphson   *'')')
  endif
  if (lrfo) then
    write(ioout,'(''*  rfo          - optimisation step to be determined by RFO method             *'')')
  endif
  if (index(keyword,'hess').ne.0) then
    if (lrfo) then
      write(ioout,'(''*  hessian      - output eigenvalues and vectors of diagonalised Hessian       *'')')
    else
      write(ioout,'(''*  hessian      - output inverse Hessian matrix when calculated exactly        *'')')
    endif
  endif
  if (llbfgs) then
    write(ioout,'(''*  lbfgs        - use limited memory BFGS update                               *'')')
  endif
  if (ldfp) then
    write(ioout,'(''*  dfp          - use Davidon-Fletcher-Powell Hessian update                   *'')')
  endif
  if (index(keyword,'stee').ne.0) then
    write(ioout,'(''*  steepest     - use steepest descents after Hessian reset                    *'')')
  endif
  if (index(keyword,'nofr').ne.0) then
    write(ioout,'(''*  nofrequency  - suppress frequency output after phonon calculation           *'')')
  endif
  if (lbroad) then
    write(ioout,'(''*  broaden_dos  - broaden density of states curves                             *'')')
  endif
  if (index(keyword,'noks').ne.0.or.lnoksym) then
    write(ioout,'(''*  noksymmetry  - do not use Brillouin zone symmetry when generating k points  *'')')
  endif
  if (index(keyword,'opera').ne.0) then
    write(ioout,'(''*  operators    - output bulk symmetry operators                               *'')')
  endif
  if (lnoexclude) then
    write(ioout,'(''*  noexclude    - do not freeze out atoms with no derivatives                  *'')')
  endif
  if (index(keyword,'noli').ne.0) then
    write(ioout,'(''*  nolist_md    - do not use list based methods for 3 & 4-body terms in MD     *'')')
  endif
  if (index(keyword,'nokp').ne.0) then
    write(ioout,'(''*  nokpoints    - do not print out list of k points                            *'')')
  endif
  if (lsave) then
    write(ioout,'(''*  save         - save defect matrices for restart                             *'')')
  endif
  if (lrest) then
    write(ioout,'(''*  restore      - restore defect matrices from disk for restart                *'')')
  endif
  if (index(keyword,'noze').ne.0) then
    write(ioout,'(''*  nozeropt     - exclude the zero point energy from phonon/free energy calcns *'')')
  endif
  if (lmarvreg2) then
    write(ioout,'(''*  marvinSE     - calculate the surface energy as per MARVIN                   *'')')
  endif
  if (index(keyword,'libff').ne.0) then
    write(ioout,'(''*  libff        - read fitting flags from library files during fitting runs    *'')')
  endif
  if (index(keyword,'libd').ne.0) then
    write(ioout,'(''*  libdump      - dump literal symbols to restart file                         *'')')
  endif
  if (index(keyword,'norep').ne.0) then
    write(ioout,'(''*  norepulsive  - no automatic cutoff to be used for repulsive exponentials    *'')')
  endif
  if (index(keyword,'nomol').ne.0) then
    write(ioout,'(''*  nomolecular..- no molecular internal kinetic energy in initial MD velocity  *'')')
  endif
  if (lPrintTwo) then
    write(ioout,'(''*  prt_two      - print real space twobody energy contributions                *'')')
  endif
  if (lPrintThree) then
    write(ioout,'(''*  prt_three    - print threebody energy contributions                         *'')')
  endif
  if (lPrintFour) then
    write(ioout,'(''*  prt_four     - print fourbody energy contributions                          *'')')
  endif
  if (lPrintSix) then
    write(ioout,'(''*  prt_six      - print sixbody energy contributions                           *'')')
  endif
  if (lPrintEAM) then
    write(ioout,'(''*  prt_eam      - print atomic densities and energies in EAM                   *'')')
  endif
  if (lPrintVar) then
    write(ioout,'(''*  prt_var      - print list of optimisation variables                         *'')')
  endif
  if (lPrintVB) then
    write(ioout,'(''*  prt_vb       - print valence bond energy contributions                      *'')')
  endif
  if (lkjmol) then
    write(ioout,'(''*  kjmol        - print energy contributions in kJ/mol                         *'')')
  endif
  if (lkcal) then
    write(ioout,'(''*  kcal         - print energy contributions in kcal                           *'')')
  endif
  if (lhideshells) then
    write(ioout,'(''*  hideshells   - hide shells in the output                                    *'')')
  endif
  if (index(keyword,' noad').ne.0.or.index(keyword,'noad').eq.1) then
    write(ioout,'(''*  noaddshells  - do not automatically add shells where they are missing       *'')')
  endif
  if (lregionforce) then
    write(ioout,'(''*  pregionforce - print out cumulative force on each region                    *'')')
  endif
  if (lstressout) then
    write(ioout,'(''*  stress_out   - print out final stress tensor components                     *'')')
  endif
  if (latomicstress) then
    write(ioout,'(''*  atomic_stress- print out final atomic stress components                     *'')')
  endif
!
!  PDF keywords (ers29)
!
  if (lmakeeigarray) then
    write(ioout,'(''*  makeEigenArrays - store all eigenvectors and frequencies after calculation  *'')')
  endif
  if (lcoreinfo) then
    write(ioout,'(''*  coreinfo     - output atomic information for "cores" used                   *'')')
    write(ioout,'(''*                 in phonon calculations                                       *'')')
  endif
  if (lpdf) then
    write(ioout,'(''*  PDF          - calculate Pair Distribution Functions                        *'')')
    if (lxray) then 
      write(ioout,'(''*  xray         - use X-ray scattering factors rather than neutron             *'')')
    endif
    if (lnowidth) then 
      write(ioout,'(''*  nowidth      - suppress output of peak widths for PDFs                      *'')')
    endif
    if (.not.(lpartial)) then 
      write(ioout,'(''*  nopartial    - suppress output of partial Pair Distribution Functions       *'')')
    endif
    if (lfreqcut) then 
      write(ioout,'(''*  PDFcut       - use frequency cut-offs in PDF as set in neutron block        *'')')
      if (lkeepcut) then
        write(ioout,'(''*  PDFkeep      - set all w below wmin / above wmax to be wmin / wmax (as appropriate) '')')
      endif
    endif
  endif
  write(ioout,'(''********************************************************************************'')')
!********************
!  Write out title  *
!********************
  if (ntitle.gt.0) then
    do i = 1,ntitle
      write(ioout,'(''* '',a76,'' *'')') titleword(i)(1:76)
    enddo
    write(ioout,'(''********************************************************************************'',/)')
  endif
!************************
!  Final keyword tasks  *
!************************
!
!  Check for contradictory keyword combinations
!
10 if (lconp.and.lcello.and.lconv) then
    nwarn = nwarn + 1
    call outwarning('Contradictory optimisation types specified - conp selected',0_i4)
    lcello = .false.
  elseif (lcello.and.lconv.or.lconv.and.lconp.or.lcello.and.lconp) then
    nwarn = nwarn + 1
    call outwarning('Contradictory optimisation types specified - conp selected',0_i4)
    lcello = .false.
  endif
  if (lrfo.and.lfree) then
    nwarn = nwarn + 1
    call outwarning('RFO is only approximate with free energy',0_i4)
  endif
  if (lrfo.and.llbfgs) then
    call outerror('RFO keyword cannot be used with conjugate gradients',0_i4)
    call stopnow('outkey')
  endif
  if (lrfo.and.lconj) then
    call outerror('RFO keyword cannot be used with conjugate gradients',0_i4)
    call stopnow('outkey')
  endif
  if (lrfo.and.lunit) then
    call outerror('RFO keyword cannot be used with a unit Hessian',0_i4)
    call stopnow('outkey')
  endif
  if (ldefect.and.lqeq) then
    call outerror('QEq keyword cannot be used in a defect calculation',0_i4)
    call stopnow('outkey')
  endif
  if (ldefect.and.lSandM) then
    call outerror('sm keyword cannot be used in a defect calculation',0_i4)
    call stopnow('outkey')
  endif
  if (ldefect.and.lpacha) then
    call outerror('Pacha keyword cannot be used in a defect calculation',0_i4)
    call stopnow('outkey')
  endif
  if (ldefect.and.leem) then
    call outerror('EEM keyword cannot be used in a defect calculation',0_i4)
    call stopnow('outkey')
  endif
  if (ldefect.and.lqtpie) then
    call outerror('QTPie keyword cannot be used in a defect calculation',0_i4)
    call stopnow('outkey')
  endif
  if ((lmolmec.or.(lmol.and..not.lmolq)).and.lqeq) then
    call outerror('QEq keyword cannot be used with Coulomb subtraction',0_i4)
    call stopnow('outkey')
  endif
  if ((lmolmec.or.(lmol.and..not.lmolq)).and.lSandM) then
    call outerror('sm keyword cannot be used with Coulomb subtraction',0_i4)
    call stopnow('outkey')
  endif
  if ((lmolmec.or.(lmol.and..not.lmolq)).and.lpacha) then
    call outerror('Pacha keyword cannot be used with Coulomb subtraction',0_i4)
    call stopnow('outkey')
  endif
  if ((lmolmec.or.(lmol.and..not.lmolq)).and.leem) then
    call outerror('EEM keyword cannot be used with Coulomb subtraction',0_i4)
    call stopnow('outkey')
  endif
!
!  Trap gfnff keyword if pGFNFF library is not linked
!
#ifdef NOGFNFF
  if (lgfnff) then
    call outerror('gfnff keyword cannot be used as pGFNFF is not linked',0_i4)
    call stopnow('outkey')
  endif
#endif
!
!  Trap TI without MD
!
  if (ltirun.and..not.lmd) then
    call outerror('Thermodynamic integration cannot be used without molecular dynamics',0_i4)
    call stopnow('outkey')
  endif
!
!  Check that molecule keyword has been specified with rigid
!
  if (lrigid.and..not.lmolmec.and..not.lmol.and..not.lmolq) then
    call outerror('rigid specified without a molecule keyword',0_i4)
    call stopnow('outkey')
  endif
!
!  If both minimum image and spatial are specified then turn off minimum image
!
  if (lspatial.and.lminimage) then
    nwarn = nwarn + 1
    call outwarning('Contradictory algorithms specified - using spatial',0_i4)
    lminimage = .false.
  endif
!
!  Check that Monte Carlo has not been requested with more than one configuration per file
!
  if (lmc.and.ncfg.gt.1) then
    call outerror('MC is only possible with a single configuraton',0_i4)
    call stopnow('outkey')
  endif
!
!  Check that restore and save are not both set
!
  if (lsave.and.lrest) then
    call outerror('save and restore cannot both be specified',0_i4)
    call stopnow('outkey')
  endif
!
!  Check that both lkcal and lkjmol haven't been specified as they are mutually exclusive
!
  if (lkcal.and.lkjmol) then
    call outerror('kcal and kjmol cannot both be specified',0_i4)
    call stopnow('outkey')
  endif
!
  return
  end
