  subroutine initmemory
!
!  Initialises all module arrays to their default minimum size
!
!   5/06 Call to changemaxeamden added
!   7/06 Sixbody memory initiated
!   9/06 Call to changemaxiltor added
!  11/06 NEB modifications added
!   5/07 Call to changemaxbondq added
!   7/07 Initialisation of metadynamics arrays added
!   7/07 Initialisation of reaxFF arrays added
!   7/07 Initialisation of plane potential arrays added
!   4/08 Call to changemaxreaxffval3 added
!   1/09 swap move added to Monte Carlo
!   2/09 changemaxdis now initialised once here
!   2/09 Argument removed from changemaxdis call
!   6/09 EVB arrays added
!   9/10 Neutron scattering modifications added
!   9/10 Call to changemaxnpwtloc/changemaxnppa added
!   9/10 Call to changemaxedipspec added
!   8/11 Call to changemaxreaxFFfixQspec added
!   9/11 Metadynamics internal code replaced with Plumed
!   1/12 Call to changemaxobsmode added
!  10/12 Initialisation of OpenKIM memory added
!  12/12 Call to changemaxtdfield added
!   3/13 New defect arrays added
!  12/13 Call to changemaxfstrain added
!   3/15 Call to changemaxbondtype added
!   3/16 Call to changemaxnboz added
!   5/16 Call to changemaxmcswaps added
!   4/17 changemaxnpts2 renamed to changemaxcosmoA
!   2/18 Trace added
!   4/18 Call to changemaxeembond added
!   5/18 Call to changemaxqrange added
!   8/18 Call to changemaxkimodel removed
!   2/19 Call to changemaxmolat added
!   3/19 Multiple temperature ramps added
!   2/20 Calls to changemaxn3a and changemaxn3ma added
!   4/20 maxfqat changed to maxfreq since it is no longer just atom based
!  11/20 Tersoff reorganised
!  12/20 Call to changemaxconnectpermol added
!   7/21 Bond valence added
!  10/21 Transform added for Monte Carlo
!  11/21 Modifications for TI added
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
!  Julian Gale, CIC, Curtin University, November 2021
!
  use control,        only : lcosmo
  use m_gfnff_nbr3,   only : changemaxn3atom, changemaxn3nbr
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
#ifdef TRACE
  call trace_in('initmemory')
#endif
!
#ifdef KLMC_DEBUG_INITMEMORY
  write(*,'(A)') "in initmemory.F90 : call changemaxallnearseg"
#endif
  call changemaxallnearseg
#ifdef KLMC_DEBUG_INITMEMORY
  write(*,'(A)') "in initmemory.F90 : call changemaxat"
#endif
  call changemaxat
#ifdef KLMC_DEBUG_INITMEMORY
  write(*,'(A)') "in initmemory.F90 : call changemaxathbH"
#endif
  call changemaxathbH
#ifdef KLMC_DEBUG_INITMEMORY
  write(*,'(A)') "in initmemory.F90 : call changemaxatxbAB"
#endif
  call changemaxatxbAB
#ifdef KLMC_DEBUG_INITMEMORY
  write(*,'(A)') "in initmemory.F90 : call changemaxatot"
#endif
  call changemaxatot
#ifdef KLMC_DEBUG_INITMEMORY
  write(*,'(A)') "in initmemory.F90 : call changemaxbond"
#endif
  call changemaxbond
#ifdef KLMC_DEBUG_INITMEMORY
  write(*,'(A)') "in initmemory.F90 : call changemaxbondq"
#endif
  call changemaxbondq
#ifdef KLMC_DEBUG_INITMEMORY
  write(*,'(A)') "in initmemory.F90 : call changemaxbondvec"
#endif
  call changemaxbondvec
#ifdef KLMC_DEBUG_INITMEMORY
  write(*,'(A)') "in initmemory.F90 : call changemaxccspec"
#endif
  call changemaxccspec
  ! wkjee
#ifdef KLMC_DEBUG_INITMEMORY
  write(*,'(A)') "in initmemory.F90 : call changemaxcfg"
#endif
  ! wkjee end
  call changemaxcfg
  call changemaxcon
  call changemaxconnect
  call changemaxconnectpermol
  call changemaxcontot
  call changemaxdef
  call changemaxdis
  call changemaxedipspec
  call changemaxndistancetotal
  call changemaxnebreplicatot
  call changemaxeamden
  call changemaxeamspec
  call changemaxeamfnspec
  call changemaxeembond
  call changemaxextrapol
  call changemaxfgrad
  call changemaxfstrain
  call changemaxfit
  call changemaxfor
  call changemaxgcmcmol
  call changemaxiltor
  call changemaxfkpt
  call changemaxfreq(6_i4)
  call changemaxkpt
  call changemaxskpt
  call changemaxkvec
  call changemaxlambda
  call changemaxlib
  call changemaxlist3
  call changemaxlist4
  call changemaxmany
  call changemaxmcswaps
  call changemaxmcswapspec
  call changemaxmctrans
  call changemaxmol
  call changemaxmolat
  call changemaxmoleqv
  call changemaxn3atom
  call changemaxn3nbr
  call changemaxnbopot
  call changemaxnboa
  call changemaxnboo
  call changemaxnbor
  call changemaxnboq
  call changemaxnboq0
  call changemaxnboz
  call changemaxbondtype
  call changemaxnobo
  call changemaxnpts
  call changemaxcosmoA
  call changemaxnptsh
  call changemaxnptstot
  call changemaxnqstepfit
  call changemaxnwstepfit
  call changemaxnset
  call changemaxobs
  call changemaxobsmode
  call changemaxone
  call changemaxpot
  call changemaxpotsites
  call changemaxppt
  call changemaxproj
  call changemaxproji
  call changemaxpts
  call changemaxqrange
  call changemaxr1at
  call changemaxr2at
  call changemaxregion
  call changemaxsasparticles
  call changemaxsasparticlespart
  call changemaxsix
  call changemaxspec
  call changemaxtdfield
  call changemaxtempramp
  call changemaxthb
  call changemaxtitle
  call changemaxtotr1at
  call changemaxtotint
  call changemaxtotvac
  call changemaxtrialatom
  call changemaxuffspec
  call changemaxint
  call changemaxvac
  call changemaxvar
  call changemaxword
  call changemaxqatoms
  call changemaxreaxFFspec
  call changemaxreaxFFfixQspec
  call changemaxreaxFFval3
  call changemaxplanepot
  call changemaxreaxffspec
  call changemaxn3a
  call changemaxn3ma
  call changemaxvalbond
  call changemaxvalener
!
!  Neutron scattering
!
  call changemaxhold
  call changemaxqvector
!
!  COSMO
!
  if (lcosmo) then
    call changemaxnppa
    call changemaxnpwtloc
  endif
#ifdef TRACE
  call trace_out('initmemory')
#endif
!
  return
  end
