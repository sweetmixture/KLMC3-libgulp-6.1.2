  subroutine gfnffbody(eqeq,zeta,lgrad1,lgrad2,nor,nor0,factor,qli,qlj)
!
!  Calculates the correction to the energy and derivative terms for the
!  GFNFF scheme due to the use of an error function at short range.
!  Must follow call to twobody as this sets the derivative terms to zero.
!
!   9/20 Created from smbody
!   2/21 Correction to second derivatives
!  10/21 lgfnff moved to control
!
!  eqeq      = correction to energy due to GFNFF
!  zeta      = exponent for erf
!  lgrad1    = .true. if first derivatives are needed
!  lgrad2    = .true. if second derivatives are needed
!  nor       = upper bound to distances in array
!  nor0      = lower bound to distances in array
!  dist      = array of distances
!  deriv     = first derivative on return
!  deriv2    = second derivative on return
!  factor    = product of occupancies / sym factor
!  qli       = charge on i
!  qlj       = charge on j
!  d1i       = derivative of first derivative with respect to qi
!  d1j       = derivative of first derivative with respect to qj
!  d2i2      = 2nd derivative of energy with respect to qi/qi
!  d2ij      = 2nd derivative of energy with respect to qi/qj
!  d2j2      = 2nd derivative of energy with respect to qj/qj
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
  use control,       only : lgfnff
  use element
  use gulp_gfnff,    only : gfnff_eeq_rad
  use realvectors
  use shells
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!  
!  Passed variables
!     
  integer(i4), intent(in)    :: nor
  integer(i4), intent(in)    :: nor0
  logical,     intent(in)    :: lgrad1
  logical,     intent(in)    :: lgrad2
  real(dp),    intent(inout) :: eqeq
  real(dp),    intent(in)    :: factor
  real(dp),    intent(in)    :: qli
  real(dp),    intent(in)    :: qlj
  real(dp),    intent(in)    :: zeta
!
!  Local variables
!
  integer(i4)                :: i
  real(dp)                   :: dgam
  real(dp)                   :: d2gamr2
  real(dp)                   :: dtrm1
  real(dp)                   :: dtrm2
  real(dp)                   :: gam
  real(dp)                   :: r
  real(dp)                   :: rrl
!
!  If this isn't a GFNFF run then subroutine shouldn't have been called
!
  if (.not.lgfnff) return
#ifdef TRACE
  call trace_in('gfnffbody')
#endif
!************************
!  Loop over distances  *
!************************
  do i = nor0,nor
    r = dist(i)
!
!  Exclude distances outside maximum cutoff and core-shell contacts
!
    if (r.lt.gfnff_eeq_rad.and.r.gt.cuts) then
!*****************************
!  GFNFF scheme corrections  *
!*****************************
      call gfnff_gamma(zeta,r,gam,dgam,d2gamr2,lgrad1,lgrad2)
      rrl = 1.0_dp/r
!
!  Energy
!
      eqeq = eqeq + qli*qlj*factor*(gam - rrl)
!
      if (lgrad1) then
!
!  First derivatives
!
        dtrm1 = qli*qlj*factor*(dgam + rrl*rrl)*rrl
        if (lgrad2) then
!
!  Second derivatives
!
          dtrm2 = qli*qlj*factor*(d2gamr2 - 2.0_dp*rrl**3)
          dtrm2 = (dtrm2 - dtrm1)*rrl**2
          deriv2(i) = deriv2(i) + dtrm2
!
!  Charge derivative terms
!
!  Standard correction for gamma-1/r
!
          d1i(i) = d1i(i) + qlj*factor*(dgam + rrl*rrl)*rrl
          d1j(i) = d1j(i) + qli*factor*(dgam + rrl*rrl)*rrl
          d2ij(i) = d2ij(i) + factor*(gam - rrl)
        endif
        deriv(i) = deriv(i) + dtrm1
      endif
    endif
!*******************************
!  End of loop over distances  *
!*******************************
  enddo
#ifdef TRACE
  call trace_out('gfnffbody')
#endif
!
  return
  end
!*********************************************************************
!  Wrapper for gfnffbody
!*********************************************************************
  subroutine gfnffbody1(eqeq,zeta,lgrad1,lgrad2,nor,nor0,dist1,deriv1,deriv21, &
    factor,qli,qlj,d1i1,d1j1,d2i21,d2ij1,d2j21)
!
!  Wrapper for call to gfnffbody when number of distances is one.
!
!   9/20 Created from smbody1
!
!  eqeq      = correction to energy due to GFNFF
!  zeta      = exponent for erf
!  lgrad1    = .true. if first derivatives are needed
!  lgrad2    = .true. if second derivatives are needed
!  nor       = upper bound to distances in array
!  nor0      = lower bound to distances in array
!  dist      = array of distances
!  deriv     = first derivative on return
!  deriv2    = second derivative on return
!  factor    = product of occupancies / sym factor
!  qli       = charge on i
!  qlj       = charge on j
!  d1i       = derivative of first derivative with respect to qi
!  d1j       = derivative of first derivative with respect to qj
!  d2i2      = 2nd derivative of energy with respect to qi/qi
!  d2ij      = 2nd derivative of energy with respect to qi/qj
!  d2j2      = 2nd derivative of energy with respect to qj/qj
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, September 2020
!
  use realvectors
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!     
!  Passed variables
!     
  integer(i4), intent(in)    :: nor
  integer(i4), intent(in)    :: nor0
  logical,     intent(in)    :: lgrad1
  logical,     intent(in)    :: lgrad2
  real(dp),    intent(inout) :: d1i1
  real(dp),    intent(inout) :: d1j1
  real(dp),    intent(inout) :: d2i21
  real(dp),    intent(inout) :: d2ij1
  real(dp),    intent(inout) :: d2j21
  real(dp),    intent(inout) :: deriv1
  real(dp),    intent(inout) :: deriv21
  real(dp),    intent(inout) :: dist1
  real(dp),    intent(inout) :: eqeq
  real(dp),    intent(in)    :: factor
  real(dp),    intent(in)    :: qli
  real(dp),    intent(in)    :: qlj
  real(dp),    intent(in)    :: zeta
#ifdef TRACE
  call trace_in('gfnffbody1')
#endif
!
!  Set variables that should be in realvectors module
!
  dist(1) = dist1
  deriv(1) = deriv1
  deriv2(1) = deriv21
  d1i(1) = d1i1
  d1j(1) = d1j1
  d2i2(1) = d2i21
  d2ij(1) = d2ij1
  d2j2(1) = d2j21
!
!  Call qeqbody
!
  call gfnffbody(eqeq,zeta,lgrad1,lgrad2,nor,nor0,factor,qli,qlj)
!
!  Set return variables from realvectors module
!
  dist1 = dist(1)
  deriv1 = deriv(1)
  deriv21 = deriv2(1)
  d1i1 = d1i(1)
  d1j1 = d1j(1)
  d2i21 = d2i2(1)
  d2ij1 = d2ij(1)
  d2j21 = d2j2(1)
#ifdef TRACE
  call trace_out('gfnffbody1')
#endif
!
  return
  end
