  subroutine d2chargepd(iloc,i,j,nor,ix,iy,iz,jx,jy,jz,xkv,ykv,zkv,dei,dej, &
                        d1ixr,d1iyr,d1izr,d1jxr,d1jyr,d1jzr,d2i2r,d2ijr, &
                        d2j2r,d2self,dei0,dej0,lreal,lDoSelf)
!
!  Calculates the contribution to the dynamical matrix due to charge derivatives from EEM/QEq. 
!  This version is for calls from distributed memory routines.
!
!  At present this is gamma point only.
!
!  lreal = .true.  => call from real space routine
!        = .false. => call from reciprocal space routine
!
!   1/22 Created from d2chargep and d2charged
!   1/22 lDoSelf added to d2chargepd arguments
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
!  Julian Gale, CIC, Curtin University, January 2022
!
  use configurations, only : nregionno
  use control
  use current
  use derivatives
  use eemdata
  use element
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_gam
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: i       ! Atom i
  integer(i4), intent(in)  :: iloc    ! Atom i number on local node
  integer(i4), intent(in)  :: ix
  integer(i4), intent(in)  :: iy
  integer(i4), intent(in)  :: iz
  integer(i4), intent(in)  :: j
  integer(i4), intent(in)  :: jx
  integer(i4), intent(in)  :: jy
  integer(i4), intent(in)  :: jz
  integer(i4), intent(in)  :: nor
  logical,     intent(in)  :: lreal
  logical,     intent(in)  :: lDoSelf
  real(dp),    intent(in)  :: dei(*)
  real(dp),    intent(in)  :: dej(*)
  real(dp),    intent(in)  :: dei0
  real(dp),    intent(in)  :: dej0
  real(dp),    intent(in)  :: d1ixr
  real(dp),    intent(in)  :: d1iyr
  real(dp),    intent(in)  :: d1izr
  real(dp),    intent(in)  :: d1jxr
  real(dp),    intent(in)  :: d1jyr
  real(dp),    intent(in)  :: d1jzr
  real(dp),    intent(in)  :: d2i2r(*)
  real(dp),    intent(in)  :: d2ijr(*)
  real(dp),    intent(in)  :: d2j2r(*)
  real(dp),    intent(in)  :: d2self
  real(dp),    intent(in)  :: xkv
  real(dp),    intent(in)  :: ykv
  real(dp),    intent(in)  :: zkv
!
!  Local variables     
!
  integer(i4)              :: iv    
  integer(i4)              :: ixf
  integer(i4)              :: iyf
  integer(i4)              :: izf
  integer(i4)              :: k
  integer(i4)              :: kx
  integer(i4)              :: ky 
  integer(i4)              :: kz  
  integer(i4)              :: m
  integer(i4)              :: mnxx
  integer(i4)              :: mnxy
  integer(i4)              :: mnxz
  integer(i4)              :: mnyy
  integer(i4)              :: mnyz
  integer(i4)              :: mnzz
  integer(i4)              :: mx
  integer(i4)              :: my
  integer(i4)              :: mz
  integer(i4)              :: nqr
  integer(i4)              :: nk
  real(dp)                 :: d1ix
  real(dp)                 :: d1iy
  real(dp)                 :: d1iz
  real(dp)                 :: d1jx
  real(dp)                 :: d1jy
  real(dp)                 :: d1jz
  real(dp)                 :: d2i2s
  real(dp)                 :: d2ijs
  real(dp)                 :: d2j2s
  real(dp)                 :: d2qk 
  real(dp)                 :: deisum
  real(dp)                 :: dejsum
  real(dp)                 :: dkix
  real(dp)                 :: dkiy
  real(dp)                 :: dkiz 
  real(dp)                 :: dkjx 
  real(dp)                 :: dkjy 
  real(dp)                 :: dkjz 
  real(dp)                 :: ock 
  real(dp)                 :: g_cpu_time
  real(dp)                 :: qlk 
  real(dp)                 :: time1
  real(dp)                 :: time2
  real(dp), parameter      :: tsqrt2pi = 0.797884560802866_dp
  real(dp)                 :: zetah0
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2chargepd')
#endif
!
  time1 = g_cpu_time()
!
  ixf = 3*(i-1) + 1
  iyf = ixf + 1
  izf = ixf + 2
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  d1ix = d1ixr
  d1iy = d1iyr
  d1iz = d1izr
  d1jx = d1jxr
  d1jy = d1jyr
  d1jz = d1jzr
  d2ijs = d2self
  d2i2s = 0.0_dp
  d2j2s = 0.0_dp
  do iv = 1,nor
    d2ijs = d2ijs + d2ijr(iv)
    d2i2s = d2i2s + d2i2r(iv)
    d2j2s = d2j2s + d2j2r(iv)
  enddo
  if ((leem.and..not.lelementOK(nat(i))).or.nregionno(nsft+nrelf2a(i)).ne.1) then
    d2ijs = 0.0_dp
    d2i2s = 0.0_dp
    d1ix = 0.0_dp
    d1iy = 0.0_dp
    d1iz = 0.0_dp
  endif
  if ((leem.and..not.lelementOK(nat(j))).or.nregionno(nsft+nrelf2a(j)).ne.1) then
    d2ijs = 0.0_dp
    d2j2s = 0.0_dp
    d1jx = 0.0_dp
    d1jy = 0.0_dp
    d1jz = 0.0_dp 
  endif
!
!  Add d1ix / d1iy / d1iz to sums for globalisation later
!
  dedqc(1,i) = dedqc(1,i) + d1ix
  dedqc(2,i) = dedqc(2,i) + d1iy
  dedqc(3,i) = dedqc(3,i) + d1iz
!
!  Add d1jx / d1iy / d1jz to sums for globalisation later
!
  d2edqc(ixf,j) = d2edqc(ixf,j) + d1jx
  d2edqc(iyf,j) = d2edqc(iyf,j) + d1jy
  d2edqc(izf,j) = d2edqc(izf,j) + d1jz
!******************
!  New algorithm  *
!******************
  d2edqdq(j,iloc) = d2edqdq(j,iloc) + d2ijs
  d2edq2(i) = d2edq2(i) + d2i2s
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lgfnff) then
    deisum = dei0
    dejsum = dej0
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
    enddo
!
!  Terms involving Q derivatives for i/j and one other atom
!
    do m = 1,nqatoms(iloc)
      k = nqatomptr(m,iloc)
      if (k.ne.i) then
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = ky + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
!
        mnxx = mx*(mx + 1)/2 
        mnyy = my*(my + 1)/2 
        mnzz = mz*(mz + 1)/2 
        mnxy = mnyy - 1
        mnxz = mnzz - 2
        mnyz = mnzz - 1
!
        derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,iloc)
        derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,iloc)
        derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,iloc)
        derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,iloc)
        derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,iloc)
        derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,iloc)
        derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,iloc)
        derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,iloc)
        derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,iloc)
      endif
    enddo
  endif
!
!  Loop over atoms to add ij-k contributions
!
  kx = - 2
  ky = - 1
  kz =   0
  do k = 1,numat
    kx = kx + 3
    ky = ky + 3
    kz = kz + 3
!
    if (lreal.and.lDoSelf.and.i.ne.j) then
!**************************************************
!  Contributions from derivatives of self energy  *
!**************************************************
      nk = nat(k)
      ock = occuf(k)
      if (leem) then
        if (lmultiqrange.and.neemrptr(k).ne.0) then
          nqr = nqrnow(neemrptr(k))
        else
          nqr = 1
        endif
        if (lqeq) then
          if (nk.ne.1) then
            d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
          else
            qlk = qf(k)
            zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
            d2qk = 2.0_dp*ock*murange(nqr,1,neemtype)*(1.0_dp+2.0_dp*qlk/zetah0)
          endif
        else
          d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
        endif
      elseif (lgfnff) then
        d2qk = ock*(tsqrt2pi/sqrt(gfnff_eeq_alp(k)) + gfnff_eeq_gam(k))*angstoev
      else
        d2qk = 0.0_dp
      endif
      dkix = dqdxyz(ixf,k)
      dkiy = dqdxyz(iyf,k)
      dkiz = dqdxyz(izf,k)
      dkjx = dqdxyz(jx,k)
      dkjy = dqdxyz(jy,k)
      dkjz = dqdxyz(jz,k)
      derv2(jx,ix) = derv2(jx,ix) + d2qk*dkjx*dkix
      derv2(jy,ix) = derv2(jy,ix) + d2qk*dkjx*dkiy
      derv2(jz,ix) = derv2(jz,ix) + d2qk*dkjx*dkiz
      derv2(jx,iy) = derv2(jx,iy) + d2qk*dkjy*dkix
      derv2(jy,iy) = derv2(jy,iy) + d2qk*dkjy*dkiy
      derv2(jz,iy) = derv2(jz,iy) + d2qk*dkjy*dkiz
      derv2(jx,iz) = derv2(jx,iz) + d2qk*dkjz*dkix
      derv2(jy,iz) = derv2(jy,iz) + d2qk*dkjz*dkiy
      derv2(jz,iz) = derv2(jz,iz) + d2qk*dkjz*dkiz
    endif
!
!  End of loop over k
!
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2chargepd')
#endif
!
  return
  end
