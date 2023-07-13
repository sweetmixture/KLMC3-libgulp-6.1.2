  subroutine d2chargep2Dd(iloc,i,j,nor,xtmp,ytmp,ix,iy,iz,jx,jy,jz,dei,dej,d1ir,d1zi,d1jr,d1zj, &
                          d2i2r,d2ijr,d2j2r,d2self,d2trm1dij,d2trm1diz,d2trm1djz,dtrm1di, &
                          dtrm1dj)
!
!  Calculates the contribution to the phonon matrices
!  due to charge derivatives from EEM/QEq. 2-D version.
!  This version is for calls from distributed memory routines.
!
!   1/22 Created from d2chargep2D and d2charged
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
  use element
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed arguments
!
  integer(i4), intent(in)    :: i       ! Atom i
  integer(i4), intent(in)    :: iloc    ! Atom i number on local node
  integer(i4), intent(in)    :: ix
  integer(i4), intent(in)    :: iy
  integer(i4), intent(in)    :: iz
  integer(i4), intent(in)    :: j
  integer(i4), intent(in)    :: jx
  integer(i4), intent(in)    :: jy
  integer(i4), intent(in)    :: jz
  integer(i4), intent(in)    :: nor
  real(dp),    intent(in)    :: d1ir(*)
  real(dp),    intent(in)    :: d1jr(*)
  real(dp),    intent(in)    :: d1zi(*)
  real(dp),    intent(in)    :: d1zj(*)
  real(dp),    intent(in)    :: dtrm1di
  real(dp),    intent(in)    :: dtrm1dj
  real(dp),    intent(in)    :: d2i2r(*)
  real(dp),    intent(in)    :: d2ijr(*)
  real(dp),    intent(in)    :: d2j2r(*)
  real(dp),    intent(in)    :: d2self
  real(dp),    intent(in)    :: d2trm1dij
  real(dp),    intent(in)    :: d2trm1diz
  real(dp),    intent(in)    :: d2trm1djz
  real(dp),    intent(in)    :: dei(*)
  real(dp),    intent(in)    :: dej(*)
  real(dp),    intent(in)    :: xtmp(*)
  real(dp),    intent(in)    :: ytmp(*)
!
!  Local variables
!
  integer(i4)                :: iv
  integer(i4)                :: ixf
  integer(i4)                :: iyf
  integer(i4)                :: izf
  integer(i4)                :: k
  integer(i4)                :: kx
  integer(i4)                :: ky
  integer(i4)                :: kz
  integer(i4)                :: m
  integer(i4)                :: mnxx
  integer(i4)                :: mnxy
  integer(i4)                :: mnxz
  integer(i4)                :: mnyy
  integer(i4)                :: mnyz
  integer(i4)                :: mnzz
  integer(i4)                :: mx
  integer(i4)                :: my
  integer(i4)                :: mz
  real(dp)                   :: d1ix
  real(dp)                   :: d1iy
  real(dp)                   :: d1iz
  real(dp)                   :: d1jx
  real(dp)                   :: d1jy
  real(dp)                   :: d1jz
  real(dp)                   :: d2ijs
  real(dp)                   :: d2i2s
  real(dp)                   :: d2j2s
  real(dp)                   :: deisum
  real(dp)                   :: dejsum
  real(dp)                   :: g_cpu_time
  real(dp)                   :: time1
  real(dp)                   :: time2
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2chargep2Dd')
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
  d1ix = 0.0_dp
  d1iy = 0.0_dp
  d1iz = 0.0_dp
  d1jx = 0.0_dp
  d1jy = 0.0_dp
  d1jz = 0.0_dp
  d2ijs = d2self
  d2i2s = 0.0_dp
  d2j2s = 0.0_dp
  do iv = 1,nor
    d1ix = d1ix + d1ir(iv)*xtmp(iv)
    d1iy = d1iy + d1ir(iv)*ytmp(iv)
    d1iz = d1iz + d1zi(iv)
    d1jx = d1jx + d1jr(iv)*xtmp(iv)
    d1jy = d1jy + d1jr(iv)*ytmp(iv)
    d1jz = d1jz + d1zj(iv)
    d2ijs = d2ijs + d2ijr(iv)
    d2i2s = d2i2s + d2i2r(iv)
    d2j2s = d2j2s + d2j2r(iv)
  enddo
  d1iz = d1iz + d2trm1diz
  d1jz = d1jz + d2trm1djz
  d2ijs = d2ijs + d2trm1dij
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
    deisum = dtrm1di
    dejsum = dtrm1dj
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
    enddo
!
!  Terms involving Q derivatives for i/j and one other atom
!
    do m = 1,nqatoms(iloc)
      k = nqatomptr(m,iloc)
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
      derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,i)
      derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,i)
      derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,i)
      derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,i)
      derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,i)
      derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,i)
      derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,i)
      derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,i)
      derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,i)
    enddo
  endif
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2chargep2Dd')
#endif
!
  return
  end
