  subroutine d2charged(iloc,i,j,nor,ix,iy,iz,jx,jy,jz,dei,dej,d1xi,d1yi,d1zi, &
                       d1xj,d1yj,d1zj,ds1i,ds1j,d2i2,d2ij,d2j2,d2self,dei0,dej0, &
                       lreal,lDoSelf)
!
!  Calculates the contribution to the second derivative matrices
!  due to charge derivatives from variable charge models.
!  This version is for calls from distributed memory routines.
!
!  lreal = .true.  => call from real space routine
!        = .false. => call from reciprocal space routine
!  NB: loldd2q = .true. not supported
!
!   1/22 Created from d2charge
!   1/22 lNonIJQDeriv removed as not currently used
!   1/22 lopi/lopj removed from arguments as freezing is not allowed for this case
!   1/22 dedqc and d2edqc added for variable charge second derivatives
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
  use optimisation
  use parallel
  use m_strain,       only : strainddetds, straindet
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
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
  logical,     intent(in)    :: lreal
  logical,     intent(in)    :: lDoSelf
  real(dp),    intent(in)    :: d1xi          ! d2E/dqi.dx
  real(dp),    intent(in)    :: d1yi          ! d2E/dqi.dy
  real(dp),    intent(in)    :: d1zi          ! d2E/dqi.dz
  real(dp),    intent(in)    :: d1xj          ! d2E/dqj.dx
  real(dp),    intent(in)    :: d1yj          ! d2E/dqj.dy
  real(dp),    intent(in)    :: d1zj          ! d2E/dqj.dz
  real(dp),    intent(in)    :: d2i2(*)       ! d2E/dqi2
  real(dp),    intent(in)    :: d2ij(*)       ! d2E/dqi.dqj
  real(dp),    intent(in)    :: d2j2(*)       ! d2E/dqj2
  real(dp),    intent(in)    :: d2self        ! Self term
  real(dp),    intent(in)    :: dei(*)        ! dE/dqi
  real(dp),    intent(in)    :: dej(*)        ! dE/dqj
  real(dp),    intent(in)    :: dei0          ! dE(selfterm)/dqi
  real(dp),    intent(in)    :: dej0          ! dE(selfterm)/dqj
  real(dp),    intent(in)    :: ds1i(*)       ! d2E/dqi.d(epsilon)
  real(dp),    intent(in)    :: ds1j(*)       ! d2E/dqj.d(epsilon)
!
!  Local variables
!
  integer(i4)                :: ind
  integer(i4)                :: indi
  integer(i4)                :: indj
  integer(i4)                :: ixf
  integer(i4)                :: iyf
  integer(i4)                :: izf
  integer(i4)                :: iv
  integer(i4)                :: k
  integer(i4)                :: kk
  integer(i4)                :: kl
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
  integer(i4)                :: nqr
  integer(i4)                :: nk
  real(dp)                   :: d1ix
  real(dp)                   :: d1iy
  real(dp)                   :: d1iz
  real(dp)                   :: d1jx
  real(dp)                   :: d1jy
  real(dp)                   :: d1jz
  real(dp)                   :: d2i2s
  real(dp)                   :: d2ijs
  real(dp)                   :: d2j2s
  real(dp)                   :: d2qk
  real(dp)                   :: d1si
  real(dp)                   :: d2si
  real(dp)                   :: d3si
  real(dp)                   :: d4si
  real(dp)                   :: d5si
  real(dp)                   :: d6si
  real(dp)                   :: d1sj
  real(dp)                   :: d2sj
  real(dp)                   :: d3sj
  real(dp)                   :: d4sj
  real(dp)                   :: d5sj
  real(dp)                   :: d6sj
  real(dp)                   :: deisum
  real(dp)                   :: dejsum
  real(dp)                   :: dikx
  real(dp)                   :: diky
  real(dp)                   :: dikz
  real(dp)                   :: dis1
  real(dp)                   :: dis2
  real(dp)                   :: dis3
  real(dp)                   :: dis4
  real(dp)                   :: dis5
  real(dp)                   :: dis6
  real(dp)                   :: djkx
  real(dp)                   :: djky
  real(dp)                   :: djkz
  real(dp)                   :: djs1
  real(dp)                   :: djs2
  real(dp)                   :: djs3
  real(dp)                   :: djs4
  real(dp)                   :: djs5
  real(dp)                   :: djs6
  real(dp)                   :: dkix
  real(dp)                   :: dkiy
  real(dp)                   :: dkiz
  real(dp)                   :: dkjx
  real(dp)                   :: dkjy
  real(dp)                   :: dkjz
  real(dp)                   :: dk1
  real(dp)                   :: dk2
  real(dp)                   :: dk3
  real(dp)                   :: dk4
  real(dp)                   :: dk5
  real(dp)                   :: dk6
  real(dp)                   :: dsi(6)
  real(dp)                   :: dsi2(6)
  real(dp)                   :: dsj(6)
  real(dp)                   :: dsj2(6)
  real(dp)                   :: g_cpu_time
  real(dp)                   :: time1
  real(dp)                   :: time2
  real(dp)                   :: ock
  real(dp)                   :: qlk
  real(dp)                   :: zetah0
  real(dp), parameter        :: tsqrt2pi = 0.797884560802866_dp
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charged')
#endif
!
  time1 = g_cpu_time()
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  indi = 3*(i - 1)
  indj = 3*(j - 1)
  ixf = indi + 1
  iyf = indi + 2
  izf = indi + 3
!
  d2ijs = d2self
  d2i2s = 0.0_dp
  d2j2s = 0.0_dp
  do iv = 1,nor
    d2ijs = d2ijs + d2ij(iv)
    d2i2s = d2i2s + d2i2(iv)
    d2j2s = d2j2s + d2j2(iv)
  enddo
  d1ix = d1xi
  d1iy = d1yi
  d1iz = d1zi
  d1jx = d1xj
  d1jy = d1yj
  d1jz = d1zj
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
  if (lstr) then
    if (ndim.eq.3) then
      dis1 = dqds(1,i)
      dis2 = dqds(2,i)
      dis3 = dqds(3,i)
      dis4 = dqds(4,i)
      dis5 = dqds(5,i)
      dis6 = dqds(6,i)
      djs1 = dqds(1,j)
      djs2 = dqds(2,j)
      djs3 = dqds(3,j)
      djs4 = dqds(4,j)
      djs5 = dqds(5,j)
      djs6 = dqds(6,j)
    elseif (ndim.eq.2) then
      dis1 = dqds(1,i)
      dis2 = dqds(2,i)
      dis3 = dqds(3,i)
      djs1 = dqds(1,j)
      djs2 = dqds(2,j)
      djs3 = dqds(3,j)
    elseif (ndim.eq.1) then
      dis1 = dqds(1,i)
      djs1 = dqds(1,j)
    endif
    if (lreal) then
!
!  Real space
!
      do kl = 1,nstrains
        dsi(kl) = ds1i(kl)
        dsj(kl) = ds1j(kl)
        dsi2(kl) = dsi(kl)
        dsj2(kl) = dsj(kl)
      enddo
    else
!
!  Reciprocal space
!
      do kl = 1,nstrains
        dsi(kl) = ds1i(kl)
        dsj(kl) = ds1j(kl)
        dsi2(kl) = dsi(kl)
        dsj2(kl) = dsj(kl)
      enddo
      deisum = 0.0_dp
      dejsum = 0.0_dp
      do iv = 1,nor
        deisum = deisum + dei(iv)
        dejsum = dejsum + dej(iv)
      enddo
!
      if (ndim.eq.3) then
        if (lfinitestrain) then
          do kl = 1,3
            dsi(kl) = dsi(kl) - deisum*strainddetds(kl)*straindet
            dsj(kl) = dsj(kl) - dejsum*strainddetds(kl)*straindet
          enddo
        else
          do kl = 1,3
            dsi(kl) = dsi(kl) - deisum
            dsj(kl) = dsj(kl) - dejsum
          enddo
        endif
      elseif (ndim.eq.2) then
        if (lfinitestrain) then
          do kl = 1,3
            dsi(kl) = dsi(kl) - deisum*strainddetds(kl)*straindet
            dsj(kl) = dsj(kl) - dejsum*strainddetds(kl)*straindet
          enddo
        else
          do kl = 1,2
            dsi(kl) = dsi(kl) - deisum
            dsj(kl) = dsj(kl) - dejsum
          enddo
        endif
      elseif (ndim.eq.1) then
        dsi(1) = dsi(1) - deisum
        dsj(1) = dsj(1) - dejsum
      endif
      if (leem.or.lgfnff) then
        do kl = 1,nstrains
          dsi2(kl) = dsi(kl)
          dsj2(kl) = dsj(kl)
        enddo
      endif
    endif
!
!  Add contributions to arrays for globalisation and use later
!
    if (i.eq.j) then
      do kl = 1,nstrains
        ds2g(kl,i)   = ds2g(kl,i)   + 2.0_dp*dsi2(kl)
        ds2gs(kl,i)  = ds2gs(kl,i)  + 2.0_dp*d2ijs*dqds(kl,j)
        d2s2gs(kl,i) = d2s2gs(kl,i) + d2i2s*dqds(kl,i)
      enddo
    else
      do kl = 1,nstrains
        ds2g(kl,i)   = ds2g(kl,i)   + dsi2(kl)
        ds2gs(kl,i)  = ds2gs(kl,i)  + d2ijs*dqds(kl,j)
        d2s2gs(kl,i) = d2s2gs(kl,i) + d2i2s*dqds(kl,i)
      enddo
    endif
  endif
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
!  Strain - strain contribution
!   
    if (lstr.and.ndim.gt.0) then
      ind = 0
      do kk = 1,nstrains
        do kl = 1,kk-1
          ind = ind + 1
          sderv2(kl,kk) = sderv2(kl,kk) + deisum*d2qds2(ind,iloc)
          sderv2(kk,kl) = sderv2(kk,kl) + deisum*d2qds2(ind,iloc)
        enddo
        ind = ind + 1
        sderv2(kk,kk) = sderv2(kk,kk) + deisum*d2qds2(ind,iloc)
      enddo
    endif
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
      derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,iloc)
      derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,iloc)
      derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,iloc)
      derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,iloc)
      derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,iloc)
      derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,iloc)
      derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,iloc)
      derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,iloc)
      derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,iloc)
!
!  Mixed strain contribution
!
      if (lstr.and.ndim.gt.0) then
        do kk = 1,nstrains
! DEBUG - k terms not computed and so may not be correct
! DEBUG - bond order charges not currently enabled and so needs checking in future
          derv3(ix,kk) = derv3(ix,kk) - deisum*d2qdxyzs(kk,mx,iloc)
          derv3(iy,kk) = derv3(iy,kk) - deisum*d2qdxyzs(kk,my,iloc)
          derv3(iz,kk) = derv3(iz,kk) - deisum*d2qdxyzs(kk,mz,iloc)
        enddo
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
    dikx = dqdxyz(kx,i)
    diky = dqdxyz(ky,i)
    dikz = dqdxyz(kz,i)
    djkx = dqdxyz(kx,j)
    djky = dqdxyz(ky,j)
    djkz = dqdxyz(kz,j)
    if (lreal.and.lDoSelf) then
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
            d2qk = 2.0_dp*ock*murange(nqr,1,neemtype)*(1.0_dp + 2.0_dp*qlk/zetah0)
          endif
        else
          d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
        endif
      elseif (lgfnff) then
        d2qk = ock*(tsqrt2pi/sqrt(gfnff_eeq_alp(k)) + gfnff_eeq_gam(k))*angstoev
      else
        d2qk = 0.0_dp
      endif
      dkix = dqdxyz(indi+1,k)
      dkiy = dqdxyz(indi+2,k)
      dkiz = dqdxyz(indi+3,k)
      dkjx = dqdxyz(indj+1,k)
      dkjy = dqdxyz(indj+2,k)
      dkjz = dqdxyz(indj+3,k)
      derv2(jx,ix) = derv2(jx,ix) + d2qk*dkjx*dkix
      derv2(jy,ix) = derv2(jy,ix) + d2qk*dkjx*dkiy
      derv2(jz,ix) = derv2(jz,ix) + d2qk*dkjx*dkiz
      derv2(jx,iy) = derv2(jx,iy) + d2qk*dkjy*dkix
      derv2(jy,iy) = derv2(jy,iy) + d2qk*dkjy*dkiy
      derv2(jz,iy) = derv2(jz,iy) + d2qk*dkjy*dkiz
      derv2(jx,iz) = derv2(jx,iz) + d2qk*dkjz*dkix
      derv2(jy,iz) = derv2(jy,iz) + d2qk*dkjz*dkiy
      derv2(jz,iz) = derv2(jz,iz) + d2qk*dkjz*dkiz
      if (lstr.and.i.eq.j) then
!
!  Mixed strain terms from self energy - only need to do once
!  for each atom - hence check on i = j
!
        do kl = 1,nstrains
          derv3(ix,kl) = derv3(ix,kl) + d2qk*dqds(kl,k)*dkix
          derv3(iy,kl) = derv3(iy,kl) + d2qk*dqds(kl,k)*dkiy
          derv3(iz,kl) = derv3(iz,kl) + d2qk*dqds(kl,k)*dkiz
        enddo
!
!  Strain-strain derivatives - only need to do if i=j=k
!
        if (i.eq.k) then
          if (ndim.eq.3) then
            dk1 = dqds(1,k)
            dk2 = dqds(2,k)
            dk3 = dqds(3,k)
            dk4 = dqds(4,k)
            dk5 = dqds(5,k)
            dk6 = dqds(6,k)
            sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
            sderv2(2,1) = sderv2(2,1) + d2qk*dk2*dk1
            sderv2(3,1) = sderv2(3,1) + d2qk*dk3*dk1
            sderv2(4,1) = sderv2(4,1) + d2qk*dk4*dk1
            sderv2(5,1) = sderv2(5,1) + d2qk*dk5*dk1
            sderv2(6,1) = sderv2(6,1) + d2qk*dk6*dk1
            sderv2(2,2) = sderv2(2,2) + d2qk*dk2*dk2
            sderv2(3,2) = sderv2(3,2) + d2qk*dk3*dk2
            sderv2(4,2) = sderv2(4,2) + d2qk*dk4*dk2
            sderv2(5,2) = sderv2(5,2) + d2qk*dk5*dk2
            sderv2(6,2) = sderv2(6,2) + d2qk*dk6*dk2
            sderv2(3,3) = sderv2(3,3) + d2qk*dk3*dk3
            sderv2(4,3) = sderv2(4,3) + d2qk*dk4*dk3
            sderv2(5,3) = sderv2(5,3) + d2qk*dk5*dk3
            sderv2(6,3) = sderv2(6,3) + d2qk*dk6*dk3
            sderv2(4,4) = sderv2(4,4) + d2qk*dk4*dk4
            sderv2(5,4) = sderv2(5,4) + d2qk*dk5*dk4
            sderv2(6,4) = sderv2(6,4) + d2qk*dk6*dk4
            sderv2(5,5) = sderv2(5,5) + d2qk*dk5*dk5
            sderv2(6,5) = sderv2(6,5) + d2qk*dk6*dk5
            sderv2(6,6) = sderv2(6,6) + d2qk*dk6*dk6
          elseif (ndim.eq.2) then
            dk1 = dqds(1,k)
            dk2 = dqds(2,k)
            dk3 = dqds(3,k)
            sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
            sderv2(2,1) = sderv2(2,1) + d2qk*dk2*dk1
            sderv2(3,1) = sderv2(3,1) + d2qk*dk3*dk1
            sderv2(2,2) = sderv2(2,2) + d2qk*dk2*dk2
            sderv2(3,2) = sderv2(3,2) + d2qk*dk3*dk2
            sderv2(3,3) = sderv2(3,3) + d2qk*dk3*dk3
          elseif (ndim.eq.1) then
            dk1 = dqds(1,k)
            sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
          endif
        endif
      endif
    endif
    if (lstr) then
!
!  Mix strain-internal contribution
!
      if (abs(d2i2s).gt.1.0d-8) then
        do kl = 1,nstrains
          derv3(ix,kl) = derv3(ix,kl) - d2i2s*dqds(kl,i)*dikx
          derv3(iy,kl) = derv3(iy,kl) - d2i2s*dqds(kl,i)*diky
          derv3(iz,kl) = derv3(iz,kl) - d2i2s*dqds(kl,i)*dikz
        enddo
      endif
    endif
!
!  End of loop over k
!
  enddo
!
!  Strain - charge derivatives terms
!
  if (lstr) then
!
!  Mixed strain-internal term
!
!  d2E/d(alpha).dq x dq/d(strain)
!
!  NB: Exclude i = j below as it doesn't cancel as it should
!
    if (i.ne.j) then
      do kl = 1,nstrains
        derv3(ix,kl) = derv3(ix,kl) - d1ix*dqds(kl,i)
        derv3(iy,kl) = derv3(iy,kl) - d1iy*dqds(kl,i)
        derv3(iz,kl) = derv3(iz,kl) - d1iz*dqds(kl,i)
        derv3(ix,kl) = derv3(ix,kl) - d1jx*dqds(kl,j)
        derv3(iy,kl) = derv3(iy,kl) - d1jy*dqds(kl,j)
        derv3(iz,kl) = derv3(iz,kl) - d1jz*dqds(kl,j)
      enddo
    endif
!
!  Strain-strain terms for charge derivatives
!
!  NB: Only add if j is less than or equal to i to avoid double counting
!
    if (j.le.i) then
      if (ndim.eq.3) then
        d1si = dsi(1)
        d2si = dsi(2)
        d3si = dsi(3)
        d4si = dsi(4)
        d5si = dsi(5)
        d6si = dsi(6)
        d1sj = dsj(1)
        d2sj = dsj(2)
        d3sj = dsj(3)
        d4sj = dsj(4)
        d5sj = dsj(5)
        d6sj = dsj(6)
        sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1)
        sderv2(2,1) = sderv2(2,1) + d1si*dis2 + d1sj*djs2
        sderv2(2,1) = sderv2(2,1) + d2si*dis1 + d2sj*djs1
        sderv2(3,1) = sderv2(3,1) + d1si*dis3 + d1sj*djs3
        sderv2(3,1) = sderv2(3,1) + d3si*dis1 + d3sj*djs1
        sderv2(4,1) = sderv2(4,1) + d1si*dis4 + d1sj*djs4
        sderv2(4,1) = sderv2(4,1) + d4si*dis1 + d4sj*djs1
        sderv2(5,1) = sderv2(5,1) + d1si*dis5 + d1sj*djs5
        sderv2(5,1) = sderv2(5,1) + d5si*dis1 + d5sj*djs1
        sderv2(6,1) = sderv2(6,1) + d1si*dis6 + d1sj*djs6
        sderv2(6,1) = sderv2(6,1) + d6si*dis1 + d6sj*djs1
        sderv2(2,2) = sderv2(2,2) + 2.0_dp*(d2si*dis2+ d2sj*djs2)
        sderv2(3,2) = sderv2(3,2) + d2si*dis3 + d2sj*djs3
        sderv2(3,2) = sderv2(3,2) + d3si*dis2 + d3sj*djs2
        sderv2(4,2) = sderv2(4,2) + d2si*dis4 + d2sj*djs4
        sderv2(4,2) = sderv2(4,2) + d4si*dis2 + d4sj*djs2
        sderv2(5,2) = sderv2(5,2) + d2si*dis5 + d2sj*djs5
        sderv2(5,2) = sderv2(5,2) + d5si*dis2 + d5sj*djs2
        sderv2(6,2) = sderv2(6,2) + d2si*dis6 + d2sj*djs6
        sderv2(6,2) = sderv2(6,2) + d6si*dis2 + d6sj*djs2
        sderv2(3,3) = sderv2(3,3) + 2.0_dp*(d3si*dis3 + d3sj*djs3)
        sderv2(4,3) = sderv2(4,3) + d3si*dis4 + d3sj*djs4
        sderv2(4,3) = sderv2(4,3) + d4si*dis3 + d4sj*djs3
        sderv2(5,3) = sderv2(5,3) + d3si*dis5 + d3sj*djs5
        sderv2(5,3) = sderv2(5,3) + d5si*dis3 + d5sj*djs3
        sderv2(6,3) = sderv2(6,3) + d3si*dis6 + d3sj*djs6
        sderv2(6,3) = sderv2(6,3) + d6si*dis3 + d6sj*djs3
        sderv2(4,4) = sderv2(4,4) + 2.0_dp*(d4si*dis4 + d4sj*djs4)
        sderv2(5,4) = sderv2(5,4) + d4si*dis5 + d4sj*djs5
        sderv2(5,4) = sderv2(5,4) + d5si*dis4 + d5sj*djs4
        sderv2(6,4) = sderv2(6,4) + d4si*dis6 + d4sj*djs6
        sderv2(6,4) = sderv2(6,4) + d6si*dis4 + d6sj*djs4
        sderv2(5,5) = sderv2(5,5) + 2.0_dp*(d5si*dis5 + d5sj*djs5)
        sderv2(6,5) = sderv2(6,5) + d5si*dis6 + d5sj*djs6
        sderv2(6,5) = sderv2(6,5) + d6si*dis5 + d6sj*djs5
        sderv2(6,6) = sderv2(6,6) + 2.0_dp*(d6si*dis6 + d6sj*djs6)
!
        sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
        sderv2(2,1) = sderv2(2,1) + d2ijs*(dis1*djs2 + dis2*djs1)
        sderv2(3,1) = sderv2(3,1) + d2ijs*(dis1*djs3 + dis3*djs1)
        sderv2(4,1) = sderv2(4,1) + d2ijs*(dis1*djs4 + dis4*djs1)
        sderv2(5,1) = sderv2(5,1) + d2ijs*(dis1*djs5 + dis5*djs1)
        sderv2(6,1) = sderv2(6,1) + d2ijs*(dis1*djs6 + dis6*djs1)
        sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2ijs*dis2*djs2
        sderv2(3,2) = sderv2(3,2) + d2ijs*(dis2*djs3 + dis3*djs2)
        sderv2(4,2) = sderv2(4,2) + d2ijs*(dis2*djs4 + dis4*djs2)
        sderv2(5,2) = sderv2(5,2) + d2ijs*(dis2*djs5 + dis5*djs2)
        sderv2(6,2) = sderv2(6,2) + d2ijs*(dis2*djs6 + dis6*djs2)
        sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2ijs*dis3*djs3
        sderv2(4,3) = sderv2(4,3) + d2ijs*(dis3*djs4 + dis4*djs3)
        sderv2(5,3) = sderv2(5,3) + d2ijs*(dis3*djs5 + dis5*djs3)
        sderv2(6,3) = sderv2(6,3) + d2ijs*(dis3*djs6 + dis6*djs3)
        sderv2(4,4) = sderv2(4,4) + 2.0_dp*d2ijs*dis4*djs4
        sderv2(5,4) = sderv2(5,4) + d2ijs*(dis4*djs5 + dis5*djs4)
        sderv2(6,4) = sderv2(6,4) + d2ijs*(dis4*djs6 + dis6*djs4)
        sderv2(5,5) = sderv2(5,5) + 2.0_dp*d2ijs*dis5*djs5
        sderv2(6,5) = sderv2(6,5) + d2ijs*(dis5*djs6 + dis6*djs5)
        sderv2(6,6) = sderv2(6,6) + 2.0_dp*d2ijs*dis6*djs6
        if (abs(d2i2s).gt.1.0d-8) then
          sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
          sderv2(2,1) = sderv2(2,1) + d2i2s*dis1*dis2
          sderv2(3,1) = sderv2(3,1) + d2i2s*dis1*dis3
          sderv2(4,1) = sderv2(4,1) + d2i2s*dis1*dis4
          sderv2(5,1) = sderv2(5,1) + d2i2s*dis1*dis5
          sderv2(6,1) = sderv2(6,1) + d2i2s*dis1*dis6
          sderv2(2,2) = sderv2(2,2) + d2i2s*dis2*dis2
          sderv2(3,2) = sderv2(3,2) + d2i2s*dis2*dis3
          sderv2(4,2) = sderv2(4,2) + d2i2s*dis2*dis4
          sderv2(5,2) = sderv2(5,2) + d2i2s*dis2*dis5
          sderv2(6,2) = sderv2(6,2) + d2i2s*dis2*dis6
          sderv2(3,3) = sderv2(3,3) + d2i2s*dis3*dis3
          sderv2(4,3) = sderv2(4,3) + d2i2s*dis3*dis4
          sderv2(5,3) = sderv2(5,3) + d2i2s*dis3*dis5
          sderv2(6,3) = sderv2(6,3) + d2i2s*dis3*dis6
          sderv2(4,4) = sderv2(4,4) + d2i2s*dis4*dis4
          sderv2(5,4) = sderv2(5,4) + d2i2s*dis4*dis5
          sderv2(6,4) = sderv2(6,4) + d2i2s*dis4*dis6
          sderv2(5,5) = sderv2(5,5) + d2i2s*dis5*dis5
          sderv2(6,5) = sderv2(6,5) + d2i2s*dis5*dis6
          sderv2(6,6) = sderv2(6,6) + d2i2s*dis6*dis6
        endif
        if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
          sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
          sderv2(2,1) = sderv2(2,1) + d2j2s*djs1*djs2
          sderv2(3,1) = sderv2(3,1) + d2j2s*djs1*djs3
          sderv2(4,1) = sderv2(4,1) + d2j2s*djs1*djs4
          sderv2(5,1) = sderv2(5,1) + d2j2s*djs1*djs5
          sderv2(6,1) = sderv2(6,1) + d2j2s*djs1*djs6
          sderv2(2,2) = sderv2(2,2) + d2j2s*djs2*djs2
          sderv2(3,2) = sderv2(3,2) + d2j2s*djs2*djs3
          sderv2(4,2) = sderv2(4,2) + d2j2s*djs2*djs4
          sderv2(5,2) = sderv2(5,2) + d2j2s*djs2*djs5
          sderv2(6,2) = sderv2(6,2) + d2j2s*djs2*djs6
          sderv2(3,3) = sderv2(3,3) + d2j2s*djs3*djs3
          sderv2(4,3) = sderv2(4,3) + d2j2s*djs3*djs4
          sderv2(5,3) = sderv2(5,3) + d2j2s*djs3*djs5
          sderv2(6,3) = sderv2(6,3) + d2j2s*djs3*djs6
          sderv2(4,4) = sderv2(4,4) + d2j2s*djs4*djs4
          sderv2(5,4) = sderv2(5,4) + d2j2s*djs4*djs5
          sderv2(6,4) = sderv2(6,4) + d2j2s*djs4*djs6
          sderv2(5,5) = sderv2(5,5) + d2j2s*djs5*djs5
          sderv2(6,5) = sderv2(6,5) + d2j2s*djs5*djs6
          sderv2(6,6) = sderv2(6,6) + d2j2s*djs6*djs6
        endif
      elseif (ndim.eq.2) then
        d1si = dsi(1)
        d2si = dsi(2)
        d3si = dsi(3)
        d1sj = dsj(1)
        d2sj = dsj(2)
        d3sj = dsj(3)
        sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1)
        sderv2(2,1) = sderv2(2,1) + d1si*dis2 + d1sj*djs2
        sderv2(2,1) = sderv2(2,1) + d2si*dis1 + d2sj*djs1
        sderv2(3,1) = sderv2(3,1) + d1si*dis3 + d1sj*djs3
        sderv2(3,1) = sderv2(3,1) + d3si*dis1 + d3sj*djs1
        sderv2(2,2) = sderv2(2,2) + 2.0_dp*(d2si*dis2 + d2sj*djs2)
        sderv2(3,2) = sderv2(3,2) + d2si*dis3 + d2sj*djs3
        sderv2(3,2) = sderv2(3,2) + d3si*dis2 + d3sj*djs2
        sderv2(3,3) = sderv2(3,3) + 2.0_dp*(d3si*dis3 + d3sj*djs3)
!
        sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
        sderv2(2,1) = sderv2(2,1) + d2ijs*(dis1*djs2 + dis2*djs1)
        sderv2(3,1) = sderv2(3,1) + d2ijs*(dis1*djs3 + dis3*djs1)
        sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2ijs*dis2*djs2
        sderv2(3,2) = sderv2(3,2) + d2ijs*(dis2*djs3 + dis3*djs2)
        sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2ijs*dis3*djs3
        if (abs(d2i2s).gt.1.0d-8) then
          sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
          sderv2(2,1) = sderv2(2,1) + d2i2s*dis1*dis2
          sderv2(3,1) = sderv2(3,1) + d2i2s*dis1*dis3
          sderv2(2,2) = sderv2(2,2) + d2i2s*dis2*dis2
          sderv2(3,2) = sderv2(3,2) + d2i2s*dis2*dis3
          sderv2(3,3) = sderv2(3,3) + d2i2s*dis3*dis3
        endif
        if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
          sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
          sderv2(2,1) = sderv2(2,1) + d2j2s*djs1*djs2
          sderv2(3,1) = sderv2(3,1) + d2j2s*djs1*djs3
          sderv2(2,2) = sderv2(2,2) + d2j2s*djs2*djs2
          sderv2(3,2) = sderv2(3,2) + d2j2s*djs2*djs3
          sderv2(3,3) = sderv2(3,3) + d2j2s*djs3*djs3
        endif
      elseif (ndim.eq.1) then
        d1si = dsi(1)
        d1sj = dsj(1)
        sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1)
        sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
        if (abs(d2i2s).gt.1.0d-8) then
          sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
        endif
        if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
          sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
        endif
      endif
    endif
  endif
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2charged')
#endif
!
  return
  end
