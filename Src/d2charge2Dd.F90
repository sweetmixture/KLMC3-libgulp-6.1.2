  subroutine d2charge2Dd(iloc,i,j,nor,ix,iy,iz,jx,jy,jz,dei,dej, &
                         d1xi,d1yi,d1zi,d1xj,d1yj,d1zj,ds1i,ds1j,d2i2,d2ij,d2j2, &
                         d2self,d2trm1dij,d2trm1diz,d2trm1djz,dtrm1di,dtrm1dj)
!
!  Calculates the contribution to the second derivative matrices
!  due to charge derivatives from EEM/QEq for the 2-D case. Only
!  needed for reciprocal space part.
!  Parallel distributed memory version.
!
!   1/22 Created from d2charge2D and d2charged
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
  use m_strain,       only : strainddetds, straindet
  use optimisation
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: iloc
  integer(i4), intent(in)  :: i
  integer(i4), intent(in)  :: ix
  integer(i4), intent(in)  :: iy
  integer(i4), intent(in)  :: iz
  integer(i4), intent(in)  :: j
  integer(i4), intent(in)  :: jx
  integer(i4), intent(in)  :: jy
  integer(i4), intent(in)  :: jz
  integer(i4), intent(in)  :: nor
  real(dp),    intent(in)  :: d1xi          ! d2E/dqi.dx
  real(dp),    intent(in)  :: d1yi          ! d2E/dqi.dy
  real(dp),    intent(in)  :: d1zi          ! d2E/dqi.dz
  real(dp),    intent(in)  :: d1xj          ! d2E/dqj.dx
  real(dp),    intent(in)  :: d1yj          ! d2E/dqj.dy
  real(dp),    intent(in)  :: d1zj          ! d2E/dqj.dz
  real(dp),    intent(in)  :: d2i2(*)
  real(dp),    intent(in)  :: d2ij(*)
  real(dp),    intent(in)  :: d2j2(*)
  real(dp),    intent(in)  :: d2self
  real(dp),    intent(in)  :: d2trm1dij
  real(dp),    intent(in)  :: d2trm1diz
  real(dp),    intent(in)  :: d2trm1djz
  real(dp),    intent(in)  :: dei(*)
  real(dp),    intent(in)  :: dej(*)
  real(dp),    intent(in)  :: ds1i(*)       ! d2E/dqi.d(epsilon)
  real(dp),    intent(in)  :: ds1j(*)       ! d2E/dqj.d(epsilon)
  real(dp),    intent(in)  :: dtrm1di
  real(dp),    intent(in)  :: dtrm1dj
!
!  Local variables
!
  integer(i4)              :: ind
  integer(i4)              :: iv
  integer(i4)              :: ixf
  integer(i4)              :: iyf
  integer(i4)              :: izf
  integer(i4)              :: k
  integer(i4)              :: kk
  integer(i4)              :: kl
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
  real(dp)                 :: d1ix
  real(dp)                 :: d1iy
  real(dp)                 :: d1iz
  real(dp)                 :: d1jx
  real(dp)                 :: d1jy
  real(dp)                 :: d1jz
  real(dp)                 :: d2i2s
  real(dp)                 :: d2ijs
  real(dp)                 :: d2j2s
  real(dp)                 :: d1si   
  real(dp)                 :: d2si   
  real(dp)                 :: d3si
  real(dp)                 :: d1sj   
  real(dp)                 :: d2sj   
  real(dp)                 :: d3sj
  real(dp)                 :: deisum
  real(dp)                 :: dejsum
  real(dp)                 :: dikx   
  real(dp)                 :: diky
  real(dp)                 :: dikz
  real(dp)                 :: dis1 
  real(dp)                 :: dis2 
  real(dp)                 :: dis3
  real(dp)                 :: djs1  
  real(dp)                 :: djs2  
  real(dp)                 :: djs3
  real(dp)                 :: dsi(3)
  real(dp)                 :: dsi2(3)
  real(dp)                 :: dsj(3)
  real(dp)                 :: dsj2(3)
  real(dp)                 :: g_cpu_time
  real(dp)                 :: time1
  real(dp)                 :: time2
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charge2Dd')
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
!
  d1iz = d1iz + d2trm1diz
  d1jz = d1jz + d2trm1djz
  d2ijs = d2ijs + d2trm1dij
!
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
    dis1 = dqds(1,i)
    dis2 = dqds(2,i)
    dis3 = dqds(3,i)
    djs1 = dqds(1,j)
    djs2 = dqds(2,j)
    djs3 = dqds(3,j)
    do kk = 1,nstrains
      dsi(kk) = ds1i(kk)
      dsj(kk) = ds1j(kk)
      dsi2(kk) = dsi(kk)
      dsj2(kk) = dsj(kk)
    enddo
    if (lfinitestrain) then
      do kk = 1,nstrains
        do iv = 1,nor
          dsi(kk) = dsi(kk) - dei(iv)*strainddetds(kk)*straindet
          dsj(kk) = dsj(kk) - dej(iv)*strainddetds(kk)*straindet
        enddo
        dsi(kk) = dsi(kk) - dtrm1di*strainddetds(kk)*straindet
        dsj(kk) = dsj(kk) - dtrm1dj*strainddetds(kk)*straindet
      enddo
    else
      do iv = 1,nor
        dsi(1) = dsi(1) - dei(iv)
        dsj(1) = dsj(1) - dej(iv)
        dsi(2) = dsi(2) - dei(iv)
        dsj(2) = dsj(2) - dej(iv)
      enddo
      dsi(1) = dsi(1) - dtrm1di
      dsj(1) = dsj(1) - dtrm1dj
      dsi(2) = dsi(2) - dtrm1di
      dsj(2) = dsj(2) - dtrm1dj
    endif
    if (leem.or.lgfnff) then
      do kk = 1,nstrains
        dsi2(kk) = dsi(kk)
        dsj2(kk) = dsj(kk)
      enddo
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
      if (i.gt.k) then
        derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,iloc)
        derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,iloc)
        derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,iloc)
        derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,iloc)
        derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,iloc)
        derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,iloc)
        derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,iloc)
        derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,iloc)
        derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,iloc)
      elseif (i.lt.k) then
        derv2(ix,kx) = derv2(ix,kx) + deisum*d2qdxyz2(mnxx,iloc)
        derv2(iy,kx) = derv2(iy,kx) + deisum*d2qdxyz2(mnxy,iloc)
        derv2(iz,kx) = derv2(iz,kx) + deisum*d2qdxyz2(mnxz,iloc)
        derv2(ix,ky) = derv2(ix,ky) + deisum*d2qdxyz2(mnxy,iloc)
        derv2(iy,ky) = derv2(iy,ky) + deisum*d2qdxyz2(mnyy,iloc)
        derv2(iz,ky) = derv2(iz,ky) + deisum*d2qdxyz2(mnyz,iloc)
        derv2(ix,kz) = derv2(ix,kz) + deisum*d2qdxyz2(mnxz,iloc)
        derv2(iy,kz) = derv2(iy,kz) + deisum*d2qdxyz2(mnyz,iloc)
        derv2(iz,kz) = derv2(iz,kz) + deisum*d2qdxyz2(mnzz,iloc)
      endif
!
!  Strain - strain contribution
!
      if (lstr.and.ndim.gt.0) then
        ind = 0
        do kk = 1,nstrains
          do kl = 1,kk-1
            ind = ind + 1
            sderv2(kk,kl) = sderv2(kk,kl) + deisum*d2qds2(ind,iloc)
          enddo
          ind = ind + 1
          sderv2(kk,kk) = sderv2(kk,kk) + deisum*d2qds2(ind,iloc)
!
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
    dikx = dqdxyz(kx,i)
    diky = dqdxyz(ky,i)
    dikz = dqdxyz(kz,i)
!
    if (lstr) then
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
    do kl = 1,nstrains
      derv3(ix,kl) = derv3(ix,kl) - d1ix*dqds(kl,i)
      derv3(iy,kl) = derv3(iy,kl) - d1iy*dqds(kl,i)
      derv3(iz,kl) = derv3(iz,kl) - d1iz*dqds(kl,i)
      derv3(ix,kl) = derv3(ix,kl) - d1jx*dqds(kl,j)
      derv3(iy,kl) = derv3(iy,kl) - d1jy*dqds(kl,j)
      derv3(iz,kl) = derv3(iz,kl) - d1jz*dqds(kl,j)
    enddo
!
!  Strain-strain terms for charge derivatives
!
!  NB: Only add if j is less than or equal to i to avoid double counting
!
    if (j.le.i) then
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
!
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
    endif
  endif
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2charge2Dd')
#endif
!
  return
  end
