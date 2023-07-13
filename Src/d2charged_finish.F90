  subroutine d2charged_finish
!
!  Calculates the second contribution to the second derivative matrices
!  due to charge derivatives from variable charge models using the pre-summed terms
!  This version is for distributed memory parallel second derivatives.
!
!  NB: It is assumed that all atoms are included in the 2nd derivatives
!      as lfreeze is incompatible with variable charges
!
!   1/22 Created from d2charge_finish
!   1/22 dedqc and d2edqc added for variable charge second derivatives
!   6/22 Modified to avoid gfortran warning
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
!  Julian Gale, CIC, Curtin University, June 2022
!
  use control
  use current
  use derivatives
  use parallel
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloc
  integer(i4)                                  :: ind
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixf
  integer(i4)                                  :: iyf
  integer(i4)                                  :: izf
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: kx
  integer(i4)                                  :: ky
  integer(i4)                                  :: kz
  integer(i4)                                  :: kxf
  integer(i4)                                  :: kyf
  integer(i4)                                  :: kzf
  integer(i4)                                  :: l
  integer(i4)                                  :: lx
  integer(i4)                                  :: ly
  integer(i4)                                  :: lz
  integer(i4)                                  :: status
  real(dp)                                     :: d1ix
  real(dp)                                     :: d1iy
  real(dp)                                     :: d1iz
  real(dp)                                     :: d1kx
  real(dp)                                     :: d1ky
  real(dp)                                     :: d1kz
  real(dp)                                     :: d2i2
  real(dp)                                     :: d2ij
  real(dp)                                     :: dikx
  real(dp)                                     :: diky
  real(dp)                                     :: dikz
  real(dp)                                     :: dilx
  real(dp)                                     :: dily
  real(dp)                                     :: dilz
  real(dp)                                     :: djix
  real(dp)                                     :: djiy
  real(dp)                                     :: djiz
  real(dp)                                     :: dkix
  real(dp)                                     :: dkiy
  real(dp)                                     :: dkiz
  real(dp),    dimension(:), allocatable, save :: dEdxyzdq
  real(dp)                                     :: dqikx
  real(dp)                                     :: dqiky
  real(dp)                                     :: dqikz
  real(dp)                                     :: dqilx
  real(dp)                                     :: dqily
  real(dp)                                     :: dqilz
  real(dp),    dimension(:), allocatable, save :: dtmp
  real(dp),    dimension(:), allocatable, save :: dtmp2
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: time1
  real(dp)                                     :: time2
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charged_finish')
#endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(dEdxyzdq(3*numat),stat=status)
  if (status/=0) call outofmemory('d2charged_finish','dEdxyzdq')
  allocate(dtmp(3*numat),stat=status)
  if (status/=0) call outofmemory('d2charged_finish','dtmp')
  allocate(dtmp2(3*numat),stat=status)
  if (status/=0) call outofmemory('d2charged_finish','dtmp2')
!
!  Sum d2edq2 over all nodes
!
  dtmp(1:numat) = d2edq2(1:numat)
  call sumall(dtmp,d2edq2,numat,"d2charged_finish","d2edq2") 
!
!  Sum dedqc over all nodes
!
  ind = 0
  do i = 1,numat
    do j = 1,3
      ind = ind + 1
      dtmp(ind) = dedqc(j,i)
    enddo
  enddo
  call sumall(dtmp,dedqc,3_i4*numat,"d2charged_finish","dedqc") 
!
!  Sum d2edqc over all nodes
!
  do i = 1,numat
    dtmp(1:3*numat) = d2edqc(1:3*numat,i)
    call sumall(dtmp,dtmp2,3_i4*numat,"d2charged_finish","d2edqc") 
    d2edqc(1:3*numat,i) = dtmp2(1:3*numat)
  enddo
  if (lstr) then
!
!  Sum ds2g / ds2gs / d2s2gs over all nodes
!
    do kl = 1,nstrains
      dtmp(1:numat) = ds2g(kl,1:numat)
      dtmp(numat+1:2*numat) = ds2gs(kl,1:numat)
      dtmp(2*numat+1:3*numat) = d2s2gs(kl,1:numat)
!
      call sumall(dtmp,dtmp2,3_i4*numat,"d2charged_finish","ds2g") 
!
      ds2g(kl,1:numat) = dtmp2(1:numat)
      ds2gs(kl,1:numat) = dtmp2(numat+1:2*numat)
      d2s2gs(kl,1:numat) = dtmp2(2*numat+1:3*numat)
    enddo
  endif
!***************
!  First pass  *
!***************
!
!  Loop over i-j
!
  do i = 1,numat
    dEdxyzdq(1:3*numat) = 0.0_dp
    iloc = atom2local(i)
    if (iloc.gt.0) then
      dtmp(1:numat) = d2edqdq(1:numat,iloc)
    endif
    call sendall(dtmp,numat,atom2node(i),"d2charged_finish","d2edqdq")
    call mpbarrier
    do j = 1,i
      kx = - 2
      ky = - 1
      kz =   0
      d2ij = dtmp(j)
      do k = 1,numat
        kx = kx + 3
        ky = ky + 3
        kz = kz + 3
        dEdxyzdq(kx) = dEdxyzdq(kx) + d2ij*dqdxyz(kx,j)
        dEdxyzdq(ky) = dEdxyzdq(ky) + d2ij*dqdxyz(ky,j)
        dEdxyzdq(kz) = dEdxyzdq(kz) + d2ij*dqdxyz(kz,j)
      enddo
    enddo
!****************
!  Second pass  *
!****************
    kx = - 2
    ky = - 1
    kz =   0
    do kk = 1,natomsonnode
      k = node2atom(kk)
      kx = kx + 3
      ky = ky + 3
      kz = kz + 3
!
      kxf = 3*(k-1) + 1
      kyf = kxf + 1
      kzf = kxf + 2
!
      dikx = dEdxyzdq(kxf)
      diky = dEdxyzdq(kyf)
      dikz = dEdxyzdq(kzf)
!
      dqikx = dqdxyz(kxf,i)
      dqiky = dqdxyz(kyf,i)
      dqikz = dqdxyz(kzf,i)
!
      lx = - 2
      ly = - 1
      lz =   0
!
      do l = 1,numat
        lx = lx + 3
        ly = ly + 3
        lz = lz + 3
!
        if (l.eq.k) cycle
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta) 
!
        dilx = dEdxyzdq(lx)
        dily = dEdxyzdq(ly)
        dilz = dEdxyzdq(lz)
!
        dqilx = dqdxyz(lx,i)
        dqily = dqdxyz(ly,i)
        dqilz = dqdxyz(lz,i)
!
        derv2(lx,kx) = derv2(lx,kx) + dqilx*dikx + dqikx*dilx
        derv2(ly,kx) = derv2(ly,kx) + dqilx*diky + dqiky*dilx
        derv2(lz,kx) = derv2(lz,kx) + dqilx*dikz + dqikz*dilx
        derv2(lx,ky) = derv2(lx,ky) + dqily*dikx + dqikx*dily
        derv2(ly,ky) = derv2(ly,ky) + dqily*diky + dqiky*dily
        derv2(lz,ky) = derv2(lz,ky) + dqily*dikz + dqikz*dily
        derv2(lx,kz) = derv2(lx,kz) + dqilz*dikx + dqikx*dilz
        derv2(ly,kz) = derv2(ly,kz) + dqilz*diky + dqiky*dilz
        derv2(lz,kz) = derv2(lz,kz) + dqilz*dikz + dqikz*dilz
!
!  End of loop over l
!
      enddo
!
!  End of loop over k
!
    enddo
!*************************
!  Extra term for QEq/H  *
!*************************
    if (abs(d2edq2(i)).gt.1.0d-8) then
      d2i2 = d2edq2(i)
!
      kx = - 2
      ky = - 1
      kz =   0
      do kk = 1,natomsonnode
        k = node2atom(kk)
        kx = kx + 3
        ky = ky + 3
        kz = kz + 3
!
        dqikx = dqdxyz(kx,i)*d2i2
        dqiky = dqdxyz(ky,i)*d2i2
        dqikz = dqdxyz(kz,i)*d2i2
!
        lx = - 2
        ly = - 1
        lz =   0
!
        do l = 1,numat
          lx = lx + 3
          ly = ly + 3
          lz = lz + 3
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta)
!
          dqilx = dqdxyz(lx,i)
          dqily = dqdxyz(ly,i)
          dqilz = dqdxyz(lz,i)
!
          derv2(lx,kx) = derv2(lx,kx) + dqilx*dqikx
          derv2(ly,kx) = derv2(ly,kx) + dqilx*dqiky
          derv2(lz,kx) = derv2(lz,kx) + dqilx*dqikz
          derv2(lx,ky) = derv2(lx,ky) + dqily*dqikx
          derv2(ly,ky) = derv2(ly,ky) + dqily*dqiky
          derv2(lz,ky) = derv2(lz,ky) + dqily*dqikz
          derv2(lx,kz) = derv2(lx,kz) + dqilz*dqikx
          derv2(ly,kz) = derv2(ly,kz) + dqilz*dqiky
          derv2(lz,kz) = derv2(lz,kz) + dqilz*dqikz
!
!  End of loop over l
!
        enddo
!
!  End of loop over k
!
      enddo
    endif
  enddo
!****************************
!  Contribution from dedqc  *
!****************************
  ix = - 2
  iy = - 1
  iz =   0
  do ii = 1,natomsonnode
    i = node2atom(ii)
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    ixf = 3*(i-1) + 1
    iyf = ixf + 1
    izf = ixf + 2
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
      dkix = dqdxyz(ixf,k)
      dkiy = dqdxyz(iyf,k)
      dkiz = dqdxyz(izf,k)
!
      d1ix = dedqc(1,i)
      d1iy = dedqc(2,i)
      d1iz = dedqc(3,i)
      d1kx = dedqc(1,k)
      d1ky = dedqc(2,k)
      d1kz = dedqc(3,k)
!
      if (i.ne.k) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : i-k
!
        if (k.gt.i) then
          derv2(kx,ix) = derv2(kx,ix) - d1ix*dikx
          derv2(ky,ix) = derv2(ky,ix) - d1ix*diky
          derv2(kz,ix) = derv2(kz,ix) - d1ix*dikz
          derv2(kx,iy) = derv2(kx,iy) - d1iy*dikx
          derv2(ky,iy) = derv2(ky,iy) - d1iy*diky
          derv2(kz,iy) = derv2(kz,iy) - d1iy*dikz
          derv2(kx,iz) = derv2(kx,iz) - d1iz*dikx
          derv2(ky,iz) = derv2(ky,iz) - d1iz*diky
          derv2(kz,iz) = derv2(kz,iz) - d1iz*dikz
!
          derv2(kx,ix) = derv2(kx,ix) - d1kx*dkix
          derv2(ky,ix) = derv2(ky,ix) - d1kx*dkiy
          derv2(kz,ix) = derv2(kz,ix) - d1kx*dkiz
          derv2(kx,iy) = derv2(kx,iy) - d1ky*dkix
          derv2(ky,iy) = derv2(ky,iy) - d1ky*dkiy
          derv2(kz,iy) = derv2(kz,iy) - d1ky*dkiz
          derv2(kx,iz) = derv2(kx,iz) - d1kz*dkix
          derv2(ky,iz) = derv2(ky,iz) - d1kz*dkiy
          derv2(kz,iz) = derv2(kz,iz) - d1kz*dkiz
        else
          derv2(kx,ix) = derv2(kx,ix) - d1ix*dikx
          derv2(ky,ix) = derv2(ky,ix) - d1iy*dikx
          derv2(kz,ix) = derv2(kz,ix) - d1iz*dikx
          derv2(kx,iy) = derv2(kx,iy) - d1ix*diky
          derv2(ky,iy) = derv2(ky,iy) - d1iy*diky
          derv2(kz,iy) = derv2(kz,iy) - d1iz*diky
          derv2(kx,iz) = derv2(kx,iz) - d1ix*dikz
          derv2(ky,iz) = derv2(ky,iz) - d1iy*dikz
          derv2(kz,iz) = derv2(kz,iz) - d1iz*dikz
!
          derv2(kx,ix) = derv2(kx,ix) - d1kx*dkix
          derv2(ky,ix) = derv2(ky,ix) - d1ky*dkix
          derv2(kz,ix) = derv2(kz,ix) - d1kz*dkix
          derv2(kx,iy) = derv2(kx,iy) - d1kx*dkiy
          derv2(ky,iy) = derv2(ky,iy) - d1ky*dkiy
          derv2(kz,iy) = derv2(kz,iy) - d1kz*dkiy
          derv2(kx,iz) = derv2(kx,iz) - d1kx*dkiz
          derv2(ky,iz) = derv2(ky,iz) - d1ky*dkiz
          derv2(kz,iz) = derv2(kz,iz) - d1kz*dkiz
        endif
      endif
    enddo
  enddo
!*****************************
!  Contribution from d2edqc  *
!*****************************
  ix = - 2
  iy = - 1
  iz =   0
  do ii = 1,natomsonnode
    i = node2atom(ii)
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
!
    ixf = 3*(i-1) + 1
    iyf = ixf + 1
    izf = ixf + 2
    do j = 1,numat
      djix = dqdxyz(ixf,j)
      djiy = dqdxyz(iyf,j)
      djiz = dqdxyz(izf,j)
      kx = -2
      ky = -1
      kz =  0
      do k = 1,numat
        kx = kx + 3
        ky = ky + 3
        kz = kz + 3
        if (k.gt.i) then
          derv2(kx,ix) = derv2(kx,ix) - d2edqc(kx,j)*djix - d2edqc(ixf,j)*dqdxyz(kx,j)
          derv2(kx,iy) = derv2(kx,iy) - d2edqc(ky,j)*djix - d2edqc(iyf,j)*dqdxyz(kx,j)
          derv2(kx,iz) = derv2(kx,iz) - d2edqc(kz,j)*djix - d2edqc(izf,j)*dqdxyz(kx,j)
          derv2(ky,ix) = derv2(ky,ix) - d2edqc(kx,j)*djiy - d2edqc(ixf,j)*dqdxyz(ky,j)
          derv2(ky,iy) = derv2(ky,iy) - d2edqc(ky,j)*djiy - d2edqc(iyf,j)*dqdxyz(ky,j)
          derv2(ky,iz) = derv2(ky,iz) - d2edqc(kz,j)*djiy - d2edqc(izf,j)*dqdxyz(ky,j)
          derv2(kz,ix) = derv2(kz,ix) - d2edqc(kx,j)*djiz - d2edqc(ixf,j)*dqdxyz(kz,j)
          derv2(kz,iy) = derv2(kz,iy) - d2edqc(ky,j)*djiz - d2edqc(iyf,j)*dqdxyz(kz,j)
          derv2(kz,iz) = derv2(kz,iz) - d2edqc(kz,j)*djiz - d2edqc(izf,j)*dqdxyz(kz,j)
        else
          derv2(kx,ix) = derv2(kx,ix) - d2edqc(kx,j)*djix - d2edqc(ixf,j)*dqdxyz(kx,j)
          derv2(kx,iy) = derv2(kx,iy) - d2edqc(kx,j)*djiy - d2edqc(ixf,j)*dqdxyz(ky,j)
          derv2(kx,iz) = derv2(kx,iz) - d2edqc(kx,j)*djiz - d2edqc(ixf,j)*dqdxyz(kz,j)
          derv2(ky,ix) = derv2(ky,ix) - d2edqc(ky,j)*djix - d2edqc(iyf,j)*dqdxyz(kx,j)
          derv2(ky,iy) = derv2(ky,iy) - d2edqc(ky,j)*djiy - d2edqc(iyf,j)*dqdxyz(ky,j)
          derv2(ky,iz) = derv2(ky,iz) - d2edqc(ky,j)*djiz - d2edqc(iyf,j)*dqdxyz(kz,j)
          derv2(kz,ix) = derv2(kz,ix) - d2edqc(kz,j)*djix - d2edqc(izf,j)*dqdxyz(kx,j)
          derv2(kz,iy) = derv2(kz,iy) - d2edqc(kz,j)*djiy - d2edqc(izf,j)*dqdxyz(ky,j)
          derv2(kz,iz) = derv2(kz,iz) - d2edqc(kz,j)*djiz - d2edqc(izf,j)*dqdxyz(kz,j)
        endif
      enddo
    enddo
  enddo
!*****************************
!  Contributions to strains  *
!*****************************
  if (lstr) then
    ix = - 2
    iy = - 1
    iz =   0
    do ii = 1,natomsonnode
      i = node2atom(ii)
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      ixf = 3*(i-1) + 1
      iyf = ixf + 1
      izf = ixf + 2
!
      do k = 1,numat
        dkix = dqdxyz(ixf,k)
        dkiy = dqdxyz(iyf,k)
        dkiz = dqdxyz(izf,k)
!
!  Mix strain-internal contribution
!
!  d2E/d(strain).dq x dq/d(alpha)
!
        do kl = 1,nstrains
          derv3(ix,kl) = derv3(ix,kl) + ds2g(kl,k)*dkix
          derv3(iy,kl) = derv3(iy,kl) + ds2g(kl,k)*dkiy
          derv3(iz,kl) = derv3(iz,kl) + ds2g(kl,k)*dkiz
        enddo
!
!  d2E/dq.dq x dq/d(alpha) x dq/d(strain)
!
        do kl = 1,nstrains
          derv3(ix,kl) = derv3(ix,kl) + ds2gs(kl,k)*dkix
          derv3(iy,kl) = derv3(iy,kl) + ds2gs(kl,k)*dkiy
          derv3(iz,kl) = derv3(iz,kl) + ds2gs(kl,k)*dkiz
        enddo
!
        do kl = 1,nstrains
          derv3(ix,kl) = derv3(ix,kl) + d2s2gs(kl,k)*dkix
          derv3(iy,kl) = derv3(iy,kl) + d2s2gs(kl,k)*dkiy
          derv3(iz,kl) = derv3(iz,kl) + d2s2gs(kl,k)*dkiz
        enddo
      enddo
    enddo
  endif
!
!  Free local memory
!
  deallocate(dtmp2,stat=status)
  if (status/=0) call deallocate_error('d2charged_finish','dtmp2')
  deallocate(dtmp,stat=status)
  if (status/=0) call deallocate_error('d2charged_finish','dtmp')
  deallocate(dEdxyzdq,stat=status)
  if (status/=0) call deallocate_error('d2charged_finish','dEdxyzdq')
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2charged_finish')
#endif
!
  return
  end
