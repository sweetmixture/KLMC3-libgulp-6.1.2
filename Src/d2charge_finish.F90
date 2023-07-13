  subroutine d2charge_finish
!
!  Calculates the second contribution to the second derivative matrices
!  due to charge derivatives from variable charge models using the pre-summed terms
!
!  NB: It is assumed that all atoms are included in the 2nd derivatives
!      as lfreeze is incompatible with variable charges
!
!   7/21 Created from d2charge
!   7/21 QEq self term corrections added
!   1/22 Symmetrisation of derv2 added
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
  use control
  use current
  use derivatives
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
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: kx
  integer(i4)                                  :: ky
  integer(i4)                                  :: kz
  integer(i4)                                  :: l
  integer(i4)                                  :: lx
  integer(i4)                                  :: ly
  integer(i4)                                  :: lz
  integer(i4)                                  :: status
  real(dp)                                     :: d2i2
  real(dp)                                     :: d2ij
  real(dp)                                     :: dikx
  real(dp)                                     :: diky
  real(dp)                                     :: dikz
  real(dp)                                     :: dilx
  real(dp)                                     :: dily
  real(dp)                                     :: dilz
  real(dp),    dimension(:), allocatable, save :: dEdxyzdq
  real(dp)                                     :: dqikx
  real(dp)                                     :: dqiky
  real(dp)                                     :: dqikz
  real(dp)                                     :: dqilx
  real(dp)                                     :: dqily
  real(dp)                                     :: dqilz
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: time1
  real(dp)                                     :: time2
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charge_finish')
#endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(dEdxyzdq(3*numat),stat=status)
  if (status/=0) call outofmemory('d2charge_finish','dEdxyzdq')
!***************
!  First pass  *
!***************
!
!  Loop over i-j
!
  do i = 1,numat
    dEdxyzdq(1:3*numat) = 0.0_dp
    do j = 1,i
      kx = - 2
      ky = - 1
      kz =   0
      d2ij = d2edqdq(j,i)
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
    do k = 1,numat
      kx = kx + 3
      ky = ky + 3
      kz = kz + 3
!
      dikx = dEdxyzdq(kx)
      diky = dEdxyzdq(ky)
      dikz = dEdxyzdq(kz)
!
      dqikx = dqdxyz(kx,i)
      dqiky = dqdxyz(ky,i)
      dqikz = dqdxyz(kz,i)
!
      lx = - 2
      ly = - 1
      lz =   0
!
      do l = 1,k-1
        lx = lx + 3
        ly = ly + 3
        lz = lz + 3
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
        if (l.lt.k) then
          derv2(lx,kx) = derv2(lx,kx) + dqilx*dikx + dqikx*dilx
          derv2(ly,kx) = derv2(ly,kx) + dqilx*diky + dqiky*dilx
          derv2(lz,kx) = derv2(lz,kx) + dqilx*dikz + dqikz*dilx
          derv2(lx,ky) = derv2(lx,ky) + dqily*dikx + dqikx*dily
          derv2(ly,ky) = derv2(ly,ky) + dqily*diky + dqiky*dily
          derv2(lz,ky) = derv2(lz,ky) + dqily*dikz + dqikz*dily
          derv2(lx,kz) = derv2(lx,kz) + dqilz*dikx + dqikx*dilz
          derv2(ly,kz) = derv2(ly,kz) + dqilz*diky + dqiky*dilz
          derv2(lz,kz) = derv2(lz,kz) + dqilz*dikz + dqikz*dilz
        else
          derv2(kx,lx) = derv2(kx,lx) + dqilx*dikx + dqikx*dilx
          derv2(ky,lx) = derv2(ky,lx) + dqily*dikx + dqikx*dily
          derv2(kz,lx) = derv2(kz,lx) + dqilz*dikx + dqikx*dilz
          derv2(kx,ly) = derv2(kx,ly) + dqilx*diky + dqiky*dilx
          derv2(ky,ly) = derv2(ky,ly) + dqily*diky + dqiky*dily
          derv2(kz,ly) = derv2(kz,ly) + dqilz*diky + dqiky*dilz
          derv2(kx,lz) = derv2(kx,lz) + dqilx*dikz + dqikz*dilx
          derv2(ky,lz) = derv2(ky,lz) + dqily*dikz + dqikz*dily
          derv2(kz,lz) = derv2(kz,lz) + dqilz*dikz + dqikz*dilz
        endif
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
      do k = 1,numat
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
        do l = 1,k-1
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
!
!  Free local memory
!
  deallocate(dEdxyzdq,stat=status)
  if (status/=0) call deallocate_error('d2charge_finish','dEdxyzdq')
!
!  Symmetrise second derivative matrix after applying charge derivatives
!
  do i = 2,3*numat
    do j = 1,i-1
      derv2(i,j) = derv2(j,i)
    enddo
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2charge_finish')
#endif
!
  return
  end
