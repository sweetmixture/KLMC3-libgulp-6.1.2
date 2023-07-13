  program gulp
!
!  General Utility Lattice Program including
!
!     (1) Electrostatic potential calculations
!     (2) Electronegativity equalisation method
!     (3) Interatomic potential fitting
!     (4) Structure optimisation
!     (5) Phonon calculations
!     (6) Defect calculations
!     (7) Molecular dynamics
!     (8) Free energy minimisation
!     (9) Calculations on 2-D and 3-D systems
!    (10) Monte Carlo calculations
!    (11) COSMO/COSMIC solvation calculations
!    (12) PDF calculations
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
!  Copyright Curtin University 2009
!
!   3/07 Main routine now called gulpmain as a subroutine called
!        from here. 
!   5/07 Option to use seed argument added (VM)
!   7/07 Modified to remove argument checking (VM)
!  10/08 COSMIC model merged into this branch
!   3/09 Communicator passed to gulpmain to aid task farming
!   6/09 PDF calculation by Elizabeth Cope merged in
!   4/17 ChemShell interaction modified
!   6/17 Old ChemShell restored as a compile option
!
!  Julian Gale, NRI, Curtin University, June 2017
!

!
! wkjee - KLMC development 07/2023
!
  use datatypes,  only : i4
#ifdef KLMC
  use klmc
  implicit none
  include 'mpif.h'
  ! able to use 'gulp_klmc_iopath' module from modules.F90 (added)
#else 
  implicit none
#endif

  integer ierr, rank, size

  integer(i4)        :: iret
#ifdef OLDCS
  integer(i4)        :: ichemsh_qm
#else
  integer(i4)        :: ichemsh_link
#endif
  integer*4          :: MPI_comm_dummy
  
  iret = 0
#ifdef OLDCS
  ichemsh_qm = -1
#else
  ! wkjee - standard gulp case
  ichemsh_link = 0
#endif
!
!  This is a dummy communicator here - for task farming gulpmain can be called
!  with different values of MPI_comm_dummy and in this case the ichemsh_qm
!  flag should be used to ensure that MPI_init is not called again later. 
!
  MPI_comm_dummy = 1
!
! wkjee - KLMC developement 07/2023
! wkjee - not touching OLDCS original but slightly modified to add KLMC


! wkjee - MPI_Comm test

  call MPI_Init(ierr)
  MPI_comm_dummy = MPI_comm_world
  write(*,'(A,I)') "in gulp : MPI_comm_dummy = ", MPI_comm_dummy
  call MPI_Comm_rank(MPI_comm_dummy,rank,ierr)
  call MPI_Comm_size(MPI_comm_dummy,size,ierr)
 
  !if ( rank .eq. 0 ) then
  write(*,'(A,I4,A4,I4)') "in gulp (main) rank / size : ", rank, "/", size
  !end if

! wkjee - MPI_Comm test end

#ifdef OLDCS
  call gulpmain(iret, ichemsh_qm, MPI_comm_dummy)
! wkjee
#elif defined KLMC
  gulp_klmc_iopath = ""
  call gulpmain(iret, ichemsh_link, MPI_comm_dummy)
! wkjee
#else
  call gulpmain(iret, ichemsh_link, MPI_comm_dummy)
#endif

!********************
!  Finish GULP run  *
!********************
  ! wkjee
  write(*,'(A)') "in gulp: calling gulpfinish"
  call gulpfinish

! wkjee - MPI_Comm test
  write(*,'(A)') "in gulp: calling MPI_Finalize()"
  call MPI_Finalize(ierr)
! wkjee - MPI_Comm test end
  write(*,'(A,I4)') "in gulp: MPI_Finalized() ierr: ", ierr
  write(*,'(A)') "in gulp: program end -- "
  end
