  subroutine mpfinish
!
!  Closes down MPI and stops program
!

#ifdef MPI
  integer ierr

#ifndef KLMC
  call MPI_Finalize(ierr)
! wkjee - check if MPI_Finalized()
! wkjee - do not finalize mpi - if KLMC is running
#ifdef KLMC_DEBUG_MPFINISH
  write(*,'(A,I4)') "in mpfinish: MPI_Finalized(): ", ierr
#endif

#else

#ifdef KLMC_DEBUG_MPFINISH
  write(*,'(A)') "in mpfinish: MPI_Finalise() is not called"
#endif

#endif 

#endif

#ifndef KLMC
! wkjee - do not stop the program
  stop
#else
  return
#endif

  end
