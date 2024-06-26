subroutine subX

  implicit none

  integer, save :: ia = 0
  integer       :: ib = 0

  write(*,'(A,I4,I4)') "ia / ib : ", ia ,ib

  ia = ia + 1
  ib = ib + 1

  return
end subroutine subX
