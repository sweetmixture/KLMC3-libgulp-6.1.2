program fsub_main
!  call fsub("ABDFAFD")
end program fsub_main

subroutine fsub(string) bind(C,name='ffoo')

  use iso_c_binding

  character(kind=c_char), dimension(*), intent(in) :: string
  integer :: len, i

  character(len=64) tmp

  character(len=256) :: env_value
  character(len=30) :: env_variable_name



  tmp = ""
  len = 0
  do
    if(string(len+1) == C_NULL_CHAR) exit
    len = len + 1
    ! write(*,'(A)') string(len)
    tmp(len:len) = string(len)
  end do
  write(*,'(A)') tmp
  write(*,'(A,I4)') "length: ", len

  write (*,'(a)') "this is fortran subroutine"

  ! Specify the name of the environment variable you want to retrieve
  env_variable_name = "HOME"
  ! env_variable_name = "asdfasdf"
  ! Call the GETENV intrinsic subroutine to get the value of the environment
  ! variable
  call GETENV(trim(env_variable_name), env_value)
  ! Print the value of the environment variable
  write(*, '(A, A)') "The value of ", trim(env_variable_name), " is ", trim(env_value)



end subroutine
