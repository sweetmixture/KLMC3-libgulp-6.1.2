program getenv_example
    implicit none

    character(len=256) :: env_value
    character(len=30) :: env_variable_name

    ! Specify the name of the environment variable you want to retrieve
    env_variable_name = "HOME"
    ! env_variable_name = "asdfasdf"

    ! Call the GETENV intrinsic subroutine to get the value of the environment
    ! variable
    call GETENV(trim(env_variable_name), env_value)

    ! Print the value of the environment variable
    write(*, '(A, A)') "The value of ", trim(env_variable_name), " is ", trim(env_value)


    call subX
    call subX
    call subX

end program getenv_example
