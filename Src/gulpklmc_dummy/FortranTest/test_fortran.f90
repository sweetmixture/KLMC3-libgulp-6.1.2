program string_test

  character(len=256) :: path = "/work/e05/e05/wkjee/Software/gulp-6.1.2/Src/Custom/path_test"
  character(len=256) :: var

  CHARACTER (LEN=16)   ENVDATA
  CALL GETENV('HOME', VALUE=ENVDATA)

  write(var,'(A)') path
  write(var,'(A)') trim(path)//"\n"

  write(*,'(A)') var

end program
