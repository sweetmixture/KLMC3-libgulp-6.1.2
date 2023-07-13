program file_read_example
    implicit none
    character(len=100) :: filename = "/work/e05/e05/wkjee/Software/gulp-6.1.2/Src/Custom/path_test/gulp_klmc  "
    character(len=100) :: line
    integer ierr
    logical lopen

    ! Open the file
    !open(1, file=filename, status='old', action='read', iostat=ierr)
    open(1, file=trim(adjustl(filename))//".gin",form='formatted',iostat=ierr)
    if (ierr /= 0) then
        write(*, '(A)', advance='no') 'Error opening file: ', trim(filename)
        stop
    endif

    ! Read the first line
    read(1, '(A)', iostat=ierr) line
    if (ierr /= 0) then
        write(*, '(A)', advance='no') 'Error reading file: ', trim(filename)
        close(1)
        stop
    endif

    ! Close the file
    close(1)

    ! Print the first line
    write(*, '(A)', advance='no') 'First line of file: ', trim(line)


    ! open(unit=ioin,file=trim(seedname)//".gin",form='formatted')
    ! open(unit=1,file=trim(adjustl(filename))//".gout",form='formatted')
    open(unit=6,file=trim(filename)//".gout",form='formatted')
    inquire(unit=6,opened=lopen)
    write(*,'(A,L)') "is opened : ", lopen
    write(unit=6,'(A)') "banna is white"
    close(6)
    inquire(unit=6,opened=lopen)
    write(*,'(A,L)') "is opened : ", lopen


end program file_read_example
