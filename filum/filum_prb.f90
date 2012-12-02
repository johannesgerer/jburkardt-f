program main

!*****************************************************************************80
!
!! MAIN is the main program for FILUM_PRB.
!
!  Discussion:
!
!    FILUM_PRB tests routines from the FILUM library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILUM_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FILUM library.'
 
  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test065 ( )
  call test07 ( )
  call test08 ( )
  call test085 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test145 ( )
  call test15 ( )
  call test16 ( )
  call test165 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )

  call test20 ( )
  call test21 ( )
  call test215 ( )
  call test22 ( )
  call test225 ( )
  call test23 ( )
  call test24 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILUM_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests FILE_APPEND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4  ) ierror
  character ( len = 255 ) filename
  integer ( kind = 4  ) ios
  integer ( kind = 4  ) iunit
  integer ( kind = 4  ) nrec
  character ( len = 255 ) old_filename

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  FILE_APPEND makes a new version of a file'
  write ( *, '(a)' ) '  which is "appendable".'
  write ( *, '(a)' ) ' '

  old_filename = 'filum_prb_test.txt'
  filename = 'filum_prb_append.txt'

  call file_copy ( old_filename, filename, ierror )

  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST001'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
    return
  end if

  write ( iunit, '(a)' ) '  This is the file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '  After the first pass, it has a total'
  write ( iunit, '(a)' ) '  of 4 lines of text.'

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Current contents of file:' 
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  call file_print ( filename )

  call file_append ( filename, ierror, iunit, nrec )

  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '  This is new information that has been'
  write ( iunit, '(a)' ) '  written to the file during a separate pass.'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '  The file should now contain a total of '
  write ( iunit, '(a)' ) '  11 lines of text.'

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Contents of file after reopening in APPEND mode:'
  write ( *, '(a)' ) ' '

  call file_print ( filename )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests FILE_CHAR_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4  ) char_num
  character ( len = 255 ) :: filename = 'filum_prb_test.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  FILE_CHAR_COUNT counts the characters in a file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  call file_char_count ( filename, char_num )
  write ( *, '(a,i8)' ) '  Number of characters: ', char_num

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests FILE_COLUMN_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt 
!
  implicit none

  integer ( kind = 4  ) column_num
  character ( len = 255 ) filename

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  FILE_COLUMN_COUNT counts the columns in a file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  It is assumed that the file contains a number of lines,'
  write ( *, '(a)' ) '  with each line containing the same number of words.'
  write ( *, '(a)' ) '  The task is to determine the number of words in a line,'
  write ( *, '(a)' ) '  that is, the number of "columns" of text.'

  filename = 'filum_prb_4by5.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining the file: ' 
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'

  call file_column_count ( filename, column_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of columns: ', column_num

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests FILE_COLUMN_COUNT and FILE_COLUMN_RANGE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8  ) col_max(5)
  real ( kind = 8  ) col_min(5)
  integer ( kind = 4  ) column_num
  character ( len = 255 ) filename
  integer ( kind = 4  ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  FILE_COLUMN_COUNT counts the columns in a file.'
  write ( *, '(a)' ) '  FILE_COLUMN_RANGE finds the range of the columns.'
  write ( *, '(a)' ) ' '

  filename = 'filum_prb_columns.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining the file: ' 
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Contents of file:'
  write ( *, '(a)' ) ' '

  call file_print ( filename )

  call file_column_count ( filename, column_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of columns in file: ', column_num

  call file_column_range ( filename, column_num, col_min, col_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Col  Minimum  Maximum'
  write ( *, '(a)' ) ' '

  do i = 1, column_num
    write ( *, '(2x,i3,2f10.4)' ) i, col_min(i), col_max(i)
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests FILE_COPY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4  ) ierror
  character ( len = 255 ) :: new_filename = 'filum_prb_test.txt'
  character ( len = 255 ) :: old_filename = 'filum_prb_copy.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  FILE_COPY makes a copy of a file.'

  call file_copy ( old_filename, new_filename, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FILE_COPY returned IERROR = ', ierror
    return
  end if
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Finished copying.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input file is  "' // trim ( old_filename ) // '".'
  write ( *, '(a)' ) '  Output file is "' // trim ( new_filename ) // '".'

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests FILE_EXIST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  logical file_exist
  character ( len = 255 ) filename

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  FILE_EXIST reports whether a file "exists".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exist?  File_name'
  write ( *, '(a)' ) ' '

  filename = 'filum_prb.f90'
  write ( *, '(7x,l1,2x,a)' ) file_exist ( filename ), trim ( filename )
  filename = 'filum.f90'
  write ( *, '(7x,l1,2x,a)' ) file_exist ( filename ), trim ( filename )
  filename = 'raisin.txt'
  write ( *, '(7x,l1,2x,a)' ) file_exist ( filename ), trim ( filename )
  filename = 'make.money.fast'
  write ( *, '(7x,l1,2x,a)' ) file_exist ( filename ), trim ( filename )

  return
end
subroutine test065 ( )

!*****************************************************************************80
!
!! TEST065 tests FILE_GET_NEXT_INTEGER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) filename
  integer ( kind = 4  ) file_status
  integer ( kind = 4  ) file_unit
  integer ( kind = 4  ) ierror
  logical more
  integer ( kind = 4  ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST065'
  write ( *, '(a)' ) '  FILE_GET_NEXT_INTEGER gets the next integer'
  write ( *, '(a)' ) '  from a file.'

  filename = 'integers.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Contents of file:'
  write ( *, '(a)' ) ' '

  call file_print ( filename )

  call get_unit ( file_unit )

  open ( unit = file_unit, file = filename, status = 'old', &
    iostat = file_status )

  if ( file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
    return
  end if

  more = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Successive Integers:'
  write ( *, '(a)' ) ' '

  do

    call file_get_next_integer ( file_unit, more, value )

    if ( .not. more ) then
      exit
    end if

    write ( *, '(2x,i8)' ) value

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  End of file!'

  close ( unit = file_unit )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests FILE_GET_NEXT_WORD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) filename
  integer ( kind = 4  ) ierror
  integer ( kind = 4  ) ios
  integer ( kind = 4  ) iunit
  integer ( kind = 4  ) num_text
  integer ( kind = 4  ) num_text_old
  integer ( kind = 4  ) num_word
  character ( len = 255 ) text
  character ( len = 255 ) word

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  FILE_GET_NEXT_WORD gets the next word'
  write ( *, '(a)' ) '  from a file.'

  filename = 'filum_prb_test.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
    return
  end if

  num_text = 0
  num_word = 0
  text = ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Line    Word  ---WORD---  ---TEXT--------------'
  write ( *, '(a)' ) ' '

  do

    num_text_old = num_text

    call file_get_next_word ( iunit, word, text, num_text, ierror )

    if ( ierror /= 0 ) then
      exit
    end if

    if ( num_text == num_text_old ) then
      num_word = num_word + 1
    else
      num_word = 1
    end if

    write ( *, '(i8,2x,i8,2x,a,2x,a)' ) num_text, num_word, word(1:10), &
      trim ( text )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  End of file!'

  close ( unit = iunit )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests FILE_LINE_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) filename
  character ( len = 255 ) line
  integer ( kind = 4  ) line_index
  integer ( kind = 4  ) line_num
  integer ( kind = 4  ) seed
  integer ( kind = 4  ) test
  integer ( kind = 4  ), parameter :: test_num = 5

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  FILE_LINE_UNIFORM selects one line of a file,'
  write ( *, '(a)' ) '  uniformly at random, reading the file one time.'
  write ( *, '(a)' ) ' '

  filename = 'filum_prb_test.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining the file: ' 
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Line    Text:'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call file_line_uniform ( filename, seed, line, line_index, line_num )
    write ( *, '(2x,i8,2x,a)' ) line_index, trim ( line )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Total number of lines in file is ', line_num

  return
end
subroutine test085 ( )

!*****************************************************************************80
!
!! TEST085 tests FILE_LINE_WIDTH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4  ) file_line_width
  character ( len = 255 ) filename
  integer ( kind = 4  ) line_width

  filename = 'filum_prb_test.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085'
  write ( *, '(a)' ) '  FILE_LINE_WIDTH counts the longest line in a file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  line_width = file_line_width ( filename )

  write ( *, '(a,i8)' ) '  Longest line width: ', line_width

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests FILE_LINES_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4  ), parameter :: n = 3

  character ( len = 255 ) filename
  integer ( kind = 4  ) j
  character ( len = 255 ) line(3)
  integer ( kind = 4  ) line_index(3)
  integer ( kind = 4  ) line_num
  integer ( kind = 4  ) seed
  integer ( kind = 4  ) test
  integer ( kind = 4  ), parameter :: test_num = 20

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  FILE_LINES_UNIFORM selects one line of a file,'
  write ( *, '(a)' ) '  uniformly at random, reading the file one time.'
  write ( *, '(a)' ) ' '

  filename = 'filum_prb_test.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining the file: ' 
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Line  Text:'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    write ( *, '(a)' ) ' '

    call file_lines_uniform ( filename, seed, n, line, line_index, line_num )

    do j = 1, n
      write ( *, '(2x,i8,2x,a)' ) line_index(j), trim ( line(j) )
    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Total number of lines in file is ', line_num

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests FILENAME_APPEND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  character ( len = 3  ) append
  character ( len = 10 ) filename
  character ( len = 10 ) string(test_num)
  integer ( kind = 4 ) test

  append = 'XYZ'

  string(1) = 'bob.for'
  string(2) = 'N.B.C.D'
  string(3) = 'Naomi.'
  string(4) = 'Arthur'
  string(5) = '.amos'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  FILENAME_APPEND appends a string to the part'
  write ( *, '(a)' ) '  of a file name that precedes the extension.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Old                  New'
  write ( *, '(a)' ) '  FILENAME   APPEND   FILENAME'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    filename = string(test)
    call filename_append ( filename, append )
    write ( *, '(2x,a10,2x,a3,3x,2x,a10)' ) string(test), append, trim ( filename )
  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests FILENAME_DEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  character ( len = 20 ) filename
  character ( len = 20 ) filename_old
  integer ( kind = 4 ) j
  character ( len = 20 ) string(test_num)
  integer ( kind = 4 ) test

  string(1) = 'file???.dat'
  string(2) = 'file076.dat'
  string(3) = '3cat3.dat'
  string(4) = 'fred03.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  FILENAME_DEC decrements a string'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Input             Output'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    filename = string(test)
    write ( *, '(a)' ) ' '
    do j = 1, 4
      filename_old = filename
      call filename_dec ( filename )
      write ( *, '(2x,a,2x,a)' ) filename_old, filename
      if ( len_trim ( filename ) <= 0 ) then
        write ( *, '(a)' ) '  (Empty output string.  Quit loop!)'
        exit
      end if
    end do
  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests FILENAME_EXT_GET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  character ( len = 10 ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character ( len = 10 ) string(test_num)
  integer ( kind = 4 ) test

  string(1) = 'bob.for'
  string(2) = 'N.B.C.D'
  string(3) = 'Naomi.'
  string(4) = 'Arthur'
  string(5) = '.amos'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  FILENAME_EXT_GET finds a file extension.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FILENAME     Begin    End'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    filename = string(test)
    call filename_ext_get ( filename, i, j )
    write ( *, '(2x,a,i3,i3)' ) string(test), i, j
  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests FILENAME_EXT_SWAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  character ( len = 3  ) ext
  character ( len = 3  ) ext_test(test_num)
  character ( len = 12 ) filename
  character ( len = 12 ) string(test_num)
  integer ( kind = 4 ) test

  string(1) = 'bob.for'
  string(2) = 'bob.bob.bob'
  string(3) = 'bob.'
  string(4) = 'bob'
  string(5) = '.oops'

  ext_test(1) = 'obj'
  ext_test(2) = 'txt'
  ext_test(3) = 'yak'
  ext_test(4) = 'ps'
  ext_test(5) = 'eek'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  FILENAME_EXT_SWAP changes a file extension.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input                 Output'
  write ( *, '(a)' ) '  FILENAME     EXT     FILENAME'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    filename = string(test)
    ext = ext_test(test)
    call filename_ext_swap ( filename, ext )
    write ( *, '(2x,a,3x,a,4x,a)' ) string(test), ext_test(test), filename
  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests FILENAME_INC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  character ( len = 20 ) filename
  character ( len = 20 ) filename_old
  integer ( kind = 4 ) j
  character ( len = 20 ) string(test_num)
  integer ( kind = 4 ) test

  string(1) = 'file???.dat'
  string(2) = 'file072.dat'
  string(3) = '2cat9.dat'
  string(4) = 'fred98.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  FILENAME_INC increments a string'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Input             Output'

  do test = 1, test_num
    filename = string(test)
    write ( *, '(a)' ) ' '
    do j = 1, 4
      filename_old = filename
      call filename_inc ( filename )
      write ( *, '(2x,a,2x,a)' ) filename_old, filename
      if ( len_trim ( filename ) <= 0 ) then
        write ( *, '(a)' ) '  (Empty output string.  Quit loop!)'
        exit
      end if
    end do
  end do

  return
end
subroutine test145 ( )

!*****************************************************************************80
!
!! TEST145 tests FILENAME_INC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 20 ) filename
  character ( len = 20 ) filename_old
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST145'
  write ( *, '(a)' ) '  FILENAME_INC increments a string.'
  write ( *, '(a)' ) '  This test checks that a file name is properly'
  write ( *, '(a)' ) '  incremented when carrying is required.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Input             Output'

  filename = 'filename_0000'
  write ( *, '(a)' ) ' '
  do test = 1, 10
    filename_old = filename
    call filename_inc ( filename )
    write ( *, '(2x,a,2x,a)' ) filename_old, filename
  end do

  filename = 'filename_0008'
  write ( *, '(a)' ) ' '
  do test = 1, 10
    filename_old = filename
    call filename_inc ( filename )
    write ( *, '(2x,a,2x,a)' ) filename_old, filename
  end do

  filename = 'filename_0096'
  write ( *, '(a)' ) ' '
  do test = 1, 10
    filename_old = filename
    call filename_inc ( filename )
    write ( *, '(2x,a,2x,a)' ) filename_old, filename
  end do

  filename = 'filename_0995'
  write ( *, '(a)' ) ' '
  do test = 1, 10
    filename_old = filename
    call filename_inc ( filename )
    write ( *, '(2x,a,2x,a)' ) filename_old, filename
  end do

  filename = 'filename_9997'
  write ( *, '(a)' ) ' '
  do test = 1, 10
    filename_old = filename
    call filename_inc ( filename )
    write ( *, '(2x,a,2x,a)' ) filename_old, filename
  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests FILENAME_INC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ios
  character ( len = 20 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  FILENAME_INC "increments" a string, such as'
  write ( *, '(a)' ) '  a file name.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, the string is a file name'
  write ( *, '(a)' ) '  of the form'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     "file_00.txt".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We know we have a sequence of files named'
  write ( *, '(a)' ) '    file_001.txt, file_002.txt, ...'
  write ( *, '(a)' ) '  and we want to generate the name of the next'
  write ( *, '(a)' ) '  file and open it.  If it doesn''t exist, exit.'

  write ( *, '(a)' ) ' '

  s = 'file_00.txt'

  do

    call filename_inc ( s )

    write ( *, '(a)' ) '  Looking for file "' // trim ( s ) // '".'

    open ( unit = 1, file = s, status = 'old', iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  File "' // trim ( s ) // '" does not exist.'
      write ( *, '(a)' ) '  Exiting...'
      exit
    end if

    close ( unit = 1 )

  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests FILENAME_INC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) sim
  character ( len = 20 ) s
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  integer ( kind = 4 ) time_step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  FILENAME_INC "increments" a string, such as'
  write ( *, '(a)' ) '  a file name.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, the string is a file name'
  write ( *, '(a)' ) '  of the form'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     "file_s00_t000.txt'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The user is going to carry out several simulations.'
  write ( *, '(a)' ) '  For each simulation, a number of time steps are done.'
  write ( *, '(a)' ) '  In the file name, the "s" file records the simulation,'
  write ( *, '(a)' ) '  and the "t" field records the time step.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  So a typical file name is "file_s05_t017.txt".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Assuming we have 5 simulations, and 4 time steps,'
  write ( *, '(a)' ) '  the following double loop will generate all the'
  write ( *, '(a)' ) '  file names, from'
  write ( *, '(a)' ) '    "file_s01_t001.txt"'
  write ( *, '(a)' ) '  to'
  write ( *, '(a)' ) '    "file_s05_t004.txt".'

  s1 = 'file_s00'
  s2 = '_t000.txt'

  do sim = 1, 5

    call filename_inc ( s1 )

    s = trim ( s1 ) // trim ( s2 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Begin simulation ', sim
    write ( *, '(a)' ) ' '

    do time_step = 1, 4

      call filename_inc ( s )

      write ( *, '(2x,a)' ) trim ( s )

    end do

  end do

  return
end
subroutine test165 ( )

!*****************************************************************************80
!
!! TEST165 tests FILENAME_INC_NOWRAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  character ( len = 20 ) filename
  character ( len = 20 ) filename_old
  integer ( kind = 4 ) j
  character ( len = 20 ) string(test_num)
  integer ( kind = 4 ) test

  string(1) = 'file???.dat'
  string(2) = 'file072.dat'
  string(3) = '2cat9.dat'
  string(4) = 'fred98.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST165'
  write ( *, '(a)' ) '  FILENAME_INC_NOWRAP increments a string'
  write ( *, '(a)' ) '  but does not allow "wrap around".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Input             Output'

  do test = 1, test_num
    filename = string(test)
    write ( *, '(a)' ) ' '
    do j = 1, 4
      filename_old = filename
      call filename_inc_nowrap ( filename )
      write ( *, '(2x,a,2x,a)' ) filename_old, filename
      if ( len_trim ( filename ) <= 0 ) then
        write ( *, '(a)' ) '  (File name not incrementable.  Quit loop!)'
        exit
      end if
    end do
  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests FILE_PARA_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) filename
  integer ( kind = 4 ) para_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  FILE_PARA_COUNT counts the paragraphs in a file.'

  filename = 'story.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  call file_para_count ( filename, para_num )
  write ( *, '(a,i8)' ) '  Number of paragraphs: ', para_num

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests FILE_PAREN_CHECK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) filename

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  FILE_PAREN_CHECK checks a file for parenthesis errors.'

  filename = 'filum_prb_parens1.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  call file_paren_check ( filename )

  filename = 'filum_prb_parens2.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  call file_paren_check ( filename )

  filename = 'filum_prb_parens3.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  call file_paren_check ( filename )

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests FILE_RENAME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  logical file_exist
  character ( len = 255 ) :: filename_new = 'fred.txt'
  character ( len = 255 ) :: filename_old = 'bob.txt'
  character ( len = 255 ) :: filename_save = 'filum_prb_test.txt'
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) line_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  FILE_RENAME renames a file'
!
!  Make a temporary file.
!
  call file_copy ( filename_save, filename_old, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FILE_COPY failed; cancelling this test.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Old file name: "' // trim ( filename_old ) // '".'
  write ( *, '(a)' ) '  New file name: "' // trim ( filename_new ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  FILE_EXIST(' // trim ( filename_old ) // ') = ', &
    file_exist ( filename_old )
  write ( *, '(a,l1)' ) '  FILE_EXIST(' // trim ( filename_new ) // ') = ', &
    file_exist ( filename_new )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now call FILE_RENAME:'

  call file_rename ( filename_old, filename_new )

  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  FILE_EXIST(' // trim ( filename_old ) // ') = ', &
    file_exist ( filename_old )
  write ( *, '(a,l1)' ) '  FILE_EXIST(' // trim ( filename_new ) // ') = ', &
    file_exist ( filename_new )
!
!  Clean up.
!
  call file_delete ( filename_new )

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests FILE_REVERSE_COLUMNS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) filename_new
  character ( len = 255 ) filename_old

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  FILE_REVERSE_COLUMNS makes a copy of a file with'
  write ( *, '(a)' ) '  each line reversed.'

  filename_old = 'filum_prb_test.txt'
  filename_new = 'filum_prb_reverse_columns.txt'

  call file_reverse_columns ( filename_old, filename_new )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Finished reversal.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input file is  "' // trim ( filename_old ) // '".'
  write ( *, '(a)' ) '  Output file is "' // trim ( filename_new ) // '".'

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests FILE_REVERSE_ROWS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) filename_new
  character ( len = 255 ) filename_old

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  FILE_REVERSE_ROWS makes a copy of a file with the'
  write ( *, '(a)' ) '  lines in reverse order.'

  filename_old = 'filum_prb_test.txt'
  filename_new = 'filum_prb_reverse_rows.txt'

  call file_reverse_rows ( filename_old, filename_new )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Finished reversal.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input file is  "' // trim ( filename_old ) // '".'
  write ( *, '(a)' ) '  Output file is "' // trim ( filename_new ) // '".'

  return
end
subroutine test215 ( )

!*****************************************************************************80
!
!! TEST215 tests FILE_ROT13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) filename1
  character ( len = 255 ) filename2
  character ( len = 255 ) filename3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST215'
  write ( *, '(a)' ) '  FILE_ROT13 makes a copy of a file with'
  write ( *, '(a)' ) '  each line encoded by ROT13.'

  filename1 = 'story.txt'
  filename2 = 'story_rot13.txt'
  filename3 = 'story_rot26.txt'

  call file_rot13 ( filename1, filename2 )
  call file_rot13 ( filename2, filename3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input file is   "' // trim ( filename1 ) // '".'
  write ( *, '(a)' ) '  ROT13 of "' // trim ( filename1 ) // &
    '" is "' // trim ( filename2 ) // '".'
  write ( *, '(a)' ) '  ROT13 of "' // trim ( filename2 ) // &
    '" is "' // trim ( filename3 ) // '".'

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests FILE_ROW_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) :: filename = 'filum_prb_test.txt'
  integer ( kind = 4  ) line_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  FILE_ROW_COUNT counts the "rows" in a file.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  call file_row_count ( filename, line_num )
  write ( *, '(a,i8)' ) '  Number of lines:      ', line_num

  return
end
subroutine test225 ( )

!*****************************************************************************80
!
!! TEST225 tests FILE_SEQUENCE_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4  ) file_dim
  character ( len = 255 ) :: filename = 'data_100.txt'
  integer ( kind = 4  ) file_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST225'
  write ( *, '(a)' ) '  FILE_SEQUENCE_SIZE "sizes" the files in a file sequence.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining files in sequence starting with:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'

  call file_sequence_size ( filename, file_dim, file_num )
  
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of files:      ', file_num
  write ( *, '(a,i8)' ) '  Number of data items: ', file_dim

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests FILE_TAG_CHECK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) filename
  logical file_tag_check
  logical value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  FILE_TAG_CHECK checks a file for tag errors.'

  filename = 'filum_prb_parens1.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  value = file_tag_check ( filename, '(', ')' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  FILE_TAG_CHECK = ', value

  filename = 'filum_prb_parens2.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  value = file_tag_check ( filename, '(', ')' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  FILE_TAG_CHECK = ', value

  filename = 'filum_prb_parens3.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  value = file_tag_check ( filename, '(', ')' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  FILE_TAG_CHECK = ', value

  filename = 'filum_prb_braces.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  value = file_tag_check ( filename, '{', '}' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  FILE_TAG_CHECK = ', value

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests FILE_WORD_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) :: filename = 'filum_prb_test.txt'
  integer ( kind = 4  ) word_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  FILE_WORD_COUNT counts the words in a file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examining file:'
  write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '

  call file_word_count ( filename, word_num )
  write ( *, '(a,i8)' ) '  Number of words:      ', word_num

  return
end
