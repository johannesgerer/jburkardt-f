subroutine file_copy ( old_file_name, new_file_name, ierror )

!*****************************************************************************80
!
!! FILE_COPY makes a copy of a file.
!
!  Discussion:
!
!    The file is assumed to be sequential access, with variable 
!    length records.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OLD_FILE_NAME, the name of the file 
!    to be copied.
!
!    Input, character ( len = * ) NEW_FILE_NAME, the name of the copy of 
!    the file.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, the file names are the same.
!    2, a free unit number could not be found for the old file.
!    3, the routine could not open the old file.
!    4, a free unit number could not be found for the new file.
!    5, the routine could not open the new file.
!
  implicit none

  logical file_exist
  logical file_is_open
  integer ierror
  integer ios
  character ( len = 256 ) line
  character ( len = * ) new_file_name
  integer new_unit
  character ( len = * ) old_file_name
  integer old_unit

  ierror = 0
!
!  Does the original file exist?
!
  if ( .not. file_exist ( old_file_name ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  The old file does not exist.'
    return
  end if
!
!  Is the original file open?
!
  if ( file_is_open ( old_file_name ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  The old file is open.'
    write ( *, '(a)' ) '  It must be closed before it can be copied.'
    return
  end if
!
!  Make sure the file names aren't the same.
!
  if ( new_file_name == old_file_name ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  The old and new file names are identical.'
    return
  end if
!
!  Does the new file exist?
!
  if ( file_exist ( new_file_name ) ) then

    if ( file_is_open ( new_file_name ) ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
      write ( *, '(a)' ) '  A file is already open with the new name.'
      write ( *, '(a)' ) '  It must be closed before it can be overwritten.'
      return
    end if
  
    call file_delete ( new_file_name )

  end if
!
!  At this point:
!
!    The old file exists, and is not open.
!    The new file does not exist, and has a different name.
!
!  Open the old file.
!
  call get_unit ( old_unit )

  if ( old_unit == 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  Could not get a unit number for the old file.'
    return
  end if

  open ( unit = old_unit, file = old_file_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  Could not open the old file:'
    write ( *, '(4x,a)' ) '"' // trim ( old_file_name ) // '".'
    return
  end if
!
!  Open the new file.
!
  call get_unit ( new_unit )

  if ( new_unit == 0 ) then
    ierror = 4
    close ( unit = old_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free unit for the copy file.'
    return
  end if

  open ( unit = new_unit, file = new_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 5
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  Could not open the new file:'
    write ( *, '(4x,a)' ) '"' // trim ( new_file_name ) // '".'
    close ( unit = old_unit )
    return
  end if
!
!  Read an old line, write a new line.
!
  do
 
    read ( old_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    write ( new_unit, '(a)' ) trim ( line )

  end do
 
  close ( unit = old_unit )

  endfile ( unit = new_unit )
  close ( unit = new_unit )

  return
end
subroutine file_delete ( file_name )

!*****************************************************************************80
!
!! FILE_DELETE deletes a named file if it exists.
!
!  Discussion:
!
!    You might want to call this routine to get rid of any old copy
!    of a file, before trying to open a new copy with the OPEN argument:
!      status = 'new'.
!
!    It's not always safe to open a file with " STATUS = 'UNKNOWN' ".
!    For instance, on the SGI, the most recent version of the FORTRAN
!    compiler seems to go crazy when I open an unformatted direct
!    access file this way.  It creates an enormous file (of somewhat
!    random size).  The problem goes away if I delete any old copy
!    using this routine, and then open a fresh copy with
!    " STATUS = 'NEW' ".  It's a scary world.
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
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
  implicit none

  logical file_exist
  logical file_is_open
  character ( len = * ) file_name
  integer ios
  integer iunit
  logical, parameter :: verbose = .false.
!
!  Does the file exist?
!
  if ( .not. file_exist ( file_name ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE - Warning!'
    write ( *, '(a)' ) '  There is no file of the given name.'
    return
  end if
!
!  Is the file open?
!
  if ( file_is_open ( file_name ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE - Warning!'
    write ( *, '(a)' ) '  The file is currently open.'
    write ( *, '(a)' ) '  It must be closed before it can be deleted.'
    return
  end if
!
!  Get a free unit number.
!
  call get_unit ( iunit )

  if ( iunit == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE: Warning!'
    write ( *, '(a)' ) '  A free FORTRAN unit could not be found.'
    return
  end if

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE:'
    write ( *, '(a)' ) '  Deleting "' // trim ( file_name ) // '".'
  end if

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE: Warning!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    return
  end if

  close ( unit = iunit, status = 'delete' )

  return
end
function file_exist ( file_name )

!*****************************************************************************80
!
!! FILE_EXIST reports whether a file exists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, logical FILE_EXIST, is TRUE if the file exists.
!
  implicit none

  character ( len = * ) file_name
  logical file_exist

  inquire ( file = file_name, exist = file_exist )

  return
end
function file_is_open ( file_name )

!*****************************************************************************80
!
!! FILE_IS_OPEN reports whether a file (specified by filename) is open.
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
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, logical FILE_IS_OPEN, is TRUE if the file is open.
!
  implicit none

  character ( len = * ) file_name
  logical file_is_open

  inquire ( file = file_name, opened = file_is_open )

  return
end
subroutine file_rename ( file_name_old, file_name_new )

!*****************************************************************************80
!
!! FILE_RENAME renames a file.
!
!  Discussion:
!
!    Actually, this routine copies the file, and deletes the original.
!    But to the user, it should look like a rename, just a little slower.
!
!    If a file already exists with the new name, it is deleted.
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
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME_OLD, the name of the original file.
!
!    Output, character ( len = * ) FILE_NAME_NEW, the name of the new file.
!
  implicit none

  logical file_exist
  logical file_is_open
  character ( len = * ) file_name_new
  character ( len = * ) file_name_old
  integer ierror
  integer ios
  integer iunit
  logical s_eqi
!
!  Does the old file exist?
!
  if ( .not. file_exist ( file_name_old ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RENAME - Error!'
    write ( *, '(a)' ) '  The original file to be renamed does not exist.'
    return
  end if
!
!  Is the old file open?
!
  if ( file_is_open ( file_name_old ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RENAME - Error!'
    write ( *, '(a)' ) '  The original file is open.'
    write ( *, '(a)' ) '  It must be closed before it can be renamed.'
    return
  end if
!
!  Does old file name = new file name?
!
  if ( s_eqi ( file_name_new, file_name_old ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RENAME: Warning!'
    write ( *, '(a)' ) '  The old and new file names are the same.'
    write ( *, '(a)' ) '  I suppose this means there is nothing to do.'
    return
  end if
!
!  Does the new file exist?
!
  if ( file_exist ( file_name_new ) ) then
!
!  Is the new file open?
!
    if ( file_is_open ( file_name_new ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_RENAME - Error!'
      write ( *, '(a)' ) '  The new file is already open.'
      write ( *, '(a)' ) '  It must be closed before it can be overwritten.'
      return
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_RENAME:'
      write ( *, '(a)' ) '  Deleting pre-existing file with new name.'
      call file_delete ( file_name_new )
    end if

  end if
!
!  Copy old into new.
!
  call file_copy ( file_name_old, file_name_new, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RENAME: Warning!'
    write ( *, '(a)' ) '  Could not copy the old file!'
    return
  end if
!
!  Delete the old file.
!
  call file_delete ( file_name_old )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine r8ad_io ( action, file_name, file_unit, record, n, x )

!*****************************************************************************80
!
!! R8AD_IO reads or writes fixed size vectors, using R8AD protocol.
!
!  Notes:
!
!    Let's cross our fingers and see whether formatted direct access
!    is possible.  That's probably in the book...
!
!    Then consider having a maximum line length of say, 5 values, so that
!    a 13 entry vector takes three lines, and you have enough intelligence
!    to jump around properly...
!
!  Discussion:
!
!    Data type: R8 ( double precision or real ( kind = 8 ) ).
!    Format:    A  ( formatted or ASCII )
!    Access:    D  ( direct )
!
!    It is assumed that the user wants to store and retrieve a number
!    of vectors.  All the vectors are of the same size, and the user
!    always specifies an index or record number when writing or retrieving
!    a particular vector.  At any time, the user can write a new vector,
!    or retrieve any vector that has been written to the file earlier.
!
!    The first call to this routine for a given file should be with
!    the 'C' action, specifying FILE_NAME and the vector size N.  
!    The output of this call is FILE_UNIT, an integer which must be 
!    included on all subsequent calls that manipulate that file.
!
!    To put information into the file, use the 'W' action, specifying
!    the FILE_NAME and FILE_UNIT, as well as the RECORD number,
!    the size N of the vector, and the vector X(1:N) itself.  (Note
!    that every vector written to the file must have the same size.)
!
!    To get information from the file, use the 'R' action, specifying
!    the FILE_NAME and FILE_UNIT, as well as the RECORD number,
!    the size N of the vector.  The routine will return the given vector
!    in X(1:N).  
!
!    When you are done with the file, you can close and save it with the
!    'F' action, if you specify the FILE_NAME and FILE_UNIT.
!
!    You can, instead, close and delete the file, using the 'D'
!    action and FILE_NAME and FILE_UNIT.
!
!
!    Things that can go wrong: 
!
!    * You should create a file before doing anything else to it.
!
!    * You should only read a record if you have already written it.
!
!    * You should always use exactly the same set of values FILE_NAME,
!      FILE_UNIT and N, when accessing the file.
!
!    * Because I want to allow multiple files, I haven't kept track
!      of properties associated with a single file.  In particular,
!      I don't keep track of the maximum record written, which would
!      make it easy for the user to write the next record, or to
!      determine the total size of the file.
!
!    Features:
!
!    * You can read or write or rewrite any record at any time.
!
!    * You can handle two or more files simultaneously, as long as
!      you use the appropriate set of FILE_NAME, FILE_UNIT and N
!      values for each.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ACTION:
!    'c' to create the file FILE_NAME with record size N;
!    'w' to write vector X(1:N) as record RECORD;
!    'r' to read vector X(1:N) which is record RECORD;
!    'f' to close and save the file FILE_NAME;
!    'd' to close and delete the file FILE_NAME.
!    's' to print some statistics about file FILE_NAME.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Input/output, integer FILE_UNIT, the unit number associated with the file.
!    For the 'c' command, this value is output.
!    For the 'w', 'r', 'f' and 'd' commands, this value is input.
!    
!    Input, integer RECORD, the index of the record to be written or read.
!    (This value is needed for the 'W' and 'R' actions only.)
!
!    Input, integer N, the size of each vector.
!    (This value is needed for the 'C', 'W' and 'R' actions only.)
!
!    Input/output, real ( kind = 8 ) X(N), the vector to be written or read.
!    (This value is input for the 'W' action, output for the 'R' action,
!    and not needed otherwise.)
!
  implicit none

  integer n

  character action
  integer, save :: call_num = 0
  integer, parameter :: char_length = 16
  logical file_exist
  logical file_is_open
  character ( len = * ) file_name
  integer file_unit
  integer ios
  integer record
  integer, save :: record_length = -1
  real ( kind = 8 ) x(n)
  character ( len = 22 ), save :: x_format = ' '
!
!  Create a file of given name.
!
  if ( action(1:1) == 'c' ) then

    if ( file_exist ( file_name ) ) then
      call file_delete ( file_name )
    end if

    if ( n <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8AD_IO - Fatal error!'
      write ( *, '(a)' ) '  The value of N is not positive.'
      write ( *, '(a)' ) '  Could not create "' // trim ( file_name ) // '".'
      stop
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Creating file "' // trim ( file_name ) // '".'

    call get_unit ( file_unit )

    record_length = n * char_length

    write ( x_format, '(a,i8,a)' ) '(', n, 'g16.6)'

    open ( file = file_name, unit = file_unit, status = 'new', &
      form = 'formatted', access = 'direct', recl = record_length, &
      iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8AD_IO - Fatal error!'
      write ( *, '(a,i8)' ) '  IO error of type ', ios
      write ( *, '(a)' ) '  Could not create "' // trim ( file_name ) // '".'
      stop
    end if

    call_num = 1
!
!  Close (if necessary) and delete the file.
!
  else if ( action(1:1) == 'd' ) then

    call_num = call_num + 1

    if ( .not. file_exist ( file_name ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8AD_IO - Warning!'
      write ( *, '(a)' ) '  The file does not exist.'
      write ( *, '(a)' ) '  Could not delete "' // trim ( file_name ) // '".'
      return
    end if

    if ( file_is_open ( file_name ) ) then
      close ( unit = file_unit )
    end if

    call file_delete ( file_name )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The file "' // trim ( file_name ) // &
      '" has been deleted.'

    file_unit = -1
    record_length = -1
!
!  Close (if necessary) and save the file.
!
  else if ( action(1:1) == 'f' ) then

    call_num = call_num + 1

    if ( .not. file_exist ( file_name ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8AD_IO - Warning!'
      write ( *, '(a)' ) '  The file does not exist.'
      write ( *, '(a)' ) '  Could not close and save "' &
        // trim ( file_name ) // '".'
      return
    end if

    if ( file_is_open ( file_name ) ) then
      close ( unit = file_unit )
    end if
!
!  Read a vector X.
!
  else if ( action(1:1) == 'r' ) then

    call_num = call_num + 1

    if ( .not. file_exist ( file_name ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8AD_IO - Fatal error!'
      write ( *, '(a)' ) '  You have asked to READ data from a file'
      write ( *, '(a)' ) '  that does not exist.'
      write ( *, '(a)' ) '  The file name is "' // trim ( file_name ) // '".'
      stop
    end if

    if ( record < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8AD_IO - Fatal error!'
      write ( *, '(a)' ) '  You have asked to read an illegal record.'
      write ( *, '(a)' ) '  Could not read from "' // trim ( file_name ) // '".'
      stop
    end if
    
    read ( unit = file_unit, fmt = x_format, rec = record,  &
      iostat = ios ) x(1:n)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8AD_IO - Fatal error!'
      write ( *, '(a,i8)' ) '  IO error of type ', ios
      write ( *, '(a)' ) '  Could not read from "' // trim ( file_name ) // '".'
      stop
    end if
!
!  Print statistics.
!
  else if ( action(1:1) == 's' ) then

    call_num = call_num + 1

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8AD_IO - File statistics:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The file name is "' // trim ( file_name ) // '".'
    write ( *, '(a,i12)' ) '  The length of each vector is N =    ', n
    write ( *, '(a,i12)' ) '  The character length of one data item =  ', &
      char_length
    write ( *, '(a,i12)' ) '  The length of each record is =      ', &
      record_length
    write ( *, '(a)' ) '  The record format is "' &
      // trim ( x_format ) // '".'
!
!  Write a vector X.
!
  else if ( action(1:1) == 'w' ) then

    call_num = call_num + 1

    if ( .not. file_exist ( file_name ) ) then

      if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8AD_IO - Fatal error!'
        write ( *, '(a)' ) '  The value of N is not positive.'
        write ( *, '(a)' ) '  Could not write to "' &
          // trim ( file_name ) // '".'
        stop
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Creating file "' // trim ( file_name ) // '".'

      call get_unit ( file_unit )

      record_length = n * char_length

      write ( x_format, '(a,i8,a)' ) '(', n, 'g16.6)'

      open ( file = file_name, unit = file_unit, status = 'new', &
        form = 'formatted', access = 'direct', recl = record_length, &
        iostat = ios )

    end if

    if ( record < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8AD_IO - Fatal error!'
      write ( *, '(a)' ) '  You have asked to write an illegal record.'
      write ( *, '(a)' ) '  Could not write to "' // trim ( file_name ) // '".'
      stop
    end if
    
    write ( unit = file_unit, fmt = x_format, rec = record, &
      iostat = ios ) x(1:n)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8AD_IO - Fatal error!'
      write ( *, '(a,i8)' ) '  IO error of type ', ios
      write ( *, '(a)' ) '  Could not write to "' // trim ( file_name ) // '".'
      stop
    end if
!
!  Unrecognized ACTION.
!
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8AD_IO - Warning!'
    write ( *, '(a)' ) '  The requested action "' // trim ( action ) // &
      '" is not recognized!'
    write ( *, '(a)' ) '  The file was "' // trim ( file_name ) // '".'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Legal actions include:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "c", create a file;'
    write ( *, '(a)' ) '  "d", delete a file;'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "w", write a vector;'
    write ( *, '(a)' ) '  "r", read a vector;'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "s", print statistics;'
    
  end if

  return
end
subroutine r8ud_io ( action, file_name, file_unit, record, n, x )

!*****************************************************************************80
!
!! R8UD_IO reads or writes fixed size vectors, using R8UD protocol.
!
!  Discussion:
!
!    Data type: R8 ( double precision or real ( kind = 8 ) ).
!    Format:    U  ( unformatted or binary )
!    Access:    D  ( direct access )
!
!    It is assumed that the user wants to store and retrieve a number
!    of vectors.  All the vectors are of the same size, and the user
!    always specifies an index or record number when writing or retrieving
!    a particular vector.  At any time, the user can write a new vector,
!    or retrieve any vector that has been written to the file earlier.
!
!    The first call to this routine for a given file should be with
!    the 'C' action, specifying FILE_NAME and the vector size N.  
!    The output of this call is FILE_UNIT, an integer which must be 
!    included on all subsequent calls that manipulate that file.
!
!    To put information into the file, use the 'W' action, specifying
!    the FILE_NAME and FILE_UNIT, as well as the RECORD number,
!    the size N of the vector, and the vector X(1:N) itself.  (Note
!    that every vector written to the file must have the same size.)
!
!    To get information from the file, use the 'R' action, specifying
!    the FILE_NAME and FILE_UNIT, as well as the RECORD number,
!    the size N of the vector.  The routine will return the given vector
!    in X(1:N).  
!
!    When you are done with the file, you can close and save it with the
!    'F' action, if you specify the FILE_NAME and FILE_UNIT.
!
!    You can, instead, close and delete the file, using the 'D'
!    action and FILE_NAME and FILE_UNIT.
!
!
!    Things that can go wrong: 
!
!    * You should create a file before doing anything else to it.
!
!    * You should only read a record if you have already written it.
!
!    * You should always use exactly the same set of values FILE_NAME,
!      FILE_UNIT and N, when accessing the file.
!
!    * Because I want to allow multiple files, I haven't kept track
!      of properties associated with a single file.  In particular,
!      I don't keep track of the maximum record written, which would
!      make it easy for the user to write the next record, or to
!      determine the total size of the file.
!
!    Features:
!
!    * You can read or write or rewrite any record at any time.
!
!    * You can handle two or more files simultaneously, as long as
!      you use the appropriate set of FILE_NAME, FILE_UNIT and N
!      values for each.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ACTION:
!    'c' to create the file FILE_NAME with record size N;
!    'w' to write vector X(1:N) as record RECORD;
!    'r' to read vector X(1:N) which is record RECORD;
!    'f' to close and save the file FILE_NAME;
!    'd' to close and delete the file FILE_NAME.
!    's' to print some statistics about file FILE_NAME.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Input/output, integer FILE_UNIT, the unit number associated with the file.
!    For the 'c' command, this value is output.
!    For the 'w', 'r', 'f' and 'd' commands, this value is input.
!    
!    Input, integer RECORD, the index of the record to be written or read.
!    (This value is needed for the 'W' and 'R' actions only.)
!
!    Input, integer N, the size of each vector.
!    (This value is needed for the 'C', 'W' and 'R' actions only.)
!
!    Input/output, real ( kind = 8 ) X(N), the vector to be written or read.
!    (This value is input for the 'W' action, output for the 'R' action,
!    and not needed otherwise.)
!
  implicit none

  integer n
  integer, parameter :: word_length = 2

  character action
  integer, parameter :: byte_length = 4 * word_length
  integer, save :: call_num = 0
  logical file_exist
  logical file_is_open
  character ( len = * ) file_name
  integer file_unit
  integer ios
  integer record
  integer, save :: record_length = -1
  real ( kind = 8 ) x(n)
!
!  Create a file of given name.
!
  if ( action(1:1) == 'c' ) then

    if ( file_exist ( file_name ) ) then
      call file_delete ( file_name )
    end if

    if ( n <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a)' ) '  The value of N is not positive.'
      write ( *, '(a)' ) '  Could not create "' // trim ( file_name ) // '".'
      stop
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Creating file "' // trim ( file_name ) // '".'

    call get_unit ( file_unit )
!
!  NOTE:
!    On some machines, the record length of a binary direct access file
!    is in WORDS, and on some machines, it is in BYTES.
!    So on some systems (Alpha's for one), use the first statement.
!    on others (Apple's, for one), use the second.
!
    record_length = n * word_length
!   record_length = n * byte_length

    open ( file = file_name, unit = file_unit, status = 'new', &
      form = 'unformatted', access = 'direct', recl = record_length, &
      iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a,i8)' ) '  IO error of type ', ios
      write ( *, '(a)' ) '  Could not create "' // trim ( file_name ) // '".'
      stop
    end if

    call_num = 1
!
!  Close (if necessary) and delete the file.
!
  else if ( action(1:1) == 'd' ) then

    call_num = call_num + 1

    if ( .not. file_exist ( file_name ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The file does not exist.'
      write ( *, '(a)' ) '  Could not delete "' // trim ( file_name ) // '".'
      return
    end if

    if ( file_is_open ( file_name ) ) then
      close ( unit = file_unit )
    end if

    call file_delete ( file_name )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The file "' // trim ( file_name ) // &
      '" has been deleted.'

    file_unit = -1
    record_length = -1
!
!  Close (if necessary) and save the file.
!
  else if ( action(1:1) == 'f' ) then

    call_num = call_num + 1

    if ( .not. file_exist ( file_name ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Warning!'
      write ( *, '(a)' ) '  The file does not exist.'
      write ( *, '(a)' ) '  Could not close and save "' &
        // trim ( file_name ) // '".'
      return
    end if

    if ( file_is_open ( file_name ) ) then
      close ( unit = file_unit )
    end if
!
!  Read a vector X.
!
  else if ( action(1:1) == 'r' ) then

    call_num = call_num + 1

    if ( .not. file_exist ( file_name ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a)' ) '  You have asked to READ data from a file'
      write ( *, '(a)' ) '  that does not exist.'
      write ( *, '(a)' ) '  The file name is "' // trim ( file_name ) // '".'
      stop
    end if

    if ( record < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a)' ) '  You have asked to read an illegal record.'
      write ( *, '(a)' ) '  Could not read from "' // trim ( file_name ) // '".'
      stop
    end if
    
    read ( unit = file_unit, rec = record, iostat = ios ) x(1:n)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a,i8)' ) '  IO error of type ', ios
      write ( *, '(a)' ) '  Could not read from "' // trim ( file_name ) // '".'
      stop
    end if
!
!  Print statistics.
!
  else if ( action(1:1) == 's' ) then

    call_num = call_num + 1

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8UD_IO - File statistics:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The file name is "' // trim ( file_name ) // '".'
    write ( *, '(a,i12)' ) '  The length of each vector is N =    ', n
    write ( *, '(a,i12)' ) '  The word length of one data item =  ', word_length
    write ( *, '(a,i12)' ) '  The byte length of one data item =  ', byte_length
    write ( *, '(a,i12)' ) '  The length of each record is =      ', &
      record_length
!
!  Write a vector X.
!
  else if ( action(1:1) == 'w' ) then

    call_num = call_num + 1

    if ( .not. file_exist ( file_name ) ) then

      if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
        write ( *, '(a)' ) '  The value of N is not positive.'
        write ( *, '(a)' ) '  Could not write to "' &
          // trim ( file_name ) // '".'
        stop
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Creating file "' // trim ( file_name ) // '".'

      call get_unit ( file_unit )

      record_length = n * byte_length

      open ( file = file_name, unit = file_unit, status = 'new', &
        form = 'unformatted', access = 'direct', recl = record_length, &
        iostat = ios )

    end if

    if ( record < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a)' ) '  You have asked to write an illegal record.'
      write ( *, '(a)' ) '  Could not write to "' // trim ( file_name ) // '".'
      stop
    end if
    
    write ( unit = file_unit, rec = record, iostat = ios ) x(1:n)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a,i8)' ) '  IO error of type ', ios
      write ( *, '(a)' ) '  Could not write to "' // trim ( file_name ) // '".'
      stop
    end if
!
!  Unrecognized ACTION.
!
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8UD_IO - Warning!'
    write ( *, '(a)' ) '  The requested action "' // trim ( action ) // &
      '" is not recognized!'
    write ( *, '(a)' ) '  The file was "' // trim ( file_name ) // '".'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Legal actions include:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "c", create a file;'
    write ( *, '(a)' ) '  "d", delete a file;'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "w", write a vector;'
    write ( *, '(a)' ) '  "r", read a vector;'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "s", print statistics;'
    
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
