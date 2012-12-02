subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine i4_to_s_zero ( intval, s )

!*****************************************************************************80
!
!! I4_TO_S_ZERO converts an integer to a string, with zero padding.
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  000001
!        -1  -00001
!         0  000000
!      1952  001952
!    123456  123456
!   1234567  ******  <-- Not enough room!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be right justified, and zero padded.
!    If there is not enough space, the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * ) s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  Working from right to left, strip off the digits of the integer
!  and place them into S(ILO:IHI).
!
  ipos = ihi

  do while ( ival /= 0 .or. ipos == ihi )

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

  end do
!
!  Fill the empties with zeroes.
!
  do i = ilo, ipos
    s(i:i) = '0'
  end do

  return
end
subroutine pfile_close ( npx, npy, base_unit )

!*****************************************************************************80
!
!! PFILE_CLOSE closes all the "parallel" files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPX, NPY, the dimension of the processor
!    grid in the X and Y directions.
!
!    Input, integer ( kind = 4 ) BASE_UNIT, the FORTRAN unit number to be assigned the
!    file associated with processor (1,1).  Subsequent files are associated
!    with subsequent unit numbers.
!    BASE_UNIT should be between 0 and 99.  Other restriction may apply.
!    An error is likely to occur if BASE_UNIT + NPX * NPY - 1 > 99.
!
  implicit none

  integer ( kind = 4 ) npx
  integer ( kind = 4 ) npy

  integer ( kind = 4 ) base_unit
  integer ( kind = 4 ) i_proc
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j_proc

  iunit = base_unit

  do j_proc = 1, npy

    do i_proc = 1, npx

      close ( unit = iunit, iostat = ios )

      iunit = iunit + 1

    end do

  end do

  return
end
subroutine pfile_name_gen ( file_name_fore, file_name_mid, file_name_post, &
  proc_no, file_name )

!*****************************************************************************80
!
!! PFILE_NAME_GEN generates a file name that contains a processor ID.
!
!  Discussion:
!
!    The filename is assumed to be made up of a beginning, middle, and
!    end.  The middle part is supplied by the user simply to specify the
!    size of the field to contain the processor number.  On output, this
!    field contains the zero-filled processor number.  Either of the
!    beginning and end parts may be blank.
!
!  Example:
!
!    Input:
!
!      FILE_NAME_FORE = 'proc'
!      FILE_NAME_MID  = '****'
!      FILE_NAME_POST = '.dat'
!      PROC_NO = 17
!
!    Output:
!
!      FILE_NAME = 'proc0017.dat'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME_FORE, FILE_NAME_MID, FILE_NAME_POST,
!    the beginning, middle, and end of the string used to construct the
!    family of file names.
!
!    Input, integer ( kind = 4 ) PROC_NO, the processor number associated
!    with the file whose name is to be constructed.
!
!    Output, character ( len = * ) FILE_NAME, the file name associated with the
!    given processor number.
!
  implicit none

  character ( len = * ) file_name
  character ( len = * ) file_name_fore
  character ( len = * ) file_name_mid
  character ( len = * ) file_name_post
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lenc1
  integer ( kind = 4 ) lenc2
  integer ( kind = 4 ) lenc3
  integer ( kind = 4 ) proc_no
  character ( len = 4 ) string4

  lenc1 = len_trim ( file_name_fore )
  lenc2 = len_trim ( file_name_mid )
  lenc3 = len_trim ( file_name_post )
!
!  Write the processor ID into the middle string as a zero-filled integer.
!
  call i4_to_s_zero ( proc_no, string4 )

  if ( 4 <= lenc2 ) then
    file_name_mid = string4
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_NAME_GEN - Fatal error!'
    write ( *, '(a)' ) '  Too much data to fit into file name.'
    stop
  end if
!
!  Make sure that FILE_NAME is big enough.
!
  lenc = len ( file_name )

  if ( lenc < lenc1 + lenc2 + lenc3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_NAME_GEN - Fatal error!'
    write ( *, '(a)' ) '  Too much data to fit into file name.'
    stop
  end if
!
!  Concatenate the non-blank information in the strings.
!
  lenc = 0
  file_name = ' '

  if ( 0 < lenc1 ) then
    file_name(1:lenc1) = file_name_fore(1:lenc1)
    lenc = lenc1
  end if

  if ( 0 < lenc2 ) then
    file_name(lenc+1:lenc+lenc2) = file_name_mid(1:lenc2)
    lenc = lenc + lenc2
  end if

  if ( 0 < lenc3 ) then
    file_name(lenc+1:lenc+lenc3) = file_name_post(1:lenc3)
    lenc = lenc + lenc3
  end if

  return
end
subroutine pfile_open_read ( npx, npy, file_name_fore, file_name_mid, &
  file_name_post, base_unit, ierror, pform )

!*****************************************************************************80
!
!! PFILE_OPEN_READ opens all the "parallel" files for reading.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPX, NPY, the dimension of the processor
!    grid in the X and Y directions.
!
!    Input, character ( len = * ) FILE_NAME_FORE, FILE_NAME_MID, FILE_NAME_POST,
!    the beginning, middle, and end of the string used to construct the
!    family of file names.
!
!    Input, integer ( kind = 4 ) BASE_UNIT, the FORTRAN unit number to be assigned the
!    file associated with processor (1,1).  Subsequent files are associated
!    with subsequent unit numbers.
!    BASE_UNIT should be between 0 and 99.  Other restriction may apply.
!    An error is likely to occur if BASE_UNIT + NPX * NPY - 1 > 99.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error;
!    1, could not open the file.
!
!    Input, character ( len = * ) PFORM, specifies the form of the parallel
!    files.
!    'formatted', the files are formatted;
!    'unformatted' or any other value, the files are unformatted.
!
  implicit none

  integer ( kind = 4 ) npx
  integer ( kind = 4 ) npy

  integer ( kind = 4 ) base_unit
  character ( len = 100 ) file_name
  character ( len = * ) file_name_fore
  character ( len = * ) file_name_mid
  character ( len = * ) file_name_post
  character ( len = * ) pform
  integer ( kind = 4 ) i_proc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j_proc
  integer ( kind = 4 ) proc_no

  ierror = 0

  iunit = base_unit

  proc_no = 0

  do j_proc = 1, npy

    do i_proc = 1, npx

      call pfile_name_gen ( file_name_fore, file_name_mid, &
        file_name_post, proc_no, file_name )

      open ( unit = iunit, file = file_name, form = pform, status = 'old', &
        iostat = ios )

      if ( ios /= 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PFILE_OPEN_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open one of the parallel files:'
        write ( *, '(a)' ) trim ( file_name )
        return
      end if

      iunit = iunit + 1
      proc_no = proc_no + 1

    end do

  end do

  return
end
subroutine pfile_open_write ( npx, npy, file_name_fore, file_name_mid, &
  file_name_post, base_unit, ierror, pform )

!*****************************************************************************80
!
!! PFILE_OPEN_WRITE opens all the "parallel" files for writing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPX, NPY, the dimension of the processor
!    grid in the X and Y directions.
!
!    Input, character ( len = * ) FILE_NAME_FORE, FILE_NAME_MID, FILE_NAME_POST,
!    the beginning, middle, and end of the string used to construct the
!    family of file names.
!
!    Input, integer ( kind = 4 ) BASE_UNIT, the FORTRAN unit number to be assigned the
!    file associated with processor (1,1).  Subsequent files are associated
!    with subsequent unit numbers.
!    BASE_UNIT should be between 0 and 99.  Other restriction may apply.
!    An error is likely to occur if BASE_UNIT + NPX * NPY - 1 > 99.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error;
!    1, could not open the file.
!
!    Input, character ( len = * ) PFORM, specifies the form of the parallel
!    files.
!    'formatted', the files are formatted;
!    'unformatted' or any other value, the files are unformatted.
!
  implicit none

  integer ( kind = 4 ) npx
  integer ( kind = 4 ) npy

  integer ( kind = 4 ) base_unit
  character ( len = 100 ) file_name
  character ( len = * ) file_name_fore
  character ( len = * ) file_name_mid
  character ( len = * ) file_name_post
  character ( len = * ) pform
  integer ( kind = 4 ) i_proc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j_proc
  integer ( kind = 4 ) proc_no

  ierror = 0

  iunit = base_unit

  proc_no = 0

  do j_proc = 1, npy

    do i_proc = 1, npx

      call pfile_name_gen ( file_name_fore, file_name_mid, &
        file_name_post, proc_no, file_name )

      open ( unit = iunit, file = file_name, form = pform, status = 'replace', &
        iostat = ios )

      if ( ios /= 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PFILE_OPEN_WRITE - Fatal error!'
        write ( *, '(a)' ) '  Could not open one of the parallel files:'
        write ( *, '(a)' ) trim ( file_name )
        return
      end if

      iunit = iunit + 1
      proc_no = proc_no + 1

    end do

  end do

  return
end
subroutine pfile_read_pdata ( npx, npy, nx_global, ny_global, nx_local, &
  ny_local, nz, base_unit, a_global, a_local, ierror, pform )

!*****************************************************************************80
!
!! PFILE_READ_PDATA reads "parallel" data from the parallel files.
!
!  Discussion:
!
!    Construct a vector A_GLOBAL(I,J,K) out of data vectors A_LOCAL(I2,J2,K2)
!    which are equal sized "parcels" of the data, divided up among
!    a grid of NPX by NPY processors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPX, NPY, the dimension of the processor
!    grid in the X and Y directions.
!
!    Input, integer ( kind = 4 ) NX_GLOBAL, NY_GLOBAL, the number of grid points in the
!    X and Y directions, for the global grid.
!
!    Input, integer ( kind = 4 ) NX_LOCAL, NY_LOCAL, the number of grid points in the
!    X and Y directions, for the local grid associated with a single
!    processor.
!
!    Input, integer ( kind = 4 ) NZ, the third dimension of the array A to be read.
!    In some cases, this is 1, because the array represents scalar data.
!
!    Input, integer ( kind = 4 ) BASE_UNIT, the FORTRAN unit number to be assigned the
!    file associated with processor (1,1).  Subsequent files are associated
!    with subsequent unit numbers.
!    BASE_UNIT should be between 0 and 99.  Other restriction may apply.
!    An error is likely to occur if BASE_UNIT + NPX * NPY - 1 > 99.
!
!    Workspace, real A_GLOBAL(NX_GLOBAL,NY_GLOBAL,NZ).
!
!    Workspace, real A_LOCAL(NX_LOCAL,NY_LOCAL,NZ).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!
!    Input, character ( len = * ) PFORM, specifies the form of the parallel
!    files.
!    'formatted', the files are formatted;
!    'unformatted' or any other value, the files are unformatted.
!
  implicit none

  integer ( kind = 4 ) npx
  integer ( kind = 4 ) npy
  integer ( kind = 4 ) nx_global
  integer ( kind = 4 ) nx_local
  integer ( kind = 4 ) ny_global
  integer ( kind = 4 ) ny_local
  integer ( kind = 4 ) nz

  real a_global(nx_global,ny_global,nz)
  real a_local(nx_local,ny_local,nz)
  integer ( kind = 4 ) base_unit
  integer ( kind = 4 ) i_global
  integer ( kind = 4 ) i_local
  integer ( kind = 4 ) i_proc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j_global
  integer ( kind = 4 ) j_local
  integer ( kind = 4 ) j_proc
  integer ( kind = 4 ) k
  character ( len = * ) pform

  ierror = 0

  iunit = base_unit

  do j_proc = 1, npy

    do i_proc = 1, npx

      if ( pform == 'formatted' ) then
        do k = 1, nz
          do j_local = 1, ny_local
            do i_local = 1, nx_local
              read ( iunit, *, iostat = ios ) a_local(i_local,j_local,k)
              if ( ios /= 0 ) then
                ierror = 2
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'PFILE_READ_PDATA - Fatal error!'
                write ( *, '(a)' ) '  Error while reading from a file.'
                return
              end if
            end do
          end do
        end do
      else
        read ( iunit, iostat = ios ) a_local
        if ( ios /= 0 ) then
          ierror = 2
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PFILE_READ_PDATA - Fatal error!'
          write ( *, '(a)' ) '  Error while reading from a file.'
          return
        end if
      end if

      do k = 1, nz
        do j_local = 1, ny_local
          j_global = ( j_proc - 1 ) * ny_local + j_local
          do i_local = 1, nx_local
            i_global = ( i_proc - 1 ) * nx_local + i_local
            a_global(i_global,j_global,k) = a_local(i_local,j_local,k)
          end do
        end do
      end do

      iunit = iunit + 1

    end do

  end do

  return
end
subroutine pfile_read_sdata ( npx, npy, base_unit, a, nval, ierror, pform )

!*****************************************************************************80
!
!! PFILE_READ_SDATA reads "scalar data" from a set of parallel files.
!
!  Discussion:
!
!    It is assumed that the next record of every parallel file contains
!    the same data, namely, NVAL real values.
!
!    The data is read from EVERY file, so that all of them are advanced
!    to the next record.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPX, NPY, the dimension of the processor
!    grid in the X and Y directions.
!
!    Input, integer ( kind = 4 ) BASE_UNIT, the FORTRAN unit number to be assigned the
!    file associated with processor (1,1).  Subsequent files are associated
!    with subsequent unit numbers.
!    BASE_UNIT should be between 0 and 99.  Other restriction may apply.
!    An error is likely to occur if BASE_UNIT + NPX * NPY - 1 > 99.
!
!    Output, real A(NVAL), the data, as read from one of the files.
!
!    Input, integer ( kind = 4 ) NVAL, the number of entries of A.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error;
!    1, end of file;
!    2, format error during a read.
!
!    Input, character ( len = * ) PFORM, specifies the form of the parallel
!    files.
!    'formatted', the files are formatted;
!    'unformatted' or any other value, the files are unformatted.
!
  implicit none

  integer ( kind = 4 ) npx
  integer ( kind = 4 ) npy
  integer ( kind = 4 ) nval

  real a(nval)
  integer ( kind = 4 ) base_unit
  character ( len = * ) pform
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_proc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j_proc

  ierror = 0

  iunit = base_unit

  do j_proc = 1, npy

    do i_proc = 1, npx

      if ( pform == 'formatted' ) then
        read ( iunit, *, iostat = ios ) a(1:nval)
      else
        read ( iunit, iostat = ios ) a
      end if

      if ( ios /= 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PFILE_READ_SDATA - Fatal error!'
        write ( *, '(a)' ) '  Error while reading from a file.'
        return
      end if

      iunit = iunit + 1

    end do

  end do

  return
end
subroutine pfile_write_pdata ( npx, npy, nx_global, ny_global, nx_local, &
  ny_local, nz, base_unit, a_global, a_local, ierror, pform )

!*****************************************************************************80
!
!! PFILE_WRITE_PDATA writes "parallel" data to the parallel files.
!
!  Discussion:
!
!    Parcel out a vector A_GLOBAL(I,J,K) into data vectors A_LOCAL(I2,J2,K2)
!    which are equal sized "parcels" of the data, divided up among
!    a grid of NPX by NPY processors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPX, NPY, the dimension of the processor
!    grid in the X and Y directions.
!
!    Input, integer ( kind = 4 ) NX_GLOBAL, NY_GLOBAL, the number of grid
!    points in the X and Y directions, for the global grid.
!
!    Input, integer ( kind = 4 ) NX_LOCAL, NY_LOCAL, the number of grid
!    points in the X and Y directions, for the local grid associated with
!    a single processor.
!
!    Input, integer ( kind = 4 ) NZ, the third dimension of the array A to be read.
!    In some cases, this is 1, because the array represents scalar data.
!
!    Input, integer ( kind = 4 ) BASE_UNIT, the FORTRAN unit number to be assigned the
!    file associated with processor (1,1).  Subsequent files are associated
!    with subsequent unit numbers.
!    BASE_UNIT should be between 0 and 99.  Other restriction may apply.
!    An error is likely to occur if BASE_UNIT + NPX * NPY - 1 > 99.
!
!    Workspace, real A_GLOBAL(NX_GLOBAL,NY_GLOBAL,NZ).
!
!    Workspace, real A_LOCAL(NX_LOCAL,NY_LOCAL,NZ).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!
!    Input, character ( len = * ) PFORM, specifies the form of the parallel
!    files.
!    'formatted', the files are formatted;
!    'unformatted' or any other value, the files are unformatted.
!
  implicit none

  integer ( kind = 4 ) npx
  integer ( kind = 4 ) npy
  integer ( kind = 4 ) nx_global
  integer ( kind = 4 ) nx_local
  integer ( kind = 4 ) ny_global
  integer ( kind = 4 ) ny_local
  integer ( kind = 4 ) nz

  real a_global(nx_global,ny_global,nz)
  real a_local(nx_local,ny_local,nz)
  integer ( kind = 4 ) base_unit
  integer ( kind = 4 ) i_global
  integer ( kind = 4 ) i_local
  integer ( kind = 4 ) i_proc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j_global
  integer ( kind = 4 ) j_local
  integer ( kind = 4 ) j_proc
  integer ( kind = 4 ) k
  character ( len = * ) pform

  ierror = 0

  iunit = base_unit

  do j_proc = 1, npy

    do i_proc = 1, npx

      do k = 1, nz
        do j_local = 1, ny_local
          j_global = ( j_proc - 1 ) * ny_local + j_local
          do i_local = 1, nx_local
            i_global = ( i_proc - 1 ) * nx_local + i_local
            a_local(i_local,j_local,k) = a_global(i_global,j_global,k)
          end do
        end do
      end do

      if ( pform == 'formatted' ) then
        do k = 1, nz
          do j_local = 1, ny_local
            do i_local = 1, nx_local
              write ( iunit, * ) a_local(i_local,j_local,k)
            end do
          end do
        end do
      else
        write ( iunit ) a_local
      end if

      iunit = iunit + 1

    end do

  end do

  return
end
subroutine pfile_write_sdata ( npx, npy, base_unit, a, nval, ierror, pform )

!*****************************************************************************80
!
!! PFILE_WRITE_SDATA writes "scalar data" to a set of parallel files.
!
!  Discussion:
!
!    It is assumed that the next record of every parallel file is to
!    contain the same data, namely, NVAL real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPX, NPY, the dimension of the processor
!    grid in the X and Y directions.
!
!    Input, integer ( kind = 4 ) BASE_UNIT, the FORTRAN unit number to be
!    assigned the file associated with processor (1,1).  Subsequent files
!    are associated with subsequent unit numbers.
!    BASE_UNIT should be between 0 and 99.  Other restriction may apply.
!    An error is likely to occur if BASE_UNIT + NPX * NPY - 1 > 99.
!
!    Input, real A(NVAL), the data, to be written to each of the files.
!
!    Input, integer ( kind = 4 ) NVAL, the number of entries of A.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!
!    Input, character ( len = * ) PFORM, specifies the form of the parallel
!    files.
!    'formatted', the files are formatted;
!    'unformatted' or any other value, the files are unformatted.
!
  implicit none

  integer ( kind = 4 ) npx
  integer ( kind = 4 ) npy
  integer ( kind = 4 ) nval

  real a(nval)
  integer ( kind = 4 ) base_unit
  character ( len = * ) pform
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_proc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j_proc

  ierror = 0

  iunit = base_unit

  do j_proc = 1, npy

    do i_proc = 1, npx

      if ( pform == 'formatted' ) then
        write ( iunit, * ) ( a(i), i = 1, nval )
      else
        write ( iunit ) a
      end if

      iunit = iunit + 1

    end do

  end do

  return
end
subroutine rejoin_save ( i_status, npx, npy, nsc, nx_global, ny_global, &
  nx_local, ny_local, time, pform, sform )

!*****************************************************************************80
!
!! REJOIN_SAVE "rejoins" parallel DNS SAVE files into one sequential file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_STATUS is:
!    0, to convert the FIELD_T_... files;
!    1, to convert the FIELD_L_... files.
!
!    Input, integer ( kind = 4 ) NPX, NPY, the dimension of the processor
!    grid in the X and Y directions.
!
!    Input, integer ( kind = 4 ) NSC, used to define the third dimension
!    of QZ1_LOCAL.
!
!    Input, integer ( kind = 4 ) NX_GLOBAL, NY_GLOBAL, the number of grid
!    points in the X and Y directions, for the global grid.
!
!    Input, integer ( kind = 4 ) NX_LOCAL, NY_LOCAL, the number of grid
!    points in the X and Y directions, for the local grid associated with
!    a single processor.
!
!    Input, real TIME, the time, which is only needed if I_STATUS = 0,
!    in which case the time value is part of the file name.
!
!    Input, character ( len = * ) PFORM, specifies the form of the parallel
!    files.
!    'formatted', the files are formatted;
!    'unformatted' or any other value, the files are unformatted.
!
!    Input, character ( len = * ) SFORM, specifies the form of the
!    sequential file.
!    'formatted', the file is formatted;
!    'unformatted' or any other value, the file is unformatted.
!
  implicit none

  integer ( kind = 4 ) nsc
  integer ( kind = 4 ) nx_global
  integer ( kind = 4 ) nx_local
  integer ( kind = 4 ) ny_global
  integer ( kind = 4 ) ny_local

  integer ( kind = 4 ), parameter :: BASE_UNIT = 10
  real dummy(1)
  character ( len = 100 ) file_name
  character ( len = 100 ) file_name_fore
  character ( len = 100 ) file_name_mid
  character ( len = 100 ) file_name_post
  integer ( kind = 4 ) i_status
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) npx
  integer ( kind = 4 ) npy
  integer ( kind = 4 ) nval
  integer ( kind = 4 ) nz
  character ( len = * ) pform
  real qz1_global(nx_global,ny_global,nsc+4)
  real qz1_local(nx_local,ny_local,nsc+4)
  character ( len = * ) sform
  real temp_global(nx_global,ny_global)
  real temp_local(nx_local,ny_local)
  real time
!
!  Construct the family of names of the parallel SAVE files.
!
  if ( i_status == 0 ) then

    file_name_fore = 'field_t_'

    if ( pform == 'formatted' ) then
      write ( file_name_post, '(''_'',1pe9.3,''.txt'')' ) time
    else
      write ( file_name_post, '(''_'',1pe9.3)' ) time
    end if

  else if ( i_status == 1 ) then

    file_name_fore = 'field_l_'

    if ( pform == 'formatted' ) then
      file_name_post = '_.txt'
    else
      file_name_post = '_'
    end if

  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REJOIN_SAVE - Fatal error!'
    write ( *, '(a,i6)' ) '  Unexpected value of I_STATUS = ', i_status
    stop
  end if

  file_name_mid  = '0000'
!
!  Open all the parallel SAVE files.
!
  call pfile_open_read ( npx, npy, file_name_fore, file_name_mid, &
    file_name_post, BASE_UNIT, ierror, pform )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REJOIN_SAVE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input parallel SAVE files.'
    stop
  end if
!
!  Open one output SAVE file.
!
  if ( i_status == 0 ) then

    if ( sform == 'formatted' ) then
      write ( file_name, '(''field_t_'',1pe9.3,''.txt'' )' ) time
    else
      write ( file_name, '(''field_t_'',1pe9.3)' ) time
    end if

  else if ( i_status == 1 ) then

    if ( sform == 'formatted' ) then
      file_name = 'field_l.txt'
    else
      file_name = 'field_l'
    end if

  end if

  iunit = BASE_UNIT - 1

  call sfile_open_write ( file_name, iunit, sform, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REJOIN_SAVE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output sequential SAVE file!'
    stop
  end if
!
!  The first record is a single value, TIME, which is the same
!  in every file.  Read this value from one file, and advance them all.
!
  nval = 1

  call pfile_read_sdata ( npx, npy, BASE_UNIT, dummy, nval, ierror, pform )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REJOIN_SAVE - Fatal error!'
    write ( *, '(a)' ) '  End of file in input parallel SAVE files!'
    stop
  end if
!
!  Write the first record.
!
  call sfile_write ( iunit, sform, nval, dummy )
!
!  The second record is an array, QZ1, which needs to be reassembled.
!
  nz = nsc + 4

  call pfile_read_pdata ( npx, npy, nx_global, ny_global, nx_local, ny_local, &
    nz, BASE_UNIT, qz1_global, qz1_local, ierror, pform )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REJOIN_SAVE - Fatal error!'
    write ( *, '(a)' ) '  End of file in input parallel SAVE files!'
    stop
  end if
!
!  Write the second record.
!
  nval = nx_global * ny_global * ( nsc + 4 )

  call sfile_write ( iunit, sform, nval, qz1_global )
!
!  The third record is an array, TEMP, which needs to be reassembled.
!
  nz = 1

  call pfile_read_pdata ( npx, npy, nx_global, ny_global, nx_local, ny_local, &
    nz, BASE_UNIT, temp_global, temp_local, ierror, pform )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REJOIN_SAVE - Fatal error!'
    write ( *, '(a)' ) '  End of file in input parallel SAVE files!'
    stop
  end if
!
!  Write the third record.
!
  nval = nx_global * ny_global

  call sfile_write ( iunit, sform, nval, temp_global )
!
!  Close all the parallel files.
!
  call pfile_close ( npx, npy, BASE_UNIT )
!
!  Close the sequential file.
!
  call sfile_close ( iunit )

  return
end
subroutine sfile_close ( iunit )

!*****************************************************************************80
!
!! SFILE_CLOSE closes a sequential file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit number associated
!    with the file.
!
  implicit none

  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit

  close ( unit = iunit, iostat = ios )

  return
end
subroutine sfile_open_read ( file_name, iunit, sform, ierror )

!*****************************************************************************80
!
!! SFILE_OPEN_READ opens a sequential file for reading.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to be opened.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit number associated
!    with the file.
!
!    Input, character ( len = * ) SFORM, specifies the form of the
!    sequential file.
!    'formatted', the file is formatted;
!    'unformatted' or any other value, the file is unformatted.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error;
!    1, could not open the file.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = * ) sform

  ierror = 0

  open ( unit = iunit, file = file_name, form = sform, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SFILE_OPEN_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the sequential file for reading:'
    write ( *, '(a)' ) trim ( file_name )
    return
  end if

  return
end
subroutine sfile_open_write ( file_name, iunit, sform, ierror )

!*****************************************************************************80
!
!! SFILE_OPEN_WRITE opens a sequential file for writing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to be opened.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit number associated
!    with the file.
!
!    Input, character ( len = * ) SFORM, specifies the form of the
!    sequential file.
!    'formatted', the file is formatted;
!    'unformatted' or any other value, the file is unformatted.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error;
!    1, could not open the file.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = * ) sform

  ierror = 0

  open ( unit = iunit, file = file_name, form = sform, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SFILE_OPEN_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the sequential file for writing:'
    write ( *, '(a)' ) trim ( file_name )
    return
  end if

  return
end
subroutine sfile_read ( iunit, sform, nval, a, ierror )

!*****************************************************************************80
!
!! SFILE_READ reads data from a sequential file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit number associated
!    with the file.
!
!    Input, character ( len = * ) SFORM, specifies the form of the
!    sequential file.
!    'formatted', the file is formatted;
!    'unformatted' or any other value, the file is unformatted.
!
!    Input, integer ( kind = 4 ) NVAL, the number of entries in A.
!
!    Output, real A(NVAL), the data read from the file.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    nonzero, an error occurred.
!
  implicit none

  integer ( kind = 4 ) nval

  real a(nval)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = * ) sform

  ierror = 0

  if ( sform == 'formatted' ) then
    read ( iunit, *, iostat = ios ) a(1:nval)
  else
    read ( iunit, iostat = ios ) a
  end if

  if ( ios /= 0 ) then
    ierror = 1
    return
  endif

  return
end
subroutine sfile_write ( iunit, sform, nval, a )

!*****************************************************************************80
!
!! SFILE_WRITE writes data to the sequential file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit number associated
!    with the file.
!
!    Input, character ( len = * ) SFORM, specifies the form of the
!    sequential file.
!    'formatted', the file is formatted;
!    'unformatted' or any other value, the file is unformatted.
!
!    Input, integer ( kind = 4 ) NVAL, the number of entries in A.
!
!    Input, real A(NVAL), the data to be written.
!
  implicit none

  integer ( kind = 4 ) nval

  real a(nval)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  character ( len = * ) sform

  if ( sform == 'formatted' ) then
    write ( iunit, * ) a(1:nval)
  else
    write ( iunit ) a
  end if

  return
end
subroutine split_save ( i_status, npx, npy, nsc, nx_global, ny_global, &
  nx_local, ny_local, time, pform, sform )

!*****************************************************************************80
!
!! SPLIT_SAVE "splits" a single sequential DNS SAVE file into parcels.
!
!  Discussion:
!
!    This routine can handle either a FIELD_T_... or a FIELD_L_... file.
!
!    The sequential DNS SAVE file may be formatted or unformatted.
!
!    The unformatted DNS SAVE file contains three records:
!
!      1) TIME
!      2) ( ( ( QZ1(I,J,K), I=1,?), J=1,?), K=1,NSC+4 )
!      3) ( ( TEMP(I,J), I = 1, ?), J=1,?)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_STATUS is:
!    0, to convert the FIELD_T_... files;
!    1, to convert the FIELD_L_... files.
!
!    Input, integer ( kind = 4 ) NPX, NPY, the dimension of the processor
!    grid in the X and Y directions.
!
!    Input, integer ( kind = 4 ) NSC, used to define the third dimension
!    of QZ1_LOCAL.
!
!    Input, integer ( kind = 4 ) NX_GLOBAL, NY_GLOBAL, the number of grid
!    points in the X and Y directions, for the global grid.
!
!    Input, integer ( kind = 4 ) NX_LOCAL, NY_LOCAL, the number of grid
!    points in the X and Y directions, for the local grid associated with a
!    single processor.
!
!    Input, real TIME, the time, which is only needed if I_STATUS = 0,
!    in which case the time value is part of the file name.
!
!    Input, character ( len = * ) PFORM, specifies the form of the parallel
!    files.
!    'formatted', the files are formatted;
!    'unformatted' or any other value, the files are unformatted.
!
!    Input, character ( len = * ) SFORM, specifies the form of the sequential
!    file.
!    'formatted', the file is formatted;
!    'unformatted' or any other value, the file is unformatted.
!
  implicit none

  integer ( kind = 4 ) nsc
  integer ( kind = 4 ) nx_global
  integer ( kind = 4 ) nx_local
  integer ( kind = 4 ) ny_global
  integer ( kind = 4 ) ny_local

  integer ( kind = 4 ), parameter :: BASE_UNIT = 10
  real dummy(1)
  character ( len = 100 ) file_name
  character ( len = 100 ) file_name_fore
  character ( len = 100 ) file_name_mid
  character ( len = 100 ) file_name_post
  integer ( kind = 4 ) i_status
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) npx
  integer ( kind = 4 ) npy
  integer ( kind = 4 ) nval
  integer ( kind = 4 ) nz
  character ( len = * ) pform
  real qz1_global(nx_global,ny_global,nsc+4)
  real qz1_local(nx_local,ny_local,nsc+4)
  character ( len = * ) sform
  real temp_global(nx_global,ny_global)
  real temp_local(nx_local,ny_local)
  real time
!
!  Construct the family of names of the parallel SAVE files.
!
  if ( i_status == 0 ) then

    file_name_fore = 'field_t_'

    if ( pform == 'formatted' ) then
      write ( file_name_post, '(''_'',1pe9.3,''.txt'')' ) time
    else
      write ( file_name_post, '(''_'',1pe9.3)' ) time
    end if

  else if ( i_status == 1 ) then

    file_name_fore = 'field_l_'

    if ( pform == 'formatted' ) then
      file_name_post = '_.txt'
    else
      file_name_post = '_'
    end if

  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLIT_SAVE - Fatal error!'
    write ( *, '(a,i6)' ) '  Unexpected value of I_STATUS = ', i_status
    stop
  end if

  file_name_mid  = '0000'
!
!  Open the sequential input SAVE file.
!
  if ( i_status == 0 ) then

    if ( sform == 'formatted' ) then
      write ( file_name, '(''field_t_'',1pe9.3,''.txt'' )' ) time
    else
      write ( file_name, '(''field_t_'',1pe9.3)' ) time
    end if

  else if ( i_status == 1 ) then

    if ( sform == 'formatted' ) then
      file_name = 'field_l.txt'
    else
      file_name = 'field_l'
    end if

  end if

  iunit = BASE_UNIT - 1

  call sfile_open_read ( file_name, iunit, sform, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLIT_SAVE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the sequential input SAVE file.'
    stop
  end if
!
!  Open all the parallel output SAVE files.
!
  call pfile_open_write ( npx, npy, file_name_fore, file_name_mid, &
    file_name_post, BASE_UNIT, ierror, pform )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLIT_SAVE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the parallel output SAVE files.'
    stop
  end if
!
!  Read the first record.
!
  nval = 1

  call sfile_read ( iunit, sform, nval, dummy, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLIT_SAVE - Fatal error!'
    write ( *, '(a)' ) '  End of file in input sequential SAVE files!'
    stop
  end if
!
!  Write a copy of the first record to every parallel file.
!
  call pfile_write_sdata ( npx, npy, BASE_UNIT, dummy, nval, ierror, pform )
!
!  Read the second record, the array QZ1.
!
  nz = nsc + 4
  nval = nx_global * ny_global * nz

  call sfile_read ( iunit, sform, nval, qz1_global, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLIT_SAVE - Fatal error!'
    write ( *, '(a)' ) '  End of file in input sequential SAVE files!'
    stop
  end if
!
!  Parcel out the second record.
!
  call pfile_write_pdata ( npx, npy, nx_global, ny_global, nx_local, &
    ny_local, nz, BASE_UNIT, qz1_global, qz1_local, ierror, pform )
!
!  Read the third record.
!
  nz = 1
  nval = nx_global * ny_global * nz

  call sfile_read ( iunit, sform, nval, temp_global, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLIT_SAVE - Fatal error!'
    write ( *, '(a)' ) '  End of file in input sequential SAVE files!'
    stop
  end if
!
!  Parcel out the third record.
!
  call pfile_write_pdata ( npx, npy, nx_global, ny_global, nx_local, &
    ny_local, nz, BASE_UNIT, temp_global, temp_local, ierror, pform )
!
!  Close the sequential file.
!
  call sfile_close ( iunit )
!
!  Close all the parallel files.
!
  call pfile_close ( npx, npy, BASE_UNIT )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
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
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine u_to_f_all ( i_status, npx, npy, nsc, nx_local, ny_local, time )

!*****************************************************************************80
!
!! U_TO_F_ALL converts all the files from unformatted to formatted form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_STATUS is:
!    0, to convert the FIELD_T_... files;
!    1, to convert the FIELD_L_... files.
!
!    Input, integer ( kind = 4 ) NPX, NPY, the dimension of the processor
!    grid in the X and Y directions.
!
!    Input, integer ( kind = 4 ) NSC, used to define the third dimension
!    of QZ1_LOCAL.
!
!    Input, integer ( kind = 4 ) NX_LOCAL, NY_LOCAL, the number of grid
!    points in the X and Y directions, for the local grid associated with a
!    single processor.
!
!    Input, real TIME, the time associated with the FIELD_T files.
!    This is only needed when I_STATUS = 0, in order to form the names
!    of the FIELD_T_ files.
!
  implicit none

  integer ( kind = 4 ) nsc
  integer ( kind = 4 ) nx_local
  integer ( kind = 4 ) ny_local

  character ( len = 100 ) file_name_fore
  character ( len = 100 ) file_name_in
  character ( len = 100 ) file_name_mid
  character ( len = 100 ) file_name_out
  character ( len = 100 ) file_name_post
  integer ( kind = 4 ) i_proc
  integer ( kind = 4 ) i_status
  integer ( kind = 4 ) j_proc
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) npx
  integer ( kind = 4 ) npy
  integer ( kind = 4 ) proc_no
  real time
!
!  Construct the family of names of the parallel SAVE files.
!
  if ( i_status == 0 ) then
    file_name_fore = 'field_t_'
    write ( file_name_post, '(''_'',1pe9.3)' ) time
  else if ( i_status == 1 ) then
    file_name_fore = 'field_l_'
    file_name_post = '_'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'U_TO_F_ALL - Fatal error!'
    write ( *, '(a,i6)' ) '  Unexpected value of I_STATUS = ', i_status
    stop
  end if

  file_name_mid  = '0000'

  proc_no = 0

  proc_no = 0

  do j_proc = 1, npy

    do i_proc = 1, npx

      call pfile_name_gen ( file_name_fore, file_name_mid, file_name_post, &
        proc_no, file_name_in )

      lenc = len_trim ( file_name_in )

      file_name_out = ' '
      file_name_out(1:lenc) = file_name_in(1:lenc)
      file_name_out(lenc+1:lenc+4) = '.txt'

      write ( *, '(a)' ) 'Converting ' // trim (  file_name_in )

      call u_to_f_one ( file_name_in, file_name_out, nsc, nx_local, &
        ny_local )

      proc_no = proc_no + 1

    end do

  end do

  return
end
subroutine u_to_f_one ( file_name_in, file_name_out, nsc, nx_local, ny_local )

!*****************************************************************************80
!
!! U_TO_F_ONE converts one file from unformatted to formatted form.
!
!  Discussion:
!
!    This routine is used for debugging.  The unformatted data created
!    by a run of the parallel version of DNS can be converted to a formatted
!    form, for visual examination, or use by a program on another
!    architecture.
!
!    The conversion from unformatted to formatted data must be
!    carried out on the same architecture used to generate the original
!    unformatted data (or else you must be lucky or clever).
!
!    The unformatted file comprises three records:
!
!    1) TIME;
!    2) QZ1_LOCAL(NX_LOCAL,NY_LOCAL,NSC+4);
!    3) TEMP_LOCAL(NX_LOCAL,NY_LOCAL).
!
!    The formatted file contains the same information,
!    but with one (scalar) value per line of the file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME_IN, FILE_NAME_OUT, the names of the
!    unformatted input file to be read, and the formatted output file
!    to be created.
!
!    Input, integer ( kind = 4 ) NSC, used to define the third dimension
!    of QZ1_LOCAL.
!
!    Input, integer ( kind = 4 ) NX_LOCAL, NY_LOCAL, the number of grid
!    points in the X and Y directions, for the local grid associated with
!    a single processor.
!
  implicit none

  integer ( kind = 4 ) nsc
  integer ( kind = 4 ) nx_local
  integer ( kind = 4 ) ny_local

  character ( len = * ) file_name_in
  character ( len = * ) file_name_out
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: INPUT_UNIT = 1
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), parameter :: OUTPUT_UNIT = 2
  real qz1_local(nx_local,ny_local,nsc+4)
  real temp_local(nx_local,ny_local)
  real time
!
!  Open the input file.
!
  open ( unit = INPUT_UNIT, file = file_name_in, form = 'unformatted', &
    status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'U_TO_F_ONE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    close ( unit = INPUT_UNIT )
    return
  end if
!
!  Open a new output file.
!
  open ( unit = OUTPUT_UNIT, file = file_name_out, form = 'formatted', &
    status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'U_TO_F_ONE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    close ( unit = INPUT_UNIT )
    close ( unit = OUTPUT_UNIT )
    return
  end if
!
!  Transfer data from the input file to the output file.
!
  read ( INPUT_UNIT, iostat = ios ) time

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'U_TO_F_ONE - Fatal error!'
    write ( *, '(a)' ) '  Unexpected end of information on input file.'
    write ( *, '(a)' ) '  The output file will be DEFECTIVE.'
    close ( unit = INPUT_UNIT )
    close ( unit = OUTPUT_UNIT )
    return
  end if

  write ( OUTPUT_UNIT, * ) time

  read ( INPUT_UNIT, iostat = ios ) qz1_local

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'U_TO_F_ONE - Fatal error!'
    write ( *, '(a)' ) '  Unexpected end of information on input file.'
    write ( *, '(a)' ) '  The output file will be DEFECTIVE.'
    close ( unit = INPUT_UNIT )
    close ( unit = OUTPUT_UNIT )
    return
  end if

  do i = 1, nx_local
    do j = 1, ny_local
      do k = 1, nsc+4
        write ( OUTPUT_UNIT, * ) qz1_local(i,j,k)
      end do
    end do
  end do

  read ( INPUT_UNIT, iostat = ios ) temp_local

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'U_TO_F_ONE - Fatal error!'
    write ( *, '(a)' ) '  Unexpected end of information on input file.'
    write ( *, '(a)' ) '  The output file will be DEFECTIVE.'
    close ( unit = INPUT_UNIT )
    close ( unit = OUTPUT_UNIT )
    return
  end if

  do i = 1, nx_local
    do j = 1, ny_local
      write ( OUTPUT_UNIT, * ) temp_local(i,j)
    end do
  end do
!
!  Close the files.
!
  close ( unit = INPUT_UNIT )
  close ( unit = OUTPUT_UNIT )

  return
end
