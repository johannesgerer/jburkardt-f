subroutine face_print ( i, j, k, nx, ny, nz, face, iunit )

!*****************************************************************************80
!
!! FACE_PRINT prints the nodes that define one face of a voxel.
!
!  Discussion:
!
!    FACE_PRINT prints a face using the format appropriate for an OBJ
!    file, and is, in fact, a utility routine for I4VOXEL_TO_OBJ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, K, the indices of a voxel.
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions
!    of the voxel data.
!
!    Input, character ( len = 2 ) FACE, specifies the voxel face.
!    '-I' is the lower I=constant side, '+I' is the upper I=constant side;
!    '-J', '+J', '-K' and '+K' are also valid.
!
!    Input, integer ( kind = 4 ) IUNIT, the output unit number.
!
  implicit none

  character ( len = 2 ) face
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n000
  integer ( kind = 4 ) n001
  integer ( kind = 4 ) n010
  integer ( kind = 4 ) n011
  integer ( kind = 4 ) n100
  integer ( kind = 4 ) n101
  integer ( kind = 4 ) n110
  integer ( kind = 4 ) n111
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  call voxel_nodes ( i, j, k, nx, ny, nz, n000, n001, n010, n011, n100, &
    n101, n110, n111 )

  if ( face == '+I' ) then
    write ( iunit, '(a,2x,i7,2x,i7,2x,i7,2x,i7)' ) 'f', n100, n110, n111, n101
  else if ( face == '-I' ) then
    write ( iunit, '(a,2x,i7,2x,i7,2x,i7,2x,i7)' ) 'f', n000, n001, n011, n010
  else if ( face == '+J' ) then
    write ( iunit, '(a,2x,i7,2x,i7,2x,i7,2x,i7)' ) 'f', n010, n011, n111, n110
  else if ( face == '-J' ) then
    write ( iunit, '(a,2x,i7,2x,i7,2x,i7,2x,i7)' ) 'f', n000, n100, n101, n001
  else if ( face == '+K' ) then
    write ( iunit, '(a,2x,i7,2x,i7,2x,i7,2x,i7)' ) 'f', n001, n101, n111, n011
  else if ( face == '-K' ) then
    write ( iunit, '(a,2x,i7,2x,i7,2x,i7,2x,i7)' ) 'f', n000, n010, n110, n100
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

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
subroutine i4voxel_bound_print ( nx, ny, nz, i4voxel, iunit, num_face )

!*****************************************************************************80
!
!! I4VOXEL_BOUND_PRINT writes bounding faces to an OBJ file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions
!    of the voxel data.
!
!    Input, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of voxel data.
!
!    Input, integer ( kind = 4 ) IUNIT, the output unit number.
!
!    Output, integer ( kind = 4 ) NUM_FACE, the number of faces printed.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) num_face

  num_face = 0
!
!  I-planes
!
  do k = 1, nz
    do j = 1, ny
      do i = 0, nx

        if ( i == 0 ) then

          if ( i4voxel(i+1,j,k) /= 0 ) then
            call face_print ( i + 1, j, k, nx, ny, nz, '-I', iunit )
            num_face = num_face + 1
          end if

        else if ( i < nx ) then

          if ( i4voxel(i,j,k) == 0 .and. i4voxel(i+1,j,k) /= 0 ) then
            call face_print ( i+1, j, k, nx, ny, nz, '-I', iunit )
            num_face = num_face + 1
          else if ( i4voxel(i,j,k) /= 0 .and. i4voxel(i+1,j,k) == 0 ) then
            call face_print ( i, j, k, nx, ny, nz, '+I', iunit )
            num_face = num_face + 1
          end if

        else if ( i == nx ) then

          if ( i4voxel(i,j,k) /= 0 ) then
            call face_print ( i, j, k, nx, ny, nz, '+I', iunit )
            num_face = num_face + 1
          end if

        end if

      end do
    end do
  end do
!
!  J-planes
!
  do k = 1, nz
    do i = 1, nx
      do j = 0, ny

        if ( j == 0 ) then

          if ( i4voxel(i,j+1,k) /= 0 ) then
            call face_print ( i, j+1, k, nx, ny, nz, '-J', iunit )
            num_face = num_face + 1
          end if

        else if ( j < ny ) then

          if ( i4voxel(i,j,k) == 0 .and. i4voxel(i,j+1,k) /= 0 ) then
            call face_print ( i, j+1, k, nx, ny, nz, '-J', iunit )
            num_face = num_face + 1
          else if ( i4voxel(i,j,k) /= 0 .and. i4voxel(i,j+1,k) == 0 ) then
            call face_print ( i, j, k, nx, ny, nz, '+J', iunit )
            num_face = num_face + 1
          end if

        else if ( j == ny ) then

          if ( i4voxel(i,j,k) /= 0 ) then
            call face_print ( i, j, k, nx, ny, nz, '+J', iunit )
            num_face = num_face + 1
          end if

        end if

      end do
    end do
  end do
!
!  K-planes
!
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz

        if ( k == 0 ) then

          if ( i4voxel(i,j,k+1) /= 0 ) then
            call face_print ( i, j, k+1, nx, ny, nz, '-K', iunit )
            num_face = num_face + 1
          end if

        else if ( k < nz ) then

          if ( i4voxel(i,j,k) == 0 .and. i4voxel(i,j,k+1) /= 0 ) then
            call face_print ( i, j, k+1, nx, ny, nz, '-K', iunit )
            num_face = num_face + 1
          else if ( i4voxel(i,j,k) /= 0 .and. i4voxel(i,j,k+1) == 0 ) then
            call face_print ( i, j, k, nx, ny, nz, '+K', iunit )
            num_face = num_face + 1
          end if

        else if ( k == nz ) then

          if ( i4voxel(i,j,k) /= 0 ) then
            call face_print ( i, j, k, nx, ny, nz, '+K', iunit )
            num_face = num_face + 1
          end if

        end if

      end do
    end do
  end do

  return
end
subroutine i4voxel_count_positive ( nx, ny, nz, i4voxel, num_pos )

!*****************************************************************************80
!
!! I4VOXEL_COUNT_POSITIVE counts the positive entries in a voxel array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions
!    of the voxel data.
!
!    Input, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of voxel data.
!
!    Output, integer ( kind = 4 ) NUM_POS, the number of strictly positive
!    entries in the voxel array.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) num_pos
  integer ( kind = 4 ) i4voxel(nx,ny,nz)

  num_pos = 0

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz

        if ( 0 < i4voxel(i,j,k) ) then
          num_pos = num_pos + 1
        end if

      end do
    end do
  end do

  return
end
subroutine i4voxel_plot ( nx, ny, nz, i4voxel )

!*****************************************************************************80
!
!! I4VOXEL_PLOT prints out a typewriter plot of the Z slices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions
!    of the voxel data.
!
!    Input, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of voxel data.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 64 ) string

  do k = 1, nz

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Z plane ', k

    do j = 1, ny

      do i = 1, nx

        if ( i4voxel(i,j,k) == 0 ) then
          string(i:i)=' '
        else if ( i4voxel(i,j,k) == 158 ) then
          string(i:i) = '.'
        else
          string(i:i) = '*'
        end if

      end do

      write ( *, '(a)' ) trim ( string )

    end do

  end do

  return
end
subroutine i4voxel_plot2 ( nx, ny, nz, i4voxel )

!*****************************************************************************80
!
!! I4VOXEL_PLOT2 prints out a typewriter plot of the 3D regions.
!
!  Discussion:
!
!    NZ individual plots are made of NX by NY images.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions
!    of the voxel data.
!
!    Input, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of voxel data.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 68 ) string

  do k = 1, nz

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Z plane ', k
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  0    5   10   15   20   25   30   35   40   45   50   55   60   65'
    write ( *, '(a)' ) &
      '  +----+----+----+----+----+----+----+----+----+----+----+----+----+'

    do j = 1, ny

      if ( mod ( j, 5 ) == 0 ) then
        string = &
          '|----+----+----+----+----+----+----+----+----+----+----+----+----|'
      else
        string = &
          '|    |    |    |    |    |    |    |    |    |    |    |    |    |'
      end if

      do i = 1, nx

        if ( i4voxel(i,j,k) /= 0 ) then
          write ( string(i:i), '(i1)' ) i4voxel(i,j,k)
        end if

      end do

      if ( mod ( j, 5 ) == 0 ) then
        write ( *, '(i2,a)' ) j, string
      else
        write ( *, '(2x,a)' ) string
      end if

    end do

    write ( *, '(a)' ) &
      '  +----+----+----+----+----+----+----+----+----+----+----+----+----+'
    write ( *, '(a)' ) &
      '  0    5   10   15   20   25   30   35   40   45   50   55   60   65'

  end do

  return
end
subroutine i4voxel_plot3 ( nx, ny, nz, i4voxel )

!*****************************************************************************80
!
!! I4VOXEL_PLOT3 prints out a typewriter plot of data/100.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions
!    of the voxel data.
!
!    Input, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of voxel data.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 64 ) string

  do k = 1, nz

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'Z plane ', k

    do j = 1, ny

      string = ' '

      do i = 1, nx

        if ( i4voxel(i,j,k) /= 0 ) then
          write ( string(i:i), '(i1)' ) i4voxel(i,j,k) / 100
        end if

      end do

      write ( *, '(a)' ) trim ( string )

    end do

  end do

  return
end
subroutine i4voxel_read ( nx, ny, nz, i4voxel, filename )

!*****************************************************************************80
!
!! I4VOXEL_READ reads MRI data from an ASCII file, one item per line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions
!    of the voxel data.
!
!    Output, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of voxel data.
!
!    Input, character ( len = * ) FILENAME, the name of the file to be read.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  character ( len = 80 ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrec
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', iostat = ios, &
    form = 'formatted', access = 'sequential' )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VOXEL_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    stop
  end if
!
!  Read the data.
!
  nrec = 0

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx

        read ( iunit, *, iostat = ios ) i4voxel(i,j,k)

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4VOXEL_READ - Fatal error!'
          write ( *, '(a,i6)' ) '  END or ERR reading record ', nrec+1
          stop
        end if

        nrec = nrec + 1

      end do
    end do
  end do
!
!  Close the file.
!
  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'I4VOXEL_READ read ', nrec, ' records.'

  return
end
subroutine i4voxel_sum ( nx, ny, nz, i4voxel, vsum )

!*****************************************************************************80
!
!! I4VOXEL_SUM sums the entries in a voxel array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions
!    of the voxel data.
!
!    Input, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of voxel data.
!
!    Output, integer ( kind = 4 ) VSUM, the sum of the entries in the
!    voxel array.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) vsum

  vsum = sum ( i4voxel(1:nx,1:ny,1:nz) )

  return
end
subroutine i4voxel_thicken ( nx, ny, nz, i4voxel )

!*****************************************************************************80
!
!! I4VOXEL_THICKEN "thickens" the voxels.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions
!    of the voxel data.
!
!    Input/output, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of
!    voxel data.  On output, the data has been "thickened".
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz

        if ( 0 < i4voxel(i,j,k) ) then

          ilo = max ( i - 1, 1 )
          ihi = min ( i + 1, nx )
          jlo = max ( j - 1, 1 )
          jhi = min ( j + 1, nx )
          klo = max ( k - 1, 1 )
          khi = min ( k + 1, nx )

          do i2 = ilo, ihi
            do j2 = jlo, jhi
              do k2 = klo, khi
                if ( i4voxel(i2,j2,k2) == 0 ) then
                  i4voxel(i2,j2,k2) = - i4voxel(i,j,k)
                end if
              end do
            end do
          end do

        end if

      end do
    end do
  end do
!
!  The new voxels were marked with negative values, so that we
!  didn't allow yet more voxels to grow from THEM.
!  Now that we're done, we can reset the values in these voxels
!  to positive values, so they join their region.
!
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz

        if ( i4voxel(i,j,k) < 0 ) then
          i4voxel(i,j,k) = - i4voxel(i,j,k)
        end if

      end do
    end do
  end do

  return
end
subroutine i4voxel_thresh ( nx, ny, nz, i4voxel, thresh )

!*****************************************************************************80
!
!! I4VOXEL_THRESH zeroes out array entries below a given threshhold.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions
!    of the voxel data.
!
!    Input/output, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array
!    of voxel data.  On output, values below the threshhold value have
!    been reset to 0.
!
!    Input, integer ( kind = 4 ) THRESH, the threshhold value.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nkeep
  integer ( kind = 4 ) nthresh
  integer ( kind = 4 ) thresh

  nkeep = 0
  nthresh = 0

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz

        if ( 0 < i4voxel(i,j,k) ) then

          if ( i4voxel(i,j,k) < thresh ) then

            nthresh = nthresh + 1
            i4voxel(i,j,k) = 0

          else

            nkeep = nkeep + 1

          end if

        end if

      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4VOXEL_THRESH'
  write ( *, '(a,i6)' ) '  Elements kept             ', nkeep
  write ( *, '(a,i6)' ) '  Elements zeroed out       ', nthresh
  write ( *, '(a,i6)' ) '  which were below THRESH = ', thresh

  return
end
subroutine i4voxel_to_obj ( nx, ny, nz, i4voxel, filename )

!*****************************************************************************80
!
!! I4VOXEL_TO_OBJ writes out an OBJ file from a voxel array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions of
!    the voxel data.
!
!    Input, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of voxel data.
!
!    Input, character ( len = * ) FILENAME, the name of the file to be created.
!    An extension of '.obj' is recommended.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  character ( len = * ) filename
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) num_face
  integer ( kind = 4 ) num_node

  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'replace', &
    form = 'formatted', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VOXEL_TO_OBJ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    stop
  end if
!
!  Print the header.
!
  write ( iunit, '(a)' ) &
    '# "' // trim ( filename ) // '", created by i4voxel_to_obj.'
!
!  Print the node position data.
!
  write ( iunit, '(a)' ) ' '
  call node_print ( nx, ny, nz, iunit )
!
!  Print the face data.
!
  write ( iunit, '(a)' ) ' '

  call i4voxel_bound_print ( nx, ny, nz, i4voxel, iunit, num_face )
!
!  Close the file.
!
  close ( unit = iunit )
!
!  Report.
!
  num_node = ( nx + 1 ) * ( ny + 1 ) * ( nz + 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4VOXEL_TO_OBJ'
  write ( *, '(a)' ) '  Created OBJ file "' // trim ( filename ) // '".'
  write ( *, '(a,i6)' ) '  Number of nodes: ', num_node
  write ( *, '(a,i6)' ) '  Number of faces: ', num_face

  return
end
subroutine i4voxel_to_region ( nx, ny, nz, i4voxel, list, maxlist, nlist, &
  nregion )

!*****************************************************************************80
!
!! I4VOXEL_TO_REGION arranges a set of voxels into contiguous regions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the number of voxels in the
!    X, Y and Z directions.
!
!    Input/output, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ).
!    * On input, I4VOXEL(I,J,K) has the values:
!      0, if the voxel is OFF;
!      anything else, if the voxel is ON.
!    * On output, I4VOXEL(I,J,K) has the values:
!      0, if the voxel is off,
!      N, if the voxel is ON, and is part of region N.
!
!    Output, integer ( kind = 4 ) LIST(MAXLIST), contains, in stack form, a list
!    of the indices of the elements in each region.
!    * The number of elements in NREGION is NELEM = LIST(NLIST).  The
!      (I,J,K) indices of the last element in this region are in
!      LIST(NLIST-3) through LIST(NLIST-1), and the first element is
!      listed in LIST(NLIST-3*NELEM), LIST(NLIST-3*NELEM+1),
!      LIST(NLIST-3*NELEM+2).
!    * The number of elements in NREGION-1 is listed in
!      LIST(NLIST-3*NELEM-1), and so on.
!
!    Input, integer ( kind = 4 ) MAXLIST, the maximum length of the array used
!    to list the elements of the regions.
!
!    Output, integer ( kind = 4 ) NLIST, the number of entries of LIST that
!    were used.  However, if MAXLIST < NLIST, then there was not enough space in
!    LIST to store the data properly, and LIST should not be used,
!    although the data in I4VOXEL should be correct.
!
!    Output, integer ( kind = 4 ) NREGION, the number of regions discovered.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxstack = 5000
  integer ( kind = 4 ) maxlist
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ibase
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jbase
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) kbase
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) list(maxlist)
  integer ( kind = 4 ) nabes
  integer ( kind = 4 ) ncan
  integer ( kind = 4 ) nelements
  integer ( kind = 4 ) nlist
  integer ( kind = 4 ) nregion
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) stack(maxstack)
!
!  Reset all nonzero entries of I4VOXEL to -1.
!
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz

        if ( i4voxel(i,j,k) /= 0 ) then
          i4voxel(i,j,k) = -1
        end if

      end do
    end do
  end do
!
!  Start the number of items in the region list at 0.
!
  nlist = 0
!
!  Start the number of regions at 0.
!
  nregion = 0
!
!  The stack begins empty.
!
  nstack = 0
!
!  Search for an unused "ON" voxel from which we can "grow" a new region.
!
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
!
!  We found a voxel that is "ON", and does not belong to any region.
!
        if ( i4voxel(i,j,k) == -1 ) then
!
!  Increase the number of regions.
!
          nregion = nregion + 1
!
!  Add this voxel to the region.
!
          i4voxel(i,j,k) = nregion
!
!  Add this voxel to the stack.
!
          if ( maxstack < nstack + 4 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'I4VOXEL_TO_REGION - Fatal error!'
            write ( *, '(a)' ) '  The internal stack overflowed.'
            write ( *, '(a)' ) '  The algorithm has failed.'
            stop
          end if

          stack(nstack+1) = i
          stack(nstack+2) = j
          stack(nstack+3) = k

          stack(nstack+4) = 1

          nstack = nstack + 4
!
!  Add this voxel to the description of the region.
!
          nelements = 1

          if ( nlist + 3 <= maxlist ) then
            list(nlist+1) = i
            list(nlist+2) = j
            list(nlist+3) = k
          end if

          nlist = nlist + 3

10            continue
!
!  Find all neighbors of BASE that are "ON" but unused.
!  Mark them as belonging to this region, and stack their indices.
!
          ibase = stack(nstack-3)
          jbase = stack(nstack-2)
          kbase = stack(nstack-1)

          ilo = max ( ibase-1, 1 )
          ihi = min ( ibase+1, nx )
          jlo = max ( jbase-1, 1 )
          jhi = min ( jbase+1, ny )
          klo = max ( kbase-1, 1 )
          khi = min ( kbase+1, nz )

          nabes = 0

          do i2 = ilo, ihi
            do j2 = jlo, jhi
              do k2 = klo, khi
!
!  We found a neighbor to our current search point, which is "ON" and unused.
!
                if ( i4voxel(i2,j2,k2) == -1 ) then
!
!  Increase the number of neighbors.
!
                  nabes = nabes + 1
!
!  Mark the neighbor as belonging to the region.
!
                  i4voxel(i2,j2,k2) = nregion
!
!  Add the neighbor to the stack.
!
                  if ( maxstack < nstack+3 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'I4VOXEL_TO_REGION - Fatal error!'
                    write ( *, '(a)' ) '  The internal stack overflowed.'
                    write ( *, '(a)' ) '  The algorithm has failed.'
                    stop
                  end if

                  stack(nstack+1) = i2
                  stack(nstack+2) = j2
                  stack(nstack+3) = k2

                  nstack = nstack+3
!
!  Add the neighbor to the description of the region.
!
                  nelements = nelements + 1

                  if ( nlist+3 <= maxlist ) then
                    list(nlist+1) = i2
                    list(nlist+2) = j2
                    list(nlist+3) = k2
                  end if

                  nlist = nlist + 3

                end if

              end do
            end do
          end do
!
!  If any new neighbors were found, take the last one as the basis
!  for a deeper search.
!
          if ( 0 < nabes ) then

            if ( maxstack < nstack+1 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'I4VOXEL_TO_REGION - Fatal error!'
              write ( *, '(a)' ) '  The internal stack overflowed.'
              write ( *, '(a)' ) '  The algorithm has failed.'
              stop
            end if

            stack(nstack+1) = nabes
            nstack = nstack + 1
            go to 10

          end if
!
!  If the current search point had no new neighbors, drop it from the stack.
!
          ncan = stack(nstack) - 1
          nstack = nstack - 3
          stack(nstack) = ncan
!
!  If there are still any unused candidates at this level, take the
!  last one as the basis for a deeper search.
!
          if ( 0 < stack(nstack) ) then
            go to 10
          end if
!
!  If there are no more unused candidates at this level, then we need
!  to back up a level in the stack.  If there are any candidates at
!  that earlier level, then we can still do more searching.
!
          nstack = nstack - 1

          if ( 0 < nstack ) then
            go to 10
          end if
!
!  If we have exhausted the stack, we have completed this region.
!  Tag the number of elements to the end of the region description list.
!
          nlist = nlist + 1
          if ( nlist <= maxlist ) then
            list(nlist) = nelements
          end if

        end if

      end do
    end do
  end do
!
!  Print some warnings.
!
  if ( maxlist < nlist ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VOXEL_TO_REGION - Warning!'
    write ( *, '(a)' ) '  MAXLIST was too small to list the regions.'
    write ( *, '(a)' ) '  Do not try to use the LIST array!'
    write ( *, '(a)' ) '  The I4VOXEL data is OK, however.'
  end if

  return
end
subroutine i4voxel_write ( nx, ny, nz, i4voxel, filename )

!*****************************************************************************80
!
!! I4VOXEL_WRITE writes MRI data to an ASCII file, one item per line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions of the
!    voxel data.
!
!    Input, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of voxel data.
!
!    Input, character ( len = * ) FILENAME, the name of the file to be created.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrec
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'replace', iostat = ios, &
    form = 'formatted', access = 'sequential' )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VOXEL_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    stop
  end if

  nrec = 0

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        write ( iunit, * ) i4voxel(i,j,k)
        nrec = nrec + 1
      end do
    end do
  end do

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4VOXEL_WRITE'
  write ( *, '(a,i6)' ) '  Number of records written: ', nrec

  return
end
subroutine node_print ( nx, ny, nz, iunit )

!*****************************************************************************80
!
!! NODE_PRINT prints the nodes that define one face of a voxel.
!
!  Discussion:
!
!    NODE_PRINT prints the node information using the format appropriate
!    for an OBJ file, and is, in fact, a utility routine for I4VOXEL_TO_OBJ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions of the
!   voxel data.
!
!    Input, integer ( kind = 4 ) IUNIT, the output unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) z
  real ( kind = 8 ) zmax

  xmax = real ( nx, kind = 8 )
  ymax = real ( ny, kind = 8 )
  zmax = real ( nz, kind = 8 )

  do k = 0, nz
    z = zmax * real ( k  ) / real ( nz, kind = 8 )
    do j = 0, ny
      y = ymax * real ( j ) / real ( ny, kind = 8 )
      do i = 0, nx
        x = xmax * real ( i ) / real ( nx, kind = 8 )
        write ( iunit, '(a,2x,f8.4,2x,f8.4,2x,f8.4)' ) 'v', x, y, z
      end do
    end do
  end do

  return
end
subroutine r8voxel_to_i4voxel ( nx, ny, nz, r8voxel, i4voxel )

!*****************************************************************************80
!
!! R8VOXEL_TO_I4VOXEL copies R8VOXEL data into I4VOXEL data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions
!    of the voxel data.
!
!    Input, real ( kind = 8 ) R8VOXEL(NX,NY,NZ), an array of real voxel data.
!
!    Output, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), a copy of the real data,
!    rounded using the nearest-integer function NINT.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  real ( kind = 8 ) r8voxel(nx,ny,nz)

  i4voxel(1:nx,1:ny,1:nz) = nint ( r8voxel(1:nx,1:ny,1:nz) )

  return
end
subroutine region_blank ( nx, ny, nz, i4voxel, center, iregion, max_region, &
  nregion )

!*****************************************************************************80
!
!! REGION_BLANK zeroes out voxels in a particular numbered region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions of
!    the voxel data.
!
!    Input, integer ( kind = 4 ) MAX_REGION, the maximum number of regions
!    for which storage has been allocated.
!
!    Input, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of voxel data.
!
  implicit none

  integer ( kind = 4 ) max_region
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) center(4,max_region)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) iregion
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nregion

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz

        if ( i4voxel(i,j,k) == iregion ) then
          i4voxel(i,j,k) = 0
        else if ( iregion < i4voxel(i,j,k) ) then
          i4voxel(i,j,k) = i4voxel(i,j,k) - 1
        end if

      end do
    end do
  end do
!
!  Shift the CENTER data down one index.
!
  do j = iregion+1, nregion
    do i = 1, 4
      center(i,j-1) = center(i,j)
    end do
  end do
!
!  Update the number of regions.
!
  nregion = nregion - 1

  return
end
subroutine region_center ( nx, ny, nz, i4voxel, max_region, nregion, center )

!*****************************************************************************80
!
!! REGION_CENTER computes the centers of mass of the regions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions of
!    the voxel data.
!
!    Input, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of voxel data.
!
!    Input, integer ( kind = 4 ) MAX_REGION, the maximum number of regions
!    for which storage has been allocated.
!
!    Input, integer ( kind = 4 ) NREGION, the number of regions.
!
!    Output, real ( kind = 8 ) CENTER(4,MAX_REGION), contains the coordinates
!    of the center of mass of region I in CENTER(1,I), CENTER(2,I), CENTER(3,I),
!    and the number of voxels in region I in CENTER(4,I).
!
  implicit none

  integer ( kind = 4 ) max_region
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) center(4,max_region)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) iregion
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nregion

  center(1:4,1:max_region) = 0
!
!  Add each (I,J,K) to its region's total.
!
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz

        iregion = i4voxel(i,j,k)

        if ( iregion /= 0 ) then

          if ( 0 < iregion .and. iregion < max_region ) then
            center(1,iregion) = center(1,iregion) + i
            center(2,iregion) = center(2,iregion) + j
            center(3,iregion) = center(3,iregion) + k
            center(4,iregion) = center(4,iregion) + 1
          end if

        end if

      end do
    end do
  end do
!
!  Now normalize each center of mass by the number of voxels in the
!  region.  We round to the nearest integer.
!
  do iregion = 1, nregion
    if ( center(4,iregion) /= 0 ) then
      do i = 1, 3
        a = real ( center(i,iregion), kind = 8 )
        b = real ( center(4,iregion), kind = 8 )
        center(i,iregion) = nint ( a / b )
      end do
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CMASS:'
  write ( *, '(a)' ) '  Region, Center(I,J,K), Voxels:'
  write ( *, '(a)' ) ' '

  do iregion = 1, nregion
    write ( *, '(i6,4i6)' ) iregion, center(1:4,iregion)
  end do

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
subroutine transport ( nx, ny, nz, rvoxel, c, center, i4voxel, max_region )

!*****************************************************************************80
!
!! TRANSPORT transports voxels to the boundary, counts intermediate hits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions of
!    the voxel data.
!
!    Output, real ( kind = 8 ) RVOXEL(NX,NY,NZ), ...
!
!    Input, integer ( kind = 4 ) C(3), ...
!
!    Input, integer ( kind = 4 ) CENTER(3), ...
!
!    Input, integer ( kind = 4 ) I4VOXEL(NX,NY,NZ), an array of voxel data.
!
!    Input, integer ( kind = 4 ) MAX_REGION, the maximum number of regions
!    for which storage has been allocated.
!
  implicit none

  integer ( kind = 4 ) max_region
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) a(3)
  integer ( kind = 4 ) b(3)
  integer ( kind = 4 ) c(3)
  integer ( kind = 4 ) center(4,max_region)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) iregion
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) percent
  real ( kind = 8 ) rvoxel(nx,ny,nz)

  rvoxel(1:nx,1:ny,1:nz) = 0.0D+00

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz

        if ( i4voxel(i,j,k) /= 0 ) then

          iregion = i4voxel(i,j,k)
          a(1) = i
          a(2) = j
          a(3) = k
          b(1) = center(1,iregion)
          b(2) = center(2,iregion)
          b(3) = center(3,iregion)
          percent = 100.0D+00 / real ( center(4,iregion), kind = 8 )

          call transvox ( nx, ny, nz, a, rvoxel, b, c, iregion, percent )

        end if

      end do
    end do
  end do

  return
end
subroutine transvox ( nx, ny, nz, a, rvoxel, b, c, iregion, percent )

!*****************************************************************************80
!
!! TRANSVOX transports one voxel to the boundary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions of
!    the voxel data.
!
!    Input, integer ( kind = 4 ) A(3), ...
!
!    Input/output, real ( kind = 8 ) RVOXEL(NX,NY,NZ), ...
!
!    Input, integer ( kind = 4 ) B(3), ...
!
!    Input, integer ( kind = 4 ) C(3), ...
!
!    Input, integer ( kind = 4 ) IREGION, ...
!
!    Input, real ( kind = 8 ) PERCENT, ...
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) a(3)
  integer ( kind = 4 ) b(3)
  integer ( kind = 4 ) c(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) iregion
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jnc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) knc
  real ( kind = 8 ) percent
  real ( kind = 8 ) rvoxel(nx,ny,nz)

  i = a(1)
  j = a(2)
  k = a(3)
!
!  Mark the image of the voxel.
!
  if ( rvoxel(i,j,k) == 0.0D+00 ) then
    rvoxel(i,j,k) = 100.0D+00 * real ( iregion, kind = 8 )
  end if

  rvoxel(i,j,k) = rvoxel(i,j,k) + percent

  inc = b(1) - c(1)
  jnc = b(2) - c(2)
  knc = b(3) - c(3)
!
!  Starting at (I,J,K), transport the voxel out to the boundary.
!
  i2 = i
  j2 = j
  k2 = k

  do

    call voxel_index_step ( i, j, k, i2, j2, k2, inc, jnc, knc )

    if ( i2 < 1 .or. nx < i2 .or. &
         j2 < 1 .or. ny < j2 .or. &
         k2 < 1 .or. nz < k2 ) then
      exit
    end if

    if ( rvoxel(i2,j2,k2) == 0.0D+00 ) then
      rvoxel(i2,j2,k2) = 100.0D+00 * real ( iregion, kind = 8 )
    end if

    rvoxel(i2,j2,k2) = rvoxel(i2,j2,k2) + percent

  end do

  return
end
subroutine voxel_index_step ( i1, j1, k1, i2, j2, k2, inc, jnc, knc )

!*****************************************************************************80
!
!! VOXEL_INDEX_STEP computes indices of voxels along a line from a given point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1, J1, K1, the coordinates of the base
!    voxel from which the line begins.
!
!    Input/output, integer ( kind = 4 ) I2, J2, K2.
!    * On input, these are the coordinates of the current voxel on
!      the line.  For the first call, these might be I1, J1 and K1.
!    * On output, these are the coordinates of the next voxel along
!      the line.
!
!    Input, integer ( kind = 4 ) INC, JNC, KNC, the increments to the voxels.
!    These values define the direction along which the line proceeds.
!    However, the voxels on the line will typically be incremented
!    by a fractional value of the vector (INC,JNC,KNC), and the
!    result is essentially rounded.
!    * If you input INC = JNC = KNC, then no movement is possible,
!      and none is made.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) alphai
  real ( kind = 8 ) alphaj
  real ( kind = 8 ) alphak
  real ( kind = 8 ), parameter :: big = 100000.0D+00
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jnc
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) knc
!
!  Assuming for the moment that (I,J,K) can take on real values,
!  points on the line have the form:
!
!    I = I1 + alpha * inc
!    J = J1 + alpha * jnc
!    K = K1 + alpha * knc
!
  if ( inc == 0 .and. jnc == 0 .and. knc == 0 ) then
    return
  end if

  alpha = 0.0D+00
!
!  Compute the smallest ALPHA that will change I2, J2 or K2 by +-0.5.
!
  if ( 0 < inc ) then
    alphai = real ( i2 - i1 + 0.5D+00, kind = 8 ) / real ( inc, kind = 8 )
  else if ( inc < 0 ) then
    alphai = real ( i2 - i1 - 0.5D+00, kind = 8 ) / real ( inc, kind = 8 )
  else
    alphai = big
  end if

  if ( 0 < jnc ) then
    alphaj = real ( j2 - j1 + 0.5D+00, kind = 8 ) / real ( jnc, kind = 8 )
  else if ( jnc < 0 ) then
    alphaj = real ( j2 - j1 - 0.5D+00, kind = 8 ) / real ( jnc, kind = 8 )
  else
    alphaj = big
  end if

  if ( 0 < knc ) then
    alphak = real ( k2 - k1 + 0.5D+00, kind = 8 ) / real ( knc, kind = 8 )
  else if ( knc < 0 ) then
    alphak = real ( k2 - k1 - 0.5D+00, kind = 8 ) / real ( knc, kind = 8 )
  else
    alphaj = big
  end if
!
!  The ALPHA of smallest positive magnitude represents the closest
!  next voxel.
!
  alpha = big

  if ( 0.0D+00 < alphai ) then
    alpha = min ( alpha, alphai )
  end if

  if ( 0.0D+00 < alphaj ) then
    alpha = min ( alpha, alphaj )
  end if

  if ( 0.0D+00 < alphak ) then
    alpha = min ( alpha, alphak )
  end if
!
!  Move to the new voxel.  Whichever index just made the half
!  step must be forced to take a whole step.
!
  if ( alpha == alphai ) then
    i2 = i2 + sign ( 1, inc )
    j2 = j1 + nint ( alpha * jnc )
    k2 = k1 + nint ( alpha * knc )
  else if ( alpha == alphaj ) then
    i2 = i1 + nint ( alpha * inc )
    j2 = j2 + sign ( 1, jnc )
    k2 = k1 + nint ( alpha * knc )
  else if ( alpha == alphak ) then
    i2 = i1 + nint ( alpha * inc )
    j2 = j1 + nint ( alpha * jnc )
    k2 = k2 + sign ( 1, knc )
  end if

  return
end
subroutine voxel_nodes ( i, j, k, nx, ny, nz, n000, n001, n010, &
  n011, n100, n101, n110, n111 )

!*****************************************************************************80
!
!! VOXEL_NODES returns the indices of the nodes of a voxel.
!
!  Diagram:
!
!                           n011-----n111
!                            /|      /|
!                           / |     / |
!   ^                     n001----n101|
!   K                      |  |    |  |
!   |                      |  |    |  |
!   T                      |  |    |  |
!   h                      |  |    |  |
!   i                      |  |    |  |
!   r                      |n010---|-n110
!   d                      | /     | /
!   |                      |/      |/
!   |  I Numbered first->-n000----n100
!   |                     /
!   |              J Numbered
!                   second
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, K, the indices of a voxel.
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the X, Y and Z dimensions
!    of the voxel data.
!
!    Output, integer ( kind = 4 ) N000, N001, N010, N011, N100, N101,
!    N110, N111, the indices of the 8 nodes associated with a voxel.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n000
  integer ( kind = 4 ) n001
  integer ( kind = 4 ) n010
  integer ( kind = 4 ) n011
  integer ( kind = 4 ) n100
  integer ( kind = 4 ) n101
  integer ( kind = 4 ) n110
  integer ( kind = 4 ) n111
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  n000 = ( k - 1 ) * ( ny + 1 ) * ( nx + 1 ) + ( j - 1 ) * ( nx + 1 ) + i

  n100 = n000 + 1
  n010 = n000 + ( nx + 1 )
  n110 = n010 + 1

  n001 = n000 + ( ny + 1 ) * ( nx + 1 )
  n101 = n001 + 1
  n011 = n001 + ( nx + 1 )
  n111 = n011 + 1

  return
end
