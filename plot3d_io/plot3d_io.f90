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
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
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
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
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
subroutine r4_b_1d_dim ( iunit, idim, ierror )

!*****************************************************************************80
!
!! R4_B_1D_DIM reads a binary 1D file for the dimension.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, the number of nodes in the X direction.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit

  ierror = 0

  read ( iunit, iostat = ios ) idim

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_1D_DIM - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_1d_dimn ( iunit, idim, nvar, ierror )

!*****************************************************************************80
!
!! R4_B_1D_DIMN reads a binary 1D file for the dimension and NVAR.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, the number of nodes in the X direction.
!
!    Output, integer ( kind = 4 ) NVAR, the number of F variables defined.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) nvar

  ierror = 0

  read ( iunit, iostat = ios ) idim, nvar

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_1D_DIMN - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_1d_f ( iunit, idim, nvar, f, ierror )

!*****************************************************************************80
!
!! R4_B_1D_F reads a binary 1D file for the F data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, the number of nodes in the X direction.
!
!    Input, integer ( kind = 4 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 4 ) F(IDIM,NVAR), the F data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) nvar

  real ( kind = 4 ) f(idim,nvar)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) l

  ierror = 0

  read ( iunit, iostat = ios ) (( f(i,l), i = 1, idim ), l = 1, nvar )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_1D_F - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_1d_ffile ( iunit, idim, nvar, f, ierror )

!*****************************************************************************80
!
!! R4_B_1D_FFILE reads a binary 1D F file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, the number of nodes in the X direction.
!
!    Output, integer ( kind = 4 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 4 ) F(IDIM,NVAR), the F data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) nvar

  real ( kind = 4 ) f(idim,nvar)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit

  ierror = 0

  call r4_b_1d_dimn ( iunit, idim, nvar, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_1d_f ( iunit, idim, nvar, f, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r4_b_1d_q ( iunit, idim, q, ierror )

!*****************************************************************************80
!
!! R4_B_1D_Q reads a binary 1D Q file for the Q data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, the number of nodes in the X direction.
!
!    Output, real ( kind = 4 ) Q(IDIM,3), the Q data, density, X momentum, and
!    stagnation energy per unit volume.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzer, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) l
  real ( kind = 4 ) q(idim,3)

  ierror = 0

  read ( iunit, iostat = ios ) (( q(i,l), i = 1, idim ), l = 1, 3 )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_1D_Q - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_1d_qfile ( iunit, idim, alpha, fsmach, re, time, q, ierror )

!*****************************************************************************80
!
!! R4_B_1D_QFILE reads a binary 1D Q file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, the number of nodes in the X direction.
!
!    Output, real ( kind = 4 ) ALPHA, the angle of attack, in degrees.
!
!    Output, real ( kind = 4 ) FSMACH, the free stream Mach number.
!
!    Output, real ( kind = 4 ) RE, the Reynolds number for the flow.
!
!    Output, real ( kind = 4 ) TIME, the time associated with the data.
!
!    Output, real ( kind = 4 ) Q(IDIM,3), the Q data, density, X momentum, and
!    stagnation energy per unit volume.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the parameters.
!    3, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim

  real ( kind = 4 ) alpha
  real ( kind = 4 ) fsmach
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 4 ) q(idim,3)
  real ( kind = 4 ) re
  real ( kind = 4 ) time

  ierror = 0

  call r4_b_1d_dim ( iunit, idim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_param ( iunit, fsmach, alpha, re, time, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  call r4_b_1d_q ( iunit, idim, q, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  return
end
subroutine r4_b_1d_x ( iunit, idim, x, ierror )

!*****************************************************************************80
!
!! R4_B_1D_X reads a binary 1D X file for the X data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, the number of nodes in the X direction.
!
!    Output, real ( kind = 4 ) X(IDIM), the X coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  real ( kind = 4 ) x(idim)

  ierror = 0

  read ( iunit, iostat = ios ) x(1:idim)

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_1D_XY - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_1d_xfile ( iunit, idim, x, ierror )

!*****************************************************************************80
!
!! R4_B_1D_XYFILE reads a binary 1D X file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, the number of nodes in the X direction.
!
!    Output, real ( kind = 4 ) X(IDIM), the X coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 4 ) x(idim)

  ierror = 0

  call r4_b_1d_dim ( iunit, idim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_1d_x ( iunit, idim, x, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r4_b_2d_dim ( iunit, idim, jdim, ierror )

!*****************************************************************************80
!
!! R4_B_2D_DIM reads a binary 2D file for the dimensions.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jdim

  ierror = 0

  read ( iunit, iostat = ios ) idim, jdim

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_2D_DIM - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_2d_dimn ( iunit, idim, jdim, nvar, ierror )

!*****************************************************************************80
!
!! R4_B_2D_DIMN reads a binary 2D file for the dimensions and NVAR.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) NVAR, the number of F variables defined.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) nvar

  ierror = 0

  read ( iunit, iostat = ios ) idim, jdim, nvar

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_2D_DIMN - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_2d_f ( iunit, idim, jdim, nvar, f, ierror )

!*****************************************************************************80
!
!! R4_B_2D_F reads a binary 2D F file for the F data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Input, integer ( kind = 4 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 4 ) F(IDIM,JDIM,NVAR), the F data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) nvar

  real ( kind = 4 ) f(idim,jdim,nvar)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  ierror = 0

  read ( iunit, iostat = ios ) ((( f(i,j,l), i = 1, idim ), j = 1, jdim ), &
    l = 1, nvar )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_2D_F - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_2d_ffile ( iunit, idim, jdim, nvar, f, ierror )

!*****************************************************************************80
!
!! R4_B_2D_FFILE reads a binary 2D F file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 4 ) F(IDIM,JDIM,NVAR), the F data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) nvar

  real ( kind = 4 ) f(idim,jdim,nvar)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit

  ierror = 0

  call r4_b_2d_dimn ( iunit, idim, jdim, nvar, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_2d_f ( iunit, idim, jdim, nvar, f, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r4_b_2d_q ( iunit, idim, jdim, q, ierror )

!*****************************************************************************80
!
!! R4_B_2D_Q reads a binary 2D Q file for the Q data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, real ( kind = 4 ) Q(IDIM,JDIM,4), the Q data, density, X and Y momentum,
!    and stagnation energy per unit volume.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  real ( kind = 4 ) q(idim,jdim,4)

  ierror = 0

  read ( iunit, iostat = ios ) ((( q(i,j,l), i = 1, idim ), j = 1, jdim ), &
    l = 1, 4 )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_2D_Q - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_2d_qfile ( iunit, idim, jdim, alpha, fsmach, re, time, q, ierror )

!*****************************************************************************80
!
!! R4_B_2D_QFILE reads a binary 2D Q file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, real ( kind = 4 ) ALPHA, the angle of attack, in degrees.
!
!    Output, real ( kind = 4 ) FSMACH, the free stream Mach number.
!
!    Output, real ( kind = 4 ) RE, the Reynolds number for the flow.
!
!    Output, real ( kind = 4 ) TIME, the time associated with the data.
!
!    Output, real ( kind = 4 ) Q(IDIM,JDIM,4), the Q data, density, X and Y momentum,
!    and stagnation energy per unit volume.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the parameters.
!    3, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim

  real ( kind = 4 ) alpha
  real ( kind = 4 ) fsmach
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 4 ) q(idim,jdim,4)
  real ( kind = 4 ) re
  real ( kind = 4 ) time

  ierror = 0

  call r4_b_2d_dim ( iunit, idim, jdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_param ( iunit, fsmach, alpha, re, time, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  call r4_b_2d_q ( iunit, idim, jdim, q, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  return
end
subroutine r4_b_2d_xy ( iunit, idim, jdim, x, y, ierror )

!*****************************************************************************80
!
!! R4_B_2D_XY reads a binary 2D XY file for the XY data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, real ( kind = 4 ) X(IDIM,JDIM), Y(IDIM,JDIM), the X and Y coordinates
!    of the nodes.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  real ( kind = 4 ) x(idim,jdim)
  real ( kind = 4 ) y(idim,jdim)

  ierror = 0

  read ( iunit, iostat = ios ) (( x(i,j), i = 1, idim ), j = 1, jdim ), &
                               (( y(i,j), i = 1, idim ), j = 1, jdim )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_2D_XY - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_2d_xyfile ( iunit, idim, jdim, x, y, ierror )

!*****************************************************************************80
!
!! R4_B_2D_XYFILE reads a binary 2D XY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, real ( kind = 4 ) X(IDIM,JDIM), Y(IDIM,JDIM), the X and Y coordinates
!    of the nodes.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 4 ) x(idim,jdim)
  real ( kind = 4 ) y(idim,jdim)

  ierror = 0

  call r4_b_2d_dim ( iunit, idim, jdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_2d_xy ( iunit, idim, jdim, x, y, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r4_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

!*****************************************************************************80
!
!! R4_B_3D_DIM reads a binary 3D file for the dimensions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  ierror = 0

  read ( iunit, iostat = ios ) idim, jdim, kdim

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_3D_DIM - Error!'
    write ( *, '(a)' ) '  A read error occurred.'
    return
  end if

  return
end
subroutine r4_b_3d_dimn ( iunit, idim, jdim, kdim, nvar, ierror )

!*****************************************************************************80
!
!! R4_B_3D_DIMN reads a binary 3D file the dimensions and NVAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, integer ( kind = 4 ) NVAR, the number of F variables defined.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim
  integer ( kind = 4 ) nvar

  ierror = 0

  read ( iunit, iostat = ios ) idim, jdim, kdim, nvar

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_3D_DIMN - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_3d_f ( iunit, idim, jdim, kdim, nvar, f, ierror )

!*****************************************************************************80
!
!! R4_B_3D_F reads a binary 3D F file for the F data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Input, integer ( kind = 4 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 4 ) F(IDIM,JDIM,KDIM,NVAR), the F data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim
  integer ( kind = 4 ) nvar

  real ( kind = 4 ) f(idim,jdim,kdim,nvar)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  ierror = 0

  read ( iunit, iostat = ios ) (((( f(i,j,k,l), i = 1, idim ), &
    j = 1, jdim ), k = 1, kdim ), l = 1, nvar )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_3D_F - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_3d_ffile ( iunit, idim, jdim, kdim, nvar, f, ierror )

!*****************************************************************************80
!
!! R4_B_3D_FFILE reads a binary 3D F file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, integer ( kind = 4 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 4 ) F(IDIM,JDIM,KDIM,NVAR), the F data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim
  integer ( kind = 4 ) nvar

  real ( kind = 4 ) f(idim,jdim,kdim,nvar)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit

  ierror = 0

  call r4_b_3d_dimn ( iunit, idim, jdim, kdim, nvar, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_3d_f ( iunit, idim, jdim, kdim, nvar, f, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r4_b_3d_q ( iunit, idim, jdim, kdim, q, ierror )

!*****************************************************************************80
!
!! R4_B_3D_Q reads a binary 3D Q file for the Q data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 4 ) Q(IDIM,JDIM,KDIM,5), the Q data, density, X, Y and Z
!    momentum, and stagnation energy per unit volume.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 4 ) q(idim,jdim,kdim,5)

  ierror = 0

  read ( iunit, iostat = ios ) (((( q(i,j,k,l), i = 1, idim ), &
    j = 1, jdim ), k = 1, kdim ), l = 1, 5 )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_3D_Q - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_3d_qfile ( iunit, idim, jdim, kdim, alpha, fsmach, re, time, &
  q, ierror )

!*****************************************************************************80
!
!! R4_B_3D_QFILE reads a binary 3D Q file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 4 ) ALPHA, the angle of attack, in degrees.
!
!    Output, real ( kind = 4 ) FSMACH, the free stream Mach number.
!
!    Output, real ( kind = 4 ) RE, the Reynolds number for the flow.
!
!    Output, real ( kind = 4 ) TIME, the time associated with the data.
!
!    Output, real ( kind = 4 ) Q(IDIM,JDIM,KDIM,5), the Q data, density, X, Y and Z
!    momentum, and stagnation energy per unit volume.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the parameters.
!    3, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  real ( kind = 4 ) alpha
  real ( kind = 4 ) fsmach
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 4 ) q(idim,jdim,kdim,5)
  real ( kind = 4 ) re
  real ( kind = 4 ) time

  ierror = 0

  call r4_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_param ( iunit, fsmach, alpha, re, time, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  call r4_b_3d_q ( iunit, idim, jdim, kdim, q, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  return
end
subroutine r4_b_3d_xyz ( iunit, idim, jdim, kdim, x, y, z, ierror )

!*****************************************************************************80
!
!! R4_B_3D_XYZ reads a binary 3D XYZ file for the XYZ data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 4 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM), Z(IDIM,JDIM,KDIM),
!    the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) x(idim,jdim,kdim)
  real ( kind = 4 ) y(idim,jdim,kdim)
  real ( kind = 4 ) z(idim,jdim,kdim)

  ierror = 0

  read ( iunit, iostat = ios ) &
    ((( x(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim ), &
    ((( y(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim ), &
    ((( z(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_3D_XYZ - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_3d_xyzb ( iunit, idim, jdim, kdim, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R4_B_3D_XYZB reads a binary 3D XYZB file for the XYZB data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 4 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM), Z(IDIM,JDIM,KDIM),
!    the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) B(IDIM,JDIM,KDIM), the blanking information.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  integer ( kind = 4 ) b(idim,jdim,kdim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) x(idim,jdim,kdim)
  real ( kind = 4 ) y(idim,jdim,kdim)
  real ( kind = 4 ) z(idim,jdim,kdim)

  ierror = 0

  read ( iunit, iostat = ios ) &
    ((( x(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim ), &
    ((( y(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim ), &
    ((( z(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim ), &
    ((( b(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim )

  return

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_3D_XYZB - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r4_b_3d_xyzbfile ( iunit, idim, jdim, kdim, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R4_B_3D_XYZBFILE reads a binary 3D XYZB file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 4 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM), Z(IDIM,JDIM,KDIM),
!    the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) B(IDIM,JDIM,KDIM), the blanking information.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  integer ( kind = 4 ) b(idim,jdim,kdim)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 4 ) x(idim,jdim,kdim)
  real ( kind = 4 ) y(idim,jdim,kdim)
  real ( kind = 4 ) z(idim,jdim,kdim)

  ierror = 0

  call r4_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_3d_xyzb ( iunit, idim, jdim, kdim, x, y, z, b, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r4_b_3d_xyzfile ( iunit, idim, jdim, kdim, x, y, z, ierror )

!*****************************************************************************80
!
!! R4_B_3D_XYZFILE reads a binary 3D XYZ file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 4 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM), Z(IDIM,JDIM,KDIM),
!    the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 4 ) x(idim,jdim,kdim)
  real ( kind = 4 ) y(idim,jdim,kdim)
  real ( kind = 4 ) z(idim,jdim,kdim)

  ierror = 0

  call r4_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_3d_xyz ( iunit, idim, jdim, kdim, x, y, z, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r4_b_3dm_dim ( iunit, idim, jdim, kdim, maxgrid, ngrid, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_DIM reads a binary 3D multiple grid file for the dimensions.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID), the number
!    of nodes in the X, Y, and Z directions for each grid.
!
!    Input, integer ( kind = 4 ) MAXGRID, the maximum value of NGRID.
!
!    Input, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) maxgrid

  integer ( kind = 4 ) idim(maxgrid)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jdim(maxgrid)
  integer ( kind = 4 ) kdim(maxgrid)
  integer ( kind = 4 ) ngrid

  ierror = 0

  read ( iunit, iostat = ios ) ( idim(igrid), jdim(igrid), &
    kdim(igrid), igrid = 1, ngrid )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_3DM_DIM - Error!'
    write ( *, '(a)' ) '  A read error occurred.'
    return
  end if

  return
end
subroutine r4_b_3dm_ngrid ( iunit, ngrid, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_NGRID reads a binary 3D multiple grid file for the number of grids.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ngrid

  ierror = 0

  read ( iunit, iostat = ios ) ngrid

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_3DM_DIM - Error!'
    write ( *, '(a)' ) '  A read error occurred.'
    return
  end if

  return
end
subroutine r4_b_3dm_q ( iunit, idim, jdim, kdim, maxi, maxj, maxk, maxgrid, &
  ngrid, fsmach, alpha, re, time, q, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_Q reads a binary 3D multiple grid Q file for the parameters and Q data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions in each grid.
!
!    Input, integer ( kind = 4 ) MAXI, MAXJ, MAXK, the maximum dimension for
!    the first three components of Q.
!
!    Input, integer ( kind = 4 ) MAXGRID, the maximum value of NGRID.
!
!    Input, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, real ( kind = 4 ) FSMACH(MAXGRID), ALPHA(MAXGRID), RE(MAXGRID),
!    TIME(MAXGRID), the values of the free stream Mach number, the angle of
!    attack, in degrees, the Reynolds number for the flow, and the time,
!    for each grid.
!
!    Output, real ( kind = 4 ) Q(MAXI,MAXJ,MAXK,5,MAXGRID),
!    the Q values of the nodes for each grid.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) maxgrid
  integer ( kind = 4 ) maxi
  integer ( kind = 4 ) maxj
  integer ( kind = 4 ) maxk

  real ( kind = 4 ) alpha(maxgrid)
  real ( kind = 4 ) fsmach(maxgrid)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idim(maxgrid)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jdim(maxgrid)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kdim(maxgrid)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ngrid
  real ( kind = 4 ) q(maxi,maxj,maxk,5,maxgrid)
  real ( kind = 4 ) re(maxgrid)
  real ( kind = 4 ) time(maxgrid)

  ierror = 0

  do igrid = 1, ngrid

    read ( iunit, iostat = ios ) fsmach(igrid), alpha(igrid), re(igrid), &
      time(igrid)

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_B_3DM_Q - Fatal error!'
      write ( *, '(a)' ) '  A read error occurred.'
      return
    end if

    read ( iunit, iostat = ios ) (((( q(i,j,k,l,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), l = 1, 5 )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_B_3DM_Q - Fatal error!'
      write ( *, '(a)' ) '  A read error occurred.'
      return
    end if

  end do

  return
end
subroutine r4_b_3dm_xyz ( iunit, idim, jdim, kdim, maxi, maxj, maxk, maxgrid, &
  ngrid, x, y, z, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_XYZ reads a binary 3D multiple grid XYZ file for the XYZ data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions in each grid.
!
!    Input, integer ( kind = 4 ) MAXI, MAXJ, MAXK, the maximum dimension for
!    the first three components of X, Y and Z.
!
!    Input, integer ( kind = 4 ) MAXGRID, the maximum value of NGRID.
!
!    Input, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, real ( kind = 4 ) X(MAXI,MAXJ,MAXK,MAXGRID),
!    Y(MAXI,MAXJ,MAXK,MAXGRID),
!    Z(MAXI,MAXJ,MAXK,MAXGRID), the X, Y and Z coordinates of the nodes
!    for each grid.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) maxgrid
  integer ( kind = 4 ) maxi
  integer ( kind = 4 ) maxj
  integer ( kind = 4 ) maxk

  integer ( kind = 4 ) i
  integer ( kind = 4 ) idim(maxgrid)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jdim(maxgrid)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kdim(maxgrid)
  integer ( kind = 4 ) ngrid
  real ( kind = 4 ) x(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) y(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) z(maxi,maxj,maxk,maxgrid)

  ierror = 0

  do igrid = 1, ngrid

    read ( iunit, iostat = ios ) ((( x(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( y(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( z(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_B_3DM_XYZ - Fatal error!'
      write ( *, '(a)' ) '  A read error occurred.'
      return
    end if

  end do

  return
end
subroutine r4_b_3dm_xyzb ( iunit, idim, jdim, kdim, maxi, maxj, maxk, maxgrid, &
  ngrid, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_XYZB reads a binary 3D multiple grid XYZB file for the XYZB data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions in each grid.
!
!    Input, integer ( kind = 4 ) MAXI, MAXJ, MAXK, the maximum dimension for
!    the first three components of X, Y and Z.
!
!    Input, integer ( kind = 4 ) MAXGRID, the maximum value of NGRID.
!
!    Input, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, real ( kind = 4 ) X(MAXI,MAXJ,MAXK,MAXGRID),
!    Y(MAXI,MAXJ,MAXK,MAXGRID),
!    Z(MAXI,MAXJ,MAXK,MAXGRID), the X, Y and Z coordinates of the nodes
!    for each grid.
!
!    Output, integer ( kind = 4 ) B(MAXI,MAXJ,MAXK,MAXGRID), the blanking
!    information.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) maxgrid
  integer ( kind = 4 ) maxi
  integer ( kind = 4 ) maxj
  integer ( kind = 4 ) maxk

  integer ( kind = 4 ) b(maxi,maxj,maxk,maxgrid)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idim(maxgrid)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jdim(maxgrid)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kdim(maxgrid)
  integer ( kind = 4 ) ngrid
  real ( kind = 4 ) x(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) y(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) z(maxi,maxj,maxk,maxgrid)

  ierror = 0

  do igrid = 1, ngrid

    read ( iunit, iostat = ios ) ((( x(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( y(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( z(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( b(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_B_3DM_XYZB - Fatal error!'
      write ( *, '(a)' ) '  A read error occurred.'
      return
    end if

  end do

  return
end
subroutine r4_b_3dm_qfile ( iunit, idim, jdim, kdim, maxi, maxj, maxk, maxgrid, &
  ngrid, fsmach, alpha, re, time, q, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_QFILE reads a binary 3D multiple grid Q file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions for each grid.
!
!    Input, integer ( kind = 4 ) MAXI, MAXJ, MAXK, the maximum dimension for
!    the first three components of Q.
!
!    Input, integer ( kind = 4 ) MAXGRID, the maximum value of NGRID.
!
!    Output, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, real ( kind = 4 ) FSMACH(MAXGRID), ALPHA(MAXGRID), RE(MAXGRID), TIME(MAXGRID),
!    the values of the free stream Mach number, the angle of attack, in degrees,
!    the Reynolds number for the flow, and the time, for each grid.
!
!    Output, real ( kind = 4 ) Q(MAXI,MAXJ,MAXK,5,MAXGRID),
!    the Q values of the nodes for each grid.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading NGRID.
!    2, MAXGRID < NGRID.
!    3, an error occurred while reading the dimensions.
!    4, some IDIM, JDIM or KDIM is greater than MAXI, MAXJ or MAXK.
!    5, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) maxgrid
  integer ( kind = 4 ) maxi
  integer ( kind = 4 ) maxj
  integer ( kind = 4 ) maxk

  real ( kind = 4 ) alpha(maxgrid)
  real ( kind = 4 ) fsmach(maxgrid)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idim(maxgrid)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jdim(maxgrid)
  integer ( kind = 4 ) kdim(maxgrid)
  integer ( kind = 4 ) ngrid
  real ( kind = 4 ) q(maxi,maxj,maxk,5,maxgrid)
  real ( kind = 4 ) re(maxgrid)
  real ( kind = 4 ) time(maxgrid)

  ierror = 0

  call r4_b_3dm_ngrid ( iunit, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  if ( maxgrid < ngrid ) then
    ierror = 2
    return
  end if

  call r4_b_3dm_dim ( iunit, idim, jdim, kdim, maxgrid, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  do i = 1, ngrid

    if ( maxi < idim(i) .or. maxj < jdim(i) .or. maxk < kdim(i) ) then
      ierror = 4
      return
    end if

  end do

  call r4_b_3dm_q ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
    maxgrid, ngrid, fsmach, alpha, re, time, q, ierror )

  if ( ierror /= 0 ) then
    ierror = 5
    return
  end if

  return
end
subroutine r4_b_3dm_xyzfile ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
  maxgrid, ngrid, x, y, z, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_XYZFILE reads a binary 3D multiple grid XYZ file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions for each grid.
!
!    Input, integer ( kind = 4 ) MAXI, MAXJ, MAXK, the maximum dimension for
!    the first three components of X, Y and Z.
!
!    Input, integer ( kind = 4 ) MAXGRID, the maximum value of NGRID.
!
!    Output, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, real ( kind = 4 ) X(MAXI,MAXJ,MAXK,MAXGRID),
!    Y(MAXI,MAXJ,MAXK,MAXGRID),
!    Z(MAXI,MAXJ,MAXK,MAXGRID), the X, Y and Z coordinates of the nodes
!    for each grid.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading NGRID.
!    2, MAXGRID < NGRID.
!    3, an error occurred while reading the dimensions.
!    4, some IDIM, JDIM or KDIM is greater than MAXI, MAXJ or MAXK.
!    5, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) maxgrid
  integer ( kind = 4 ) maxi
  integer ( kind = 4 ) maxj
  integer ( kind = 4 ) maxk

  integer ( kind = 4 ) i
  integer ( kind = 4 ) idim(maxgrid)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jdim(maxgrid)
  integer ( kind = 4 ) kdim(maxgrid)
  integer ( kind = 4 ) ngrid
  real ( kind = 4 ) x(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) y(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) z(maxi,maxj,maxk,maxgrid)

  ierror = 0

  call r4_b_3dm_ngrid ( iunit, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  if ( maxgrid < ngrid ) then
    ierror = 2
    return
  end if

  call r4_b_3dm_dim ( iunit, idim, jdim, kdim, maxgrid, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  do i = 1, ngrid

    if ( maxi < idim(i) .or. maxj < jdim(i) .or. maxk < kdim(i) ) then
      ierror = 4
      return
    end if

  end do

  call r4_b_3dm_xyz ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
    maxgrid, ngrid, x, y, z, ierror )

  if ( ierror /= 0 ) then
    ierror = 5
    return
  end if

  return
end
subroutine r4_b_3dm_xyzbfile ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
  maxgrid, ngrid, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_XYZBFILE reads a binary 3D multiple grid XYZB file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions for each grid.
!
!    Input, integer ( kind = 4 ) MAXI, MAXJ, MAXK, the maximum dimension for the first three
!    components of X, Y and Z.
!
!    Input, integer ( kind = 4 ) MAXGRID, the maximum value of NGRID.
!
!    Output, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, real ( kind = 4 ) X(MAXI,MAXJ,MAXK,MAXGRID),
!    Y(MAXI,MAXJ,MAXK,MAXGRID),
!    Z(MAXI,MAXJ,MAXK,MAXGRID), the X, Y and Z coordinates of the nodes
!    for each grid.
!
!    Output, integer ( kind = 4 ) B(MAXI,MAXJ,MAXK,MAXGRID), the blanking information.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading NGRID.
!    2, MAXGRID < NGRID.
!    3, an error occurred while reading the dimensions.
!    4, some IDIM, JDIM or KDIM is greater than MAXI, MAXJ or MAXK.
!    5, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) maxgrid
  integer ( kind = 4 ) maxi
  integer ( kind = 4 ) maxj
  integer ( kind = 4 ) maxk

  integer ( kind = 4 ) b(maxi,maxj,maxk,maxgrid)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idim(maxgrid)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jdim(maxgrid)
  integer ( kind = 4 ) kdim(maxgrid)
  integer ( kind = 4 ) ngrid
  real ( kind = 4 ) x(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) y(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) z(maxi,maxj,maxk,maxgrid)

  ierror = 0

  call r4_b_3dm_ngrid ( iunit, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  if ( maxgrid < ngrid ) then
    ierror = 2
    return
  end if

  call r4_b_3dm_dim ( iunit, idim, jdim, kdim, maxgrid, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  do i = 1, ngrid

    if ( maxi < idim(i) .or. maxj < jdim(i) .or. maxk < kdim(i) ) then
      ierror = 4
      return
    end if

  end do

  call r4_b_3dm_xyzb ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
    maxgrid, ngrid, x, y, z, b, ierror )

  if ( ierror /= 0 ) then
    ierror = 5
    return
  end if

  return
end
subroutine r4_b_3dp_f ( iunit, idim, jdim, kdim, nvar, f, ierror )

!*****************************************************************************80
!
!! R4_B_3DP_F reads a binary 3D plane F file for the F data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Input, integer ( kind = 4 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 4 ) F(IDIM,JDIM,KDIM,NVAR), the F data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim
  integer ( kind = 4 ) nvar

  real ( kind = 4 ) f(idim,jdim,kdim,nvar)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  ierror = 0

  do k = 1, kdim

    read ( iunit, iostat = ios ) &
      ((( f(i,j,k,l), i = 1, idim ), j = 1, jdim ), l = 1, nvar )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_B_3DP_F - Error!'
      write ( *, '(a)' ) '  An end-of-file condition occurred.'
      return
    end if

  end do

  return
end
subroutine r4_b_3dp_ffile ( iunit, idim, jdim, kdim, nvar, f, ierror )

!*****************************************************************************80
!
!! R4_B_3DP_FFILE reads a binary 3D plane F file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, integer ( kind = 4 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 4 ) F(IDIM,JDIM,KDIM,NVAR), the F data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim
  integer ( kind = 4 ) nvar

  real ( kind = 4 ) f(idim,jdim,kdim,nvar)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit

  ierror = 0

  call r4_b_3d_dimn ( iunit, idim, jdim, kdim, nvar, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_3dp_f ( iunit, idim, jdim, kdim, nvar, f, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r4_b_3dp_q ( iunit, idim, jdim, kdim, q, ierror )

!*****************************************************************************80
!
!! R4_B_3DP_Q reads a binary 3D plane Q file for the Q data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 4 ) Q(IDIM,JDIM,KDIM,5), the Q data, density, X, Y and Z
!    momentum, and stagnation energy per unit volume.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 4 ) q(idim,jdim,kdim,5)

  ierror = 0

  do k = 1, kdim

    read ( iunit, iostat = ios ) &
      ((( q(i,j,k,l), i = 1, idim ), j = 1, jdim ), l = 1, 5 )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_B_3DP_Q - Error!'
      write ( *, '(a)' ) '  An end-of-file condition occurred.'
      return
    end if

  end do

  return
end
subroutine r4_b_3dp_qfile ( iunit, idim, jdim, kdim, alpha, fsmach, re, time, &
  q, ierror )

!*****************************************************************************80
!
!! R4_B_3DP_QFILE reads a binary 3D plane Q file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 4 ) ALPHA, the angle of attack, in degrees.
!
!    Output, real ( kind = 4 ) FSMACH, the free stream Mach number.
!
!    Output, real ( kind = 4 ) RE, the Reynolds number for the flow.
!
!    Output, real ( kind = 4 ) TIME, the time associated with the data.
!
!    Output, real ( kind = 4 ) Q(IDIM,JDIM,KDIM,5), the Q data, density, X, Y and Z
!    momentum, and stagnation energy per unit volume.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the parameters.
!    3, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  real ( kind = 4 ) alpha
  real ( kind = 4 ) fsmach
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 4 ) q(idim,jdim,kdim,5)
  real ( kind = 4 ) re
  real ( kind = 4 ) time

  ierror = 0

  call r4_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_param ( iunit, fsmach, alpha, re, time, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  call r4_b_3dp_q ( iunit, idim, jdim, kdim, q, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  return
end
subroutine r4_b_3dp_xyz ( iunit, idim, jdim, kdim, x, y, z, ierror )

!*****************************************************************************80
!
!! R4_B_3DP_XYZ reads a binary 3D plane XYZ file for the XYZ data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 4 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM),
!    Z(IDIM,JDIM,KDIM), the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) x(idim,jdim,kdim)
  real ( kind = 4 ) y(idim,jdim,kdim)
  real ( kind = 4 ) z(idim,jdim,kdim)

  ierror = 0

  do k = 1, kdim

    read ( iunit, iostat = ios ) (( x(i,j,k), i = 1, idim ), j = 1, jdim ), &
      (( y(i,j,k), i = 1, idim ), j = 1, jdim ), &
      (( z(i,j,k), i = 1, idim ), j = 1, jdim )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_B_3DP_XYZ - Error!'
      write ( *, '(a)' ) '  An end-of-file condition occurred.'
      return
    end if

  end do

  return
end
subroutine r4_b_3dp_xyzb ( iunit, idim, jdim, kdim, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R4_B_3DP_XYZB reads a binary 3D plane XYZB file for the XYZB data.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 4 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM),
!    Z(IDIM,JDIM,KDIM), the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) B(IDIM,JDIM,KDIM), the blanking information.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  integer ( kind = 4 ) b(idim,jdim,kdim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) x(idim,jdim,kdim)
  real ( kind = 4 ) y(idim,jdim,kdim)
  real ( kind = 4 ) z(idim,jdim,kdim)

  ierror = 0

  do k = 1, kdim

    read ( iunit, iostat = ios ) (( x(i,j,k), i = 1, idim ), j = 1, jdim ), &
      (( y(i,j,k), i = 1, idim ), j = 1, jdim ), &
      (( z(i,j,k), i = 1, idim ), j = 1, jdim ), &
      (( b(i,j,k), i = 1, idim ), j = 1, jdim )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_B_3DP_XYZB - Error!'
      write ( *, '(a)' ) '  An end-of-file condition occurred.'
    end if

  end do

  return
end
subroutine r4_b_3dp_xyzbfile ( iunit, idim, jdim, kdim, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R4_B_3DP_XYZBFILE reads a binary 3D plane XYZB file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 4 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM),
!    Z(IDIM,JDIM,KDIM), the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) B(IDIM,JDIM,KDIM), the blanking information.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  integer ( kind = 4 ) b(idim,jdim,kdim)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 4 ) x(idim,jdim,kdim)
  real ( kind = 4 ) y(idim,jdim,kdim)
  real ( kind = 4 ) z(idim,jdim,kdim)

  ierror = 0

  call r4_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_3dp_xyzb ( iunit, idim, jdim, kdim, x, y, z, b, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r4_b_3dp_xyzfile ( iunit, idim, jdim, kdim, x, y, z, ierror )

!*****************************************************************************80
!
!! R4_B_3DP_XYZFILE reads a binary 3D plane XYZ file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 4 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM),
!    Z(IDIM,JDIM,KDIM), the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) idim
  integer ( kind = 4 ) jdim
  integer ( kind = 4 ) kdim

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 4 ) x(idim,jdim,kdim)
  real ( kind = 4 ) y(idim,jdim,kdim)
  real ( kind = 4 ) z(idim,jdim,kdim)

  ierror = 0

  call r4_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r4_b_3dp_xyz ( iunit, idim, jdim, kdim, x, y, z, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r4_b_param ( iunit, fsmach, alpha, re, time, ierror )

!*****************************************************************************80
!
!! R4_B_PARAM reads a binary Q file for the flow parameters.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, real ( kind = 4 ) FSMACH, the free stream Mach number.
!
!    Output, real ( kind = 4 ) ALPHA, the angle of attack, in degrees.
!
!    Output, real ( kind = 4 ) RE, the Reynolds number for the flow.
!
!    Output, real ( kind = 4 ) TIME, the time associated with the data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  real ( kind = 4 ) alpha
  real ( kind = 4 ) fsmach
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  real ( kind = 4 ) re
  real ( kind = 4 ) time

  ierror = 0

  read ( iunit, iostat = ios ) fsmach, alpha, re, time

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_PARAM - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_1d_dim ( iunit, idim, ierror )

!*****************************************************************************80
!
!! R8_B_1D_DIM reads a binary 1D file for the dimension.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, the number of nodes in the X direction.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit

  ierror = 0

  read ( iunit, iostat = ios ) idim

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_1D_DIM - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_1d_dimn ( iunit, idim, nvar, ierror )

!*****************************************************************************80
!
!! R8_B_1D_DIMN reads a binary 1D file for the dimension and NVAR.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, the number of nodes in the X direction.
!
!    Output, integer ( kind = 8 ) NVAR, the number of F variables defined.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) nvar

  ierror = 0

  read ( iunit, iostat = ios ) idim, nvar

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_1D_DIMN - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_1d_f ( iunit, idim, nvar, f, ierror )

!*****************************************************************************80
!
!! R8_B_1D_F reads a binary 1D file for the F data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, the number of nodes in the X direction.
!
!    Input, integer ( kind = 8 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 8 ) F(IDIM,NVAR), the F data.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) nvar

  real ( kind = 8 ) f(idim,nvar)
  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) l

  ierror = 0

  read ( iunit, iostat = ios ) (( f(i,l), i = 1, idim ), l = 1, nvar )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_1D_F - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_1d_ffile ( iunit, idim, nvar, f, ierror )

!*****************************************************************************80
!
!! R8_B_1D_FFILE reads a binary 1D F file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, the number of nodes in the X direction.
!
!    Output, integer ( kind = 8 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 8 ) F(IDIM,NVAR), the F data.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) nvar

  real ( kind = 8 ) f(idim,nvar)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit

  ierror = 0

  call r8_b_1d_dimn ( iunit, idim, nvar, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_1d_f ( iunit, idim, nvar, f, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r8_b_1d_q ( iunit, idim, q, ierror )

!*****************************************************************************80
!
!! R8_B_1D_Q reads a binary 1D Q file for the Q data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, the number of nodes in the X direction.
!
!    Output, real ( kind = 8 ) Q(IDIM,3), the Q data, density, X momentum, and
!    stagnation energy per unit volume.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzer, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim

  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) l
  real ( kind = 8 ) q(idim,3)

  ierror = 0

  read ( iunit, iostat = ios ) (( q(i,l), i = 1, idim ), l = 1, 3 )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_1D_Q - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_1d_qfile ( iunit, idim, alpha, fsmach, re, time, q, ierror )

!*****************************************************************************80
!
!! R8_B_1D_QFILE reads a binary 1D Q file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, the number of nodes in the X direction.
!
!    Output, real ( kind = 8 ) ALPHA, the angle of attack, in degrees.
!
!    Output, real ( kind = 8 ) FSMACH, the free stream Mach number.
!
!    Output, real ( kind = 8 ) RE, the Reynolds number for the flow.
!
!    Output, real ( kind = 8 ) TIME, the time associated with the data.
!
!    Output, real ( kind = 8 ) Q(IDIM,3), the Q data, density, X momentum, and
!    stagnation energy per unit volume.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the parameters.
!    3, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim

  real ( kind = 8 ) alpha
  real ( kind = 8 ) fsmach
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  real ( kind = 8 ) q(idim,3)
  real ( kind = 8 ) re
  real ( kind = 8 ) time

  ierror = 0

  call r8_b_1d_dim ( iunit, idim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_param ( iunit, fsmach, alpha, re, time, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  call r8_b_1d_q ( iunit, idim, q, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  return
end
subroutine r8_b_1d_x ( iunit, idim, x, ierror )

!*****************************************************************************80
!
!! R8_B_1D_X reads a binary 1D X file for the X data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, the number of nodes in the X direction.
!
!    Output, real ( kind = 8 ) X(IDIM), the X coordinates of the nodes.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim

  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  real ( kind = 8 ) x(idim)

  ierror = 0

  read ( iunit, iostat = ios ) x(1:idim)

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_1D_XY - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_1d_xfile ( iunit, idim, x, ierror )

!*****************************************************************************80
!
!! R8_B_1D_XYFILE reads a binary 1D X file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, the number of nodes in the X direction.
!
!    Output, real ( kind = 8 ) X(IDIM), the X coordinates of the nodes.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim

  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  real ( kind = 8 ) x(idim)

  ierror = 0

  call r8_b_1d_dim ( iunit, idim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_1d_x ( iunit, idim, x, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r8_b_2d_dim ( iunit, idim, jdim, ierror )

!*****************************************************************************80
!
!! R8_B_2D_DIM reads a binary 2D file for the dimensions.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) jdim

  ierror = 0

  read ( iunit, iostat = ios ) idim, jdim

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_2D_DIM - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_2d_dimn ( iunit, idim, jdim, nvar, ierror )

!*****************************************************************************80
!
!! R8_B_2D_DIMN reads a binary 2D file for the dimensions and NVAR.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, integer ( kind = 8 ) NVAR, the number of F variables defined.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) nvar

  ierror = 0

  read ( iunit, iostat = ios ) idim, jdim, nvar

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_2D_DIMN - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_2d_f ( iunit, idim, jdim, nvar, f, ierror )

!*****************************************************************************80
!
!! R8_B_2D_F reads a binary 2D F file for the F data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Input, integer ( kind = 8 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 8 ) F(IDIM,JDIM,NVAR), the F data.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) nvar

  real ( kind = 8 ) f(idim,jdim,nvar)
  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) l

  ierror = 0

  read ( iunit, iostat = ios ) ((( f(i,j,l), i = 1, idim ), j = 1, jdim ), &
    l = 1, nvar )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_2D_F - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_2d_ffile ( iunit, idim, jdim, nvar, f, ierror )

!*****************************************************************************80
!
!! R8_B_2D_FFILE reads a binary 2D F file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, integer ( kind = 8 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 8 ) F(IDIM,JDIM,NVAR), the F data.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) nvar

  real ( kind = 8 ) f(idim,jdim,nvar)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit

  ierror = 0

  call r8_b_2d_dimn ( iunit, idim, jdim, nvar, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_2d_f ( iunit, idim, jdim, nvar, f, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r8_b_2d_q ( iunit, idim, jdim, q, ierror )

!*****************************************************************************80
!
!! R8_B_2D_Q reads a binary 2D Q file for the Q data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, real ( kind = 8 ) Q(IDIM,JDIM,4), the Q data, density, X and
!    Y momentum, and stagnation energy per unit volume.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim

  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) l
  real ( kind = 8 ) q(idim,jdim,4)

  ierror = 0

  read ( iunit, iostat = ios ) ((( q(i,j,l), i = 1, idim ), j = 1, jdim ), &
    l = 1, 4 )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_2D_Q - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_2d_qfile ( iunit, idim, jdim, alpha, fsmach, re, time, q, ierror )

!*****************************************************************************80
!
!! R8_B_2D_QFILE reads a binary 2D Q file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, real ( kind = 8 ) ALPHA, the angle of attack, in degrees.
!
!    Output, real ( kind = 8 ) FSMACH, the free stream Mach number.
!
!    Output, real ( kind = 8 ) RE, the Reynolds number for the flow.
!
!    Output, real ( kind = 8 ) TIME, the time associated with the data.
!
!    Output, real ( kind = 8 ) Q(IDIM,JDIM,4), the Q data, density, X and
!    Y momentum, and stagnation energy per unit volume.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the parameters.
!    3, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim

  real ( kind = 8 ) alpha
  real ( kind = 8 ) fsmach
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  real ( kind = 8 ) q(idim,jdim,4)
  real ( kind = 8 ) re
  real ( kind = 8 ) time

  ierror = 0

  call r8_b_2d_dim ( iunit, idim, jdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_param ( iunit, fsmach, alpha, re, time, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  call r8_b_2d_q ( iunit, idim, jdim, q, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  return
end
subroutine r8_b_2d_xy ( iunit, idim, jdim, x, y, ierror )

!*****************************************************************************80
!
!! R8_B_2D_XY reads a binary 2D XY file for the XY data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, real ( kind = 8 ) X(IDIM,JDIM), Y(IDIM,JDIM), the X and Y coordinates
!    of the nodes.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim

  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  real ( kind = 8 ) x(idim,jdim)
  real ( kind = 8 ) y(idim,jdim)

  ierror = 0

  read ( iunit, iostat = ios ) (( x(i,j), i = 1, idim ), j = 1, jdim ), &
                               (( y(i,j), i = 1, idim ), j = 1, jdim )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_2D_XY - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_2d_xyfile ( iunit, idim, jdim, x, y, ierror )

!*****************************************************************************80
!
!! R8_B_2D_XYFILE reads a binary 2D XY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, the number of nodes in the
!    X and Y directions.
!
!    Output, real ( kind = 8 ) X(IDIM,JDIM), Y(IDIM,JDIM), the X and Y coordinates
!    of the nodes.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim

  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  real ( kind = 8 ) x(idim,jdim)
  real ( kind = 8 ) y(idim,jdim)

  ierror = 0

  call r8_b_2d_dim ( iunit, idim, jdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_2d_xy ( iunit, idim, jdim, x, y, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r8_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

!*****************************************************************************80
!
!! R8_B_3D_DIM reads a binary 3D file for the dimensions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  ierror = 0

  read ( iunit, iostat = ios ) idim, jdim, kdim

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_3D_DIM - Error!'
    write ( *, '(a)' ) '  A read error occurred.'
    return
  end if

  return
end
subroutine r8_b_3d_dimn ( iunit, idim, jdim, kdim, nvar, ierror )

!*****************************************************************************80
!
!! R8_B_3D_DIMN reads a binary 3D file the dimensions and NVAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, integer ( kind = 8 ) NVAR, the number of F variables defined.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim
  integer ( kind = 8 ) nvar

  ierror = 0

  read ( iunit, iostat = ios ) idim, jdim, kdim, nvar

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_3D_DIMN - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_3d_f ( iunit, idim, jdim, kdim, nvar, f, ierror )

!*****************************************************************************80
!
!! R8_B_3D_F reads a binary 3D F file for the F data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Input, integer ( kind = 8 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 8 ) F(IDIM,JDIM,KDIM,NVAR), the F data.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim
  integer ( kind = 8 ) nvar

  real ( kind = 8 ) f(idim,jdim,kdim,nvar)
  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  integer ( kind = 8 ) l

  ierror = 0

  read ( iunit, iostat = ios ) (((( f(i,j,k,l), i = 1, idim ), &
    j = 1, jdim ), k = 1, kdim ), l = 1, nvar )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_3D_F - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_3d_ffile ( iunit, idim, jdim, kdim, nvar, f, ierror )

!*****************************************************************************80
!
!! R8_B_3D_FFILE reads a binary 3D F file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, integer ( kind = 8 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 8 ) F(IDIM,JDIM,KDIM,NVAR), the F data.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim
  integer ( kind = 8 ) nvar

  real ( kind = 8 ) f(idim,jdim,kdim,nvar)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit

  ierror = 0

  call r8_b_3d_dimn ( iunit, idim, jdim, kdim, nvar, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_3d_f ( iunit, idim, jdim, kdim, nvar, f, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r8_b_3d_q ( iunit, idim, jdim, kdim, q, ierror )

!*****************************************************************************80
!
!! R8_B_3D_Q reads a binary 3D Q file for the Q data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 8 ) Q(IDIM,JDIM,KDIM,5), the Q data, density, X, Y and Z
!    momentum, and stagnation energy per unit volume.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  integer ( kind = 8 ) l
  real ( kind = 8 ) q(idim,jdim,kdim,5)

  ierror = 0

  read ( iunit, iostat = ios ) (((( q(i,j,k,l), i = 1, idim ), &
    j = 1, jdim ), k = 1, kdim ), l = 1, 5 )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_3D_Q - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_3d_qfile ( iunit, idim, jdim, kdim, alpha, fsmach, re, time, &
  q, ierror )

!*****************************************************************************80
!
!! R8_B_3D_QFILE reads a binary 3D Q file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 8 ) ALPHA, the angle of attack, in degrees.
!
!    Output, real ( kind = 8 ) FSMACH, the free stream Mach number.
!
!    Output, real ( kind = 8 ) RE, the Reynolds number for the flow.
!
!    Output, real ( kind = 8 ) TIME, the time associated with the data.
!
!    Output, real ( kind = 8 ) Q(IDIM,JDIM,KDIM,5), the Q data, density, X, Y and Z
!    momentum, and stagnation energy per unit volume.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the parameters.
!    3, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  real ( kind = 8 ) alpha
  real ( kind = 8 ) fsmach
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  real ( kind = 8 ) q(idim,jdim,kdim,5)
  real ( kind = 8 ) re
  real ( kind = 8 ) time

  ierror = 0

  call r8_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_param ( iunit, fsmach, alpha, re, time, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  call r8_b_3d_q ( iunit, idim, jdim, kdim, q, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  return
end
subroutine r8_b_3d_xyz ( iunit, idim, jdim, kdim, x, y, z, ierror )

!*****************************************************************************80
!
!! R8_B_3D_XYZ reads a binary 3D XYZ file for the XYZ data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 8 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM),
!    Z(IDIM,JDIM,KDIM), the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  real ( kind = 8 ) x(idim,jdim,kdim)
  real ( kind = 8 ) y(idim,jdim,kdim)
  real ( kind = 8 ) z(idim,jdim,kdim)

  ierror = 0

  read ( iunit, iostat = ios ) &
    ((( x(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim ), &
    ((( y(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim ), &
    ((( z(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_3D_XYZ - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_3d_xyzb ( iunit, idim, jdim, kdim, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R8_B_3D_XYZB reads a binary 3D XYZB file for the XYZB data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 8 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM),
!    Z(IDIM,JDIM,KDIM), the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 8 ) B(IDIM,JDIM,KDIM), the blanking information.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  integer ( kind = 8 ) b(idim,jdim,kdim)
  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  real ( kind = 8 ) x(idim,jdim,kdim)
  real ( kind = 8 ) y(idim,jdim,kdim)
  real ( kind = 8 ) z(idim,jdim,kdim)

  ierror = 0

  read ( iunit, iostat = ios ) &
    ((( x(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim ), &
    ((( y(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim ), &
    ((( z(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim ), &
    ((( b(i,j,k), i = 1, idim ), j = 1, jdim ), k = 1, kdim )

  return

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_3D_XYZB - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
  end if

  return
end
subroutine r8_b_3d_xyzbfile ( iunit, idim, jdim, kdim, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R8_B_3D_XYZBFILE reads a binary 3D XYZB file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 8 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM),
!    Z(IDIM,JDIM,KDIM), the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 8 ) B(IDIM,JDIM,KDIM), the blanking information.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  integer ( kind = 8 ) b(idim,jdim,kdim)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  real ( kind = 8 ) x(idim,jdim,kdim)
  real ( kind = 8 ) y(idim,jdim,kdim)
  real ( kind = 8 ) z(idim,jdim,kdim)

  ierror = 0

  call r8_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_3d_xyzb ( iunit, idim, jdim, kdim, x, y, z, b, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r8_b_3d_xyzfile ( iunit, idim, jdim, kdim, x, y, z, ierror )

!*****************************************************************************80
!
!! R8_B_3D_XYZFILE reads a binary 3D XYZ file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 8 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM),
!    Z(IDIM,JDIM,KDIM), the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  real ( kind = 8 ) x(idim,jdim,kdim)
  real ( kind = 8 ) y(idim,jdim,kdim)
  real ( kind = 8 ) z(idim,jdim,kdim)

  ierror = 0

  call r8_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_3d_xyz ( iunit, idim, jdim, kdim, x, y, z, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r8_b_3dm_dim ( iunit, idim, jdim, kdim, maxgrid, ngrid, ierror )

!*****************************************************************************80
!
!! R8_B_3DM_DIM reads a binary 3D multiple grid file for the dimensions.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions for each grid.
!
!    Input, integer ( kind = 8 ) MAXGRID, the maximum value of NGRID.
!
!    Input, integer ( kind = 8 ) NGRID, the number of grids.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) maxgrid

  integer ( kind = 8 ) idim(maxgrid)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) igrid
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) jdim(maxgrid)
  integer ( kind = 8 ) kdim(maxgrid)
  integer ( kind = 8 ) ngrid

  ierror = 0

  read ( iunit, iostat = ios ) ( idim(igrid), jdim(igrid), &
    kdim(igrid), igrid = 1, ngrid )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_3DM_DIM - Error!'
    write ( *, '(a)' ) '  A read error occurred.'
    return
  end if

  return
end
subroutine r8_b_3dm_ngrid ( iunit, ngrid, ierror )

!*****************************************************************************80
!
!! R8_B_3DM_NGRID reads a binary 3D multiple grid file for the number of grids.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) NGRID, the number of grids.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) ngrid

  ierror = 0

  read ( iunit, iostat = ios ) ngrid

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_3DM_DIM - Error!'
    write ( *, '(a)' ) '  A read error occurred.'
    return
  end if

  return
end
subroutine r8_b_3dm_q ( iunit, idim, jdim, kdim, maxi, maxj, maxk, maxgrid, &
  ngrid, fsmach, alpha, re, time, q, ierror )

!*****************************************************************************80
!
!! R8_B_3DM_Q reads a binary 3D multiple grid Q file for the parameters and Q data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions in each grid.
!
!    Input, integer ( kind = 8 ) MAXI, MAXJ, MAXK, the maximum dimension for the first three
!    components of Q.
!
!    Input, integer ( kind = 8 ) MAXGRID, the maximum value of NGRID.
!
!    Input, integer ( kind = 8 ) NGRID, the number of grids.
!
!    Output, real ( kind = 8 ) FSMACH(MAXGRID), ALPHA(MAXGRID), RE(MAXGRID), TIME(MAXGRID),
!    the values of the free stream Mach number, the angle of attack, in degrees,
!    the Reynolds number for the flow, and the time, for each grid.
!
!    Output, real ( kind = 8 ) Q(MAXI,MAXJ,MAXK,5,MAXGRID),
!    the Q values of the nodes for each grid.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) maxgrid
  integer ( kind = 8 ) maxi
  integer ( kind = 8 ) maxj
  integer ( kind = 8 ) maxk

  real ( kind = 8 ) alpha(maxgrid)
  real ( kind = 8 ) fsmach(maxgrid)
  integer ( kind = 8 ) i
  integer ( kind = 8 ) idim(maxgrid)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) igrid
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jdim(maxgrid)
  integer ( kind = 8 ) k
  integer ( kind = 8 ) kdim(maxgrid)
  integer ( kind = 8 ) l
  integer ( kind = 8 ) ngrid
  real ( kind = 8 ) q(maxi,maxj,maxk,5,maxgrid)
  real ( kind = 8 ) re(maxgrid)
  real ( kind = 8 ) time(maxgrid)

  ierror = 0

  do igrid = 1, ngrid

    read ( iunit, iostat = ios ) fsmach(igrid), alpha(igrid), re(igrid), &
      time(igrid)

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_B_3DM_Q - Fatal error!'
      write ( *, '(a)' ) '  A read error occurred.'
      return
    end if

    read ( iunit, iostat = ios ) (((( q(i,j,k,l,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), l = 1, 5 )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_B_3DM_Q - Fatal error!'
      write ( *, '(a)' ) '  A read error occurred.'
      return
    end if

  end do

  return
end
subroutine r8_b_3dm_xyz ( iunit, idim, jdim, kdim, maxi, maxj, maxk, maxgrid, &
  ngrid, x, y, z, ierror )

!*****************************************************************************80
!
!! R8_B_3DM_XYZ reads a binary 3D multiple grid XYZ file for the XYZ data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions in each grid.
!
!    Input, integer ( kind = 8 ) MAXI, MAXJ, MAXK, the maximum dimension for the first three
!    components of X, Y and Z.
!
!    Input, integer ( kind = 8 ) MAXGRID, the maximum value of NGRID.
!
!    Input, integer ( kind = 8 ) NGRID, the number of grids.
!
!    Output, real ( kind = 8 ) X(MAXI,MAXJ,MAXK,MAXGRID),
!    Y(MAXI,MAXJ,MAXK,MAXGRID),
!    Z(MAXI,MAXJ,MAXK,MAXGRID), the X, Y and Z coordinates of the nodes
!    for each grid.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) maxgrid
  integer ( kind = 8 ) maxi
  integer ( kind = 8 ) maxj
  integer ( kind = 8 ) maxk

  integer ( kind = 8 ) i
  integer ( kind = 8 ) idim(maxgrid)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) igrid
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jdim(maxgrid)
  integer ( kind = 8 ) k
  integer ( kind = 8 ) kdim(maxgrid)
  integer ( kind = 8 ) ngrid
  real ( kind = 8 ) x(maxi,maxj,maxk,maxgrid)
  real ( kind = 8 ) y(maxi,maxj,maxk,maxgrid)
  real ( kind = 8 ) z(maxi,maxj,maxk,maxgrid)

  ierror = 0

  do igrid = 1, ngrid

    read ( iunit, iostat = ios ) ((( x(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( y(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( z(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_B_3DM_XYZ - Fatal error!'
      write ( *, '(a)' ) '  A read error occurred.'
      return
    end if

  end do

  return
end
subroutine r8_b_3dm_xyzb ( iunit, idim, jdim, kdim, maxi, maxj, maxk, maxgrid, &
  ngrid, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R8_B_3DM_XYZB reads a binary 3D multiple grid XYZB file for the XYZB data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions in each grid.
!
!    Input, integer ( kind = 8 ) MAXI, MAXJ, MAXK, the maximum dimension for the first three
!    components of X, Y and Z.
!
!    Input, integer ( kind = 8 ) MAXGRID, the maximum value of NGRID.
!
!    Input, integer ( kind = 8 ) NGRID, the number of grids.
!
!    Output, real ( kind = 8 ) X(MAXI,MAXJ,MAXK,MAXGRID),
!    Y(MAXI,MAXJ,MAXK,MAXGRID),
!    Z(MAXI,MAXJ,MAXK,MAXGRID), the X, Y and Z coordinates of the nodes
!    for each grid.
!
!    Output, integer ( kind = 8 ) B(MAXI,MAXJ,MAXK,MAXGRID), the blanking information.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) maxgrid
  integer ( kind = 8 ) maxi
  integer ( kind = 8 ) maxj
  integer ( kind = 8 ) maxk

  integer ( kind = 8 ) b(maxi,maxj,maxk,maxgrid)
  integer ( kind = 8 ) i
  integer ( kind = 8 ) idim(maxgrid)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) igrid
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jdim(maxgrid)
  integer ( kind = 8 ) k
  integer ( kind = 8 ) kdim(maxgrid)
  integer ( kind = 8 ) ngrid
  real ( kind = 8 ) x(maxi,maxj,maxk,maxgrid)
  real ( kind = 8 ) y(maxi,maxj,maxk,maxgrid)
  real ( kind = 8 ) z(maxi,maxj,maxk,maxgrid)

  ierror = 0

  do igrid = 1, ngrid

    read ( iunit, iostat = ios ) ((( x(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( y(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( z(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( b(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_B_3DM_XYZB - Fatal error!'
      write ( *, '(a)' ) '  A read error occurred.'
      return
    end if

  end do

  return
end
subroutine r8_b_3dm_qfile ( iunit, idim, jdim, kdim, maxi, maxj, maxk, maxgrid, &
  ngrid, fsmach, alpha, re, time, q, ierror )

!*****************************************************************************80
!
!! R8_B_3DM_QFILE reads a binary 3D multiple grid Q file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions for each grid.
!
!    Input, integer ( kind = 8 ) MAXI, MAXJ, MAXK, the maximum dimension for the first three
!    components of Q.
!
!    Input, integer ( kind = 8 ) MAXGRID, the maximum value of NGRID.
!
!    Output, integer ( kind = 8 ) NGRID, the number of grids.
!
!    Output, real ( kind = 8 ) FSMACH(MAXGRID), ALPHA(MAXGRID), RE(MAXGRID), TIME(MAXGRID),
!    the values of the free stream Mach number, the angle of attack, in degrees,
!    the Reynolds number for the flow, and the time, for each grid.
!
!    Output, real ( kind = 8 ) Q(MAXI,MAXJ,MAXK,5,MAXGRID),
!    the Q values of the nodes for each grid.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading NGRID.
!    2, MAXGRID < NGRID.
!    3, an error occurred while reading the dimensions.
!    4, some IDIM, JDIM or KDIM is greater than MAXI, MAXJ or MAXK.
!    5, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) maxgrid
  integer ( kind = 8 ) maxi
  integer ( kind = 8 ) maxj
  integer ( kind = 8 ) maxk

  real ( kind = 8 ) alpha(maxgrid)
  real ( kind = 8 ) fsmach(maxgrid)
  integer ( kind = 8 ) i
  integer ( kind = 8 ) idim(maxgrid)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) jdim(maxgrid)
  integer ( kind = 8 ) kdim(maxgrid)
  integer ( kind = 8 ) ngrid
  real ( kind = 8 ) q(maxi,maxj,maxk,5,maxgrid)
  real ( kind = 8 ) re(maxgrid)
  real ( kind = 8 ) time(maxgrid)

  ierror = 0

  call r8_b_3dm_ngrid ( iunit, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  if ( maxgrid < ngrid ) then
    ierror = 2
    return
  end if

  call r8_b_3dm_dim ( iunit, idim, jdim, kdim, maxgrid, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  do i = 1, ngrid

    if ( maxi < idim(i) .or. maxj < jdim(i) .or. maxk < kdim(i) ) then
      ierror = 4
      return
    end if

  end do

  call r8_b_3dm_q ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
    maxgrid, ngrid, fsmach, alpha, re, time, q, ierror )

  if ( ierror /= 0 ) then
    ierror = 5
    return
  end if

  return
end
subroutine r8_b_3dm_xyzfile ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
  maxgrid, ngrid, x, y, z, ierror )

!*****************************************************************************80
!
!! R8_B_3DM_XYZFILE reads a binary 3D multiple grid XYZ file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions for each grid.
!
!    Input, integer ( kind = 8 ) MAXI, MAXJ, MAXK, the maximum dimension for the first three
!    components of X, Y and Z.
!
!    Input, integer ( kind = 8 ) MAXGRID, the maximum value of NGRID.
!
!    Output, integer ( kind = 8 ) NGRID, the number of grids.
!
!    Output, real ( kind = 8 ) X(MAXI,MAXJ,MAXK,MAXGRID),
!    Y(MAXI,MAXJ,MAXK,MAXGRID),
!    Z(MAXI,MAXJ,MAXK,MAXGRID), the X, Y and Z coordinates of the nodes
!    for each grid.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading NGRID.
!    2, MAXGRID < NGRID.
!    3, an error occurred while reading the dimensions.
!    4, some IDIM, JDIM or KDIM is greater than MAXI, MAXJ or MAXK.
!    5, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) maxgrid
  integer ( kind = 8 ) maxi
  integer ( kind = 8 ) maxj
  integer ( kind = 8 ) maxk

  integer ( kind = 8 ) i
  integer ( kind = 8 ) idim(maxgrid)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) jdim(maxgrid)
  integer ( kind = 8 ) kdim(maxgrid)
  integer ( kind = 8 ) ngrid
  real ( kind = 8 ) x(maxi,maxj,maxk,maxgrid)
  real ( kind = 8 ) y(maxi,maxj,maxk,maxgrid)
  real ( kind = 8 ) z(maxi,maxj,maxk,maxgrid)

  ierror = 0

  call r8_b_3dm_ngrid ( iunit, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  if ( maxgrid < ngrid ) then
    ierror = 2
    return
  end if

  call r8_b_3dm_dim ( iunit, idim, jdim, kdim, maxgrid, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  do i = 1, ngrid

    if ( maxi < idim(i) .or. maxj < jdim(i) .or. maxk < kdim(i) ) then
      ierror = 4
      return
    end if

  end do

  call r8_b_3dm_xyz ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
    maxgrid, ngrid, x, y, z, ierror )

  if ( ierror /= 0 ) then
    ierror = 5
    return
  end if

  return
end
subroutine r8_b_3dm_xyzbfile ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
  maxgrid, ngrid, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R8_B_3DM_XYZBFILE reads a binary 3D multiple grid XYZB file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID),
!    the number of nodes in the X, Y, and Z directions for each grid.
!
!    Input, integer ( kind = 8 ) MAXI, MAXJ, MAXK, the maximum dimension for the first three
!    components of X, Y and Z.
!
!    Input, integer ( kind = 8 ) MAXGRID, the maximum value of NGRID.
!
!    Output, integer ( kind = 8 ) NGRID, the number of grids.
!
!    Output, real ( kind = 8 ) X(MAXI,MAXJ,MAXK,MAXGRID),
!    Y(MAXI,MAXJ,MAXK,MAXGRID),
!    Z(MAXI,MAXJ,MAXK,MAXGRID), the X, Y and Z coordinates of the nodes
!    for each grid.
!
!    Output, integer ( kind = 8 ) B(MAXI,MAXJ,MAXK,MAXGRID), the blanking information.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading NGRID.
!    2, MAXGRID < NGRID.
!    3, an error occurred while reading the dimensions.
!    4, some IDIM, JDIM or KDIM is greater than MAXI, MAXJ or MAXK.
!    5, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) maxgrid
  integer ( kind = 8 ) maxi
  integer ( kind = 8 ) maxj
  integer ( kind = 8 ) maxk

  integer ( kind = 8 ) b(maxi,maxj,maxk,maxgrid)
  integer ( kind = 8 ) i
  integer ( kind = 8 ) idim(maxgrid)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) jdim(maxgrid)
  integer ( kind = 8 ) kdim(maxgrid)
  integer ( kind = 8 ) ngrid
  real ( kind = 8 ) x(maxi,maxj,maxk,maxgrid)
  real ( kind = 8 ) y(maxi,maxj,maxk,maxgrid)
  real ( kind = 8 ) z(maxi,maxj,maxk,maxgrid)

  ierror = 0

  call r8_b_3dm_ngrid ( iunit, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  if ( maxgrid < ngrid ) then
    ierror = 2
    return
  end if

  call r8_b_3dm_dim ( iunit, idim, jdim, kdim, maxgrid, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  do i = 1, ngrid

    if ( maxi < idim(i) .or. maxj < jdim(i) .or. maxk < kdim(i) ) then
      ierror = 4
      return
    end if

  end do

  call r8_b_3dm_xyzb ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
    maxgrid, ngrid, x, y, z, b, ierror )

  if ( ierror /= 0 ) then
    ierror = 5
    return
  end if

  return
end
subroutine r8_b_3dp_f ( iunit, idim, jdim, kdim, nvar, f, ierror )

!*****************************************************************************80
!
!! R8_B_3DP_F reads a binary 3D plane F file for the F data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Input, integer ( kind = 8 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 8 ) F(IDIM,JDIM,KDIM,NVAR), the F data.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim
  integer ( kind = 8 ) nvar

  real ( kind = 8 ) f(idim,jdim,kdim,nvar)
  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  integer ( kind = 8 ) l

  ierror = 0

  do k = 1, kdim

    read ( iunit, iostat = ios ) &
      ((( f(i,j,k,l), i = 1, idim ), j = 1, jdim ), l = 1, nvar )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_B_3DP_F - Error!'
      write ( *, '(a)' ) '  An end-of-file condition occurred.'
      return
    end if

  end do

  return
end
subroutine r8_b_3dp_ffile ( iunit, idim, jdim, kdim, nvar, f, ierror )

!*****************************************************************************80
!
!! R8_B_3DP_FFILE reads a binary 3D plane F file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, integer ( kind = 8 ) NVAR, the number of F variables defined.
!
!    Output, real ( kind = 8 ) F(IDIM,JDIM,KDIM,NVAR), the F data.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim
  integer ( kind = 8 ) nvar

  real ( kind = 8 ) f(idim,jdim,kdim,nvar)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit

  ierror = 0

  call r8_b_3d_dimn ( iunit, idim, jdim, kdim, nvar, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_3dp_f ( iunit, idim, jdim, kdim, nvar, f, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r8_b_3dp_q ( iunit, idim, jdim, kdim, q, ierror )

!*****************************************************************************80
!
!! R8_B_3DP_Q reads a binary 3D plane Q file for the Q data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 8 ) Q(IDIM,JDIM,KDIM,5), the Q data, density, X, Y and Z
!    momentum, and stagnation energy per unit volume.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  integer ( kind = 8 ) l
  real ( kind = 8 ) q(idim,jdim,kdim,5)

  ierror = 0

  do k = 1, kdim

    read ( iunit, iostat = ios ) &
      ((( q(i,j,k,l), i = 1, idim ), j = 1, jdim ), l = 1, 5 )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_B_3DP_Q - Error!'
      write ( *, '(a)' ) '  An end-of-file condition occurred.'
      return
    end if

  end do

  return
end
subroutine r8_b_3dp_qfile ( iunit, idim, jdim, kdim, alpha, fsmach, re, time, &
  q, ierror )

!*****************************************************************************80
!
!! R8_B_3DP_QFILE reads a binary 3D plane Q file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 8 ) ALPHA, the angle of attack, in degrees.
!
!    Output, real ( kind = 8 ) FSMACH, the free stream Mach number.
!
!    Output, real ( kind = 8 ) RE, the Reynolds number for the flow.
!
!    Output, real ( kind = 8 ) TIME, the time associated with the data.
!
!    Output, real ( kind = 8 ) Q(IDIM,JDIM,KDIM,5), the Q data, density, X, Y and Z
!    momentum, and stagnation energy per unit volume.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the parameters.
!    3, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  real ( kind = 8 ) alpha
  real ( kind = 8 ) fsmach
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  real ( kind = 8 ) q(idim,jdim,kdim,5)
  real ( kind = 8 ) re
  real ( kind = 8 ) time

  ierror = 0

  call r8_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_param ( iunit, fsmach, alpha, re, time, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  call r8_b_3dp_q ( iunit, idim, jdim, kdim, q, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  return
end
subroutine r8_b_3dp_xyz ( iunit, idim, jdim, kdim, x, y, z, ierror )

!*****************************************************************************80
!
!! R8_B_3DP_XYZ reads a binary 3D plane XYZ file for the XYZ data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 8 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM),
!    Z(IDIM,JDIM,KDIM), the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  real ( kind = 8 ) x(idim,jdim,kdim)
  real ( kind = 8 ) y(idim,jdim,kdim)
  real ( kind = 8 ) z(idim,jdim,kdim)

  ierror = 0

  do k = 1, kdim

    read ( iunit, iostat = ios ) (( x(i,j,k), i = 1, idim ), j = 1, jdim ), &
      (( y(i,j,k), i = 1, idim ), j = 1, jdim ), &
      (( z(i,j,k), i = 1, idim ), j = 1, jdim )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_B_3DP_XYZ - Error!'
      write ( *, '(a)' ) '  An end-of-file condition occurred.'
      return
    end if

  end do

  return
end
subroutine r8_b_3dp_xyzb ( iunit, idim, jdim, kdim, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R8_B_3DP_XYZB reads a binary 3D plane XYZB file for the XYZB data.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 8 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM), Z(IDIM,JDIM,KDIM),
!    the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 8 ) B(IDIM,JDIM,KDIM), the blanking information.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  integer ( kind = 8 ) b(idim,jdim,kdim)
  integer ( kind = 8 ) i
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  real ( kind = 8 ) x(idim,jdim,kdim)
  real ( kind = 8 ) y(idim,jdim,kdim)
  real ( kind = 8 ) z(idim,jdim,kdim)

  ierror = 0

  do k = 1, kdim

    read ( iunit, iostat = ios ) (( x(i,j,k), i = 1, idim ), j = 1, jdim ), &
      (( y(i,j,k), i = 1, idim ), j = 1, jdim ), &
      (( z(i,j,k), i = 1, idim ), j = 1, jdim ), &
      (( b(i,j,k), i = 1, idim ), j = 1, jdim )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_B_3DP_XYZB - Error!'
      write ( *, '(a)' ) '  An end-of-file condition occurred.'
    end if

  end do

  return
end
subroutine r8_b_3dp_xyzbfile ( iunit, idim, jdim, kdim, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R8_B_3DP_XYZBFILE reads a binary 3D plane XYZB file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 8 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM), Z(IDIM,JDIM,KDIM),
!    the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 8 ) B(IDIM,JDIM,KDIM), the blanking information.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  integer ( kind = 8 ) b(idim,jdim,kdim)
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  real ( kind = 8 ) x(idim,jdim,kdim)
  real ( kind = 8 ) y(idim,jdim,kdim)
  real ( kind = 8 ) z(idim,jdim,kdim)

  ierror = 0

  call r8_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_3dp_xyzb ( iunit, idim, jdim, kdim, x, y, z, b, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r8_b_3dp_xyzfile ( iunit, idim, jdim, kdim, x, y, z, ierror )

!*****************************************************************************80
!
!! R8_B_3DP_XYZFILE reads a binary 3D plane XYZ file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 8 ) IDIM, JDIM, KDIM, the number of nodes in the
!    X, Y, and Z directions.
!
!    Output, real ( kind = 8 ) X(IDIM,JDIM,KDIM), Y(IDIM,JDIM,KDIM), Z(IDIM,JDIM,KDIM),
!    the X, Y and Z coordinates of the nodes.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the dimensions.
!    2, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 8 ) idim
  integer ( kind = 8 ) jdim
  integer ( kind = 8 ) kdim

  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) iunit
  real ( kind = 8 ) x(idim,jdim,kdim)
  real ( kind = 8 ) y(idim,jdim,kdim)
  real ( kind = 8 ) z(idim,jdim,kdim)

  ierror = 0

  call r8_b_3d_dim ( iunit, idim, jdim, kdim, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  call r8_b_3dp_xyz ( iunit, idim, jdim, kdim, x, y, z, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    return
  end if

  return
end
subroutine r8_b_param ( iunit, fsmach, alpha, re, time, ierror )

!*****************************************************************************80
!
!! R8_B_PARAM reads a binary Q file for the flow parameters.
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
!    Input, integer ( kind = 8 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, real ( kind = 8 ) FSMACH, the free stream Mach number.
!
!    Output, real ( kind = 8 ) ALPHA, the angle of attack, in degrees.
!
!    Output, real ( kind = 8 ) RE, the Reynolds number for the flow.
!
!    Output, real ( kind = 8 ) TIME, the time associated with the data.
!
!    Output, integer ( kind = 8 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) fsmach
  integer ( kind = 8 ) ierror
  integer ( kind = 8 ) ios
  integer ( kind = 8 ) iunit
  real ( kind = 8 ) re
  real ( kind = 8 ) time

  ierror = 0

  read ( iunit, iostat = ios ) fsmach, alpha, re, time

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B_PARAM - Error!'
    write ( *, '(a)' ) '  An end-of-file condition occurred.'
    return
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
