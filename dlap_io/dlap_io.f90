subroutine dlap_file_print ( n, nelt, isym, irhs, isoln, ia, ja, a, rhs, &
  soln, title )

!*****************************************************************************80
!
!! DLAP_FILE_PRINT prints a DLAP linear system that was stored in a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NELT, the number of non-zeros stored in A.
!
!    Input, integer ( kind = 4 ) ISYM, a flag to indicate symmetric storage format.
!    * 0, all nonzero entries of the matrix are stored.
!    * 1, the matrix is symmetric, and only the lower triangle of the
!    matrix is stored.
!
!    Input, integer ( kind = 4 ) IRHS, is 1 if a right hand side vector is included.
!
!    Input, integer ( kind = 4 ) ISOLN, is 1 if a solution vector is included.
!
!    Input, integer ( kind = 4 ) IA(NELT), integer ( kind = 4 ) JA(NELT),
!    real ( kind = 8 ) A(NELT), the DLAP triad matrix description.
!
!    Input, real ( kind = 8 ) RHS(N), the right hand side vector.
!
!    Input, real ( kind = 8 ) SOLN(N), the solution to the linear system.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt

  real ( kind = 8 ) a(nelt)
  integer ( kind = 4 ) ia(nelt)
  integer ( kind = 4 ) irhs
  integer ( kind = 4 ) isoln
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) ja(nelt)
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) soln(n)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  call dlap_header_print ( n, nelt, isym, irhs, isoln )
!
!  Write out the matrix.
!
  call ds3_print ( n, n, nelt, isym, ia, ja, a, '  The sparse matrix' )
!
!  Write the right hand side.
!
  if ( irhs == 1 ) then
    call dlap_rhs_print ( n, rhs )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  No right hand side vector was supplied.'
  end if
!
!  Write the solution.
!
  if ( isoln == 1 ) then
    call dlap_soln_print ( n, soln )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  No solution vector was supplied.'
  end if

  return
end
subroutine dlap_file_read ( n_max, nelt_max, n, nelt, isym, irhs, isoln, &
  ia, ja, a, rhs, soln, iunit )

!*****************************************************************************80
!
!! DLAP_FILE_READ reads in a DLAP matrix contained in a file.
!
!  Discussion:
!
!    This routine reads in a DLAP Triad Format Linear System,
!    including the matrix, right hand side, and solution, if known.
!
!
!    The original version of this program seems to have a minor
!    logical flaw.  If the user requests the solution but not
!    the right hand side, and the file contains both, the original
!    program would not correctly read past the right hand side
!    to get to the solution.  The current version should fix
!    that flaw.
!
!
!    The expected format of the file is as follows.  On the first line
!    are counters and flags: N, NELT, ISYM, IRHS, ISOLN.  N, NELT
!    and ISYM are described below.  IRHS is a flag indicating if
!    the RHS was written out (1 is yes, 0 is  no).  ISOLN is a
!    flag indicating if the SOLN was written out  (1 is yes, 0 is
!    no).  The format for the first line is: 5i10.  Then comes the
!    NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format
!    for these lines is   :  1X,I5,1X,I5,1X,E16.7.   Then comes
!    RHS(I), I = 1, N, if IRHS = 1.  Then comes SOLN(I), I  = 1,
!    N, if ISOLN = 1.  The format for these lines is: 1X,E16.7.
!
!
!    This routine requires that the  matrix A be stored in the
!    DLAP Triad format.  In this format only the non-zeros  are
!    stored.  They may appear in ANY order.  The user supplies
!    three arrays of length NELT, where NELT is the number of
!    non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
!    each non-zero the user puts the row and column index of that
!    matrix element in the IA and JA arrays.  The value of the
!    non-zero matrix element is placed in the corresponding
!    location of the A array.   This is an extremely easy data
!    structure to generate.  On the other hand it is not too
!    efficient on vector computers for the iterative solution of
!    linear systems.  Hence,  DLAP changes this input data
!    structure to the DLAP Column format for the iteration (but
!    does not change it back).
!
!    Here is an example of the DLAP Triad storage format for a
!    5x5 Matrix.  Recall that the entries may appear in any order.
!
!        5x5 Matrix       DLAP Triad format for 5x5 matrix on left.
!                              1  2  3  4  5  6  7  8  9 10 11
!    |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
!    |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
!    | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
!    | 0  0  0 44  0|
!    |51  0 53  0 55|
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_MAX, the maximum value of N for which storage
!    has been allocated.
!
!    Input, integer ( kind = 4 ) NELT_MAX, the maximum value of NELT for which
!    storage has been allocated.
!
!    Output, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) NELT, the number of non-zeros stored in A.
!
!    Output, integer ( kind = 4 ) ISYM, a flag to indicate symmetric storage 
!    format.
!    * 0, all nonzero entries of the matrix are stored.
!    * 1, the matrix is symmetric, and only the lower triangle of the
!    matrix is stored.
!
!    Output, integer ( kind = 4 ) IRHS, is 1 if a right hand side vector is included.
!
!    Output, integer ( kind = 4 ) ISOLN, is 1 if a solution vector is included.
!
!    Output, integer ( kind = 4 ) IA(NELT), integer ( kind = 4 ) JA(NELT),
!    real ( kind = 8 ) A(NELT).  On output these arrays hold the matrix A in the
!    DLAP Triad format.
!
!    Output, real ( kind = 8 ) RHS(N), the right hand side vector.
!
!    Output, real ( kind = 8 ) SOLN(N), the solution to the linear system, 
!    if present.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN device unit from which the
!    matrix is to be read.
!
  implicit none

  integer ( kind = 4 ) n_max
  integer ( kind = 4 ) nelt_max

  real ( kind = 8 ) a(nelt_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nelt_max)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) irhs
  integer ( kind = 4 ) isoln
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ja(nelt_max)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt
  real ( kind = 8 ) rhs(n_max)
  real ( kind = 8 ) soln(n_max)
!
!  Read the header line.
!
  call dlap_header_read ( iunit, n, nelt, isym, irhs, isoln, ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SLAP_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  Error while reading header line of DLAP file.'
    return
  end if

  if ( nelt_max < nelt ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SLAP_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  NELT_MAX < NELT.'
    write ( *, '(a,i8)' ) '  NELT_MAX = ', nelt_max
    write ( *, '(a,i8)' ) '  NELT     = ', nelt
    stop
  end if
!
!  Read the nonzero matrix entries in DLAP Triad format.
!
  do i = 1, nelt

    read ( iunit, '(1x,i5,1x,i5,1x,e16.7)', iostat = ios ) ia(i), ja(i), a(i)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SLAP_FILE_READ - Fatal error!'
      write ( *, '(a,i8)' ) '  Error while reading matrix element ', i+1
      return
    end if

  end do
!
!  If a value for RHS is available in the file, read it in.
!
  if ( irhs == 1 ) then

    if ( n_max < n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SLAP_FILE_READ - Fatal error!'
      write ( *, '(a)' ) '  N_MAX < N.'
      write ( *, '(a,i8)' ) '  N_MAX = ', n_max
      write ( *, '(a,i8)' ) '  N     = ', n
      stop
    end if

    call dlap_rhs_read ( iunit, n, rhs, ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SLAP_FILE_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading RHS from DLAP file.'
      return
    end if

  end if
!
!  If a value of SOLN is available in the file, read it.
!
  if ( isoln == 1 ) then

    if ( n_max < n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SLAP_FILE_READ - Fatal error!'
      write ( *, '(a)' ) '  N_MAX < N.'
      write ( *, '(a,i8)' ) '  N_MAX = ', n_max
      write ( *, '(a,i8)' ) '  N     = ', n
      stop
    end if

    call dlap_soln_read ( iunit, n, soln, ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SLAP_FILE_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading SOLN from DLAP file.'
      return
    end if

  end if

  return
end
subroutine dlap_file_write ( n, nelt, isym, irhs, isoln, ia, ja, a, rhs, &
  soln, iunit )

!*****************************************************************************80
!
!! DLAP_FILE_WRITE writes out DLAP Triad Format Linear System.
!
!  Discussion:
!
!    This routine writes out a DLAP Triad Format Linear System.
!    including the matrix, right hand side, and solution to the
!    system, if known.
!
!
!    The format for the output is as follows.  On  the first line
!    are counters and flags:
!
!      N, NELT, ISYM, IRHS, ISOLN.
!
!    N, NELT and ISYM are described below.  IRHS is a flag indicating if
!    the RHS was written out (1 is  yes, 0 is  no).  ISOLN  is a
!    flag indicating if the SOLN was written out  (1 is yes, 0 is
!    no).  The format for the first line is: 5i10.  Then comes the
!    NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format
!    for  these lines is   :  1X,I5,1X,I5,1X,E16.7.   Then comes
!    RHS(I), I = 1, N, if IRHS = 1.  Then comes SOLN(I), I  = 1,
!    N, if ISOLN = 1.  The format for these lines is: 1X,E16.7.
!
!
!    This routine requires that the  matrix A be stored in the
!    DLAP Triad format.  In this format only the non-zeros  are
!    stored.  They may appear in ANY order.  The user supplies
!    three arrays of length NELT, where NELT is the number of
!    non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
!    each non-zero the user puts the row and column index of that
!    matrix element in the IA and JA arrays.  The value of the
!    non-zero matrix element is placed in the corresponding
!    location of the A array.   This is an extremely easy data
!    structure to generate.  On the other hand it is not too
!    efficient on vector computers for the iterative solution of
!    linear systems.  Hence,  DLAP changes this input data
!    structure to the DLAP Column format for the iteration (but
!    does not change it back).
!
!    Here is an example of the DLAP Triad storage format for a
!    5x5 Matrix.  Recall that the entries may appear in any order.
!
!        5x5 Matrix       DLAP Triad format for 5x5 matrix on left.
!                              1  2  3  4  5  6  7  8  9 10 11
!    |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
!    |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
!    | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
!    | 0  0  0 44  0|
!    |51  0 53  0 55|
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NELT, the number of non-zeros stored in A.
!
!    Input, integer ( kind = 4 ) ISYM, indicates symmetric storage format.
!    * 0, all nonzero entries of the matrix are stored.
!    * 1, the matrix is symmetric, and only the lower triangle of
!      the matrix is stored.
!
!    Input, integer ( kind = 4 ) IRHS, is 1 if a right hand side vector is included.
!
!    Input, integer ( kind = 4 ) ISOLN, is 1 if a solution vector is included.
!
!    Input, integer ( kind = 4 ) IA(NELT), integer ( kind = 4 ) JA(NELT),
!    real ( kind = 8 ) A(NELT), the DLAP triad matrix description.
!
!    Input, real ( kind = 8 ) RHS(N), the right hand side vector.  This array is
!    accessed if JOB is set to print it out.
!
!    Input, real ( kind = 8 ) SOLN(N), the solution to the linear system, if known.
!    This array is accessed if and only if JOB is set to print it out.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt

  real ( kind = 8 ) a(nelt)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nelt)
  integer ( kind = 4 ) irhs
  integer ( kind = 4 ) isoln
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ja(nelt)
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) soln(n)

  call dlap_header_write ( iunit, n, nelt, isym, irhs, isoln )
!
!  Write the matrix non-zeros in Triad format.
!
  do i = 1, nelt
    write ( iunit, '(1x,i5,1x,i5,1x,e16.7)' ) ia(i), ja(i), a(i)
  end do
!
!  Write the right hand side.
!
  if ( irhs == 1 ) then
    call dlap_rhs_write ( iunit, n, rhs )
  end if
!
!  Write the solution.
!
  if ( isoln == 1 ) then
    call dlap_soln_write ( iunit, n, soln )
  end if

  return
end
subroutine dlap_header_print ( n, nelt, isym, irhs, isoln )

!*****************************************************************************80
!
!! DLAP_HEADER_PRINT prints the header line of a DLAP file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NELT, the number of non-zeros stored in A.
!
!    Input, integer ( kind = 4 ) ISYM, indicates symmetric storage format.
!    * 0, all nonzero entries of the matrix are stored.
!    * 1, the matrix is symmetric, and only the lower triangle of
!      the matrix is stored.
!
!    Input, integer ( kind = 4 ) IRHS, is 1 if a right hand side vector 
!    is included.
!
!    Input, integer ( kind = 4 ) ISOLN, is 1 if a solution vector is included.
!
  implicit none

  integer ( kind = 4 ) isoln
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) irhs
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SLAP Sparse Matrix File Header:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N,     the matrix order =                ', n
  write ( *, '(a,i8)' ) '  NELT,  the number of nonzeros stored =   ', nelt
  write ( *, '(a,i8)' ) '  ISYM,  1 if symmetric storage used =     ', isym
  write ( *, '(a,i8)' ) '  IRHS,  1 if a right hand side included = ', irhs
  write ( *, '(a,i8)' ) '  ISOLN, 1 if a solution vector included = ', isoln

  return
end
subroutine dlap_header_read ( iunit, n, nelt, isym, irhs, isoln, ios )

!*****************************************************************************80
!
!! DLAP_HEADER_READ reads the header line from a DLAP file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
!    Output, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) NELT, the number of non-zeros stored in A.
!
!    Output, integer ( kind = 4 ) ISYM, indicates symmetric storage format.
!    * 0, all nonzero entries of the matrix are stored.
!    * 1, the matrix is symmetric, and only the lower triangle of
!      the matrix is stored.
!
!    Output, integer ( kind = 4 ) IRHS, is 1 if a right hand side vector
!    is included.
!
!    Output, integer ( kind = 4 ) ISOLN, is 1 if a solution vector is included.
!
!    Output, integer ( kind = 4 ) IOS, the I/O status variable, which is 0 if
!    no I/O error occurred.
!
  implicit none

  integer ( kind = 4 ) ios
  integer ( kind = 4 ) isoln
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) irhs
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt

  read ( iunit, *, iostat = ios ) n, nelt, isym, irhs, isoln

  return
end
subroutine dlap_header_write ( iunit, n, nelt, isym, irhs, isoln )

!*****************************************************************************80
!
!! DLAP_HEADER_WRITE writes the header line to a DLAP file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NELT, the number of non-zeros stored in A.
!
!    Input, integer ( kind = 4 ) ISYM, indicates symmetric storage format.
!    * 0, all nonzero entries of the matrix are stored.
!    * 1, the matrix is symmetric, and only the lower triangle of
!      the matrix is stored.
!
!    Input, integer ( kind = 4 ) IRHS, is 1 if a right hand side vector is included.
!
!    Input, integer ( kind = 4 ) ISOLN, is 1 if a solution vector is included.
!
  implicit none

  integer ( kind = 4 ) isoln
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) irhs
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt

  write ( iunit, '(5i10)' ) n, nelt, isym, irhs, isoln

  return
end
subroutine dlap_rhs_print ( n, rhs )

!*****************************************************************************80
!
!! DLAP_RHS_PRINT prints the right hand side vector from a DLAP file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) RHS(N), the right hand side vector to be written.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) rhs(n)

  call r8vec_print ( n, rhs, '  DLAP right hand side vector:' )

  return
end
subroutine dlap_rhs_read ( iunit, n, rhs, ios )

!*****************************************************************************80
!
!! DLAP_RHS_READ reads the right hand side vector from a DLAP file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) RHS(N), the right hand side vector.
!
!    Output, integer ( kind = 4 ) IOS, the I/O status variable, which is 0 if
!    no I/O error occurred.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) rhs(n)

  read ( iunit, '(1x,e16.7)', iostat = ios ) rhs(1:n)

  return
end
subroutine dlap_rhs_write ( iunit, n, rhs )

!*****************************************************************************80
!
!! DLAP_RHS_WRITE writes a right hand side vector to a DLAP file.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) RHS(N), the right hand side vector to be written.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iunit
  real ( kind = 8 ) rhs(n)

  write ( iunit, '(1x,e16.7)' ) rhs(1:n)

  return
end
subroutine dlap_soln_print ( n, soln )

!*****************************************************************************80
!
!! DLAP_SOLN_PRINT prints the solution vector from a DLAP file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) SOLN(N), the solution vector to be written.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) soln(n)

  call r8vec_print ( n, soln, '  DLAP solution vector:' )

  return
end
subroutine dlap_soln_read ( iunit, n, soln, ios )

!*****************************************************************************80
!
!! DLAP_SOLN_READ reads the solution vector from a DLAP file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) SOLN(N), the solution vector.
!
!    Output, integer ( kind = 4 ) IOS, the I/O status variable, which is 0 if
!    no I/O error occurred.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) soln(n)

  read ( iunit, '(1x,e16.7)', iostat = ios ) soln(1:n)

  return
end
subroutine dlap_soln_write ( iunit, n, soln )

!*****************************************************************************80
!
!! DLAP_SOLN_WRITE writes a solution vector to a DLAP file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) SOLN(N), the solution vector to be written.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iunit
  real ( kind = 8 ) soln(n)

  write ( iunit, '(2x,e16.7)' ) soln(1:n)

  return
end
subroutine ds3_print ( m, n, nz_num, isym, row, col, a, title )

!*****************************************************************************80
!
!! DS3_PRINT prints a DS3 matrix.
!
!  Discussion:
!
!    The DS3 storage format corresponds to the DLAP Triad format.
!
!    The DS3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ( kind = 4 ) ISYM, is 0 if the matrix is not symmetric, and 1
!    if the matrix is symmetric.  The symmetric case only makes sense
!    if the matrix is also square, that is, M = N.  In this case, only
!    the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  call ds3_print_some ( m, n, nz_num, isym, row, col, a, 1, 1, m, n )

  return
end
subroutine ds3_print_some ( m, n, nz_num, isym, row, col, a, ilo, jlo, ihi, &
  jhi )

!*****************************************************************************80
!
!! DS3_PRINT_SOME prints some of a DS3 matrix.
!
!  Discussion:
!
!    The DS3 storage format corresponds to the DLAP Triad format.
!
!    The DS3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ( kind = 4 ) ISYM, is 0 if the matrix is not symmetric, and 1
!    if the matrix is symmetric.  The symmetric case only makes sense
!    if the matrix is also square, that is, M = N.  In this case, only
!    the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column 
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) aij
  integer ( kind = 4 ) col(nz_num)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  logical nonzero
  integer ( kind = 4 ) row(nz_num)
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      nonzero = .false.

      aij = 0.0D+00
      do j2 = 1, inc
        write ( ctemp(j2), '(f8.0,6x)' ) aij
      end do

      do k = 1, nz_num

        if ( i == row(k) .and. j2lo <= col(k) .and. col(k) <= j2hi ) then

          j2 = col(k) - j2lo + 1
          aij = a(k)

          if ( aij == 0.0D+00 ) then
            cycle
          end if

          nonzero = .true.

          write ( ctemp(j2), '(g14.6)' ) aij


        else if ( isym == 1 .and. m == n .and. &
          i == col(k) .and. j2lo <= row(k) .and. row(k) <= j2hi ) then

          j2 = row(k) - j2lo + 1
          aij = a(k)

          if ( aij == 0.0D+00 ) then
            cycle
          end if

          nonzero = .true.

          write ( ctemp(j2), '(g14.6)' ) aij

        end if

      end do

      if ( nonzero ) then
        write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )
      end if

    end do

  end do

  write ( *, '(a)' ) ' '

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
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,g16.8)' ) i, a(i)
  end do

  return
end
subroutine r8vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print-2
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
    end do
    write ( *, '(a)' ) '  ......  ..............'
    i = n
    write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,2x,g14.6,2x,a)' ) i, a(i), '...more entries...'

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
!    26 February 2005
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
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
