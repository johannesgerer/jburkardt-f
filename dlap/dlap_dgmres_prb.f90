program main

!*****************************************************************************80
!
!! DLAP_DGMRES_PRB tests the DLAP GMRES solver.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ligw = 20
  integer ( kind = 4 ), parameter :: maxl = 20
  integer ( kind = 4 ), parameter :: n = 100

  integer ( kind = 4 ), parameter :: lrgw = 1 + n * ( maxl + 6 ) + maxl * ( maxl + 3 )

  real ( kind = 8 ), allocatable, dimension ( : ) :: a
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) err
  integer ( kind = 4 ), allocatable, dimension ( : ) :: ia
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) igwk(ligw)
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) itmax
  integer ( kind = 4 ) itol
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) iwork(1)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: ja
  external matvec_triad
  integer ( kind = 4 ) mode
  external msolve_identity
  integer ( kind = 4 ) nelt
  real ( kind = 8 ) rgwk(lrgw)
  real ( kind = 8 ) rwork(1)
  real ( kind = 8 ) sb(n)
  real ( kind = 8 ) sx(n)
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DLAP_DGMRES_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the DLAP GMRES solver.'

  isym = 0
  itol = 0
  tol = 0.00001D+00
  itmax = 500
  iunit = 0
  sb(1:n) = 1.0D+00
  sx(1:n) = 1.0D+00

  igwk(1) = maxl
  igwk(2) = maxl
  igwk(3) = 0
  igwk(4) = 0
  igwk(5) = 60
!
!  Determine the amount of storage needed for the matrix.
!
  mode = 0
  call matset_triad ( n, mode, nelt, ia, ja, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Problem dimension N = ', n
  write ( *, '(a,i6)' ) '  Matrix storage NELT = ', nelt

  allocate ( ia(1:nelt) )
  allocate ( ja(1:nelt) )
  allocate (  a(1:nelt) )
!
!  Assign the values for the matrix.
!
  mode = 1
  call matset_triad ( n, mode, nelt, ia, ja, a )
!
!  Set the solution X.
!
  call r8vec_indicator ( n, x )

  call r8vec_print_some ( n, x, 10, '  Part of the exact solution X:' )
!
!  Set the right hand side B!
!
  call matvec_triad ( n, x, b, nelt, ia, ja, a, isym )

  call r8vec_print_some ( n, b, 10, '  Part of the right hand side B:' )
!
!  Zero out X, because it will be used as an initial guess.
!
  x(1:n) = 0.0D+00
!
!  Call DGMRES to solve the system..
!
  call dgmres ( n, b, x, nelt, ia, ja, a, isym, matvec_triad, &
    msolve_identity, itol, tol, itmax, iter, err, ierr, iunit, sb, &
    sx, rgwk, lrgw, igwk, ligw, rwork, iwork )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of iterations:  ', iter
  write ( *, '(a,g14.6)' ) '  Convergence measure is ', rgwk(1)
  write ( *, '(a,g14.6)' ) '  Error estimate ', err
  write ( *, '(a,i6)' ) '  Error code is ', ierr

  call r8vec_print_some ( n, x, 10, '  Part of the computed solution X:' )
!
!  Free memory.
!
  deallocate ( ia )
  deallocate ( ja )
  deallocate ( a )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DLAP_DGMRES_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine matset_triad ( n, mode, nelt, ia, ja, a )

!*****************************************************************************80
!
!! MATSET_TRIAD sets the data structure for a sparse matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vectors.
!
!    Input, integer MODE.
!    0, setup mode.  Count the number of nonzero elements.
!    1, assignment mode.  NELT is set.  Set the matrix entries.
!
!    Input/output, integer NELT, the number of nonzero entries.
!    On call with MODE = 0, this is an output quantity.
!    On call with MODE = 1, this is an input quantity.
!
!    Output, integer IA(NELT), JA(NELT), real ( kind = 8 ) A(NELT), the data
!    structure storing the sparse matrix.  These values are only
!    set on a call with MODE = 1.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt

  real ( kind = 8 ) a(nelt)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nelt)
  integer ( kind = 4 ) ja(nelt)
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) k

  k = 0
  do i = 1, n

    if ( 1 < i ) then
      k = k + 1
      if ( mode == 1 ) then
        ia(k) = i
        ja(k) = i-1
        a(k) = -1.0D+00
      end if
    end if

    k = k + 1
    if ( mode == 1 ) then
      ia(k) = i
      ja(k) = i
      a(k) = 2.0D+00
    end if

    if ( i < n ) then
      k = k + 1
      if ( mode == 1 ) then
        ia(k) = i
        ja(k) = i+1
        a(k) = -1.0D+00
      end if
    end if

  end do

  if ( mode == 0 ) then
    nelt = k
  end if

  return
end
subroutine matvec_triad ( n, x, y, nelt, ia, ja, a, isym )

!*****************************************************************************80
!
!! MATVEC_TRIAD computes A*X for a matrix A stored in SLAP Triad form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vectors.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) Y(N), the product A * X.
!
!    Input, integer NELT, the number of nonzero entries in A.
!
!    Input, integer IA(NELT), JA(NELT), real ( kind = 8 ) A(NELT), the data
!    structure storing the sparse matrix.
!
!    Input, integer ISYM, is 0 if all nonzero entries of the matrix
!    are stored, and 1 if only the diagonal and upper or lower triangle
!    are stored.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt

  real ( kind = 8 ) a(nelt)
  integer ( kind = 4 ) ia(nelt)
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) ja(nelt)
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  y(1:n) = 0.0D+00

  do k = 1, nelt
    y(ia(k)) = y(ia(k)) + a(k) * x(ja(k))
  end do

  return
end
subroutine msolve_identity ( n, r, z, nelt, ia, ja, a, isym, rwork, iwork )

!*****************************************************************************80
!
!! MSOLVE_IDENTITY applies the identity matrix preconditioner.
!
!  Discussion:
!
!    Most SLAP solver routines require a preconditioner routine
!    that can solve M * Z = R.  If no preconditioning is required,
!    then you can simply behave as though the preconditioning matrix
!    M was the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the vectors.
!
!    Input, real ( kind = 8 ) R(N), the right hand side.
!
!    Output, real ( kind = 8 ) Z(N), the solution of M * Z = R.
!
!    Input, integer ( kind = 4 ) NELT, the number of nonzero entries in A.
!
!    Input, integer ( kind = 4 ) IA(NELT), JA(NELT), real ( kind = 8 ) A(NELT), 
!    the data structure storing the sparse matrix.
!
!    Input, integer ( kind = 4 ) ISYM, is 0 if all nonzero entries of the matrix
!    are stored, and 1 if only the diagonal and upper or lower triangle
!    are stored.
!
!    Input, real ( kind = 8 ) RWORK(*), a real array that
!    can be used to pass information to the preconditioner.
!
!    Input, integer ( kind = 4 ) IWORK(*), an integer array that
!    can be used to pass information to the preconditioner.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt

  real ( kind = 8 ) a(nelt)
  integer ( kind = 4 ) ia(nelt)
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) iwork(*)
  integer ( kind = 4 ) ja(nelt)
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rwork(*)
  real ( kind = 8 ) z(n)

  z(1:n) = r(1:n)

  return
end
subroutine r8vec_indicator ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR sets a real vector to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, real ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = real ( i, kind = 8 )
  end do

  return
end
subroutine r8vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of a real vector.
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
!    16 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
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

    if ( all ( a(1:n) == aint ( a(1:n) ) ) ) then
      do i = 1, n
        write ( *, '(i6,2x,i6)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:n) ) < 1000000.0D+00 ) ) then
      do i = 1, n
        write ( *, '(i6,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, n
        write ( *, '(i6,2x,g14.6)' ) i, a(i)
      end do
    end if

  else if ( 3 <= max_print ) then

    if ( all ( a(1:max_print-2) == aint ( a(1:max_print-2) ) ) ) then
      do i = 1, max_print-2
        write ( *, '(i6,2x,i6)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-2) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print-2
        write ( *, '(i6,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print-2
        write ( *, '(i6,2x,g14.6)' ) i, a(i)
      end do
    end if

    write ( *, '(a)' ) '......  ..............'
    i = n

    if ( a(i) == aint ( a(i) ) ) then
      write ( *, '(i6,2x,i6)' ) i, int ( a(i) )
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i6,2x,f14.6)' ) i, a(i)
    else
      write ( *, '(i6,2x,g14.6)' ) i, a(i)
    end if

  else

    if ( all ( a(1:max_print-1) == aint ( a(1:max_print-1) ) ) ) then
      do i = 1, max_print-1
        write ( *, '(i6,2x,i6)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-1) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print-1
        write ( *, '(i6,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print-1
        write ( *, '(i6,2x,g14.6)' ) i, a(i)
      end do
    end if

    i = max_print

    if ( a(i) == aint ( a(i) ) ) then
      write ( *, '(i6,2x,i6,a)' ) i, int ( a(i) ), '...more entries...'
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i6,2x,f14.6,a)' ) i, a(i), '...more entries...'
    else
      write ( *, '(i6,2x,g14.6,a)' ) i, a(i), '...more entries...'
    end if

  end if

  return
end
