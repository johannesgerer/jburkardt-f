subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine toep_cholesky_lower ( n, g, l )

!*****************************************************************************80
!
!! TOEP_CHOLESKY_LOWER: lower Cholesky factor of a Toeplitz matrix.
!
!  Discussion:
!
!    The Toeplitz matrix A is supplied in a compressed form G.
!
!    The Toeplitz matrix must be positive semi-definite.
!
!    After factorization, A = L * L'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Michael Stewart,
!    Cholesky factorization of semi-definite Toeplitz matrices.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) G(2,N), the compressed Toeplitz matrix.
!    G(1,1:N) contains the first row.
!    G(2,2:N) contains the first column.
!
!    Output, real ( kind = 8 ) L(N,N), the lower Cholesky factor.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(2,n)
  real ( kind = 8 ) h(2,2)
  integer ( kind = 4 ) i
  real ( kind = 8 ) l(n,n)
  real ( kind = 8 ) rho

  l(1:n,1:n) = 0.0D+00

  l(1:n,1) = g(1,1:n)
  g(1,2:n) = g(1,1:n-1)
  g(1,1) = 0.0D+00
  do i = 2, n
    rho = - g(2,i) / g(1,i)
    h = reshape ( (/ 1.0D+00, rho, rho, 1.0D+00 /), (/ 2, 2 /) )
    g(1:2,i:n) = matmul ( h, g(1:2,i:n) ) &
      / sqrt ( ( 1.0D+00 - rho ) * ( 1.0D+00 + rho ) )
    l(i:n,i) = g(1,i:n)
    g(1,i+1:n) = g(1,i:n-1)
    g(1,i) = 0.0D+00
  end do

  return
end
subroutine toep_cholesky_upper ( n, g, r )

!*****************************************************************************80
!
!! TOEP_CHOLESKY_UPPER: upper Cholesky factor of a Toeplitz matrix.
!
!  Discussion:
!
!    The Toeplitz matrix A is supplied in a compressed form G.
!
!    The Toeplitz matrix must be positive semi-definite.
!
!    After factorization, A = R' * R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Michael Stewart,
!    Cholesky factorization of semi-definite Toeplitz matrices.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) G(2,N), the compressed Toeplitz matrix.
!    G(1,1:N) contains the first row.
!    G(2,2:N) contains the first column.
!
!    Output, real ( kind = 8 ) R(N,N), the upper Cholesky factor.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) g(2,n)
  real ( kind = 8 ) h(2,2)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) rho

  r(1:n,1:n) = 0.0D+00

  r(1,1:n) = g(1,1:n)
  g(1,2:n) = g(1,1:n-1)
  g(1,1) = 0.0D+00

  do i = 2, n
    rho = - g(2,i) / g(1,i)
    h = reshape ( (/ 1.0D+00, rho, rho, 1.0D+00 /), (/ 2, 2 /) )
    g(1:2,i:n) = matmul ( h, g(1:2,i:n) ) &
      / sqrt ( ( 1.0D+00 - rho ) * ( 1.0D+00 + rho ) )
    r(i,i:n) = g(1,i:n)
    g(1,i+1:n) = g(1,i:n-1)
    g(1,i) = 0.0D+00
  end do

  return
end
subroutine toeplitz_cholesky_lower ( n, a, l )

!*****************************************************************************80
!
!! TOEPLITZ_CHOLESKY_LOWER: lower Cholesky factor of a Toeplitz matrix.
!
!  Discussion:
!
!    The Toeplitz matrix must be positive semi-definite.
!
!    After factorization, A = L * L'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Michael Stewart,
!    Cholesky factorization of semi-definite Toeplitz matrices.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the Toeplitz matrix.
!
!    Output, real ( kind = 8 ) L(N,N), the lower Cholesky factor.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) g(2,n)
  real ( kind = 8 ) h(2,2)
  integer ( kind = 4 ) i
  real ( kind = 8 ) l(n,n)
  real ( kind = 8 ) rho

  l(1:n,1:n) = 0.0D+00

  g(1,1:n) = a(1,1:n)
  g(2,1) = 0.0D+00
  g(2,2:n) = a(2:n,1) 

  l(1:n,1) = g(1,1:n)
  g(1,2:n) = g(1,1:n-1)
  g(1,1) = 0.0D+00
  do i = 2, n
    rho = - g(2,i) / g(1,i)
    h = reshape ( (/ 1.0D+00, rho, rho, 1.0D+00 /), (/ 2, 2 /) )
    g(1:2,i:n) = matmul ( h, g(1:2,i:n) ) &
      / sqrt ( ( 1 - rho ) * ( 1 + rho ) )
    l(i:n,i) = g(1,i:n)
    g(1,i+1:n) = g(1,i:n-1)
    g(1,i) = 0.0D+00
  end do

  return
end
subroutine toeplitz_cholesky_upper ( n, a, r )

!*****************************************************************************80
!
!! TOEPLITZ_CHOLESKY_UPPER: upper Cholesky factor of a Toeplitz matrix.
!
!  Discussion:
!
!    The Toeplitz matrix must be positive semi-definite.
!
!    After factorization, A = R' * R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Michael Stewart,
!    Cholesky factorization of semi-definite Toeplitz matrices.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the Toeplitz matrix.
!
!    Output, real ( kind = 8 ) R(N,N), the upper Cholesky factor.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) g(2,n)
  real ( kind = 8 ) h(2,2)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) rho

  r(1:n,1:n) = 0.0D+00

  g(1,1:n) = a(1,1:n)
  g(2,1) = 0.0D+00
  g(2,2:n) = a(2:n,1) 

  r(1,1:n) = g(1,1:n)
  g(1,2:n) = g(1,1:n-1)
  g(1,1) = 0.0D+00
  do i = 2, n
    rho = - g(2,i) / g(1,i)
    h = reshape ( (/ 1.0D+00, rho, rho, 1.0D+00 /), (/ 2, 2 /) )
    g(1:2,i:n) = matmul ( h, g(1:2,i:n) ) &
      / sqrt ( ( 1.0D+00 - rho ) * ( 1.0D+00 + rho ) )
    r(i,i:n) = g(1,i:n)
    g(1,i+1:n) = g(1,i:n-1)
    g(1,i) = 0.0D+00
  end do

  return
end

