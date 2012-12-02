subroutine p00_a ( prob, m, n, a )

!*****************************************************************************80
!
!! P00_A returns the matrix A for any least squares problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) prob

  if ( prob == 1 ) then
    call p01_a ( m, n, a )
  else if ( prob == 2 ) then
    call p02_a ( m, n, a )
  else if ( prob == 3 ) then
    call p03_a ( m, n, a )
  else if ( prob == 4 ) then
    call p04_a ( m, n, a )
  else if ( prob == 5 ) then
    call p05_a ( m, n, a )
  else if ( prob == 6 ) then
    call p06_a ( m, n, a )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_A - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of PROB = ', prob
    stop
  end if

  return
end
subroutine p00_b ( prob, m, b )

!*****************************************************************************80
!
!! P00_B returns the right hand side B for any least squares problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Output, real ( kind = 8 ) B(M), the right hand side.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) prob

  if ( prob == 1 ) then
    call p01_b ( m, b )
  else if ( prob == 2 ) then
    call p02_b ( m, b )
  else if ( prob == 3 ) then
    call p03_b ( m, b )
  else if ( prob == 4 ) then
    call p04_b ( m, b )
  else if ( prob == 5 ) then
    call p05_b ( m, b )
  else if ( prob == 6 ) then
    call p06_b ( m, b )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_B - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of PROB = ', prob
    stop
  end if

  return
end
subroutine p00_m ( prob, m )

!*****************************************************************************80
!
!! P00_M returns the number of equations M for any least squares problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Output, integer ( kind = 4 ) M, the number of equations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) prob

  if ( prob == 1 ) then
    call p01_m ( m )
  else if ( prob == 2 ) then
    call p02_m ( m )
  else if ( prob == 3 ) then
    call p03_m ( m )
  else if ( prob == 4 ) then
    call p04_m ( m )
  else if ( prob == 5 ) then
    call p05_m ( m )
  else if ( prob == 6 ) then
    call p06_m ( m )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_M - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of PROB = ', prob
    stop
  end if

  return
end
subroutine p00_n ( prob, n )

!*****************************************************************************80
!
!! P00_N returns the number of variables N for any least squares problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Output, integer ( kind = 4 ) N, the number of variables.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) prob

  if ( prob == 1 ) then
    call p01_n ( n )
  else if ( prob == 2 ) then
    call p02_n ( n )
  else if ( prob == 3 ) then
    call p03_n ( n )
  else if ( prob == 4 ) then
    call p04_n ( n )
  else if ( prob == 5 ) then
    call p05_n ( n )
  else if ( prob == 6 ) then
    call p06_n ( n )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_N - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of PROB = ', prob
    stop
  end if

  return
end
subroutine p00_prob_num ( prob_num )

!*****************************************************************************80
!
!! P00_PROB_NUM returns the number of least squares problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) PROB_NUM, the number of problems.
!
  implicit none

  integer ( kind = 4 ) prob_num

  prob_num = 6

  return
end
subroutine p00_x ( prob, n, x )

!*****************************************************************************80
!
!! P00_X returns the least squares solution X for any least squares problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) prob
  real ( kind = 8 ) x(n)

  if ( prob == 1 ) then
    call p01_x ( n, x )
  else if ( prob == 2 ) then
    call p02_x ( n, x )
  else if ( prob == 3 ) then
    call p03_x ( n, x )
  else if ( prob == 4 ) then
    call p04_x ( n, x )
  else if ( prob == 5 ) then
    call p05_x ( n, x )
  else if ( prob == 6 ) then
    call p06_x ( n, x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_X - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of PROB = ', prob
    stop
  end if

  return
end
subroutine p01_a ( m, n, a )

!*****************************************************************************80
!
!! P01_A returns the matrix A for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m
    a(i,1) = 1.0D+00
    do j = 2, n
      a(i,j) = a(i,j-1) * real ( i, kind = 8 )
    end do
  end do

  return
end
subroutine p01_b ( m, b )

!*****************************************************************************80
!
!! P01_B returns the right hand side B for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Output, real ( kind = 8 ) B(M), the right hand side.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) b(m)

  b(1:m) = (/ 1.0D+00, 2.3D+00, 4.6D+00, 3.1D+00, 1.2D+00 /)

  return
end
subroutine p01_m ( m )

!*****************************************************************************80
!
!! P01_M returns the number of equations M for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the number of equations.
!
  implicit none

  integer ( kind = 4 ) m

  m = 5

  return
end
subroutine p01_n ( n )

!*****************************************************************************80
!
!! P01_N returns the number of variables N for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N, the number of variables.
!
  implicit none

  integer ( kind = 4 ) n

  n = 3

  return
end
subroutine p01_x ( n, x )

!*****************************************************************************80
!
!! P01_X returns the least squares solution X for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ -3.0200000D+00, 4.4914286D+00, -0.72857143D+00 /)

  return
end
subroutine p02_a ( m, n, a )

!*****************************************************************************80
!
!! P02_A returns the matrix A for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cleve Moler,
!    Numerical Computing with MATLAB,
!    SIAM, 2004,
!    ISBN13: 978-0-898716-60-3,
!    LC: QA297.M625,
!    ebook: http://www.mathworks.com/moler/chapters.html
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m
    a(i,n) = 1.0D+00
    do j = n - 1, 1, -1
      a(i,j) = a(i,j+1) * real ( i - 1, kind = 8 ) / 5.0D+00
    end do
  end do

  return
end
subroutine p02_b ( m, b )

!*****************************************************************************80
!
!! P02_B returns the right hand side B for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Output, real ( kind = 8 ) B(M), the right hand side.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) b(m)

  b(1:m) = (/ 150.697D+00, 179.323D+00, 203.212D+00, 226.505D+00, 249.633D+00, &
    281.422D+00 /)

  return
end
subroutine p02_m ( m )

!*****************************************************************************80
!
!! P02_M returns the number of equations M for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the number of equations.
!
  implicit none

  integer ( kind = 4 ) m

  m = 6

  return
end
subroutine p02_n ( n )

!*****************************************************************************80
!
!! P02_N returns the number of variables N for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N, the number of variables.
!
  implicit none

  integer ( kind = 4 ) n

  n = 3

  return
end
subroutine p02_x ( n, x )

!*****************************************************************************80
!
!! P02_X returns the least squares solution X for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ 5.7013D+00, 121.1341D+00, 152.4745D+00 /)

  return
end
subroutine p03_a ( m, n, a )

!*****************************************************************************80
!
!! P03_A returns the matrix A for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cleve Moler,
!    Numerical Computing with MATLAB,
!    SIAM, 2004,
!    ISBN13: 978-0-898716-60-3,
!    LC: QA297.M625,
!    ebook: http://www.mathworks.com/moler/chapters.html
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ), dimension ( 5, 3 ) :: a_save = reshape ( (/ &
    1.0D+00, 4.0D+00, 7.0D+00, 10.0D+00, 13.0D+00, &
    2.0D+00, 5.0D+00, 8.0D+00, 11.0D+00, 14.0D+00, &
    3.0D+00, 6.0D+00, 9.0D+00, 12.0D+00, 15.0D+00 /), (/ 5, 3 /) )

  call r8mat_copy ( m, n, a_save, a )

  return
end
subroutine p03_b ( m, b )

!*****************************************************************************80
!
!! P03_B returns the right hand side B for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Output, real ( kind = 8 ) B(M), the right hand side.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) b(m)

  b(1:m) = (/ 16.0D+00, 17.0D+00, 18.0D+00, 19.0D+00, 20.0D+00 /)

  return
end
subroutine p03_m ( m )

!*****************************************************************************80
!
!! P03_M returns the number of equations M for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the number of equations.
!
  implicit none

  integer ( kind = 4 ) m

  m = 5

  return
end
subroutine p03_n ( n )

!*****************************************************************************80
!
!! P03_N returns the number of variables N for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N, the number of variables.
!
  implicit none

  integer ( kind = 4 ) n

  n = 3

  return
end
subroutine p03_x ( n, x )

!*****************************************************************************80
!
!! P03_X returns the least squares solution X for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ -7.5555556D+00, 0.1111111D+00, 7.7777778D+00 /)

  return
end
subroutine p04_a ( m, n, a )

!*****************************************************************************80
!
!! P04_A returns the matrix A for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n
    do i = 1, m
      a(i,j) = real ( j ** ( i - 1 ), kind = 8 )
    end do
  end do

  return
end
subroutine p04_b ( m, b )

!*****************************************************************************80
!
!! P04_B returns the right hand side B for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Output, real ( kind = 8 ) B(M), the right hand side.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) b(m)

  b(1:m) = (/ 15.0D+00, 55.0D+00, 225.0D+00 /)

  return
end
subroutine p04_m ( m )

!*****************************************************************************80
!
!! P04_M returns the number of equations M for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the number of equations.
!
  implicit none

  integer ( kind = 4 ) m

  m = 3

  return
end
subroutine p04_n ( n )

!*****************************************************************************80
!
!! P04_N returns the number of variables N for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N, the number of variables.
!
  implicit none

  integer ( kind = 4 ) n

  n = 5

  return
end
subroutine p04_x ( n, x )

!*****************************************************************************80
!
!! P04_X returns the least squares solution X for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = (/ 1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /)

  return
end
subroutine p05_a ( m, n, a )

!*****************************************************************************80
!
!! P05_A returns the matrix A for problem 5.
!
!  Discussion:
!
!    A is the Hilbert matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n
    do i = 1, m
      a(i,j) = 1.0D+00 / real ( i + j - 1, kind = 8 )
    end do
  end do

  return
end
subroutine p05_b ( m, b )

!*****************************************************************************80
!
!! P05_B returns the right hand side B for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Output, real ( kind = 8 ) B(M), the right hand side.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) b(m)

  b(1) = 1.0D+00
  b(2:m) = 0.0D+00

  return
end
subroutine p05_m ( m )

!*****************************************************************************80
!
!! P05_M returns the number of equations M for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the number of equations.
!
  implicit none

  integer ( kind = 4 ) m

  m = 10

  return
end
subroutine p05_n ( n )

!*****************************************************************************80
!
!! P05_N returns the number of variables N for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N, the number of variables.
!
  implicit none

  integer ( kind = 4 ) n

  n = 10

  return
end
subroutine p05_x ( n, x )

!*****************************************************************************80
!
!! P05_X returns the least squares solution X for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = r8_mop ( i + 1 ) * real ( i, kind = 8 ) &
      * r8_choose ( n + i - 1, n - 1 ) * r8_choose ( n, n - i )
  end do

  return
end
subroutine p06_a ( m, n, a )

!*****************************************************************************80
!
!! P06_A returns the matrix A for problem 6.
!
!  Discussion:
!
!    A is a symmetric, orthogonal matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  do i = 1, m
    do j = 1, n
      angle = real ( i * j, kind = 8 ) * pi / real ( n + 1, kind = 8 )
      a(i,j) = sin ( angle )
    end do
  end do

  a(1:m,1:n) = a(1:m,1:n) * sqrt ( 2.0D+00 / real ( n + 1, kind = 8 ) )

  return
end
subroutine p06_b ( m, b )

!*****************************************************************************80
!
!! P06_B returns the right hand side B for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of equations.
!
!    Output, real ( kind = 8 ) B(M), the right hand side.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) b(m)

  b(1) = 1.0D+00
  b(2:m) = 0.0D+00

  return
end
subroutine p06_m ( m )

!*****************************************************************************80
!
!! P06_M returns the number of equations M for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the number of equations.
!
  implicit none

  integer ( kind = 4 ) m

  m = 10

  return
end
subroutine p06_n ( n )

!*****************************************************************************80
!
!! P06_N returns the number of variables N for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N, the number of variables.
!
  implicit none

  integer ( kind = 4 ) n

  n = 10

  return
end
subroutine p06_x ( n, x )

!*****************************************************************************80
!
!! P06_X returns the least squares solution X for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  do i = 1, n
    angle = real ( i, kind = 8 ) * pi / real ( n + 1, kind = 8 )
    x(i) = sin ( angle )
  end do

  x(1:n) = x(1:n) * sqrt ( 2.0D+00 / real ( n + 1, kind = 8 ) )

  return
end
