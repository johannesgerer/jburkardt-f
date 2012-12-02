function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  real ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = 8 )
  end do

  return
end
subroutine r8mat_det ( n, a, det )

!*****************************************************************************80
!
!! R8MAT_DET computes the determinant of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Helmut Spaeth.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 125-127.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) piv(1)
  real ( kind = 8 ) t

  b(1:n,1:n) = a(1:n,1:n)

  det = 1.0D+00

  do k = 1, n

    piv = maxloc ( abs ( b(k:n,k) ) )

    m = piv(1) + k - 1

    if ( m /= k ) then
      det = - det
      t      = b(m,k)
      b(m,k) = b(k,k)
      b(k,k) = t
    end if

    det = det * b(k,k)

    if ( b(k,k) /= 0.0D+00 ) then

      b(k+1:n,k) = -b(k+1:n,k) / b(k,k)

      do j = k + 1, n
        if ( m /= k ) then
          t      = b(m,j)
          b(m,j) = b(k,j)
          b(k,j) = t
        end if
        b(k+1:n,j) = b(k+1:n,j) + b(k+1:n,k) * b(k,j)
      end do

    end if

  end do

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 June 2004
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
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
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine simplex_coordinates1 ( n, x )

!*****************************************************************************80
!
!! SIMPLEX_COORDINATES1 computes the Cartesian coordinates of simplex vertices.
!
!  Discussion:
!
!    The simplex will have its centroid at 0;
!
!    The sum of the vertices will be zero.
!
!    The distance of each vertex from the origin will be 1.
!
!    The length of each edge will be constant.
!
!    The dot product of the vectors defining any two vertices will be - 1 / N.
!    This also means the angle subtended by the vectors from the origin
!    to any two distinct vertices will be arccos ( - 1 / N ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, real ( kind = 8 ) X(N,N+1), the coordinates of the vertices
!    of a simplex in N dimensions.  
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n,n+1)

  x(1:n,1:n+1) = 0.0D+00

  do i = 1, n
!
!  Set X(I,I) so that sum ( X(1:I,I)**2 ) = 1.
!
    x(i,i) = sqrt ( 1.0D+00 - sum ( x(1:i-1,i)**2 ) )
!
!  Set X(I,J) for J = I+1 to N+1 by using the fact that XI dot XJ = - 1 / N 
!
    do j = i + 1, n + 1
      x(i,j) = ( - 1.0D+00 / real ( n, kind = 8 ) &
        - dot_product ( x(1:i-1,i), x(1:i-1,j) ) ) / x(i,i)
    end do

  end do

  return
end
subroutine simplex_coordinates2 ( n, x )

!*****************************************************************************80
!
!! SIMPLEX_COORDINATES2 computes the Cartesian coordinates of simplex vertices.
!
!  Discussion:
!
!    This routine uses a simple approach to determining the coordinates of
!    the vertices of a regular simplex in n dimensions.
!
!    We want the vertices of the simplex to satisfy the following conditions:
!
!    1) The centroid, or average of the vertices, is 0.
!    2) The distance of each vertex from the centroid is 1.
!       By 1), this is equivalent to requiring that the sum of the squares
!       of the coordinates of any vertex be 1.
!    3) The distance between any pair of vertices is equal (and is not zero.)
!    4) The dot product of any two coordinate vectors for distinct vertices
!       is -1/N; equivalently, the angle subtended by two distinct vertices
!       from the centroid is arccos ( -1/N).
!
!    Note that if we choose the first N vertices to be the columns of the
!    NxN identity matrix, we are almost there.  By symmetry, the last column
!    must have all entries equal to some value A.  Because the square of the
!    distance between the last column and any other column must be 2 (because
!    that's the distance between any pair of columns), we deduce that
!    (A-1)^2 + (N-1)*A^2 = 2, hence A = (1-sqrt(1+N))/N.  Now compute the 
!    centroid C of the vertices, and subtract that, to center the simplex 
!    around the origin.  Finally, compute the norm of one column, and rescale 
!    the matrix of coordinates so each vertex has unit distance from the origin.
!
!    This approach devised by John Burkardt, 19 September 2010.  What,
!    I'm not the first?
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, real ( kind = 8 ) X(N,N+1), the coordinates of the vertices
!    of a simplex in N dimensions.  
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) s
  real ( kind = 8 ) x(n,n+1)

  x(1:n,1:n+1) = 0.0D+00

  do j = 1, n
    x(j,j) = 1.0D+00
  end do

  a = ( 1.0D+00 - sqrt ( 1.0D+00 + real ( n, kind = 8 ) ) ) &
    / real ( n, kind = 8 )

  x(1:n,n+1) = a
!
!  Now adjust coordinates so the centroid is at zero.
!
  c(1:n) = sum ( x(1:n,1:n+1), dim = 2 ) / real ( n + 1, kind = 8 )

  do j = 1, n + 1
    x(1:n,j) = x(1:n,j) - c(1:n)
  end do
!
!  Now scale so each column has norm 1.
!
  s = sqrt ( sum ( x(1:n,1)**2 ) )

  x(1:n,1:n+1) = x(1:n,1:n+1) / s

  return
end
subroutine simplex_volume ( n, x, volume )

!*****************************************************************************80
!
!! SIMPLEX_VOLUME computes the volume of a simplex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(N,N+1), the coordinates of the vertices
!    of a simplex in N dimensions.  
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the simplex.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) volume
  real ( kind = 8 ) x(n,n+1)

  a(1:n,1:n) = x(1:n,1:n)
  do j = 1, n
    a(1:n,j) = a(1:n,j) - x(1:n,n+1)
  end do

  call r8mat_det ( n, a, det )

  volume = abs ( det )
  do i = 1, n
    volume = volume / real ( i, kind = 8 )
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
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
