subroutine fem_basis_1d ( i, j, x, lij )

!*****************************************************************************80
!
!! FEM_BASIS_1D evaluates an arbitrary 1D basis function.
!
!  Discussion:
!
!    Given the maximum degree D for the polynomial basis defined
!    on a reference interval, we have D + 1 monomials
!    of degree at most D.  In each barycentric coordinate, we define
!    D+1 points, so that 0 <= I, J <= D and I+J = D, with
!    (I,J) corresponding to 
!    * the basis point X(I,J) = ( I/D );
!    * the basis monomial P(I,J)(X) = X^I.
!
!    For example, with D = 2, we have simply:
!
!      A---B---C
!
!    with 
!
!       I J    X      P(I,J)(X) 
!
!    A (0 2) ( 0.0 )  1
!    B (1 1) ( 0.5 )  x
!    C (2 0) ( 1.0 )  x^2
!
!    Now instead of the monomials P(I,J)(X), we want a set of
!    polynomials L(I,J)(X) which span the same space, but have
!    the Lagrange property, namely L(I,J) (X) is 1 if X is
!    equal to X(I,J), and 0 if X is equal to any other 
!    of the basis points.
!    
!    This is easily arranged.  Given an index (I,J), we compute
!    1) I factors of the form (   X -0/D) * (   X -1/D) * ... * (   X -(I-1)/D);
!    2) J factors of the form ((1-X)-0/D) * ((1-X)-1/D) * ... * ((1-X)-(J-1)/D).
!
!    This results in the product of I+J linear factors, in other words,
!    a polynomial of degree D.  This polynomial is 0 at all basis points
!    except X(I,J).  If we divide this polynomial by its value at
!    the basis point, we arrive at the desired Lagrange polynomial
!    L(I,J)(X). 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the integer barycentric coordinates of
!    the basis function, 0 <= I, J.  The polynomial degree D = I + J.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) LIJ, the value of the basis function at X.
!
  implicit none

  real ( kind = 8 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lij
  integer ( kind = 4 ) p
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  
  d = i + j
  lij = 1.0D+00
  c = 1.0D+00
  do p = 0, i - 1
    lij = lij * ( d * x - p )
    c = c     * (     i - p )
  end do
  w = 1.0D+00 - x
  do p = 0, j - 1
    lij = lij * ( d * w - p )
    c = c     * (     j - p )
  end do
  
  lij = lij / c

  return
end
subroutine fem_basis_2d ( i, j, k, x, y, lijk )

!*****************************************************************************80
!
!! FEM_BASIS_2D evaluates an arbitrary triangular basis function.
!
!  Discussion:
!
!    Given the maximum degree D for the polynomial basis defined
!    on a reference triangle, we have ( ( D + 1 ) * ( D + 2 ) ) / 2 monomials
!    of degree at most D.  In each barycentric coordinate, we define
!    D+1 planes, so that 0 <= I, J, K <= D and I+J+K = D, with
!    (I,J,K) corresponding to 
!    * the basis point (X,Y)(I,J,K) = ( I/D, J/D );
!    * the basis monomial P(I,J,K)(X,Y) = X^I Y^J.
!
!    For example, with D = 2, we have simply:
!
!    F
!    |\
!    C-E
!    |\|\
!    A-B-D
!
!    with 
!
!       I J K    X    Y    P(I,J,K)(X,Y) 
!
!    A (0 0 2) (0.0, 0.0)  1
!    B (1 0 1) (0.5, 0.0)  x
!    C (0 1 1) (0.0, 0.5)  y
!    D (2 0 0) (1.0, 0.0)  x^2
!    E (1 1 0) (0.5, 0.5)  x y
!    F (0 2 0) (0.0, 1.0)  y^2
!
!    Now instead of the monomials P(I,J,K)(X,Y), we want a set of
!    polynomials L(I,J,K)(X,Y) which span the same space, but have
!    the Lagrange property, namely L(I,J,K) (X,Y) is 1 if (X,Y) is
!    equal to (X,Y)(I,J,K), and 0 if (X,Y) is equal to any other 
!    of the basis points.
!    
!    This is easily arranged.  Given an index (I,J,K), we compute
!    1) I factors of the form (X-0)   * (X-1/D)   * ... * (X-(I-1)/D);
!    2) J factors of the form (Y-0)   * (Y-1/D)   * ... * (Y-(J-1)/D);
!    3) K factors of the form ((1-X-Y)-0/D) * ((1-X-Y)-1/D) * ... 
!       * ((1-X-Y)-(K-1)/D).
!
!    This results in the product of I+J+K linear factors, in other words,
!    a polynomial of degree D.  This polynomial is 0 at all basis points
!    except (X,Y)(I,J,K).  If we divide this polynomial by its value at
!    the basis point, we arrive at the desired Lagrange polynomial
!    L(I,J,K)(X,Y). 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, K, the integer barycentric coordinates of
!    the basis function, 0 <= I, J, K.  The polynomial degree D = I + J + K.
!
!    Input, real ( kind = 8 ) X, Y, the evaluation point.
!
!    Output, real ( kind = 8 ) LIJK, the value of the basis function at (X,Y).
!
  implicit none

  real ( kind = 8 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) lijk
  integer ( kind = 4 ) p
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  
  d = i + j + k
  lijk = 1.0D+00
  c = 1.0D+00
  do p = 0, i - 1
    lijk = lijk * ( d * x - p )
    c = c       * (     i - p )
  end do
  do p = 0, j - 1
    lijk = lijk * ( d * y - p )
    c = c       * (     j - p )
  end do
  w = 1.0D+00 - x - y
  do p = 0, k - 1
    lijk = lijk * ( d * w - p )
    c = c       * (     k - p )
  end do
  
  lijk = lijk / c

  return
end
subroutine fem_basis_3d ( i, j, k, l, x, y, z, lijkl )

!*****************************************************************************80
!
!! FEM_BASIS_3D evaluates an arbitrary tetrahedral basis function.
!
!  Discussion:
!
!    Given the maximum degree D for the polynomial basis defined
!    on a reference tetrahedron, we have 
!    ( D + 1 ) * ( D + 2 ) * ( D + 3 ) / 6 monomials
!    of degree at most D.  In each barycentric coordinate, we define
!    D+1 planes, so that 0 <= I, J, K, L <= D and I+J+K+L = D, with
!    (I,J,K,L) corresponding to 
!    * the basis point (X,Y,Z)(I,J,K,L) = ( I/D, J/D, K/D );
!    * the basis monomial P(I,J,K,L)(X,Y,Z) = X^I Y^J Z^K.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, K, L, the integer barycentric 
!    coordinates of the basis function, 0 <= I, J, K, L. 
!    The polynomial degree D = I + J + K + L.
!
!    Input, real ( kind = 8 ) X, Y, Z, the evaluation point.
!
!    Output, real ( kind = 8 ) LIJKL, the value of the basis function 
!    at (X,Y,Z).
!
  implicit none

  real ( kind = 8 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) lijkl
  integer ( kind = 4 ) p
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  
  d = i + j + k + l
  lijkl = 1.0D+00
  c = 1.0D+00
  do p = 0, i - 1
    lijkl = lijkl * ( d * x - p )
    c = c         * (     i - p )
  end do
  do p = 0, j - 1
    lijkl = lijkl * ( d * y - p )
    c = c         * (     j - p )
  end do
  do p = 0, k - 1
    lijkl = lijkl * ( d * z - p )
    c = c         * (     k - p )
  end do
  w = 1.0D+00 - x - y - z
  do p = 0, l - 1
    lijkl = lijkl * ( d * w - p )
    c = c         * (     l - p )
  end do
  
  lijkl = lijkl / c

  return
end
subroutine fem_basis_md ( m, i, x, l )

!*****************************************************************************80
!
!! FEM_BASIS_MD evaluates an arbitrary M-dimensional basis function.
!
!  Discussion:
!
!    This routine evaluates the generalization of the formula used for
!    the 1D, 2D and 3D cases.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) I(M+1), the integer barycentric 
!    coordinates of the basis function, 0 <= I(1:M+1). 
!    The polynomial degree D = sum(I(1:M+1)).
!
!    Input, real ( kind = 8 ) X(M), the evaluation point.
!
!    Output, real ( kind = 8 ) L, the value of the basis function at X.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i(m+1)
  real ( kind = 8 ) l
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  real ( kind = 8 ) w
  real ( kind = 8 ) x(m)
  
  d = sum ( i(1:m+1) )

  l = 1.0D+00
  c = 1.0D+00

  do q = 1, m
    do p = 0, i(q) - 1
      l = l * ( d * x(q) - p )
      c = c * (     i(q) - p )
    end do
  end do

  w = 1.0D+00 - sum ( x(1:m) )

  do p = 0, i(m+1) - 1
    l = l * ( d * w      - p )
    c = c * (     i(m+1) - p )
  end do
  
  l = l / c

  return
end
subroutine fem_basis_prism_triangle ( i, j, xyz, b )

!*****************************************************************************80
!
!! FEM_BASIS_PRISM_TRIANGLE evaluates a triangular prism basis function.
!
!  Discussion:
!
!    The element is a 3D prism, formed from a triangular base in the
!    XY plane that is extended vertically in the Z direction.
!
!    I(1:3) are the integer barycentric coordinates of a point in the
!    triangle.  I(1) + I(2) + I(3) = DI, the degree of the triangular
!    basis function BI.  X = I(1) / DI, Y = I(2) / DI.
!    The triangle is assumed to be the unit reference
!    triangle 0 <= X <= 1, 0 <= Y <= 1, 0 <= X + Y <= 1.
!
!    J(1:2) are the integer barycentric coordinates of a point in the
!    line segment.  J(1) + J(2) = DJ, the degree of the linear basis 
!    function BJ.  Z = J(1) / DJ.  
!    The line is assumed to be the unit line 0 <= Z <= 1.
!
!    The degree of the basis function B = BI * BJ is D = DI + DJ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I(3), the integer barycentric coordinates of
!    the triangular basis function, 0 <= I(*).  
!    The polynomial degree DI = I(1) + I(2) + I(3).
!
!    Input, integer ( kind = 4 ) J(2), the integer barycentric coordinates of
!    the linear basis function, 0 <= J(*).  
!    The polynomial degree DJ = J(1) + J(2).
!
!    Input, real ( kind = 8 ) XYZ(3), the evaluation point.
!
!    Output, real ( kind = 8 ) B, the value of the basis function at XYZ.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) bi
  real ( kind = 8 ) bj
  integer ( kind = 4 ) i(3)
  integer ( kind = 4 ) j(2)
  real ( kind = 8 ) xyz(3)

  call fem_basis_2d ( i(1), i(2), i(3), xyz(1), xyz(2), bi )

  call fem_basis_1d ( j(1), j(2), xyz(3), bj )
  
  b = bi * bj

  return
end
function r8_fraction ( i, j )

!*****************************************************************************80
!
!! R8_FRACTION uses real arithmetic on an integer ratio.
!
!  Discussion:
!
!    Given integer variables I and J, both FORTRAN and C will evaluate 
!    an expression such as "I/J" using what is called "integer division",
!    with the result being an integer.  It is often convenient to express
!    the parts of a fraction as integers but expect the result to be computed
!    using real arithmetic.  This function carries out that operation.
!
!  Example:
!
!       I     J   I/J  R8_FRACTION
!
!       1     2     0  0.5
!       7     4     1  1.75
!       8     4     2  2.00
!       9     4     2  2.25
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the arguments.
!
!    Output, real ( kind = 8 ) R8_FRACTION, the value of the ratio.
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_fraction

  r8_fraction = real ( i, kind = 8 ) / real ( j, kind = 8 )

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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
