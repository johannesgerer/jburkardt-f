function flin ( n, j, l, f, x, nf, v, q0, q1 )

!*****************************************************************************80
!
!! FLIN is the function of one variable to be minimized by MINNY.
!
!  Discussion:
!
!    In fact, what is happening is that the scalar function F(X),
!    where X is an N dimensional vector, is being minimized along a 
!    fixed line.
!
!  Modified:
!
!    22 May 2006
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, integer ( kind = 4 ) J, indicates the kind of search.
!    If J is nonzero, then the search is linear in the direction of V(*,J).
!    If J is zero, then the search is parabolic, based on X, Q0 and Q1.
!
!    Input, real ( kind = 8 ) L, is the parameter determining the particular
!    point at which F is to be evaluated.  
!    For a linear search ( J is nonzero ), L is the size of the step
!    along the direction V(*,J).
!    For a quadratic search ( J is zero ), L is a parameter which specifies
!    a point in the plane of X, Q0 and Q1.
!
!    Input, external F, is the name of the function to be minimized.
!    The function should have the form 
!      function f(x,n)
!      integer ( kind = 4 ) ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!    and accepts X and N as input, returning in F the function value.
!
!    Input, real ( kind = 8 ) X(N), the base point of the search.
!
!    Input/output, integer ( kind = 4 ) NF, the function evaluation counter.
!
!    Input, real ( kind = 8 ) V(N,N), a matrix whose columns constitute 
!    search directions.  If J is nonzero, then a linear search is being 
!    carried out in the direction of V(*,J).
!
!    Input, real ( kind = 8 ) Q0(N), Q1(N), two auxiliary points used to
!    determine the plane when a quadratic search is performed.
!
!    Output, real ( kind = 8 ) FLIN, the value of the function at the 
!    given point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), external :: f
  real ( kind = 8 ) flin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) l
  integer ( kind = 4 ) nf
  real ( kind = 8 ) q0(n)
  real ( kind = 8 ) q1(n)
  real ( kind = 8 ) qa
  real ( kind = 8 ) qb
  real ( kind = 8 ) qc
  real ( kind = 8 ) qd0
  real ( kind = 8 ) qd1
  real ( kind = 8 ) qf1
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) x(n)

  common /q/ qa, qb, qc, qd0, qd1, qf1

  if ( j /= 0 ) then
!
!  The search is linear.
!
    t(1:n) = x(1:n) + l * v(1:n,j)

  else
!
!  The search is along a parabolic space curve.
!
    qa = ( l * ( l - qd1 ) ) / ( qd0 * ( qd0 + qd1 ) )
    qb = ( ( l + qd0 ) * ( qd1 - l ) ) / ( qd0 * qd1 )
    qc = ( l * ( l + qd0 ) ) / ( qd1 * ( qd0 + qd1 ) )

    t(1:n) = qa * q0(1:n) + qb * x(1:n) + qc * q1(1:n)

  end if
!
!  The function evaluation counter NF is incremented.
!
  nf = nf + 1
!
!  Evaluate the function.
!
  flin = f ( t, n )

  return
end
subroutine minfit ( m, n, tol, ab, q )

!*****************************************************************************80
!
!! MINFIT computes the singular value decomposition of an N by N array.
!
!  Discussion:
!
!    This is an improved version of the EISPACK routine MINFIT
!    restricted to the case M = N and P = 0.
!
!    The singular values of the array AB are returned in Q.  AB is
!    overwritten with the orthogonal matrix V such that u.diag(q) = ab.v,
!    where U is another orthogonal matrix.
!
!    Thanks to Andreas Zuend for pointing out a potential for overflow in
!    the computation z = sqrt ( f*f + 1 ), 22 March 2012.
!
!  Modified:
!
!    22 March 2012
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, Jack Dongarra, Burton Garbow, Y Ikebe, 
!    Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the leading dimension of AB, which must be
!    at least N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix AB.
!
!    Input, real ( kind = 8 ) TOL, a tolerance which determines when a vector
!    (a column or part of a column of the matrix) may be considered
!    "essentially" equal to zero.
!
!    Input/output, real ( kind = 8 ) AB(M,N).  On input, an N by N array whose
!    singular value decomposition is desired.  On output, ?
!
!    Input, real ( kind = 8 ) Q(N), ?
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) ab(m,n)
  real ( kind = 8 ) c
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) eps
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kt
  integer ( kind = 4 ), parameter :: kt_max = 30
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) ll2
  real ( kind = 8 ) q(n)
  real ( kind = 8 ) r8_hypot
  real ( kind = 8 ) s
  real ( kind = 8 ) temp
  real ( kind = 8 ) tol
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!
!  Householder's reduction to bidiagonal form.
!
  if ( n == 1 ) then
    q(1) = ab(1,1)
    ab(1,1) = 1.0D+00
    return
  end if

  eps = epsilon ( eps )
  g = 0.0D+00
  x = 0.0D+00

  do i = 1, n

    e(i) = g
    l = i + 1

    s = sum ( ab(1:n,i)**2 )

    g = 0.0D+00

    if ( tol <= s ) then

      f = ab(i,i)

      g = sqrt ( s )
      if ( 0.0D+00 <= f ) then
        g = -g
      end if

      h = f * g - s
      ab(i,i) = f - g

      do j = l, n

        f = dot_product ( ab(i:n,i), ab(i:n,j) ) / h

        ab(i:n,j) = ab(i:n,j) + f * ab(i:n,i)

      end do 

    end if

    q(i) = g

    s = sum ( ab(i,l:n)**2 )

    g = 0.0D+00

    if ( tol <= s ) then

      if ( i /= n ) then
        f = ab(i,i+1)
      end if

      g = sqrt ( s )
      if ( 0.0D+00 <= f ) then
        g = - g
      end if

      h = f * g - s

      if ( i /= n ) then

        ab(i,i+1) = f - g
        e(l:n) = ab(i,l:n) / h

        do j = l, n

          s = dot_product ( ab(j,l:n), ab(i,l:n) )

          ab(j,l:n) = ab(j,l:n) + s * e(l:n)

        end do

      end if

    end if

    y = abs ( q(i) ) + abs ( e(i) )

  end do

  x = max ( x, y )
!
!  Accumulation of right-hand transformations.
!
  ab(n,n) = 1.0D+00
  g = e(n)
  l = n

  do ii = 2, n

    i = n - ii + 1

    if ( g /= 0.0D+00 ) then

      h = ab(i,i+1) * g

      ab(l:n,i) = ab(i,l:n) / h

      do j = l, n

        s = dot_product ( ab(i,l:n), ab(l:n,j) )

        ab(l:n,j) = ab(l:n,j) + s * ab(l:n,i)

      end do

    end if

    ab(i,l:n) = 0.0D+00
    ab(l:n,i) = 0.0D+00
    ab(i,i) = 1.0D+00

    g = e(i)

  end do

  l = i
!
!  Diagonalization of the bidiagonal form.
!
  eps = eps * x

  do kk = 1, n

    k = n - kk + 1
    kt = 0

10  continue

    kt = kt + 1

    if ( kt_max < kt ) then
      e(k) = 0.0D+00
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MINFIT - Warning!'
      write ( *, '(a)' ) '  The QR algorithm failed to converge.'
    end if

    do ll2 = 1, k

      l2 = k - ll2 + 1
      l = l2

      if ( abs ( e(l) ) <= eps ) then
        go to 20
      end if

      if ( l /= 1 ) then
        if ( abs ( q(l-1) ) <= eps ) then
          exit
        end if
      end if

    end do
!
!  Cancellation of E(L) if 1 < L.
!
    c = 0.0D+00
    s = 1.0D+00

    do i = l, k

      f = s * e(i)
      e(i) = c * e(i)
      if ( abs ( f ) <= eps ) then
        exit
      end if
      g = q(i)
!
!  q(i) = h = sqrt(g*g + f*f).
!
      h = r8_hypot ( f, g )

      q(i) = h

      if ( h == 0.0D+00 ) then
        g = 1.0D+00
        h = 1.0D+00
      end if

      c = g / h
      s = - f / h

    end do
!
!  Test for convergence.
!
20  continue

    z = q(k)

    if ( l == k ) then
      if ( z < 0.0D+00 ) then
        q(k) = - z
        ab(1:n,k) = - ab(1:n,k)
      end if
      cycle
    end if
!
!  Shift from bottom 2*2 minor.
!
    x = q(l)
    y = q(k-1)
    g = e(k-1)
    h = e(k)
    f = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) / ( 2.0D+00 * h * y )

    g = r8_hypot ( f, 1.0D+00 )

    if ( f < 0.0D+00 ) then
      temp = f - g
    else
      temp = f + g
    end if

    f = ( ( x - z ) * ( x + z ) + h * ( y / temp - h ) ) / x
!
!  Next QR transformation.
!
    c = 1.0D+00
    s = 1.0D+00

    do i = l + 1, k

      g = e(i)
      y = q(i)
      h = s * g
      g = g * c

      z = r8_hypot ( f, h )

      e(i-1) = z

      if ( z == 0.0D+00 ) then
        f = 1.0D+00
        z = 1.0D+00
      end if

      c = f / z
      s = h / z
      f =   x * c + g * s
      g = - x * s + g * c
      h = y * s
      y = y * c

      do j = 1, n
        x = ab(j,i-1)
        z = ab(j,i)
        ab(j,i-1) = x * c + z * s
        ab(j,i) = - x * s + z * c
      end do

      z = r8_hypot ( f, h )

      q(i-1) = z

      if ( z == 0.0D+00 ) then
        f = 1.0D+00
        z = 1.0D+00
      end if

      c = f / z
      s = h / z
      f = c * g + s * y
      x = - s * g + c * y

    end do

    e(l) = 0.0D+00
    e(k) = f
    q(k) = x
    go to 10

  end do

  return
end
subroutine minny ( n, j, nits, d2, x1, f1, fk, f, x, t, h, v, q0, q1 )

!*****************************************************************************80
!
!! MINNY minimizes a scalar function of N variables along a line.
!
!  Discussion:
!
!    MINNY minimizes F along the line from X in the direction V(*,J) unless
!    J is less than 1, when a quadratic search is made in the plane
!    defined by q0,q1,x.
!
!    If fk = .true., then f1 is flin(x1).  Otherwise x1 and f1 are ignored
!    on entry unless final fx is greater than f1.
!
!  Modified:
!
!    22 May 2006
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, integer ( kind = 4 ) J, indicates the kind of search.
!    If J is positive, then the search is linear in the direction of V(*,J).
!    If J is zero, then the search is parabolic, based on X, Q0 and Q1.
!
!    Input, integer ( kind = 4 ) NITS, the maximum number of times the interval 
!    may be halved to retry the calculation.
!
!    Input, real ( kind = 8 ) D2, is either zero, or an approximation to 
!    the value of (1/2) times the second derivative of F.
!
!    Input/output, real ( kind = 8 ) X1, on entry, an estimate of the 
!    distance from X to the minimum along V(*,J), or, if J = 0, a curve.  
!    On output, the distance between X and the minimizer that was found.
!
!    Input/output, real ( kind = 8 ) F1, ?
!
!    Input, logical FK; if FK is TRUE, then on input F1 contains 
!    the value FLIN(X1).
!
!    Input, external real ( kind = 8 ) F, is the name of the function to 
!    be minimized.  The function should have the form 
!      function f(x,n)
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!    and accepts X and N as input, returning in F the function value.
!
!    ?, real ( kind = 8 ) X(N), ?
!
!    ?, real ( kind = 8 ) T, ?
!
!    ?, real ( kind = 8 ) H, ?
!
!    Input, real ( kind = 8 ) V(N,N), a matrix whose columns are direction
!    vectors along which the function may be minimized.
!
!    ?, real ( kind = 8 ) Q0(N), ?
!
!    ?, real ( kind = 8 ) Q1(N), ?
!
  implicit none

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dmin
  logical dz
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  logical fk
  real ( kind = 8 ) flin
  real ( kind = 8 ) fm
  real ( kind = 8 ) fx
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) ldt
  real ( kind = 8 ) m2
  real ( kind = 8 ) m4
  real ( kind = 8 ) machep
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nits
  integer ( kind = 4 ) nl
  real ( kind = 8 ) q0(n)
  real ( kind = 8 ) q1(n)
  real ( kind = 8 ) qa
  real ( kind = 8 ) qb
  real ( kind = 8 ) qc
  real ( kind = 8 ) qd0
  real ( kind = 8 ) qd1
  real ( kind = 8 ) qf1
  real ( kind = 8 ) s
  real ( kind = 8 ) sf1
  real ( kind = 8 ) small
  real ( kind = 8 ) sx1
  real ( kind = 8 ) t
  real ( kind = 8 ) t2
  real ( kind = 8 ) temp
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xm

  common /global/ fx, ldt, dmin, nf, nl
  common /q/ qa, qb, qc, qd0, qd1, qf1

  machep = epsilon ( machep )
  small = machep**2
  m2 = sqrt ( machep )
  m4 = sqrt ( m2 )
  sf1 = f1
  sx1 = x1
  k = 0
  xm = 0.0D+00
  fm = fx
  f0 = fx
  dz = ( d2 < machep )
!
!  Find the step size.
! 
  s = sqrt ( sum ( x(1:n)**2 ) )

  if ( dz ) then
    temp = dmin
  else
    temp = d2
  end if

  t2 = m4 * sqrt ( abs ( fx ) / temp + s * ldt ) + m2 * ldt
  s = m4 * s + t
  if ( dz .and. s < t2 ) then
    t2 = s
  end if

  t2 = max ( t2, small )
  t2 = min ( t2, 0.01D+00 * h )

  if ( fk .and. f1 <= fm ) then
    xm = x1
    fm = f1
  end if

  if ( .not. fk .or. abs ( x1 ) < t2 ) then

    if ( 0.0D+00 <= x1 ) then
      temp = 1.0D+00
    else
      temp = - 1.0D+00
    end if

    x1 = temp * t2
    f1 = flin ( n, j, x1, f, x, nf, v, q0, q1 )

  end if

  if ( f1 <= fm ) then
    xm = x1
    fm = f1
  end if
!
!  Evaluate FLIN at another point and estimate the second derivative.
!
10 continue

  if ( dz ) then

    if ( f1 <= f0 ) then
      x2 = 2.0D+00 * x1
    else
      x2 = - x1
    end if

    f2 = flin ( n, j, x2, f, x, nf, v, q0, q1 )

    if ( f2 <= fm ) then
      xm = x2
      fm = f2
    end if

    d2 = ( x2 * ( f1 - f0 ) - x1 * ( f2 - f0 ) ) / ( ( x1 * x2 ) * ( x1 - x2 ) )

  end if
!
!  Estimate the first derivative at 0.
!
  d1 = ( f1 - f0 ) / x1 - x1 * d2
  dz = .true.
!
!  Predict the minimum.
!
  if ( d2 <= small ) then

    if ( 0.0D+00 <= d1 ) then
      x2 = - h
    else
      x2 = h
    end if

  else

    x2 = ( - 0.5D+00 * d1 ) / d2

  end if

  if ( h < abs ( x2 ) ) then

    if ( x2 <= 0.0D+00 ) then
      x2 = - h
    else
      x2 = h
    end if

  end if
!
!  Evaluate F at the predicted minimum.
!
  do

    f2 = flin ( n, j, x2, f, x, nf, v, q0, q1 )

    if ( nits <= k .or. f2 <= f0 ) then
      exit
    end if

    k = k + 1

    if ( f0 < f1 .and. 0.0D+00 < x1 * x2 ) then
      go to 10
    end if

    x2 = 0.5D+00 * x2

  end do
!
!  Increment the one-dimensional search counter.
!
  nl = nl + 1

  if ( fm < f2 ) then
    x2 = xm
  else
    fm = f2
  end if
!
!  Get a new estimate of the second derivative.
!
  if ( small < abs ( x2 * ( x2 - x1 ) ) ) then
    d2 = ( x2 * ( f1 - f0 ) - x1 * ( fm - f0 ) ) / ( ( x1 * x2 ) * ( x1 - x2 ) )
  else
    if ( 0 < k ) then
      d2 = 0.0D+00
    end if
  end if

  d2 = max ( d2, small )

  x1 = x2
  fx = fm

  if ( sf1 < fx ) then
    fx = sf1
    x1 = sx1
  end if
!
!  Update X for linear but not parabolic search.
!
  if ( j /= 0 ) then

    x(1:n) = x(1:n) + x1 * v(1:n,j)

  end if

  return
end
function praxis ( t0, h0, n, prin, x, f )

!*****************************************************************************80
!
!! PRAXIS seeks an N-dimensional minimizer X of a scalar function F(X).
!
!  Discussion:
!
!    PRAXIS returns the minimum of the function F(X,N) of N variables
!    using the principal axis method.  The gradient of the function is
!    not required.
!
!    The approximating quadratic form is
!
!      Q(x') = F(x,n) + (1/2) * (x'-x)' * A * (x'-x)
!
!    where X is the best estimate of the minimum and 
!
!      A = inverse(V') * D * inverse(V)
!
!    V(*,*) is the matrix of search directions; 
!    D(*) is the array of second differences.  
!
!    If F(X) has continuous second derivatives near X0, then A will tend 
!    to the hessian of F at X0 as X approaches X0.
!
!    Thanks to Andreas Zuend for pointing out an error in the form of the
!    call to the routine r8mat_print (), 22 March 2012.
!
!  Modified:
!
!    22 March 2012
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T0, is a tolerance.  PRAXIS attempts to return 
!    praxis = f(x) such that if X0 is the true local minimum near X, then
!    norm ( x - x0 ) < T0 + sqrt ( EPSILON ( X ) ) * norm ( X ),
!    where EPSILON ( X ) is the machine precision for X.
!
!    Input, real ( kind = 8 ) H0, is the maximum step size.  H0 should be 
!    set to about the maximum distance from the initial guess to the minimum.
!    If H0 is set too large or too small, the initial rate of
!    convergence may be slow.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, integer ( kind = 4 ) PRIN, controls printing intermediate results.
!    0, nothing is printed.
!    1, F is printed after every n+1 or n+2 linear minimizations.  
!       final X is printed, but intermediate X is printed only 
!       if N is at most 4.
!    2, the scale factors and the principal values of the approximating 
!       quadratic form are also printed.
!    3, X is also printed after every few linear minimizations.
!    4, the principal vectors of the approximating quadratic form are 
!       also printed.
!
!    Input/output, real ( kind = 8 ) X(N), is an array containing on entry a
!    guess of the point of minimum, on return the estimated point of minimum.
!
!    Input, external real ( kind = 8 ) F, is the name of the function to be
!    minimized.  The function should have the form 
!      function f(x,n)
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!    and accepts X and N as input, returning in F the function value.
!
!    Output, real ( kind = 8 ) PRAXIS, the function value at the minimizer.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) df
  real ( kind = 8 ) dmin
  real ( kind = 8 ) dn
  real ( kind = 8 ) dni
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) f1
  real ( kind = 8 ) fx
  real ( kind = 8 ) h
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  logical illc
  integer ( kind = 4 ), save :: iseed = 1234567
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) kl
  integer ( kind = 4 ) kt
  integer ( kind = 4 ) ktm
  real ( kind = 8 ) large
  real ( kind = 8 ) ldfac
  real ( kind = 8 ) lds
  real ( kind = 8 ) ldt
  real ( kind = 8 ) m2
  real ( kind = 8 ) m4
  real ( kind = 8 ) machep
  integer ( kind = 4 ) nits
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nf
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) q0(n)
  real ( kind = 8 ) q1(n)
  real ( kind = 8 ) qa
  real ( kind = 8 ) qb
  real ( kind = 8 ) qc
  real ( kind = 8 ) qd0
  real ( kind = 8 ) qd1
  real ( kind = 8 ) qf1
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) scbd
  real ( kind = 8 ) sf
  real ( kind = 8 ) sl
  real ( kind = 8 ) small
  real ( kind = 8 ) t
  real ( kind = 8 ) t0
  real ( kind = 8 ) t2
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) value
  real ( kind = 8 ) vlarge
  real ( kind = 8 ) vsmall
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  common /global/ fx, ldt, dmin, nf, nl
  common /q/ qa, qb, qc, qd0, qd1, qf1
!
!  Initialization.
!
  machep = epsilon ( machep )
  small = machep * machep
  vsmall = small * small
  large = 1.0D+00 / small
  vlarge = 1.0D+00 / vsmall
  m2 = sqrt ( machep )
  m4 = sqrt ( m2 )
!
!  Heuristic numbers:
!
!  If the axes may be badly scaled (which is to be avoided if
!  possible), then set SCBD = 10.  Otherwise set SCBD = 1.
!
!  If the problem is known to be ill-conditioned, initialize ILLC = true.
!
!  KTM is the number of iterations without improvement before the
!  algorithm terminates.  KTM = 4 is very cautious; usually KTM = 1
!  is satisfactory.
!
  scbd = 1.0D+00
  illc = .false.
  ktm = 1

  if ( illc ) then
    ldfac = 0.1D+00
  else
    ldfac = 0.01D+00
  end if

  kt = 0
  nl = 0
  nf = 1
  fx = f(x,n)
  qf1 = fx
  t = small + abs ( t0 )
  t2 = t
  dmin = small
  h = h0
  h = max ( h, 100.0D+00 * t )
  ldt = h
!
!  The initial set of search directions V is the identity matrix.
!
  v(1:n,1:n) = 0.0D+00
  do i = 1, n
    v(i,i) = 1.0D+00
  end do

  d(1) = 0.0D+00
  qd0 = 0.0D+00
  q0(1:n) = x(1:n)
  q1(1:n) = x(1:n)

  if ( 0 < prin ) then
    call print2 ( n, x, prin, fx, nf, nl )
  end if
!
!  The main loop starts here.
!
  do

    sf = d(1)
    d(1) = 0.0D+00
    s = 0.0D+00
!
!  Minimize along the first direction v(*,1).
!
    nits = 2
    value = fx

    call minny ( n, 1, nits, d(1), s, value, .false., f, x, t, h, v, q0, q1 )

    if ( s <= 0.0D+00 ) then
      v(1:n,1) = - v(1:n,1)
    end if

    if ( sf <= 0.9D+00 * d(1) .or. d(1) <= 0.9D+00 * sf ) then
      d(2:n) = 0.0D+00
    end if
!
!  The inner loop starts here.
!
    do k = 2, n

      y(1:n) = x(1:n)

      sf = fx

      if ( 0 < kt ) then
        illc = .true.
      end if

10    continue

      kl = k
      df = 0.0D+00
!
!  A random step follows (to avoid resolution valleys).
!  PRAXIS assumes that the random number generator returns a random 
!  number uniformly distributed in (0,1).
!
      if ( illc ) then

        do i = 1, n
          call random_number ( harvest = r )
          s = ( 0.1D+00 * ldt + t2 * 10.0D+00**kt ) * ( r - 0.5D+00 )
          z(i) = s
          x(1:n) = x(1:n) + s * v(1:n,i)
        end do

        fx = f(x,n)
        nf = nf + 1

      end if
!
!  Minimize along the "non-conjugate" directions V(*,K),...,V(*,N).
!
      do k2 = k, n

        sl = fx
        s = 0.0D+00
        nits = 2
        value = fx

        call minny ( n, k2, nits, d(k2), s, value, .false., f, x, t, h, v, &
          q0, q1 )

        if ( illc ) then
          s = d(k2) * ( ( s + z(k2) )**2 )
        else
          s = sl - fx
        end if

        if ( df <= s ) then
          df = s
          kl = k2
        end if

      end do
!
!  If there was not much improvement on the first try, set
!  ILLC = true and start the inner loop again.
!
      if ( .not. illc ) then
        if ( df < abs ( 100.0D+00 * machep * fx ) ) then
          illc = .true.
          go to 10
        end if
      end if

      if ( k == 2 .and. 1 < prin ) then
        call r8vec_print ( n, d, '  The second difference array' )
      end if
!
!  Minimize along the "conjugate" directions V(*,1),...,V(*,K-1).
!
      do k2 = 1, k-1

        s = 0.0D+00
        nits = 2
        value = fx

        call minny ( n, k2, nits, d(k2), s, value, .false., f, x, t, &
          h, v, q0, q1 )

      end do

      f1 = fx
      fx = sf
      lds = 0

      do i = 1, n
        sl = x(i)
        x(i) = y(i)
        sl = sl - y(i)
        y(i) = sl
        lds = lds + sl**2
      end do

      lds = sqrt ( lds )
!
!  Discard direction V(*,kl).
!
!  If no random step was taken, V(*,KL) is the "non-conjugate"
!  direction along which the greatest improvement was made.
!
      if ( small < lds ) then

        do ii = 1, kl - k
          i = kl - ii
          v(1:n,i+1) = v(1:n,i)
          d(i+1) = d(i)
        end do

        d(k) = 0.0D+00
        v(1:n,k) = y(1:n) / lds
!
!  Minimize along the new "conjugate" direction V(*,k), which is
!  the normalized vector:  (new x) - (old x).
!
        nits = 4
        value = f1

        call minny ( n, k, nits, d(k), lds, value, .true., f, x, t, &
          h, v, q0, q1 )

        if ( lds <= 0.0D+00 ) then
          lds = - lds
          v(1:n,k) = - v(1:n,k)
        end if

      end if

      ldt = ldfac * ldt
      ldt = max ( ldt, lds )

      if ( 0 < prin ) then
        call print2 ( n, x, prin, fx, nf, nl )
      end if

      t2 = m2 * sqrt ( sum ( x(1:n)**2 ) ) + t
!
!  See whether the length of the step taken since starting the
!  inner loop exceeds half the tolerance.
!
      if ( 0.5D+00 * t2 < ldt ) then
        kt = -1
      end if

      kt = kt + 1

      if ( ktm < kt ) then
        go to 20
      end if

    end do
!
!  The inner loop ends here.
!
!  Try quadratic extrapolation in case we are in a curved valley.
!
    call quad ( n, f, x, t, h, v, q0, q1 )

    d(1:n) = 1.0D+00 / sqrt ( d(1:n) )

    dn = maxval ( d(1:n) )

    if ( 3 < prin ) then
      call r8mat_print ( n, n, v, '  The new direction vectors' )
    end if

    do j = 1, n
      v(1:n,j) = ( d(j) / dn ) * v(1:n,j)
    end do
!
!  Scale the axes to try to reduce the condition number.
!
    if ( 1.0D+00 < scbd ) then

      do i = 1, n
        z(i) = sqrt ( sum ( v(i,1:n)**2 ) )
        z(i) = max ( z(i), m4 )
      end do

      s = minval ( z(1:n) )

      do i = 1, n

        sl = s / z(i)
        z(i) = 1.0D+00 / sl

        if ( scbd < z(i) ) then
          sl = 1.0D+00 / scbd
          z(i) = scbd
        end if

        v(i,1:n) = sl * v(i,1:n)

      end do

    end if
!
!  Calculate a new set of orthogonal directions before repeating
!  the main loop.
!
!  Transpose V for MINFIT:
!
    v(1:n,1:n) = transpose ( v(1:n,1:n) )
!
!  Call MINFIT to find the singular value decomposition of V.
!
!  This gives the principal values and principal directions of the
!  approximating quadratic form without squaring the condition number.
!
    call minfit ( n, n, vsmall, v, d )
!
!  Unscale the axes.
!
    if ( 1.0D+00 < scbd ) then

      do i = 1, n
        v(i,1:n) = z(i) * v(i,1:n)
      end do

      do i = 1, n

        s = sqrt ( sum ( v(1:n,i)**2 ) )

        d(i) = s * d(i)
        v(1:n,i) = v(1:n,i) / s

      end do

    end if

    do i = 1, n

      dni = dn * d(i)

      if ( large < dni ) then
        d(i) = vsmall
      else if ( dni < small ) then
        d(i) = vlarge
      else
        d(i) = 1.0D+00 / dni**2
      end if

    end do
!
!  Sort the eigenvalues and eigenvectors.
!
    call sort ( n, n, d, v )
!
!  Determine the smallest eigenvalue.
!
    dmin = max ( d(n), small )
!
!  The ratio of the smallest to largest eigenvalue determines whether
!  the system is ill conditioned.
!
    if ( dmin < m2 * d(1) ) then
      illc = .true.
    else
      illc = .false.
    end if

    if ( 1 < prin ) then

      if ( 1.0D+00 < scbd ) then
        call r8vec_print ( n, z, '  The scale factors' )
      end if 

      call r8vec_print ( n, d, '  Principal values of the quadratic form' )

    end if

    if ( 3 < prin ) then
      call r8mat_print ( n, n, v, '  The principal axes:' )
    end if
!
!  The main loop ends here.
!
  end do

20 continue

  if ( 0 < prin ) then
    call r8vec_print ( n, x, '  X:' )
  end if

  praxis = fx

  return
end
subroutine print2 ( n, x, prin, fx, nf, nl )

!*****************************************************************************80
!
!! PRINT2 prints certain data about the progress of the iteration.
!
!  Modified:
!
!    22 May 2006
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the current estimate of the minimizer.
!
!    Input, integer ( kind = 4 ) PRIN, the user-specifed print level.
!    0, nothing is printed.
!    1, F is printed after every n+1 or n+2 linear minimizations.  
!       final X is printed, but intermediate X is printed only 
!       if N is at most 4.
!    2, the scale factors and the principal values of the approximating 
!       quadratic form are also printed.
!    3, X is also printed after every few linear minimizations.
!    4, the principal vectors of the approximating quadratic form are 
!       also printed.
!
!    Input, real ( kind = 8 ) FX, the smallest value of F(X) found so far.
!
!    Input, integer ( kind = 4 ) NF, the number of function evaluations.
!
!    Input, integer ( kind = 4 ) NL, the number of linear searches.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) prin
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Linear searches      ', nl
  write ( *, '(a,i8)' ) '  Function evaluations ', nf 

  if ( n <= 4 .or. 2 < prin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'X:'
    write ( *, '(5g14.6)' ) x(1:n)
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) 'FX: ', fx

  return
end
subroutine quad ( n, f, x, t, h, v, q0, q1 )

!*****************************************************************************80
!
!! QUAD seeks to minimize the scalar function F along a particular curve.
!
!  Discussion:
!
!    The minimizer to be sought is required to lie on a curve defined
!    by Q0, Q1 and X.
!
!  Modified:
!
!    22 May 2006
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, external real ( kind = 8 ) F, is the name of the function to 
!    be minimized.  The function should have the form 
!      function f(x,n)
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!    and accepts X and N as input, returning in F the function value.
!
!    Input, real ( kind = 8 ) X(N), ?
!
!    Input, real ( kind = 8 ) T, ?
!
!    Input, real ( kind = 8 ) H, ?
!
!    Input, real ( kind = 8 ) V(N,N), the matrix of search directions.
!
!    Input, real ( kind = 8 ) Q0(N), Q1(N), two auxiliary points used to define
!    a curve through X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dmin
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  real ( kind = 8 ) l
  real ( kind = 8 ) ldt
  real ( kind = 8 ) machep
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nits
  integer ( kind = 4 ) nl
  real ( kind = 8 ) q0(n)
  real ( kind = 8 ) q1(n)
  real ( kind = 8 ) qa
  real ( kind = 8 ) qb
  real ( kind = 8 ) qc
  real ( kind = 8 ) qd0
  real ( kind = 8 ) qd1
  real ( kind = 8 ) qf1
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) value
  real ( kind = 8 ) x(n)

  common /global/ fx, ldt, dmin, nf, nl
  common /q/ qa, qb, qc, qd0, qd1, qf1

  machep = epsilon ( machep )

  temp = fx
  fx   = qf1
  qf1  = temp

  call r8vec_swap ( n, x, q1 )

  qd1 = sqrt ( sum ( ( x(1:n) - q1(1:n) )**2 ) )

  l = qd1
  s = 0.0D+00

  if ( qd0 <= 0.0D+00 .or. qd1 <= 0.0D+00 .or. nl < 3 * n**2 ) then

    fx = qf1
    qa = 0.0D+00
    qb = 0.0D+00
    qc = 1.0D+00

  else

    nits = 2
    value = qf1

    call minny ( n, 0, nits, s, l, value, .true., f, x, t, &
      h, v, q0, q1 )

    qa = ( l * ( l - qd1 ) ) / ( qd0 * ( qd0 + qd1 ) )
    qb = ( ( l + qd0 ) * ( qd1 - l ) ) / ( qd0 * qd1 )
    qc = ( l * ( l + qd0 ) ) / ( qd1 * ( qd0 + qd1 ) )

  end if

  qd0 = qd1

  do i = 1, n
    s = q0(i)
    q0(i) = x(i)
    x(i) = ( qa * s + qb * x(i) ) + qc * q1(i)
  end do

  return
end
function r8_hypot ( x, y )

!*****************************************************************************80
!
!! R8_HYPOT returns the value of sqrt ( X^2 + Y^2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the arguments.
!
!    Output, real ( kind = 8 ) R8_HYPOT, the value of sqrt ( X^2 + Y^2 ).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) r8_hypot
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( abs ( x ) < abs ( y ) ) then
    a = abs ( y )
    b = abs ( x )
  else
    a = abs ( x )
    b = abs ( y )
  end if
!
!  A contains the larger value.
!
  if ( a == 0.0D+00 ) then
    c = 0.0D+00
  else
    c = a * sqrt ( 1.0D+00 + ( b / a )**2 )
  end if

  r8_hypot = c

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
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
!    An R8MAT is an array of R8 values.
!
!  Modified:
!
!    26 March 2005
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

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
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
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end
subroutine r8vec_swap ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_SWAP swaps the entries of two R8VECs.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Modified:
!
!    04 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the arrays.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N), the vectors to swap.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a3(n)

  a3(1:n) = a1(1:n)
  a1(1:n) = a2(1:n)
  a2(1:n) = a3(1:n)

  return
end
subroutine sort ( m, n, d, v ) 

!*****************************************************************************80
!
!! SORT sorts a vector D and adjusts the corresponding columns of a matrix V.
!
!  Discussion:
!
!    A simple bubble sort is used on D.
!
!    In our application, D contains eigenvalues, and the columns of V are
!    the corresponding eigenvectors.
!
!  Modified:
!
!    25 February 2002
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the row dimension of V, which must be 
!    at least N.
!
!    Input, integer ( kind = 4 ) N, the length of D, and the order of V.
!
!    Input/output, real ( kind = 8 ) D(N), the vector to be sorted.  
!    On output, the entries of D are in descending order.
!
!    Input/output, real ( kind = 8 ) V(M,N), an N by N array to be adjusted 
!    as D is sorted.  In particular, if the value that was in D(I) on input is
!    moved to D(J) on output, then the input column V(*,I) is moved to
!    the output column V(*,J).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k(1)
  real ( kind = 8 ) temp
  real ( kind = 8 ) v(m,n)

  do i = 1, n - 1
!
!  Find K, the index of the largest entry in D(I:N).
!  MAXLOC apparently requires its output to be an array.
!
    k = maxloc ( d(i:n) )
!
!  If I < K, swap D(K) and D(I), and columns K and I of V.
!
    if ( i < k(1) ) then

      temp    = d(i)
      d(i)    = d(k(1))
      d(k(1)) = temp

      call r8vec_swap ( n, v(1:n,i), v(1:n,k(1)) )

    end if

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
