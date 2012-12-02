subroutine cascade ( n, t_length, t, c_length, c, s )

!*****************************************************************************80
!
!! CASCADE carries out the cascade algorithm.
!
!  Discussion:
!
!    The value of T3 computed by
!
!      call cascade ( 3, t_length, t0, c_length, c, t3 )
!
!    will be the same if computed in three steps by:
!
!      call cascade ( 1, t_length, t0, c_length, c, t1 );
!      call cascade ( 1, t_length, t1, c_length, c, t2 );
!      call cascade ( 1, t_length, t2, c_length, c, t3 );
!
!    If C represents a vector of Daubechies filter coefficients, then
!
!      call cascade ( 5, c_length, c, c_length, c, c5 );
!
!    computes an approximation to the corresponding scaling function, and
!
!      call r8vec_conjugate ( c_length, c, w )
!      call cascade ( 5, c_length, w, c_length, c, w5 );
!
!    computes an approximation to the corresponding wavelet.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of iterations to carry out.
!    0 <= N.
!
!    Input, integer ( kind = 4 ) T_LENGTH, the length of T.
!
!    Input, real ( kind = 8 ) T(T_LENGTH), the initial value of the quantity, 
!    or the value of the quantity at the integers 0 through P-1.
!
!    Input, integer ( kind = 4 ) C_LENGTH, the number of transform coefficients.
!
!    Input, real ( kind = 8 ) C(C_LENGTH), the transform coefficients.
!
!    Output, real ( kind = 8 ) S(2^N * T_LENGTH + (2^N-1)*C_LENGTH - 2*(2^N-1)),
!    the values of  the function.
!
  implicit none

  integer ( kind = 4 ) c_length
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t_length

  real ( kind = 8 ) c(c_length)
  integer ( kind = 4 ) i
  real ( kind = 8 ) s(*)
  integer ( kind = 4 ) s_length
  real ( kind = 8 ) t(t_length)
  integer ( kind = 4 ) x_length
  real ( kind = 8 ), allocatable :: x(:)

  s_length = t_length

  s(1:t_length) = t(1:t_length)

  do i = 1, n

    x_length = s_length * 2 - 1

    allocate ( x(1:x_length) )

    x(1:x_length  :2) = s(1:s_length)
    x(2:x_length-1:2) = 0.0D+00

    call r8vec_convolution ( x_length, x, c_length, c, s )

    deallocate ( x )

    s_length = x_length + c_length - 1

  end do

  return
end
subroutine daub_coefficients ( n, c )

!*****************************************************************************80
!
!! DAUB_COEFFICIENTS returns a set of Daubechies coefficients.
!
!  Discussion:
!
!    Often, the uses to which these coefficients are applied require that they
!    be rescaled, by being multiplied by sqrt ( 2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the coefficient set.
!    2 <= N <= 20, and N must be even.
!
!    Output, real ( kind = 8 ) C(N), the coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)

  if ( n == 2 ) then

    c = (/ &
      7.071067811865475D-01, &
      7.071067811865475D-01 /)

  else if ( n == 4 ) then

    c = (/ &
      0.4829629131445341D+00, &
      0.8365163037378079D+00, &
      0.2241438680420133D+00, &
    - 0.1294095225512603D+00 /)

  else if ( n == 6 ) then

    c = (/ &
      0.3326705529500826D+00, &
      0.8068915093110925D+00, &
      0.4598775021184915D+00, &
    - 0.1350110200102545D+00, &
    - 0.8544127388202666D-01, &
      0.3522629188570953D-01 /)

  else if ( n == 8 ) then

    c = (/ &
      0.2303778133088965D+00, &
      0.7148465705529156D+00, &
      0.6308807679298589D+00, &
     -0.2798376941685985D-01, &
     -0.1870348117190930D+00, &
      0.3084138183556076D-01, &
      0.3288301166688519D-01, &
     -0.1059740178506903D-01 /)

  else if ( n == 10 ) then

    c = (/ &
      0.1601023979741929D+00, &
      0.6038292697971896D+00, &
      0.7243085284377729D+00, &
      0.1384281459013207D+00, &
     -0.2422948870663820D+00, &
     -0.3224486958463837D-01, &
      0.7757149384004571D-01, &
     -0.6241490212798274D-02, &
     -0.1258075199908199D-01, &
      0.3335725285473771D-02 /);

  else if ( n == 12 ) then

    c = (/ &
      0.1115407433501094D+00, &
      0.4946238903984530D+00, &
      0.7511339080210953D+00, &
      0.3152503517091976D+00, &
     -0.2262646939654398D+00, &
     -0.1297668675672619D+00, &
      0.9750160558732304D-01, &
      0.2752286553030572D-01, &
     -0.3158203931748602D-01, &
      0.5538422011614961D-03, &
      0.4777257510945510D-02, &
     -0.1077301085308479D-02 /)

  else if ( n == 14 ) then

    c = (/ &
       7.785205408500917D-02, &
       3.965393194819173D-01, &
       7.291320908462351D-01, &
       4.697822874051931D-01, &
      -1.439060039285649D-01, &
      -2.240361849938749D-01, &
       7.130921926683026D-02, &
       8.061260915108307D-02, &
      -3.802993693501441D-02, &
      -1.657454163066688D-02, &
       1.255099855609984D-02, &
       4.295779729213665D-04, &
      -1.801640704047490D-03, &
       3.537137999745202D-04 /)

  else if ( n == 16 ) then

    c = (/ &
       5.441584224310400D-02, &
       3.128715909142999D-01, &
       6.756307362972898D-01, &
       5.853546836542067D-01, &
      -1.582910525634930D-02, &
      -2.840155429615469D-01, &
       4.724845739132827D-04, &
       1.287474266204784D-01, &
      -1.736930100180754D-02, &
      -4.408825393079475D-02, &
       1.398102791739828D-02, &
       8.746094047405776D-03, &
      -4.870352993451574D-03, &
      -3.917403733769470D-04, &
       6.754494064505693D-04, &
      -1.174767841247695D-04 /)

  else if ( n == 18 ) then

    c = (/ &
       3.807794736387834D-02, &
       2.438346746125903D-01, &
       6.048231236901111D-01, &
       6.572880780513005D-01, &
       1.331973858250075D-01, &
      -2.932737832791749D-01, &
      -9.684078322297646D-02, &
       1.485407493381063D-01, &
       3.072568147933337D-02, &
      -6.763282906132997D-02, &
       2.509471148314519D-04, &
       2.236166212367909D-02, &
      -4.723204757751397D-03, &
      -4.281503682463429D-03, &
       1.847646883056226D-03, &
       2.303857635231959D-04, &
      -2.519631889427101D-04, &
       3.934732031627159D-05 /)

  else if ( n == 20 ) then

    c = (/ &
       2.667005790055555D-02, &
       1.881768000776914D-01, &
       5.272011889317255D-01, &
       6.884590394536035D-01, &
       2.811723436605774D-01, &
      -2.498464243273153D-01, &
      -1.959462743773770D-01, &
       1.273693403357932D-01, &
       9.305736460357235D-02, &
      -7.139414716639708D-02, &
      -2.945753682187581D-02, &
       3.321267405934100D-02, &
       3.606553566956169D-03, &
      -1.073317548333057D-02, &
       1.395351747052901D-03, &
       1.992405295185056D-03, &
      -6.858566949597116D-04, &
      -1.164668551292854D-04, &
       9.358867032006959D-05, &
      -1.326420289452124D-05 /)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DAUB_COEFFICIENTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Value of N = ', n
    write ( *, '(a)' ) '  Legal values are 2, 4, 6, 8, 10, 12, 14, 16, 18, 20.'
    stop

  end if

  return
end
subroutine daub2_matrix ( n, a )

!*****************************************************************************80
!
!! DAUB2_MATRIX returns the DAUB2 matrix.
!
!  Discussion:
!
!    The DAUB2 matrix is the Daubechies wavelet transformation matrix 
!    with 2 coefficients.
!
!    The DAUB2 matrix is also known as the Haar matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2 and a multiple of 2.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 1

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ), dimension(0:p) :: c = (/ &
    7.071067811865475D-01, &
    7.071067811865475D-01 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) m

  if ( n < 2 .or. mod ( n, 2 ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DAUB2_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  Order N must be at least 2 and a multiple of 2.'
    stop
  end if

  a(1:n,1:n) = 0.0D+00

  m = n / 2

  do i = 1, n - 1, 2

    a(i,i)     =   c(0)
    a(i,i+1)   =   c(1)

    a(i+1,i)   =   c(1)
    a(i+1,i+1) = - c(0)

  end do

  return
end
recursive function daub2_scale ( n, x ) result ( y )

!*****************************************************************************80
!
!! DAUB2_SCALE recursively evaluates the DAUB2 scaling function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the recursion level.
!
!    Input, real ( kind = 8 ) X, the point at which the function is to 
!    be evaluated.
!
!    Output, real ( kind = 8 ) Y, the estimated value of the function.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension(0:1) :: c = (/ &
    7.071067811865475D-01, &
    7.071067811865475D-01 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( 0 < n ) then
    y = sqrt ( 2.0D+00 ) * &
        ( c(0) * daub2_scale ( n - 1, 2.0D+00 * x           ) &
        + c(1) * daub2_scale ( n - 1, 2.0D+00 * x - 1.0D+00 ) )
  else
    if ( 0.0D+00 <= x .and. x < 1.0D+00 ) then
      y = 1.0D+00
    else
      y = 0.0D+00
    end if
  end if

  return
end
subroutine daub2_transform ( n, x, y )

!*****************************************************************************80
!
!! DAUB2_TRANSFORM computes the DAUB2 transform of a vector.
!
!  Discussion:
!
!    DAUB2 is better known as the Haar transform.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2.
!
!    Input, real ( kind = 8 ) X(N), the vector to be transformed. 
!
!    Output, real ( kind = 8 ) Y(N), the transformed vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 1

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
    7.071067811865475D-01, &
    7.071067811865475D-01 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  y(1:n) = x(1:n)
  z(1:n) = 0.0D+00

  m = n

  do while ( 2 <= m )
  
    m = m / 2

    do i = 1, m
      z(i)   = c(0) * ( y(2*i-1) + y(2*i) )
      z(i+m) = c(1) * ( y(2*i-1) - y(2*i) )
    end do

    y(1:2*m) = z(1:2*m)

  end do

  return
end
subroutine daub2_transform_inverse ( n, y, x )

!*****************************************************************************80
!
!! DAUB2_TRANSFORM_INVERSE inverts the DAUB2 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2.
!
!    Input, real ( kind = 8 ) Y(N), the transformed vector.
!
!    Output, real ( kind = 8 ) X(N), the original vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 1

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
    7.071067811865475D-01, &
    7.071067811865475D-01 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  x(1:n) = y(1:n)
  z(1:n) = 0.0D+00
 
  m = 1

  do while ( m * 2 <= n )

    do i = 1, m
      z(2*i-1) = c(0) * ( x(i) + x(i+m) )
      z(2*i)   = c(1) * ( x(i) - x(i+m) )
    end do

    x(1:2*m) = z(1:2*m)

    m = m * 2

  end do

  return
end
subroutine daub4_matrix ( n, a )

!*****************************************************************************80
!
!! DAUB4_MATRIX returns the DAUB4 matrix.
!
!  Discussion:
!
!    The DAUB4 matrix is the Daubechies wavelet transformation matrix 
!    with 4 coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 4 and a multiple of 2.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 3

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     0.4829629131445341D+00, &
     0.8365163037378079D+00, &
     0.2241438680420133D+00, &
    -0.1294095225512603D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) m

  if ( n < 4 .or. mod ( n, 2 ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DAUB4_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  Order N must be at least 4 and a multiple of 2.'
    stop
  end if

  a(1:n,1:n) = 0.0D+00

  do i = 1, n - 1, 2

    a(i,i)                  =   c(0)
    a(i,i+1)                =   c(1)
    a(i,i4_wrap(i+2,1,n))   =   c(2)
    a(i,i4_wrap(i+3,1,n))   =   c(3)

    a(i+1,i)                =   c(3)
    a(i+1,i+1)              = - c(2)
    a(i+1,i4_wrap(i+2,1,n)) =   c(1)
    a(i+1,i4_wrap(i+3,1,n)) = - c(0)

  end do

  return
end
recursive function daub4_scale ( n, x ) result ( y )

!*****************************************************************************80
!
!! DAUB4_SCALE recursively evaluates the DAUB4 scaling function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the recursion level.
!
!    Input, real ( kind = 8 ) X, the point at which the function is to 
!    be evaluated.
!
!    Output, real ( kind = 8 ) Y, the estimated value of the function.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension(0:3) :: c = (/ &
     0.4829629131445341D+00, &
     0.8365163037378079D+00, &
     0.2241438680420133D+00, &
    -0.1294095225512603D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( 0 < n ) then
    y = sqrt ( 2.0D+00 ) * &
        ( c(0) * daub4_scale ( n - 1, 2.0D+00 * x           ) &
        + c(1) * daub4_scale ( n - 1, 2.0D+00 * x - 1.0D+00 ) &
        + c(2) * daub4_scale ( n - 1, 2.0D+00 * x - 2.0D+00 ) &
        + c(3) * daub4_scale ( n - 1, 2.0D+00 * x - 3.0D+00 ) )
  else
    if ( 0.0D+00 <= x .and. x < 1.0D+00 ) then
      y = 1.0D+00
    else
      y = 0.0D+00
    end if
  end if

  return
end
subroutine daub4_transform ( n, x, y )

!*****************************************************************************80
!
!! DAUB4_TRANSFORM computes the DAUB4 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) X(N), the vector to be transformed. 
!
!    Output, real ( kind = 8 ) Y(N), the transformed vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 3

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     0.4829629131445341D+00, &
     0.8365163037378079D+00, &
     0.2241438680420133D+00, &
    -0.1294095225512603D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  y(1:n) = x(1:n)
  z(1:n) = 0.0D+00

  m = n

  do while ( 4 <= m )
  
    i = 1

    do j = 1, m - 1, 2

      j0 = i4_wrap ( j,     1, m )
      j1 = i4_wrap ( j + 1, 1, m )
      j2 = i4_wrap ( j + 2, 1, m )
      j3 = i4_wrap ( j + 3, 1, m )

      z(i)     = c(0) * y(j0) + c(1) * y(j1) &
               + c(2) * y(j2) + c(3) * y(j3)

      z(i+m/2) = c(3) * y(j0) - c(2) * y(j1) &
               + c(1) * y(j2) - c(0) * y(j3)

      i = i + 1

    end do

    y(1:m) = z(1:m)

    m = m / 2

  end do

  return
end
subroutine daub4_transform_inverse ( n, y, x )

!*****************************************************************************80
!
!! DAUB4_TRANSFORM_INVERSE inverts the DAUB4 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) Y(N), the transformed vector. 
!
!    Output, real ( kind = 8 ) X(N), the original vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 3

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     0.4829629131445341D+00, &
     0.8365163037378079D+00, &
     0.2241438680420133D+00, &
    -0.1294095225512603D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  x(1:n) = y(1:n)
  z(1:n) = 0.0D+00

  m = 4

  do while ( m <= n )
  
    j = 1

    do i = 0, m / 2 - 1
      
      i0 = i4_wrap ( i,                 1,         m / 2 )
      i2 = i4_wrap ( i + 1,             1,         m / 2 )

      i1 = i4_wrap ( i + m / 2,         m / 2 + 1, m     )
      i3 = i4_wrap ( i + m / 2 + 1,     m / 2 + 1, m     )

      z(j)   = c(2) * x(i0) + c(1) * x(i1) &
             + c(0) * x(i2) + c(3) * x(i3)

      z(j+1) = c(3) * x(i0) - c(0) * x(i1) &
             + c(1) * x(i2) - c(2) * x(i3)

      j = j + 2

    end do

    x(1:m) = z(1:m)

    m = m * 2

  end do

  return
end
subroutine daub6_matrix ( n, a )

!*****************************************************************************80
!
!! DAUB6_MATRIX returns the DAUB6 matrix.
!
!  Discussion:
!
!    The DAUB6 matrix is the Daubechies wavelet transformation matrix 
!    with 6 coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 6, and a multiple of 2.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     0.3326705529500826D+00, &
     0.8068915093110925D+00, &
     0.4598775021184915D+00, &
   - 0.1350110200102545D+00, &
   - 0.08544127388202666D+00, &
     0.03522629188570953D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) m

  if ( n < 6 .or. mod ( n, 2 ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DAUB6_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 6 and a multiple of 2.'
    stop
  end if

  a(1:n,1:n) = 0.0D+00

  do i = 1, n - 1, 2

    a(i,i)                  =   c(0)
    a(i,i+1)                =   c(1)
    a(i,i4_wrap(i+2,1,n))   =   c(2)
    a(i,i4_wrap(i+3,1,n))   =   c(3)
    a(i,i4_wrap(i+4,1,n))   =   c(4)
    a(i,i4_wrap(i+5,1,n))   =   c(5)

    a(i+1,i)                =   c(5)
    a(i+1,i+1)              = - c(4)
    a(i+1,i4_wrap(i+2,1,n)) =   c(3)
    a(i+1,i4_wrap(i+3,1,n)) = - c(2)
    a(i+1,i4_wrap(i+4,1,n)) =   c(1)
    a(i+1,i4_wrap(i+5,1,n)) = - c(0)

  end do

  return
end
recursive function daub6_scale ( n, x ) result ( y )

!*****************************************************************************80
!
!! DAUB6_SCALE recursively evaluates the DAUB6 scaling function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the recursion level.
!
!    Input, real ( kind = 8 ) X, the point at which the function is to 
!    be evaluated.
!
!    Output, real ( kind = 8 ) Y, the estimated value of the function.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension(0:5) :: c = (/ &
     0.3326705529500826D+00, &
     0.8068915093110925D+00, &
     0.4598775021184915D+00, &
   - 0.1350110200102545D+00, &
   - 0.08544127388202666D+00, &
     0.03522629188570953D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( 0 < n ) then
    y = sqrt ( 2.0D+00 ) * &
        ( c(0) * daub6_scale ( n - 1, 2.0D+00 * x           ) &
        + c(1) * daub6_scale ( n - 1, 2.0D+00 * x - 1.0D+00 ) &
        + c(2) * daub6_scale ( n - 1, 2.0D+00 * x - 2.0D+00 ) &
        + c(3) * daub6_scale ( n - 1, 2.0D+00 * x - 3.0D+00 ) &
        + c(4) * daub6_scale ( n - 1, 2.0D+00 * x - 4.0D+00 ) &
        + c(5) * daub6_scale ( n - 1, 2.0D+00 * x - 5.0D+00 ) )
  else
    if ( 0.0D+00 <= x .and. x < 1.0D+00 ) then
      y = 1.0D+00
    else
      y = 0.0D+00
    end if
  end if

  return
end
subroutine daub6_transform ( n, x, y )

!*****************************************************************************80
!
!! DAUB6_TRANSFORM computes the DAUB6 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) X(N), the vector to be transformed. 
!
!    Output, real ( kind = 8 ) Y(N), the transformed vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 5

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     0.3326705529500826D+00, &
     0.8068915093110925D+00, &
     0.4598775021184915D+00, &
   - 0.1350110200102545D+00, &
   - 0.08544127388202666D+00, &
     0.03522629188570953D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  y(1:n) = x(1:n)

  m = n
  q = ( p - 1 ) / 2

  do while ( 4 <= m )
  
    i = 1
    z(1:m) = 0.0D+00

    do j = 1, m - 1, 2

      do k = 0, p - 1, 2
        j0 = i4_wrap ( j + k,     1, m )
        j1 = i4_wrap ( j + k + 1, 1, m )
        z(i)     = z(i)     + c(  k) * y(j0) + c(  k+1) * y(j1)
        z(i+m/2) = z(i+m/2) + c(p-k) * y(j0) - c(p-k-1) * y(j1)
      end do

      i = i + 1

    end do

    y(1:m) = z(1:m)

    m = m / 2

  end do

  return
end
subroutine daub6_transform_inverse ( n, y, x )

!*****************************************************************************80
!
!! DAUB6_TRANSFORM_INVERSE inverts the DAUB6 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) Y(N), the transformed vector. 
!
!    Output, real ( kind = 8 ) X(N), the original vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 5

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     0.3326705529500826D+00, &
     0.8068915093110925D+00, &
     0.4598775021184915D+00, &
   - 0.1350110200102545D+00, &
   - 0.08544127388202666D+00, &
     0.03522629188570953D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  x(1:n) = y(1:n)

  m = 4
  q = ( p - 1 ) / 2

  do while ( m <= n )
  
    z(1:m) = 0.0D+00

    j = 1

    do i = - q + 1, m / 2 - q
      
      do k = 0, p - 1, 2
        i0 = i4_wrap ( i         + k / 2,     1,         m / 2 )
        i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2 + 1, m     )
        z(j)   = z(j)   + c(p-k-1) * x(i0) + c(k+1) * x(i1)
        z(j+1) = z(j+1) + c(p-k)   * x(i0) - c(k)   * x(i1)
      end do

      j = j + 2

    end do

    x(1:m) = z(1:m)

    m = m * 2

  end do

  return
end
subroutine daub8_matrix ( n, a )

!*****************************************************************************80
!
!! DAUB8_MATRIX returns the DAUB8 matrix.
!
!  Discussion:
!
!    The DAUB8 matrix is the Daubechies wavelet transformation matrix 
!    with 8 coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gilbert Strang, Truong Nguyen,
!    Wavelets and Filter Banks,
!    Wellesley-Cambridge Press, 1997,
!    ISBN: 0-9614088-7-1,
!    LC: TK7872.F5S79 / QA403.3.S87
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 8, and a multiple of 2.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ), dimension ( 0 : 7 ) :: c = (/ &
    0.2303778133088964D+00, &
    0.7148465705529154D+00, &
    0.6308807679298587D+00, &
   -0.0279837694168599D+00, &
   -0.1870348117190931D+00, &
    0.0308413818355607D+00, &
    0.0328830116668852D+00, &
   -0.0105974017850690D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) m

  if ( n < 8 .or. mod ( n, 2 ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DAUB8_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 8 and a multiple of 2.'
    stop
  end if

  a(1:n,1:n) = 0.0D+00

  do i = 1, n - 1, 2

    a(i,i)                  =   c(0)
    a(i,i+1)                =   c(1)
    a(i,i4_wrap(i+2,1,n))   =   c(2)
    a(i,i4_wrap(i+3,1,n))   =   c(3)
    a(i,i4_wrap(i+4,1,n))   =   c(4)
    a(i,i4_wrap(i+5,1,n))   =   c(5)
    a(i,i4_wrap(i+6,1,n))   =   c(6)
    a(i,i4_wrap(i+7,1,n))   =   c(7)

    a(i+1,i)                =   c(7)
    a(i+1,i+1)              = - c(6)
    a(i+1,i4_wrap(i+2,1,n)) =   c(5)
    a(i+1,i4_wrap(i+3,1,n)) = - c(4)
    a(i+1,i4_wrap(i+4,1,n)) =   c(3)
    a(i+1,i4_wrap(i+5,1,n)) = - c(2)
    a(i+1,i4_wrap(i+6,1,n)) =   c(1)
    a(i+1,i4_wrap(i+7,1,n)) = - c(0)

  end do

  return
end
recursive function daub8_scale ( n, x ) result ( y )

!*****************************************************************************80
!
!! DAUB8_SCALE recursively evaluates the DAUB8 scaling function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the recursion level.
!
!    Input, real ( kind = 8 ) X, the point at which the function is to 
!    be evaluated.
!
!    Output, real ( kind = 8 ) Y, the estimated value of the function.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension(0:7) :: c = (/ &
    0.2303778133088964D+00, &
    0.7148465705529154D+00, &
    0.6308807679298587D+00, &
   -0.0279837694168599D+00, &
   -0.1870348117190931D+00, &
    0.0308413818355607D+00, &
    0.0328830116668852D+00, &
   -0.0105974017850690D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( 0 < n ) then
    y = sqrt ( 2.0D+00 ) * &
        ( c(0) * daub8_scale ( n - 1, 2.0D+00 * x           ) &
        + c(1) * daub8_scale ( n - 1, 2.0D+00 * x - 1.0D+00 ) &
        + c(2) * daub8_scale ( n - 1, 2.0D+00 * x - 2.0D+00 ) &
        + c(3) * daub8_scale ( n - 1, 2.0D+00 * x - 3.0D+00 ) &
        + c(4) * daub8_scale ( n - 1, 2.0D+00 * x - 4.0D+00 ) &
        + c(5) * daub8_scale ( n - 1, 2.0D+00 * x - 5.0D+00 ) &
        + c(6) * daub8_scale ( n - 1, 2.0D+00 * x - 6.0D+00 ) &
        + c(7) * daub8_scale ( n - 1, 2.0D+00 * x - 7.0D+00 ) )
  else
    if ( 0.0D+00 <= x .and. x < 1.0D+00 ) then
      y = 1.0D+00
    else
      y = 0.0D+00
    end if
  end if

  return
end
subroutine daub8_transform ( n, x, y )

!*****************************************************************************80
!
!! DAUB8_TRANSFORM computes the DAUB8 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) X(N), the vector to be transformed. 
!
!    Output, real ( kind = 8 ) Y(N), the transformed vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 7

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     0.2303778133088964D+00, &
     0.7148465705529154D+00, &
     0.6308807679298587D+00, &
    -0.02798376941685985D+00, &
    -0.1870348117190931D+00, &
     0.03084138183556076D+00, &
     0.03288301166688519D+00, &
    -0.01059740178506903D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  y(1:n) = x(1:n)

  m = n
  q = ( p - 1 ) / 2

  do while ( 4 <= m )
  
    i = 1
    z(1:m) = 0.0D+00

    do j = 1, m - 1, 2

      do k = 0, p - 1, 2
        j0 = i4_wrap ( j + k,     1, m )
        j1 = i4_wrap ( j + k + 1, 1, m )
        z(i)     = z(i)     + c(  k) * y(j0) + c(  k+1) * y(j1)
        z(i+m/2) = z(i+m/2) + c(p-k) * y(j0) - c(p-k-1) * y(j1)
      end do

      i = i + 1

    end do

    y(1:m) = z(1:m)

    m = m / 2

  end do

  return
end
subroutine daub8_transform_inverse ( n, y, x )

!*****************************************************************************80
!
!! DAUB8_TRANSFORM_INVERSE inverts the DAUB8 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) Y(N), the transformed vector. 
!
!    Output, real ( kind = 8 ) X(N), the original vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 7

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     0.2303778133088964D+00, &
     0.7148465705529154D+00, &
     0.6308807679298587D+00, &
    -0.02798376941685985D+00, &
    -0.1870348117190931D+00, &
     0.03084138183556076D+00, &
     0.03288301166688519D+00, &
    -0.01059740178506903D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  x(1:n) = y(1:n)

  m = 4
  q = ( p - 1 ) / 2

  do while ( m <= n )
  
    z(1:m) = 0.0D+00

    j = 1

    do i = - q + 1, m / 2 - q
      
      do k = 0, p - 1, 2
        i0 = i4_wrap ( i         + k / 2,     1,         m / 2 )
        i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2 + 1, m     )
        z(j)   = z(j)   + c(p-k-1) * x(i0) + c(k+1) * x(i1)
        z(j+1) = z(j+1) + c(p-k)   * x(i0) - c(k)   * x(i1)
      end do

      j = j + 2

    end do

    x(1:m) = z(1:m)

    m = m * 2

  end do

  return
end
subroutine daub10_matrix ( n, a )

!*****************************************************************************80
!
!! DAUB10_MATRIX returns the DAUB10 matrix.
!
!  Discussion:
!
!    The DAUB10 matrix is the Daubechies wavelet transformation matrix 
!    with 10 coefficients.
!
!    Note that in the reference, the coefficient 0.0775714938400459D+00 
!    is given incorrectly, with the "8" misrepresented as a "0".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gilbert Strang, Truong Nguyen,
!    Wavelets and Filter Banks,
!    Wellesley-Cambridge Press, 1997,
!    ISBN: 0-9614088-7-1,
!    LC: TK7872.F5S79 / QA403.3.S87
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 10, and a multiple of 2.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ), dimension ( 0 : 9 ) :: c = (/ &
    0.1601023979741929D+00, &
    0.6038292697971895D+00, &
    0.7243085284377726D+00, &
    0.1384281459013203D+00, &
   -0.2422948870663823D+00, &
   -0.0322448695846381D+00, &
    0.0775714938400459D+00, &
   -0.0062414902127983D+00, &
   -0.0125807519990820D+00, &
    0.0033357252854738D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) m

  if ( n < 10 .or. mod ( n, 2 ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DAUB10_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 10 and a multiple of 2.'
    stop
  end if

  a(1:n,1:n) = 0.0D+00

  do i = 1, n - 1, 2

    a(i,i)                  =   c(0)
    a(i,i+1)                =   c(1)
    a(i,i4_wrap(i+2,1,n))   =   c(2)
    a(i,i4_wrap(i+3,1,n))   =   c(3)
    a(i,i4_wrap(i+4,1,n))   =   c(4)
    a(i,i4_wrap(i+5,1,n))   =   c(5)
    a(i,i4_wrap(i+6,1,n))   =   c(6)
    a(i,i4_wrap(i+7,1,n))   =   c(7)
    a(i,i4_wrap(i+8,1,n))   =   c(8)
    a(i,i4_wrap(i+9,1,n))   =   c(9)

    a(i+1,i)                =   c(9)
    a(i+1,i+1)              = - c(8)
    a(i+1,i4_wrap(i+2,1,n)) =   c(7)
    a(i+1,i4_wrap(i+3,1,n)) = - c(6)
    a(i+1,i4_wrap(i+4,1,n)) =   c(5)
    a(i+1,i4_wrap(i+5,1,n)) = - c(4)
    a(i+1,i4_wrap(i+6,1,n)) =   c(3)
    a(i+1,i4_wrap(i+7,1,n)) = - c(2)
    a(i+1,i4_wrap(i+8,1,n)) =   c(1)
    a(i+1,i4_wrap(i+9,1,n)) = - c(0)

  end do

  return
end
recursive function daub10_scale ( n, x ) result ( y )

!*****************************************************************************80
!
!! DAUB10_SCALE recursively evaluates the DAUB10 scaling function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the recursion level.
!
!    Input, real ( kind = 8 ) X, the point at which the function is to 
!    be evaluated.
!
!    Output, real ( kind = 8 ) Y, the estimated value of the function.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension(0:9) :: c = (/ &
    0.1601023979741929D+00, &
    0.6038292697971895D+00, &
    0.7243085284377726D+00, &
    0.1384281459013203D+00, &
   -0.2422948870663823D+00, &
   -0.0322448695846381D+00, &
    0.0775714938400459D+00, &
   -0.0062414902127983D+00, &
   -0.0125807519990820D+00, &
    0.0033357252854738D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( 0 < n ) then
    y = sqrt ( 2.0D+00 ) * &
        ( c(0) * daub10_scale ( n - 1, 2.0D+00 * x           ) &
        + c(1) * daub10_scale ( n - 1, 2.0D+00 * x - 1.0D+00 ) &
        + c(2) * daub10_scale ( n - 1, 2.0D+00 * x - 2.0D+00 ) &
        + c(3) * daub10_scale ( n - 1, 2.0D+00 * x - 3.0D+00 ) &
        + c(4) * daub10_scale ( n - 1, 2.0D+00 * x - 4.0D+00 ) &
        + c(5) * daub10_scale ( n - 1, 2.0D+00 * x - 5.0D+00 ) &
        + c(6) * daub10_scale ( n - 1, 2.0D+00 * x - 6.0D+00 ) &
        + c(7) * daub10_scale ( n - 1, 2.0D+00 * x - 7.0D+00 ) &
        + c(8) * daub10_scale ( n - 1, 2.0D+00 * x - 8.0D+00 ) &
        + c(9) * daub10_scale ( n - 1, 2.0D+00 * x - 9.0D+00 ) )
  else
    if ( 0.0D+00 <= x .and. x < 1.0D+00 ) then
      y = 1.0D+00
    else
      y = 0.0D+00
    end if
  end if

  return
end
subroutine daub10_transform ( n, x, y )

!*****************************************************************************80
!
!! DAUB10_TRANSFORM computes the DAUB10 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) X(N), the vector to be transformed. 
!
!    Output, real ( kind = 8 ) Y(N), the transformed vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 9

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     1.601023979741929D-01, &
     6.038292697971896D-01, &
     7.243085284377729D-01, &
     1.384281459013207D-01, &
    -2.422948870663820D-01, &
    -3.224486958463837D-02, &
     7.757149384004571D-02, &
    -6.241490212798274D-03, &
    -1.258075199908199D-02, &
     3.335725285473771D-03 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  y(1:n) = x(1:n)

  m = n
  q = ( p - 1 ) / 2

  do while ( 4 <= m )
  
    i = 1
    z(1:m) = 0.0D+00

    do j = 1, m - 1, 2

      do k = 0, p - 1, 2
        j0 = i4_wrap ( j + k,     1, m )
        j1 = i4_wrap ( j + k + 1, 1, m )
        z(i)     = z(i)     + c(  k) * y(j0) + c(  k+1) * y(j1)
        z(i+m/2) = z(i+m/2) + c(p-k) * y(j0) - c(p-k-1) * y(j1)
      end do

      i = i + 1

    end do

    y(1:m) = z(1:m)

    m = m / 2

  end do

  return
end
subroutine daub10_transform_inverse ( n, y, x )

!*****************************************************************************80
!
!! DAUB10_TRANSFORM_INVERSE inverts the DAUB10 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) Y(N), the transformed vector. 
!
!    Output, real ( kind = 8 ) X(N), the original vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 9

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     1.601023979741929D-01, &
     6.038292697971896D-01, &
     7.243085284377729D-01, &
     1.384281459013207D-01, &
    -2.422948870663820D-01, &
    -3.224486958463837D-02, &
     7.757149384004571D-02, &
    -6.241490212798274D-03, &
    -1.258075199908199D-02, &
     3.335725285473771D-03 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  x(1:n) = y(1:n)

  m = 4
  q = ( p - 1 ) / 2

  do while ( m <= n )
  
    z(1:m) = 0.0D+00

    j = 1

    do i = - q + 1, m / 2 - q
      
      do k = 0, p - 1, 2
        i0 = i4_wrap ( i         + k / 2,     1,         m / 2 )
        i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2 + 1, m     )
        z(j)   = z(j)   + c(p-k-1) * x(i0) + c(k+1) * x(i1)
        z(j+1) = z(j+1) + c(p-k)   * x(i0) - c(k)   * x(i1)
      end do

      j = j + 2

    end do

    x(1:m) = z(1:m)

    m = m * 2

  end do

  return
end
subroutine daub12_matrix ( n, a )

!*****************************************************************************80
!
!! DAUB12_MATRIX returns the DAUB12 matrix.
!
!  Discussion:
!
!    The DAUB12 matrix is the Daubechies wavelet transformation matrix 
!    with 12 coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gilbert Strang, Truong Nguyen,
!    Wavelets and Filter Banks,
!    Wellesley-Cambridge Press, 1997,
!    ISBN: 0-9614088-7-1,
!    LC: TK7872.F5S79 / QA403.3.S87
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 12, and a multiple of 2.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 11

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ), dimension ( 0 : p ) :: c = (/ &
    0.1115407433501095D+00, &
    0.4946238903984533D+00, &
    0.7511339080210959D+00, &
    0.3152503517091982D+00, &
   -0.2262646939654400D+00, &
   -0.1297668675672625D+00, &
    0.0975016055873225D+00, &
    0.0275228655303053D+00, &
   -0.0315820393174862D+00, &
    0.0005538422011614D+00, &
    0.0047772575109455D+00, &
   -0.0010773010853085D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) q

  if ( n < 12 .or. mod ( n, 2 ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DAUB12_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 12 and a multiple of 2.'
    stop
  end if

  q = ( p - 1 ) / 2

  a(1:n,1:n) = 0.0D+00

  do i = 1, n - 1, 2

    a(i,i)                   =    c(0)
    a(i,i+1)                 =    c(1)
    a(i,i4_wrap(i+ 2,1,n))   =    c(2)
    a(i,i4_wrap(i+ 3,1,n))   =    c(3)
    a(i,i4_wrap(i+ 4,1,n))   =    c(4)
    a(i,i4_wrap(i+ 5,1,n))   =    c(5)
    a(i,i4_wrap(i+ 6,1,n))   =    c(6)
    a(i,i4_wrap(i+ 7,1,n))   =    c(7)
    a(i,i4_wrap(i+ 8,1,n))   =    c(8)
    a(i,i4_wrap(i+ 9,1,n))   =    c(9)
    a(i,i4_wrap(i+10,1,n))   =   c(10)
    a(i,i4_wrap(i+11,1,n))   =   c(11)

    a(i+1,i)                 =   c(11)
    a(i+1,i+1)               = - c(10)
    a(i+1,i4_wrap(i+ 2,1,n)) =    c(9)
    a(i+1,i4_wrap(i+ 3,1,n)) =  - c(8)
    a(i+1,i4_wrap(i+ 4,1,n)) =    c(7)
    a(i+1,i4_wrap(i+ 5,1,n)) =  - c(6)
    a(i+1,i4_wrap(i+ 6,1,n)) =    c(5)
    a(i+1,i4_wrap(i+ 7,1,n)) =  - c(4)
    a(i+1,i4_wrap(i+ 8,1,n)) =    c(3)
    a(i+1,i4_wrap(i+ 9,1,n)) =  - c(2)
    a(i+1,i4_wrap(i+10,1,n)) =    c(1)
    a(i+1,i4_wrap(i+11,1,n)) =  - c(0)

  end do

  return
end
subroutine daub12_transform ( n, x, y )

!*****************************************************************************80
!
!! DAUB12_TRANSFORM computes the DAUB12 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) X(N), the vector to be transformed. 
!
!    Output, real ( kind = 8 ) Y(N), the transformed vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 11

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
    0.1115407433501095D+00, &
    0.4946238903984533D+00, &
    0.7511339080210959D+00, &
    0.3152503517091982D+00, &
   -0.2262646939654400D+00, &
   -0.1297668675672625D+00, &
    0.0975016055873225D+00, &
    0.0275228655303053D+00, &
   -0.0315820393174862D+00, &
    0.0005538422011614D+00, &
    0.0047772575109455D+00, &
   -0.0010773010853085D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  y(1:n) = x(1:n)

  m = n
  q = ( p - 1 ) / 2

  do while ( 4 <= m )
  
    i = 1
    z(1:m) = 0.0D+00

    do j = 1, m - 1, 2

      do k = 0, p - 1, 2
        j0 = i4_wrap ( j + k,     1, m )
        j1 = i4_wrap ( j + k + 1, 1, m )
        z(i)     = z(i)     + c(  k) * y(j0) + c(  k+1) * y(j1)
        z(i+m/2) = z(i+m/2) + c(p-k) * y(j0) - c(p-k-1) * y(j1)
      end do

      i = i + 1

    end do

    y(1:m) = z(1:m)

    m = m / 2

  end do

  return
end
subroutine daub12_transform_inverse ( n, y, x )

!*****************************************************************************80
!
!! DAUB12_TRANSFORM_INVERSE inverts the DAUB12 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) Y(N), the transformed vector. 
!
!    Output, real ( kind = 8 ) X(N), the original vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 11

  real ( kind = 8 ), dimension ( 0:p ) :: c = (/ &
    0.1115407433501095D+00, &
    0.4946238903984533D+00, &
    0.7511339080210959D+00, &
    0.3152503517091982D+00, &
   -0.2262646939654400D+00, &
   -0.1297668675672625D+00, &
    0.0975016055873225D+00, &
    0.0275228655303053D+00, &
   -0.0315820393174862D+00, &
    0.0005538422011614D+00, &
    0.0047772575109455D+00, &
   -0.0010773010853085D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  x(1:n) = y(1:n)

  m = 4
  q = ( p - 1 ) / 2

  do while ( m <= n )
  
    z(1:m) = 0.0D+00

    j = 1

    do i = - q + 1, m / 2 - q
      
      do k = 0, p - 1, 2
        i0 = i4_wrap ( i         + k / 2,     1,         m / 2 )
        i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2 + 1, m     )
        z(j)   = z(j)   + c(p-k-1) * x(i0) + c(k+1) * x(i1)
        z(j+1) = z(j+1) + c(p-k)   * x(i0) - c(k)   * x(i1)
      end do

      j = j + 2

    end do

    x(1:m) = z(1:m)

    m = m * 2

  end do

  return
end
subroutine daub14_transform ( n, x, y )

!*****************************************************************************80
!
!! DAUB14_TRANSFORM computes the DAUB14 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) X(N), the vector to be transformed. 
!
!    Output, real ( kind = 8 ) Y(N), the transformed vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 13

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     7.785205408500917D-02, &
     3.965393194819173D-01, &
     7.291320908462351D-01, &
     4.697822874051931D-01, &
    -1.439060039285649D-01, &
    -2.240361849938749D-01, &
     7.130921926683026D-02, &
     8.061260915108307D-02, &
    -3.802993693501441D-02, &
    -1.657454163066688D-02, &
     1.255099855609984D-02, &
     4.295779729213665D-04, &
    -1.801640704047490D-03, &
     3.537137999745202D-04 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  y(1:n) = x(1:n)

  m = n
  q = ( p - 1 ) / 2

  do while ( 4 <= m )
  
    i = 1
    z(1:m) = 0.0D+00

    do j = 1, m - 1, 2

      do k = 0, p - 1, 2
        j0 = i4_wrap ( j + k,     1, m )
        j1 = i4_wrap ( j + k + 1, 1, m )
        z(i)     = z(i)     + c(  k) * y(j0) + c(  k+1) * y(j1)
        z(i+m/2) = z(i+m/2) + c(p-k) * y(j0) - c(p-k-1) * y(j1)
      end do

      i = i + 1

    end do

    y(1:m) = z(1:m)

    m = m / 2

  end do

  return
end
subroutine daub14_transform_inverse ( n, y, x )

!*****************************************************************************80
!
!! DAUB14_TRANSFORM_INVERSE inverts the DAUB14 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) Y(N), the transformed vector. 
!
!    Output, real ( kind = 8 ) X(N), the original vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 13

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     7.785205408500917D-02, &
     3.965393194819173D-01, &
     7.291320908462351D-01, &
     4.697822874051931D-01, &
    -1.439060039285649D-01, &
    -2.240361849938749D-01, &
     7.130921926683026D-02, &
     8.061260915108307D-02, &
    -3.802993693501441D-02, &
    -1.657454163066688D-02, &
     1.255099855609984D-02, &
     4.295779729213665D-04, &
    -1.801640704047490D-03, &
     3.537137999745202D-04 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  x(1:n) = y(1:n)

  m = 4
  q = ( p - 1 ) / 2

  do while ( m <= n )
  
    z(1:m) = 0.0D+00

    j = 1

    do i = - q + 1, m / 2 - q
      
      do k = 0, p - 1, 2
        i0 = i4_wrap ( i         + k / 2,     1,         m / 2 )
        i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2 + 1, m     )
        z(j)   = z(j)   + c(p-k-1) * x(i0) + c(k+1) * x(i1)
        z(j+1) = z(j+1) + c(p-k)   * x(i0) - c(k)   * x(i1)
      end do

      j = j + 2

    end do

    x(1:m) = z(1:m)

    m = m * 2

  end do

  return
end
subroutine daub16_transform ( n, x, y )

!*****************************************************************************80
!
!! DAUB16_TRANSFORM computes the DAUB16 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) X(N), the vector to be transformed. 
!
!    Output, real ( kind = 8 ) Y(N), the transformed vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 15

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     5.441584224310400D-02, &
     3.128715909142999D-01, &
     6.756307362972898D-01, &
     5.853546836542067D-01, &
    -1.582910525634930D-02, &
    -2.840155429615469D-01, &
     4.724845739132827D-04, &
     1.287474266204784D-01, &
    -1.736930100180754D-02, &
    -4.408825393079475D-02, &
     1.398102791739828D-02, &
     8.746094047405776D-03, &
    -4.870352993451574D-03, &
    -3.917403733769470D-04, &
     6.754494064505693D-04, &
    -1.174767841247695D-04 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  y(1:n) = x(1:n)

  m = n
  q = ( p - 1 ) / 2

  do while ( 4 <= m )
  
    i = 1
    z(1:m) = 0.0D+00

    do j = 1, m - 1, 2

      do k = 0, p - 1, 2
        j0 = i4_wrap ( j + k,     1, m )
        j1 = i4_wrap ( j + k + 1, 1, m )
        z(i)     = z(i)     + c(  k) * y(j0) + c(  k+1) * y(j1)
        z(i+m/2) = z(i+m/2) + c(p-k) * y(j0) - c(p-k-1) * y(j1)
      end do

      i = i + 1

    end do

    y(1:m) = z(1:m)

    m = m / 2

  end do

  return
end
subroutine daub16_transform_inverse ( n, y, x )

!*****************************************************************************80
!
!! DAUB16_TRANSFORM_INVERSE inverts the DAUB16 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) Y(N), the transformed vector. 
!
!    Output, real ( kind = 8 ) X(N), the original vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 15

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     5.441584224310400D-02, &
     3.128715909142999D-01, &
     6.756307362972898D-01, &
     5.853546836542067D-01, &
    -1.582910525634930D-02, &
    -2.840155429615469D-01, &
     4.724845739132827D-04, &
     1.287474266204784D-01, &
    -1.736930100180754D-02, &
    -4.408825393079475D-02, &
     1.398102791739828D-02, &
     8.746094047405776D-03, &
    -4.870352993451574D-03, &
    -3.917403733769470D-04, &
     6.754494064505693D-04, &
    -1.174767841247695D-04 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  x(1:n) = y(1:n)

  m = 4
  q = ( p - 1 ) / 2

  do while ( m <= n )
  
    z(1:m) = 0.0D+00

    j = 1

    do i = - q + 1, m / 2 - q
      
      do k = 0, p - 1, 2
        i0 = i4_wrap ( i         + k / 2,     1,         m / 2 )
        i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2 + 1, m     )
        z(j)   = z(j)   + c(p-k-1) * x(i0) + c(k+1) * x(i1)
        z(j+1) = z(j+1) + c(p-k)   * x(i0) - c(k)   * x(i1)
      end do

      j = j + 2

    end do

    x(1:m) = z(1:m)

    m = m * 2

  end do

  return
end
subroutine daub18_transform ( n, x, y )

!*****************************************************************************80
!
!! DAUB18_TRANSFORM computes the DAUB18 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) X(N), the vector to be transformed. 
!
!    Output, real ( kind = 8 ) Y(N), the transformed vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 17

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     3.807794736387834D-02, &
     2.438346746125903D-01, &
     6.048231236901111D-01, &
     6.572880780513005D-01, &
     1.331973858250075D-01, &
    -2.932737832791749D-01, &
    -9.684078322297646D-02, &
     1.485407493381063D-01, &
     3.072568147933337D-02, &
    -6.763282906132997D-02, &
     2.509471148314519D-04, &
     2.236166212367909D-02, &
    -4.723204757751397D-03, &
    -4.281503682463429D-03, &
     1.847646883056226D-03, &
     2.303857635231959D-04, &
    -2.519631889427101D-04, &
     3.934732031627159D-05 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  y(1:n) = x(1:n)

  m = n
  q = ( p - 1 ) / 2

  do while ( 4 <= m )
  
    i = 1
    z(1:m) = 0.0D+00

    do j = 1, m - 1, 2

      do k = 0, p - 1, 2
        j0 = i4_wrap ( j + k,     1, m )
        j1 = i4_wrap ( j + k + 1, 1, m )
        z(i)     = z(i)     + c(  k) * y(j0) + c(  k+1) * y(j1)
        z(i+m/2) = z(i+m/2) + c(p-k) * y(j0) - c(p-k-1) * y(j1)
      end do

      i = i + 1

    end do

    y(1:m) = z(1:m)

    m = m / 2

  end do

  return
end
subroutine daub18_transform_inverse ( n, y, x )

!*****************************************************************************80
!
!! DAUB18_TRANSFORM_INVERSE inverts the DAUB18 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) Y(N), the transformed vector. 
!
!    Output, real ( kind = 8 ) X(N), the original vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 17

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     3.807794736387834D-02, &
     2.438346746125903D-01, &
     6.048231236901111D-01, &
     6.572880780513005D-01, &
     1.331973858250075D-01, &
    -2.932737832791749D-01, &
    -9.684078322297646D-02, &
     1.485407493381063D-01, &
     3.072568147933337D-02, &
    -6.763282906132997D-02, &
     2.509471148314519D-04, &
     2.236166212367909D-02, &
    -4.723204757751397D-03, &
    -4.281503682463429D-03, &
     1.847646883056226D-03, &
     2.303857635231959D-04, &
    -2.519631889427101D-04, &
     3.934732031627159D-05 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  x(1:n) = y(1:n)

  m = 4
  q = ( p - 1 ) / 2

  do while ( m <= n )
  
    z(1:m) = 0.0D+00

    j = 1

    do i = - q + 1, m / 2 - q
      
      do k = 0, p - 1, 2
        i0 = i4_wrap ( i         + k / 2,     1,         m / 2 )
        i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2 + 1, m     )
        z(j)   = z(j)   + c(p-k-1) * x(i0) + c(k+1) * x(i1)
        z(j+1) = z(j+1) + c(p-k)   * x(i0) - c(k)   * x(i1)
      end do

      j = j + 2

    end do

    x(1:m) = z(1:m)

    m = m * 2

  end do

  return
end
subroutine daub20_transform ( n, x, y )

!*****************************************************************************80
!
!! DAUB20_TRANSFORM computes the DAUB20 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) X(N), the vector to be transformed. 
!
!    Output, real ( kind = 8 ) Y(N), the transformed vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 19

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     2.667005790055555D-02, &
     1.881768000776914D-01, &
     5.272011889317255D-01, &
     6.884590394536035D-01, &
     2.811723436605774D-01, &
    -2.498464243273153D-01, &
    -1.959462743773770D-01, &
     1.273693403357932D-01, &
     9.305736460357235D-02, &
    -7.139414716639708D-02, &
    -2.945753682187581D-02, &
     3.321267405934100D-02, &
     3.606553566956169D-03, &
    -1.073317548333057D-02, &
     1.395351747052901D-03, &
     1.992405295185056D-03, &
    -6.858566949597116D-04, &
    -1.164668551292854D-04, &
     9.358867032006959D-05, &
    -1.326420289452124D-05 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  y(1:n) = x(1:n)

  m = n
  q = ( p - 1 ) / 2

  do while ( 4 <= m )
  
    i = 1
    z(1:m) = 0.0D+00

    do j = 1, m - 1, 2

      do k = 0, p - 1, 2
        j0 = i4_wrap ( j + k,     1, m )
        j1 = i4_wrap ( j + k + 1, 1, m )
        z(i)     = z(i)     + c(  k) * y(j0) + c(  k+1) * y(j1)
        z(i+m/2) = z(i+m/2) + c(p-k) * y(j0) - c(p-k-1) * y(j1)
      end do

      i = i + 1

    end do

    y(1:m) = z(1:m)

    m = m / 2

  end do

  return
end
subroutine daub20_transform_inverse ( n, y, x )

!*****************************************************************************80
!
!! DAUB20_TRANSFORM_INVERSE inverts the DAUB20 transform of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!    N must be a power of 2 and at least 4.
!
!    Input, real ( kind = 8 ) Y(N), the transformed vector. 
!
!    Output, real ( kind = 8 ) X(N), the original vector.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: p = 19

  real ( kind = 8 ), dimension(0:p) :: c = (/ &
     2.667005790055555D-02, &
     1.881768000776914D-01, &
     5.272011889317255D-01, &
     6.884590394536035D-01, &
     2.811723436605774D-01, &
    -2.498464243273153D-01, &
    -1.959462743773770D-01, &
     1.273693403357932D-01, &
     9.305736460357235D-02, &
    -7.139414716639708D-02, &
    -2.945753682187581D-02, &
     3.321267405934100D-02, &
     3.606553566956169D-03, &
    -1.073317548333057D-02, &
     1.395351747052901D-03, &
     1.992405295185056D-03, &
    -6.858566949597116D-04, &
    -1.164668551292854D-04, &
     9.358867032006959D-05, &
    -1.326420289452124D-05 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  x(1:n) = y(1:n)

  m = 4
  q = ( p - 1 ) / 2

  do while ( m <= n )
  
    z(1:m) = 0.0D+00

    j = 1

    do i = - q + 1, m / 2 - q
      
      do k = 0, p - 1, 2
        i0 = i4_wrap ( i         + k / 2,     1,         m / 2 )
        i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2 + 1, m     )
        z(j)   = z(j)   + c(p-k-1) * x(i0) + c(k+1) * x(i1)
        z(j+1) = z(j+1) + c(p-k)   * x(i0) - c(k)   * x(i1)
      end do

      j = j + 2

    end do

    x(1:m) = z(1:m)

    m = m * 2

  end do

  return
end
function i4_is_power_of_2 ( n )

!*****************************************************************************80
!
!! I4_IS_POWER_OF_2 reports whether an I4 is a power of 2.
!
!  Discussion:
!
!    The powers of 2 are 1, 2, 4, 8, 16, and so on.
!
!    An I4 is an integer ( kind = 4 ) value.
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
!    Input, integer ( kind = 4 ) N, the integer to be tested.
!
!    Output, logical I4_IS_POWER_OF_2, is TRUE if N is a power of 2.
!
  implicit none

  logical              i4_is_power_of_2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_copy

  n_copy = n
  i4_is_power_of_2 = .false.

  if ( n_copy <= 0 ) then
    return
  end if

  do while ( n_copy /= 1 )

    if ( mod ( n_copy, 2 ) == 1 ) then
      return
    end if

    n_copy = n_copy / 2

  end do

  i4_is_power_of_2 = .true.

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
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
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, a value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) value
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8mat_is_identity ( n, a, error_frobenius )

!*****************************************************************************80
!
!! R8MAT_IS_IDENTITY determines if an R8MAT is the identity.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = 8 ) values.
!
!    The routine returns the Frobenius norm of A - I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) ERROR_FROBENIUS, the Frobenius norm
!    of the difference matrix A - I, which would be exactly zero
!    if A were the identity matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) error_frobenius
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) value

  error_frobenius = 0.0D+00

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        error_frobenius = error_frobenius + ( a(i,j) - 1.0D+00 )**2
      else
        error_frobenius = error_frobenius + a(i,j)**2
      end if
    end do 
  end do

  error_frobenius = sqrt ( error_frobenius )

  return
end
subroutine r8vec_conjugate ( n, c, d )

!*****************************************************************************80
!
!! R8VEC_CONJUGATE reverses a vector and negates even-indexed entries.
!
!  Discussion:
!
!    There are many times in wavelet computations when such an operation
!    is invoked.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, real ( kind = 8 ) C(N), the input vector.
!
!    Output, real ( kind = 8 ) D(N), the "conjugated" vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)

  d(1:n) = c(n:-1:1)

  d(2:2:n) = - d(2:2:n)

  return
end
subroutine r8vec_convolution ( m, x, n, y, z )

!*****************************************************************************80
!
!! R8VEC_CONVOLUTION returns the convolution of two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The I-th entry of the convolution can be formed by summing the products 
!    that lie along the I-th diagonal of the following table:
!
!    Y3 | 3   4   5   6   7
!    Y2 | 2   3   4   5   6
!    Y1 | 1   2   3   4   5
!       +------------------
!        X1  X2  X3  X4  X5
!
!    which will result in:
!
!    Z = ( X1 * Y1,
!          X1 * Y2 + X2 * Y1,
!          X1 * Y3 + X2 * Y2 + X3 * Y1,
!                    X2 * Y3 + X3 * Y2 + X4 * Y1,
!                              X3 * Y3 + X4 * Y2 + X5 * Y1,
!                                        X4 * Y3 + X5 * Y2,
!                                                  X5 * Y3 )
!            
!  Example:
!
!    Input:
!
!      X = (/ 1, 2, 3, 4 /)
!      Y = (/ -1, 5, 3 /)
!
!    Output:
!
!      Z = (/ -1, 3, 10, 17, 29, 12 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of X.
!
!    Input, real ( kind = 8 ) X(M), the first vector to be convolved.
!
!    Input, integer ( kind = 4 ) N, the dimension of Y.
!
!    Input, real ( kind = 8 ) Y(N), the second vector to be convolved.
!
!    Output, real ( kind = 8 ) Z(M+N-1), the convolution of X and Y.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(m+n-1)

  z(1:m+n-1) = 0.0D+00

  do j = 1, n
    z(j:j+m-1) = z(j:j+m-1) + x(1:m) * y(j)
  end do

  return
end
subroutine r8vec_linspace ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = 8 ) A(N), a vector of linearly spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_first
  real ( kind = 8 ) a_last
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) / 2.0D+00

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * a_first &
             + real (     i - 1, kind = 8 ) * a_last ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
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

  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

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
