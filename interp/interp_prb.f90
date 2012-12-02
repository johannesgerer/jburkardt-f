program main

!*****************************************************************************80
!
!! MAIN is the main program for INTERP_PRB.
!
!  Discussion:
!
!    INTERP_PRB calls the INTERP tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) data_num

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTERP_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the INTERP library.'

  call test01 ( )

  call test02 ( )

  data_num = 6
  call test03 ( data_num )

  data_num = 11
  call test03 ( data_num )

  data_num = 6
  call test04 ( data_num )

  data_num = 11
  call test04 ( data_num )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTERP_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests INTERP_NEAREST on 1-dimensional data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: data_num = 11
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ) after
  integer ( kind = 4 ) before
  integer ( kind = 4 ) fat
  integer ( kind = 4 ) i
  integer ( kind = 4 ) interp
  integer ( kind = 4 ) interp_num
  real ( kind = 8 ) p
  real ( kind = 8 ) p_data(m,data_num)
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p_interp
  real ( kind = 8 ), allocatable, dimension ( : ) :: p_value
  real ( kind = 8 ) t
  real ( kind = 8 ) t_data(data_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: t_interp
  real ( kind = 8 ) t_max
  real ( kind = 8 ) t_min

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  INTERP_NEAREST evaluates a nearest-neighbor interpolant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, the function we are interpolating is'
  write ( *, '(a)' ) '  Runge''s function, with Chebyshev knots.'

  t_min = -1.0D+00
  t_max = +1.0D+00

  call cc_abscissas_ab ( t_min, t_max, data_num, t_data )

  call f_runge ( m, data_num, t_data, p_data )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension =     ', m
  write ( *, '(a,i8)' ) '  Number of data values = ', data_num
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T_data        P_data'
  write ( *, '(a)'    ) ' '
  do i = 1, data_num
    write ( *, '(2x,2g14.6)' ) t_data(i), p_data(1,i)
  end do
!
!  Our interpolation values will include the original T values, plus
!  3 new values in between each pair of original values.
!
  before = 4
  fat = 3
  after = 2

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after

  allocate ( t_interp(interp_num) )
  allocate ( p_interp(m,interp_num) )

  call r8vec_expand_linear2 ( data_num, t_data, before, fat, after, t_interp )

  call interp_nearest ( m, data_num, t_data, p_data, interp_num, &
    t_interp, p_interp )

  allocate ( p_value(1:interp_num) )

  call f_runge ( m, interp_num, t_interp, p_value )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  Interpolation:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '    T_interp      P_interp        P_exact        Error'
  write ( *, '(a)'    ) ' '

  do interp = 1, interp_num

    write ( *, '(2x,f10.4,2x,g14.6,2x,g14.6,2x,g10.2)' ) &
      t_interp(interp), p_interp(1,interp), p_value(interp), &
      p_interp(1,interp) - p_value(interp)

  end do

  deallocate ( p_interp )
  deallocate ( p_value )
  deallocate ( t_interp )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests INTERP_LINEAR on 1-dimensional data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: data_num = 11
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ) after
  integer ( kind = 4 ) before
  integer ( kind = 4 ) fat
  integer ( kind = 4 ) i
  integer ( kind = 4 ) interp
  integer ( kind = 4 ) interp_num
  real ( kind = 8 ) p
  real ( kind = 8 ) p_data(m,data_num)
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p_interp
  real ( kind = 8 ), allocatable, dimension ( : ) :: p_value
  real ( kind = 8 ) t
  real ( kind = 8 ) t_data(data_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: t_interp
  real ( kind = 8 ) t_max
  real ( kind = 8 ) t_min

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  INTERP_LINEAR evaluates a piecewise linear spline.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, the function we are interpolating is'
  write ( *, '(a)' ) '  Runge''s function, with evenly spaced knots.'

  t_min = -1.0D+00
  t_max = +1.0D+00

  call ncc_abscissas_ab ( t_min, t_max, data_num, t_data )

  call f_runge ( m, data_num, t_data, p_data )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension =     ', m
  write ( *, '(a,i8)' ) '  Number of data values = ', data_num
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T_data        P_data'
  write ( *, '(a)'    ) ' '
  do i = 1, data_num
    write ( *, '(2x,2g14.6)' ) t_data(i), p_data(1,i)
  end do
!
!  Our interpolation values will include the original T values, plus
!  3 new values in between each pair of original values.
!
  before = 4
  fat = 3
  after = 2

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after

  allocate ( t_interp(interp_num) )
  allocate ( p_interp(m,interp_num) )

  call r8vec_expand_linear2 ( data_num, t_data, before, fat, after, t_interp )

  call interp_linear ( m, data_num, t_data, p_data, interp_num, &
    t_interp, p_interp )

  allocate ( p_value(1:interp_num) )

  call f_runge ( m, interp_num, t_interp, p_value )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  Interpolation:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '    T_interp      P_interp        P_exact        Error'
  write ( *, '(a)'    ) ' '

  do interp = 1, interp_num

    write ( *, '(2x,f10.4,2x,g14.6,2x,g14.6,2x,g10.2)' ) &
      t_interp(interp), p_interp(1,interp), p_value(interp), &
      p_interp(1,interp) - p_value(interp)

  end do

  deallocate ( p_interp )
  deallocate ( p_value )
  deallocate ( t_interp )

  return
end
subroutine test03 ( data_num )

!*****************************************************************************80
!
!! TEST03 tests INTERP_LAGRANGE on 1-dimensional data, equally spaced data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data values.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ) after
  integer ( kind = 4 ) before
  integer ( kind = 4 ) fat
  integer ( kind = 4 ) i
  integer ( kind = 4 ) interp
  integer ( kind = 4 ) interp_num
  real ( kind = 8 ) p
  real ( kind = 8 ) p_data(m,data_num)
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p_interp
  real ( kind = 8 ), allocatable, dimension ( : ) :: p_value
  real ( kind = 8 ) t
  real ( kind = 8 ) t_data(data_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: t_interp
  real ( kind = 8 ) t_max
  real ( kind = 8 ) t_min

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  INTERP_LAGRANGE evaluates a polynomial interpolant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, the function we are interpolating is'
  write ( *, '(a)' ) '  Runge''s function, with evenly spaced knots.'

  t_min = -1.0D+00
  t_max = +1.0D+00

  call ncc_abscissas_ab ( t_min, t_max, data_num, t_data )

  call f_runge ( m, data_num, t_data, p_data )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension =     ', m
  write ( *, '(a,i8)' ) '  Number of data values = ', data_num
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T_data        P_data'
  write ( *, '(a)'    ) ' '
  do i = 1, data_num
    write ( *, '(2x,2g14.6)' ) t_data(i), p_data(1,i)
  end do
!
!  Our interpolation values will include the original T values, plus
!  3 new values in between each pair of original values.
!
  before = 4
  fat = 3
  after = 2

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after

  allocate ( t_interp(interp_num) )
  allocate ( p_interp(m,interp_num) )

  call r8vec_expand_linear2 ( data_num, t_data, before, fat, after, t_interp )

  call interp_lagrange ( m, data_num, t_data, p_data, interp_num, &
    t_interp, p_interp )

  allocate ( p_value(1:interp_num) )

  call f_runge ( m, interp_num, t_interp, p_value )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  Interpolation:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '    T_interp      P_interp        P_exact        Error'
  write ( *, '(a)'    ) ' '

  do interp = 1, interp_num

    write ( *, '(2x,f10.4,2x,g14.6,2x,g14.6,2x,g10.2)' ) &
      t_interp(interp), p_interp(1,interp), p_value(interp), &
      p_interp(1,interp) - p_value(interp)

  end do

  deallocate ( p_interp )
  deallocate ( p_value )
  deallocate ( t_interp )

  return
end
subroutine test04 ( data_num )

!*****************************************************************************80
!
!! TEST04 tests INTERP_LAGRANGE on 1-dimensional data, Clenshaw-Curtis data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data values.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ) after
  integer ( kind = 4 ) before
  integer ( kind = 4 ) fat
  integer ( kind = 4 ) i
  integer ( kind = 4 ) interp
  integer ( kind = 4 ) interp_num
  real ( kind = 8 ) p
  real ( kind = 8 ) p_data(m,data_num)
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p_interp
  real ( kind = 8 ), allocatable, dimension ( : ) :: p_value
  real ( kind = 8 ) t
  real ( kind = 8 ) t_data(data_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: t_interp
  real ( kind = 8 ) t_max
  real ( kind = 8 ) t_min

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  INTERP_LAGRANGE evaluates a polynomial interpolant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, the function we are interpolating is'
  write ( *, '(a)' ) '  Runge''s function, with Clenshaw Curtis knots.'

  t_min = -1.0D+00
  t_max = +1.0D+00

  call cc_abscissas_ab ( t_min, t_max, data_num, t_data )

  call f_runge ( m, data_num, t_data, p_data )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension =     ', m
  write ( *, '(a,i8)' ) '  Number of data values = ', data_num
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T_data        P_data'
  write ( *, '(a)'    ) ' '
  do i = 1, data_num
    write ( *, '(2x,2g14.6)' ) t_data(i), p_data(1,i)
  end do
!
!  Our interpolation values will include the original T values, plus
!  3 new values in between each pair of original values.
!
  before = 4
  fat = 3
  after = 2

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after

  allocate ( t_interp(interp_num) )
  allocate ( p_interp(m,interp_num) )

  call r8vec_expand_linear2 ( data_num, t_data, before, fat, after, t_interp )

  call interp_lagrange ( m, data_num, t_data, p_data, interp_num, &
    t_interp, p_interp )

  allocate ( p_value(1:interp_num) )

  call f_runge ( m, interp_num, t_interp, p_value )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  Interpolation:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '    T_interp      P_interp        P_exact        Error'
  write ( *, '(a)'    ) ' '

  do interp = 1, interp_num

    write ( *, '(2x,f10.4,2x,g14.6,2x,g14.6,2x,g10.2)' ) &
      t_interp(interp), p_interp(1,interp), p_value(interp), &
      p_interp(1,interp) - p_value(interp)

  end do

  deallocate ( p_interp )
  deallocate ( p_value )
  deallocate ( t_interp )

  return
end
subroutine f_runge ( m, n, x, f )

!*****************************************************************************80
!
!! F_RUNGE evaluates the Runge function.
!
!  Discussion:
!
!    Interpolation of the Runge function at evenly spaced points in [-1,1]
!    is a common test.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(m,n)

  f(1:n) = 1.0D+00 / ( 1.0D+00 + 25.0D+00 * sum ( x(1:m,1:n)**2, 1 ) )

  return
end
