program main

!*****************************************************************************80
!
!! SIMPLEX_COORDINATES_PRB tests SIMPLEX_COORDINATES.
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
  implicit none

  integer ( kind = 4 ) n

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIMPLEX_COORDINATES_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SIMPLEX_COORDINATES library.'

  n = 3
  call test01 ( n )
  call test02 ( n )

  n = 4
  call test01 ( n )
  call test02 ( n )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIMPLEX_COORDINATES_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( n )

!*****************************************************************************80
!
!! TEST01 calls SIMPLEX_COORDINATES1.
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
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) side
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume2
  real ( kind = 8 ) x(n,n+1)
  real ( kind = 8 ) xtx(n+1,n+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Call SIMPLEX_COORDINATES1'

  call simplex_coordinates1 ( n, x )

  call r8mat_transpose_print ( n, n + 1, x, '  Simplex vertex coordinates:' )

  side = sqrt ( sum ( ( x(1:n,1) - x(1:n,2) )**2 ) )

  call simplex_volume ( n, x, volume )

  volume2 = sqrt ( real ( n + 1, kind = 8 ) ) / r8_factorial ( n ) &
    / sqrt ( 2.0D+00**n ) * side**n

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Side length =     ', side
  write ( *, '(a,g14.6)' ) '  Volume =          ', volume
  write ( *, '(a,g14.6)' ) '  Expected volume = ', volume2

  xtx = matmul ( transpose ( x ), x )

  call r8mat_transpose_print ( n + 1, n + 1, xtx, '  Dot product matrix:' )

  return
end
subroutine test02 ( n )

!*****************************************************************************80
!
!! TEST02 calls SIMPLEX_COORDINATES2.
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
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) side
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume2
  real ( kind = 8 ) x(n,n+1)
  real ( kind = 8 ) xtx(n+1,n+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Call SIMPLEX_COORDINATES2'

  call simplex_coordinates2 ( n, x )

  call r8mat_transpose_print ( n, n + 1, x, '  Simplex vertex coordinates:' )

  side = sqrt ( sum ( ( x(1:n,1) - x(1:n,2) )**2 ) )

  call simplex_volume ( n, x, volume )

  volume2 = sqrt ( real ( n + 1, kind = 8 ) ) / r8_factorial ( n ) &
    / sqrt ( 2.0D+00**n ) * side**n

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Side length =     ', side
  write ( *, '(a,g14.6)' ) '  Volume =          ', volume
  write ( *, '(a,g14.6)' ) '  Expected volume = ', volume2

  xtx = matmul ( transpose ( x ), x )

  call r8mat_transpose_print ( n + 1, n + 1, xtx, '  Dot product matrix:' )

  return
end
