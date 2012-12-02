program main

!*****************************************************************************80
!
!! MAIN is the main program for RANDOM_DATA_PRB.
!
!  Discussion:
!
!    RANDOM_DATA_PRB demonstrates the routines in RANDOM_DATA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_DATA_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the RANDOM_DATA library.'

  call test005 ( )
  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test115 ( )
  call test12 ( )
  call test125 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )

  call test20 ( )
  call test205 ( )
  call test21 ( )
  call test22 ( )
  call test23 ( )
  call test24 ( )
  call test245 ( )
  call test25 ( )
  call test26 ( )
  call test264 ( )
  call test265 ( )
  call test267 ( )
  call test27 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_DATA_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests BAD_IN_SIMPLEX01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n
  character ( len = 255 ) output_filename
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable ::  x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005:'
  write ( *, '(a)' ) '  BAD_IN_SIMPLEX01 is a "bad" sampling technique'
  write ( *, '(a)' ) '  for the unit simplex.'

  do dim_num = 2, 3

    seed = 123456789

    n = 10000

    if ( dim_num == 2 ) then
      output_filename = 'bad_in_triangle.txt'
    else if ( dim_num == 3 ) then
      output_filename = 'bad_in_tetrahedron.txt'
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
    write ( *, '(a,i8)'  ) '  Number of points N =          ', n
    write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

    allocate ( x(1:dim_num,1:n) )

    call bad_in_simplex01 ( dim_num, n, seed, x )

    call r8mat_write ( output_filename, dim_num, n, x )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Data written to "' // trim ( output_filename ) // '".'

    deallocate ( x )

  end do

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BROWNIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 100

  character ( len = 80 ) :: output_filename = 'brownian.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  BROWNIAN generates Brownian motion points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call brownian ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call scale_to_block01 ( dim_num, n, x )

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to "' // trim ( output_filename ) // '".'

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests R8_NORMAL_01
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) seed_in
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  R8_NORMAL_01 generates a single normal'
  write ( *, '(a)' ) '  pseudorandom value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Seed          Seed       R8_NORMAL_01'
  write ( *, '(a)' ) '    (Input)       (Output)'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    seed_in = seed
    x = r8_normal_01 ( seed )

    write ( *, '(2x,i12,2x,i12,2x,f12.8)' ) seed_in, seed, x

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests R8_UNIFORM_01
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) seed_in
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  R8_UNIFORM_01 generates a single uniform'
  write ( *, '(a)' ) '  pseudorandom value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Seed          Seed       R8_UNIFORM_01'
  write ( *, '(a)' ) '    (Input)       (Output)'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    seed_in = seed
    x = r8_uniform_01 ( seed )

    write ( *, '(2x,i12,2x,i12,2x,f12.8)' ) seed_in, seed, x

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests GRID_IN_CUBE01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 85

  integer ( kind = 4 ), parameter :: center = 1
  character ( len = 80 ) :: output_filename = 'grid_in_cube01.txt'
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(dim_num,n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  GRID_IN_CUBE01 generates grid points'
  write ( *, '(a)' ) '  in the unit hypercube.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i8)'  ) '  CENTER option =               ', center
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call grid_in_cube01 ( dim_num, n, center, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests HALTON_IN_CIRCLE01_ACCEPT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 400

  character ( len = 80 ) :: output_filename = 'halton_in_circle01_accept.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  HALTON_IN_CIRCLE01_ACCEPT generates'
  write ( *, '(a)' ) '  Halton points in a unit circle by acceptance/rejection.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call halton_in_circle01_accept ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests HALTON_IN_CIRCLE01_MAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 400

  character ( len = 80 ) :: output_filename = 'halton_in_circle01_map.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  HALTON_IN_CIRCLE01_MAP maps'
  write ( *, '(a)' ) '  Halton points into a unit circle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call halton_in_circle01_map ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests HALTON_IN_CUBE01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 510

  character ( len = 80 ) :: output_filename = 'halton_in_cube01.txt'
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(dim_num,n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  HALTON_IN_CUBE01 generates Halton points'
  write ( *, '(a)' ) '  in the unit hypercube.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call halton_in_cube01 ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests HAMMERSLEY_IN_CUBE01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 100

  character ( len = 80 ) :: output_filename = 'hammersley_in_cube01.txt'
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(dim_num,n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  HAMMERSLEY_IN_CUBE01 generates Hammersley points'
  write ( *, '(a)' ) '  in the unit hypercube.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call hammersley_in_cube01 ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests NORMAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 1000

  character ( len = 80 ) :: output_filename = 'normal.txt'
  integer ( kind = 4 ) info
  real ( kind = 8 ), dimension ( dim_num ) :: mu = (/ 6.0D+00, 100.0D+00 /)
  real ( kind = 8 ) r(dim_num,dim_num)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ), dimension ( dim_num, dim_num ) :: v = reshape ( (/ &
    1.0D+00, 0.3D+00, 0.3D+00, 1.0D+00 /), (/ 2, 2 /) )
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  NORMAL generates normal points'
  write ( *, '(a)' ) '    in M dimensions, using a nonzero mean, and with'
  write ( *, '(a)' ) '    user-specified variance-covariance matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call r8vec_print ( dim_num, mu, '  Mean vector MU:' )

  call r8mat_print ( dim_num, dim_num, v, '  Variance-covariance matrix V:' )
 
  r(1:dim_num,1:dim_num) = v(1:dim_num,1:dim_num)

  call dpo_fa ( dim_num, r, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04 - Fatal error!'
    write ( *, '(a)' ) '  Variance-covariance matrix factorization failed.'
    write ( *, '(a,i8)' ) '  INFO = ', info
    stop
  end if

  call r8mat_print ( dim_num, dim_num, r, '  Cholesky factor R:' )

  call normal ( dim_num, n, r, mu, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call scale_to_block01 ( dim_num, n, x )

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests NORMAL_CIRCULAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 2000

  character ( len = 80 ) :: output_filename = 'normal_circular.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  NORMAL_CIRCULAR generates points in 2D'
  write ( *, '(a)' ) '    distributed according to a circular normal.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call normal_circular ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call scale_to_block01 ( dim_num, n, x )

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests NORMAL_SIMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 1000

  character ( len = 80 ) :: output_filename = 'normal_simple.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  NORMAL_SIMPLE generates normal points'
  write ( *, '(a)' ) '    in M dimensions, using a zero mean, and with'
  write ( *, '(a)' ) '    the identity as the variance-covariance matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call normal_simple ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call scale_to_block01 ( dim_num, n, x )

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test115 ( )

!*****************************************************************************80
!
!! TEST115 tests UNIFORM_IN_ANNULUS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 400

  character ( len = 80 ) :: output_filename = 'uniform_in_annulus.txt'
  real ( kind = 8 ), dimension ( 2 ) :: pc = (/ 10.0D+00, 5.0D+00 /)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 3.0D+00
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST115'
  write ( *, '(a)' ) '  UNIFORM_IN_ANNULUS generates uniform '
  write ( *, '(a)' ) '  points in an annulus by mapping.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'     ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'     ) '  Number of points N =          ', n
  write ( *, '(a,2g14.6)' ) '  Center PC(1:2) =              ', pc(1:2)
  write ( *, '(a,g14.6)'  ) '  Inner radius is R1 =          ', r1
  write ( *, '(a,g14.6)'  ) '  Outer radius is R2 =          ', r2
  write ( *, '(a,i12)'    ) '  Initial random number SEED =  ', seed
 
  call uniform_in_annulus ( pc, r1, r2, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests UNIFORM_IN_ANNULUS_ACCEPT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 400

  character ( len = 80 ) :: output_filename = 'uniform_in_annulus_accept.txt'
  real ( kind = 8 ), dimension ( 2 ) :: pc = (/ 10.0D+00, 5.0D+00 /)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 3.0D+00
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  UNIFORM_IN_ANNULUS_ACCEPT generates uniform '
  write ( *, '(a)' ) '  points in an annulus by acceptance/rejection.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'     ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'     ) '  Number of points N =          ', n
  write ( *, '(a,2g14.6)' ) '  Center PC(1:2) =              ', pc(1:2)
  write ( *, '(a,g14.6)'  ) '  Inner radius is R1 =          ', r1
  write ( *, '(a,g14.6)'  ) '  Outer radius is R2 =          ', r2
  write ( *, '(a,i12)'    ) '  Initial random number SEED =  ', seed
 
  call uniform_in_annulus_accept ( pc, r1, r2, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test125 ( )

!*****************************************************************************80
!
!! TEST125 tests UNIFORM_IN_ANNULUS_SECTOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 400

  character ( len = 80 ) :: output_filename = 'uniform_in_annulus_sector.txt'
  real ( kind = 8 ), dimension ( 2 ) :: pc = (/ 10.0D+00, 5.0D+00 /)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 3.0D+00
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ), parameter :: theta1 = 0.0D+00
  real ( kind = 8 ), parameter :: theta2 = 1.5707964D+00
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST125'
  write ( *, '(a)' ) '  UNIFORM_IN_ANNULUS_SECTOR generates uniform '
  write ( *, '(a)' ) '  points in an annular sector by mapping.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'     ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'     ) '  Number of points N =          ', n
  write ( *, '(a,2g14.6)' ) '  Center PC(1:2) =              ', pc(1:2)
  write ( *, '(a,g14.6)'  ) '  Inner radius is R1 =          ', r1
  write ( *, '(a,g14.6)'  ) '  Outer radius is R2 =          ', r2
  write ( *, '(a,g14.6)'  ) '  THETA1 =                      ', theta1
  write ( *, '(a,g14.6)'  ) '  THETA2 =                      ', theta2
  write ( *, '(a,i12)'    ) '  Initial random number SEED =  ', seed
 
  call uniform_in_annulus_sector ( pc, r1, r2, theta1, theta2, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests UNIFORM_IN_CIRCLE01_MAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2005
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 400

  character ( len = 80 ) :: output_filename = 'uniform_in_circle01_map.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  UNIFORM_IN_CIRCLE01_MAP maps uniform '
  write ( *, '(a)' ) '  points into a unit circle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call uniform_in_circle01_map ( n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests UNIFORM_IN_CUBE01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 1000

  character ( len = 80 ) :: output_filename = 'uniform_in_cube01.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  UNIFORM_IN_CUBE01 generates uniform '
  write ( *, '(a)' ) '  points in the unit hypercube.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call uniform_in_cube01 ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests UNIFORM_IN_ELLIPSOID_MAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 1000

  real ( kind = 8 ), dimension(dim_num,dim_num) :: a = reshape ( (/ &
    3.0D+00, 1.0D+00, 1.0D+00, 2.0D+00 /), (/ 2, 2 /) )
  integer ( kind = 4 ) fail_num
  character ( len = 80 ) :: output_filename = 'uniform_in_ellipisoid_map.txt'
  integer ( kind = 4 ) j
  real ( kind = 8 ) :: r = 1.0D+00
  real ( kind = 8 ) r2
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) success_num
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  UNIFORM_IN_ELLIPSOID_MAP maps uniform '
  write ( *, '(a)' ) '  points into an ellipsoid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call uniform_in_ellipsoid_map ( dim_num, n, a, r, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'
!
!  Test the data.
!
  fail_num = 0
  success_num = 0

  do j = 1, n

    r2 = sqrt ( dot_product ( x(1:dim_num,j), &
                matmul ( a(1:dim_num,1:dim_num), x(1:dim_num,j) ) ) )

    if ( r < r2 ) then
      fail_num = fail_num + 1
    else
      success_num = success_num + 1
    end if

  end do
  
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8,a)' ) fail_num, '  points failed the ellipsoid test.'
  write ( *, '(2x,i8,a)' ) success_num, ' points satisfy the ellipsoid test.'

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests UNIFORM_IN_PARALLELOGRAM_MAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 1000

  character ( len = 80 ) :: output_filename = 'uniform_in_parallelogram_map.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ), dimension ( dim_num ) :: v1 = (/ &
    0.75D+00, 0.90D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: v2 = (/ &
    0.00D+00, 0.20D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: v3 = (/ &
    1.10D+00, 0.65D+00 /)
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  UNIFORM_IN_PARALLELOGRAM_MAP maps uniform'
  write ( *, '(a)' ) '  points into an arbitrary parallelogram.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a,2f8.2)' ) '  V1 = ', v1(1:2)
  write ( *, '(a,2f8.2)' ) '  V2 = ', v2(1:2)
  write ( *, '(a,2f8.2)' ) '  V3 = ', v3(1:2)
  write ( *, '(a,2f8.2)' ) '  V4 = ', v3(1:2)+v2(1:2)-v1(1:2)

  call uniform_in_parallelogram_map ( v1, v2, v3, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call scale_to_block01 ( dim_num, n, x )

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests UNIFORM_IN_POLYGON_MAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 1000
  integer ( kind = 4 ), parameter :: nv = 10

  character ( len = 80 ) :: output_filename = 'uniform_in_polygon_map.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ), dimension ( dim_num, nv ) :: v = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.5D+00, 0.3D+00, &
    1.0D+00, 0.0D+00, &
    0.7D+00, 0.4D+00, &
    1.0D+00, 0.6D+00, &
    0.6D+00, 0.6D+00, &
    0.5D+00, 1.0D+00, &
    0.4D+00, 0.6D+00, &
    0.0D+00, 0.6D+00, &
    0.3D+00, 0.4D+00 /), (/ dim_num, nv /) )
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  UNIFORM_IN_POLYGON_MAP maps uniform '
  write ( *, '(a)' ) '  points into a polygon.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call r8mat_print ( dim_num, nv, v, '  Polygonal vertices:' )

  call uniform_in_polygon_map ( nv, v, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( 'polygon_vertices.txt', dim_num, nv, v )

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests UNIFORM_IN_SECTOR_MAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 300

  character ( len = 80 ) :: output_filename = 'uniform_in_sector_map.txt'
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 2.0D+00
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ), parameter :: t1 = 0.78D+00
  real ( kind = 8 ), parameter :: t2 = 2.35D+00
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  UNIFORM_IN_SECTOR_MAP maps uniform '
  write ( *, '(a)' ) '  points into a circular sector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R1 = ', r1
  write ( *, '(a,g14.6)' ) '  R2 = ', r2
  write ( *, '(a,g14.6)' ) '  T1 = ', t1
  write ( *, '(a,g14.6)' ) '  T2 = ', t2
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call uniform_in_sector_map ( r1, r2, t1, t2, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests UNIFORM_IN_SIMPLEX01_MAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 10000

  character ( len = 80 ) :: output_filename = 'uniform_in_simplex01_map.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  UNIFORM_IN_SIMPLEX01_MAP maps uniform '
  write ( *, '(a)' ) '  points into the unit simplex.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call uniform_in_simplex01_map ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to "' // trim ( output_filename ) // '".'

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests UNIFORM_IN_SPHERE01_MAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 1000

  character ( len = 80 ) :: output_filename = 'uniform_in_sphere01_map.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  UNIFORM_IN_SPHERE01_MAP maps uniform '
  write ( *, '(a)' ) '  points into the unit sphere.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call uniform_in_sphere01_map ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test205 ( )

!*****************************************************************************80
!
!! TEST205 tests UNIFORM_IN_TETRAHEDRON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 1000

  character ( len = 80 ) :: output_filename = 'uniform_in_tetrahedron.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ), dimension ( 3,4 ) :: v = reshape ( (/ &
    1.0D+00,  2.0D+00,  3.0D+00, &
    4.0D+00,  1.0D+00,  2.0D+00, &
    2.0D+00,  4.0D+00,  4.0D+00, &
    3.0D+00,  2.0D+00,  5.0D+00 /), (/ 3, 4 /) )
  real ( kind = 8 ) x(3,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST205'
  write ( *, '(a)' ) '  UNIFORM_IN_TETRAHEDRON returns uniform '
  write ( *, '(a)' ) '  points from a tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', 3
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call r8mat_print ( 3, 4, v, '  Tetrahedron vertices:' )

  call uniform_in_tetrahedron ( v, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, 3, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests UNIFORM_IN_TRIANGLE_MAP1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 1000

  character ( len = 80 ) :: output_filename = 'uniform_in_triangle_map1.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ), dimension ( dim_num ) :: v1 = (/ &
    0.75D+00, 0.90D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: v2 = (/ &
    0.00D+00, 0.20D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: v3 = (/ &
    0.95D+00, 0.65D+00 /)
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  UNIFORM_IN_TRIANGLE_MAP1 maps uniform '
  write ( *, '(a)' ) '  points into a triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a,2f8.2)' ) '  V1 = ', v1(1:2)
  write ( *, '(a,2f8.2)' ) '  V2 = ', v2(1:2)
  write ( *, '(a,2f8.2)' ) '  V3 = ', v3(1:2)

  call uniform_in_triangle_map1 ( v1, v2, v3, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests UNIFORM_IN_TRIANGLE_MAP2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 1000

  character ( len = 80 ) :: output_filename = 'uniform_in_triangle_map2.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ), dimension ( dim_num ) :: v1 = (/ &
    0.75D+00, 0.90D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: v2 = (/ &
    0.00D+00, 0.20D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: v3 = (/ &
    0.95D+00, 0.65D+00 /)
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  UNIFORM_IN_TRIANGLE_MAP maps uniform '
  write ( *, '(a)' ) '  points into a triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a,2f8.2)' ) '  V1 = ', v1(1:2)
  write ( *, '(a,2f8.2)' ) '  V2 = ', v2(1:2)
  write ( *, '(a,2f8.2)' ) '  V3 = ', v3(1:2)

  call uniform_in_triangle_map2 ( v1, v2, v3, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests UNIFORM_IN_TRIANGLE01_MAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 2000

  character ( len = 80 ) :: output_filename = 'uniform_in_triangle01_map.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  UNIFORM_IN_TRIANGLE01_MAP maps uniform '
  write ( *, '(a)' ) '  points into the unit triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call uniform_in_triangle01_map ( n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests UNIFORM_ON_ELLIPSOID_MAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 200

  character ( len = 80 ) :: output_filename = 'uniform_on_ellipsoid_map.txt'
  real ( kind = 8 ) a(dim_num,dim_num)
  real ( kind = 8 ) r
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  UNIFORM_ON_ELLIPSOID_MAP maps uniform '
  write ( *, '(a)' ) '  points onto an ellipsoid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  a(1:dim_num,1:dim_num) = reshape ( (/ &
    3.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00 /), (/ 2, 2 /) )
  r = 1.0D+00

  call uniform_on_ellipsoid_map ( dim_num, n, a, r, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test245 ( )

!*****************************************************************************80
!
!! TEST245 tests UNIFORM_ON_HEMISPHERE01_PHONG.
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
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: n = 50

  character ( len = 80 ) :: output_filename = 'uniform_on_hemisphere01_phong.txt'
  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST245'
  write ( *, '(a)' ) '  UNIFORM_ON_HEMISPHERE01_PHONG maps uniform '
  write ( *, '(a)' ) '  points onto the unit hemisphere with Phong density.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i8)'  ) '  Phong exponent M =            ', m
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call uniform_on_hemisphere01_phong ( n, m, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests UNIFORM_ON_SIMPLEX01_MAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 50

  character ( len = 80 ) :: output_filename = 'uniform_on_simplex01_map.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  UNIFORM_ON_SIMPLEX01_MAP maps uniform '
  write ( *, '(a)' ) '  points onto the unit simplex.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call uniform_on_simplex01_map ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26 tests UNIFORM_ON_SPHERE01_MAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 50

  character ( len = 80 ) :: output_filename = 'uniform_on_sphere01_map.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  UNIFORM_ON_SPHERE01_MAP maps uniform '
  write ( *, '(a)' ) '  points onto the unit sphere, in any dimension.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call uniform_on_sphere01_map ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test264 ( )

!*****************************************************************************80
!
!! TEST264 tests UNIFORM_ON_SPHERE01_PATCH_TP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5000

  character ( len = 80 ) :: output_filename = 'uniform_on_sphere01_patch_tp.txt'
  real ( kind = 8 ) phi1
  real ( kind = 8 ) phi2
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) tp(2,n)

  phi1 = 0.0D+00 * ( pi / 180.0D+00 )
  phi2 = 180.0D+00 * ( pi / 180.0D+00 )
  theta1 =  0.0D+00 * ( pi / 360.0D+00 )
  theta2 = 30.0D+00 * ( pi / 360.0D+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST264'
  write ( *, '(a)' ) '  UNIFORM_ON_SPHERE01_PATCH_TP maps uniform '
  write ( *, '(a)' ) '  points onto a TP (THETA,PHI) patch of the unit sphere.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Spatial dimension DIM_NUM =   ', 3
  write ( *, '(a,i8)'    ) '  Data dimension =              ', 2
  write ( *, '(a,i8)'    ) '  Number of points N =          ', n
  write ( *, '(a,g14.6)' ) '  Latitudinal angle PHI1 =      ', phi1
  write ( *, '(a,g14.6)' ) '  Latitudinal angle PHI2 =      ', phi2
  write ( *, '(a,g14.6)' ) '  Longitudinal angle THETA1 =   ', theta1
  write ( *, '(a,g14.6)' ) '  Longitudinal angle THETA2 =   ', theta2
  write ( *, '(a,i12)'   ) '  Initial random number SEED =  ', seed

  call uniform_on_sphere01_patch_tp ( n, phi1, phi2, theta1, theta2, seed, tp )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, 2, n, tp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test265 ( )

!*****************************************************************************80
!
!! TEST265 tests UNIFORM_ON_SPHERE01_PATCH_XYZ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: n = 50

  character ( len = 80 ) :: output_filename = 'uniform_on_sphere01_patch_xyz.txt'
  real ( kind = 8 ) phi1
  real ( kind = 8 ) phi2
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) x(dim_num,n)

  phi1 = 75.0D+00 * ( pi / 180.0D+00 )
  phi2 = 90.0D+00 * ( pi / 180.0D+00 )
  theta1 =  0.0D+00 * ( pi / 360.0D+00 )
  theta2 = 30.0D+00 * ( pi / 360.0D+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST265'
  write ( *, '(a)' ) '  UNIFORM_ON_SPHERE01_PATCH_XYZ maps uniform '
  write ( *, '(a)' ) '  points onto an XYZ patch of the unit sphere.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'    ) '  Number of points N =          ', n
  write ( *, '(a,g14.6)' ) '  Latitudinal angle PHI1 =      ', phi1
  write ( *, '(a,g14.6)' ) '  Latitudinal angle PHI2 =      ', phi2
  write ( *, '(a,g14.6)' ) '  Longitudinal angle THETA1 =   ', theta1
  write ( *, '(a,g14.6)' ) '  Longitudinal angle THETA2 =   ', theta2
  write ( *, '(a,i12)'   ) '  Initial random number SEED =  ', seed

  call uniform_on_sphere01_patch_xyz ( n, phi1, phi2, theta1, theta2, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test267 ( )

!*****************************************************************************80
!
!! TEST267 tests UNIFORM_ON_SPHERE01_TRIANGLE_XYZ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 500

  character ( len = 80 ) :: output_filename = 'uniform_on_sphere01_triangle_xyz.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)
  real ( kind = 8 ) x(3,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST267'
  write ( *, '(a)' ) '  UNIFORM_ON_SPHERE01_TRIANGLE_XYZ maps uniform '
  write ( *, '(a)' ) '  points onto a spherical triangle using XYZ coordinates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Spatial dimension DIM_NUM =   ', 3
  write ( *, '(a,i8)'    ) '  Number of points N =          ', n
  write ( *, '(a,i12)'   ) '  Initial random number SEED =  ', seed

  if ( .true. ) then

    call uniform_on_sphere01_map ( 3, 1, seed, v1 )
    call uniform_on_sphere01_map ( 3, 1, seed, v2 )
    call uniform_on_sphere01_map ( 3, 1, seed, v3 )

  else

    v1 = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)
    v2 = (/ 0.0D+00, 1.0D+00, 0.0D+00 /)
    v3 = (/ 0.0D+00, 0.0D+00, 1.0D+00 /)

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertices of spherical triangle:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  V1:', v1(1:3)
  write ( *, '(a,3g14.6)' ) '  V2:', v2(1:3)
  write ( *, '(a,3g14.6)' ) '  V3:', v3(1:3)

  call uniform_on_sphere01_triangle_xyz ( n, v1, v2, v3, seed, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call r8mat_write ( output_filename, 3, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests UNIFORM_WALK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 400

  character ( len = 80 ) :: output_filename = 'uniform_walk.txt'
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27:'
  write ( *, '(a)' ) '  UNIFORM_WALK generates points on a uniform random walk'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension DIM_NUM =   ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points N =          ', n
  write ( *, '(a,i12)' ) '  Initial random number SEED =  ', seed

  call uniform_walk ( dim_num, n, seed, x )

  write ( *, '(a,i12)' ) '  Final random number SEED =    ', seed

  call scale_to_block01 ( dim_num, n, x )

  call r8mat_write ( output_filename, dim_num, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( output_filename ) &
    // '".'

  return
end
