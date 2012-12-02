program main

!*****************************************************************************80
!
!! MAIN is the main program for CVT_PRB.
!
!  Discussion:
!
!    CVT_PRB calls a set of problems for CVT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CVT library.'
 
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
  call test12 ( )
  call test13 ( )
  call test14 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CVT with uniform initialization and uniform sampling.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) batch
  real ( kind = 8 ) energy
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'

  batch = 1000
  init = 0
  init_string = 'uniform'
  it_max = 40
  it_fixed = 1
  sample = 0
  sample_num = 10000
  sample_string = 'uniform'
  seed = 123456789

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  Generators (rows):' )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 repeats test 1, but uses twice as many iterations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) batch
  real ( kind = 8 ) energy
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) '  Repeat test 1, but with twice the number of iterations.'

  batch = 1000
  init = 0
  init_string = 'uniform'
  it_max = 80
  it_fixed = 1
  sample = 0
  sample_num = 10000
  sample_string = 'uniform'
  seed = 123456789

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  Generators (rows):' )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 repeats test 1 but uses 100 times as many sample points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) batch
  real ( kind = 8 ) energy
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) '  Repeat test 1, but with 100 times the sample points.'

  batch = 1000
  init = 0
  init_string = 'uniform'
  it_max = 40
  it_fixed = 1
  sample = 0
  sample_num = 1000000
  sample_string = 'uniform'
  seed = 123456789

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  Generators (rows):' )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 repeats test 1 with uniform initialization and Halton sampling.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) batch
  real ( kind = 8 ) energy
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) '  Repeat test 1, but with Halton sampling.'

  batch = 1000
  init = 0
  init_string = 'uniform'
  it_max = 40
  it_fixed = 1
  sample = 1
  sample_num = 10000
  sample_string = 'halton'
  seed = 123456789

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  Generators (rows):' )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 repeats test 1 with uniform initialization and grid sampling.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) batch
  real ( kind = 8 ) energy
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) '  Repeat test 1, but with grid sampling.'

  batch = 1000
  init = 0
  init_string = 'uniform'
  it_max = 40
  it_fixed = 1
  sample = 2
  sample_num = 10000
  sample_string = 'grid'
  seed = 123456789

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  Generators (rows):' )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 repeats test 1 with uniform initialization and RANDOM_NUMBER sampling.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) batch
  real ( kind = 8 ) energy
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) '  Repeat test 1, with FORTRAN90 RANDOM_NUMBER sampling.'

  batch = 1000
  init = 0
  init_string = 'uniform'
  it_max = 40
  it_fixed = 1
  sample = -1
  sample_num = 10000
  sample_string = 'RANDOM'
  seed = 123456789

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  Generators (rows):' )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 repeats test 1 with a different seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) batch
  real ( kind = 8 ) energy
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) '  Repeat test 1, but with a different seed.'

  batch = 1000
  init = 0
  init_string = 'uniform'
  it_max = 40
  it_fixed = 1
  sample = 0
  sample_num = 10000
  sample_string = 'uniform'
  seed = 987654321

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  Generators (rows):' )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 repeats test 1 with a different batch size.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) batch
  real ( kind = 8 ) energy
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) '  Repeat test 1 with a different batch size.'

  batch = 5
  init = 0
  init_string = 'uniform'
  it_max = 40
  it_fixed = 1
  sample = 0
  sample_num = 10000
  sample_string = 'uniform'
  seed = 123456789

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  Generators (rows):' )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests CVT with IT_FIXED = IT_MAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) batch
  real ( kind = 8 ) energy
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) '  Repeat test 1 with a fixed sample set.'
  write ( *, '(a)' ) '  (IT_FIXED = IT_MAX)'

  batch = 1000
  init = 0
  init_string = 'uniform'
  it_max = 40
  it_fixed = it_max
  sample = 0
  sample_num = 10000
  sample_string = 'uniform'
  seed = 123456789

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  Generators (rows):' )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 generates 100 points in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) batch
  real ( kind = 8 ) energy
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) '  Compute 100 points in 3D.'

  batch = 1000
  init = 0
  init_string = 'uniform'
  it_max = 40
  it_fixed = 1
  sample = 0
  sample_num = 10000
  sample_string = 'uniform'
  seed = 123456789

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print_some ( dim_num, n, r, 1, 1, dim_num, 10, &
    '  First 10 generators (rows):' )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests CVT.
!
!  Discussion:
!
!    In this test, we initialize the generators to grid points; this is 
!    an unstable CVT solution.  The data would "prefer" to be in a
!    different form.  However, even if we take 2000 steps of CVT iteration,
!    the data is still only slowly progressing towards that other 
!    configuration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 16
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) batch
  real ( kind = 8 ) energy
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) ngrid
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init
  integer ( kind = 4 ) tuple(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we initialize the generators to'
  write ( *, '(a)' ) '  grid points; this is an unstable CVT solution.'

  batch = 1000
  init = 4
  init_string = 'user initialization'
  it_max = 40
  it_fixed = 1
  sample = 0
  sample_num = 1000
  sample_string = 'uniform'
  seed = 123456789

  seed_init = seed
!
!  Initialize the tuple generator.
!
  rank = -1
  ngrid = 4
  call tuple_next_fast ( ngrid, dim_num, rank, tuple )
!
!  Pick points on a grid.
!
  do rank = 0, n-1
    call tuple_next_fast ( ngrid, dim_num, rank, tuple )
    r(1:dim_num,rank+1) = real ( 2 * tuple(1:dim_num) - 1, kind = 8 ) &
                     / real ( 2 * ngrid, kind = 8 ) 
  end do

  call r8mat_transpose_print ( dim_num, n, r, '  Initial generators (rows):' )

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )

  call r8mat_transpose_print ( dim_num, n, r, '  Final generators (rows):' )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests CVT with 'RANDOM' initialization.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) batch
  real ( kind = 8 ) energy
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  The "random" initialization option calls the'
  write ( *, '(a)' ) '  system random number generator.  There is some'
  write ( *, '(a)' ) '  question about whether this works correctly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The test is as follows:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CVT call #1:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    DIM_NUM   =      2'
  write ( *, '(a)' ) '    N         =     10'
  write ( *, '(a)' ) '    INIT      =     -1'
  write ( *, '(a)' ) '    IT_MAX    =      0'
  write ( *, '(a)' ) '    SEED      = 100000'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Print output values of SEED and R #1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CVT call #2: (jump SEED)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    DIM_NUM   =      2'
  write ( *, '(a)' ) '    N         =     10'
  write ( *, '(a)' ) '    INIT      =     -1'
  write ( *, '(a)' ) '    IT_MAX    =      0'
  write ( *, '(a)' ) '    SEED      = 200000.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Print output values of SEED and R #2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CVT call #3: (restore SEED)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    DIM_NUM   =      2'
  write ( *, '(a)' ) '    N         =     10'
  write ( *, '(a)' ) '    INIT      =     -1'
  write ( *, '(a)' ) '    IT_MAX    =      0'
  write ( *, '(a)' ) '    SEED_INIT = 100000'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Print output values of SEED and R #3.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We expect that:'
  write ( *, '(a)' ) '  * the values of R #1 and R #2 differ;'
  write ( *, '(a)' ) '  AND'
  write ( *, '(a)' ) '  * the values of R #1 and R #3 agree.'
!
!  Run #1.
!
  batch = 1000
  init = -1
  init_string = 'random'
  it_max = 0
  it_fixed = 1
  sample = 0
  sample_num = 10000
  sample_string = 'uniform'
  seed = 100000

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  R #1:' )
!
!  Run #2.
!
  batch = 1000
  init = -1
  init_string = 'random'
  it_max = 0
  it_fixed = 1
  sample = 0
  sample_num = 10000
  sample_string = 'uniform'
  seed = 200000

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  R #2:' )
!
!  Run #3.
!
  batch = 1000
  init = -1
  init_string = 'random'
  it_max = 0
  it_fixed = 1
  sample = 0
  sample_num = 10000
  sample_string = 'uniform'
  seed = 100000

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  R #3:' )

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests CVT with the "user" routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) batch
  logical comment
  real ( kind = 8 ) energy
  character ( len = 80 ) :: file_out_name = 'cvt_circle.txt'
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) '  In this example, we call the "USER" routine,'
  write ( *, '(a)' ) '  which allows the user to define the geometry and'
  write ( *, '(a)' ) '  density implicitly, by returning sample points.'

  batch = 1000
  init = 3
  init_string = 'user'
  it_max = 40
  it_fixed = 1
  sample = 3
  sample_num = 10000
  sample_string = 'user'
  seed = 123456789

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_write ( file_out_name, dim_num, n, r )

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 generates a 10 point CVT on [0,1].
!
!  Discussion:
!
!    Generate 10 CVT points on the interval [0,1].
!    We expect them to be at 1/20, 3/20, 5/20, ..., 19/20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 1

  integer ( kind = 4 ) batch
  real    ( kind = 8 ) energy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) init
  character ( len = 80 ) init_string
  real    ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real    ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  Generate a CVT in the interval [0,1] using 10 points.'

  batch = 10000
  init = 0
!
!  TEMP
!
  if ( .false. ) then
    init = 4
    do i = 1, n
      r(1,i) = real ( 2 * i - 1, kind = 8 ) / real ( 2 * n, kind = 8 )
    end do
  end if

  init_string = 'uniform'
  it_max = 40
  it_fixed = 1
  sample = 0
  sample_num = 100000
  sample_string = 'uniform'
  seed = 123456789

  seed_init = seed

  call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, &
    it_fixed, seed, r, it_num, it_diff, energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)'   ) '  Number of points N =     ', n
  write ( *, '(a,i12)'   ) '  Initial SEED =           ', seed_init
  write ( *, '(a,i12)'   ) '  Current SEED =           ', seed
  write ( *, '(a)'       ) '  INIT =                   "' &
    // trim ( init_string ) // '".'
  write ( *, '(a,i12)'   ) '  Max iterations IT_MAX =  ', it_max
  write ( *, '(a,i12)'   ) '  IT_FIXED (fixed samples) ', it_fixed
  write ( *, '(a,i12)'   ) '  Iterations IT_NUM =      ', it_num
  write ( *, '(a,g14.6)' ) '  Difference IT_DIFF =     ', it_diff
  write ( *, '(a,g14.6)' ) '  CVT ENERGY =             ', energy
  write ( *, '(a)'       ) '  SAMPLE =                 "' &
    // trim ( sample_string ) // '".'
  write ( *, '(a,i12)'   ) '  Samples SAMPLE_NUM    =  ', &
    sample_num
  write ( *, '(a,i12)'   ) '  Sampling BATCH size =    ', batch
  write ( *, '(a,g14.6)' ) '  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  
  call r8mat_transpose_print ( dim_num, n, r, '  Generators (rows):' )

  return
end
