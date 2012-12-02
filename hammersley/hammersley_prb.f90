program main

!*****************************************************************************80
!
!! MAIN is the main program for HAMMERSLEY_PRB.
!
!  Discussion:
!
!    HAMMERSLEY_PRB calls a set of problems for HAMMERSLEY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HAMMERSLEY_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the HAMMERSLEY library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HAMMERSLEY_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests I4_TO_HAMMERSLEY_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 4

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  integer ( kind = 4 ), parameter :: nmax = 1000
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of '
  write ( *, '(a)' ) '  a Hammersley sequence on a single call.'
  write ( *, '(a)' ) '  All arguments are specified explicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we compute the first 10 elements'
  write ( *, '(a)' ) '  of a "classical" Hammersley sequence, and then'
  write ( *, '(a)' ) '  the "last" 10 elements.'

  step = 1
  seed(1:dim_num) = 0
  leap(1:dim_num) = 1
  base(1) = -nmax
  do i = 2, dim_num
    base(i) = prime ( i - 1 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Hammersley'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,2x,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We can jump ahead in the sequence by changing STEP:'

  step = nmax - n + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  STEP = ', step

  call i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Hammersley'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,2x,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests I4_TO_HAMMERSLEY_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ), parameter :: dim_num = 4

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of '
  write ( *, '(a)' ) '  a Hammersley sequence on a single call.'
  write ( *, '(a)' ) '  All arguments are specified explicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We are free to choose the values of BASE.'
  write ( *, '(a)' ) '  Any negative value indicates a sequence of'
  write ( *, '(a)' ) '  J/(-BASE) in that coordinate.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, that is the only kind of base we use.'
  step = 0
  seed(1:dim_num) = 0
  leap(1:dim_num) = 1
  do i = 1, dim_num
    base(i) = -(10**i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Hammersley'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,2x,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests I4_TO_HAMMERSLEY_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ), parameter :: dim_num = 4

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of '
  write ( *, '(a)' ) '  a Hammersley sequence on a single call.'
  write ( *, '(a)' ) '  All arguments are specified explicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The SEED vector allows us to define the zeroth'
  write ( *, '(a)' ) '  element of the coordinate subsequence.'
  write ( *, '(a)' ) '  That is, if we ask for the STEP=0 entry of the'
  write ( *, '(a)' ) '  subsequence, we will get the SEED(I)th entry'
  write ( *, '(a)' ) '  of the full sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we use a fixed base for simplicity.'

  step = 0
  do i = 1, dim_num
    seed(i) = 10 * i
  end do
  leap(1:dim_num) = 1
  base(1:dim_num) = -100

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Hammersley'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,2x,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests I4_TO_HAMMERSLEY_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ), parameter :: dim_num = 4

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of '
  write ( *, '(a)' ) '  a Hammersley sequence on a single call.'
  write ( *, '(a)' ) '  All arguments are specified explicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The LEAP vector allows us to define the distance'
  write ( *, '(a)' ) '  (in the original sequence) between successive'
  write ( *, '(a)' ) '  subsequence elements.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A LEAP of 1 means that, once we start sampling'
  write ( *, '(a)' ) '  the sequence, we are taking every element.'
  write ( *, '(a)' ) '  A LEAP of 2 takes every other sequence element,'
  write ( *, '(a)' ) '  and so on.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we use a fixed base for simplicity.'

  step = 0
  seed(1:dim_num) = 0
  do i = 1, dim_num
    leap(i) = 2**(i-1)
  end do
  base(1:dim_num) = -100

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Hammersley'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,2x,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests I4_TO_HAMMERSLEY_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: dim_num = 4

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) leap(dim_num)
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of '
  write ( *, '(a)' ) '  a Hammersley sequence on a single call.'
  write ( *, '(a)' ) '  All arguments are specified explicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Any entry of the Hammersley sequence can be computed'
  write ( *, '(a)' ) '  immediately, without having to compute the previous'
  write ( *, '(a)' ) '  entries.  This is also true of the entries of the'
  write ( *, '(a)' ) '  leaped Hammersley subsequences we generate.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The value of a component of the Hammersley sequence'
  write ( *, '(a)' ) '  is computed directly from its index.  But there'
  write ( *, '(a)' ) '  should not be much difficulty handling indices'
  write ( *, '(a)' ) '  that go as high as a million or a billion.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we look at high index entries,'
  write ( *, '(a)' ) '  attained by large values of STEP, or SEED or LEAP.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we use the default bases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BIG VALUES OF STEP:'
  write ( *, '(a)' ) ' '

  do k = 1, 4

    step = 100**k
    do i = 1, dim_num
      seed(i) = 0
    end do
    leap(1:dim_num) = 1
    base(1) = -(step+n-1)
    do i = 2, dim_num
      base(i) = prime ( i - 1 )
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
    write ( *, '(a,i12)' ) '  N =    ', n
    write ( *, '(a,i12)' ) '  STEP = ', step
    call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
    call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
    call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

    call i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '          STEP   Hammersley'
    write ( *, '(a)' ) ' '
    do j = 1, n
      write ( *, '(2x,i12,2x,4f14.9)' ) step+j-1, r(1:dim_num,j)
    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BIG VALUES OF SEED:'
  write ( *, '(a)' ) ' '

  step = 0
  do i = 1, dim_num
    seed(i) = 100**i
  end do
  leap(1:dim_num) = 1
  base(1) = -(100+n-1)
  do i = 2, dim_num
    base(i) = prime ( i - 1 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Hammersley'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,2x,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BIG VALUES OF LEAP:'
  write ( *, '(a)' ) ' '

  step = 0
  seed(1:dim_num) = 0
  do i = 1, dim_num
    leap(i) = 100**i
  end do
  base(1) = -( n * 100 )
  do i = 2, dim_num
    base(i) = prime ( i - 1 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Hammersley'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,2x,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests HAMMERSLEY_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 4

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  integer ( kind = 4 ), parameter :: nmax = 1000
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  HAMMERSLEY_SEQUENCE computes N elements of '
  write ( *, '(a)' ) '  a Hammersley sequence on a single call.'
  write ( *, '(a)' ) '  All arguments are specified externally, by calling'
  write ( *, '(a)' ) '  various setup routines.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we compute the first 10 elements'
  write ( *, '(a)' ) '  of a "classical" Hammersley sequence, and then'
  write ( *, '(a)' ) '  the "last" 10 elements.'

  step = 1
  seed(1:dim_num) = 0
  leap(1:dim_num) = 1
  base(1) = -nmax
  do i = 2, dim_num
    base(i) = prime ( i - 1 )
  end do

  call hammersley_dim_num_set ( dim_num )
  call hammersley_step_set ( step )
  call hammersley_seed_set ( seed )
  call hammersley_leap_set ( leap )
  call hammersley_base_set ( base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call hammersley_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Hammersley'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,2x,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We can jump ahead in the sequence by changing STEP:'

  step = nmax - n + 1
  call hammersley_step_set ( step )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  STEP = ', step

  call hammersley_sequence ( dim_num, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Hammersley'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,2x,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests HALHAM_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: dim_num = 4

  integer ( kind = 4 ) base(dim_num)
  character ( len = 80 ) :: file_name = 'hammersley_04_00010.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  HALHAM_WRITE writes a Halton or Hammersley '
  write ( *, '(a)' ) '  dataset to a file'

  step = 0
  seed(1) = 1
  seed(2:dim_num) = 0
  leap(1:dim_num) = 1
  base(1) = -n
  do i = 2, dim_num
    base(i) = prime ( i - 1 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Hammersley'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,2x,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  call halham_write ( dim_num, n, step, seed, leap, base, r, file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data was written to "' // trim ( file_name ) // '".'

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests I4_TO_HAMMERSLEY_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20
  integer ( kind = 4 ), parameter :: dim_num = 4

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  integer ( kind = 4 ), parameter :: nmax = 10
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of '
  write ( *, '(a)' ) '  a Hammersley sequence on a single call.'
  write ( *, '(a)' ) '  All arguments are specified explicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we demonstrate that any coordinate of'
  write ( *, '(a)' ) '  the generalized Hammersley sequence that is generated'
  write ( *, '(a)' ) '  as a fractional sequence J/|BASE(I)| will'
  write ( *, '(a)' ) '  "wrap around".'

  step = 1
  seed(1:dim_num) = 0
  leap(1:dim_num) = 1
  base(1) = -nmax
  do i = 2, dim_num
    base(i) = prime ( i - 1 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
  write ( *, '(a,i12)' ) '  N =    ', n
  write ( *, '(a,i12)' ) '  STEP = ', step
  call i4vec_transpose_print ( dim_num, seed, '  SEED = ' )
  call i4vec_transpose_print ( dim_num, leap, '  LEAP = ' )
  call i4vec_transpose_print ( dim_num, base, '  BASE = ' )

  call i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STEP   Hammersley'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i6,2x,4g14.6)' ) step+j-1, r(1:dim_num,j)
  end do

  return
end
