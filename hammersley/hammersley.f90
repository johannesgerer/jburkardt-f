function arc_cosine ( c )

!*****************************************************************************80
!
!! ARC_COSINE computes the arc cosine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ACOS routine with an input argument that is
!    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
!    This routine truncates arguments outside the range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, the argument.
!
!    Output, real ( kind = 8 ) ARC_COSINE, an angle whose cosine is C.
!
  implicit none

  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) c
  real ( kind = 8 ) c2

  c2 = c
  c2 = max ( c2, real ( -1.0, kind = 8 ) )
  c2 = min ( c2, real ( +1.0, kind = 8 ) )

  arc_cosine = acos ( c2 )

  return
end
function atan4 ( y, x )

!*****************************************************************************80
!
!! ATAN4 computes the inverse tangent of the ratio Y / X.
!
!  Discussion:
!
!    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
!    the built in functions ATAN and ATAN2 already do.
!
!    However:
!
!    * ATAN4 always returns a positive angle, between 0 and 2 PI,
!      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
!      and [-PI,+PI] respectively;
!
!    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
!     function by contrast always returns an angle in the first or fourth
!     quadrants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Y, X, two quantities which represent the
!    tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
!
!    Output, real ( kind = 8 ) ATAN4, an angle between 0 and 2 * PI,
!    whose tangent is (Y/X), and which lies in the appropriate quadrant
!    so that the signs of its cosine and sine match those of X and Y.
!
  implicit none

  real ( kind = 8 ) abs_x
  real ( kind = 8 ) abs_y
  real ( kind = 8 ) atan4
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta_0
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Special cases:
!
  if ( x == 0.0 ) then

    if ( 0.0 < y ) then
      theta = pi / real ( 2.0, kind = 8 )
    else if ( y < 0.0 ) then
      theta = real ( 3.0, kind = 8 ) * pi / real ( 2.0, kind = 8 )
    else if ( y == 0.0 ) then
      theta = 0.0
    end if

  else if ( y == 0.0 ) then

    if ( 0.0D+00 < x ) then
      theta = 0.0D+00
    else if ( x < 0.0D+00 ) then
      theta = pi
    end if
!
!  We assume that ATAN2 is correct when both arguments are positive.
!
  else

    abs_y = abs ( y )
    abs_x = abs ( x )

    theta_0 = atan2 ( abs_y, abs_x )

    if ( 0.0D+00 < x .and. 0.0D+00 < y ) then
      theta = theta_0
    else if ( x < 0.0D+00 .and. 0.0D+00 < y ) then
      theta = pi - theta_0
    else if ( x < 0.0D+00 .and. y < 0.0D+00 ) then
      theta = pi + theta_0
    else if ( 0.0D+00 < x .and. y < 0.0D+00 ) then
      theta = 2.0D+00 * pi - theta_0
    end if

  end if

  atan4 = theta

  return
end
subroutine get_seed ( seed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
!
  implicit none

  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 ) today
  integer ( kind = 4 ) values(8)
  character ( len = 5 ) zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) / real (  11.0D+00, kind = 8 )
  temp = temp + real ( values(3) - 1, kind = 8 ) / real (  30.0D+00, kind = 8 )
  temp = temp + real ( values(5),     kind = 8 ) / real (  23.0D+00, kind = 8 )
  temp = temp + real ( values(6),     kind = 8 ) / real (  59.0D+00, kind = 8 )
  temp = temp + real ( values(7),     kind = 8 ) / real (  59.0D+00, kind = 8 )
  temp = temp + real ( values(8),     kind = 8 ) / real ( 999.0D+00, kind = 8 )
  temp = temp / real ( 6.0D+00, kind = 8 )

  if ( temp <= 0.0D+00 ) then
    temp = real ( 1.0D+00, kind = 8 ) / real ( 3.0D+00, kind = 8 )
  else if ( 1.0 <= temp ) then
    temp = real ( 2.0D+00, kind = 8 ) / real ( 3.0D+00, kind = 8 )
  end if

  seed = int ( real ( huge ( 1 ), kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == huge ( 1 ) ) then
    seed = seed - 1
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
function halham_leap_check ( dim_num, leap )

!*****************************************************************************80
!
!! HALHAM_LEAP_CHECK checks LEAP for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the leap vector.
!
!    Output, logical, HALHAM_LEAP_CHECK, is true if LEAP is legal.
!
  implicit none

  integer ( kind = 4 ) dim_num

  logical halham_leap_check
  integer ( kind = 4 ) leap(dim_num)

  if ( any ( leap(1:dim_num) < 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_LEAP_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some entry of LEAP < 1!'
    write ( *, '(a)' ) ' '
    call i4vec_transpose_print ( dim_num, leap, 'LEAP:  ' )
    halham_leap_check = .false.
  else
    halham_leap_check = .true.
  end if

  return
end
function halham_n_check ( n )

!*****************************************************************************80
!
!! HALHAM_N_CHECK checks N for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, logical HALHAM_N_CHECK, is true if N is legal.
!
  implicit none

  logical halham_n_check
  integer ( kind = 4 ) n

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_N_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i12)' ) '  N = ', n
    halham_n_check = .false.
  else
    halham_n_check = .true.
  end if

  return
end
function halham_dim_num_check ( dim_num )

!*****************************************************************************80
!
!! HALHAM_DIM_NUM_CHECK checks DIM_NUM for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, logical HALHAM_DIM_NUM_CHECK, is true if DIM_NUM is legal.
!
  implicit none

  logical halham_dim_num_check
  integer ( kind = 4 ) dim_num

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_DIM_NUM_CHECK - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
    halham_dim_num_check = .false.
  else
    halham_dim_num_check = .true.
  end if

  return
end
function halham_seed_check ( dim_num, seed )

!*****************************************************************************80
!
!! HALHAM_SEED_CHECK checks SEED for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the seed vector.
!
!    Output, logical, HALHAM_SEED_CHECK, is true if SEED is legal.
!
  implicit none

  integer ( kind = 4 ) dim_num

  logical halham_seed_check
  integer ( kind = 4 ) seed(dim_num)

  if ( any ( seed(1:dim_num) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_SEED_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some entry of SEED < 0!'
    write ( *, '(a)' ) ' '
    call i4vec_transpose_print ( dim_num, seed, 'SEED:  ' )
    halham_seed_check = .false.
  else
    halham_seed_check = .true.
  end if

  return
end
function halham_step_check ( step )

!*****************************************************************************80
!
!! HALHAM_STEP_CHECK checks STEP for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!
!    Output, logical HALHAM_STEP_CHECK, is true if STEP is legal.
!
  implicit none

  logical halham_step_check
  integer ( kind = 4 ) step

  if ( step < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_STEP_CHECK - Fatal error!'
    write ( *, '(a)' ) '  STEP < 0.'
    write ( *, '(a,i12)' ) '  STEP = ', step
    halham_step_check = .false.
  else
    halham_step_check = .true.
  end if

  return
end
subroutine halham_write ( dim_num, n, step, seed, leap, base, r, file_out_name )

!*****************************************************************************80
!
!! HALHAM_WRITE writes a Halton or Hammersley subsequence to a file.
!
!  Discussion:
!
!    The initial lines of the file are comments, which begin with a
!    '#' character.
!
!    Thereafter, each line of the file contains the DIM_NUM-dimensional
!    components of the next entry of the dataset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the subsequence.
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!    0 <= STEP is required.
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the sequence index for STEP = 0.
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the successive jumps in
!    the sequence.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases.
!
!    Input, real ( kind = 8 ) R(DIM_NUM,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the output file name.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) base(dim_num)
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  integer ( kind = 4 ) mhi
  integer ( kind = 4 ) mlo
  real    ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step
  character ( len = 80 ) string

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file:'
    write ( *, '(a)' ) '  "' // trim ( file_out_name ) // '".'
    stop
  end if

  write ( file_out_unit, '(a)'       ) '#  ' // trim ( file_out_name )
  write ( file_out_unit, '(a)'       ) '#  created by HALHAM_WRITE.F90'
  write ( file_out_unit, '(a)'       ) '#'
  write ( file_out_unit, '(a)'       ) '#'
  write ( file_out_unit, '(a,i12)'    ) '#  DIM_NUM = ', dim_num
  write ( file_out_unit, '(a,i12)'    ) '#  N =    ', n
  write ( file_out_unit, '(a,i12)'    ) '#  STEP = ', step
  do mlo = 1, dim_num, 5
    mhi = min ( mlo + 5 - 1, dim_num )
    if ( mlo == 1 ) then
      write ( file_out_unit, '(a,5i12)' ) '#  SEED = ', seed(mlo:mhi)
    else
      write ( file_out_unit, '(a,5i12)' ) '#         ', seed(mlo:mhi)
    end if
  end do
  do mlo = 1, dim_num, 5
    mhi = min ( mlo + 5 - 1, dim_num )
    if ( mlo == 1 ) then
      write ( file_out_unit, '(a,5i12)' ) '#  LEAP = ', leap(mlo:mhi)
    else
      write ( file_out_unit, '(a,5i12)' ) '#         ', leap(mlo:mhi)
    end if
  end do
  do mlo = 1, dim_num, 5
    mhi = min ( mlo + 5 - 1, dim_num )
    if ( mlo == 1 ) then
      write ( file_out_unit, '(a,5i12)' ) '#  BASE = ', base(mlo:mhi)
    else
      write ( file_out_unit, '(a,5i12)' ) '#         ', base(mlo:mhi)
    end if
  end do
  write ( file_out_unit, '(a,g14.6)' ) '#  EPSILON (unit roundoff ) = ', &
    epsilon ( r(1,1) )
  write ( file_out_unit, '(a)'      ) '#'

  write ( string, '(a,i3,a)' ) '(', dim_num, '(2x,f10.6))'

  do j = 1, n
    write ( file_out_unit, string ) r(1:dim_num,j)
  end do

  close ( unit = file_out_unit )

  return
end
subroutine hammersley ( dim_num, r )

!*****************************************************************************80
!
!! HAMMERSLEY computes the next element in a leaped Hammersley subsequence.
!
!  Discussion:
!
!    The DIM_NUM-dimensional Hammersley sequence is really DIM_NUM separate
!    sequences, each generated by a particular base.  If the base is
!    greater than 1, a standard 1-dimensional
!    van der Corput sequence is generated.  But if the base is
!    negative, this is a signal that the much simpler sequence J/(-BASE)
!    is to be generated.  For the standard Hammersley sequence, the
!    first spatial coordinate uses a base of (-N), and subsequent
!    coordinates use bases of successive primes (2, 3, 5, 7, 11, ...).
!    This program allows the user to specify any combination of bases,
!    included nonprimes and repeated values.
!
!    This routine selects elements of a "leaped" subsequence of the
!    Hammersley sequence.  The subsequence elements are indexed by a
!    quantity called STEP, which starts at 0.  The STEP-th subsequence
!    element is simply element
!
!      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
!
!    of the original Hammersley sequence.
!
!
!    This routine "hides" a number of input arguments.  To specify these
!    arguments explicitly, use I4_TO_HAMMERSLEY instead.
!
!    All the arguments have default values.  However, if you want to
!    examine or change them, you may call the appropriate routine first.
!
!    * DIM_NUM, the spatial dimension,
!      Default: DIM_NUM = 1;
!      Required: 1 <= DIM_NUM is required.
!
!    * STEP, the subsequence index.
!      Default: STEP = 0.
!      Required: 0 <= STEP.
!
!    * SEED(1:DIM_NUM), the Hammersley sequence element for  STEP = 0.
!      Default SEED = (0, 0, ... 0).
!      Required: 0 <= SEED(1:DIM_NUM).
!
!    * LEAP(1:DIM_NUM), the succesive jumps in the sequence.
!      Default: LEAP = (1, 1, ..., 1).
!      Required: 1 <= LEAP(1:DIM_NUM).
!
!    * BASE(1:DIM_NUM), the bases.
!      Default: BASE = (2, 3, 5, 7, 11, ... ) or ( -N, 2, 3, 5, 7, 11,...)
!      if N is known.
!      Required: 0, 1 /= BASE(1:DIM_NUM).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J M Hammersley,
!    Monte Carlo methods for solving multivariable problems,
!    Proceedings of the New York Academy of Science,
!    Volume 86, 1960, pages 844-874.
!
!    Ladislav Kocis and William Whiten,
!    Computational Investigations of Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 2, 1997, pages 266-294.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) R(DIM_NUM), the next element of the
!    leaped Hammersley subsequence.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) leap(dim_num)
  real    ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step
  integer ( kind = 4 ) value(1)

  value(1) = dim_num
  call hammersley_memory ( 'SET', 'DIM_NUM', 1, value )
  call hammersley_memory ( 'GET', 'STEP', 1, value )
  step = value(1)
  call hammersley_memory ( 'GET', 'SEED', dim_num, seed )
  call hammersley_memory ( 'GET', 'LEAP', dim_num, leap )
  call hammersley_memory ( 'GET', 'BASE', dim_num, base )

  call i4_to_hammersley ( dim_num, step, seed, leap, base, r )

  value(1) = 1
  call hammersley_memory ( 'INC', 'STEP', 1, value )

  return
end
function hammersley_base_check ( dim_num, base )

!*****************************************************************************80
!
!! HAMMERSLEY_BASE_CHECK checks BASE for a Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases.
!
!    Output, logical, HAMMERSLEY_BASE_CHECK, is true if BASE is legal.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) base(dim_num)
  logical hammersley_base_check

  if ( any ( base(1:dim_num) == 0 ) .or. any ( base(1:dim_num) == 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HAMMERSLEY_BASE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some entry of BASE is 0 or 1!'
    write ( *, '(a)' ) ' '
    call i4vec_transpose_print ( dim_num, base, 'BASE:  ' )
    hammersley_base_check = .false.
  else
    hammersley_base_check = .true.
  end if

  return
end
subroutine hammersley_base_get ( base )

!*****************************************************************************80
!
!! HAMMERSLEY_BASE_GET gets the base vector for a leaped Hammersley subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BASE(DIM_NUM), the bases.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) base(*)
  integer ( kind = 4 ) value(1)

  call hammersley_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  call hammersley_memory ( 'GET', 'BASE', dim_num, base )

  return
end
subroutine hammersley_base_set ( base )

!*****************************************************************************80
!
!! HAMMERSLEY_BASE_SET sets the base vector for a leaped Hammersley subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases.
!
  implicit none

  integer ( kind = 4 ) base(*)
  logical hammersley_base_check
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) value(1)

  call hammersley_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  if ( .not. hammersley_base_check ( dim_num, base ) ) then
    stop
  end if

  call hammersley_memory ( 'SET', 'BASE', dim_num, base )

  return
end
subroutine hammersley_leap_get ( leap )

!*****************************************************************************80
!
!! HAMMERSLEY_LEAP_GET gets the leap vector for a leaped Hammersley subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) LEAP(DIM_NUM), the successive jumps in
!    the sequence.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) leap(*)
  integer ( kind = 4 ) value(1)

  call hammersley_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  call hammersley_memory ( 'GET', 'LEAP', dim_num, leap )

  return
end
subroutine hammersley_leap_set ( leap )

!*****************************************************************************80
!
!! HAMMERSLEY_LEAP_SET sets the leap vector for a leaped Hammersley subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the successive jumps in
!    the sequence.
!
  implicit none

  logical halham_leap_check
  integer ( kind = 4 ) leap(*)
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) value(1)

  call hammersley_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  if ( .not. halham_leap_check ( dim_num, leap ) ) then
    stop
  end if

  call hammersley_memory ( 'SET', 'LEAP', dim_num, leap )

  return
end
subroutine hammersley_memory ( action, name, dim_num, value )

!*****************************************************************************80
!
!! HAMMERSLEY_MEMORY holds data associated with a leaped Hammersley subsequence.
!
!  Discussion:
!
!    If you're going to define a new problem, it's important that
!    you set the value of DIM_NUM before setting the values of BASE,
!    LEAP or SEED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the desired action.
!    'GET' means get the value of a particular quantity.
!    'SET' means set the value of a particular quantity.
!    'INC' means increment the value of a particular quantity.
!          (Only SEED and STEP can be incremented.)
!
!    Input, character ( len = * ) NAME, the name of the quantity.
!    'BASE' means the base vector.
!    'LEAP' means the leap vector.
!    'DIM_NUM' means the spatial dimension.
!    'SEED' means the seed vector.
!    'STEP' means the step.
!
!    Input/output, integer ( kind = 4 ) DIM_NUM, the dimension of the quantity.
!    If ACTION is 'SET' and NAME is 'BASE', then DIM_NUM is input, and
!    is the number of entries in VALUE to be put into BASE.
!
!    Input/output, integer ( kind = 4 ) VALUE(DIM_NUM), contains a value.
!    If ACTION is 'SET', then on input, VALUE contains values to be assigned
!    to the internal variable.
!    If ACTION is 'GET', then on output, VALUE contains the values of
!    the specified internal variable.
!    If ACTION is 'INC', then on input, VALUE contains the increment to
!    be added to the specified internal variable.
!
  implicit none

  character ( len = * ) action
  integer ( kind = 4 ), allocatable, save, dimension ( : ) :: base
  logical, save :: first_call = .true.
  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable, save, dimension ( : ) :: leap
  character ( len = * ) name
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), save :: dim_num_save = 0
  integer ( kind = 4 ) prime
  integer ( kind = 4 ), allocatable, save, dimension ( : ) :: seed
  integer ( kind = 4 ), save :: step = 0
  integer ( kind = 4 ) value(*)

  if ( first_call ) then
    dim_num_save = 1
    allocate ( base(dim_num_save) )
    allocate ( leap(dim_num_save) )
    allocate ( seed(dim_num_save) )
    base(1) = 2
    leap(1) = 1
    seed(1) = 0
    step = 0
    first_call = .false.
  end if
!
!  If this is a SET DIM_NUM call, and the input value of DIM_NUM
!  differs from the internal value, discard all old information.
!
  if ( action(1:1) == 'S' .or. action(1:1) == 's') then
    if ( name == 'DIM_NUM' .or. name == 'dim_num' ) then
      if ( dim_num_save /= value(1) ) then
        deallocate ( base )
        deallocate ( leap )
        deallocate ( seed )
        dim_num_save = value(1)
        allocate ( base(dim_num_save) )
        allocate ( leap(dim_num_save) )
        allocate ( seed(dim_num_save) )
        do i = 1, dim_num_save
          base(i) = prime ( i )
        end do
        leap(1:dim_num_save) = 1
        seed(1:dim_num_save) = 0
      end if
    end if
  end if
!
!  Set
!
  if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'BASE' .or. name == 'base' ) then

      if ( dim_num_save /= dim_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HAMMERSLEY_MEMORY - Fatal error!'
        write ( *, '(a)' ) '  Internal and input values of DIM_NUM disagree'
        write ( *, '(a)' ) '  while setting BASE.'
        stop
      end if

      base(1:dim_num) = value(1:dim_num)

    else if ( name == 'LEAP' .or. name == 'leap' ) then

      if ( dim_num_save /= dim_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HAMMERSLEY_MEMORY - Fatal error!'
        write ( *, '(a)' ) '  Internal and input values of DIM_NUM disagree'
        write ( *, '(a)' ) '  while setting LEAP.'
        stop
      end if

      leap(1:dim_num) = value(1:dim_num)

    else if ( name == 'DIM_NUM' .or. name == 'dim_num' ) then

      dim_num_save = value(1)

    else if ( name == 'SEED' .or. name == 'seed' ) then

      if ( dim_num_save /= dim_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HAMMERSLEY_MEMORY - Fatal error!'
        write ( *, '(a)' ) '  Internal and input values of DIM_NUM disagree'
        write ( *, '(a)' ) '  while setting SEED.'
        stop
      end if

      seed(1:dim_num) = value(1:dim_num)

    else if ( name == 'STEP' .or. name == 'step' ) then

      if ( value(1) < 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HAMMERSLEY_MEMORY - Fatal error!'
        write ( *, '(a)' ) '  Input value of STEP < 0.'
        stop
      end if

      step = value(1)

    end if
!
!  Get
!
  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'BASE' .or. name == 'base' ) then

      value(1:dim_num_save) = base(1:dim_num_save)

    else if ( name == 'LEAP' .or. name == 'leap' ) then

      value(1:dim_num_save) = leap(1:dim_num_save)

    else if ( name == 'DIM_NUM' .or. name == 'dim_num' ) then

      value(1) = dim_num_save

    else if ( name == 'SEED' .or. name == 'seed' ) then

      value(1:dim_num_save) = seed(1:dim_num_save)

    else if ( name == 'STEP' .or. name == 'step' ) then

      value(1) = step

    end if
!
!  Increment
!
  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == 'SEED' .or. name == 'seed' ) then
      if ( dim_num == 1 ) then
        seed(1:dim_num_save) = seed(1:dim_num_save) + value(1)
      else
        seed(1:dim_num_save) = seed(1:dim_num_save) + value(1:dim_num_save)
      end if
    else if ( name == 'STEP' .or. name == 'step' ) then
      step = step + value(1)
    end if

  end if

  return
end
subroutine hammersley_dim_num_get ( dim_num )

!*****************************************************************************80
!
!! HAMMERSLEY_DIM_NUM_GET: spatial dimension, leaped Hammersley subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) value(1)

  call hammersley_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  return
end
subroutine hammersley_dim_num_set ( dim_num )

!*****************************************************************************80
!
!! HAMMERSLEY_DIM_NUM_SET sets spatial dimension, leaped Hammersley subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!    1 <= DIM_NUM is required.
!
  implicit none

  logical halham_dim_num_check
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) value(1)

  if ( .not. halham_dim_num_check ( dim_num ) ) then
    stop
  end if

  value(1) = dim_num
  call hammersley_memory ( 'SET', 'DIM_NUM', 1, value )

  return
end
subroutine hammersley_seed_get ( seed )

!*****************************************************************************80
!
!! HAMMERSLEY_SEED_GET gets the seed vector for a leaped Hammersley subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED(DIM_NUM), the Hammersley sequence
!    index corresponding to STEP = 0.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) seed(*)
  integer ( kind = 4 ) value(1)

  call hammersley_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)
  call hammersley_memory ( 'GET', 'SEED', dim_num, seed )

  return
end
subroutine hammersley_seed_set ( seed )

!*****************************************************************************80
!
!! HAMMERSLEY_SEED_SET sets the seed vector for a leaped Hammersley subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the Hammersley sequence index
!    corresponding to STEP = 0.
!
  implicit none

  logical halham_seed_check
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) seed(*)
  integer ( kind = 4 ) value(1)

  call hammersley_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  if ( .not. halham_seed_check ( dim_num, seed ) ) then
    stop
  end if

  call hammersley_memory ( 'SET', 'SEED', dim_num, seed )

  return
end
subroutine hammersley_sequence ( dim_num, n, r )

!*****************************************************************************80
!
!! HAMMERSLEY_SEQUENCE computes N elements of a leaped Hammersley subsequence.
!
!  Discussion:
!
!    The DIM_NUM-dimensional Hammersley sequence is really DIM_NUM separate
!    sequences, each generated by a particular base.  If the base is
!    greater than 1, a standard 1-dimensional
!    van der Corput sequence is generated.  But if the base is
!    negative, this is a signal that the much simpler sequence J/(-BASE)
!    is to be generated.  For the standard Hammersley sequence, the
!    first spatial coordinate uses a base of (-N), and subsequent
!    coordinates use bases of successive primes (2, 3, 5, 7, 11, ...).
!    This program allows the user to specify any combination of bases,
!    included nonprimes and repeated values.
!
!    This routine selects elements of a "leaped" subsequence of the
!    Hammersley sequence.  The subsequence elements are indexed by a
!    quantity called STEP, which starts at 0.  The STEP-th subsequence
!    element is simply element
!
!      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
!
!    of the original Hammersley sequence.
!
!
!    This routine "hides" a number of input arguments.  To specify these
!    arguments explicitly, use I4_TO_HAMMERSLEY_SEQUENCE instead.
!
!    All the arguments have default values.  However, if you want to
!    examine or change them, you may call the appropriate routine first.
!
!    The arguments that the user may set include:
!
!    * DIM_NUM, the spatial dimension,
!      Default: DIM_NUM = 1;
!      Required: 1 <= DIM_NUM is required.
!
!    * STEP, the subsequence index.
!      Default: STEP = 0.
!      Required: 0 <= STEP.
!
!    * SEED(1:DIM_NUM), the sequence element corresponding to STEP = 0.
!      Default SEED = (0, 0, ... 0).
!      Required: 0 <= SEED(1:DIM_NUM).
!
!    * LEAP(1:DIM_NUM), the succesive jumps in the sequence.
!      Default: LEAP = (1, 1, ..., 1).
!      Required: 1 <= LEAP(1:DIM_NUM).
!
!    * BASE(1:DIM_NUM), the bases.
!      Default: BASE = (2, 3, 5, 7, 11, ... ) or ( -N, 2, 3, 5, 7, 11,...)
!      if N is known.
!      Required: 0, 1 /= BASE(1:DIM_NUM).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J M Hammersley,
!    Monte Carlo methods for solving multivariable problems,
!    Proceedings of the New York Academy of Science,
!    Volume 86, 1960, pages 844-874.
!
!    Ladislav Kocis and William Whiten,
!    Computational Investigations of Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 2, 1997, pages 266-294.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of elements desired.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the next N elements of the
!    leaped Hammersley subsequence.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) leap(dim_num)
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step
  integer ( kind = 4 ) value(1)

  value(1) = dim_num
  call hammersley_memory ( 'SET', 'DIM_NUM', 1, value )
  call hammersley_memory ( 'GET', 'STEP', 1, value )
  step = value(1)
  call hammersley_memory ( 'GET', 'SEED', dim_num, seed )
  call hammersley_memory ( 'GET', 'LEAP', dim_num, leap )
  call hammersley_memory ( 'GET', 'BASE', dim_num, base )

  call i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

  value(1) = n
  call hammersley_memory ( 'INC', 'STEP', 1, value )

  return
end
subroutine hammersley_step_get ( step )

!*****************************************************************************80
!
!! HAMMERSLEY_STEP_GET gets the "step" for a leaped Hammersley subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) STEP, the index of the subsequence element.
!
  implicit none

  integer ( kind = 4 ) step
  integer ( kind = 4 ) value(1)

  call hammersley_memory ( 'GET', 'STEP', 1, value )
  step = value(1)

  return
end
subroutine hammersley_step_set ( step )

!*****************************************************************************80
!
!! HAMMERSLEY_STEP_SET sets the "step" for a leaped Hammersley subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!    0 <= STEP is required.
!
  implicit none

  logical halham_step_check
  integer ( kind = 4 ) step
  integer ( kind = 4 ) value(1)

  if ( .not. halham_step_check ( step ) ) then
    stop
  end if

  value(1) = step
  call hammersley_memory ( 'SET', 'STEP', 1, value )

  return
end
subroutine i4_to_hammersley ( dim_num, step, seed, leap, base, r )

!*****************************************************************************80
!
!! I4_TO_HAMMERSLEY computes one element of a leaped Hammersley subsequence.
!
!  Discussion:
!
!    The DIM_NUM-dimensional Hammersley sequence is really DIM_NUM separate
!    sequences, each generated by a particular base.  If the base is
!    greater than 1, a standard 1-dimensional
!    van der Corput sequence is generated.  But if the base is
!    negative, this is a signal that the much simpler sequence J/(-BASE)
!    is to be generated.  For the standard Hammersley sequence, the
!    first spatial coordinate uses a base of (-N), and subsequent
!    coordinates use bases of successive primes (2, 3, 5, 7, 11, ...).
!    This program allows the user to specify any combination of bases,
!    included nonprimes and repeated values.
!
!    This routine selects elements of a "leaped" subsequence of the
!    Hammersley sequence.  The subsequence elements are indexed by a
!    quantity called STEP, which starts at 0.  The STEP-th subsequence
!    element is simply element
!
!      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
!
!    of the original Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J M Hammersley,
!    Monte Carlo methods for solving multivariable problems,
!    Proceedings of the New York Academy of Science,
!    Volume 86, 1960, pages 844-874.
!
!    Ladislav Kocis and William Whiten,
!    Computational Investigations of Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 2, 1997, pages 266-294.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!    1 <= DIM_NUM is required.
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!    0 <= STEP is required.
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the sequence index corresponding
!    to STEP = 0.
!    0 <= SEED(1:DIM_NUM) is required.
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the successive jumps in
!    the sequence.
!    1 <= LEAP(1:DIM_NUM) is required.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases.
!
!    Output, real ( kind = 8 ) R(DIM_NUM), the STEP-th element of the leaped
!    Hammersley subsequence.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) base(dim_num)
  real ( kind = 8 ) base_inv
  integer ( kind = 4 ) digit
  real ( kind = 8 ) :: fiddle = 1.0D+00
  logical halham_leap_check
  logical halham_dim_num_check
  logical halham_seed_check
  logical halham_step_check
  logical hammersley_base_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) leap(dim_num)
  real ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) seed2
  integer ( kind = 4 ) step
!
!  Check the input.
!
  if ( .not. halham_dim_num_check ( dim_num ) ) then
    stop
  end if

  if ( .not. halham_step_check ( step ) ) then
    stop
  end if

  if ( .not. halham_seed_check ( dim_num, seed ) ) then
    stop
  end if

  if ( .not. halham_leap_check ( dim_num, leap ) ) then
    stop
  end if

  if ( .not. hammersley_base_check ( dim_num, base ) ) then
    stop
  end if
!
!  Calculate the data.
!
  do i = 1, dim_num

    if ( 1 < base(i) ) then

      seed2 = seed(i) + step * leap(i)

      r(i) = 0.0D+00

      base_inv = real ( 1.0D+00, kind = 8 ) / real ( base(i), kind = 8 )

      do while ( seed2 /= 0 )
        digit = mod ( seed2, base(i) )
        r(i) = r(i) + real ( digit, kind = 8 ) * base_inv
        base_inv = base_inv / real ( base(i), kind = 8 )
        seed2 = seed2 / base(i)
      end do
!
!  In the following computation, the value of FIDDLE can be:
!
!    0,   for the sequence 0/N, 1/N, ..., N-1/N
!    1,   for the sequence 1/N, 2/N, ..., N/N
!    1/2, for the sequence 1/(2N), 3/(2N), ..., (2*N-1)/(2N)
!
    else if ( base(i) <= -1 ) then

      seed2 = seed(i) + step * leap(i)

      seed2 = mod ( seed2, abs ( base(i) ) )

      r(i) = ( real ( seed2, kind = 8 ) + fiddle ) &
             / real ( -base(i), kind = 8 )

    end if

  end do

  return
end
subroutine i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

!*****************************************************************************80
!
!! I4_TO_HAMMERSLEY_SEQUENCE: N elements of a leaped Hammersley subsequence.
!
!  Discussion:
!
!    The DIM_NUM-dimensional Hammersley sequence is really DIM_NUM separate
!    sequences, each generated by a particular base.  If the base is
!    greater than 1, a standard 1-dimensional
!    van der Corput sequence is generated.  But if the base is
!    negative, this is a signal that the much simpler sequence J/(-BASE)
!    is to be generated.  For the standard Hammersley sequence, the
!    first spatial coordinate uses a base of (-N), and subsequent
!    coordinates use bases of successive primes (2, 3, 5, 7, 11, ...).
!    This program allows the user to specify any combination of bases,
!    included nonprimes and repeated values.
!
!    This routine selects elements of a "leaped" subsequence of the
!    Hammersley sequence.  The subsequence elements are indexed by a
!    quantity called STEP, which starts at 0.  The STEP-th subsequence
!    element is simply element
!
!      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
!
!    of the original Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J M Hammersley,
!    Monte Carlo methods for solving multivariable problems,
!    Proceedings of the New York Academy of Science,
!    Volume 86, 1960, pages 844-874.
!
!    Ladislav Kocis and William Whiten,
!    Computational Investigations of Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 2, 1997, pages 266-294.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!    1 <= DIM_NUM is required.
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence.
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!    0 <= STEP is required.
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the sequence index corresponding
!    to STEP = 0.
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the succesive jumps in
!    the sequence.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the next N elements of the
!    leaped Hammersley subsequence, beginning with element STEP.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) base(dim_num)
  real ( kind = 8 ) base_inv
  integer ( kind = 4 ) digit(n)
  real ( kind = 8 ) :: fiddle = 1.0D+00
  logical halham_leap_check
  logical halham_dim_num_check
  logical halham_seed_check
  logical halham_step_check
  logical hammersley_base_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) seed2(n)
  integer ( kind = 4 ) step
!
!  Check the input.
!
  if ( .not. halham_dim_num_check ( dim_num ) ) then
    stop
  end if

  if ( .not. halham_step_check ( step ) ) then
    stop
  end if

  if ( .not. halham_seed_check ( dim_num, seed ) ) then
    stop
  end if

  if ( .not. halham_leap_check ( dim_num, leap ) ) then
    stop
  end if

  if ( .not. hammersley_base_check ( dim_num, base ) ) then
    stop
  end if
!
!  Calculate the data.
!
  do i = 1, dim_num

    if ( 1 < base(i) ) then

      do j = 1, n
        seed2(j) = seed(i) + ( step + j - 1 ) * leap(i)
      end do

      r(i,1:n) = 0.0D+00

      base_inv = real ( 1.0D+00, kind = 8 ) / real ( base(i), kind = 8 )

      do while ( any ( seed2(1:n) /= 0 ) )
        digit(1:n) = mod ( seed2(1:n), base(i) )
        r(i,1:n) = r(i,1:n) + real ( digit(1:n), kind = 8 ) * base_inv
        base_inv = base_inv / real ( base(i), kind = 8 )
        seed2(1:n) = seed2(1:n) / base(i)
      end do
!
!  In the following computation, the value of FIDDLE can be:
!
!    0,   for the sequence 0/N, 1/N, ..., N-1/N
!    1,   for the sequence 1/N, 2/N, ..., N/N
!    1/2, for the sequence 1/(2N), 3/(2N), ..., (2*N-1)/(2N)
!
    else if ( base(i) <= -1 ) then

      do j = 1, n
        seed2(j) = seed(i) + ( step + j - 1 ) * leap(i)
      end do

      seed2(1:n) = mod ( seed2(1:n), abs ( base(i) ) )

      r(i,1:n) = ( real ( seed2(1:n), kind = 8 ) + fiddle ) &
                 / real ( -base(i), kind = 8 )

    end if

  end do

  return
end
subroutine i4vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
!
!  Example:
!
!    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
!    TITLE = 'My vector:  '
!
!    My vector:      1    2    3    4    5
!                    6    7    8    9   10
!                   11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = 11 ) string
  character ( len = * ) title
  integer ( kind = 4 ) title_len

  if ( 0 < len ( title ) ) then

    title_len = len ( title )

    write ( string, '(a,i3,a)' ) '(', title_len, 'x,5i12)'

    do ilo = 1, n, 5
      ihi = min ( ilo + 5 - 1, n )
      if ( ilo == 1 ) then
        write ( *, '(a, 5i12)' ) title, a(ilo:ihi)
      else
        write ( *, string      )        a(ilo:ihi)
      end if
    end do

  else

    do ilo = 1, n, 5
      ihi = min ( ilo + 5 - 1, n )
      write ( *, '(5i12)' ) a(ilo:ihi)
    end do

  end if

  return
end
function prime ( n )

!*****************************************************************************80
!
!! PRIME returns any of the first PRIME_MAX prime numbers.
!
!  Discussion:
!
!    PRIME_MAX is 1600, and the largest prime stored is 13499.
!
!    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964, pages 870-873.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 95-98.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the desired prime number.
!    In general, is should be true that 0 <= N <= PRIME_MAX.
!    N = -1 returns PRIME_MAX, the index of the largest prime available.
!    N = 0 is legal, returning PRIME = 1.
!
!    Output, integer ( kind = 4 ) PRIME, the N-th prime.  If N is out of range,
!    PRIME is returned as -1.
!
  implicit none

  integer ( kind = 4 ), parameter :: prime_max = 1600

  integer ( kind = 4 ), save :: icall = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( prime_max ) :: npvec
  integer ( kind = 4 ) prime

  if ( icall == 0 ) then

    icall = 1

    npvec(1:100) = (/ &
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /)

    npvec(101:200) = (/ &
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /)

    npvec(201:300) = (/ &
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /)

    npvec(301:400) = (/ &
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /)

    npvec(401:500) = (/ &
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /)

    npvec(501:600) = (/ &
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /)

    npvec(601:700) = (/ &
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /)

    npvec(701:800) = (/ &
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /)

    npvec(801:900) = (/ &
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /)

    npvec(901:1000) = (/ &
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)

    npvec(1001:1100) = (/ &
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /)

    npvec(1101:1200) = (/ &
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /)

    npvec(1201:1300) = (/ &
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /)

    npvec(1301:1400) = (/ &
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949, &
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /)

    npvec(1401:1500) = (/ &
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, &
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, &
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, &
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, &
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, &
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, &
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, &
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, &
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, &
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /)

   npvec(1501:1600) = (/ &
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, &
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, &
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, &
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, &
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, &
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, &
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, &
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, &
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, &
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 /)

  end if

  if ( n == -1 ) then
    prime = prime_max
  else if ( n == 0 ) then
    prime = 1
  else if ( n <= prime_max ) then
    prime = npvec(n)
  else
    prime = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRIME - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i8)' ) '  N should be between 1 and PRIME_MAX =', prime_max
    stop
  end if

  return
end
subroutine u1_to_sphere_unit_2d ( u, x )

!*****************************************************************************80
!
!! U1_TO_SPHERE_UNIT_2D maps a point in the unit interval to the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U, a point in the unit interval.
!
!    Output, real ( kind = 8 ) X(2), the corresponding point on the circle.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) x(2)

  angle = real ( 2.0D+00, kind = 8 ) * pi * u

  x(1) = cos ( angle )
  x(2) = sin ( angle )

  return
end
subroutine u2_to_ball_unit_2d ( u, x )

!*****************************************************************************80
!
!! U2_TO_BALL_UNIT_2D maps points from the unit box to the unit ball in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U(2), a point in the unit square.
!
!    Output, real ( kind = 8 ) X(2), the corresponding point in the
!    unit ball.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) theta
  real ( kind = 8 ) u(2)
  real ( kind = 8 ) x(2)

  r = sqrt ( u(1) )
  theta = real ( 2.0D+00, kind = 8 ) * pi * u(2)

  x(1) = r * cos ( theta )
  x(2) = r * sin ( theta )

  return
end
subroutine u2_to_sphere_unit_3d ( u, x )

!*****************************************************************************80
!
!! U2_TO_SPHERE_UNIT_3D maps a point in the unit box onto the unit sphere in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U(2), the point in the unit box.
!
!    Output, real ( kind = 8 ) X(3), the corresponding point on the unit sphere.
!
  implicit none

  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) phi
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) u(2)
  real ( kind = 8 ) vdot
  real ( kind = 8 ) x(3)
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  vdot = real ( 2.0D+00, kind = 8 ) * u(1) - real ( 1.0D+00, kind = 8 )

  phi = arc_cosine ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  theta = real ( 2.0D+00, kind = 8 ) * pi * u(2)

  x(1) = cos ( theta ) * sin ( phi )
  x(2) = sin ( theta ) * sin ( phi )
  x(3) = cos ( phi )

  return
end
subroutine u3_to_ball_unit_3d ( u, x )

!*****************************************************************************80
!
!! U3_TO_BALL_UNIT_3D maps points from the unit box to the unit ball in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U(3), a point in the unit box in 3D.
!
!    Output, real ( kind = 8 ) X(3), the corresponding point in the
!    unit ball in 3D.
!
  implicit none

  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) phi
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) theta
  real ( kind = 8 ) u(3)
  real ( kind = 8 ) vdot
  real ( kind = 8 ) x(3)
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  vdot = real ( 2.0D+00, kind = 8 ) * u(1) - real ( 1.0D+00, kind = 8 )

  phi = arc_cosine ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  theta = real ( 2.0D+00, kind = 8 ) * pi * u(2)
!
!  Pick a random radius R.
!
  r = u(3)**( real ( 1.0D+00, kind = 8 ) / real ( 3.0D+00, kind = 8 ) )

  x(1) = r * cos ( theta ) * sin ( phi )
  x(2) = r * sin ( theta ) * sin ( phi )
  x(3) = r * cos ( phi )

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
