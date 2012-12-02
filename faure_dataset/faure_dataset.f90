program main

!*****************************************************************************80
!
!! MAIN is the main program for FAURE_DATASET.
!
!  Discussion:
!
!    FAURE_DATASET generates a Faure dataset and writes it to a file.
!
!  Usage:
!
!    faure_dataset m n skip
!
!    where
!
!    * M, the spatial dimension,
!    * N, the number of points to generate,
!    * SKIP, the number of initial values to skip.
!
!    creates an M by N dataset and writes it to the
!    file "faure_M_N.txt".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 )  arg_num
  integer   ( kind = 4 )  base
  integer   ( kind = 4 )  iarg
  integer   ( kind = 4 )  iargc
  integer   ( kind = 4 )  ierror
  integer   ( kind = 4 )  ios
  integer   ( kind = 4 )  last
  integer   ( kind = 4 )  m
  integer   ( kind = 4 )  n
  character ( len = 255 ) output_filename
  integer   ( kind = 4 )  prime_ge
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r
  integer   ( kind = 4 )  skip
  integer   ( kind = 4 )  skip_suggest
  character ( len = 255 ) string

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FAURE_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Generate a Faure dataset.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get M.
!
  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, m, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter M, the spatial dimension:'
    read ( *, *, iostat = ios ) m
  end if
!
!  Get N.
!
  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, n, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter N, the number of points to generate:'

    read ( *, *, iostat = ios ) n

  end if
!
!  Get the SKIP.
!
  base = prime_ge ( m )
  skip_suggest = base**4 - 1

  if ( 3 <= arg_num ) then
    iarg = 3
    call getarg ( iarg, string )
    call s_to_i4 ( string, skip, ierror, last )
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter SKIP, the number of initial points to skip:'
    write ( *, '(a)' ) '  (0 is a common choice).'
    write ( *, '(a,i6,a)' ) '  (', skip_suggest, ' is a recommended choice).'

    read ( *, *, iostat = ios ) skip

  end if
!
!  Compute the data.
!
  allocate ( r(1:m,1:n) )

  call faure_generate ( m, n, skip, base, r )
!
!  Write the data to a file.
!
  write ( output_filename, '(a,i2.2,a,i5.5,a)' ) &
    'faure_', m, '_', n, '.txt'

  call r8mat_write ( output_filename, m, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Faure data was written to the file "' &
     // trim ( output_filename ) // '".'
!
!  Free memory.
!
  deallocate ( r )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FAURE_DATASET:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine binomial_table ( qs, m, n, coef )

!*****************************************************************************80
!
!! BINOMIAL_TABLE computes a table of bionomial coefficients MOD QS.
!
!  Discussion:
!
!    Thanks to Michael Baudin for pointing out an error in a previous
!    version of this function, 07 December 2009.
!
!  Modified:
!
!    07 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) QS, the base for the MOD operation.
!
!    Input, integer ( kind = 4 ) M, N, the limits of the binomial table.
!
!    Output, integer ( kind = 4 ) COEF(0:M,0:N), the table of binomial 
!    coefficients modulo QS.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) coef(0:m,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) qs

  coef(0:m,0:n) = 0

  coef(0:m,0) = 1

  do j = 1, min ( m, n )
    coef(j,j) = 1
  end do

  do j = 1, n
    do i = j + 1, m
      coef(i,j) = mod ( coef(i-1,j) + coef(i-1,j-1), qs )
    end do
  end do

  return
end
subroutine faure ( dim_num, seed, quasi )

!*****************************************************************************80
!
!! FAURE generates a new quasirandom Faure vector with each call.
!
!  Discussion:
!
!    This routine implements Faure's method of computing
!    quasirandom numbers.  It is a merging and adaptation of
!    Bennett Fox's routines INFAUR and GOFAUR from ACM TOMS 
!    Algorithm 647.
!
!  Modified:
!
!    16 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henri Faure,
!    Discrepance de suites associees a un systeme de numeration
!    (en dimension s),
!    Acta Arithmetica,
!    Volume 41, 1982, pages 337-351, especially page 342.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom 
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension, which should be
!    at least 2.
!
!    Input/output, integer ( kind = 4 ) SEED, the seed, which can be used to index
!    the values.  On first call, set the input value of SEED to be 0
!    or negative.  The routine will automatically initialize data,
!    and set SEED to a new value.  Thereafter, to compute successive
!    entries of the sequence, simply call again without changing
!    SEED.  On the first call, if SEED is negative, it will be set
!    to a positive value that "skips over" an early part of the sequence
!    (This is recommended for better results).
!
!    Output, real ( kind = 8 ) QUASI(DIM_NUM), the next quasirandom vector.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ), save, allocatable, dimension ( :, : ) :: coef
  integer ( kind = 4 ) hisum
  integer ( kind = 4 ), save :: hisum_save = -1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_i4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ktemp
  integer ( kind = 4 ) ltemp
  integer ( kind = 4 ) mtemp
  integer ( kind = 4 ) prime_ge
  integer ( kind = 4 ), save :: qs = -1
  real ( kind = 8 ) quasi(dim_num)
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save, allocatable, dimension ( : ) :: ytemp
  integer ( kind = 4 ) ztemp
!
!  Initialization required or requested?
!
  if ( qs <= 0 .or. seed <= 0 ) then

    qs = prime_ge ( dim_num )

    if ( qs < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAURE - Fatal error!'
      write ( *, '(a)' ) '  PRIME_GE failed.'
      stop
    end if

    hisum_save = -1

  end if
!
!  If SEED < 0, reset for recommended initial skip.
!
  if ( seed < 0 ) then

    hisum = 3
    seed = qs**( hisum + 1 ) - 1

  elseif ( seed == 0 ) then

    hisum = 0

  else

    hisum = i4_log_i4 ( seed, qs )

  end if
!
!  Is it necessary to recompute the coefficient table?
!
  if ( hisum_save < hisum ) then

    if ( allocated ( coef ) ) then
      deallocate ( coef )
    end if

    if ( allocated ( ytemp ) ) then
      deallocate ( ytemp )
    end if

    hisum_save = hisum

    allocate ( coef(0:hisum,0:hisum) )
    allocate ( ytemp(0:hisum) )

    call binomial_table ( qs, hisum, hisum, coef )

  end if
!
!  Find QUASI(1) using the method of Faure.
!
!  SEED has a representation in base QS of the form: 
!
!    Sum ( 0 <= J <= HISUM ) YTEMP(J) * QS**J
!
!  We now compute the YTEMP(J)'s.
!
  ktemp = qs**( hisum + 1 )
  ltemp = seed
  do i = hisum, 0, -1
    ktemp = ktemp / qs
    mtemp = mod ( ltemp, ktemp )
    ytemp(i) = ( ltemp - mtemp ) / ktemp
    ltemp = mtemp
  end do
!
!  QUASI(K) has the form
!
!    Sum ( 0 <= J <= HISUM ) YTEMP(J) / QS**(J+1)
!
!  Compute QUASI(1) using nested multiplication.
!
  r = real ( ytemp(hisum), kind = 8 )
  do i = hisum - 1, 0, -1
    r = real ( ytemp(i), kind = 8 ) + r / real ( qs, kind = 8 )
  end do

  quasi(1) = r / real ( qs, kind = 8 )
!
!  Find components QUASI(2:DIM_NUM) using the Faure method.
!
  do k = 2, dim_num

    quasi(k) = 0.0D+00
    r = 1.0D+00 / real ( qs, kind = 8 )

    do j = 0, hisum

      ztemp = dot_product ( ytemp(j:hisum), coef(j:hisum,j) )
!
!  New YTEMP(J) is:
!
!    Sum ( J <= I <= HISUM ) ( old ytemp(i) * binom(i,j) ) mod QS.
!
      ytemp(j) = mod ( ztemp, qs )
      quasi(k) = quasi(k) + real ( ytemp(j), kind = 8 ) * r
      r = r / real ( qs, kind = 8 )

    end do

  end do
!
!  Update SEED.
!
  seed = seed + 1

  return
end
subroutine faure_generate ( dim_num, n, skip, base, r )

!*****************************************************************************80
!
!! FAURE_GENERATE generates a Faure dataset.
!
!  Modified:
!
!    19 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points to generate.
!
!    Input, integer ( kind = 4 ) SKIP, the number of initial points to skip.
!
!    Output, integer ( kind = 4 ) BASE, the base used for the sequence.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) j
  integer ( kind = 4 ) prime_ge
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) skip

  base = prime_ge ( dim_num )

  do j = 1, n
    seed = skip + j - 1
    call faure ( dim_num, seed, r(1:dim_num,j) )
  end do

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
  logical              lopen

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
function i4_log_i4 ( i4, j4 )

!*****************************************************************************80
!
!! I4_LOG_I4 returns the logarithm of an I4 to an I4 base.
!
!  Discussion:
!
!    Only the integer part of the logarithm is returned.
!
!    If
!
!      K4 = I4_LOG_J4 ( I4, J4 ),
!
!    then we ordinarily have
!
!      J4^(K4-1) < I4 <= J4^K4.
!
!    The base J4 should be positive, and at least 2.  If J4 is negative,
!    a computation is made using the absolute value of J4.  If J4 is
!    -1, 0, or 1, the logarithm is returned as 0.
!
!    The number I4 should be positive and at least 2.  If I4 is negative,
!    a computation is made using the absolute value of I4.  If I4 is
!    -1, 0, or 1, then the logarithm is returned as 0.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!    I4  J4  K4
!
!     0   3   0
!     1   3   0
!     2   3   0
!     3   3   1
!     4   3   1
!     8   3   1
!     9   3   2
!    10   3   2
!
!  Modified:
!
!    09 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I4, the number whose logarithm is desired.
!
!    Input, integer ( kind = 4 ) J4, the base of the logarithms.
!
!    Output, integer ( kind = 4 ) I4_LOG_I4, the integer part of the logarithm
!    base abs(J4) of abs(I4).
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_abs
  integer ( kind = 4 ) i4_log_i4
  integer ( kind = 4 ) j4
  integer ( kind = 4 ) j4_abs
  integer ( kind = 4 ) value

  value = 0

  i4_abs = abs ( i4 )

  if ( 2 <= i4_abs ) then

    j4_abs = abs ( j4 )

    if ( 2 <= j4_abs ) then

      do while ( j4_abs <= i4_abs )
        i4_abs = i4_abs / j4_abs
        value = value + 1
      end do

    end if

  end if

  i4_log_i4 = value

  return
end
function prime_ge ( n )

!*****************************************************************************80
!
!! PRIME_GE returns the smallest prime greater than or equal to N.
!
!  Example:
!
!    N     PRIME_GE
!
!    -10    2
!      1    2
!      2    2
!      3    3
!      4    5
!      5    5
!      6    7
!      7    7
!      8   11
!      9   11
!     10   11
!
!  Modified:
!
!    09 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number to be bounded.
!
!    Output, integer ( kind = 4 ) PRIME_GE, the smallest prime number that is 
!    greater than or equal to N.  However, if N is larger than 12,553, (the
!    largest prime stored), then PRIME_GE is returned as -1.
!
  implicit none

  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  integer ( kind = 4 ) i_mid
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p_hi
  integer ( kind = 4 ) p_lo
  integer ( kind = 4 ) p_mid
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) prime_ge

  if ( n <= 2 ) then
    prime_ge = 2
    return
  end if

  i_lo = 1
  p_lo = prime(i_lo)
  i_hi = prime(-1)
  p_hi = prime(i_hi)

  if ( p_hi < n ) then
    prime_ge = -p_hi
    return
  end if

  do

    if ( i_lo + 1 == i_hi ) then
      prime_ge = p_hi
      return
    end if

    i_mid = ( i_lo + i_hi ) / 2
    p_mid = prime(i_mid)

    if ( p_mid < n ) then
      i_lo = i_mid
      p_lo = p_mid
    else if ( n <= p_mid ) then
      i_hi = i_mid
      p_hi = p_mid
    end if

  end do

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
    write ( *, '(a,i6)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i6)' ) '  N should be between 1 and PRIME_MAX =', prime_max
    stop
  end if

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
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
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_filename
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real      ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine s_to_i4 ( s, value, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Licensing:
!
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) VALUE, the integer value read from the string.
!    If the string is blank, then VALUE will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters 
!    of S used to make the integer.
!
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) length
  character ( len = * )  s
  integer   ( kind = 4 ) state
  integer   ( kind = 4 ) value

  value = 0
  ierror = 0
  length = 0

  state = 0
  isgn = 1

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  STATE = 0, haven't read anything.
!
    if ( state == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        state = 1
        isgn = -1
      else if ( c == '+' ) then
        state = 1
        isgn = +1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 2, have read at least one digit, expecting more.
!
    else if ( state == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + iachar ( c ) - iachar ( '0' )

      else

        value = isgn * value
        ierror = 0
        length = i - 1
        return

      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( state == 2 ) then

    value = isgn * value
    ierror = 0
    length = len_trim ( s )

  else

    value = 0
    ierror = 1
    length = 0

  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  character ( len = 8 )  date
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  character ( len = 10 ) time
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y
  character ( len = 5 )  zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
