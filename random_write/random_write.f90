program main

!*****************************************************************************80
!
!! MAIN is the main program for RANDOM_WRITE.
!
!  Discussion:
!
!    This program assumes that the following functions are available:
!
!    AND ( )
!    IRAND ( )
!    LSHIFT ( )
!    OR ( )
!    RSHIFT ( )
!    XOR ( )
!
!  Modified:
!
!    10 February 2008
!
  implicit none

  integer jch

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_WRITE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program makes a file of random integers '
  write ( *, '(a)' ) '  for input to the DIEHARD tests.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Select a sequence from this list:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  1.  A multiply-with-carry (MWC) generator'
  write ( *, '(a)' ) '      x(n) = a * x(n-1) + carry mod 2^32'
  write ( *, '(a)' ) '  2.  A MWC generator on pairs of 16 bits'
  write ( *, '(a)' ) '  3.  The "Mother of all random number generators"'
  write ( *, '(a)' ) '  4.  The KISS generator'
  write ( *, '(a)' ) '  5.  The simple but very good generator COMBO'
  write ( *, '(a)' ) '  6.  The lagged Fibonacci-MWC combination ULTRA'
  write ( *, '(a)' ) '  7.  A combination MWC/subtract-with-borrow '
  write ( *, '(a)' ) '      (SWB) generator, period ~ 10^364'
  write ( *, '(a)' ) '  8.  An extended congruential generator'
  write ( *, '(a)' ) '  9.  The Super-Duper generator'
  write ( *, '(a)' ) ' 10.  A subtract-with-borrow generator'
  write ( *, '(a)' ) ' 11.  Any specified congruential generator'
  write ( *, '(a)' ) ' 12.  The 31-bit generator ran2 from Numerical Recipes'
  write ( *, '(a)' ) ' 13.  Any specified shift-register generator,'
  write ( *, '(a)' ) '      31 or 32 bits'
  write ( *, '(a)' ) ' 14.  The system generator in Sun Fortran f77'
  write ( *, '(a)' ) ' 15.  Any lagged-Fibonacci generator,  '
  write ( *, '(a)' ) '      x(n) = x(n-r) op x(n-s)'
  write ( *, '(a)' ) ' 16.  An inverse congruential generator'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_WRITE:'
  write ( *, '(a)' ) '  Enter your choice, 1 to 16:'

  read ( *, * ) jch

  if ( jch == 1 ) then
    call make_mwc1
  else if ( jch == 2 ) then
    call make_1616
  else if ( jch == 3 ) then
    call make_mthr
  else if ( jch == 4 ) then
    call make_kiss
  else if ( jch == 5 ) then
    call make_cmbo
  else if ( jch == 6 ) then
    call make_ltra
  else if ( jch == 7 ) then
    call make_sbmc
  else if ( jch == 8 ) then
    call make_xcng
  else if ( jch == 9 ) then
    call make_supr
  else if ( jch == 10 ) then
    call make_swb
  else if ( jch == 11 ) then
    call make_cong
  else if ( jch == 12 ) then
    call make_ran2
  else if ( jch == 13 ) then
    call make_shrg
  else if ( jch == 14 ) then
    call make_sunr
  else if ( jch == 15 ) then
    call make_fibo
  else if ( jch == 16 ) then
    call make_invc
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_WRITE - Fatal error!'
    write ( *, '(a,i6)' ) '  Unrecognized choice = ', jch
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_WRITE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function iran2 ( idum )

!*****************************************************************************80
!
!! IRAN2 is a random number generator from Numerical Recipes.
!
  implicit none

  real ( kind = 8 ), parameter :: eps = 1.2D-07
  integer, parameter :: im1 = 2147483563
  integer, parameter :: ntab = 32

  real ( kind = 8 ), parameter :: am = 1.0D+00 / real ( im1, kind = 8 )
  integer, parameter :: ia1 = 40014
  integer, parameter :: ia2 = 40692
  integer idum
  integer, save :: idum2 = 123456789
  integer, parameter :: im2 = 2147483399
  integer, parameter :: imm1 = im1 - 1
  integer, parameter :: iq1 = 53668
  integer, parameter :: iq2 = 52774
  integer, parameter :: ir1 = 12211
  integer, parameter :: ir2 = 3791
  integer iran2
  integer, save, dimension ( ntab ) :: iv = (/ &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0 /)
  integer, save :: iy = 0
  integer j
  integer k
  integer, parameter :: ndiv = 1 + ( imm1 / ntab )
  real ( kind = 8 ), parameter :: rnmx = 1.0D+00 - eps

  if ( idum <= 0 ) then

    idum = max ( -idum, 1 )
    idum2 = idum

    do j = ntab + 8, 1, -1

      k = idum / iq1
      idum = ia1 * ( idum - k * iq1 ) - k * IR1

      if ( idum < 0 ) then
        idum = idum + im1
      end if

      if ( j <= ntab ) then
        iv(j) = idum
      end if

    end do

    iy = iv(1)

  end if

  k = idum / iq1
  idum = ia1 * ( idum - k * iq1 ) - k * ir1

  if ( idum < 0 ) then
    idum = idum + im1
  end if

  k = idum2 / iq2
  idum2 = ia2 * ( idum2 - k * iq2 ) - k * ir2

  if ( idum2 < 0 ) then
    idum2 = idum2 + im2
  end if

  j = 1 + iy / ndiv
  iy = iv(j) - idum2
  iv(j) = idum

  if ( iy < 1 ) then
    iy = iy + imm1
  end if

  iran2 = iy

  return
end
subroutine ks_test ( y, n, p )

!*****************************************************************************80
!
!! KS_TEST applies the Kolmorogov-Smirnov test.
!
!  Discussion:
!
!    The test is based on the distance between the empirical
!    and theoretical distribution functions.
!
!  Modified:
!
!    07 July 2003
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Y(N), an array of numbers, supposed to be drawn
!    from a uniform distribution in [0,1].
!
!    Input, integer N, the number of elements of Y.
!
!    Output, real ( kind = 8 ) P, the probability associated with the
!    observed value of the Anderson-Darling statistic, which is N times
!    the integral of ( FN(X) - X )**2 / ( X * ( 1 - X ) )
!
  implicit none

  integer n

  real ( kind = 8 ) a
  real ( kind = 8 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension ( 8, 10 ) :: l
  integer ( kind = 4 ) m
  real ( kind = 8 ) p
  real ( kind = 8 ) sp
  real ( kind = 8 ) t
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z

  data l /40,46,37,34,27,24,20,20,88,59,43,37,29,27,20, &
  22,92,63,48,41,30,30,25,24,82,59,42,37,26,28,26,22,62,48,33,30 &
  ,23,23,22,18,49,34,22,20,16,17,17,12,17,17,7,8,4,7,5,1,40,18, &
  19,14,16,13,10,9,59,20,10,4,1,1,0,-1,41,43,36,112,15,95,32,58/
!
!  Ascending sort the Y array.
!
  call r8vec_sort_heap_a ( n, y )

  z = - real ( n, kind = 8 ) * real ( n, kind = 8 )

  do i = 1, n

    t = y(i) * ( 1.0D+00 - y(n+1-i) )

    if ( t < 1.0D-20 ) then
      t = 1.0D-20
    end if

    z = z - real ( i + i - 1, kind = 8 ) * log ( t )

  end do

  z = z / real ( n, kind = 8 )

  if ( z < 0.01D+00 ) then

    p = 0.0D+00

  else if ( z <= 2.0D+00 ) then

    p = 2.0D+00 * exp ( -1.2337D+00 / Z ) &
      * ( 1.0D+00 + z / 8.0D+00 &
      - 0.04958D+00 * z * z / ( 1.325D+00 + z ) ) / sqrt ( z )

  else if ( 4.0D+00 < z ) then

    p = 1.0D+00 - 0.4938691D+00 * exp ( -1.050321D+00 * z ) &
      - 0.5946335D+00 * exp ( -1.527198D+00 * z )

  else

    p = 1.0D+00 - 0.6621361D+00 * exp ( -1.091638D+00 * z ) &
      - 0.95059D+00 * exp ( -2.005138D+00 * z )

  end if

  m = min ( n - 2, 8 )

  e = 0.0D+00
  do j = 1, 10
    e = e + l(m,j) * sp ( p, j ) * 0.0001D+00
  end do

  if ( 10 < n ) then
    e = 10.0D+00 * e / real ( n, kind = 8 )
  end if

  a = p + e

  return
end
subroutine make_1616

!*****************************************************************************80
!
!! MAKE_1616 generates data from two 16 bit multiply-with-carry generators.
!
!  Modified:
!
!    10 February 2008
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n(4096)
  character ( len = 80 ) :: output_file = '1616.32'
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MAKE_1616:'
  write ( *, '(a)' ) '  This program creates the binary file 1616.32, '
  write ( *, '(a)' ) '  containing  11+ megabytes of integers made by '
  write ( *, '(a)' ) '  concatenating two 16-bit multiply-with-carry '
  write ( *, '(a)' ) '  generators.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The two generators have the form '
  write ( *, '(a)' ) '    x(n)=a*x(n-1)+carry mod 2^16      and '
  write ( *, '(a)' ) '    y(n)=b*y(n-1)+carry mod 2^16, '
  write ( *, '(a)' ) '  with suggested choices for multipliers A and B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The carry C works as follows:  If a and x are '
  write ( *, '(a)' ) '  16-bit  and C at most 14 bits, then forming a*x+c '
  write ( *, '(a)' ) '  produces an at-most 31-bit result.  That result '
  write ( *, '(a)' ) '  mod 2^16 (the rightmost 16 bits) is the new x '
  write ( *, '(a)' ) '  and the topmost 16 bits the new carry c.  '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sequence of resulting x''s has period the '
  write ( *, '(a)' ) '  order of 2^16 in the group of residues relatively'
  write ( *, '(a)' ) '  prime to m=a*2^16-1, which will be a*2^15-1 for '
  write ( *, '(a)' ) '  the multipliers suggested here.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  You will be prompted to choose a and b and two '
  write ( *, '(a)' ) '  seeds.  Output is a 32-bit integer, the pair x,y '
  write ( *, '(a)' ) '  side by side.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This multiply-with-carry generator is best done '
  write ( *, '(a)' ) '  in assembler, where it takes about 200 nanoseconds '
  write ( *, '(a)' ) '  with a Pentium 120.  A Fortran version takes about '
  write ( *, '(a)' ) '  300 ns.  It seems to pass all tests and is highly '
  write ( *, '(a)' ) '  recommended for speed and simplicity. '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The essence of a version in C requires only two '
  write ( *, '(a)' ) '  statements:'
  write ( *, '(a)' ) '    x = a * ( x & 65535 ) + ( x >> 16 ); '
  write ( *, '(a)' ) '    y = b * ( y & 65535 ) + ( y >> 16 ); '
  write ( *, '(a)' ) '  if x and y are 32-bit integers with carry in '
  write ( *, '(a)' ) '  the top and output in the bottom half.  The '
  write ( *, '(a)' ) '  32-bit integer returned is (x<<16)+(y&65525)'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Select multipliers a and b, a /= b, from this list:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  -----------------------------------------------------------'
  write ( *, '(a)' ) '  18000 18030 18273 18513 18879 19074 19098 19164 19215 19584'
  write ( *, '(a)' ) '  19599 19950 20088 20508 20544 20664 20814 20970 21153 21243'
  write ( *, '(a)' ) '  21423 21723 21954 22125 22188 22293 22860 22938 22965 22974'
  write ( *, '(a)' ) '  23109 23124 23163 23208 23508 23520 23553 23658 23865 24114'
  write ( *, '(a)' ) '  24219 24660 24699 24864 24948 25023 25308 25443 26004 26088'
  write ( *, '(a)' ) '  26154 26550 26679 26838 27183 27258 27753 27795 27810 27834'
  write ( *, '(a)' ) '  27960 28320 28380 28689 28710 28794 28854 28959 28980 29013'
  write ( *, '(a)' ) '  29379 29889 30135 30345 30459 30714 30903 30963 31059 31083'
  write ( *, '(a)' ) '  -----------------------------------------------------------'


  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter a and b: (my favorites are 18000 and 30903)'
  read ( *, * ) a, b

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter two (<=31 bit) seed integers, not zero'
  read ( *, * ) x, y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 2, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  jk = 0

  do i = 1, 700

    jk = jk + 1

    do j = 1, 4096
      x = a * and ( x, 65535 ) + rshift ( x, 16 )
      y = b * and ( y, 65535 ) + rshift ( y, 16 )
      n(j) = lshift ( x, 16 ) + and ( y, 65535 )
    end do

    write(2,rec=jk) n(1:4096)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MAKE_1616:'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  2,867,200 32-bit random integers (11,468,800 bytes)'
  write ( *, '(a)' ) '  have been written to the file.'
  write ( *, '(a)' ) '  They are the concatenated output of two 16-bit'
  write ( *, '(a)' ) '  multiply-with-carry generators.  The period is'

  write ( *, '(f40.0)' ) (a*2.0D+00**15-1.0D+00)*(b*2.0D+00**15-1.0D+00)
  write ( *, '(a)' ) '  (the last 4 digits may be lost to roundoff).'
  k = mod(mod(a*32768-1,10000)*mod(b*32768-1,10000),10000)
  write ( *, '(a,i8)' ) '  The correct last 4 digits are ', k
  return
end
subroutine make_cmbo

!*****************************************************************************80
!
!! MAKE_CMBO generates data from the COMBO generator.
!
!  Simple combo, period> 2^60.5
!  x(n)=x(n-1)*x(n-2) mod 2^32 added to
!  period of x(n)=x(n-1)*x(n-2) is 3*2^29 if seeds odd, and one is +or-3 mod
!  easy to ensure: replace seed x with 3*x*x.
!  mwc z=30903*iand(z,65535)+ishft(z,-16)
!
  implicit none

  integer ( kind = 4 ) b(4096)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  character ( len = 80 ) :: output_file = 'cmbo.32'
  integer ( kind = 4 ) v
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
  integer ( kind = 4 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates a binary file'
  write ( *, '(a)' ) 'containing 11+ megabytes of integers from a simple '
  write ( *, '(a)' ) 'but very good combination generator.  It combines, '
  write ( *, '(a)' ) 'by addition mod 2^32,'
  write ( *, '(a)' ) '  x(n)=x(n-1)*x(n-2) mod 2^32'
  write ( *, '(a)' ) 'and'
  write ( *, '(a)' ) '  y(n)=30903*y(n-1) + carry mod 2^16'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The period of the first is 3*2^29, on odd integers, '
  write ( *, '(a)' ) 'and the period of the second, a multiply-with-carry '
  write ( *, '(a)' ) 'generator, is 30903*2^15-1=1012629503, so the period '
  write ( *, '(a)' ) 'of combo exceeds 2^60. '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This generator is simple to program in Fortran or C '
  write ( *, '(a)' ) 'and quite fast.   It seems to pass all tests in DIEHARD.'
  write ( *, '(a)' ) 'Try it.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'You will be prompted for three seed integers.  The x''s '
  write ( *, '(a)' ) 'of the seeds x1, x2, y must be 3 or 5 mod 8, which is '
  write ( *, '(a)' ) 'ensured by replacing  x1 by 3*(x1+x1+1)^2 and x2 by 2*x2+1.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter three seed integers:'
  read ( *, * ) x, y, z

  x = x + x + 1
  x = 3 * x * x
  y = y + y + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 1, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  jk = 0

  do i = 1, 700

    jk = jk + 1

    do j = 1, 4096
      v = x * y
      x = y
      y = v
      z = 30903 * and ( z, 65535 ) + rshift ( z, 16 )
      b(j) = y + z
    end do

    write(1,rec=jk) b(1:4096)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  2,867,200 32-bit random integers (11,468,800 bytes)'
  write ( *, '(a)' ) '  have been written.'

  return
end
subroutine make_cong

!*****************************************************************************80
!
!! MAKE_CONG generates data from a congruential generator.
!
  implicit none

! implicit real*8(a-h,o-z)

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bb(4096)
  real ( kind = 8 ) da
  real ( kind = 8 ) db
  real ( kind = 8 ) dj
  real ( kind = 8 ) dl
  real ( kind = 8 ) dm
  real ( kind = 8 ) dr
  real ( kind = 8 ) dx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jlat
  integer ( kind = 4 ) jseed
  integer ( kind = 4 ) k
  integer ( kind = 4 ) klim
  integer ( kind = 4 ) kount
  character op
  character ( len = 80 ) :: output_file = 'cong.32'
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  real ( kind = 8  )x(2000)
  real ( kind = 8 ) y(2000)
  real ( kind = 8 ) yy

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates a binary file containing 11+ '
  write ( *, '(a)' ) 'megabytes  of 32-bit random integers from a congruential'
  write ( *, '(a)' ) 'generator, x(n)=a*x(n-1)+b mod m'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'You will be prompted to choose a,b and m, the latter'
  write ( *, '(a)' ) 'in the form m=2^r+s.   If r<=31 then r-bit integers will'
  write ( *, '(a)' ) 'be left-justified (shifted left)  to meet DIEHARD '
  write ( *, '(a)' ) 'requirements.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'If 32 < r then the recursion is carried out in double'
  write ( *, '(a)' ) 'precision and the result sent to the file as the '
  write ( *, '(a)' ) 'integer part of c*x, where c is the ratio 2.^32/m.  '
  write ( *, '(a)' ) 'In addition, this program will display a (line-printer)'
  write ( *, '(a)' ) 'plot of the 2-lattice of the chosen congruential'
  write ( *, '(a)' ) 'generator.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter a,b,r and s for the generator'
  write ( *, '(a)' ) '    x(n)=a*x(n-1)+b mod 2^r+s'
  write ( *, '(a)' ) '  To avoid overflow, make sure a*2^r < 2^53) '
  read ( *, * ) a, b, r, s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter a seed integer:'
  read ( *, * ) jseed

  da = a
  db = b
  dm = 2.0D+00**r + s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If you want the 2-lattice without generating'
  write ( *, '(a)' ) '  the binary file, enter 0, else enter 1:'
  write ( *, '(a)' ) '  Create binary file? 0 for NO, 1 for YES:'
  read ( *, * ) jlat

  if ( jlat == 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Computation is beginning.'
    write ( *, '(a)' ) '  Please be patient!'

    open ( 1, file = output_file, form = 'unformatted', access = 'direct', &
      recl = 16384 )

    dj = jseed
    dr = 2.0D+00**( 32 - r )
    jk = 0

    do i = 1, 700

      jk = jk + 1

      do j = 1, 4096

        dj = mod ( da * dj + db, dm )
        dx = dr * dj

        if ( dx < 2.0D+00**31 ) then
          bb(j) = dx
        else
          bb(j) = dx - 2.0D+00**32
        end if

      end do

      write(1,rec=jk) bb(1:4096)

    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FINISHED'
    write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
    write ( *, '(a)' ) '  2,867,200 32-bit random integers (11,468,800 bytes)'
    write ( *, '(a)' ) '  have been written.'

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  To display the 2-lattice of this generator, '
  write ( *, '(a)' ) '  hit any letter or number key:'
  read ( *, '(a)' ) op

  dl = 2.0D+00**( 20 + r - 32 )
  kount = 0
  klim = da * dl / dm

  do k = 1, klim

    j1 = max ( 0, int ( ( real ( k, kind = 8 ) * dm - b ) / da ) )
    j2 = int ( ( real ( k, kind = 8 ) * dm - b + dl ) / da )

    do j = j1, j2

      yy = mod ( da * j + db, dm )

      if ( yy < dl ) then
        kount = kount + 1
        x(kount) = j / dm
        y(kount) = yy / dm
      end if

    end do

  end do

  call plot1 ( x, y, kount, 54 )

  return
end
subroutine make_fibo

!*****************************************************************************80
!
!! MAKE_FIBO generates data from a lagged Fibonacci generator.
!
  implicit none

  integer ( kind = 4 ) bb(4096)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ijk
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) irand
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jjn
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) js
  integer ( kind = 4 ) ju(607)
  integer ( kind = 4 ) juni
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) op
  character ( len = 80 ) :: output_file = 'fibo.32'
  real ( kind = 8 ) period
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates a binary file of integers from'
  write ( *, '(a)' ) 'a specified lagged Fibonacci generator, '
  write ( *, '(a)' ) '  x(n)=x(n-r) op x(n-s), '
  write ( *, '(a)' ) 'with op one of the four operations +, -, *, xor.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The user is prompted to enter integers R and S '
  write ( *, '(a)' ) 'from a list of lags that provide long periods, '
  write ( *, '(a)' ) 'and a choice of the operation +, -, *, xor.  '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'For simpler program logic, choice of operation is '
  write ( *, '(a)' ) 'indicated by an integer:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  1 for +'
  write ( *, '(a)' ) '  2 for -'
  write ( *, '(a)' ) '  3 for *'
  write ( *, '(a)' ) '  4 for xor'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter lags r and s from this list:'
  write ( *, '(a)' ) '    17,5  33,13  39,14  55,24  63,31  73,25  97,33  607,273'

  read ( *, * ) r, s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter op code: '
  write ( *, '(a)' ) '  1  for +'
  write ( *, '(a)' ) '  2  for -'
  write ( *, '(a)' ) '  3  for *'
  write ( *, '(a)' ) '  4  for xor:'

  read ( *, * ) op

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter four positive integers for seeds:'
  read ( *, * ) i, j, k, l

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 1, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  do ijk = 1, i
    jjn = irand(1)
  end do

  do ii = 1, 97

    js = 0

    do jj = 1, 32

      m = mod ( mod ( i * j, 179 ) * k, 179 )
      i = j
      j = k
      k = m
      l = mod ( 53 * l + 1, 169 )
      js = 2 * js

      if ( 32 <= mod ( l * m, 64 ) ) then
        js = js + 1
      end if

    end do

    if ( op == 3 ) then
      js = xor ( js, 1 )
    end if

    ju(ii) = js

  end do

  ip = r
  jp = s
  jk = 0

  do i = 1, 700

    jk = jk + 1

    do j = 1, 4096

      if ( op == 1 ) then
        juni = ju(ip) + ju(jp)
      else if ( op == 2 ) then
        juni = ju(ip) - ju(jp)
      else if ( op == 3 ) then
        juni = ju(ip) * ju(jp)
      else
        juni = xor ( ju(ip), ju(jp) )
      end if

      ju(ip) = juni

      ip = ip - 1
      if ( ip == 0 ) then
        ip = 97
      end if

      jp = jp - 1
      if ( jp == 0 ) then
        jp = 97
      end if

      bb(j) = juni

    end do

    write ( 1, rec=jk ) bb(1:4096)

  end do

  period = ( 2.0D+00**r - 1.0D+00 ) * 2.0D+00**32

  if ( op == 3 ) then
    period = 0.25D+00 * period
  else if ( op == 4 ) then
    period = 2.0D+00**r - 1.0D+00
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  has been created with 32-bit integers'
  write ( *, '(a)' ) '  from the lagged-Fibonacci sequence'

  if ( op == 1 ) then
    write ( *, 351 ) r,s
351     format(15x,'x(n) = x(n-',i3,') + x(n-',i3,') mod 2^32')
  else if ( op == 2 ) then
    write ( *, 352 )r,s
352     format(15x,'x(n) = x(n-',i3,') - x(n-',i3,') mod 2^32')
  else if ( op == 3 ) then
    write ( *, 353 )r,s
353     format(15x,'x(n) = x(n-',i3,') * x(n-',i3,') mod 2^32')
  else
    write ( *, 354 )r,s
354     format(15x,'x(n) = x(n-',i3,') xor x(n-',i3,')')
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g18.10)' ) '   The period of the sequence is ', period

  return
end
subroutine make_invc

!*****************************************************************************80
!
!! MAKE_INVC generates data from an inverse congruential generator.
!
  implicit none

  integer ( kind = 4 ) :: a = 69069
  real ( kind = 8 ) a0
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  integer ( kind = 4 ) :: b = 362436069
  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) :: dz = 123456789.0D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icong
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) n(4096)
  character ( len = 80 ) :: output_file = 'invc.32'

  b0 = 4294967296.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates a file, INVC.32, with 11+ '
  write ( *, '(a)' ) 'megabytes of 32-bit integers from an inverse '
  write ( *, '(a)' ) 'congruential generator,  described by '
  write ( *, '(a)' ) 'Eichenauer-Herrmann in International Statistics'
  write ( *, '(a)' ) 'Review v60, 167-176 and earlier papers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The attraction of this generator seems to be the '
  write ( *, '(a)' ) 'theory behind it; as a random number generator it '
  write ( *, '(a)' ) 'is not very good!'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'It takes from twenty to fifty times as long as much'
  write ( *, '(a)' ) 'better generators and fails many of the tests in '
  write ( *, '(a)' ) 'DIEHARD.  The first 12-16 bits seem to be good, but '
  write ( *, '(a)' ) 'trailing bits are as bad or worse than those from '
  write ( *, '(a)' ) 'an ordinary congruential RNG mod 2^32.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Try it yourself.  You will be prompted for '
  write ( *, '(a)' ) 'parameters A and B and a seed.  Make sure A mod 4 is 1'
  write ( *, '(a)' ) 'and that B is odd.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter a and b and seed integer, free format:'
  read ( *, * ) a, b, dz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This may take a long time, as this is a very'
  write ( *, '(a)' ) '  very slow generator!'

  open ( 1, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  jk = 0

  do ij = 1, 700

    jk = jk + 1

    do i = 1, 4096

      b1 = dz
      a0 = 0
      a1 = 1
      b0 = 4294967296.0D+00

      do while ( 0 < b1 )
        b2 = mod ( b0, b1 )
        a2 = a0 - a1 * ( ( b0 - b2 ) / b1 )
        b0 = b1
        b1 = b2
        a0 = a1
        a1 = a2
      end do

      dz = mod ( a * b0 * a0 + b, 4294967296.0D+00 )

      if ( dz < 0 ) then
        dz = dz + 4294967296.0D+00
      end if

      if ( 2.0D+00**31 < dz ) then
        icong = dz - 4294967296.0D+00
      else
        icong = dz
      end if

      n(i) = icong

    end do

    write(1,rec=jk) n(1:4096)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED:'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  The file contains 11,428,800 bytes '
  write ( *, '(a)' ) '  of 32 bit integers from the inverse congruential'
  write ( *, '(a)' ) '  generator.'

  return
end
subroutine make_kiss

!*****************************************************************************80
!
!! MAKE_KISS generates data from the KISS generator.
!
  implicit none

  integer ( kind = 4 ) b(4096)
  integer ( kind = 4 ) carry
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) k
  character ( len = 80 ) :: output_file = 'kiss.32'
  integer ( kind = 4 ) m
  integer ( kind = 4 ) w
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
  integer ( kind = 4 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates the binary file kiss.32,'
  write ( *, '(a)' ) 'containing 11+ megabytes of integers from the '
  write ( *, '(a)' ) 'generator KISS, which combines three simple generators.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The acronym KISS means'
  write ( *, '(a)' ) '  Keep It Simple, Stupid!'
  write ( *, '(a)' ) 'and the idea is to use simple, fast, individually '
  write ( *, '(a)' ) 'promising generators to get a composite that will be '
  write ( *, '(a)' ) 'fast, easy to code, have a very long period and pass'
  write ( *, '(a)' ) 'all the tests put to it.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The three components of KISS are '
  write ( *, '(a)' ) '  x(n)=a*x(n-1)+1 mod 2^32'
  write ( *, '(a)' ) '  y(n)=y(n-1)(I+L^13)(I+R^17)(I+L^5),'
  write ( *, '(a)' ) '  z(n)=2*z(n-1)+z(n-2) +carry mod 2^32'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The y''s are a shift register sequence on 32bit binary '
  write ( *, '(a)' ) 'vectors period 2^32-1; see the description in executing'
  write ( *, '(a)' ) 'makesupr.exe.  The z''s are a simple multiply-with-carry '
  write ( *, '(a)' ) 'sequence with period 2^63+2^32-1.  The period of KISS '
  write ( *, '(a)' ) 'is thus 2^32*(2^32-1)*(2^63+2^32-1) > 2^127.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KISS is particularly well suited for assembler '
  write ( *, '(a)' ) 'programming, where it takes about 200 nanoseconds with'
  write ( *, '(a)' ) 'a Pentium 120.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'It seems to pass all tests and is highly recommended'
  write ( *, '(a)' ) 'for speed and simplicity (for generators with that '
  write ( *, '(a)' ) 'long a period)'

  write ( *, '(a)' ) ' '
  write ( *, * ) '  Enter four seed integers, not zero'
  read ( *, * ) x, y, z, w
  carry = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 2, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  jk = 0

  do i = 1, 700

    jk = jk + 1

    do jj = 1, 4096
      x = 69069 * x + 1
      y = xor ( y, lshift ( y, 13 ) )
      y = xor ( y, rshift ( y, 17 ) )
      y = xor ( y, lshift ( y, 5 ) )
      k = rshift ( z, 2 ) + rshift ( w, 3 ) + rshift ( carry, 2 )
      m = w + w + z + carry
      z = w
      w = m
      carry = rshift ( k, 30 )
      b(jj) = x + y + w
    end do

    write(2,rec=jk) b(1:4096)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  2,867,200 32-bit random integers (11,468,800 bytes)'
  write ( *, '(a)' ) '  have been written.'

  return
end
subroutine make_ltra

!*****************************************************************************80
!
!! MAKE_LTRA generates data from the ULTRA generator.
!
  implicit none

  integer ( kind = 4 ) bb(4096)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) js
  integer ( kind = 4 ) ju(99)
  integer ( kind = 4 ) juni
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mwcs
  character ( len = 80 ) :: output_file = 'ltra.32'
  real ( kind = 8 ) px
  integer ( kind = 4 ), parameter :: r = 99
  integer ( kind = 4 ), parameter :: s = 33
  real ( kind = 8 ) xx(607)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates a binary file of integers from '
  write ( *, '(a)' ) 'a version of ULTRA,  a generator we posted on the net '
  write ( *, '(a)' ) 'a few years ago.  It combines the lagged Fibonacci '
  write ( *, '(a)' ) 'generator'
  write ( *, '(a)' ) '  x(n)=x(n-99)*x(n-33)  mod 2^32, x''s odd '
  write ( *, '(a)' ) 'with the multiply-with-carry generator '
  write ( *, '(a)' ) '  y(n)=30903*y(n-1)+carry mod 2^16,'
  write ( *, '(a)' ) 'returning x(n)+y(n) mod 2^32.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'By itself, the lagged Fibonacci generator using '
  write ( *, '(a)' ) 'multiplication passes all tests except those dependent '
  write ( *, '(a)' ) 'on the last bit, since the output integers are always '
  write ( *, '(a)' ) 'odd.  Adding the MWC sequence provides a proper mix of '
  write ( *, '(a)' ) 'trailing bits.  The resulting combination in ULTRA '
  write ( *, '(a)' ) 'seems to pass all tests and has a very long period.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'That period is  3*2^96*(30903*2^15-1),  about 2^127.5'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter four positive integers for seeds:'
  read ( *, * ) i, j, k, l

  mwcs = i + j + k + l

  do ii = 1, r
    i = 18273 * and ( i, 65535 ) + rshift ( i, 16 )
    j = 23163 * and ( j, 65535 ) + rshift ( j, 16 )
    k = 24984 * and ( k, 65535 ) + rshift ( k, 16 )
    l = 28854 * and ( l, 65535 ) + rshift ( l, 16 )
    js = lshift ( i, 16 ) + and ( j, 65535 ) + lshift ( k, 16 ) + and ( l, 65535 )
    js = or ( js, 1 )
    xx(ii) = 0.5 + js * 0.5**32
    ju(ii) = js
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The seed table:'
  write ( *, '(a)' ) ' '
  write ( *, '(5i12)' ) ju(1:r)

  call ks_test ( xx, r, px )

  write ( *, '(a,f8.6)' ) &
    '  P-value for Kolmogorov-Smirnov test on uniformity of seeds: ', px

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 1, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  ip = r
  jp = s
  jk = 0

  do i = 1, 700

    jk = jk + 1

    do j = 1, 4096

      juni = ju(ip) * ju(jp)
      ju(ip) = juni

      ip = ip - 1
      if ( ip == 0 ) then
        ip = 97
      end if

      jp = jp - 1
      if ( jp == 0 ) then
        jp = 97
      end if

      mwcs = 30903 * and ( mwcs, 65535 ) + rshift ( mwcs, 16 )
      bb(j) = juni + mwcs

    end do

    write(1,rec=jk) bb(1:4096)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  2,867,200  32-bit random integers (11,468,800 bytes)'

  return
end
subroutine make_mthr

!*****************************************************************************80
!
!! MAKE_MTHR generates data from the "mother of all generators" generator.
!
  implicit none

  integer ( kind = 4 ), parameter :: a = 2111111111
  integer ( kind = 4 ), parameter :: b = 1492
  integer ( kind = 4 ), parameter :: c = 1776
  integer ( kind = 4 ) carry
  integer ( kind = 4 ), parameter :: d = 5115
  real ( kind = 8 ), parameter :: da = 2111111111.0D+00
  real ( kind = 8 ), parameter :: db = 1492.0D+00
  real ( kind = 8 ), parameter :: dc = 1776.0D+00
  real ( kind = 8 ), parameter :: dd = 5115.0D+00
  real ( kind = 8 ), parameter :: dm = 2.0D+000**32
  real ( kind = 8 ) dv
  real ( kind = 8 ) dw
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) dz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) n(4096)
  character ( len = 80 ) :: output_file = 'mthr.32'
  integer ( kind = 4 ) v
  integer ( kind = 4 ) w
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
  integer ( kind = 4 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MAKE_MTHR:'
  write ( *, '(a)' ) 'This program creates a binary file containing '
  write ( *, '(a)' ) '11+ megabytes of 32-bit integers from the '
  write ( *, '(a)' ) 'multiply-with-carry generator'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  x(n) = 2111111111x(n-4)+1492x(n-3)'
  write ( *, '(a)' ) '        +1776x(n-2)+5115x(n-1)+carry mod 2^32.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The period of this generator is  about 2^160.  '
  write ( *, '(a)' ) 'It is one of what I called "The Mother of All '
  write ( *, '(a)' ) 'Random Number Generators", a few years ago when '
  write ( *, '(a)' ) 'use of Mother of All... was topical and could '
  write ( *, '(a)' ) 'be used for showing bombast, defiance or derision.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'All apply to the usage here.   The carry part, c, '
  write ( *, '(a)' ) 'is the multiple of the modulus b=2^32  dropped in '
  write ( *, '(a)' ) 'the reduction; for example, if the linear '
  write ( *, '(a)' ) 'combination with the current four x''s and carry '
  write ( *, '(a)' ) 'c produced the result  125*b+3621, then the new x'
  write ( *, '(a)' ) 'becomes 3621 and the new carry 125. '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The big advantage of this and other '
  write ( *, '(a)' ) 'multiply-with-carry generators is that they allow'
  write ( *, '(a)' ) 'use of modulus 2^16 & 2^32 without the trailing-bits'
  write ( *, '(a)' ) 'regularities encountered in congruential sequences'
  write ( *, '(a)' ) 'for such moduli.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'But that advantage has to be gained through '
  write ( *, '(a)' ) 'assembly language, if b=2^32, as no common high '
  write ( *, '(a)' ) 'level language seems to allow access to the top '
  write ( *, '(a)' ) '32 bits of the 64-bit product of two 32-bit'
  write ( *, '(a)' ) 'integers. See also the file make1616.exe and '
  write ( *, '(a)' ) 'makemwc1.exe'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter four seed integers:'
  read ( *, * ) x, y, z, w

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 2, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  jk = 0

  do ij = 1, 700

    jk = jk + 1

    do i = 1, 4096

      if ( x < 0 ) then
        dx = x + dm
      else
        dx = x
      end if

      if ( y < 0 ) then
        dy = y + dm
      else
        dy = y
      end if

      if ( z < 0 ) then
        dz = z + dm
      else
        dz = z
      end if

      if ( w < 0 ) then
        dw = w + dm
      else
        dw = w
      end if

      v = a * x + b * y + c * z + d * w + carry
      x = y
      y = z
      z = w
      w = v
      dv = da * dx + db * dy + dc * dz + dd * dw + carry
      carry = int ( dv / dm )

      n(i) = w

    end do

    write ( 2, rec = jk ) n(1:4096)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  2,867,200 32-bit random integers (11,468,800 bytes)'
  write ( *, '(a)' ) '  have been written to the file'
  write ( *, '(a)' ) '  The period is about 2^158.97'

  return
end
subroutine make_mwc1

!*****************************************************************************80
!
!! MAKE_MWC1 generates data from a multiply-with-carry generator.
!
  implicit none

  integer ( kind = 4 ) a(2)
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c(4)
  integer ( kind = 4 ) cc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) n(4096)
  character op
  character ( len = 80 ) :: output_file = 'mwc1.32'
  real ( kind = 8 ) period
  integer ( kind = 4 ) t
  integer ( kind = 4 ) w(4)
  integer ( kind = 4 ) x(2)
  integer ( kind = 4 ) xx
  integer ( kind = 4 ) z(4)

  t ( xx ) = rshift ( xx, 16 )
  b ( xx ) = and ( xx, 65535 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates a binary file containing 11+ '
  write ( *, '(a)' ) 'megabytes of 32-bit integers from a multiply-with-carry'
  write ( *, '(a)' ) 'generator'
  write ( *, '(a)' ) '  x(n)=a*x(n-1)+carry mod 2^32'
  write ( *, '(a)' ) 'You choose the multiplier from a list.   The period '
  write ( *, '(a)' ) 'of the generator will be a*2^31-1.  This class of '
  write ( *, '(a)' ) 'generators is particularly well suited for '
  write ( *, '(a)' ) 'implementation in machine language, and I predict'
  write ( *, '(a)' ) 'that many system generators in the future will be '
  write ( *, '(a)' ) 'of this class rather than the linear congruential '
  write ( *, '(a)' ) 'generators for modulus 2^32 that are common today.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'To illustrate how the carry works, suppose from'
  write ( *, '(a)' ) 'the current (32-bit) x and (30 bit) c, one forms '
  write ( *, '(a)' ) 'a*x+c.  This may be done in a 64-(or double 32-) '
  write ( *, '(a)' ) 'bit register in most modern CPU''s.  Then the new '
  write ( *, '(a)' ) 'random x is the lower 32 bits, the new carry the '
  write ( *, '(a)' ) 'upper 32.  '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'To see how well such a simple and fast generator '
  write ( *, '(a)' ) 'performs on tests of randomness, this program makes'
  write ( *, '(a)' ) 'a large file with the multiply-with-carry generator'
  write ( *, '(a)' ) 'implemented in 16-bit integer arithmetic.  Those '
  write ( *, '(a)' ) 'finding it suitable may wish to make an assembler '
  write ( *, '(a)' ) 'version for their system.'
  write ( *, '(a)' ) 'It seems to pass all tests.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Select multiplier "a" from this list:'
  write ( *, '(a)' ) '-------------------------------------------------------'
  write ( *, '(a)' ) ' 1791398085 1929682203 1683268614 1965537969 1675393560'
  write ( *, '(a)' ) ' 1967773755 1517746329 1447497129 1655692410 1606218150'
  write ( *, '(a)' ) ' 2051013963 1075433238 1557985959 1781943330 1893513180'
  write ( *, '(a)' ) ' 1631296680 2131995753 2083801278 1873196400 1554115554'
  write ( *, '(a)' ) '-------------------------------------------------------'
  write ( *, '(a)' ) '  (In general, for any choice of `a`, let m=a*2^32-1. '
  write ( *, '(a)' ) '  If both m and (m-1)/2 are prime, then the period will'
  write ( *, '(a)' ) '  be (m-1)/2).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter multiplier A:'
  read ( *, * ) aa
  write ( *, '(a)' ) '  Enter a seed integer x and initial carry c:'
  read ( *, * ) xx, cc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 2, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  x(1) = b(xx)
  x(2) = t(xx)
  c(1) = b(cc)
  c(2) = t(cc)
  c(3) = 0
  c(4) = 0
  a(1) = b(aa)
  a(2) = t(aa)

  jk = 0

  do i = 1, 700

    jk = jk + 1

    do j = 1, 4096
      call prod ( x, a, z )
      call summer ( z, c, w )
      x(1) = w(1)
      x(2) = w(2)
      c(1) = w(3)
      c(2) = w(4)
      n(j)= lshift ( x(2), 16 ) + x(1)
    end do

    write ( 2, rec = jk ) n(1:4096)

  end do

  period = real ( aa, kind = 8 ) * 2.0d+00**31 - 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  2,867,200 32-bit random integers (11,468,800 bytes)'
  write ( *, '(a)' ) '  have been written.'
  write ( *, '(a,g14.6)' ) '  The period is ', period

  return
end
subroutine make_ran2

!*****************************************************************************80
!
!! MAKE_RAN2 generates data from a Numerical Recipes generator.
!
  implicit none

  integer ( kind = 4 ) b(4096)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idum
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) iran2
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jseed
  character ( len = 80 ) :: output_file = 'ran2.32'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates a binary file containing'
  write ( *, '(a)' ) 'some 11 megabytes of integers from the generator RAN2'
  write ( *, '(a)' ) 'in Numerical Recipes.  See Press and Teukolsky, '
  write ( *, '(a)' ) '"Portable Random Number Generators", Computers in '
  write ( *, '(a)' ) 'Physics, v 6, n 2, 522:524 1992.'
  write ( *, '(a)' ) 'See also "Some portable very-long-period random number'
  write ( *, '(a)' ) 'generators", v 8 n 1, 1994 by Marsaglia and Zaman, '
  write ( *, '(a)' ) 'that points out shortcomings, suggested improvements'
  write ( *, '(a)' ) 'and alternatives to RAN2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The Numerical Recipes version produces only 31 bit'
  write ( *, '(a)' ) 'integers.  They are left-adjusted, that is, shifted '
  write ( *, '(a)' ) 'left 1, before being written to the binary file '
  write ( *, '(a)' ) 'as several DIEHARD tests emphasize the leading bits.'
  write ( *, '(a)' ) 'But that means that all integers have trailing bit=0,'
  write ( *, '(a)' ) 'so that tests that use bit 32 (DIEHARD numbers bits'
  write ( *, '(a)' ) 'from left to right, 1 to 32) will yield spectacular'
  write ( *, '(a)' ) 'failures.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter a seed integer'
  read ( *, * ) jseed

  idum = -abs ( jseed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 2, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  jk = 0

  do ij = 1, 700

    jk = jk + 1

    do i = 1, 4096
      b(i) = lshift ( iran2 ( idum ), 1 )
    end do

    write ( 2, rec = jk ) b(1:4096)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data was written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  Now 2,867,200 left-adjusted integers from ran2'
  write ( *, '(a)' ) '  (11,486,800 bytes) have been written.'

  return
end
subroutine make_sbmc

!*****************************************************************************80
!
!! MAKE_SBMC generates data from a subtract/borrow/multiply/carry generator.
!
  implicit none

  integer ( kind = 4 ) bb(4096)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) js
  integer ( kind = 4 ) ju(37)
  integer ( kind = 4 ) juni
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mwcs
  character ( len = 80 ) :: output_file = 'sbmc.32'
  real ( kind = 8 ) px
  integer ( kind = 4 ), parameter :: r = 37
  integer ( kind = 4 ), parameter :: s = 24
  integer ( kind = 4 ) x
  real ( kind = 8 ) xx(37)
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates a binary file of integers '
  write ( *, '(a)' ) 'from the subtract-with-borrow random number '
  write ( *, '(a)' ) 'generator:'
  write ( *, '(a)' ) '  x(n)=x(n-24)-x(n-37) - borrow mod 2^32 '
  write ( *, '(a)' ) 'combined with the multiply-with-carry generator'
  write ( *, '(a)' ) '  y(n)=30903*y(n-1) + carry  mod 2^16.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The period of the composite is '
  write ( *, '(a)' ) '(2^1178-2^762)(30903*2^15-1),  about 2^1208 or 10^364.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter four positive integers for seeds:'
  read ( *, * ) i, j, k, l

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 1, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  mwcs = i + j + k + l

  do ii = 1, r

    js = 0

    do jj = 1, 32

      m = mod ( mod ( i * j, 179 ) * k, 179 )
      i = j
      j = k
      k = m
      l = mod ( 53 * l + 1, 169 )
      js = 2 * js

      if ( 32 <= mod ( l * m, 64 ) ) then
        js = js + 1
      end if

    end do

    xx(ii) = 0.5D+00 + 0.5D+00**32 * js
    ju(ii) = js

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The seed table:'
  write ( *, '(a)' ) ' '

  write ( *, '(5i12)' ) ju(1:r)

  call ks_test ( xx, r, px )

  write ( *, '(a,f8.6)' ) &
    '  P-value for Kolmogorov-Smirnov test on uniformity of seeds: ', px

  ip = r
  jp = s
  jk = 0

  do i = 1, 700

    jk = jk + 1

    do j = 1, 4096

      x = ju(ip)
      y = ju(jp)
      juni = x - y

      if ( rshift ( y, 1 ) < rshift ( x, 1 ) ) then
        juni = juni - 1
      end if

      ju(ip) = juni

      ip = ip - 1
      if ( ip == 0 ) then
        ip = r
      end if

      jp = jp - 1
      if ( jp == 0 ) then
        jp = r
      end if

      mwcs = 30903 * and ( mwcs, 65535 ) + rshift ( mwcs, 16 )
      bb(j) = mwcs + juni

    end do

    write(1,rec=jk) bb(1:4096)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  from the subtract-with-borrow sequence'
  write ( *, '(a)' ) '  x(n)=x(n-24)-x(n-37)-borrow mod 2^32'
  write ( *, '(a)' ) '  combined with the multiply-with-carry sequence'
  write ( *, '(a)' ) '  y(n)=30903*y(n-1)+carry mod 2^16,'
  write ( *, '(a)' ) '  the overall sequence having period about 10^364.'

  return
end
subroutine make_shrg

!*****************************************************************************80
!
!! MAKE_SHRG generates data from a specified shift register generator.
!
!  creates shift register random number file
!  use 13,18 or 7,24 or 6,25 or 3,28 or reverse for 31 bits
!  use 15,17 or 13,17,5 for 32 bits
!
  implicit none

  integer ( kind = 4 ) b(4096)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mask
  integer ( kind = 4 ) mr
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbits
  character ( len = 80 ) :: output_file = 'shrg.32'
  integer ( kind = 4 ) r

  m(k,n) = xor(k,lshift(k,n))
  mr(k,n) = xor(k,rshift(k,n))

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates a binary file with 11+ megabytes'
  write ( *, '(a)' ) 'of 32-bit integers from a specified shift register '
  write ( *, '(a)' ) 'generator.  You will be prompted to choose the shifts'
  write ( *, '(a)' ) 'and whether to generate 31- or 32-bit integers.  '
  write ( *, '(a)' ) 'If 31 bits are chosen, the resulting integers will '
  write ( *, '(a)' ) 'be left-justified (shifted left 1) before being '
  write ( *, '(a)' ) 'written to the output file you name, since many of '
  write ( *, '(a)' ) 'DIEHARD''s tests emphasize leading bits of random '
  write ( *, '(a)' ) 'integers.  That means, of course, that the generator'
  write ( *, '(a)' ) 'is likely to fail tests that depend, in any significant'
  write ( *, '(a)' ) 'way, on the rightmost bit of a 32-bit random integer.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter number of bits, 31 or 32:'
  read ( *, * ) nbits

  jk = 0

  if ( nbits == 31 ) then

    write ( *, '(a)' ) '  Choose Left,Right shifts from these choices:'
    write ( *, '(a)' ) '    13,18  18,13  24,7  7,24 6,25  25,6  28,3  3,28'
    write ( *, '(a)' ) '  Enter two integers, L and R:'
    read ( *, * ) l, r

    mask = huge ( mask )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter seed integer, not zero:'
    read ( *, * ) j

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Computation is beginning.'
    write ( *, '(a)' ) '  Please be patient!'

    open ( 1, file = output_file, form = 'unformatted', access = 'direct', &
      recl = 16384 )

    do i1 = 1, 700

      jk = jk + 1

      do i = 1, 4096
        j = mr ( and ( m(j,l), mask ), r )
        b(i) = j + j
      end do

      write(1,rec=jk) b(1:4096)

    end do

  end if

  if ( nbits == 32 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  For 32 bit integers, you have three choices.'
    write ( *, '(a)' ) '  two shifts: (L,R)=(17,15) or (15,17) and'
    write ( *, '(a)' ) '  three shifts (L1,R,L2)=(13,17,5).'
    write ( *, '(a)' ) '  Enter three integers, 17 15 0 or 15 17 0 or 13 17 5:'
    read ( *, * ) l, r, l2

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter seed integer, not zero:'
    read ( *, * ) j

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Computation is beginning.'
    write ( *, '(a)' ) '  Please be patient!'

    open ( 1, file = output_file, form = 'unformatted', access = 'direct', &
      recl = 16384 )

    if ( l2 == 0 ) then

      do i1 = 1, 700

        jk = jk + 1

        do i = 1, 4096
          j = mr(m(j,l),r)
          b(i) = j
        end do

        write(1,rec=jk) b(1:4096)

      end do

    else if ( l2 == 5 ) then

      do i1 = 1, 700

        jk = jk + 1

        do i = 1, 4096
          j = m ( mr ( m(j,L), R ), L2 )
          b(i) = j
        end do

        write(1,rec=jk) b(1:4096)

      end do

    end if

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  2,867,200 32-bit random integers (11,468,800 bytes)'
  write ( *, '(a)' ) '  have been written to the file.'

  if ( nbits == 31 .or. l2 == 5 ) then
    write ( *, '(a,i2,a)' ) ' The period of your generator is 2^', nbits ,'-1.'
  end if

  if ( l == 15 .or. l == 17 ) then
    write ( *, '(a)' ) ' The period of your generator is 2^32-2^21-2^11+1.'
  end if

  return
end
subroutine make_sunr

!*****************************************************************************80
!
!! MAKE_SUNR generates data from the Sun f77 generator.
!
  implicit none

  integer ( kind = 4 ) b(4096)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ijk
  integer ( kind = 4 ) irand
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jseed
  character ( len = 80 ) :: output_file = 'sunran.32'

  write ( *, '(a)' ) 'The system generator in Sun Fortran  f77'
  write ( *, '(a)' ) 'This program creates the binary file sunr.32,'
  write ( *, '(a)' ) 'containing 11+ megabytes of 31-bit integers made '
  write ( *, '(a)' ) 'from the generator in Sun Fortran f77.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'It is not a congruential generator.'
  write ( *, '(a)' ) 'I do not have a manual that tells what it is.  '
  write ( *, '(a)' ) 'But whatever it is, it is not a very good generator.'
  write ( *, '(a)' ) 'Try it for yourself.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'A cryptic f77 manual describes its use. It has a '
  write ( *, '(a)' ) 'desirable feature of providing either reals on [0,1)'
  write ( *, '(a)' ) 'or 32-bit integers (32-bit would be preferable). '
  write ( *, '(a)' ) 'The Fortran assignment x=rand(0) will provide the '
  write ( *, '(a)' ) 'next real in the sequence, which can be reseeded '
  write ( *, '(a)' ) 'with x=rand(iseed).  Similarly, 31-bit integers come'
  write ( *, '(a)' ) 'from successive assignments j=irand(0), with a new '
  write ( *, '(a)' ) 'seed from j=irand(iseed).  That calling procedure '
  write ( *, '(a)' ) 'seems better than that of Microsoft Fortran, which uses '
  write ( *, '(a)' ) 'call random(x) to get a random real x;  a nuisance '
  write ( *, '(a)' ) 'for those wanting to use a random variable in an '
  write ( *, '(a)' ) 'expression.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter one seed integer:'
  read ( *, * ) jseed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 2, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  ijk = irand ( jseed )
  ijk = 1
  jk = 0

  do j = 1, 700

    jk = jk + 1

    do i = 1, 4096
      b(i) = irand ( 1 )
    end do

    write(2,rec=jk) b(1:4096)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  2,867,200  32-bit random integers (11,468,800 bytes)'
  write ( *, '(a)' ) '  have been written to the file.'
  write ( *, '(a)' ) '  Note: the rightmost bit is always 0'

  return
end
subroutine make_supr

!*****************************************************************************80
!
!! MAKE_SUPR generates data from the superduper generator.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) n(4096)
  character op
  character ( len = 80 ) :: output_file = 'supr.32'
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The random number generator SUPERDUPER'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates the binary file SUPR.32,'
  write ( *, '(a)' ) 'containing 11+ megabytes of integers made by adding '
  write ( *, '(a)' ) 'the results of two random number generators, congruential '
  write ( *, '(a)' ) 'and shift register.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The two generators are'
  write ( *, '(a)' ) '  x <--- 69069*x+oddconstant mod 2^32'
  write ( *, '(a)' ) '  y <--- y(I+L^13)(I+R^17)(I+L^5),'
  write ( *, '(a)' ) 'where y is viewed as a binary vector. The transformation is'
  write ( *, '(a)' ) '  y <-- yT, with T the 32x32 binary matrix'
  write ( *, '(a)' ) 'T=(I+L^13)(I+R^17)(I+L^5) and L, R are matrices that '
  write ( *, '(a)' ) 'effect a shift of 1 left or 1 right.  The transformation '
  write ( *, '(a)' ) 'is readily done with shifts and exclusive or''s.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The period of the x''s is 2^32 for any initial seed and '
  write ( *, '(a)' ) 'that for the y''s is 2^32-1 for any seed not zero.  '
  write ( *, '(a)' ) 'So the period of  SUPERDUPER is '
  write ( *, '(a)' ) '2^64-2^32 = 18,446,744,069,414,584,320.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'You will be prompted for two seeds.'
  write ( *, '(a)' ) 'You will also be prompted to choose the method for combining'
  write ( *, '(a)' ) 'the two sequences: addition (enter +), exclusive-or (enter x) '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter two seed integers, the second not zero:'
  read ( *, * ) x, y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Emter your method of combination:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  +  for addition, '
  write ( *, '(a)' ) '  x  for exclusive-or, in column 1:'
  read ( *, '(a)' ) op

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 2, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  jk = 0

  if ( op == '+' ) then

    do i = 1, 700

      jk = jk + 1

      do j = 1, 4096
        x = 69069 * x + 1
        y = xor ( y, lshift ( y, 13 ) )
        y = xor ( y, rshift ( y, 17 ) )
        y = xor ( y, lshift ( y, 5 ) )
        n(j) = x + y
      end do

      write(2,rec=jk) n(1:4096)

    end do

  end if

  if ( op == 'x' ) then

    do i = 1, 700

      jk = jk + 1

      do j = 1, 4096
        x = 69069 * x + 1
        y = xor ( y, lshift ( y, 13 ) )
        y = xor ( y, rshift ( y, 17 ) )
        y = xor ( y, lshift ( y, 5 ) )
        n(j) = xor ( x, y )
      end do

      write(2,rec=jk) n(1:4096)

    end do

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  2,867,200 32-bit random integers (11,468,800 bytes)'
  write ( *, '(a)' ) '  have been written.'

  return
end
subroutine make_swb

!*****************************************************************************80
!
!! MAKE_SWB generates data from a subtract-with-borrow generator.
!
  implicit none

  integer ( kind = 4 ) bb(4096)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ijk
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) irand
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jchoice
  integer ( kind = 4 ), dimension ( 5 ) :: jix = (/ -6, 0, 0, 0, 2147483647 /)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jjn
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) js
  integer ( kind = 4 ) ju(48)
  integer ( kind = 4 ) juni
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  character ( len = 80 ) :: output_file = 'swb.32'
  integer ( kind = 4 ) r
  integer ( kind = 4 ), dimension ( 5 ) :: rc = (/ 43, 37, 24, 21, 48 /)
  integer ( kind = 4 ) s
  integer ( kind = 4 ), dimension ( 5 ) :: sc = (/ 22, 24, 19, 6, 8 /)
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y

  character*42 period(5)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates a binary file of integers from '
  write ( *, '(a)' ) 'a specified subtract-with-borrow random number generator.'
  write ( *, '(a)' ) 'Structures of such generators are much like those of '
  write ( *, '(a)' ) 'lagged Fibonacci generators x(n)=x(n-r)-x(n-s) mod m, '
  write ( *, '(a)' ) 'except that one forms x(n)=x(n-r)-x(n-s)-c mod m, '
  write ( *, '(a)' ) 'where c is the borrow: 0 if the subtraction produced '
  write ( *, '(a)' ) 'a positive result, and 1 if an m had to be borrowed '
  write ( *, '(a)' ) 'and added to the difference to make a positive result.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'With proper choice of the lags, tremendously long periods'
  write ( *, '(a)' ) 'can be attained.  But performance on tests is much the '
  write ( *, '(a)' ) 'same as for lagged Fibonacci, as SWB sequences behave '
  write ( *, '(a)' ) 'locally much like lagged Fibonacci using subtraction.  '
  write ( *, '(a)' ) 'There are also add-with-carry generators '
  write ( *, '(a)' ) '  x(n)=x(n-r)+x(n-s)+carry mod m'
  write ( *, '(a)' ) 'that use addition, with the carry set to 1 or 0, '
  write ( *, '(a)' ) 'depending on whether m had to be subtracted to give a '
  write ( *, '(a)' ) 'positive result mod m.  For a full description, see'
  write ( *, '(a)' ) 'Marsaglia and Zaman, '
  write ( *, '(a)' ) 'Annals of Applied Probability, '
  write ( *, '(a)' ) 'Volume 1, Number 3, 1991.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Choose your generator from this list:'
  write ( *, '(a)' ) '   1:   x(n)=x(n-43)-x(n-22)-c mod 2^32-5'
  write ( *, '(a)' ) '   2:   x(n)=x(n-37)-x(n-24)-c mod 2^32  '
  write ( *, '(a)' ) '   3:   x(n)=x(n-24)-x(n-19)-c mod 2^32  '
  write ( *, '(a)' ) '   4:   x(n)=x(n-21)-x(n- 6)-c mod 2^32  '
  write ( *, '(a)' ) '   5:   x(n)=x(n-48)-x(n- 8)-c mod 2^31  '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  (Choice 5 provides 31-bit integers that are left-justified)'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter 1,2,3,4 or 5:'
  read ( *, * ) jchoice

  r = rc(jchoice)
  s = sc(jchoice)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter four positive integers for seeds:'
  read ( *, * ) i, j, k, l

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 1, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  do ijk = 1, i
    jjn = irand ( 1 )
  end do

  do ii = 1, r

    js = 0

    do jj = 1, 32

      m = mod ( mod ( I * J, 179 ) * K, 179 )
      i = j
      j = k
      k = m
      l = mod ( 53 * l + 1, 169 )
      js = 2 * js

      if ( 32 <= mod ( l * m, 64 ) ) then
        js = js + 1
      end if

    end do

    if ( jchoice == 5 ) then
      js = xor ( js, 1 )
    end if

    ju(ii) = js

  end do

  ip = r
  jp = s
  jk = 0

  do i = 1, 700

    jk = jk + 1

    do j = 1, 4096

      x = ju(ip)
      y = ju(jp)

      juni = x - y
      if ( rshift ( x, 1 ) < rshift ( y, 1 ) ) then
        juni = juni + jix(jchoice)
      end if

      ju(ip) = juni

      ip = ip - 1
      if ( ip == 0 ) then
        ip = r
      end if

      jp = jp - 1
      if ( jp == 0 ) then
        jp = r
      end if

      if ( jchoice == 5 ) then
        juni = juni + juni
      end if

      bb(j) = juni

    end do

    write(1,rec=jk) bb(1:4096)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  from the subtract-with-borrow sequence'

  write ( *, 350 ) period(jchoice)
350   format(15x, 'The period is',/,a41)
  data period/ &
  '  (2^32-5)^43 - (2^32-5)^22, about 10^414', &
  '  2^1178 - 2^762, about 10^354           ', &
  '  (2^759 - 2^599)/3, about 10^228        ', &
  '  (2^666 - 2^186)/3, about 10^200        ', &
  '  (2^1478 - 2^247)/105, about 10^445     '/

  return
end
subroutine make_xcng

!*****************************************************************************80
!
!! MAKE_XCNG generates data from an extended congruential generator.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jch
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) n(4096)
  character ( len = 80 ) :: output_file = 'xcng.32'
  real ( kind = 8 ), parameter :: p = 2.0D+00**32-5.0D+00
  real ( kind = 8 ) r3
  real ( kind = 8 ) r4
  real ( kind = 8 ), parameter :: s = 2.0D+00**31
  real ( kind = 8 ) x
  character ( len = 80 ) text(17),dum

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program creates the binary file XCNG.32 '
  write ( *, '(a)' ) 'containing 11+ megabytes of 32-bit integers from'
  write ( *, '(a)' ) 'an extended congruential generator.  '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'You have your choice of four:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  1:  x(n)=65065*x(n-1)+67067*x(n-2)+69069*x(n-3) mod 2^32-5'
  write ( *, '(a)' ) '  2:  x(n)=2**10*[x(n-1)+x(n-2)+x(n-3)] mod 2^32-5'
  write ( *, '(a)' ) '  3:  x(n)=2000*x(n-1)+1950*x(n-2)+1900*x(n-3) mod 2^32-209'
  write ( *, '(a)' ) '  4:  x(n)=2**20*[x(n-1)+x(n-2)+x(n-3)] mod 2^32-209 '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The period of these generators is the modulus '
  write ( *, '(a)' ) 'cubed minus 1, about 2^96, for any 3 seeds not all '
  write ( *, '(a)' ) 'zero.   The recursion is implemented in double '
  write ( *, '(a)' ) 'precision, with the result converted to a 32-bit '
  write ( *, '(a)' ) 'integer.  Notice that choices 2 and 4 are well suited '
  write ( *, '(a)' ) 'for implementations that avoid multiplication: adding '
  write ( *, '(a)' ) '10 or 20 to the exponent of a double precision sum.   '
  write ( *, '(a)' ) 'Since  all four choices seem to pass all tests, 2) and '
  write ( *, '(a)' ) '4) may be preferable.  Try them yourself.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the generator, 1, 2, 3 or 4:'
  read ( *, * ) jch

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter three seed integers, not all zero:'
  read ( *, * ) a, b, c

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computation is beginning.'
  write ( *, '(a)' ) '  Please be patient!'

  open ( 2, file = output_file, form = 'unformatted', access = 'direct', &
    recl = 16384 )

  r3 = 2.0D+00**32 / 34359737519.0D+00
  r4 = 2.0D+00**32 / 34359736739.0D+00

  jk = 0

  do ij = 1, 700

    jk = jk + 1

    do i = 1, 4096

      if ( jch == 1 ) then
        x = mod ( 1776.0D+00 * a + 1476.0D+00 * b + 1176.0D+00 * c, &
          4294967291.0D+00 )
      else if ( jch == 2 ) then
        x = mod ( 8192 * ( a + b + c ), 4294967291.0D+00 )
      else if ( jch == 3 ) then
        x = mod ( 2001.0D+00 * a + 1998.0D+00 * b + 1995.0D+00 * c, &
          34359737519.0D+00 )
        x = x * r3
      else if ( jch == 4 ) then
        x = mod ( 524288.0D+00 * ( a + b + c ), 34359736739.0D+00 )
        x = x * r4
      end if

      a = b
      b = c
      c = x

      if ( x < s ) then
        n(i) = x
      else
        n(i) = x - 2.0D+00**32
      end if

    end do

    write(2,rec=jk) n(1:4096)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FINISHED'
  write ( *, '(a)' ) '  The data has been written to "' // trim ( output_file ) // '".'
  write ( *, '(a)' ) '  2,867,200 32-bit random integers (11,468,800 bytes)'
  write ( *, '(a)' ) '  have been written.'

  return
end
subroutine plot1 ( x, y, n, nc )

!*****************************************************************************80
!
!! PLOT1 makes a printer plot of data.
!
  implicit none

  integer n

  character, parameter :: blank = ' '
  character c(90,54)
  character f(4860)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) nr
  logical p2
  character, parameter :: plus = '+'
  character, parameter :: star = '*'
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xd
  real ( kind = 8 ) xl
  real ( kind = 8 ) xr
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yb
  real ( kind = 8 ) yd
  real ( kind = 8 ) yt
  real ( kind = 8 ) z(n)
  character, parameter :: zero = '0'

  character*4860 D,E

  equIVALENCE (C,D),(E,F)
  DATA F/4860*' '/

  P2=.FALSE.

  nr = 48 * nc / 100

  if ( 60 < nc ) then
    nr = 6 * nc / 10
  end if

  xl = x(1)
  xr = xl
  yb = y(1)
  yt = yb

  do i = 1, n

    xl = min ( xl, X(I) )
    XR = max ( XR, X(I) )

    if ( P2 ) then
      YB = min ( YB, Z(I) )
      YT = max ( YT, Z(I) )
    end if

    yb = min ( yb, Y(I) )
    YT = max ( YT, Y(I) )

  end do

  PRINT 22,XL,XR,YB,YT
 22   FORMAT('   X RANGE:',G12.3,' TO ',G12.3,'  Y,Z RANGE:',G12.3,'  TO',G12.3)

  xd = 1.0D+00
  if ( xl /= xr ) then
    xd = nc / ( xr - xl )
  end if

  yd = 1.0D+00
  if ( yb /= yt ) then
    yd = nr / ( yt - yb )
  end if

  d = e

  do m = 1, n

    i = int ( 1.0D+00 + ( X(M) - XL ) * XD )
    j = int ( 1.0D+00 + ( Y(M) - YB ) * YD )
    c(i,j) = '*'

    if ( P2 ) then

      k = int ( 1.0D+00 + ( Z(M) - YB ) * YD )
      c(i,k) = '+'

      if ( j == k ) then
        c(i,k) = '0'
      end if

    end if

  end do

  do j = nr, 1, -1
    write ( *, '(10x,90a1)' ) c(1:nc,j)
  end do

  return
end
subroutine prod ( x, y, z )

!*****************************************************************************80
!
!! PROD ???
!
  implicit none

  integer ( kind = 4 ) b
  integer ( kind = 4 ) d
  integer ( kind = 4 ) r(4)
  integer ( kind = 4 ) s(4)
  integer ( kind = 4 ) t
  integer ( kind = 4 ) x(2)
  integer ( kind = 4 ) y(2)
  integer ( kind = 4 ) z(4)

  t(d) = rshift ( d, 16 )
  b(d) = and ( d, 65535 )

  d = x(1) * y(1)
  z(1) = b(d)
  d = t(d) + x(1) * y(2)
  r(2) = b(d)
  r(3) = t(d)
  d = x(2) * y(1)
  s(2) = b(d)
  d = t(d) + x(2) * y(2)
  s(3) = b(d)
  s(4) = t(d)
  d = r(2) + s(2)
  z(2) = b(d)
  d = t(d) + r(3) + s(3)
  z(3) = b(d)
  z(4) = t(d) + s(4)

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r8vec_heap_d ( n, a )

!*****************************************************************************80
!
!! R8VEC_HEAP_D reorders an R8VEC into an descending heap.
!
!  Definition:
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Modified:
!
!    07 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the size of the input array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  real ( kind = 8 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine r8vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_A ascending sorts an R8VEC using heap sort.
!
!  Modified:
!
!    07 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call r8vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  call r8_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n-1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call r8vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    call r8_swap ( a(1), a(n1) )

  end do

  return
end
function sp ( x, i )

!*****************************************************************************80
!
!! SP is called by KS_TEST.
!
!  Modified:
!
!    13 November 2011
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) sp
  real ( kind = 8 ) t
  real ( kind = 8 ) x

  if ( i == 8 ) then

    if ( x <= 0.8D+00 ) then
      sp = 0.0D+00
    else if ( x < 1.0D+00 ) then
      sp = 100.0D+00 * ( x - 0.9D+00 )**2 - 1.0D+00
    else
      sp = 0.0D+00
    end if

  else if ( i == 9 ) then

    if ( x <= 0.0D+00 ) then
      sp = 0.0D+00
    else if ( x <= 0.01D+00 ) then
      sp = -100.0D+00 * x
    else if ( x < 0.05D+00 ) then
      sp = 25.0D+00 * ( x - 0.05D+00 )
    else
      sp = 0.0D+00
    end if

  else if ( i == 10 ) then

    if ( x <= 0.98D+00 ) then
      sp = 0.0D+00
    else if ( x < 1.0D+00 ) then
      sp = 0.1D+00 - 10.0D+00 * abs ( x - 0.99D+00 )
    else
      sp = 0.0D+00
    end if

  else

    t = abs ( 10.0D+00 * x - 0.5D+00 - real ( i, kind = 8 ) )

    if ( t <= 0.5D+00 ) then
      sp = 1.5D+00 - 2.0D+00 * t * t
    else if ( t <= 1.5D+00 ) then
      sp = 2.25D+00 - t * ( 3.0D+00 - t )
    else
      sp = 0.0D+00
    end if

  end if

  return
end
subroutine summer ( x, y, z )

!*****************************************************************************80
!
!! SUMMER computes X(1:4) + Y(1:4).
!
  implicit none

  integer ( kind = 4 ) b
  integer ( kind = 4 ) d
  integer ( kind = 4 ) t
  integer ( kind = 4 ) x(4)
  integer ( kind = 4 ) y(4)
  integer ( kind = 4 ) z(4)

  t(d) = rshift ( d, 16 )
  b(d) = and ( d, 65535 )

  d = x(1) + y(1)
  z(1) = b(d)
  d = t(d) + x(2) + y(2)
  z(2) = b(d)
  d = t(d) + x(3) + y(3)
  z(3) = b(d)
  d = t(d) + x(4) + y(4)
  z(4) = b(d)

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
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
