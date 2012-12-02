subroutine ccu ( n, x, w )

!*****************************************************************************80
!
!! CCU computes a Clenshaw Curtis quadrature rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [0,1].
!
!    The integral to approximate:
!
!      Integral ( 0 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCU - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop
  end if

  if ( n == 1 ) then

    x(1) = 0.0D+00
    w(1) = 2.0D+00

  else

    do i = 1, n
      x(i) = cos ( real ( n - i, kind = 8 ) * pi &
                 / real ( n - 1, kind = 8 ) )
    end do

    x(1) = -1.0D+00
    if ( mod ( n, 2 ) == 1 ) then
      x((n+1)/2) = 0.0D+00
    end if
    x(n) = +1.0D+00

    do i = 1, n

      theta = real ( i - 1, kind = 8 ) * pi &
            / real ( n - 1, kind = 8 )

      w(i) = 1.0D+00

      do j = 1, ( n - 1 ) / 2

        if ( 2 * j == ( n - 1 ) ) then
          b = 1.0D+00
        else
          b = 2.0D+00
        end if

        w(i) = w(i) - b * cos ( 2.0D+00 * real ( j, kind = 8 ) * theta ) &
             / real ( 4 * j * j - 1, kind = 8 )

      end do

    end do

    w(1)     =           w(1)     / real ( n - 1, kind = 8 )
    w(2:n-1) = 2.0D+00 * w(2:n-1) / real ( n - 1, kind = 8 )
    w(n)     =           w(n)     / real ( n - 1, kind = 8 )

  end if
!
!  Transform from [-1,+1] to [0,1].
!
  x(1:n) = 0.5D+00 * ( x(1:n) + 1.0D+00 )
  w(1:n) = 0.5D+00 * w(1:n)

  return
end
subroutine ccu_order ( l, n )

!*****************************************************************************80
!
!! CCU_ORDER computes the order of a CCU rule from the level.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, the level of the rule.  
!    1 <= L.
!
!    Output, integer ( kind = 4 ) N, the order of the rule.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  if ( l < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCU_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop
  else if ( l == 1 ) then
    n = 1
  else
    n = 2 ** ( l - 1 ) + 1
  end if

  return
end
function fn_integral ( d )

!*****************************************************************************80
!
!! FN_INTEGRAL is the integral of the Hermite test function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) D, the spatial dimension.
!
!    Output, real ( kind = 8 ) FN_INTEGRAL, the integral value.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ), parameter :: exponent = 6
  real ( kind = 8 ) fn_integral
  integer ( kind = 4 ) i4_factorial2

  fn_integral = real ( i4_factorial2 ( exponent - 1 ), kind = 8 )

  return
end
subroutine fn_value ( d, n, x, fx )

!*****************************************************************************80
!
!! FN_VALUE is a Hermite test function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) D, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(D,N), the points.
!
!    Output, real ( kind = 8 ) FX(N), the function values.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: exponent = 6
  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) x(d,n)

  fx(1:n) = x(1,1:n)**exponent

  return
end
function fu_integral ( d )

!*****************************************************************************80
!
!! FU_INTEGRAL is the integral of the test function for the [0,1]^D interval.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) D, the spatial dimension.
!
!    Output, real ( kind = 8 ) FU_INTEGRAL, the integral value.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) fu_integral

  fu_integral = ( 0.5D+00 * erf ( 0.5D+00 / sqrt ( 2.0D+00 ) ) ) ** d

  return
end
subroutine fu_value ( d, n, x, fx )

!*****************************************************************************80
!
!! FU_VALUE is a sample function for the [0,1]^D interval.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) D, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(D,N), the points.
!
!    Output, real ( kind = 8 ) FX(N), the function values.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(d,n)

  fx(1:n) = 1.0D+00

  do i = 1, d
    fx(1:n) = fx(1:n) * exp ( - ( x(i,1:n) / 2.0D+00 )**2 / 2.0D+00 ) &
      / 2.0D+00 / sqrt ( 2.0D+00 * pi )
  end do

  return
end
subroutine get_seq ( d, norm, seq_num, fs )

!*****************************************************************************80
!
!! GET_SEQ generates all positive integer D-vectors that sum to NORM.
!
!  Discussion:
!
!    This function computes a list, in reverse dictionary order, of
!    all D-vectors of positive values that sum to NORM.
!
!    For example, call get_seq ( 3, 5, 6, fs ) returns
!
!      3  1  1
!      2  2  1
!      2  1  2
!      1  3  1
!      1  2  2
!      1  1  3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 May 2010
!
!  Author:
!
!    Original MATLAB version by Florian Heiss, Viktor Winschel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) D, the dimension.
!    1 <= D.
!
!    Input, integer ( kind = 4 ) NORM, the value that each row must sum to.
!    D <= NORM.
!
!    Input, integer ( kind = 4 ) SEQ_NUM, the number of rows of FS.
!
!    Output, integer ( kind = 4 ) FS(SEQ_NUM,D).  Each row of FS represents 
!    one vector with all elements positive and summing to NORM.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) seq_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) c
  integer ( kind = 4 ) fs(seq_num,d)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) norm
  integer ( kind = 4 ) row
  integer ( kind = 4 ) seq(d)

  if ( norm < d ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GET_SEQ - Fatal error!'
    write ( *, '(a,i4,a,i4)' ) '  NORM = ', norm, ' < D = ', d
    stop
  end if

  seq(1:d) = 0
!
!  The algorithm is written to work with vectors whose minimum value is
!  allowed to be zero.  So we subtract D from NORM at the beginning and
!  then increment the result vectors by 1 at the end!
!
  a = norm - d
  seq(1) = a

  row = 1
  fs(row,1:d) = seq(1:d) + 1
  c = 1

  do while ( seq ( d ) < a )

    if ( c == d ) then
      do i = c - 1, 1, -1
        c = i
        if ( seq(i) /= 0 ) then
          exit
        end if
      end do
    end if

    seq(c) = seq(c) - 1
    c = c + 1
    seq(c) = a - sum ( seq(1:(c-1)) )

    if ( c < d ) then
      seq((c+1):d) = 0
    end if

    row = row + 1
    fs(row,1:d) = seq(1:d) + 1

  end do

  return
end
subroutine gqn ( n, x, w )

!*****************************************************************************80
!
!! GQN provides data for Gauss quadrature with a normal weight.
!
!  Discussion:
!
!    This data assumes integration over the interval (-oo,+oo) with 
!    weight function w(x) = exp(-x*x/2)/sqrt(2*pi).
!
!    Because the rule is symmetric, for each N only (N+1)/2 values are returned.
!
!    For N odd, the actual rule has points and weights:
!      X(1) with weight W(1),
!      X(2) and -X(2), with weight W(2),
!      X(3) and -X(3), with weight W(3), and so on.
!
!    For N even, the actual rule has points and weights:
!      X(1) and -X(1), with weight W(1),
!      X(2) and -X(2), with weight W(2), and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2012
!
!  Author:
!
!    Original MATLAB version by Florian Heiss, Viktor Winschel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points and weights.
!    1 <= N <= 25.
!
!    Output, real ( kind = 8 ) X((N+1)/2), the nodes.
!
!    Output, real ( kind = 8 ) W((N+1)/2), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) nhalf
  real ( kind = 8 ) w(*)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ), dimension ( 1 ) :: x01 = (/ &
     0.0000000000000000D+00 /)
  real ( kind = 8 ), dimension ( 1 ) :: w01 = (/ &
    1.0000000000000000D+00 /)
  real ( kind = 8 ), dimension ( 1 ) :: x02 = (/ &
    1.0000000000000002D+00 /)
  real ( kind = 8 ), dimension ( 1 ) :: w02 = (/ &
    5.0000000000000000D-01 /)
  real ( kind = 8 ), dimension ( 2 ) :: x03 = (/ &
    0.0000000000000000D+00, 1.7320508075688772D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: w03 = (/ &
    6.6666666666666663D-01, 1.6666666666666674D-01 /)
  real ( kind = 8 ), dimension ( 2 ) :: x04 = (/ &
    7.4196378430272591D-01, 2.3344142183389773D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: w04 = (/ &
    4.5412414523193145D-01, 4.5875854768068498D-02 /)
  real ( kind = 8 ), dimension ( 3 ) :: x05 = (/ &
    0.0000000000000000D+00, 1.3556261799742659D+00, 2.8569700138728056D+00 /)
  real ( kind = 8 ), dimension ( 3 ) :: w05 = (/ &
    5.3333333333333344D-01, 2.2207592200561263D-01, 1.1257411327720691D-02 /)
  real ( kind = 8 ), dimension ( 3 ) :: x06 = (/ &
    6.1670659019259422D-01, 1.8891758777537109D+00, 3.3242574335521193D+00 /)
  real ( kind = 8 ), dimension ( 3 ) :: w06 = (/ &
    4.0882846955602919D-01, 8.8615746041914523D-02, 2.5557844020562431D-03 /)
  real ( kind = 8 ), dimension ( 4 ) :: x07 = (/ &
    0.0000000000000000D+00, 1.1544053947399682D+00, 2.3667594107345415D+00, &
    3.7504397177257425D+00 /)
  real ( kind = 8 ), dimension ( 4 ) :: w07 = (/ &
    4.5714285714285757D-01, 2.4012317860501250D-01, 3.0757123967586491D-02, &
    5.4826885597221875D-04 /)
  real ( kind = 8 ), dimension ( 4 ) :: x08 = (/ &
    5.3907981135137517D-01, 1.6365190424351082D+00, 2.8024858612875416D+00, &
    4.1445471861258945D+00 /)
  real ( kind = 8 ), dimension ( 4 ) :: w08 = (/ &
    3.7301225767907736D-01, 1.1723990766175897D-01, 9.6352201207882630D-03, &
    1.1261453837536784D-04 /)
  real ( kind = 8 ), dimension ( 5 ) :: x09 = (/ &
    0.0000000000000000D+00, 1.0232556637891326D+00, 2.0768479786778302D+00, &
    3.2054290028564703D+00, 4.5127458633997826D+00 /)
  real ( kind = 8 ), dimension ( 5 ) :: w09 = (/ &
    4.0634920634920685D-01, 2.4409750289493909D-01, 4.9916406765217969D-02, &
    2.7891413212317675D-03, 2.2345844007746563D-05 /)
  real ( kind = 8 ), dimension ( 5 ) :: x10 = (/ &
    4.8493570751549764D-01, 1.4659890943911582D+00, 2.4843258416389546D+00, &
    3.5818234835519269D+00, 4.8594628283323127D+00 /)
  real ( kind = 8 ), dimension ( 5 ) :: w10 = (/ &
    3.4464233493201940D-01, 1.3548370298026730D-01, 1.9111580500770317D-02, &
    7.5807093431221972D-04, 4.3106526307183106D-06 /)
  real ( kind = 8 ), dimension ( 6 ) :: x11 = (/ &
    0.0000000000000000D+00, 9.2886899738106388D-01, 1.8760350201548459D+00, &
    2.8651231606436447D+00, 3.9361666071299775D+00, 5.1880012243748714D+00 /)
  real ( kind = 8 ), dimension ( 6 ) :: w11 = (/ &
    3.6940836940836957D-01, 2.4224029987397003D-01, 6.6138746071057644D-02, &
    6.7202852355372697D-03, 1.9567193027122324D-04, 8.1218497902149036D-07 /)
  real ( kind = 8 ), dimension ( 6 ) :: x12 = (/ &
    4.4440300194413901D-01, 1.3403751971516167D+00, 2.2594644510007993D+00, &
    3.2237098287700974D+00, 4.2718258479322815D+00, 5.5009017044677480D+00 /)
  real ( kind = 8 ), dimension ( 6 ) :: w12 = (/ &
    3.2166436151283007D-01, 1.4696704804532995D-01, 2.9116687912364138D-02, &
    2.2033806875331849D-03, 4.8371849225906076D-05, 1.4999271676371597D-07 /)
  real ( kind = 8 ), dimension ( 7 ) :: x13 = (/ &
    0.0000000000000000D+00, 8.5667949351945005D-01, 1.7254183795882394D+00, &
    2.6206899734322149D+00, 3.5634443802816347D+00, 4.5913984489365207D+00, &
    5.8001672523865011D+00 /)
  real ( kind = 8 ), dimension ( 7 ) :: w13 = (/ &
    3.4099234099234149D-01, 2.3787152296413588D-01, 7.9168955860450141D-02, &
    1.1770560505996543D-02, 6.8123635044292619D-04, 1.1526596527333885D-05, &
    2.7226276428059039D-08 /)
  real ( kind = 8 ), dimension ( 7 ) :: x14 = (/ &
    4.1259045795460181D-01, 1.2426889554854643D+00, 2.0883447457019444D+00, &
    2.9630365798386675D+00, 3.8869245750597696D+00, 4.8969363973455646D+00, &
    6.0874095469012914D+00 /)
  real ( kind = 8 ), dimension ( 7 ) :: w14 = (/ &
    3.0263462681301945D-01, 1.5408333984251366D-01, 3.8650108824253432D-02, &
    4.4289191069474066D-03, 2.0033955376074381D-04, 2.6609913440676334D-06, &
    4.8681612577483872D-09 /)
  real ( kind = 8 ), dimension ( 8 ) :: x15 = (/ &
    0.0000000000000000D+00, 7.9912906832454811D-01, 1.6067100690287301D+00, &
    2.4324368270097581D+00, 3.2890824243987664D+00, 4.1962077112690155D+00, &
    5.1900935913047821D+00, 6.3639478888298378D+00 /)
  real ( kind = 8 ), dimension ( 8 ) :: w15 = (/ &
    3.1825951825951820D-01, 2.3246229360973222D-01, 8.9417795399844444D-02, &
    1.7365774492137616D-02, 1.5673575035499571D-03, 5.6421464051890157D-05, &
    5.9754195979205961D-07, 8.5896498996331805D-10 /)
  real ( kind = 8 ), dimension ( 8 ) :: x16 = (/ &
    3.8676060450055738D-01, 1.1638291005549648D+00, 1.9519803457163336D+00, &
    2.7602450476307019D+00, 3.6008736241715487D+00, 4.4929553025200120D+00, &
    5.4722257059493433D+00, 6.6308781983931295D+00 /)
  real ( kind = 8 ), dimension ( 8 ) :: w16 = (/ &
    2.8656852123801241D-01, 1.5833837275094925D-01, 4.7284752354014067D-02, &
    7.2669376011847411D-03, 5.2598492657390979D-04, 1.5300032162487286D-05, &
    1.3094732162868203D-07, 1.4978147231618314D-10 /)
  real ( kind = 8 ), dimension ( 9 ) :: x17 = (/ &
    0.0000000000000000D+00, 7.5184260070389630D-01, 1.5098833077967408D+00, &
    2.2810194402529889D+00, 3.0737971753281941D+00, 3.9000657171980104D+00, &
    4.7785315896299840D+00, 5.7444600786594071D+00, 6.8891224398953330D+00 /)
  real ( kind = 8 ), dimension ( 9 ) :: w17 = (/ &
    2.9953837012660756D-01, 2.2670630846897877D-01, 9.7406371162718081D-02, &
    2.3086657025711152D-02, 2.8589460622846499D-03, 1.6849143155133945D-04, &
    4.0126794479798725D-06, 2.8080161179305783D-08, 2.5843149193749151D-11 /)
  real ( kind = 8 ), dimension ( 9 ) :: x18 = (/ &
    3.6524575550769767D-01, 1.0983955180915013D+00, 1.8397799215086457D+00, &
    2.5958336889112403D+00, 3.3747365357780907D+00, 4.1880202316294044D+00, &
    5.0540726854427405D+00, 6.0077459113595975D+00, 7.1394648491464796D+00 /)
  real ( kind = 8 ), dimension ( 9 ) :: w18 = (/ &
    2.7278323465428789D-01, 1.6068530389351263D-01, 5.4896632480222654D-02, &
    1.0516517751941352D-02, 1.0654847962916496D-03, 5.1798961441161962D-05, &
    1.0215523976369816D-06, 5.9054884788365484D-09, 4.4165887693587078D-12 /)
  real ( kind = 8 ), dimension ( 10 ) :: x19 = (/ &
    0.0000000000000000D+00, 7.1208504404237993D-01, 1.4288766760783731D+00, &
    2.1555027613169351D+00, 2.8980512765157536D+00, 3.6644165474506383D+00, &
    4.4658726268310316D+00, 5.3205363773360386D+00, 6.2628911565132519D+00, &
    7.3825790240304316D+00 /)
  real ( kind = 8 ), dimension ( 10 ) :: w19 = (/ &
    2.8377319275152108D-01, 2.2094171219914366D-01, 1.0360365727614400D-01, &
    2.8666691030118496D-02, 4.5072354203420355D-03, 3.7850210941426759D-04, &
    1.5351145954666744D-05, 2.5322200320928681D-07, 1.2203708484474786D-09, &
    7.4828300540572308D-13 /)
  real ( kind = 8 ), dimension ( 10 ) :: x20 = (/ &
    3.4696415708135592D-01, 1.0429453488027509D+00, 1.7452473208141270D+00, &
    2.4586636111723679D+00, 3.1890148165533900D+00, 3.9439673506573163D+00, &
    4.7345813340460552D+00, 5.5787388058932015D+00, 6.5105901570136551D+00, &
    7.6190485416797591D+00 /)
  real ( kind = 8 ), dimension ( 10 ) :: w20 = (/ &
    2.6079306344955544D-01, 1.6173933398400026D-01, 6.1506372063976029D-02, &
    1.3997837447101043D-02, 1.8301031310804918D-03, 1.2882627996192898D-04, &
    4.4021210902308646D-06, 6.1274902599829597D-08, 2.4820623623151838D-10, &
    1.2578006724379305D-13 /)
  real ( kind = 8 ), dimension ( 11 ) :: x21 = (/ &
    0.0000000000000000D+00, 6.7804569244064405D-01, 1.3597658232112304D+00, &
    2.0491024682571628D+00, 2.7505929810523733D+00, 3.4698466904753764D+00, &
    4.2143439816884216D+00, 4.9949639447820253D+00, 5.8293820073044706D+00, &
    6.7514447187174609D+00, 7.8493828951138225D+00 /)
  real ( kind = 8 ), dimension ( 11 ) :: w21 = (/ &
    2.7026018357287707D-01, 2.1533371569505982D-01, 1.0839228562641938D-01, &
    3.3952729786542839D-02, 6.4396970514087768D-03, 7.0804779548153736D-04, &
    4.2192347425515866D-05, 1.2253548361482522D-06, 1.4506612844930740D-08, &
    4.9753686041217464D-11, 2.0989912195656652D-14 /)
  real ( kind = 8 ), dimension ( 11 ) :: x22 = (/ &
    3.3117931571527381D-01, 9.9516242227121554D-01, 1.6641248391179071D+00, &
    2.3417599962877080D+00, 3.0324042278316763D+00, 3.7414963502665177D+00, &
    4.4763619773108685D+00, 5.2477244337144251D+00, 6.0730749511228979D+00, &
    6.9859804240188152D+00, 8.0740299840217116D+00 /)
  real ( kind = 8 ), dimension ( 11 ) :: w22 = (/ &
    2.5024359658693501D-01, 1.6190629341367538D-01, 6.7196311428889891D-02, &
    1.7569072880805774D-02, 2.8087610475772107D-03, 2.6228330325596416D-04, &
    1.3345977126808712D-05, 3.3198537498140043D-07, 3.3665141594582109D-09, &
    9.8413789823460105D-12, 3.4794606478771428D-15 /)
  real ( kind = 8 ), dimension ( 12 ) :: x23 = (/ &
    0.0000000000000000D+00, 6.4847115353449580D-01, 1.2998764683039790D+00, &
    1.9573275529334242D+00, 2.6243236340591820D+00, 3.3050400217529652D+00, &
    4.0047753217333044D+00, 4.7307241974514733D+00, 5.4934739864717947D+00, &
    6.3103498544483996D+00, 7.2146594350518622D+00, 8.2933860274173536D+00 /)
  real ( kind = 8 ), dimension ( 12 ) :: w23 = (/ &
    2.5850974080883904D-01, 2.0995966957754261D-01, 1.1207338260262091D-01, &
    3.8867183703480947D-02, 8.5796783914656640D-03, 1.1676286374978613D-03, &
    9.3408186090312983D-05, 4.0899772449921549D-06, 8.7750624838617161D-08, &
    7.6708888623999076D-10, 1.9229353115677913D-12, 5.7323831678020873D-16 /)
  real ( kind = 8 ), dimension ( 12 ) :: x24 = (/ &
    3.1737009662945231D-01, 9.5342192293210926D-01, 1.5934804298164202D+00, &
    2.2404678516917524D+00, 2.8977286432233140D+00, 3.5693067640735610D+00, &
    4.2603836050199053D+00, 4.9780413746391208D+00, 5.7327471752512009D+00, &
    6.5416750050986341D+00, 7.4378906660216630D+00, 8.5078035191952583D+00 /)
  real ( kind = 8 ), dimension ( 12 ) :: w24 = (/ &
    2.4087011554664056D-01, 1.6145951286700025D-01, 7.2069364017178436D-02, &
    2.1126344408967029D-02, 3.9766089291813113D-03, 4.6471871877939763D-04, &
    3.2095005652745989D-05, 1.2176597454425830D-06, 2.2674616734804651D-08, &
    1.7186649279648690D-10, 3.7149741527624159D-13, 9.3901936890419202D-17 /)
  real ( kind = 8 ), dimension ( 13 ) :: x25 = (/ &
    0.0000000000000000D+00, 6.2246227918607611D-01, 1.2473119756167892D+00, &
    1.8770583699478387D+00, 2.5144733039522058D+00, 3.1627756793881927D+00, &
    3.8259005699724917D+00, 4.5089299229672850D+00, 5.2188480936442794D+00, &
    5.9660146906067020D+00, 6.7674649638097168D+00, 7.6560379553930762D+00, &
    8.7175976783995885D+00 /)
  real ( kind = 8 ), dimension ( 13 ) :: w25 = (/ &
    2.4816935117648548D-01, 2.0485102565034041D-01, 1.1488092430395164D-01, &
    4.3379970167644971D-02, 1.0856755991462316D-02, 1.7578504052637961D-03, &
    1.7776690692652660D-04, 1.0672194905202536D-05, 3.5301525602454978D-07, &
    5.7380238688993763D-09, 3.7911500004771871D-11, 7.1021030370039253D-14, &
    1.5300389979986825D-17 /)

  nhalf = ( n + 1 ) / 2

  if ( n == 1 ) then
    call r8vec_copy ( nhalf, x01, x )
    call r8vec_copy ( nhalf, w01, w )
  else if ( n == 2 ) then
    call r8vec_copy ( nhalf, x02, x )
    call r8vec_copy ( nhalf, w02, w )
  else if ( n == 3 ) then
    call r8vec_copy ( nhalf, x03, x )
    call r8vec_copy ( nhalf, w03, w )
  else if ( n == 4 ) then
    call r8vec_copy ( nhalf, x04, x )
    call r8vec_copy ( nhalf, w04, w )
  else if ( n == 5 ) then
    call r8vec_copy ( nhalf, x05, x )
    call r8vec_copy ( nhalf, w05, w )
  else if ( n == 6 ) then
    call r8vec_copy ( nhalf, x06, x )
    call r8vec_copy ( nhalf, w06, w )
  else if ( n == 7 ) then
    call r8vec_copy ( nhalf, x07, x )
    call r8vec_copy ( nhalf, w07, w )
  else if ( n == 8 ) then
    call r8vec_copy ( nhalf, x08, x )
    call r8vec_copy ( nhalf, w08, w )
  else if ( n == 9 ) then
    call r8vec_copy ( nhalf, x09, x )
    call r8vec_copy ( nhalf, w09, w )
  else if ( n == 10 ) then
    call r8vec_copy ( nhalf, x10, x )
    call r8vec_copy ( nhalf, w10, w )
  else if ( n == 11 ) then
    call r8vec_copy ( nhalf, x11, x )
    call r8vec_copy ( nhalf, w11, w )
  else if ( n == 12 ) then
    call r8vec_copy ( nhalf, x12, x )
    call r8vec_copy ( nhalf, w12, w )
  else if ( n == 13 ) then
    call r8vec_copy ( nhalf, x13, x )
    call r8vec_copy ( nhalf, w13, w )
  else if ( n == 14 ) then
    call r8vec_copy ( nhalf, x14, x )
    call r8vec_copy ( nhalf, w14, w )
  else if ( n == 15 ) then
    call r8vec_copy ( nhalf, x15, x )
    call r8vec_copy ( nhalf, w15, w )
  else if ( n == 16 ) then
    call r8vec_copy ( nhalf, x16, x )
    call r8vec_copy ( nhalf, w16, w )
  else if ( n == 17 ) then
    call r8vec_copy ( nhalf, x17, x )
    call r8vec_copy ( nhalf, w17, w )
  else if ( n == 18 ) then
    call r8vec_copy ( nhalf, x18, x )
    call r8vec_copy ( nhalf, w18, w )
  else if ( n == 19 ) then
    call r8vec_copy ( nhalf, x19, x )
    call r8vec_copy ( nhalf, w19, w )
  else if ( n == 20 ) then
    call r8vec_copy ( nhalf, x20, x )
    call r8vec_copy ( nhalf, w20, w )
  else if ( n == 21 ) then
    call r8vec_copy ( nhalf, x21, x )
    call r8vec_copy ( nhalf, w21, w )
  else if ( n == 22 ) then
    call r8vec_copy ( nhalf, x22, x )
    call r8vec_copy ( nhalf, w22, w )
  else if ( n == 23 ) then
    call r8vec_copy ( nhalf, x23, x )
    call r8vec_copy ( nhalf, w23, w )
  else if ( n == 24 ) then
    call r8vec_copy ( nhalf, x24, x )
    call r8vec_copy ( nhalf, w24, w )
  else if ( n == 25 ) then
    call r8vec_copy ( nhalf, x25, x )
    call r8vec_copy ( nhalf, w25, w )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQN - Fatal error!'
    write ( *, '(a)' ) '  Value of N must be between 1 and 25.'
    stop
  end if

  return
end
subroutine gqn_order ( l, n )

!*****************************************************************************80
!
!! GQN_ORDER computes the order of a GQN rule from the level.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, the level of the rule.  
!    1 <= L.
!
!    Output, integer ( kind = 4 ) N, the order of the rule.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  if ( l < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQN_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop
  else if ( l <= 25 ) then
    n = l
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQN_ORDER - Fatal error!'
    write ( *, '(a)' ) '  L <= 25 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop
  end if

  return
end
subroutine gqu ( n, x, w )

!*****************************************************************************80
!
!! GQU provides data for Gauss quadrature with a uniform weight.
!
!  Discussion:
!
!    This data assumes integration over the interval [0,1] with 
!    weight function w(x) = 1.
!
!    Because the rule is symmetric, for each N only (N+1)/2 values are returned.
!
!    For N odd, the actual rule has points and weights:
!      X(1) with weight W(1),
!      X(2) and -X(2), with weight W(2),
!      X(3) and -X(3), with weight W(3), and so on.
!
!    For N even, the actual rule has points and weights:
!      X(1) and -X(1), with weight W(1),
!      X(2) and -X(2), with weight W(2), and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2012
!
!  Author:
!
!    Original MATLAB version by Florian Heiss, Viktor Winschel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points and weights.
!    1 <= N <= 25.
!
!    Output, real ( kind = 8 ) X((N+1)/2), the nodes.
!
!    Output, real ( kind = 8 ) W((N+1)/2), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) nhalf
  real ( kind = 8 ) w(*)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ), dimension ( 1 ) :: x01 = (/ &
    5.0000000000000000D-01 /)
  real ( kind = 8 ), dimension ( 1 ) :: w01 = (/ &
    1.0000000000000000D+00 /)
  real ( kind = 8 ), dimension ( 1 ) :: x02 = (/ &
    7.8867513459481287D-01 /)
  real ( kind = 8 ), dimension ( 1 ) :: w02 = (/ &
    5.0000000000000000D-01 /)
  real ( kind = 8 ), dimension ( 2 ) :: x03 = (/ &
    5.0000000000000000D-01, 8.8729833462074170D-01 /)
  real ( kind = 8 ), dimension ( 2 ) :: w03 = (/ &
    4.4444444444444570D-01, 2.7777777777777712D-01 /)
  real ( kind = 8 ), dimension ( 2 ) :: x04 = (/ &
    6.6999052179242813D-01, 9.3056815579702623D-01 /)
  real ( kind = 8 ), dimension ( 2 ) :: w04 = (/ &
    3.2607257743127516D-01, 1.7392742256872484D-01 /)
  real ( kind = 8 ), dimension ( 3 ) :: x05 = (/ &
    5.0000000000000000D-01, 7.6923465505284150D-01, 9.5308992296933193D-01 /)
  real ( kind = 8 ), dimension ( 3 ) :: w05 = (/ &
    2.8444444444444655D-01, 2.3931433524968501D-01, 1.1846344252809174D-01 /)
  real ( kind = 8 ), dimension ( 3 ) :: x06 = (/ &
    6.1930959304159849D-01, 8.3060469323313235D-01, 9.6623475710157603D-01 /)
  real ( kind = 8 ), dimension ( 3 ) :: w06 = (/ &
    2.3395696728634746D-01, 1.8038078652407072D-01, 8.5662246189581834D-02 /)
  real ( kind = 8 ), dimension ( 4 ) :: x07 = (/ &
    5.0000000000000000D-01, 7.0292257568869854D-01, 8.7076559279969723D-01, &
    9.7455395617137919D-01 /)
  real ( kind = 8 ), dimension ( 4 ) :: w07 = (/ &
    2.0897959183673620D-01, 1.9091502525256090D-01, 1.3985269574463935D-01, &
    6.4742483084431701D-02 /)
  real ( kind = 8 ), dimension ( 4 ) :: x08 = (/ &
    5.9171732124782495D-01, 7.6276620495816450D-01, 8.9833323870681348D-01, &
    9.8014492824876809D-01 /)
  real ( kind = 8 ), dimension ( 4 ) :: w08 = (/ &
    1.8134189168918213D-01, 1.5685332293894469D-01, 1.1119051722668793D-01, &
    5.0614268145185180D-02 /)
  real ( kind = 8 ), dimension ( 5 ) :: x09 = (/ &
    5.0000000000000000D-01, 6.6212671170190451D-01, 8.0668571635029518D-01, &
    9.1801555366331788D-01, 9.8408011975381304D-01 /)
  real ( kind = 8 ), dimension ( 5 ) :: w09 = (/ &
    1.6511967750063075D-01, 1.5617353852000226D-01, 1.3030534820146844D-01, &
    9.0324080347429253D-02, 4.0637194180784583D-02 /)
  real ( kind = 8 ), dimension ( 5 ) :: x10 = (/ &
    5.7443716949081558D-01, 7.1669769706462361D-01, 8.3970478414951222D-01, &
    9.3253168334449232D-01, 9.8695326425858587D-01 /)
  real ( kind = 8 ), dimension ( 5 ) :: w10 = (/ &
    1.4776211235737713D-01, 1.3463335965499873D-01, 1.0954318125799158D-01, &
    7.4725674575290599D-02, 3.3335672154342001D-02 /)
  real ( kind = 8 ), dimension ( 6 ) :: x11 = (/ &
    5.0000000000000000D-01, 6.3477157797617245D-01, 7.5954806460340585D-01, &
    8.6507600278702468D-01, 9.4353129988404771D-01, 9.8911432907302843D-01 /)
  real ( kind = 8 ), dimension ( 6 ) :: w11 = (/ &
    1.3646254338895086D-01, 1.3140227225512388D-01, 1.1659688229599563D-01, &
    9.3145105463867520D-02, 6.2790184732452625D-02, 2.7834283558084916D-02 /)
  real ( kind = 8 ), dimension ( 6 ) :: x12 = (/ &
    5.6261670425573451D-01, 6.8391574949909006D-01, 7.9365897714330869D-01, &
    8.8495133709715235D-01, 9.5205862818523745D-01, 9.9078031712335957D-01 /)
  real ( kind = 8 ), dimension ( 6 ) :: w12 = (/ &
    1.2457352290670189D-01, 1.1674626826917781D-01, 1.0158371336153328D-01, &
    8.0039164271673444D-02, 5.3469662997659276D-02, 2.3587668193254314D-02 /)
  real ( kind = 8 ), dimension ( 7 ) :: x13 = (/ &
    5.0000000000000000D-01, 6.1522915797756739D-01, 7.2424637551822335D-01, &
    8.2117466972017006D-01, 9.0078904536665494D-01, 9.5879919961148907D-01, &
    9.9209152735929407D-01 /)
  real ( kind = 8 ), dimension ( 7 ) :: w13 = (/ &
    1.1627577661543741D-01, 1.1314159013144903D-01, 1.0390802376844462D-01, &
    8.9072990380973202D-02, 6.9436755109893875D-02, 4.6060749918864378D-02, &
    2.0242002382656228D-02 /)
  real ( kind = 8 ), dimension ( 7 ) :: x14 = (/ &
    5.5402747435367183D-01, 6.5955618446394482D-01, 7.5762431817907705D-01, &
    8.4364645240584268D-01, 9.1360065753488251D-01, 9.6421744183178681D-01, &
    9.9314190434840621D-01 /)
  real ( kind = 8 ), dimension ( 7 ) :: w14 = (/ &
    1.0763192673157916D-01, 1.0259923186064811D-01, 9.2769198738969161D-02, &
    7.8601583579096995D-02, 6.0759285343951711D-02, 4.0079043579880291D-02, &
    1.7559730165874574D-02 /)
  real ( kind = 8 ), dimension ( 8 ) :: x15 = (/ &
    5.0000000000000000D-01, 6.0059704699871730D-01, 6.9707567353878175D-01, &
    7.8548608630426942D-01, 8.6220886568008503D-01, 9.2410329170521366D-01, &
    9.6863669620035298D-01, 9.9399625901024269D-01 /)
  real ( kind = 8 ), dimension ( 8 ) :: w15 = (/ &
    1.0128912096278091D-01, 9.9215742663556039D-02, 9.3080500007781286D-02, &
    8.3134602908497196D-02, 6.9785338963077315D-02, 5.3579610233586157D-02, &
    3.5183023744054159D-02, 1.5376620998057434D-02 /)
  real ( kind = 8 ), dimension ( 8 ) :: x16 = (/ &
    5.4750625491881877D-01, 6.4080177538962946D-01, 7.2900838882861363D-01, &
    8.0893812220132189D-01, 8.7770220417750155D-01, 9.3281560119391593D-01, &
    9.7228751153661630D-01, 9.9470046749582497D-01 /)
  real ( kind = 8 ), dimension ( 8 ) :: w16 = (/ &
    9.4725305227534431D-02, 9.1301707522462000D-02, 8.4578259697501462D-02, &
    7.4797994408288562D-02, 6.2314485627767105D-02, 4.7579255841246545D-02, &
    3.1126761969323954D-02, 1.3576229705875955D-02 /)
  real ( kind = 8 ), dimension ( 9 ) :: x17 = (/ &
    5.0000000000000000D-01, 5.8924209074792389D-01, 6.7561588172693821D-01, &
    7.5634526854323847D-01, 8.2883557960834531D-01, 8.9075700194840068D-01, &
    9.4011957686349290D-01, 9.7533776088438384D-01, 9.9528773765720868D-01 /)
  real ( kind = 8 ), dimension ( 9 ) :: w17 = (/ &
    8.9723235178103419D-02, 8.8281352683496447D-02, 8.4002051078225143D-02, &
    7.7022880538405308D-02, 6.7568184234262890D-02, 5.5941923596702053D-02, &
    4.2518074158589644D-02, 2.7729764686993612D-02, 1.2074151434273140D-02 /)
  real ( kind = 8 ), dimension ( 9 ) :: x18 = (/ &
    5.4238750652086765D-01, 6.2594311284575277D-01, 7.0587558073142131D-01, &
    7.7988541553697377D-01, 8.4584352153017661D-01, 9.0185247948626157D-01, &
    9.4630123324877791D-01, 9.7791197478569880D-01, 9.9578258421046550D-01 /)
  real ( kind = 8 ), dimension ( 9 ) :: w18 = (/ &
    8.4571191481571939D-02, 8.2138241872916504D-02, 7.7342337563132801D-02, &
    7.0321457335325452D-02, 6.1277603355739306D-02, 5.0471022053143716D-02, &
    3.8212865127444665D-02, 2.4857274447484968D-02, 1.0808006763240719D-02 /)
  real ( kind = 8 ), dimension ( 10 ) :: x19 = (/ &
    5.0000000000000000D-01, 5.8017932282011264D-01, 6.5828204998181494D-01, &
    7.3228537068798050D-01, 8.0027265233084055D-01, 8.6048308866761469D-01, &
    9.1135732826857141D-01, 9.5157795180740901D-01, 9.8010407606741501D-01, &
    9.9620342192179212D-01 /)
  real ( kind = 8 ), dimension ( 10 ) :: w19 = (/ &
    8.0527224924391946D-02, 7.9484421696977337D-02, 7.6383021032929960D-02, &
    7.1303351086803413D-02, 6.4376981269668232D-02, 5.5783322773667113D-02, &
    4.5745010811225124D-02, 3.4522271368820669D-02, 2.2407113382849821D-02, &
    9.7308941148624341D-03 /)
  real ( kind = 8 ), dimension ( 10 ) :: x20 = (/ &
    5.3826326056674867D-01, 6.1389292557082253D-01, 6.8685304435770977D-01, &
    7.5543350097541362D-01, 8.1802684036325757D-01, 8.7316595323007540D-01, &
    9.1955848591110945D-01, 9.5611721412566297D-01, 9.8198596363895696D-01, &
    9.9656429959254744D-01 /)
  real ( kind = 8 ), dimension ( 10 ) :: w20 = (/ &
    7.6376693565363113D-02, 7.4586493236301996D-02, 7.1048054659191187D-02, &
    6.5844319224588346D-02, 5.9097265980759248D-02, 5.0965059908620318D-02, &
    4.1638370788352433D-02, 3.1336024167054569D-02, 2.0300714900193556D-02, &
    8.8070035695753026D-03 /)
  real ( kind = 8 ), dimension ( 11 ) :: x21 = (/ &
    5.0000000000000000D-01, 5.7278092708044759D-01, 6.4401065840120053D-01, &
    7.1217106010371944D-01, 7.7580941794360991D-01, 8.3356940209870611D-01, &
    8.8421998173783889D-01, 9.2668168229165859D-01, 9.6004966707520034D-01, &
    9.8361341928315316D-01, 9.9687608531019478D-01 /)
  real ( kind = 8 ), dimension ( 11 ) :: w21 = (/ &
    7.3040566824845346D-02, 7.2262201994985134D-02, 6.9943697395536658D-02, &
    6.6134469316668845D-02, 6.0915708026864350D-02, 5.4398649583574356D-02, &
    4.6722211728016994D-02, 3.8050056814189707D-02, 2.8567212713428641D-02, &
    1.8476894885426285D-02, 8.0086141288864491D-03 /)
  real ( kind = 8 ), dimension ( 11 ) :: x22 = (/ &
    5.3486963665986109D-01, 6.0393021334411068D-01, 6.7096791044604209D-01, &
    7.3467791899337853D-01, 7.9382020175345580D-01, 8.4724363159334137D-01, &
    8.9390840298960406D-01, 9.3290628886015003D-01, 9.6347838609358694D-01, &
    9.8503024891771429D-01, 9.9714729274119962D-01 /)
  real ( kind = 8 ), dimension ( 11 ) :: w22 = (/ &
    6.9625936427816129D-02, 6.8270749173007697D-02, 6.5586752393531317D-02, &
    6.1626188405256251D-02, 5.6466148040269712D-02, 5.0207072221440600D-02, &
    4.2970803108533975D-02, 3.4898234212260300D-02, 2.6146667576341692D-02, &
    1.6887450792407110D-02, 7.3139976491353280D-03 /)
  real ( kind = 8 ), dimension ( 12 ) :: x23 = (/ &
    5.0000000000000000D-01, 5.6662841214923310D-01, 6.3206784048517251D-01, &
    6.9515051901514546D-01, 7.5475073892300371D-01, 8.0980493788182306D-01, &
    8.5933068156597514D-01, 9.0244420080942001D-01, 9.3837617913522076D-01, &
    9.6648554341300807D-01, 9.8627123560905761D-01, 9.9738466749877608D-01 /)
  real ( kind = 8 ), dimension ( 12 ) :: w23 = (/ &
    6.6827286093053176D-02, 6.6231019702348404D-02, 6.4452861094041150D-02, &
    6.1524542153364815D-02, 5.7498320111205814D-02, 5.2446045732270824D-02, &
    4.6457883030017563D-02, 3.9640705888359551D-02, 3.2116210704262994D-02, &
    2.4018835865542369D-02, 1.5494002928489686D-02, 6.7059297435702412D-03 /)
  real ( kind = 8 ), dimension ( 12 ) :: x24 = (/ &
    5.3202844643130276D-01, 5.9555943373680820D-01, 6.5752133984808170D-01, &
    7.1689675381302254D-01, 7.7271073569441984D-01, 8.2404682596848777D-01, &
    8.7006209578927718D-01, 9.1000099298695147D-01, 9.4320776350220048D-01, &
    9.6913727600136634D-01, 9.8736427798565474D-01, 9.9759360999851066D-01 /)
  real ( kind = 8 ), dimension ( 12 ) :: w24 = (/ &
    6.3969097673376246D-02, 6.2918728173414318D-02, 6.0835236463901793D-02, &
    5.7752834026862883D-02, 5.3722135057982914D-02, 4.8809326052057039D-02, &
    4.3095080765976693D-02, 3.6673240705540205D-02, 2.9649292457718385D-02, &
    2.2138719408709880D-02, 1.4265694314466934D-02, 6.1706148999928351D-03 /)
  real ( kind = 8 ), dimension ( 13 ) :: x25 = (/ &
    5.0000000000000000D-01, 5.6143234630535521D-01, 6.2193344186049426D-01, &
    6.8058615290469393D-01, 7.3650136572285752D-01, 7.8883146512061142D-01, &
    8.3678318423673415D-01, 8.7962963151867890D-01, 9.1672131438041693D-01, &
    9.4749599893913761D-01, 9.7148728561448716D-01, 9.8833196072975871D-01, &
    9.9777848489524912D-01 /)
  real ( kind = 8 ), dimension ( 13 ) :: w25 = (/ &
    6.1588026863357799D-02, 6.1121221495155122D-02, 5.9727881767892461D-02, &
    5.7429129572855862D-02, 5.4259812237131867D-02, 5.0267974533525363D-02, &
    4.5514130991481903D-02, 4.0070350167500532D-02, 3.4019166906178545D-02, &
    2.7452347987917691D-02, 2.0469578350653148D-02, 1.3177493307516108D-02, &
    5.6968992505125535D-03 /)

  nhalf = ( n + 1 ) / 2

  if ( n == 1 ) then
    call r8vec_copy ( nhalf, x01, x )
    call r8vec_copy ( nhalf, w01, w )
  else if ( n == 2 ) then
    call r8vec_copy ( nhalf, x02, x )
    call r8vec_copy ( nhalf, w02, w )
  else if ( n == 3 ) then
    call r8vec_copy ( nhalf, x03, x )
    call r8vec_copy ( nhalf, w03, w )
  else if ( n == 4 ) then
    call r8vec_copy ( nhalf, x04, x )
    call r8vec_copy ( nhalf, w04, w )
  else if ( n == 5 ) then
    call r8vec_copy ( nhalf, x05, x )
    call r8vec_copy ( nhalf, w05, w )
  else if ( n == 6 ) then
    call r8vec_copy ( nhalf, x06, x )
    call r8vec_copy ( nhalf, w06, w )
  else if ( n == 7 ) then
    call r8vec_copy ( nhalf, x07, x )
    call r8vec_copy ( nhalf, w07, w )
  else if ( n == 8 ) then
    call r8vec_copy ( nhalf, x08, x )
    call r8vec_copy ( nhalf, w08, w )
  else if ( n == 9 ) then
    call r8vec_copy ( nhalf, x09, x )
    call r8vec_copy ( nhalf, w09, w )
  else if ( n == 10 ) then
    call r8vec_copy ( nhalf, x10, x )
    call r8vec_copy ( nhalf, w10, w )
  else if ( n == 11 ) then
    call r8vec_copy ( nhalf, x11, x )
    call r8vec_copy ( nhalf, w11, w )
  else if ( n == 12 ) then
    call r8vec_copy ( nhalf, x12, x )
    call r8vec_copy ( nhalf, w12, w )
  else if ( n == 13 ) then
    call r8vec_copy ( nhalf, x13, x )
    call r8vec_copy ( nhalf, w13, w )
  else if ( n == 14 ) then
    call r8vec_copy ( nhalf, x14, x )
    call r8vec_copy ( nhalf, w14, w )
  else if ( n == 15 ) then
    call r8vec_copy ( nhalf, x15, x )
    call r8vec_copy ( nhalf, w15, w )
  else if ( n == 16 ) then
    call r8vec_copy ( nhalf, x16, x )
    call r8vec_copy ( nhalf, w16, w )
  else if ( n == 17 ) then
    call r8vec_copy ( nhalf, x17, x )
    call r8vec_copy ( nhalf, w17, w )
  else if ( n == 18 ) then
    call r8vec_copy ( nhalf, x18, x )
    call r8vec_copy ( nhalf, w18, w )
  else if ( n == 19 ) then
    call r8vec_copy ( nhalf, x19, x )
    call r8vec_copy ( nhalf, w19, w )
  else if ( n == 20 ) then
    call r8vec_copy ( nhalf, x20, x )
    call r8vec_copy ( nhalf, w20, w )
  else if ( n == 21 ) then
    call r8vec_copy ( nhalf, x21, x )
    call r8vec_copy ( nhalf, w21, w )
  else if ( n == 22 ) then
    call r8vec_copy ( nhalf, x22, x )
    call r8vec_copy ( nhalf, w22, w )
  else if ( n == 23 ) then
    call r8vec_copy ( nhalf, x23, x )
    call r8vec_copy ( nhalf, w23, w )
  else if ( n == 24 ) then
    call r8vec_copy ( nhalf, x24, x )
    call r8vec_copy ( nhalf, w24, w )
  else if ( n == 25 ) then
    call r8vec_copy ( nhalf, x25, x )
    call r8vec_copy ( nhalf, w25, w )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQU - Fatal error!'
    write ( *, '(a)' ) '  Value of N must be between 1 and 25.'
    stop
  end if

  return
end
subroutine gqu_order ( l, n )

!*****************************************************************************80
!
!! GQU_ORDER computes the order of a GQU rule from the level.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, the level of the rule.  
!    1 <= L <= 25.
!
!    Output, integer ( kind = 4 ) N, the order of the rule.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  if ( l < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQU_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop
  else if ( l <= 25 ) then
    n = l
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQU_ORDER - Fatal error!'
    write ( *, '(a)' ) '  L <= 25 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop
  end if

  return
end
function i4_choose ( n, k )

!*****************************************************************************80
!
!! I4_CHOOSE computes the binomial coefficient C(N,K) as an I4.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, integer ( kind = 4 ) I4_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    mx = max ( k, n - k )
    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  i4_choose = value

  return
end
function i4_factorial2 ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    The formula is:
!
!      FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                      = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    Factorial2(N)
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the double factorial 
!    function.  If N is less than 1, I4_FACTORIAL2 is returned as 1.
!
!    Output, integer ( kind = 4 ) I4_FACTORIAL2, the value of the function..
!
  implicit none

  integer ( kind = 4 ) i4_factorial2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_copy

  if ( n < 1 ) then
    i4_factorial2 = 1
    return
  end if

  n_copy = n
  i4_factorial2 = 1

  do while ( 1 < n_copy )
    i4_factorial2 = i4_factorial2 * n_copy
    n_copy = n_copy - 2
  end do

  return
end
function i4_mop ( i )

!*****************************************************************************80
!
!! I4_MOP returns the I-th power of -1 as an I4 value.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the power of -1.
!
!    Output, integer ( kind = 4 ) I4_MOP, the I-th power of -1.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_mop

  if ( mod ( i, 2 ) == 0 ) then
    i4_mop = 1
  else
    i4_mop = -1
  end if

  return
end
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2003
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
!    Input, integer ( kind = 4 ) A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

  return
end
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 8 ) ctemp(incx)
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

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8)' ) j
    end do

    write ( *, '(''  Col '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine i4vec_cum0 ( n, a, a_cum )

!*****************************************************************************80
!
!! I4VEC_CUM0 computes the cumulutive sum of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    This routine returns a vector of length N+1, with the first value
!    being 0.
!
!  Example:
!
!    Input:
!
!      A = (/ 1, 2, 3, 4 /)
!
!    Output:
!
!      A_CUM = (/ 0, 1, 3, 6, 10 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be summed.
!
!    Output, integer ( kind = 4 ) A_CUM(0:N), the cumulative sum of the
!    entries of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_cum(0:n)
  integer ( kind = 4 ) i

  a_cum(0) = 0

  do i = 1, n
    a_cum(i) = a_cum(i-1) + a(i)
  end do

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2010
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
function i4vec_product ( n, a )

!*****************************************************************************80
!
!! I4VEC_PRODUCT returns the product of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    In FORTRAN90, this facility is offered by the built in
!    PRODUCT function:
!
!      I4VEC_PRODUCT ( N, A ) = PRODUCT ( A(1:N) )
!
!    In MATLAB, this facility is offered by the built in
!    PROD function:
!
!      I4VEC_PRODUCT ( N, A ) = PROD ( A(1:N) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) I4VEC_PRODUCT, the product of the entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i4vec_product

  i4vec_product = product ( a(1:n) )

  return
end
function i4vec_sum ( n, a )

!*****************************************************************************80
!
!! I4VEC_SUM returns the sum of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    In FORTRAN90, this facility is offered by the built in
!    SUM function:
!
!      I4VEC_SUM ( N, A ) = SUM ( A(1:N) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) I4VEC_SUM, the sum of the entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i4vec_sum

  i4vec_sum = sum ( a(1:n) )

  return
end
subroutine i4vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
!    TITLE = 'My vector:  '
!
!    My vector:
!
!        1    2    3    4    5
!        6    7    8    9   10
!       11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2005
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5i12)' ) a(ilo:ihi)
  end do

  return
end
subroutine kpn ( n, x, w )

!*****************************************************************************80
!
!! KPN provides data for Kronrod-Patterson quadrature with a normal weight.
!
!  Discussion:
!
!    This data assumes integration over the interval (-oo,+oo) with 
!    weight function w(x) = exp(-x*x/2)/sqrt(2*pi).
!
!    For all orders N, the rule is formed by
!      X(1) with weight W(1),
!      X(2) and -X(2) with weight W(2),
!      X(3) and -X(3) with weight W(3) and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    Original MATLAB version by Florian Heiss, Viktor Winschel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309.
!
!    Thomas Patterson,
!    The optimal addition of points to quadrature formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) X((N+1)/2), the nodes.
!
!    Output, real ( kind = 8 ) W((N+1)/2), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) nhalf
  real ( kind = 8 ) w(*)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ), dimension ( 1 ) :: x01 = (/ &
    0.0000000000000000D+00 /)
  real ( kind = 8 ), dimension ( 1 ) :: w01 = (/ &
    1.0000000000000000D+00/)
  real ( kind = 8 ), dimension ( 2 ) :: x03 = (/ &
    0.0000000000000000D+00, 1.7320508075688772D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: w03 = (/ &
    6.6666666666666663D-01, 1.6666666666666666D-01/)
  real ( kind = 8 ), dimension ( 4 ) :: x07 = (/ &
    0.0000000000000000D+00, 7.4109534999454085D-01, 1.7320508075688772D+00, &
    4.1849560176727323D+00 /)
  real ( kind = 8 ), dimension ( 4 ) :: w07 = (/ &
    4.5874486825749189D-01, 1.3137860698313561D-01, 1.3855327472974924D-01, &
    6.9568415836913987D-04 /)
  real ( kind = 8 ), dimension ( 5 ) :: x09 = (/ &
    0.0000000000000000D+00, 7.4109534999454085D-01, 1.7320508075688772D+00, &
    2.8612795760570582D+00, 4.1849560176727323D+00 /)
  real ( kind = 8 ), dimension ( 5 ) :: w09 = (/ &
    2.5396825396825407D-01, 2.7007432957793776D-01, 9.4850948509485125D-02, &
    7.9963254708935293D-03, 9.4269457556517470D-05 /)
  real ( kind = 8 ), dimension ( 9 ) :: x17 = (/ &
    0.0000000000000000D+00, 7.4109534999454085D-01, 1.2304236340273060D+00, &
    1.7320508075688772D+00, 2.5960831150492023D+00, 2.8612795760570582D+00, &
    4.1849560176727323D+00, 5.1870160399136562D+00, 6.3633944943363696D+00 /)
  real ( kind = 8 ), dimension ( 9 ) :: w17 = (/ &
    2.6692223033505302D-01, 2.5456123204171222D-01, 1.4192654826449365D-02, &
    8.8681002152028010D-02, 1.9656770938777492D-03, 7.0334802378279075D-03, &
    1.0563783615416941D-04, -8.2049207541509217D-07, 2.1136499505424257D-08 /)
  real ( kind = 8 ), dimension ( 10 ) :: x19 = (/ &
    0.0000000000000000D+00, 7.4109534999454085D-01, 1.2304236340273060D+00, &
    1.7320508075688772D+00, 2.5960831150492023D+00, 2.8612795760570582D+00, &
    3.2053337944991944D+00, 4.1849560176727323D+00, 5.1870160399136562D+00, &
    6.3633944943363696D+00 /)
  real ( kind = 8 ), dimension ( 10 ) :: w19 = (/ &
    3.0346719985420623D-01, 2.0832499164960877D-01, 6.1151730125247716D-02, &
    6.4096054686807610D-02, 1.8085234254798462D-02, -6.3372247933737571D-03, &
    2.8848804365067559D-03, 6.0123369459847997D-05, 6.0948087314689840D-07, &
    8.6296846022298632D-10 /)
  real ( kind = 8 ), dimension ( 16 ) :: x31 = (/ &
    0.0000000000000000D+00, 2.4899229757996061D-01, 7.4109534999454085D-01, &
    1.2304236340273060D+00, 1.7320508075688772D+00, 2.2336260616769419D+00, &
    2.5960831150492023D+00, 2.8612795760570582D+00, 3.2053337944991944D+00, &
    3.6353185190372783D+00, 4.1849560176727323D+00, 5.1870160399136562D+00, &
    6.3633944943363696D+00, 7.1221067008046166D+00, 7.9807717985905606D+00, &
    9.0169397898903032D+00 /)
  real ( kind = 8 ), dimension ( 16 ) :: w31 = (/ &
    2.5890005324151566D-01, 2.8128101540033167D-02, 1.9968863511734550D-01, &
    6.5417392836092561D-02, 6.1718532565867179D-02, 1.7608475581318002D-03, &
    1.6592492698936010D-02, -5.5610063068358157D-03, 2.7298430467334002D-03, &
    1.5044205390914219D-05, 5.9474961163931621D-05, 6.1435843232617913D-07, &
    7.9298267864869338D-10, 5.1158053105504208D-12, -1.4840835740298868D-13, &
    1.2618464280815118D-15 /)
  real ( kind = 8 ), dimension ( 17 ) :: x33 = (/ &
    0.0000000000000000D+00, 2.4899229757996061D-01, 7.4109534999454085D-01, &
    1.2304236340273060D+00, 1.7320508075688772D+00, 2.2336260616769419D+00, &
    2.5960831150492023D+00, 2.8612795760570582D+00, 3.2053337944991944D+00, &
    3.6353185190372783D+00, 4.1849560176727323D+00, 5.1870160399136562D+00, &
    5.6981777684881099D+00, 6.3633944943363696D+00, 7.1221067008046166D+00, &
    7.9807717985905606D+00, 9.0169397898903032D+00 /)
  real ( kind = 8 ), dimension ( 17 ) :: w33 = (/ &
    1.3911022236338039D-01, 1.0387687125574284D-01, 1.7607598741571459D-01, &
    7.7443602746299481D-02, 5.4677556143463042D-02, 7.3530110204955076D-03, &
    1.1529247065398790D-02, -2.7712189007789243D-03, 2.1202259559596325D-03, &
    8.3236045295766745D-05, 5.5691158981081479D-05, 6.9086261179113738D-07, &
   -1.3486017348542930D-08, 1.5542195992782658D-09, -1.9341305000880955D-11, &
    2.6640625166231651D-13, -9.9313913286822465D-16 /)
  real ( kind = 8 ), dimension ( 18 ) :: x35 = (/ &
    0.0000000000000000D+00, 2.4899229757996061D-01, 7.4109534999454085D-01, &
    1.2304236340273060D+00, 1.7320508075688772D+00, 2.2336260616769419D+00, &
    2.5960831150492023D+00, 2.8612795760570582D+00, 3.2053337944991944D+00, &
    3.6353185190372783D+00, 4.1849560176727323D+00, 4.7364330859522967D+00, &
    5.1870160399136562D+00, 5.6981777684881099D+00, 6.3633944943363696D+00, &
    7.1221067008046166D+00, 7.9807717985905606D+00, 9.0169397898903032D+00 /)
  real ( kind = 8 ), dimension ( 18 ) :: w35 = (/ &
    5.1489450806921377D-04, 1.9176011588804434D-01, 1.4807083115521585D-01, &
    9.2364726716986353D-02, 4.5273685465150391D-02, 1.5673473751851151D-02, &
    3.1554462691875513D-03, 2.3113452403522071D-03, 8.1895392750226735D-04, &
    2.7524214116785131D-04, 3.5729348198975332D-05, 2.7342206801187888D-06, &
    2.4676421345798140D-07, 2.1394194479561062D-08, 4.6011760348655917D-10, &
    3.0972223576062995D-12, 5.4500412650638128D-15, 1.0541326582334014D-18 /)

  nhalf = ( n + 1 ) / 2

  if ( n == 1 ) then
    call r8vec_copy ( nhalf, x01, x )
    call r8vec_copy ( nhalf, w01, w )
  else if ( n == 3 ) then
    call r8vec_copy ( nhalf, x03, x )
    call r8vec_copy ( nhalf, w03, w )
  else if ( n == 7 ) then
    call r8vec_copy ( nhalf, x07, x )
    call r8vec_copy ( nhalf, w07, w )
  else if ( n == 9 ) then
    call r8vec_copy ( nhalf, x09, x )
    call r8vec_copy ( nhalf, w09, w )
  else if ( n == 17 ) then
    call r8vec_copy ( nhalf, x17, x )
    call r8vec_copy ( nhalf, w17, w )
  else if ( n == 19 ) then
    call r8vec_copy ( nhalf, x19, x )
    call r8vec_copy ( nhalf, w19, w )
  else if ( n == 31 ) then
    call r8vec_copy ( nhalf, x31, x )
    call r8vec_copy ( nhalf, w31, w )
  else if ( n == 33 ) then
    call r8vec_copy ( nhalf, x33, x )
    call r8vec_copy ( nhalf, w33, w )
  else if ( n == 35 ) then
    call r8vec_copy ( nhalf, x35, x )
    call r8vec_copy ( nhalf, w35, w )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KPN - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of N.'
    stop
  end if

  return
end
subroutine kpn_order ( l, n )

!*****************************************************************************80
!
!! KPN_ORDER computes the order of a KPN rule from the level.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, the level of the rule.  
!    1 <= L <= 25
!
!    Output, integer ( kind = 4 ) N, the order of the rule.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  if ( l < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KPN_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L <= 25 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop
  else if ( l == 1 ) then
    n = 1
  else if ( l <= 3 ) then
    n = 3
  else if ( l == 4 ) then
    n = 7
  else if ( l <= 8 ) then
    n = 9
  else if ( l == 9 ) then
    n = 17
  else if ( l <= 15 ) then
    n = 19
  else if ( l == 16 ) then
    n = 31
  else if ( l == 17 ) then
    n = 33
  else if ( l <= 25 ) then
    n = 35
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KPN_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L <= 25 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop
  end if

  return
end
subroutine kpu ( n, x, w )

!*****************************************************************************80
!
!! KPU provides data for Kronrod-Patterson quadrature with a uniform weight.
!
!  Discussion:
!
!    This data assumes integration over the interval [0,1] with 
!    weight function w(x) = 1.
!
!    For all orders N, the rule is formed by
!      X(1) with weight W(1),
!      X(2) and -X(2) with weight W(2),
!      X(3) and -X(3) with weight W(3) and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    Original MATLAB version by Florian Heiss, Viktor Winschel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309.
!
!    Thomas Patterson,
!    The optimal addition of points to quadrature formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) X((N+1)/2), the nodes.
!
!    Output, real ( kind = 8 ) W((N+1)/2), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) nhalf
  real ( kind = 8 ) w(*)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ), dimension ( 1 ) :: x01 = (/ &
    5.0000000000000000D-01 /)
  real ( kind = 8 ), dimension ( 1 ) :: w01 = (/ &
    1.0000000000000000D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: x03 = (/ &
    5.0000000000000000D-01, 8.8729829999999998D-01 /)
  real ( kind = 8 ), dimension ( 2 ) :: w03 = (/ &
    4.4444440000000002D-01, 2.7777780000000002D-01 /)
  real ( kind = 8 ), dimension ( 4 ) :: x07 = (/ &
    5.0000000000000000D-01, 7.1712189999999998D-01, 8.8729829999999998D-01, &
    9.8024560000000005D-01 /)
  real ( kind = 8 ), dimension ( 4 ) :: w07 = (/ &
    2.2545832254583223D-01, 2.0069872006987199D-01, 1.3424401342440134D-01, &
    5.2328105232810521D-02 /)
  real ( kind = 8 ), dimension ( 8 ) :: x15 = (/ &
    5.0000000000000000D-01, 6.1169330000000000D-01, 7.1712189999999998D-01, &
    8.1055149999999998D-01, 8.8729829999999998D-01, 9.4422960000000000D-01, &
    9.8024560000000005D-01, 9.9691600000000002D-01 /)
  real ( kind = 8 ), dimension ( 8 ) :: w15 = (/ &
    1.1275520000000000D-01, 1.0957840000000001D-01, 1.0031430000000000D-01, &
    8.5755999999999999D-02, 6.7207600000000006D-02, 4.6463600000000001D-02, &
    2.5801600000000001D-02, 8.5009000000000005D-03 /)
  real ( kind = 8 ), dimension ( 16 ) :: x31 = (/ &
    5.0000000000000000D-01, 5.5624450000000003D-01, 6.1169330000000000D-01, &
    6.6556769999999998D-01, 7.1712189999999998D-01, 7.6565989999999995D-01, &
    8.1055149999999998D-01, 8.5124809999999995D-01, 8.8729829999999998D-01, &
    9.1836300000000004D-01, 9.4422960000000000D-01, 9.6482740000000000D-01, &
    9.8024560000000005D-01, 9.9076560000000002D-01, 9.9691600000000002D-01, &
    9.9954909999999997D-01 /)
  real ( kind = 8 ), dimension ( 16 ) :: w31 = (/ &
    5.6377600000000014D-02, 5.5978400000000011D-02, 5.4789200000000017D-02, &
    5.2834900000000011D-02, 5.0157100000000017D-02, 4.6813600000000004D-02, &
    4.2878000000000006D-02, 3.8439800000000010D-02, 3.3603900000000006D-02, &
    2.8489800000000006D-02, 2.3231400000000003D-02, 1.7978600000000004D-02, &
    1.2903800000000003D-02, 8.2230000000000011D-03, 4.2173000000000011D-03, &
    1.2724000000000001D-03 /)
  real ( kind = 8 ), dimension ( 32 ) :: x63 = (/ &
    5.0000000000000000D-01, 5.2817210000000003D-01, 5.5624450000000003D-01, &
    5.8411769999999996D-01, 6.1169330000000000D-01, 6.3887490000000002D-01, &
    6.6556769999999998D-01, 6.9167970000000001D-01, 7.1712189999999998D-01, &
    7.4180900000000005D-01, 7.6565989999999995D-01, 7.8859789999999996D-01, &
    8.1055149999999998D-01, 8.3145480000000005D-01, 8.5124809999999995D-01, &
    8.6987800000000004D-01, 8.8729829999999998D-01, 9.0347029999999995D-01, &
    9.1836300000000004D-01, 9.3195399999999995D-01, 9.4422960000000000D-01, &
    9.5518559999999997D-01, 9.6482740000000000D-01, 9.7317140000000002D-01, &
    9.8024560000000005D-01, 9.8609139999999995D-01, 9.9076560000000002D-01, &
    9.9434239999999996D-01, 9.9691600000000002D-01, 9.9860309999999997D-01, &
    9.9954909999999997D-01, 9.9993650000000001D-01 /)
  real ( kind = 8 ), dimension ( 32 ) :: w63 = (/ &
    2.8188799999999993D-02, 2.8138799999999992D-02, 2.7989199999999992D-02, &
    2.7740699999999993D-02, 2.7394599999999995D-02, 2.6952699999999993D-02, &
    2.6417499999999993D-02, 2.5791599999999994D-02, 2.5078599999999993D-02, &
    2.4282199999999993D-02, 2.3406799999999995D-02, 2.2457299999999996D-02, &
    2.1438999999999996D-02, 2.0357799999999995D-02, 1.9219899999999998D-02, &
    1.8032199999999998D-02, 1.6801899999999998D-02, 1.5536799999999996D-02, &
    1.4244899999999996D-02, 1.2934799999999996D-02, 1.1615699999999998D-02, &
    1.0297099999999998D-02, 8.9892999999999987D-03, 7.7033999999999983D-03, &
    6.4518999999999983D-03, 5.2490999999999987D-03, 4.1114999999999988D-03, &
    3.0577999999999990D-03, 2.1087999999999997D-03, 1.2894999999999998D-03, &
    6.3259999999999987D-04, 1.8159999999999997D-04 /)

  nhalf = ( n + 1 ) / 2

  if ( n == 1 ) then
    call r8vec_copy ( nhalf, x01, x )
    call r8vec_copy ( nhalf, w01, w )
  else if ( n == 3 ) then
    call r8vec_copy ( nhalf, x03, x )
    call r8vec_copy ( nhalf, w03, w )
  else if ( n == 7 ) then
    call r8vec_copy ( nhalf, x07, x )
    call r8vec_copy ( nhalf, w07, w )
  else if ( n == 15 ) then
    call r8vec_copy ( nhalf, x15, x )
    call r8vec_copy ( nhalf, w15, w )
  else if ( n == 31 ) then
    call r8vec_copy ( nhalf, x31, x )
    call r8vec_copy ( nhalf, w31, w )
  else if ( n == 63 ) then
    call r8vec_copy ( nhalf, x63, x )
    call r8vec_copy ( nhalf, w63, w )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KPU - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of N.'
    stop
  end if

  return
end
subroutine kpu_order ( l, n )

!*****************************************************************************80
!
!! KPU_ORDER computes the order of a KPU rule from the level.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, the level of the rule.  
!    1 <= L <= 25
!
!    Output, integer ( kind = 4 ) N, the order of the rule.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  if ( l < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KPU_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L <= 25 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop
  else if ( l == 1 ) then
    n = 1
  else if ( l <= 3 ) then
    n = 3
  else if ( l <= 6 ) then
    n = 7
  else if ( l <= 12 ) then
    n = 15
  else if ( l <= 24 ) then
    n = 31
  else if ( l <= 25 ) then
    n = 63
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KPU_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L <= 25 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop
  end if

  return
end
subroutine num_seq ( n, k, seq_num )

!*****************************************************************************80
!
!! NUM_SEQ returns the number of compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The 28 compositions of 6 into three parts are:
!
!      6 0 0,  5 1 0,  5 0 1,  4 2 0,  4 1 1,  4 0 2,
!      3 3 0,  3 2 1,  3 1 2,  3 0 3,  2 4 0,  2 3 1,
!      2 2 2,  2 1 3,  2 0 4,  1 5 0,  1 4 1,  1 3 2,
!      1 2 3,  1 1 4,  1 0 5,  0 6 0,  0 5 1,  0 4 2,
!      0 3 3,  0 2 4,  0 1 5,  0 0 6.
!
!    The formula for the number of compositions of N into K parts is
!
!      Number = ( N + K - 1 )! / ( N! * ( K - 1 )! )
!
!    Describe the composition using N '1's and K-1 dividing lines '|'.
!    The number of distinct permutations of these symbols is the number
!    of compositions.  This is equal to the number of permutations of
!    N+K-1 things, with N identical of one kind and K-1 identical of another.
!
!    Thus, for the above example, we have:
!
!      Number = ( 6 + 3 - 1 )! / ( 6! * (3-1)! ) = 28
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer whose compositions are desired.
!
!    Input, integer ( kind = 4 ) K, the number of parts in the composition.
!
!    Output, integer ( kind = 4 ) SEQ_NUM, the number of compositions of N
!    into K parts.
!
  implicit none

  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seq_num

  seq_num = i4_choose ( n + k - 1, n )

  return
end
subroutine nwspgr_size ( rule_order, dim, k, r_size )

!*****************************************************************************80
!
!! NWSPGR_SIZE determines the size of a sparse grid rule.
!
!  Discussion:
!
!    This routine does a "raw" count, that is, it does not notice that many
!    points may be repeated, in which case, the size of the rule could be
!    reduced by merging repeated points and combining the corresponding weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
!
!    Input, external RULE_ORDER ( l, n ), the name of a subroutine which
!    is given the level L and returns the order N of the corresponding rule.
!
!    Input, integer ( kind = 4 ) DIM, the dimension of the integration problem.
!
!    Input, integer ( kind = 4 ) K, the level.  When using the built in 1D 
!    rules, the resulting sparse grid will be exact for polynomials up to total
!    order 2*K-1.  When using the built-in quadrature rules, the maximum value 
!    of K that is available is 25.
!
!    Output, integer ( kind = 4 ) R_SIZE, the "size" of the rule, that is,
!    the number of weights and multidimensional quadrature points that will
!    be needed.  The size of the rule will be reduced when duplicate points
!    are merged.
!
  implicit none

  integer ( kind = 4 ) dim

  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable :: is(:,:)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) level
  integer ( kind = 4 ) maxq
  integer ( kind = 4 ) minq
  integer ( kind = 4 ) n
  integer ( kind = 4 ), allocatable :: n1d(:)
  integer ( kind = 4 ) n1d_total
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r_size
  integer ( kind = 4 ), allocatable :: rq(:)
  external rule_order
  integer ( kind = 4 ) seq_num
!
!  Determine the size of each 1D rule.
!
  allocate ( n1d(k) )

  do level = 1, k

    call rule_order ( level, n )
    n1d(level) = n

  end do

  n1d_total = sum ( n1d(1:k) )
!
!  Go through the motions of generating the rules.
!
  minq = max ( 0, k - dim )
  maxq = k - 1
  r_size = 0

  do q = minq, maxq
!
!  Compute the D-dimensional vectors that sum to Q+DIM.
!
    call num_seq ( q, dim, seq_num )

    allocate ( is(1:seq_num,1:dim) )

    call get_seq ( dim, q + dim, seq_num, is )
!
!  Determine the size of each rule.
!
    allocate ( rq(1:seq_num) )

    do i = 1, seq_num
      rq(i) = product ( n1d(is(i,1:dim)) )
    end do
!
!  Add the sizes to the total.
!
    r_size = r_size + sum ( rq )

    deallocate ( is )
    deallocate ( rq )

  end do

  return
end
subroutine quad_rule_print ( m, n, x, w, title )

!*****************************************************************************80
!
!! QUAD_RULE_PRINT prints a multidimensional quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the abscissas.
!
!    Input, real ( kind = 8 ) W(N), the weights.
!
!    Input, character ( len = * ) TITLE, a title for the rule.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character ( len = * ) title
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do j = 1, n
    write ( *, '(2x,i2,2x,f10.6,a)', advance = 'no' ) j, w(j), ' * f ('
    do i = 1, m
      write ( *, '(f10.6)', advance = 'no' ) x(i,j)
      if ( i < m ) then
        write ( *, '('','')', advance = 'no' )
      else
        write ( *, '('' )'')', advance = 'yes' )
      end if
    end do
  end do

  return
end
subroutine r8cvv_offset ( m, nr, roff )

!*****************************************************************************80
!
!! R8CVV_OFFSET determines the row offsets of an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) NR(M), the row sizes.
!
!    Output, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) roff(m+1)
  integer ( kind = 4 ) nr(m)

  roff(1) = 0
  do i = 1, m
    roff(i+1) = roff(i) + nr(i)
  end do

  return
end
subroutine r8cvv_print ( mn, a, m, roff, title )

!*****************************************************************************80
!
!! R8CVV_PRINT prints an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  real ( kind = 8 ) a(mn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) roff(m+1)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, m

    k1 = roff(i) + 1
    k2 = roff(i+1)

    do klo = k1, k2, 5
      khi = min ( klo + 5 - 1, k2 )
      if ( klo == k1 ) then
        write ( *, '(i5,2x, 5g14.6)' ) i, a(klo:khi)
      else
        write ( *, '(5x,2x, 5g14.6)' )    a(klo:khi)
      end if
    end do

  end do

  return
end
subroutine r8cvv_rget ( mn, a, m, roff, i, ai )

!*****************************************************************************80
!
!! R8CVV_RGET gets row I from an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row.
!    1 <= I <= M.
!
!    Output, real ( kind = 8 ) AI(NR(I)), the value of A(I,*).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  real ( kind = 8 ) a(mn)
  real ( kind = 8 ) ai(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  ai(1:nv) = a(k1:k2)

  return
end
subroutine r8cvv_rset ( mn, a, m, roff, i, ai )

!*****************************************************************************80
!
!! R8CVV_RSET sets row I from an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row.
!    1 <= I <= M.
!
!    Input, real ( kind = 8 ) AI(NR(I)), the new value of A(I,*).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  real ( kind = 8 ) a(mn)
  real ( kind = 8 ) ai(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  a(k1:k2) = ai(1:nv)

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
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
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
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
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8vec_copy ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_COPY copies an R8VEC.
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
!    17 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), the vector to be copied.
!
!    Output, real ( kind = 8 ) A2(N), a copy of A1.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)

  a2(1:n) = a1(1:n)

  return
end
subroutine r8vec_direct_product ( factor_index, factor_order, factor_value, &
  factor_num, point_num, x )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out that task for the abscissas X.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      X(1:4) = ( 1, 2, 3, 4 )
!
!    Rule 2:
!      Order = 3
!      X(1:3) = ( 10, 20, 30 )
!
!    Rule 3:
!      Order = 2
!      X(1:2) = ( 100, 200 )
!
!    Product Rule:
!      Order = 24
!      X(1:24) =
!        ( 1, 10, 100 )
!        ( 2, 10, 100 )
!        ( 3, 10, 100 )
!        ( 4, 10, 100 )
!        ( 1, 20, 100 )
!        ( 2, 20, 100 )
!        ( 3, 20, 100 )
!        ( 4, 20, 100 )
!        ( 1, 30, 100 )
!        ( 2, 30, 100 )
!        ( 3, 30, 100 )
!        ( 4, 30, 100 )
!        ( 1, 10, 200 )
!        ( 2, 10, 200 )
!        ( 3, 10, 200 )
!        ( 4, 10, 200 )
!        ( 1, 20, 200 )
!        ( 2, 20, 200 )
!        ( 3, 20, 200 )
!        ( 4, 20, 200 )
!        ( 1, 30, 200 )
!        ( 2, 30, 200 )
!        ( 3, 30, 200 )
!        ( 4, 30, 200 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind = 8 ) X(FACTOR_NUM,POINT_NUM), the elements of
!    the direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer START, the first location of a block of values to set.
!
!    Local, integer CONTIG, the number of consecutive values to set.
!
!    Local, integer SKIP, the distance from the current value of START
!    to the next location of a block of values to set.
!
!    Local, integer REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  real ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real ( kind = 8 ) x(factor_num,point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    x(1:factor_num,1:point_num) = 0.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      x(factor_index,start:start+contig-1) = factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
subroutine r8vec_direct_product2 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, w )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out the task involving the weights W.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      W(1:4) = ( 2, 3, 5, 7 )
!
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind = 8 ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of values
!    to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values
!    to set.
!
!    Local, integer SKIP, the distance from the current value of START
!    to the next location of a block of values to set.
!
!    Local, integer REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  real ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real ( kind = 8 ) w(point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    w(1:point_num) = 1.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

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

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vecs_print ( m, nvec, na, a, title )

!*****************************************************************************80
!
!! R8VECS_PRINT prints a packed vector.
!
!  Example:
!
!    M = 5
!    NVEC = (/ 1, 4, 6, 11, 13, 14 /)
!    A = (/ 11, 12, 13, 21, 22, 31, 32, 33, 34, 35, 41, 42, 51 /)
!
!    11 12 13
!    21 22
!    31 32 33 34 35
!    41 42
!    51
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of vectors packed into A.
!
!    Input, integer ( kind = 4 ) NVEC(M+1), pointers to the first entry 
!    in each vector.
!
!    Input, integer ( kind = 4 ) NA, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(NA), the packed vector.  The I-th vector
!    extends from A(NVEC(I)) to A(NVEC(I+1)-1).
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) na

  real ( kind = 8 ) a(na)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nvec(m+1)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  khi = 0

  do i = 1, m

    n = nvec(i+1) - nvec(i)

    do jlo = 1, n, 5
      jhi = min ( jlo + 5 - 1, n )
      klo = khi + 1
      khi = klo + ( jhi - jlo )
      if ( jlo == 1 ) then
        write ( *, '(2x,i3,2x,5g14.6)' ) i, a(klo:khi)
      else
        write ( *, '(7x,5g14.6)' ) a(klo:khi)
      end if
    end do

  end do

  return
end
subroutine rules_1d_set ( k, sym, f_size, f_set, r1d_size, r1d, x1d, w1d )

!*****************************************************************************80
!
!! RULES_1D_SIZE determines the size of the packed RULES_1D array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the maximum level.
!    1 <= K.
!
!    Input, logical SYM, is TRUE if the rule is symmetric, so that symmetric
!    storage should be applied.
!
!    Input, external subroutine F_SIZE ( K, ORDER ), for any level K,
!    returns in ORDER the size of the corresponding 1D rule.
!
!    Input, external subroutine F_SET ( ORDER, X, W ), for any level K,
!    returns the ORDER values of nodes and weights.
!
!    Input, integer ( kind = 4 ) R1D_SIZE, the dimension required for the
!    X1D and W1D arrays.
!
!    Input, integer ( kind = 4 ) R1D(K+1), a pointer array.  
!    For 1 <= LEVEL <= K, the X values for rule I are stored in 
!    X1D(R1D(LEVEL)) through X1D(R1D(LEVEL+1)-1), and similarly for W.
!
!    Output, real ( kind = 8 ) X1D(R1D_SIZE), the nodes for levels 1 to K.
!
!    Output, real ( kind = 8 ) W1D(R1D_SIZE), the weights for levels 1 to K.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) r1d_size

  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) level
  integer ( kind = 4 ) order
  integer ( kind = 4 ) r1d(k+1)
  logical sym
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ) w1d(r1d_size)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) x1d(r1d_size)

  ihi = 0

  do level = 1, k

    call f_size ( level, order )

    allocate ( x(1:order) )
    allocate ( w(1:order) )

    call f_set ( order, x, w )

    if ( sym ) then
      jlo = ( order + 2 ) / 2
      jhi = order
      order = ( order + 1 ) / 2
    else
      jlo = 1
      jhi = order
    end if

    jhi = jlo + order - 1

    ilo = ihi + 1
    ihi = ilo + order - 1

    x1d(ilo:ihi) = x(jlo:jhi)
    w1d(ilo:ihi) = w(jlo:jhi)

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine rules_1d_size ( k, sym, f_size, r1d_size, r1d )

!*****************************************************************************80
!
!! RULES_1D_SIZE determines the size of the packed RULES_1D array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the maximum level.
!    1 <= K.
!
!    Input, logical SYM, is TRUE if the rule is symmetric, so that symmetric
!    storage should be applied.
!
!    Input, external subroutine F_SIZE ( K, ORDER ), for any level K,
!    returns in ORDER the size of the corresponding 1D rule.
!
!    Output, integer ( kind = 4 ) R1D_SIZE, the dimension required for the
!    X1D and W1D arrays.
!
!    Output, integer ( kind = 4 ) R1D(K+1), a pointer array.  
!    For 1 <= LEVEL <= K, the X values for rule I are stored in 
!    X1D(R1D(LEVEL)) through X1D(R1D(LEVEL+1)-1), and similarly for W.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) level
  integer ( kind = 4 ) order
  integer ( kind = 4 ) r1d(k+1)
  integer ( kind = 4 ) r1d_size
  logical sym

  r1d_size = 0
  r1d(1) = 1

  do level = 1, k

    call f_size ( level, order )

    if ( sym ) then
      order = ( order + 1 ) / 2
    end if

    r1d_size = r1d_size + order
    r1d(level+1) = r1d(level) + order

  end do

  return
end
subroutine symmetric_sparse_size ( nr, dim, nodes, x0, nr2 )

!*****************************************************************************80
!
!! SYMMETRIC_SPARSE_SIZE sizes a symmetric sparse rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
! 
!    Input, integer ( kind = 4 ) DIM, the dimension.
!
!    Input, integer ( kind = 4 ) NR, the dimension of the rule in the 
!    positive orthant.
!
!    Input, real ( kind = 8 ) NODES(NR,DIM), the nodes for the positive orthant.
!
!    Input, real ( kind = 8 ) X0, the point of symmetry for the 1D rule, 
!    typically 0.
!
!    Output, integer NR2, the dimension of the rule when "unfolded" to the 
!    full space.
!
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) nr

  integer ( kind = 4 ) count
  integer ( kind = 4 ) j
  real ( kind = 8 ) nodes(nr,dim)
  integer ( kind = 4 ) nr2
  integer ( kind = 4 ) r
  real ( kind = 8 ) x0
!
!  Count the size of the full rule.
!
    nr2 = 0

    do r = 1, nr
      count = 1
      do j = 1, dim
        if ( nodes(r,j) /= x0 ) then
          count = 2 * count
        end if
      end do
      nr2 = nr2 + count
    end do

  return
end
subroutine tensor_product ( d, order1d, n1d, x1d, w1d, n, xnd, wnd )

!*****************************************************************************80
!
!! TENSOR_PRODUCT generates a tensor product quadrature rule.
!
!  Discussion:
!
!    The Kronecker product of an K by L matrix A and an M by N matrix B
!    is the K*M by L*N matrix formed by
!
!      a(1,1) * B,  a(1,2) * B,  ..., a(1,l) * B
!      a(2,1) * B,  a(2,2) * B,  ..., a(2,l) * B
!      ..........   ..........   .... ..........
!      a(k,1) * B,  a(k,2) * B,  ..., a(k,l) * B
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2012
!
!  Author:
!
!    Original MATLAB version by Florian Heiss, Viktor Winschel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) D, the spatial dimension.
!
!    Input, integer ( kind = 4 ) ORDER1D(D), the order of each 1D rule.
!
!    Input, integer ( kind = 4 ) N1D, the number of 1D items.
!
!    Input, real ( kind = 8 ) X1D(N1D), the 1D nodes.
!
!    Input, real ( kind = 8 ) W1D(N1D), the 1D weights.
!
!    Input, integer ( kind = 4 ) N, the number of N-dimensional items.
!
!    Output, real ( kind = 8 ) XND(D,N), the nodes.
!
!    Output, real ( kind = 8 ) WND(N), the weights.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n1d

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) order1d(d)
  real ( kind = 8 ) w1d(n1d)
  real ( kind = 8 ) wnd(n)
  real ( kind = 8 ) x1d(n1d)
  real ( kind = 8 ) xnd(d,n)
!
!  Compute the weights.
!
  i2 = 0
  do i = 1, d
    i1 = i2 + 1
    i2 = i2 + order1d(i)
    call r8vec_direct_product2 ( i, order1d(i), w1d(i1:i2), d, n, wnd )
  end do
!
!  Compute the points.
!
  i2 = 0
  do i = 1, d
    i1 = i2 + 1
    i2 = i2 + order1d(i)
    call r8vec_direct_product ( i, order1d(i), x1d(i1:i2), d, n, xnd )
  end do

  return
end
subroutine tensor_product_cell ( nc, xc, wc, dim, nr, roff, np, xp, wp )

!*****************************************************************************80
!
!! TENSOR_PRODUCT_CELL generates a tensor product quadrature rule.
!
!  Discussion:
!
!    The Kronecker product of an K by L matrix A and an M by N matrix B
!    is the K*M by L*N matrix formed by
!
!      a(1,1) * B,  a(1,2) * B,  ..., a(1,l) * B
!      a(2,1) * B,  a(2,2) * B,  ..., a(2,l) * B
!      ..........   ..........   .... ..........
!      a(k,1) * B,  a(k,2) * B,  ..., a(k,l) * B
!
!    The 1D factors are stored in a kind of cell array structure.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NC, the number of items in the cell arrays.
!
!    Input, real ( kind = 8 ) XC(NC), a cell array containing points for 
!    1D rules.
!
!    Input, real ( kind = 8 ) WC(NC), a cell array containing weights for
!    1D rules.
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NR(DIM), the length of each row of the 
!    cell array.
!
!    Input, integer ( kind = 4 ) ROFF(DIM+1), offsets for the cell arrays.
!
!    Input, integer ( kind = 4 ) NP, the number of points in the product rule.
!
!    Output, real ( kind = 8 ) XP(DIM,NP), the nodes.
!
!    Output, real ( kind = 8 ) WP(NP), the weights.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) np

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n1d
  integer ( kind = 4 ) nr(dim)
  integer ( kind = 4 ) roff(dim+1)
  real ( kind = 8 ), allocatable :: w1d(:)
  real ( kind = 8 ) wc(nc)
  real ( kind = 8 ) wp(np)
  real ( kind = 8 ), allocatable :: x1d(:)
  real ( kind = 8 ) xc(nc)
  real ( kind = 8 ) xp(dim,np)
!
!  Compute the weights.
!
  do i = 1, dim
    n1d = nr(i)
    allocate ( w1d(1:n1d) )
    call r8cvv_rget ( nc, wc, dim, roff, i, w1d )
    call r8vec_direct_product2 ( i, n1d, w1d, dim, np, wp )
    deallocate ( w1d )
  end do
!
!  Compute the points.
!
  do i = 1, dim
    n1d = nr(i)
    allocate ( x1d(1:n1d) )
    call r8cvv_rget ( nc, xc, dim, roff, i, x1d )
    call r8vec_direct_product ( i, n1d, x1d, dim, np, xp )
    deallocate ( x1d )
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

  character ( len = 8 ) ampm
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
