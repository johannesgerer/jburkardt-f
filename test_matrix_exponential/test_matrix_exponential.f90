subroutine mexp_a ( test, n, a )

!*****************************************************************************80
!
!! MEXP_A returns the matrix for a given test.
!
!  Discussion:
!
!     1) Diagonal example
!     2) Symmetric example
!     3) Laub
!     4) Moler and Van Loan
!     5) Moler and Van Loan
!     6) Moler and Van Loan
!     7) Moler and Van Loan
!     8) Wikipedia example
!     9) NAG F01ECF
!    10) Ward #1
!    11) Ward #2
!    12) Ward #3
!    13) Ward #4
!    14) Moler example
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Laub,
!    Review of "Linear System Theory" by Joao Hespanha,
!    SIAM Review,
!    Volume 52, Number 4, December 2010, page 779-781.
!
!    Cleve Moler, Charles VanLoan,
!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!    Twenty-Five Years Later,
!    SIAM Review,
!    Volume 45, Number 1, March 2003, pages 3-49.
!
!    Cleve Moler,
!    Cleve's Corner: A Balancing Act for the Matrix Exponential,
!    July 23rd, 2012.
!
!    Robert Ward,
!    Numerical computation of the matrix exponential with accuracy estimate,
!    SIAM Journal on Numerical Analysis,
!    Volume 14, Number 4, September 1977, pages 600-610.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test case.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    a = reshape ( (/ &
      1.0D+00, 0.0D+00, &
      0.0D+00, 2.0D+00 /), (/ 2, 2 /) )
  else if ( test == 2 ) then
    a = reshape ( (/ &
      1.0D+00, 3.0D+00, &
      3.0D+00, 2.0D+00 /), (/ 2, 2 /) )
  else if ( test == 3 ) then
    a = reshape ( (/ &
      0.0D+00, -39.0D+00, &
      1.0D+00, -40.0D+00 /), (/ 2, 2 /) )
  else if ( test == 4 ) then
    a = reshape ( (/ &
      -49.0D+00, -64.0D+00, &
       24.0D+00,  31.0D+00 /), (/ 2, 2 /) )
  else if ( test == 5 ) then
    a = reshape ( (/ &
      0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
      6.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
      0.0D+00, 6.0D+00, 0.0D+00, 0.0D+00, &
      0.0D+00, 0.0D+00, 6.0D+00, 0.0D+00 /), (/ 4, 4 /) )
  else if ( test == 6 ) then
    a = reshape ( (/ &
      1.0D+00, 0.0D+00, &
      1.0D+00, 1.0D+00 /), (/ 2, 2 /) )
  else if ( test == 7 ) then
    a = reshape ( (/ &
      1.0D+00 + epsilon ( 1.0D+00 ), 0.0D+00, &
      1.0D+00,                       1.0D+00 - epsilon ( 1.0D+00 ) /), &
      (/ 2, 2 /) )
  else if ( test == 8 ) then
    a = reshape ( (/ &
      21.0D+00,  -5.0D+00,   4.0D+00, &
      17.0D+00,  -1.0D+00,   4.0D+00, &
       6.0D+00,  -6.0D+00,  16.0D+00 /), (/ 3, 3 /) )
  else if ( test == 9 ) then
    a = reshape ( (/ &
      1.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, &
      2.0D+00, 1.0D+00, 2.0D+00, 3.0D+00, &
      2.0D+00, 1.0D+00, 1.0D+00, 3.0D+00, &
      2.0D+00, 2.0D+00, 2.0D+00, 1.0D+00 /), (/ 4, 4 /) )
  elseif ( test == 10 ) then
    a = reshape ( (/ &
      4.0D+00, 1.0D+00, 1.0D+00, &
      2.0D+00, 4.0D+00, 1.0D+00, &
      0.0D+00, 1.0D+00, 4.0D+00 /), (/ 3, 3 /) )
  elseif ( test == 11 ) then
    a = reshape ( (/ &
      29.87942128909879D+00, &
       0.7815750847907159D+00, &
      -2.289519314033932D+00, &
       0.7815750847907159D+00, & 
      25.72656945571064D+00, &
       8.680737820540137D+00, &
      -2.289519314033932D+00, &
       8.680737820540137D+00, &
      34.39400925519054D+00 /), (/ 3, 3 /) )
  elseif ( test == 12 ) then
    a = reshape ( (/ &
      -131.0D+00, -390.0D+00, -387.0D+00, &
        19.0D+00,   56.0D+00,   57.0D+00, &
        18.0D+00,   54.0D+00,   52.0D+00 /), (/ 3, 3 /) )
  elseif ( test == 13 ) then
    a(1:n,1:n) = 0.0D+00
    do i = 1, n - 1
      a(i,i+1) = 1.0D+00
    end do
    a(n,1) = 1.0D-10
  elseif ( test == 14 ) then
    a(1,1) = 0.0D+00
    a(1,2) = 1.0D-08
    a(1,3) = 0.0D+00
    a(2,1) = - 2.0D+10 - 2.0D+08 / 3.0D+00
    a(2,2) = - 3.0D+00
    a(2,3) = 2.0D+10
    a(3,1) = 200.0D+00 / 3.0D+00
    a(3,2) = 0.0D+00
    a(3,3) = - 200.0D+00 / 3.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MEXP_A - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of TEST = ', test
  end if

  return
end
subroutine mexp_expa ( test, n, expa )

!*****************************************************************************80
!
!! MEXP_EXPA returns the "exact" exponential matrix for a given test.
!
!  Discussion:
!
!    In some cases, the "exact" value is given to six significant digits.
!
!     1) Diagonal example
!     2) Symmetric example
!     3) Laub
!     4) Moler and Van Loan
!     5) Moler and Van Loan
!     6) Moler and Van Loan
!     7) Moler and Van Loan
!     8) Wikipedia example
!     9) NAG F01ECF
!    10) Ward #1
!    11) Ward #2
!    12) Ward #3
!    13) Ward #4
!    14) Moler example
!
!    Thanks to Alex Griffing for correcting the value of matrix 3,
!    17 October 2012.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Laub,
!    Review of "Linear System Theory" by Joao Hespanha,
!    SIAM Review,
!    Volume 52, Number 4, December 2010, page 779-781.
!
!    Cleve Moler, Charles VanLoan,
!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!    Twenty-Five Years Later,
!    SIAM Review,
!    Volume 45, Number 1, March 2003, pages 3-49.
!
!    Cleve Moler,
!    Cleve's Corner: A Balancing Act for the Matrix Exponential,
!    July 23rd, 2012.
!
!    Robert Ward,
!    Numerical computation of the matrix exponential with accuracy estimate,
!    SIAM Journal on Numerical Analysis,
!    Volume 14, Number 4, September 1977, pages 600-610.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test case.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) EXPA(N,N), the exponential of the test matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) exp16
  real ( kind = 8 ) exp4
  real ( kind = 8 ) expa(n,n)
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    expa = reshape ( (/ &
      2.718281828459046D+00, 0.0D+00, &
      0.0D+00,               7.389056098930650D+00 /), (/ 2, 2 /) )
  else if ( test == 2 ) then
    expa = reshape ( (/ &
      39.322809708033859D+00,  46.166301438885753D+00, &
      46.166301438885768D+00,  54.711576854329110D+00 /), (/ 2, 2 /) )
  else if ( test == 3 ) then
    expa = reshape ( (/ &
       0.37756048D+00,  0.00968104D+00,&
      -0.37756048D+00, -0.00968104D+00 /), (/ 2, 2 /) )
  else if ( test == 4 ) then
    expa = reshape ( (/ &
      -0.735759D+00, -1.471518D+00, &
       0.551819D+00,  1.103638D+00 /), (/ 2, 2 /) )
  else if ( test == 5 ) then
    expa = reshape ( (/ &
      1.0D+00,  0.0D+00, 0.0D+00, 0.0D+00, &
      6.0D+00,  1.0D+00, 0.0D+00, 0.0D+00, &
     18.0D+00,  6.0D+00, 1.0D+00, 0.0D+00, &
     36.0D+00, 18.0D+00, 6.0D+00, 1.0D+00 /), (/ 4, 4 /) )
  else if ( test == 6 ) then
    expa = reshape ( (/ &
      2.718281828459046D+00, 0.0D+00, &
      2.718281828459046D+00, 2.718281828459046D+00 /), (/ 2, 2 /) )
  else if ( test == 7 ) then
    expa = reshape ( (/ &
      2.718309D+00, 0.0D+00, &
      2.718282D+00, 2.718255D+00 /), (/ 2, 2 /) )
  else if ( test == 8 ) then
    exp16 = exp ( 16.0D+00 )
    exp4 = exp ( 4.0D+00 )
    expa = 0.25D+00 * reshape ( (/ &
      13.0D+00 * exp16 -           exp4, &
      -9.0D+00 * exp16           + exp4, &
      16.0D+00 * exp16, &
      13.0D+00 * exp16 - 5.0D+00 * exp4, &
      -9.0D+00 * exp16 + 5.0D+00 * exp4, &
      16.0D+00 * exp16, &
       2.0D+00 * exp16 - 2.0D+00 * exp4, &
      -2.0D+00 * exp16 + 2.0D+00 * exp4, &
       4.0D+00 * exp16 /), (/ 3, 3 /) )
  else if ( test == 9 ) then
    expa = reshape ( (/ &
      740.7038D+00, 731.2510D+00, 823.7630D+00, 998.4355D+00, &
      610.8500D+00, 603.5524D+00, 679.4257D+00, 823.7630D+00, &
      542.2743D+00, 535.0884D+00, 603.5524D+00, 731.2510D+00, &
      549.1753D+00, 542.2743D+00, 610.8500D+00, 740.7038D+00 /), (/ 4, 4 /) )
  else if ( test == 10 ) then
    expa = reshape ( (/ &
      147.8666224463699D+00, &
      127.7810855231823D+00, &
      127.7810855231824D+00, &
      183.7651386463682D+00, &
      183.7651386463682D+00, &
      163.6796017231806D+00, &
      71.79703239999647D+00, &
      91.88256932318415D+00, &
     111.9681062463718D+00 /), (/ 3, 3 /) )
  else if ( test == 11 ) then
    expa = reshape ( (/ &
     5.496313853692378D+15, &
    -1.823188097200899D+16, &
    -3.047577080858001D+16, &
    -1.823188097200898D+16, &
     6.060522870222108D+16, &
     1.012918429302482D+17, &
    -3.047577080858001D+16, &
     1.012918429302482D+17, &
     1.692944112408493D+17 /), (/ 3, 3 /) )
  else if ( test == 12 ) then
    expa = reshape ( (/ &
    -1.509644158793135D+00, &
    -5.632570799891469D+00, &
    -4.934938326088363D+00, &
     0.3678794391096522D+00, &
     1.471517758499875D+00, &
     1.103638317328798D+00, &
     0.1353352811751005D+00, &
     0.4060058435250609D+00, &
     0.5413411267617766D+00 /), (/ 3, 3 /) )
  else if ( test == 13 ) then
    expa(1:n,1:n) = 0.0D+00
  else if ( test == 14 ) then
    expa = reshape ( (/ &
    4.468494682831735D-01, &
   -5.743067779479621D+06, &
    4.477229778494929D-01, &
    1.540441573839520D-09, &
   -1.528300386868247D-02, &
    1.542704845195912D-09, &
    4.628114535587735D-01, &
   -4.526542712784168D+06, &
    4.634806488376499D-01 /), (/ 3, 3 /) )    
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MEXP_EXPA - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of TEST = ', test
  end if

  return
end
subroutine mexp_n ( test, n )

!*****************************************************************************80
!
!! MEXP_N returns the matrix order for a given test.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test case.
!
!    Output, integer ( kind = 4 ) N, the order of the matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    n = 2
  else if ( test == 2 ) then
    n = 2
  else if ( test == 3 ) then
    n = 2
  else if ( test == 4 ) then
    n = 2
  else if ( test == 5 ) then
    n = 4
  else if ( test == 6 ) then
    n = 2
  else if ( test == 7 ) then
    n = 2
  else if ( test == 8 ) then
    n = 3
  else if ( test == 9 ) then
    n = 4
  else if ( test == 10 ) then
    n = 3
  else if ( test == 11 ) then
    n = 3
  else if ( test == 12 ) then
    n = 3
  else if ( test == 13 ) then
    n = 10
  else if ( test == 14 ) then
    n = 3
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MEXP_N - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of TEST = ', test
  end if

  return
end
subroutine mexp_story ( test )

!*****************************************************************************80
!
!! MEXP_STORY prints explanatory text for each problem.
!
!  Discussion:
!
!     1) Diagonal example
!     2) Symmetric example
!     3) Laub
!     4) Moler and Van Loan
!     5) Moler and Van Loan
!     6) Moler and Van Loan
!     7) Moler and Van Loan
!     8) Wikipedia example
!     9) NAG F01ECF
!    10) Ward #1
!    11) Ward #2
!    12) Ward #3
!    13) Ward #4
!    14) Moler example
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Laub,
!    Review of "Linear System Theory" by Joao Hespanha,
!    SIAM Review,
!    Volume 52, Number 4, December 2010, page 779-781.
!
!    Cleve Moler, Charles VanLoan,
!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!    Twenty-Five Years Later,
!    SIAM Review,
!    Volume 45, Number 1, March 2003, pages 3-49.
!
!    Cleve Moler,
!    Cleve's Corner: A Balancing Act for the Matrix Exponential,
!    July 23rd, 2012.
!
!    Robert Ward,
!    Numerical computation of the matrix exponential with accuracy estimate,
!    SIAM Journal on Numerical Analysis,
!    Volume 14, Number 4, September 1977, pages 600-610.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test case.
!
  implicit none

  integer ( kind = 4 ) test

  if ( test == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This matrix is diagonal.'
    write ( *, '(a)' ) '  The calculation of the matrix exponential is simple.'
  else if ( test == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This matrix is symmetric.'
    write ( *, '(a)' ) &
     '  The calculation of the matrix exponential is straightforward.'
  else if ( test == 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This example is due to Laub.'
    write ( *, '(a)' ) &
      '  This matrix is ill-suited for the Taylor series approach.'
    write ( *, '(a)' ) &
      '  As powers of A are computed, the entries blow up too quickly.'
  else if ( test == 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This example is due to Moler and Van Loan.'
    write ( *, '(a)' ) &
    '  The example will cause problems for the series summation approach,'
    write ( *, '(a)' ) '  as well as for diagonal Pade approximations.'
  else if ( test == 5 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This example is due to Moler and Van Loan.'
    write ( *, '(a)' ) '  This matrix is strictly upper triangular'
    write ( *, '(a)' ) '  All powers of A are zero beyond some (low) limit.' 
    write ( *, '(a)' ) &
      '  This example will cause problems for Pade approximations.'
  else if ( test == 6 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This example is due to Moler and Van Loan.'
    write ( *, '(a)' ) &
      '  This matrix does not have a complete set of eigenvectors.'
    write ( *, '(a)' ) '  That means the eigenvector approach will fail.'
  else if ( test == 7 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This example is due to Moler and Van Loan.'
    write ( *, '(a)' ) '  This matrix is very close to example 5.'
    write ( *, '(a)' ) &
      '  Mathematically, it has a complete set of eigenvectors.'
    write ( *, '(a)' ) &
      '  Numerically, however, the calculation will be suspect.'
  else if ( test == 8 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This matrix was an example in Wikipedia.'
  else if ( test == 9 ) then
    write ( *, '(a)' ) ' ' 
    write ( *, '(a)' ) '  This matrix is due to the NAG Library.'
    write ( *, '(a)' ) '  It is an example for function F01ECF.'
  else if ( test == 10 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This is Ward''s example #1.'
    write ( *, '(a)' ) '  It is defective and nonderogatory.'
    write ( *, '(a)' ) '  The eigenvalues are 3, 3 and 6.'
  else if ( test == 11 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This is Ward''s example #2.'
    write ( *, '(a)' ) '  It is a symmetric matrix.'
    write ( *, '(a)' ) '  The eigenvalues are 20, 30, 40.'
  else if ( test == 12 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This is Ward''s example #3.'
    write ( *, '(a)' ) &
      '  Ward''s algorithm has difficulty estimating the accuracy'
    write ( *, '(a)' ) '  of its results.  The eigenvalues are -1, -2, -20.'
  else if ( test == 13 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This is Ward''s example #4.'
    write ( *, '(a)' ) '  This is a version of the Forsythe matrix.'
    write ( *, '(a)' ) '  The eigenvector problem is badly conditioned.'
    write ( *, '(a)' ) &
      '  Ward''s algorithm has difficulty estimating the accuracy'
    write ( *, '(a)' ) '  of its results for this problem.'
  else if ( test == 14 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This is Moler''s example.'
    write ( *, '(a)' ) &
      '  This badly scaled matrix caused problems for MATLAB''s expm().'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MEXP_STORY - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of TEST = ', test
  end if

  return
end
subroutine mexp_test_num ( test_num )

!*****************************************************************************80
!
!! MEXP_TEST_NUM returns the number of matrix exponential tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) TEST_NUM, the number of tests.
!
  implicit none

  integer ( kind = 4 ) test_num

  test_num = 14

  return
end
