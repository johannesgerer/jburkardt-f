program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA205_PRB.
!
!  Discussion:
!
!    ASA205_PRB tests ASA205.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA205_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA205 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA205_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 examines a simple case with no repeated sum values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  external eval01
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ), dimension ( n ) :: colsum = (/ 2, 3, 1, 4 /)
  integer ( kind = 4 ), dimension ( m ) :: rowsum = (/ 5, 3, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  The tables will not have any multiplicities.'

  call i4vec_print ( m, rowsum, '  The row sums:' )

  call i4vec_print ( n, colsum, '  The column sums:' )

  call enum ( m, n, rowsum, colsum, eval01, ifault )

  if ( ifault /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ENUM returned error flag IFAULT = ', ifault
  end if

  return
end
subroutine eval01 ( iflag, table, m, n, rowsum, colsum, prob, mult )

!*****************************************************************************80
!
!! EVAL01 is called by ENUM every time a new contingency table is determined.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IFLAG, input flag.
!    1, this is the first call.  No table is input.
!    2, this is a call with a new table.
!    3, this is the last call.  No table is input.
!
!    Input, integer ( kind = 4 ) TABLE(M,N), the current contingency table.
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, integer ( kind = 4 ) ROWSUM(M), the row sums.
!
!    Input, integer ( kind = 4 ) COLSUM(N), the column sums.
!
!    Input, real PROB, the logarithm of the hypergeometric probability
!    of this table.
!
!    Input, integer ( kind = 4 ) MULT, the multiplicity of this table, that is,
!    the number of different tables that still have the same set of
!    entries, but differ by a permutation of some rows and columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) colsum(n)
  integer ( kind = 4 ), save :: count1 = 0
  integer ( kind = 4 ), save :: count2 = 0
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) mult
  real ( kind = 8 ) prob
  real ( kind = 8 ), save :: psum = 0.0D+00
  integer ( kind = 4 ) rowsum(m)
  integer ( kind = 4 ) table(m,n)
!
!  First call, no table, initialize.
!
  if ( iflag == 1 ) then

    count1 = 0
    count2 = 0
    psum = 0.0D+00

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVAL01'
    write ( *, '(a)' ) '  Only first ten cases will be printed.'
    write ( *, '(a)' ) ' '
!
!  Call with a new table.
!
  else if ( iflag == 2 ) then

    count1 = count1 + 1
    count2 = count2 + mult
    psum = psum + mult * exp ( prob )

    if ( count1 <= 10 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EVAL01:'
      write ( *, '(i3,i3,g14.6)' ) count1, mult, prob

      call i4mat_print ( m, n, table, '  Table' )

    end if
!
!  Last call, no table.
!
  else if ( iflag == 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVAL01 summary'
    write ( *, '(a,i8)' ) '  Number of cases (ignoring multiplicity):', count1
    write ( *, '(a,i8)' ) '  Number of cases (allowing multiplicity):', count2
    write ( *, '(a,g14.6)' ) '  Probability sum = ', psum

  end if

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 examines a case where a sum value is repeated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n = 3

  external eval02
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ), dimension ( n ) :: colsum = (/ 1, 2, 1 /)
  integer ( kind = 4 ), dimension ( m ) :: rowsum = (/ 3, 1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  The data will have multiplicities.'

  call i4vec_print ( m, rowsum, '  The row sums:' )

  call i4vec_print ( n, colsum, '  The column sums:' )

  call enum ( m, n, rowsum, colsum, eval02, ifault )

  if ( ifault /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ENUM returned error flag IFAULT = ', ifault
  end if

  return
end
subroutine eval02 ( iflag, table, m, n, rowsum, colsum, prob, mult )

!*****************************************************************************80
!
!! EVAL02 is called by ENUM every time a new contingency table is determined.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IFLAG, input flag.
!    1, this is the first call.  No table is input.
!    2, this is a call with a new table.
!    3, this is the last call.  No table is input.
!
!    Input, integer ( kind = 4 ) TABLE(M,N), the current contingency table.
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, integer ( kind = 4 ) ROWSUM(M), the row sums.
!
!    Input, integer ( kind = 4 ) COLSUM(N), the column sums.
!
!    Input, real PROB, the logarithm of the hypergeometric probability
!    of this table.
!
!    Input, integer ( kind = 4 ) MULT, the multiplicity of this table, that is,
!    the number of different tables that still have the same set of
!    entries, but differ by a permutation of some rows and columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) colsum(n)
  integer ( kind = 4 ), save :: count1 = 0
  integer ( kind = 4 ), save :: count2 = 0
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) mult
  real ( kind = 8 ) prob
  real ( kind = 8 ), save :: psum = 0.0D+00
  integer ( kind = 4 ) rowsum(m)
  integer ( kind = 4 ) table(m,n)
!
!  First call, no table, initialize.
!
  if ( iflag == 1 ) then

    count1 = 0
    count2 = 0
    psum = 0.0D+00
!
!  Call with a new table.
!
  else if ( iflag == 2 ) then

    count1 = count1 + 1
    count2 = count2 + mult
    psum = psum + mult * exp ( prob )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVAL02:'
    write ( *, '(i3,i3,g14.6)' ) count1, mult, prob

    call i4mat_print ( m, n, table, '  Table' )
!
!  Last call, no table.
!
  else if ( iflag == 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVAL02 summary'
    write ( *, '(a,i8)' ) '  Number of cases (ignoring multiplicity):', count1
    write ( *, '(a,i8)' ) '  Number of cases (allowing multiplicity):', count2
    write ( *, '(a,g14.6)' ) '  Probability sum = ', psum

  end if

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 examines a test case from the paper, with 118489 tables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 3

  external eval03
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ), dimension ( n ) :: colsum = (/ 4, 57, 59 /)
  integer ( kind = 4 ), dimension ( m ) :: rowsum = (/ 3, 38, 39, 40 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Big problem test from the paper.'

  call i4vec_print ( m, rowsum, '  The row sums:' )

  call i4vec_print ( n, colsum, '  The column sums:' )

  call enum ( m, n, rowsum, colsum, eval03, ifault )

  if ( ifault /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ENUM returned error flag IFAULT = ', ifault
  end if

  return
end
subroutine eval03 ( iflag, table, m, n, rowsum, colsum, prob, mult )

!*****************************************************************************80
!
!! EVAL03 is called by ENUM every time a new contingency table is determined.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IFLAG, input flag.
!    1, this is the first call.  No table is input.
!    2, this is a call with a new table.
!    3, this is the last call.  No table is input.
!
!    Input, integer ( kind = 4 ) TABLE(M,N), the current contingency table.
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, integer ( kind = 4 ) ROWSUM(M), the row sums.
!
!    Input, integer ( kind = 4 ) COLSUM(N), the column sums.
!
!    Input, real PROB, the logarithm of the hypergeometric probability
!    of this table.
!
!    Input, integer ( kind = 4 ) MULT, the multiplicity of this table, that is,
!    the number of different tables that still have the same set of
!    entries, but differ by a permutation of some rows and columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) colsum(n)
  integer ( kind = 4 ), save :: count1 = 0
  integer ( kind = 4 ), save :: count2 = 0
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) mult
  real ( kind = 8 ) prob
  real ( kind = 8 ), save :: psum = 0.0D+00
  integer ( kind = 4 ) rowsum(m)
  integer ( kind = 4 ) table(m,n)
!
!  First call, no table, initialize.
!
  if ( iflag == 1 ) then

    count1 = 0
    count2 = 0
    psum = 0.0D+00
!
!  Call with a new table.
!
  else if ( iflag == 2 ) then

    count1 = count1 + 1
    count2 = count2 + mult
    psum = psum + mult * exp ( prob )
!
!  Last call, no table.
!
  else if ( iflag == 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVAL03 summary'
    write ( *, '(a,i8)' ) '  Number of cases (ignoring multiplicity): ', count1
    write ( *, '(a,i8)' ) '  Number of cases (allowing multiplicity): ', count2
    write ( *, '(a,g14.6)' ) '  Probability sum = ', psum
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Result from paper:'
    write ( *, '(a,i8)' ) '  Number of cases (ignoring multiplicity): ', 118489

  end if

  return
end
