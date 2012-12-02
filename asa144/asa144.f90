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
!    04 November 2003
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

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
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

      write ( *, '(i5,1x,10a8)' ) i, ( ctemp(j), j = 1, inc )

    end do

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
!    An I4VEC is a vector of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
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
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i12)' ) i, a(i)
  end do

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
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
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine rcont ( nrow, ncol, nrowt, ncolt, nsubt, matrix, key, ifault )

!*****************************************************************************80
!
!! RCONT generates a random two-way table with given marginal totals.
!
!  Discussion:
!
!    Each time the program is called, another table will be randomly
!    generated.
!
!    Note that it should be the case that the sum of the row totals
!    is equal to the sum of the column totals.  However, this program
!    does not check for that condition.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    James Boyett
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    James Boyett,
!    Algorithm AS 144:
!    Random R x C Tables with Given Row and Column Totals,
!    Applied Statistics,
!    Volume 28, Number 3, pages 329-332, 1979.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the number of rows in the observed
!    matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns in the observed
!    matrix.
!
!    Input, integer ( kind = 4 ) NROWT(NROW), the row totals of the observed
!    matrix.
!
!    Input, integer ( kind = 4 ) NCOLT(NCOL), the column totals of the
!    observed matrix.
!
!    Input/output, integer ( kind = 4 ) NSUBT(NCOL), used by RCONT for partial
!    column sums.  Must not be changed by the calling program.
!
!    Output, integer ( kind = 4 ) MATRIX(NROW,NCOL), the random matrix.
!
!    Input/output, logical KEY, should be set to FALSE by the user before
!    the initial call.  RCONT will reset it to TRUE, and it should be left
!    at that value for subsequent calls in which the same values of NROW,
!    NCOL, NROWT and NCOLT are being used.
!
!    Output, integer ( kind = 4 ) IFAULT, fault indicator.
!    0, no error occured.
!    1, NROW <= 0.
!    2, NCOL <= 1.
!    3, some entry of NROWT is less than 0.
!    4, some entry of NCOLT is less than 0.
!    5, the sample size (sum of the column totals) is too large.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ), parameter :: nvec_max = 200

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical key
  integer ( kind = 4 ) limit
  integer ( kind = 4 ) matrix(nrow,ncol)
  integer ( kind = 4 ) ncolt(ncol)
  integer ( kind = 4 ) nnvect(nvec_max)
  integer ( kind = 4 ) noct
  integer ( kind = 4 ) nrowt(nrow)
  integer ( kind = 4 ) nsubt(ncol)
  integer ( kind = 4 ) ntemp
  integer ( kind = 4 ), save :: ntotal
  integer ( kind = 4 ), save, dimension ( nvec_max ) :: nvect
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ), save :: seed = 0

  ifault = 0

  if ( .not. key ) then
!
!  Set KEY for subsequent calls.
!
    key = .true.
    seed = 123456789
!
!  Check for faults and prepare for future calls.
!
    if ( nrow <= 0 ) then
      ifault = 1
      return
    end if

    if ( ncol <= 1 ) then
      ifault = 2
      return
    end if

    do i = 1, nrow
      if ( nrowt(i) <= 0 ) then
        ifault = 3
        return
      end if
    end do

    if ( ncolt(1) <= 0 ) then
      ifault = 4
      return
    end if

    nsubt(1) = ncolt(1)

    do j = 2, ncol

      if ( ncolt(j) <= 0 ) then
        ifault = 4
        return
      end if

      nsubt(j) = nsubt(j-1) + ncolt(j)

    end do

    ntotal = nsubt(ncol)

    if ( nvec_max < ntotal ) then
      ifault = 5
      return
    end if
!
!  Initialize vector to be permuted.
!
    do i = 1, ntotal
      nvect(i) = i
    end do

  end if
!
!  Initialize vector to be permuted.
!
  nnvect(1:ntotal) = nvect(1:ntotal)
!
!  Permute the vector.
!
  ntemp = ntotal
  do i = 1, ntotal
    noct = int ( r8_uniform_01 ( seed ) * real ( ntemp, kind = 8 ) + 1.0D+00 )
    nvect(i) = nnvect(noct)
    nnvect(noct) = nnvect(ntemp)
    ntemp = ntemp - 1
  end do
!
!  Construct the random matrix.
!
  matrix(1:nrow,1:ncol) = 0

  ii = 1

  do i = 1, nrow

    limit = nrowt(i)

    do k = 1, limit

      do j = 1, ncol
        if ( nvect(ii) <= nsubt(j) ) then
          ii = ii + 1
          matrix(i,j) = matrix(i,j) + 1
          exit
        end if
      end do

    end do

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
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
