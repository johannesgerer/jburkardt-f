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
    write ( *, '(2x,i8,2x,i12)' ) i, a(i)
  end do

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
subroutine rcont2 ( nrow, ncol, nrowt, ncolt, key, seed, matrix, ierror )

!*****************************************************************************80
!
!! RCONT2 constructs a random two-way contingency table with given sums.
!
!  Discussion:
!
!    It is possible to specify row and column sum vectors which
!    correspond to no table at all.  As far as I can see, this routine does
!    not detect such a case.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by WM Patefield.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    WM Patefield,
!    Algorithm AS 159:
!    An Efficient Method of Generating RXC Tables with
!    Given Row and Column Totals,
!    Applied Statistics,
!    Volume 30, Number 1, 1981, pages 91-97.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns
!    in the table.  NROW and NCOL must each be at least 2.
!
!    Input, integer ( kind = 4 ) NROWT(NROW), NCOLT(NCOL), the row and column
!    sums.  Each entry must be positive.
!
!    Input/output, logical KEY, a flag that indicates whether data has
!    been initialized for this problem.  Set KEY = .FALSE. before the first
!    call.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) MATRIX(NROW,NCOL), the matrix.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is returned
!    as 0 if no error occurred.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtot = 5000

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  logical              done1
  logical              done2
  real ( kind = 8 ), save, dimension ( maxtot+1 ) :: fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) iap
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) id
  integer ( kind = 4 ) idp
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) igp
  integer ( kind = 4 ) ihp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) jwork(ncol)
  logical              key
  integer ( kind = 4 ) l
  logical              lsm
  logical              lsp
  integer ( kind = 4 ) m
  integer ( kind = 4 ) matrix(nrow,ncol)
  integer ( kind = 4 ) ncolt(ncol)
  integer ( kind = 4 ) nll
  integer ( kind = 4 ) nlm
  integer ( kind = 4 ) nlmp
  integer ( kind = 4 ) nrowt(nrow)
  integer ( kind = 4 ) nrowtl
  integer ( kind = 4 ), save :: ntotal = 0
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sumprb
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  ierror = 0
!
!  On user's signal, set up the factorial table.
!
  if ( .not. key ) then

    key = .true.

    if ( nrow <= 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RCONT - Fatal error!'
      write ( *, '(a)' ) '  Input number of rows is less than 2.'
      ierror = 1
      return
    end if

    if ( ncol <= 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RCONT - Fatal error!'
      write ( *, '(a)' ) '  The number of columns is less than 2.'
      ierror = 2
      return
    end if

    do i = 1, nrow
      if ( nrowt(i) <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RCONT - Fatal error!'
        write ( *, '(a)' ) '  An entry in the row sum vector is not positive.'
        ierror = 3
        return
      end if
    end do

    do j = 1, ncol
      if ( ncolt(j) <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RCONT - Fatal error!'
        write ( *, '(a)' ) &
         '  An entry in the column sum vector is not positive.'
        ierror = 4
        return
      end if
    end do

    if ( sum ( ncolt(1:ncol) ) /= sum ( nrowt(1:nrow) ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RCONT - Fatal error!'
      write ( *, '(a)' ) &
        '  The row and column sum vectors do not have the same sum.'
      ierror = 6
      return
    end if

    ntotal = sum ( ncolt(1:ncol) )

    if ( maxtot < ntotal ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RCONT - Fatal error!'
      write ( *, '(a)' ) &
        '  The sum of the column sum vector entries is too large.'
      ierror = 5
      return
    end if
!
!  Calculate log-factorials.
!
    x = 0.0D+00
    fact(1) = 0.0D+00
    do i = 1, ntotal
      x = x + log ( real ( i, kind = 8 ) )
      fact(i+1) = x
    end do

  end if
!
!  Construct a random matrix.
!
  jwork(1:ncol-1) = ncolt(1:ncol-1)

  jc = ntotal

  do l = 1, nrow - 1

    nrowtl = nrowt(l)
    ia = nrowtl
    ic = jc
    jc = jc - nrowtl

    do m = 1, ncol - 1

      id = jwork(m)
      ie = ic
      ic = ic - id
      ib = ie - ia
      ii = ib - id
!
!  Test for zero entries in matrix.
!
      if ( ie == 0 ) then
        ia = 0
        matrix(l,m:ncol) = 0
        exit
      end if
!
!  Generate a pseudo-random number.
!
      r = r8_uniform_01 ( seed )
!
!  Compute the conditional expected value of MATRIX(L,M).
!
      done1 = .false.

      do

        nlm = int ( &
          real ( ia * id, kind = 8 ) / real ( ie , kind = 8 ) + 0.5D+00 )

        iap = ia + 1
        idp = id + 1
        igp = idp - nlm
        ihp = iap - nlm
        nlmp = nlm + 1
        iip = ii + nlmp
        x = exp ( fact(iap) + fact(ib+1) + fact(ic+1) + fact(idp) - &
          fact(ie+1) - fact(nlmp) - fact(igp) - fact(ihp) - fact(iip) )

        if ( r <= x ) then
          exit
        end if

        sumprb = x
        y = x
        nll = nlm
        lsp = .false.
        lsm = .false.
!
!  Increment entry in row L, column M.
!
        do while ( .not. lsp )

          j = ( id - nlm ) * ( ia - nlm )

          if ( j == 0 ) then

            lsp = .true.

          else

            nlm = nlm + 1
            x = x * real ( j, kind = 8 ) &
              / real ( nlm * ( ii + nlm ), kind = 8 )
            sumprb = sumprb + x

            if ( r <= sumprb ) then
              done1 = .true.
              exit
            end if

          end if

          done2 = .false.

          do while ( .not. lsm )
!
!  Decrement the entry in row L, column M.
!
            j = nll * ( ii + nll )

            if ( j == 0 ) then
              lsm = .true.
              exit
            end if

            nll = nll - 1
            y = y * real ( j, kind = 8 ) &
              / real ( ( id - nll ) * ( ia - nll ), kind = 8 )
            sumprb = sumprb + y

            if ( r <= sumprb ) then
              nlm = nll
              done2 = .true.
              exit
            end if

            if ( .not. lsp ) then
              exit
            end if

          end do

          if ( done2 ) then
            exit
          end if

        end do

        if ( done1 ) then
          exit
        end if

        if ( done2 ) then
          exit
        end if

        r = r8_uniform_01 ( seed )
        r = sumprb * r

      end do

      matrix(l,m) = nlm
      ia = ia - nlm
      jwork(m) = jwork(m) - nlm

    end do

    matrix(l,ncol) = ia

  end do
!
!  Compute the last row.
!
  matrix(nrow,1:ncol-1) = jwork(1:ncol-1)
  matrix(nrow,ncol) = ib - matrix(nrow,ncol-1)

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
