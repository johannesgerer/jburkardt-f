subroutine enum ( m, n, rowsum, colsum, eval, ifault )

!*****************************************************************************80
!
!! ENUM generates contingency tables with given shape and row and column sums.
!
!  Discussion:
!
!   The routine enumerates all M by N contingency tables with given row 
!   and column totals, and calculates the hypergeometric probability of 
!   each table.
!
!   For tables having two or more row sums repeated, equivalent
!   tables differing only by a row permutation are not separately
!   enumerated.  A representative of each equivalence class is enumerated
!   and the multiplicity of each class calculated.
!
!   For each table enumerated, subroutine EVAL is called to carry out
!   calculations on the table.
!
!   Note that the entries in the column sum and row sum vectors will
!   be (implicitly) sorted into ascending order, and results will be 
!   returned as though these orderings were being used!
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
!    Original FORTRAN77 version by Ian Saunders.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Ian Saunders,
!    Algorithm AS 205,
!    Enumeration of R x C Tables with Repeated Row Totals,
!    Applied Statistics,
!    Volume 33, Number 3, 1984, pages 340-352. 
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, integer ( kind = 4 ) ROWSUM(M), the row sums.
!
!    Input, integer ( kind = 4 ) COLSUM(N), the column sums.
!
!    Input, external EVAL, the name of the user supplied routine
!    that will be called each time a new table is determined.
!    The routine has the form:
!      subroutine eval ( iflag, table, m, n, rowsum, colsum, prob, mult )
!    See the dummy version of EVAL for an explanation of the meaning
!    of the arguments.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) BOUND(M_MAX,M_MAX).  BOUND(I,J) is the
!    current upper bound on TABLE2(I,J) to satisfy row and column totals.
!
!    Local, integer ( kind = 4 ) M_MAX, the maximum dimension of TABLE2.
!
!    Local, integer ( kind = 4 ) MULT2(M_MAX), MULT2(I) is the maximum number 
!    of equivalent tables given the first I rows.
!
!    Local, integer ( kind = 4 ) NTOTAL, the total of the table entries.
!
!    Local, real ( kind = 8 ) PROB2(M_MAX,M_MAX), the partial sum of terms 
!    in log(p).
!
!    Local, integer ( kind = 4 ) REPS(M_MAX), the number of previous rows equal 
!    to row I.
!
!    Local, logical REPT(M_MAX), is true if row totals ROWSUM(I) and N(I-1) 
!    are equal.
!
!    Local, integer ( kind = 4 ) TABLE2(M_MAX,M_MAX), the current table.
!
!    Local, integer ( kind = 4 ) TOTAL_MAX, the maximum number of observations.
!
!    Local, integer ( kind = 4 ) Z(M_MAX), the lower bound on the sum of 
!    entries used by algorithm C.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: m_max = 10
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: total_max = 201

  integer ( kind = 4 ) bound(m_max,m_max)
  integer ( kind = 4 ) colsum(n)
  integer ( kind = 4 ) colsum2(m_max)
  external             eval
  real ( kind = 8 ) factlm(total_max+1)
  real ( kind = 8 ) flogm(total_max+1)
  integer ( kind = 4 ) i
  logical              ieqim
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iim
  integer ( kind = 4 ) iip
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jjj
  integer ( kind = 4 ) jjm
  integer ( kind = 4 ) jkeep
  integer ( kind = 4 ) jnext
  integer ( kind = 4 ) k
  integer ( kind = 4 ) left
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) maxrc
  integer ( kind = 4 ) mtotal
  integer ( kind = 4 ) mult
  integer ( kind = 4 ) mult2(m_max)
  integer ( kind = 4 ) multc
  integer ( kind = 4 ) multr
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) ntotal
  real ( kind = 8 ) prob
  real ( kind = 8 ) prob0
  real ( kind = 8 ) prob2(m_max,m_max)
  integer ( kind = 4 ) reps(m_max)
  integer ( kind = 4 ) repsc
  integer ( kind = 4 ) repsr
  logical              rept(m_max)
  logical              reptc(m_max)
  integer ( kind = 4 ) rowbnd
  integer ( kind = 4 ) rsum
  integer ( kind = 4 ) rowsum(m)
  integer ( kind = 4 ) rowsum2(m_max)
  integer ( kind = 4 ) table(m,n)
  integer ( kind = 4 ) table2(m_max,m_max)
  integer ( kind = 4 ) z(m_max)
!
!  Check the input values.
!
  ifault = 0

  if ( m <= 0 .or. m_max < m ) then
    ifault = 1
    return
  end if

  if ( n <= 0 .or. m_max < n ) then
    ifault = 1
    return
  end if

  ntotal = 0
  do i = 1, m
    if ( rowsum(i) <= 0 ) then
      ifault = 3
      return
    end if
    ntotal = ntotal + rowsum(i)
  end do

  mtotal = 0
  do j = 1, n
    if ( colsum(j) <= 0 ) then
      ifault = 3
      return
    end if
    mtotal = mtotal + colsum(j)
  end do

  if ( mtotal /= ntotal ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ENUM - Fatal error!'
    write ( *, '(a)' ) '  Row and column sums are not equal.'
    write ( *, '(a,i8)' ) '  Row sum =    ', ntotal
    write ( *, '(a,i8)' ) '  Column sum = ', mtotal
    ifault = 2
    return
  end if

  if ( total_max < ntotal ) then
    ifault = 4
    return
  end if
!
!  Copy the parameters.
!
  m2 = m
  n2 = n
  rowsum2(1:m) = rowsum(1:m)
  colsum2(1:n) = colsum(1:n)
!
!  The first call to EVAL is simply to allow EVAL to initialize.
!
  iflag = 1
  prob = 0.0D+00
  mult = 0
  table(1:m,1:n) = 0
  call eval ( iflag, table, m, n, rowsum, colsum, prob, mult ) 
!
!  Initialize 
!    flogm(k) = log(k-1), 
!    factlm(k) = log( (k-1)! )
!
  flogm(1) = 0.0D+00
  factlm(1) = 0.0D+00
  do k = 1, ntotal
    flogm(k+1) = log ( real ( k, kind = 8 ) )
    factlm(k+1) = factlm(k) + flogm(k+1)
  end do
!
!  Sort the rows and columns into ascending order.
!
  do i = 1, m2-1
    do ii = i+1, m2
      if ( rowsum2(ii) < rowsum2(i) ) then
        call i4_swap ( rowsum2(i), rowsum2(ii) )
      end if
    end do
  end do

  do j = 1, n2-1
    do jj = j+1, n2
      if ( colsum2(jj) < colsum2(j) ) then
        call i4_swap ( colsum2(j), colsum2(jj) )
      end if
    end do
  end do
!
!  Calculate the multiplicities of rows and columns.
!
!  reptc(j) = .true. if columns J and J-1 have the same total.
!  rept(i)  = .true. if rows I and I-1 have the same total.
!
  multc = 1.0D+00
  repsc = 1.0D+00
  reptc(1) = .false.

  do j = 2, n2
    reptc(j) = ( colsum2(j) == colsum2(j-1) )
    if ( reptc(j) ) then
      repsc = repsc + 1.0D+00
      multc = multc * repsc
    else
      repsc = 1.0D+00
    end if
  end do

  multr = 1.0D+00
  repsr = 1.0D+00
  rept(1) = .false.

  do i = 2, m2

    rept(i) = ( rowsum2(i) == rowsum2(i-1) )

    if ( rept(i) ) then
      repsr = repsr + 1.0D+00
      multr = multr * repsr
    else
      repsr = 1.0D+00
    end if

  end do
!
!  If column multiplicity exceeds row multiplicity, transpose the table.
!
  if ( multr < multc ) then

    maxrc = max ( m2, n2 )
    do ij = 1, maxrc
      call i4_swap ( rowsum2(ij), colsum2(ij) )
    end do

    call i4_swap ( m2, n2 )
    rept(1:m2) = reptc(1:m2)
    multr = multc

  end if
!
!  Set up the initial table.
!
!  Maximum multiplicity.
!
  mult2(1) = multr
  reps(1) = 1.0D+00
!
!  Constant term in probability.
!
  prob0 = - factlm(ntotal+1)

  do i = 1, m2
    ii = rowsum2(i)
    prob0 = prob0 + factlm(ii+1)
  end do

  do j = 1, n2
    jj = colsum2(j)
    prob0 = prob0 + factlm(jj+1)
  end do
!
!  Calculate bounds on row 1.
!
  bound(1,1:n2) = colsum2(1:n2)
!
!  For each I, find the greatest I-th row satisfying bounds.
!
  do i = 1, m2

    if ( i /= 1 ) then
      prob0 = prob2(i-1,n2)
    end if

    left = rowsum2(i)
!
!  Elements of row I:
!
    ieqim = rept(i)

    do j = 1, n2 - 1

      ij = min ( left, bound(i,j) )
      table2(i,j) = ij

      if ( j == 1 ) then
        prob2(i,j) = prob0 - factlm(ij+1)
      else
        prob2(i,j) = prob2(i,j-1) - factlm(ij+1)
      end if

      left = left - table2(i,j)

      if ( i < m2 ) then
        bound(i+1,j) = bound(i,j) - table2(i,j)
      end if

      if ( left == 0 ) then

        do jj = j+1, n2
          table2(i,jj) = 0
          prob2(i,jj) = prob2(i,jj-1)
          bound(i+1,jj) = bound(i,jj)
        end do

        exit

      end if

      if ( ieqim ) then
        ieqim = table2(i,j) == table2(i-1,j)
      end if

    end do

    if ( left /= 0 ) then
      table2(i,n2) = left
      prob2(i,n2) = prob2(i,n2-1) - factlm(left+1)
      if ( i < m2 ) then
         bound(i+1,n2) = bound(i,n2) - left
      end if

    end if

    if ( 1 < i ) then

      mult2(i) = mult2(i-1)
      reps(i) = 1.0D+00

      if ( ieqim ) then
        reps(i) = reps(i-1) + 1.0D+00
        mult2(i) = mult2(i) / reps(i)
      end if

    end if

  end do
!
!  Call EVAL for the first table.
!
  iflag = 2
  prob = prob2(m2,n2)
  mult = mult2(m2)

  if ( m == m2 .and. n == n2 ) then
    table(1:m,1:n) = table2(1:m2,1:n2)
  else
    table(1:m,1:n) = transpose ( table2(1:m2,1:n2) )
  end if

  call eval ( iflag, table, m, n, rowsum, colsum, prob, mult ) 
!
!  Enumerate the remaining tables.
!
!  Start of main loop.
!
  do
   
    i = m2

210 continue

    i = i - 1
!
!  If I = 0, no more tables are possible.
!
    if ( i == 0 ) then
      exit
    end if

    j = n2 - 1
    left = table2(i,n2)
    rowbnd = bound(i,n2)
!
!  Try to decrease element (I,J).
!
    do
!
!  Decrease element (I,J).
!
      if ( 0 < table2(i,j) .and. left < rowbnd ) then

        ij = table2(i,j)
        prob2(i,j) = prob2(i,j) + flogm(ij+1)
        table2(i,j) = table2(i,j) - 1
        bound(i+1,j) = bound(i+1,j) + 1
!
!  If row I was the same as row I-1, it is no longer.
!
        if ( reps(i) /= 1.0D+00 ) then
          reps(i) = 1.0D+00
          mult2(i) = mult2(i-1)
        end if
!
!  Complete row I with the largest possible values.
!
        ii = i
        iip = ii + 1
        iim = ii - 1
        jnext = j + 1
        left = left + 1
        go to 380
!
!  Element (I,J) cannot be decreased.  Try (I,J-1).
!
      else if ( j <= 1 ) then

        exit

      else

        left = left + table2(i,j)
        rowbnd = rowbnd + bound(i,j)
        j = j - 1

      end if

    end do

    go to 210
!
!  Fill up the remaining rows.
!
300 continue

    ii = ii + 1
!
!  The last row is treated separately.
!
    if ( ii == m2 ) then
      go to 400
    end if

    iip = ii + 1
    iim = ii - 1
!
!  Row total ROWSUM2(II) is not a repeat.  Make row II as large as possible.
!
    if ( .not. rept(ii) ) then
      left = rowsum2(ii)
      jnext = 1
      go to 380
    end if
!
!  Repeated row totals.
!
!   (i) if row II-1 satisfies the bounds on row II repeat it
!
310 continue
  
    do j = 1, n2

      if ( bound(ii,j) < table2(iim,j) ) then
        go to 330
      end if

      ij = table2(iim,j)
      table2(ii,j) = ij
      bound(iip,j) = bound(ii,j) - table2(ii,j)

      if ( j == 1 ) then
        prob2(ii,j) = prob2(iim,n2) - factlm(ij+1)
      else
        prob2(ii,j) = prob2(ii,j-1) - factlm(ij+1)
      end if

    end do
!
!  Row II is a repeat of row II-1.
!
    reps(ii) = reps(iim) + 1.0D+00
    mult2(ii) = mult2(iim) / reps(ii)
    go to 300
!
!  Element J of row II-1 was too big.
!
!  Construct the sequence z(j) of lower bounds
!
330 continue
!
!  If J=1 the bounds are satisfied automatically.
!
    if ( j == 1 ) then
      ij = bound(ii,1)
      table2(ii,1) = ij
      prob2(ii,1) = prob2(iim,n2) - factlm(ij+1)
      jnext = 2
      left = rowsum2(ii) - table2(ii,1)
      bound(iip,1) = 0
      go to 380
    end if

    z(j) = rowsum2(ii)

    do jj = j+1, n2
      z(j) = z(j) - bound(ii,jj)
    end do

    do jjm = 1, j-1
      jj = j - jjm
      z(jj) = z(jj+1) - bound(ii,jj+1)
    end do
!
!  (ii) if the cumulative totals of row II-1 all exceed the bounds Z(J),
!  make element (II,J) equal to its bound.
!
    rsum = 0
    jkeep = 0

    do jj = 1, j-1

      rsum = rsum + table2(iim,jj)
!
!  (III) the cumulative sums violate the bounds.
!  If no element of row II-1 can be changed to satisfy the bounds,
!  then no suitable row II is possible.
!  In that case go back and try decreasing row II-1.
!
      if ( rsum < z(jj) ) then

        if ( jkeep == 0 ) then
          i = ii
          go to 210
        end if
!
!  Element (II,JKEEP) can be decreased.
!
        bound(iip,jkeep) = bound(iip,jkeep) + 1
        ij = table2(ii,jkeep)
        prob2(ii,jkeep) = prob2(ii,jkeep) + flogm(ij+1)
        table2(ii,jkeep) = table2(ii,jkeep) - 1
!
!  Complete the row.
!
        jnext = jkeep + 1
        left = rowsum2(ii)

        do jjj = 1, jkeep
          left = left - table2(ii,jjj)
        end do

        go to 380

      end  if

      if ( z(jj) < rsum .and. 0 < table2(iim,jj) ) then
        jkeep = jj
      end if

    end do

    table2(ii,j) = bound(ii,j)
    bound(iip,j) = 0
    ij = table2(ii,j)
    prob2(ii,j) = prob2(ii,j-1) - factlm(ij+1)
    reps(ii) = 1.0D+00
    mult2(ii) = mult2(iim)
!
!  Complete row II with the largest possible elements.
!
    jnext = j + 1
    left = rowsum2(ii)
    do jj = 1, j
      left = left - table2(ii,jj)
    end do
!
!  Row II is complete up to element JNEXT-1.
!  Make the remaining elements as large as possible.
!  This section of code is used for every row, repeated or not.
!
380 continue

    do j = jnext, n2 - 1

      table2(ii,j) = min ( left, bound(ii,j) )
      left = left - table2(ii,j)
      bound(iip,j) = bound(ii,j) - table2(ii,j)
      ij = table2(ii,j)

      if ( j == 1 ) then
        prob2(ii,j) = prob2(iim,n2) - factlm(ij+1)
      else
        prob2(ii,j) = prob2(ii,j-1) - factlm(ij+1)
      end if

      if ( left == 0 ) then
        do jj = j+1, n2
          table2(ii,jj) = 0
          prob2(ii,jj) = prob2(ii,jj-1)
          bound(iip,jj) = bound(ii,jj)
        end do
        exit
      end if

    end do

    if ( left .ne. 0 ) then
      table2(ii,n2) = left
      prob2(ii,n2) = prob2(ii,n2-1) - factlm(left+1)
      bound(iip,n2) = bound(ii,n2) - left
    end if

    reps(ii) = 1.0D+00

    if ( 1 < ii ) then
      mult2(ii) = mult2(iim)
    end if

    go to 300
!
!  The final row.
!
400 continue
!
!  If not a repeat, set row M2 equal to its bounds.
!
    if ( .not. rept(m2) ) then

      ij = bound(m2,1)
      table2(m2,1) = ij
      prob2(m2,1) = prob2(m2-1,n2) - factlm(ij+1)

      do j = 2, n2
        ij = bound(m2,j)
        table2(m2,j) = ij
        prob2(m2,j) = prob2(m2,j-1) - factlm(ij+1)
      end do

      mult2(m2) = mult2(m2-1)
      go to 500

    end if
!
!  Row total M2 is a repeat - ensure that it is less than row M2-1.
!
    do j = 1, n2
!
!  If row M2 would be bigger than row M2-1, go back and try
!  decreasing row M2-2.
!
      if ( table2(m2-1,j) < bound(m2,j) ) then
        i = m2 - 1
        go to 210
      end if

      ij = bound(m2,j)
      table2(m2,j) = ij

      if ( j == 1 ) then
        prob2(m2,j) = prob2(m2-1,n2) - factlm(ij+1)
      else
        prob2(m2,j) = prob2(m2,j-1) - factlm(ij+1)
      end if
!
!  Row M2 is already less then row M2-1, so no more checks are needed.
!
      if ( table2(m2,j) /= table2(m2-1,j) ) then

        do jj = j + 1, n2
          ij = bound(m2,jj)
          table2(m2,jj) = ij
          prob2(m2,jj) = prob2(m2,jj-1) - factlm(ij+1)
        end do

        mult2(m2) = mult2(m2-1)
        go to 500

      end if

    end do
!
!  Row M2 is a repeat of row M2-1.
!
    reps(m2) = reps(m2-1) + 1.0D+00
    mult2(m2) = mult2(m2-1) / reps(m2)
!
!  The table is complete.
!
500 continue

    iflag = 2
    prob = prob2(m2,n2)
    mult = mult2(m2)

    if ( m == m2 .and. n == n2 ) then
      table(1:m,1:n) = table2(1:m2,1:n2)
    else
      table(1:m,1:n) = transpose ( table2(1:m2,1:n2) )
    end if

    call eval ( iflag, table, m, n, rowsum, colsum, prob, mult )

  end do

  iflag = 3
  prob = 0.0D+00
  mult = 0
  table(1:m,1:n) = 0
  call eval ( iflag, table, m, n, rowsum, colsum, prob, mult ) 

  return
end
subroutine eval ( iflag, table, m, n, rowsum, colsum, prob, mult )

!*****************************************************************************80
!
!! EVAL is called by ENUM every time a new contingency table is determined.
!
!  Discussion:
!
!    This is a dummy version of the routine.
!
!    The user might wish to print out each contingency table, or collect
!    some statistics.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Saunders,
!    Algorithm AS 205,
!    Enumeration of R x C Tables with Repeated Row Totals,
!    Applied Statistics,
!    Volume 33, Number 3, pages 340-352, 1984. 
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
!    Input, real ( kind = 8 ) PROB, the logarithm of the hypergeometric 
!    probability of this table.
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
    psum = psum + real ( mult, kind = 8 ) * exp ( prob )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVAL:'
    write ( *, '(i3,i3,g14.6)' ) count1, mult, prob

    call i4mat_print ( m, n, table, '  Table' )
!
!  Last call, no table.
!
  else if ( iflag == 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVAL summary'
    write ( *, '(a,i8)' ) '  Number of cases (ignoring multiplicity):', count1
    write ( *, '(a,i8)' ) '  Number of cases (allowing multiplicity):', count2
    write ( *, '(a,g14.6)' ) '  Probability sum = ', psum

  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 May 2000
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
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do jlo = 1, n, 10
    jhi = min ( jlo + 9, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,10(i7))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,10i7)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 1999
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
    write ( *, '(i6,i10)' ) i, a(i)
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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

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
