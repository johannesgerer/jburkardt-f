subroutine addcom ( isub, jsub, value, perm_inv, diag, xlnz, ixlnz, nzsub, &
  xnzsub, n )

!*****************************************************************************80
!
!! ADDCOM adds values to a matrix stored in compressed storage scheme.
!
!  Discussion:
!
!    The routine is called once the equations and variables have
!    been reordered, to add numbers to the matrix.
!
!    The arrays DIAG and XLNZ, which are used as storage for array entries,
!    should be zeroed out before this routine is first called for a
!    particular problem.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ISUB, JSUB, the row and column of the matrix to
!    which the value is to be added.  Note that the matrix is assumed
!    to be symmetric, and only the lower triangle is stored.  Hence,
!    if ISUB < JSUB, the routine ignores the request.
!
!    Input, real ( kind = 8 ) VALUE, the quantity to be added to A(ISUB,JSUB).
!
!    Input, integer ( kind = 4 ) PERM_INV(N), the inverse ordering, which
!    should have been created by calling PERM_INVERSE once the reordering
!    is set.
!
!    Input/output, real ( kind = 8 ) DIAG(N), the diagonal elements of the
!    matrix.
!
!    Input/output, real ( kind = 8 ) XLNZ(*), the nonzero subdiagonal elements
!    of the matrix.
!
!    Input, integer ( kind = 4 ) IXLNZ(N+1), NZSUB(*), XNZSUB(N), data
!    structures which define the compressed storage scheme.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) ixlnz(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jsub
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kstop
  integer ( kind = 4 ) kstrt
  integer ( kind = 4 ) ksub
  integer ( kind = 4 ) nzsub(*)
  integer ( kind = 4 ) perm_inv(n)
  real ( kind = 8 ) value
  real ( kind = 8 ) xlnz(*)
  integer ( kind = 4 ) xnzsub(n)
!
!  Figure out the current locations of the given row and column.
!
  i = perm_inv(isub)
  j = perm_inv(jsub)
!
!  If the entry is on the diagonal, update DIAG.
!
  if ( i == j ) then
    diag(i) = diag(i) + value
    return
  end if
!
!  If the entry is above the diagonal, don't store it at all.
!
  if ( i < j ) then
    return
  end if

  kstrt = ixlnz(j)
  kstop = ixlnz(j+1) - 1

  if ( kstop < kstrt ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADDCOM - Fatal error!'
    write ( *, '(a)' ) '  The IXLNZ array is incorrect.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ISUB = ', isub
    write ( *, '(a,i8)' ) '  JSUB = ', jsub
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  I = ', i
    write ( *, '(a,i8)' ) '  J = ', j
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  IXLNZ(J) =     ', ixlnz(j)
    write ( *, '(a,i8)' ) '  IXLNZ(J+1)-1 = ', ixlnz(j+1) - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  However, we must have that'
    write ( *, '(a)' ) '  IXLNZ(J) <= IXLNZ(J+1)-1.'
    stop
  end if

  ksub = xnzsub(j)

  do k = kstrt, kstop

    if ( nzsub(ksub) == i ) then
      xlnz(k) = xlnz(k) + value
      return
    end if

    ksub = ksub + 1

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ADDCOM - Fatal error!'
  write ( *, '(a)' ) '  No storage was set aside for entry ISUB, JSUB'
  write ( *, '(a,i8)' ) '  ISUB = ', isub
  write ( *, '(a,i8)' ) '  JSUB = ', jsub

  stop
end
subroutine addrcm ( isub, jsub, value, perm_inv, diag, xenv, env, n )

!*****************************************************************************80
!
!! ADDRCM adds values to a matrix stored in the RCM scheme.
!
!  Discussion:
!
!    Since this routine only adds VALUE to the current matrix entry, the
!    matrix should be zeroed out first.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ISUB, JSUB, the row and column of the matrix to
!    which the value is to be added.  Note that the matrix is assumed
!    to be symmetric, and only the lower triangle is stored.  Hence,
!    if ISUB < JSUB, the routine ignores the request.
!
!    Input, real ( kind = 8 ) VALUE, the number to be added.
!
!    Input, integer ( kind = 4 ) PERM_INV(N), the inverse variable ordering.
!
!    Input/output, real ( kind = 8 ) DIAG(N), the diagonal of the matrix.
!
!    Input, integer ( kind = 4 ) XENV(N+1), describes the envelope structure.
!
!    Input/output, real ( kind = 8 ) ENV(*), the nonzeros of the matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) env(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jsub
  integer ( kind = 4 ) k
  integer ( kind = 4 ) perm_inv(n)
  real ( kind = 8 ) value
  integer ( kind = 4 ) xenv(n+1)

  i = perm_inv(isub)
  j = perm_inv(jsub)

  if ( i < j ) then
    return
  end if

  if ( i == j ) then
    diag(i) = diag(i) + value
    return
  end if

  k = xenv(i+1) - i + j

  if ( k < xenv(i) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADDRCM - Fatal error!'
    write ( *, '(a)' ) '  Indices outside of envelope.'
    write ( *, '(a,i8)' ) '  ISUB = ', isub
    write ( *, '(a,i8)' ) '  JSUB = ', jsub
    stop
  end if

  env(k) = env(k) + value

  return
end
subroutine addrhs ( perm_inv, isub, n, rhs, value )

!*****************************************************************************80
!
!! ADDRHS adds a quantity to a specific entry of the right hand side.
!
!  Discussion:
!
!    After the equations have been reordered, it may be desired to
!    alter one of the entries of the right hand side.  This routine
!    carries out this operation after adjusting for the reordering.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PERM_INV(N), the inverse permutation, constructed
!    by PERM_INVERSE.
!
!    Input, integer ( kind = 4 ) ISUB, the index of the entry of RHS to be modified.
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input/output, real ( kind = 8 ) RHS(N), the right hand side.  This vector
!    has been permuted.
!
!    Input, real ( kind = 8 ) VALUE, the quantity to be added to the right hand side.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) perm_inv(n)
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) value

  i = perm_inv(isub)

  if ( 1 <= i .and. i <= n ) then
    rhs(i) = rhs(i) + value
  end if

  return
end
subroutine addrqt ( isub, jsub, value, perm_inv, diag, xenv, env, xnonz, nonz, &
  nzsubs, n )

!*****************************************************************************80
!
!! ADDRQT adds values to a matrix stored in the implicit block storage scheme.
!
!  Discussion:
!
!    Since the routine only adds new values to those currently in storage, the
!    space used to store the matrix must be initialized to 0 before numerical
!    values are supplied.  This routine can be used with the RQT and 1WD
!    methods.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ISUB, JSUB, the row and column of the matrix to
!    which the value is to be added.  Note that the matrix is assumed
!    to be symmetric, and only the lower triangle is stored.  Hence,
!    if ISUB < JSUB, the routine ignores the request.
!
!    Input, real ( kind = 8 ) VALUE, the number to be added.
!
!    Input, integer ( kind = 4 ) PERM_INV(N), the inverse variable ordering.
!
!    Input/output, real ( kind = 8 ) DIAG(N), the diagonal of the matrix.
!
!    Input, integer ( kind = 4 ) XENV(N+1), describes the envelope structure
!    of the diagonal blocks.
!
!    Input/output, real ( kind = 8 ) ENV(*), the nonzeros of the matrix.
!
!    Input/output, integer ( kind = 4 ) XNONZ(N+1), real ( kind = 8 ) NONZ(*),
!    integer ( kind = 4 ) NZSUBS(*), level structure containing the off-block
!    diagonal parts of the rows of the lower triangle of the original matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) env(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jsub
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kstop
  integer ( kind = 4 ) kstrt
  real ( kind = 8 ) nonz(*)
  integer ( kind = 4 ) nzsubs(*)
  integer ( kind = 4 ) perm_inv(n)
  real ( kind = 8 ) value
  integer ( kind = 4 ) xenv(n+1)
  integer ( kind = 4 ) xnonz(n+1)

  i = perm_inv(isub)
  j = perm_inv(jsub)
!
!  Ignore superdiagonal entries.
!
  if ( i < j ) then
    return
  end if
!
!  Diagonal entries.
!
  if ( i == j ) then
    diag(i) = diag(i) + value
    return
  end if
!
!  Entries within the diagonal envelope.
!
  k = xenv(i+1) - i + j

  if ( xenv(i) <= k ) then
    env(k) = env(k) + value
    return
  end if
!
!  The value goes outside the diagonal blocks.
!
  kstrt = xnonz(i)
  kstop = xnonz(i+1) - 1

  do k = kstrt, kstop

    if ( nzsubs(k) == j ) then
      nonz(k) = nonz(k) + value
      return
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ADDRQT - Fatal error!'
  write ( *, '(a)' ) '  Lack of storage!'
  write ( *, '(a,i8)' ) '  ISUB = ', isub
  write ( *, '(a,i8)' ) '  JSUB = ', jsub

  stop
end
subroutine adj_env_size ( n, adj_row, adj_num, adj, perm, perm_inv, env_size )

!*****************************************************************************80
!
!! ADJ_ENV_SIZE computes the envelope size for an adjacency structure.
!
!  Discussion:
!
!    The matrix is assumed to be symmetric.
!
!    The variables are assumed to have been permuted.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) PERM(N), PERM_INV(N), the permutation
!    and inverse permutation.
!
!    Output, integer ( kind = 4 ) ENV_SIZE, the number of cells in the envelope.
!    You could think of this number as the sum of the
!    bandwidths of each row.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) add
  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) env_size
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  integer ( kind = 4 ) row

  env_size = 0

  do i = 1, n

    row = perm(i)

    add = 0

    do j = adj_row(row), adj_row(row+1) - 1
      col = perm_inv(adj(j))
      if ( col < i ) then
        add = max ( add, i - col )
      end if
    end do

    env_size = env_size + add

  end do

  return
end
subroutine adj_print ( n, adj_num, adj_row, adj )

!*****************************************************************************80
!
!! ADJ_PRINT prints the adjacency information stored in ADJ_ROW and ADJ.
!
!  Discussion:
!
!    The list has the form:
!
!    Row   Nonzeros
!
!    1       2   5   9
!    2       7   8   9   15   78   79   81  86  91  99
!          100 103
!    3      48  49  53
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the dimension of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1), organizes the entries of ADJ
!    into rows.  The entries for row I are in entries ADJ_ROW(I)
!    through ADJ_ROW(I+1)-1.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ADJ_PRINT'
  write ( *, '(a)' ) '  Show adjacency structure of sparse matrix.'
  write ( *, '(a,i8)' ) '  The matrix order is ', n
  write ( *, '(a,i8)' ) '  The number of entries is ', adj_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Row         Nonzeros '
  write ( *, '(a)' ) ' '

  do i = 1, n

    jmin = adj_row(i)
    jmax = adj_row(i+1) - 1
    do jlo = jmin, jmax, 10
      jhi = min ( jlo+9, jmax )
      if ( jlo == jmin ) then
        write ( *, '(i6,6x,10i6)' ) i, adj(jlo:jhi)
      else
        write ( *, '(6x,6x,10i6)' ) adj(jlo:jhi)
      end if
    end do

  end do

  return
end
subroutine adj_set ( adj, irow, jcol, adj_max, adj_num, n, adj_row )

!*****************************************************************************80
!
!! ADJ_SET sets up the adjacency information ADJ_ROW and ADJ.
!
!  Discussion:
!
!    The routine records the locations of each nonzero element,
!    one at a time.
!
!    The first call for a given problem should be with IROW or ICOL
!    negative.  This is a signal indicating the data structure should
!    be initialized.
!
!    Then, for each case in which A(IROW,JCOL) is nonzero, or
!    in which IROW is adjacent to JCOL, call this routine once
!    to record that fact.
!
!    Diagonal entries are not to be stored.
!
!    The matrix is assumed to be symmetric, so setting I adjacent to J
!    will also set J adjacent to I.
!
!    Repeated calls with the same values of IROW and JCOL do not
!    actually hurt.  No extra storage will be allocated.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input/output, integer ( kind = 4 ) ADJ(ADJ_MAX), stores adjacency information.
!
!    Input, integer ( kind = 4 ) IROW, JCOL, the row and column indices of a nonzero
!    entry of the matrix.
!
!    Input, integer ( kind = 4 ) ADJ_MAX, the maximum dimension of the adjacency array.
!
!    Input/output, integer ( kind = 4 ) ADJ_NUM, the number of adjaceny entries.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
  implicit none

  integer ( kind = 4 ) adj_max
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(adj_max)
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kback
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
!
!  Negative IROW or JCOL indicates the data structure should be initialized.
!
  if ( irow < 0 .or. jcol < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SET - Note:'
    write ( *, '(a)') '  Initializing adjacency information.'
    write ( *, '(a,i8)' ) '  Number of equations N = ', n
    write ( *, '(a,i8)' ) '  Maximum adjacency ADJ_MAX = ', adj_max

    adj_num = 0
    adj_row(1:n+1) = 1
    adj(1:adj_max) = 0

    return

  end if
!
!  Symmetric entries are not stored.
!
  if ( irow == jcol ) then
    return
  end if

  if ( n < irow ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SET - Fatal error!'
    write ( *, '(a)' ) '  N < IROW.'
    write ( *, '(a,i8)' ) '  IROW = ', irow
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  else if ( irow < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SET - Fatal error!'
    write ( *, '(a)' ) '  IROW < 1.'
    write ( *, '(a,i8)' ) '  IROW = ', irow
    stop
  else if ( n < jcol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SET - Fatal error!'
    write ( *, '(a)' ) '  N < JCOL.'
    write ( *, '(a,i8)' ) '  JCOL = ', jcol
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  else if ( jcol < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SET - Fatal error!'
    write ( *, '(a)' ) '  JCOL < 1.'
    write ( *, '(a,i8)' ) '  JCOL = ', jcol
    stop
  end if

  i = irow
  j = jcol

   20 continue
!
!  Search the adjacency entries already stored for row I,
!  to see if J has already been stored.
!
  klo = adj_row(i)
  khi = adj_row(i+1) - 1

  do k = klo, khi

    if ( adj(k) == j ) then

      if ( i == irow ) then
        i = jcol
        j = irow
        go to 20
      end if

      return

    end if

  end do
!
!  A new adjacency entry must be made.
!  Check that we're not exceeding the storage allocation for ADJ.
!
  if ( adj_max < adj_num + 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SET - Fatal error!'
    write ( *, '(a)' ) '  All available storage has been used.'
    write ( *, '(a)' ) '  No more information can be stored!'
    write ( *, '(a)' ) '  This error occurred for '
    write ( *, '(a,i8)' ) '  Row IROW = ',irow
    write ( *, '(a,i8)' ) '  Column JCOL = ',jcol
    stop
  end if
!
!  Shift the later adjacency entries up one.
!
  do k = adj_row(i+1), adj_row(n+1)
    kback = adj_row(n+1) + adj_row(i+1) - k
    adj(kback+1) = adj(kback)
  end do
!
!  Insert the new entry.
!
  adj(adj_row(i+1)) = j
!
!  Update the ADJ_ROW pointers.
!
  do k = i + 1, n + 1
    adj_row(k) = adj_row(k) + 1
  end do

  adj_num = adj_row(n+1) - 1

  if ( i == irow ) then
    i = jcol
    j = irow
    go to 20
  end if

  return
end
subroutine adj_show ( adj, iband, perm_inv, adj_num, n, perm, adj_row )

!*****************************************************************************80
!
!! ADJ_SHOW displays a symbolic picture of a matrix.
!
!  Discussion:
!
!    The matrix is defined by the adjacency information in adj_row and adj,
!    with a possible permutation through PERM and PERM_INV.  The routine
!    also computes the bandwidth and the size of the envelope.
!
!    If no permutation has been done, you must set PERM_INV(I) = PERM(I) = I
!    before calling this routine.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) IBAND, the bandwidth of the matrix.
!
!    Input, integer ( kind = 4 ) PERM(N), PERM_INV(N), the permutation
!    and inverse permutation.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries in ADJ.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
  implicit none

  integer, parameter :: n_max = 100

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_row(n+1)
  character band(n_max)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iband
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nonz
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)

  iband = 0
  nonz = 0

  if ( n_max < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SHOW - Fatal error!'
    write ( *, '(a)' ) '  N is too large!'
    write ( *, '(a,i8)' ) '  Maximum legal value is ', n_max
    write ( *, '(a,i8)' ) '  Your input value was ', n
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ADJ_SHOW'
  write ( *, '(a)' ) '  Display nonzero structure of matrix.'
  write ( *, '(a)' ) ' '

  do i = 1, n

    do k = 1, n
      band(k) = ' '
    end do

    band(i) = 'X'

    do j = adj_row(perm(i)), adj_row(perm(i)+1) - 1
      col = perm_inv(adj(j))
      if ( col < i ) then
        nonz = nonz + 1
      end if
      iband = max ( iband, i-col )
      band(col) = 'X'
    end do

    write ( *, '(i6,1x,100a1)' ) i, band(1:n)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Lower bandwidth = ', iband
  write ( *, '(a,i8,a)' ) '  Lower envelope contains ', nonz, ' nonzeros.'

  return
end
subroutine block_shuffle ( adj_row, adj, perm, nblks, xblk, n )

!*****************************************************************************80
!
!! BLOCK_SHUFFLE renumbers the nodes of each block to reduce its envelope.
!
!  Discussion:
!
!    Nodes in a block with no neighbors in previous blocks are renumbered
!    by RCM_SUB before the others.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input/output, integer ( kind = 4 ) PERM(N), the permutation vector, which is
!    updated on output.
!
!    Input, integer ( kind = 4 ) NBLKS, XBLK(NBLKS+1), the tree partitioning.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) BNUM(N) stores the block number of each variable.
!
!    Local, integer ( kind = 4 ) MASK(N) is used to mask the subgraph.
!
!    Local, integer ( kind = 4 ) SUBG(N) is used to prescribe a subgraph.
!
  implicit none

  integer ( kind = 4 ) nblks
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) bnum(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) istrt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) nabor
  integer ( kind = 4 ) nbrblk
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nsubg
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) subg(n)
  integer ( kind = 4 ) xblk(nblks+1)
!
!  Find the block number for each variable and initialize MASK.
!
  do k = 1, nblks

    do i = xblk(k), xblk(k+1) - 1
      node = perm(i)
      bnum(node) = k
      mask(node) = 0
    end do

  end do
!
!  For each block, find those nodes with no neighbors
!  in previous blocks and accumulate them in SUBG.
!  They will be renumbered before others in the block.
!
  do k = 1, nblks

    istrt = xblk(k)
    istop = xblk(k+1) - 1
    nsubg = 0

    do i = istrt, istop

      node = perm(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1

      if ( jstrt <= jstop ) then

        do j = jstrt, jstop
          nabor = adj(j)
          nbrblk = bnum(nabor)
          if ( nbrblk < k ) then
            go to 40
          end if
        end do

        nsubg = nsubg + 1
        subg(nsubg) = node
        ip = istrt + nsubg - 1
        perm(i) = perm(ip)

      end if

40    continue

    end do
!
!  RCM_SUB renumbers the subgraph stored in (NSUBG, SUBG).
!
    if ( 0 < nsubg ) then
      call rcm_sub ( adj_row, adj, nsubg, subg, perm(istrt), n, mask )
    end if

  end do

  return
end
subroutine degree ( root, adj_row, adj, mask, deg, ccsize, ls, n )

!*****************************************************************************80
!
!! DEGREE computes node degrees in a connected component, for the RCM method.
!
!  Discussion:
!
!    The connected component is specified by MASK and ROOT.
!    Nodes for which MASK is zero are ignored.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the node that defines the component.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) MASK(N), specifies a section subgraph.
!
!    Output, integer ( kind = 4 ) DEG(N), the degrees of the nodes in the component.
!
!    Output, integer ( kind = 4 ) CCSIZE, the size of the component specified by
!    MASK and ROOT.
!
!    Output, integer ( kind = 4 ) LS(N), stores the nodes of the component level by level.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) ccsize
  integer ( kind = 4 ) deg(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ideg
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) lbegin
  integer ( kind = 4 ) ls(n)
  integer ( kind = 4 ) lvlend
  integer ( kind = 4 ) lvsize
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) node
  integer ( kind = 4 ) root
!
!  The array ADJ_ROW is used as a temporary marker to
!  indicate which nodes have been considered so far.
!
  ls(1) = root
  adj_row(root) = -adj_row(root)
  lvlend = 0
  ccsize = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
  do

    lbegin = lvlend + 1
    lvlend = ccsize
!
!  Find the degrees of nodes in the current level,
!  and at the same time, generate the next level.
!
    do i = lbegin, lvlend

      node = ls(i)
      jstrt = -adj_row(node)
      jstop = abs ( adj_row(node+1) ) - 1
      ideg = 0

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then

          ideg = ideg + 1

          if ( 0 <= adj_row(nbr) ) then
            adj_row(nbr) = -adj_row(nbr)
            ccsize = ccsize + 1
            ls(ccsize) = nbr
          end if

        end if

      end do

      deg(node) = ideg

    end do
!
!  Compute the current level width.
!
    lvsize = ccsize - lvlend
!
!  If the current level width is nonzero, generate another level.
!
    if ( lvsize == 0 ) then
      exit
    end if

  end do
!
!  Reset ADJ_ROW to its correct sign and return.
!
  do i = 1, ccsize
    node = ls(i)
    adj_row(node) = -adj_row(node)
  end do

  return
end
subroutine el_solve ( n, xenv, env, diag, rhs )

!*****************************************************************************80
!
!! EL_SOLVE solves a lower triangular system stored in the envelope format.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, integer ( kind = 4 ) XENV(N+1), describes the envelope structure.
!
!    Input, real ( kind = 8 ) ENV(*), the nonzeros of the matrix.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal of the matrix.
!
!    Input/output, real ( kind = 8 ) RHS(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) env(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iband
  integer ( kind = 4 ) ifirst
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kstop
  integer ( kind = 4 ) kstrt
  integer ( kind = 4 ) l
  integer ( kind = 4 ) last
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) s
  integer ( kind = 4 ) xenv(n+1)
!
!  Find IFIRST, the index of the first nonzero in RHS.
!
  ifirst = 0

  do

    ifirst = ifirst + 1

    if ( rhs(ifirst) /= 0.0D+00 ) then
      exit
    end if

    if ( n <= ifirst ) then
      return
    end if

  end do
!
!  LAST contains the position of the most recently
!  computed nonzero component of the solution.
!
  last = 0

  do i = ifirst, n

    iband = xenv(i+1) - xenv(i)
    iband = min ( iband, i - 1 )

    s = rhs(i)
    l = i - iband
    rhs(i) = 0.0D+00
!
!  A row of the envelope is empty, or the corresponding
!  components of the solution are all zeros.
!
    if ( iband /= 0 .and. 1 <= last ) then

      kstrt = xenv(i+1) - iband
      kstop = xenv(i+1) - 1

      do k = kstrt, kstop
        s = s - env(k) * rhs(l)
        l = l + 1
      end do

    end if

    if ( s /= 0.0D+00 ) then
      rhs(i) = s / diag(i)
      last = i
    end if

  end do

  return
end
subroutine es_factor ( n, xenv, env, diag, ierror )

!*****************************************************************************80
!
!! ES_FACTOR factors a positive definite envelope matrix into L * L'.
!
!  Discussion:
!
!    The algorithm used is the standard bordering method.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) XENV(N+1), the envelope index vector.
!
!    Input/output, real ( kind = 8 ) ENV(*), on input, the envelope of A,
!    and on output the envelope of the factor L.
!
!    Input/output, real ( kind = 8 ) DIAG(N), on input the diagonal of A,
!    and on output the diagonal of L.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error, the factorization was carried out.
!    1, the matrix is not positive definite.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) env(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iband
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifirst
  integer ( kind = 4 ) ixenv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  real ( kind = 8 ) temp
  integer ( kind = 4 ) xenv(n+1)

  ierror = 0

  if ( diag(1) <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ES_FACTOR - Fatal error!'
    write ( *, '(a)' ) '  The matrix is not positive definite.'
    ierror = 1
    return
  end if

  diag(1) = sqrt ( diag(1) )
!
!  Loop over rows 2:N of the matrix
!
  do i = 2, n

    ixenv = xenv(i)
    iband = xenv(i+1) - ixenv
!
!  Compute row I of the triangular factor.
!
    temp = diag(i)

    if ( iband /= 0 ) then

      ifirst = i - iband

      call el_solve ( iband, xenv(ifirst), env, diag(ifirst), env(ixenv) )

      jstop = xenv(i+1) - 1

      do j = ixenv, jstop
        temp = temp - env(j)**2
      end do

    end if

    if ( temp <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ES_FACTOR - Fatal error!'
      write ( *, '(a)' ) '  The matrix is not positive definite.'
      ierror = 1
      return
    end if

    diag(i) = sqrt ( temp )

  end do

  return
end
subroutine eu_solve ( n, xenv, env, diag, rhs )

!*****************************************************************************80
!
!! EU_SOLVE solves an upper triangular system stored in the envelope format.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) XENV(N+1), describes the envelope structure.
!
!    Input, real ( kind = 8 ) ENV(*), the nonzeros of the matrix.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal of U.
!
!    Input/output, real ( kind = 8 ) RHS(N).  On input, the right hand side.
!    On output, the solution vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) env(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iband
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kstop
  integer ( kind = 4 ) kstrt
  integer ( kind = 4 ) l
  real ( kind = 8 ) rhs(n)
  integer ( kind = 4 ) xenv(n+1)

  i = n + 1

  do

    i = i - 1

    if ( i == 0 ) then
      exit
    end if

    if ( rhs(i) == 0.0D+00 ) then
      cycle
    end if

    rhs(i) = rhs(i) / diag(i)
    iband = xenv(i+1) - xenv(i)

    if ( i <= iband ) then
      iband = i - 1
    end if

    if ( iband == 0 ) then
      cycle
    end if

    kstrt = i - iband
    kstop = i - 1
    l = xenv(i+1) - iband

    do k = kstrt, kstop
      rhs(k) = rhs(k) - env(l) * rhs(i)
      l = l + 1
    end do

  end do

  return
end
subroutine fnbenv ( adj_row, adj, perm, perm_inv, nblks, xblk, xenv, env_size, &
  marker, rchset, n )

!*****************************************************************************80
!
!! FNBENV finds the envelope of the diagonal blocks of a Cholesky factor.
!
!  Discussion:
!
!    Specifically, this routine finds the exact envelope structure of
!    the diagonal blocks of the Cholesky factor of a permuted
!    partitioned matrix.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) PERM(N), PERM_INV(N), the permutation vector,
!    and its inverse.
!
!    Input, integer ( kind = 4 ) NBLKS, XBLK(NBLKS+1), the partitioning.
!
!    Output, integer ( kind = 4 ) XENV(N+1), the envelope index vector.
!
!    Output, integer ( kind = 4 ) ENV_SIZE, the size of the envelope.
!
!    Workspace, integer ( kind = 4 ) MASK(N), marks nodes that have been considered.
!
!    Workspace, integer ( kind = 4 ) MARKER(N), used by routine REACH.
!
!    Workspace, integer ( kind = 4 ) RCHSET(N), used by REACH.
!    Stores both reachable and neighborhood sets.
!
  implicit none

  integer ( kind = 4 ) nblks
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) blkbeg
  integer ( kind = 4 ) blkend
  integer ( kind = 4 ) i
  integer ( kind = 4 ) env_size
  integer ( kind = 4 ) ifirst
  integer ( kind = 4 ) inhd
  integer ( kind = 4 ) k
  integer ( kind = 4 ) marker(n)
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) newnhd
  integer ( kind = 4 ) nhdsze
  integer ( kind = 4 ) node
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  integer ( kind = 4 ) rchset(n)
  integer ( kind = 4 ) rchsze
  integer ( kind = 4 ) xblk(nblks+1)
  integer ( kind = 4 ) xenv(n+1)

  env_size = 1
!
!  MASK will mark nodes that have been considered.
!
  mask(1:n) = 0
  marker(1:n) = 1
!
!  Loop over the blocks.
!
  do k = 1, nblks

    nhdsze = 0
    blkbeg = xblk(k)
    blkend = xblk(k+1) - 1

    do i = xblk(k), xblk(k+1) - 1
      node = perm(i)
      marker(node) = 0
    end do
!
!  Loop through the nodes in the current block.
!
    do i = xblk(k), xblk(k+1) - 1
      node = perm(i)
      call reach ( node, adj_row, adj, mask, marker, rchsze, rchset(blkbeg), &
        newnhd, rchset(nhdsze+1), n )
      nhdsze = nhdsze + newnhd
      ifirst = marker(node)
      ifirst = perm_inv(ifirst)
      xenv(i) = env_size
      env_size = env_size + i - ifirst
    end do
!
!  Reset marker values of nodes in nbrhd set.
!
    do inhd = 1, nhdsze
      node = rchset(inhd)
      marker(node) = 0
    end do
!
!  Reset the marker and mask values of nodes in the current block.
!
    do i = blkbeg, blkend
      node = perm(i)
      marker(node) = 0
      mask(node) = 1
    end do

  end do

  xenv(n+1) = env_size
  env_size = env_size - 1

  return
end
subroutine fndsep ( root, adj_row, adj, mask, nsep, sep, n )

!*****************************************************************************80
!
!! FNDSEP finds a small separator for a connected component in a graph.
!
!  Discussion:
!
!    The connected component is specified by the MASK array.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the node that determines the masked component.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) NSEP, the number of variables in the separator.
!
!    Output, integer ( kind = 4 ) SEP(*), the separator nodes.
!
!    Input/output, integer ( kind = 4 ) MASK(N), nodes in the separator have their mask
!    values set to zero.
!
!  Local parameters:
!
!    Workspace, integer ( kind = 4 ) XLS(N+1),  LS(N), level structure pair
!    for level structure found by ROOT_FIND.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) ls(n)
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) midbeg
  integer ( kind = 4 ) midend
  integer ( kind = 4 ) midlvl
  integer ( kind = 4 ) mp1beg
  integer ( kind = 4 ) mp1end
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nlvl
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nsep
  integer ( kind = 4 ) root
  integer ( kind = 4 ) sep(*)
  integer ( kind = 4 ) xls(n+1)
!
!  Determine the level structure associated with ROOT.
!
  call root_find ( root, adj_row, adj, mask, nlvl, xls, ls, n )
!
!  If the number of levels is less than 3, return the whole component
!  as the separator.
!
  if ( nlvl < 3 ) then

    nsep = xls(nlvl+1) - 1

    do i = 1, nsep
      node = ls(i)
      sep(i) = node
      mask(node) = 0
    end do

    return

  end if
!
!  Find the middle level of the rooted level structure.
!
  midlvl = ( nlvl + 2 ) / 2
  midbeg = xls(midlvl)
  mp1beg = xls(midlvl+1)
  midend = mp1beg - 1
  mp1end = xls(midlvl+2) - 1
!
!  The separator is obtained by including only those
!  middle-level nodes with neighbors in the MIDDLE + 1
!  level.  ADJ_ROW is used temporarily to mark those
!  nodes in the MIDDLE + 1 level.
!
  do i = mp1beg, mp1end
    node = ls(i)
    adj_row(node) = -adj_row(node)
  end do

  nsep = 0

  do i = midbeg, midend

    node = ls(i)
    jstrt = adj_row(node)
    jstop = abs ( adj_row(node+1) ) - 1

    do j = jstrt, jstop

      nbr = adj(j)

      if ( adj_row(nbr) <= 0 ) then
        nsep = nsep + 1
        sep(nsep) = node
        mask(node) = 0
        exit
      end if

    end do

  end do
!
!  Reset ADJ_ROW to its correct sign.
!
  do i = mp1beg, mp1end
    node = ls(i)
    adj_row(node) = -adj_row(node)
  end do

  return
end
subroutine fnenv ( n, adj_row, adj, perm, perm_inv, xenv, env_size, iband )

!*****************************************************************************80
!
!! FNENV finds the envelope structure of a permuted matrix.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) PERM(N), PERM_INV(N), the permutation and inverse
!    permutation vectors.
!
!    Output, integer ( kind = 4 ) XENV(N+1), the index vector for the level structure
!    to be used to store the lower (or upper) envelope of the reordered matrix.
!
!    Output, integer ( kind = 4 ) ENV_SIZE, is equal to XENV(N+1)-1.
!
!    Output, integer ( kind = 4 ) IBAND, the bandwidth of the reordered matrix.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iband
  integer ( kind = 4 ) env_size
  integer ( kind = 4 ) ifirst
  integer ( kind = 4 ) iperm
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jband
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) nabor
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  integer ( kind = 4 ) xenv(n+1)

  iband = 0
  env_size = 1

  do i = 1, n

    xenv(i) = env_size
    iperm = perm(i)
    jstrt = adj_row(iperm)
    jstop = adj_row(iperm+1) - 1
!
!  Find the first nonzero in row I.
!
    if ( jstrt <= jstop ) then

      ifirst = i

      do j = jstrt, jstop

        nabor = adj(j)
        nabor = perm_inv(nabor)
        ifirst = min ( ifirst, nabor )

      end do

      jband = i - ifirst
      env_size = env_size + jband
      iband = max ( iband, jband )

    end if

  end do

  xenv(n+1) = env_size
  env_size = env_size - 1

  return
end
subroutine fnlvls ( root, adj_row, adj, nodlvl, nlvl, xls, ls, n )

!*****************************************************************************80
!
!! FNLVLS generates a rooted level structure, as part of the RQT method.
!
!  Discussion:
!
!    FNLVLS generates a rooted level structure for a masked connected
!    subgraph.  The structure is rooted at a pseudo-peripheral node.
!    The level numbers are recorded.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ROOT.  On input, with the array NODLVL, ROOT
!    specifies the component whose pseudo-peripheral node is to be found.
!    On output, ROOT is the value of that pseudo-peripheral node.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) NLVL, the number of levels in the level structure.
!
!    Output, integer ( kind = 4 ) XLS(N+1), LS(N), the level structure returned.
!
!    Input/output, integer ( kind = 4 ) NODLVL(N).  On input, it specifies a section
!    subgraph.  On return, the node level numbers.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ls(n)
  integer ( kind = 4 ) lvl
  integer ( kind = 4 ) nlvl
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nodlvl(n)
  integer ( kind = 4 ) root
  integer ( kind = 4 ) xls(n+1)

  call root_find ( root, adj_row, adj, nodlvl, nlvl, xls, ls, n )

  do lvl = 1, nlvl

    do j = xls(lvl), xls(lvl+1) - 1
      node = ls(j)
      nodlvl(node) = lvl
    end do

  end do

  return
end
subroutine fnofnz ( adj_row, adj, perm, perm_inv, nblks, xblk, xnonz, nzsubs, &
  maxnz, n )

!*****************************************************************************80
!
!! FNOFNZ finds columns of off-block-diagonal nonzeros.
!
!  Discussion:
!
!    It looks for these items in the lower triangle of a partitioned matrix.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) PERM(N), PERM_INV(N), the permutation and inverse
!    permutation vectors.
!
!    Input, integer ( kind = 4 ) NBLKS, XBLK(NBLKS+1), the block partitioning.
!
!    Output, integer ( kind = 4 ) XNONZ(N+1), NZSUBS(*), the column subscripts of
!    the nonzeros of the matrix to the left of the diagonal blocks are
!    stored row by row in contiguous locations in NZSUBS, and  XNONZ is
!    the index vector to NZSUBS.
!
!    Input/output, integer ( kind = 4 ) MAXNZ.  On input, the size of the vector
!    NZSUBS; and on output, the number of nonzeros found.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) nblks
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) blkbeg
  integer ( kind = 4 ) blkend
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jperm
  integer ( kind = 4 ) jxnonz
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kstop
  integer ( kind = 4 ) kstrt
  integer ( kind = 4 ) maxnz
  integer ( kind = 4 ) nabor
  integer ( kind = 4 ) nzcnt
  integer ( kind = 4 ) nzsubs(*)
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  integer ( kind = 4 ) xblk(nblks+1)
  integer ( kind = 4 ) xnonz(n+1)

  nzcnt = 1

  if ( nblks <= 0 ) then
    maxnz = 0
    return
  end if
!
!  Loop over the blocks.
!
  do i = 1, nblks

    blkbeg = xblk(i)
    blkend = xblk(i+1) - 1
!
!  Loop over the rows of the I-th block.
!
    do j = blkbeg, blkend

      xnonz(j) = nzcnt
      jperm = perm(j)
      kstrt = adj_row(jperm)
      kstop = adj_row(jperm+1) - 1
!
!  Loop over the nonzeros of row J.
!
      if ( kstrt <= kstop ) then

        do k = kstrt, kstop

          nabor = adj(k)
          nabor = perm_inv(nabor)
!
!  Check to see if it is to the left of the I-th diagonal block.
!
          if ( nabor < blkbeg ) then

            if ( maxnz < nzcnt ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'FNOFNZ - Fatal error!'
              write ( *, '(a)' ) '  NZCNT exceeds MAXNZ'
              write ( *, '(a,i8)' ) '  NZCNT = ', nzcnt
              write ( *, '(a,i8)' ) '  MAXNZ = ', maxnz
              stop
            end if

            nzsubs(nzcnt) = nabor
            nzcnt = nzcnt + 1

          end if

        end do
!
!  Sort the subscripts of row J.
!
        jxnonz = xnonz(j)
        if ( jxnonz < nzcnt ) then
          call i4vec_sort_insert_a ( nzcnt-jxnonz, nzsubs(jxnonz) )
        end if

      end if

    end do

  end do

  xnonz(blkend+1) = nzcnt

  maxnz = nzcnt - 1

  return
end
subroutine fnspan ( adj_row, adj, nodlvl, nspan, set, level, nadjs, adjs, &
  leaf, n )

!*****************************************************************************80
!
!! FNSPAN finds the span of a level subgraph subset, as part of the RQT method.
!
!  Discussion:
!
!    The adjacent set of the span in the lower level is also determined.
!    If the span has an unnumbered node in the higher level, an unnumbered
!    leaf node (i.e. one with no neighbor in next level) will be returned.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) LEVEL, the level number of the current set.
!
!    Input/output, integer ( kind = 4 ) NSPAN, SET(*), the input set.  On return,
!    the resulting span set.
!
!    Input/output, integer ( kind = 4 ) NODLVL(N), the level number vector.  Nodes
!    considered will have their NODLVL changed to zero.
!
!    Output, integer ( kind = 4 ) NADJS, ADJS(*), the adjacency set of the span
!    in the lower level.
!
!    Output, integer ( kind = 4 ) LEAF.  If the span has an unnumbered higher level node,
!    LEAF returns an unnumbered leaf node in the level
!    structure, otherwise, LEAF is zero.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) adjs(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) leaf
  integer ( kind = 4 ) level
  integer ( kind = 4 ) lvl
  integer ( kind = 4 ) lvlm1
  integer ( kind = 4 ) nadjs
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nbrlvl
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nodlvl(n)
  integer ( kind = 4 ) nspan
  integer ( kind = 4 ) set(*)
  integer ( kind = 4 ) setptr
!
!  Initialization.
!
  leaf = 0
  nadjs = 0
  setptr = 0

   10 continue

  setptr = setptr + 1

  if ( nspan < setptr ) then
    return
  end if
!
!  For each node in the partially spanned set...
!
  node = set(setptr)
  jstrt = adj_row(node)
  jstop = adj_row(node+1) - 1

  if ( jstop < jstrt ) then
    go to 10
  end if
!
!  For each neighbor of node, test its NODLVL value.
!
  do j = jstrt, jstop

    nbr = adj(j)
    nbrlvl = nodlvl(nbr)

    if ( 0 < nbrlvl ) then
!
!  Node NBR is in the level.  Add it to the span set.
!
      if ( level == nbrlvl ) then

        nspan = nspan + 1
        set(nspan) = nbr

      else if ( level < nbrlvl ) then

        go to 60
!
!  Node NBR is in LEVEL - 1.  Add it to ADJS.
!
      else

        nadjs = nadjs + 1
        adjs(nadjs) = nbr

      end if

      nodlvl(nbr) = 0

    end if

  end do

  go to 10
!
!  Node NBR is in LEVEL + 1.
!
!  Find an unnumbered leaf node by tracing a path up the level structure.
!
!  Then reset the NODLVL values of nodes in adjs.
!
   60 continue

  leaf = nbr
  lvl = level + 1

   70 continue

  jstrt = adj_row(leaf)
  jstop = adj_row(leaf+1) - 1

  do j = jstrt, jstop

    nbr = adj(j)

    if ( lvl < nodlvl(nbr) ) then
      leaf = nbr
      lvl = lvl + 1
      go to 70
    end if

  end do

  if ( nadjs <= 0 ) then
    return
  end if

  lvlm1 = level - 1

  do i = 1, nadjs
    node = adjs(i)
    nodlvl(node) = lvlm1
  end do

  return
end
subroutine fntadj ( adj_row, adj, perm, nblks, xblk, father, n )

!*****************************************************************************80
!
!! FNTADJ determines the quotient tree adjacency structure for a graph.
!
!  Discussion:
!
!    The graph is defined by the matrix that is to be factored.
!
!    The structure is represented by the father vector.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) PERM(N), the permutation vector.
!
!    Input, integer ( kind = 4 ) NBLKS, XBLK(NBLKS+1), the tree partitioning.
!
!    Output, integer ( kind = 4 ) FATHER(N), the father vector of the quotient tree.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) nblks
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) father(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) istrt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) nabor
  integer ( kind = 4 ) nbrblk
  integer ( kind = 4 ) node
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) xblk(nblks+1)
!
!  Initialize the block number vector.
!
  do k = 1, nblks

    do i = xblk(k), xblk(k+1) - 1
      node = perm(i)
      mask(node) = k
    end do

  end do
!
!  For each block, find its father block in the tree structure.
!
  father(nblks) = 0

  do k = 1, nblks - 1

    istrt = xblk(k)
    istop = xblk(k+1) - 1
    father(k) = 0

    do i = istrt, istop

      node = perm(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1

      do j = jstrt, jstop

        nabor = adj(j)
        nbrblk = mask(nabor)

        if ( k < nbrblk ) then
          father(k) = nbrblk
          go to 10
        end if

      end do

    end do

10  continue

  end do

  return
end
subroutine fntenv ( adj_row, adj, perm, perm_inv, nblks, xblk, xenv, &
  env_size, n )

!*****************************************************************************80
!
!! FNTENV determines the envelope index vector.
!
!  Discussion:
!
!    FNTENV determines the envelope index vector for the envelope of the
!    diagonal blocks of a tree partitioned system.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) PERM(N), the permutation vector.
!
!    Input, integer ( kind = 4 ) PERM_INV(N), the inverse permutation vector.
!
!    Input, integer ( kind = 4 ) NBLKS, XBLK(NBLKS+1), the tree partitioning.
!
!    Output, integer ( kind = 4 ) XENV(N+1), the envelope index vector.
!
!    Output, integer ( kind = 4 ) ENV_SIZE, the size of the envelope found.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) nblks
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) blkbeg
  integer ( kind = 4 ) blkend
  integer ( kind = 4 ) i
  integer ( kind = 4 ) env_size
  integer ( kind = 4 ) ifirst
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kfirst
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) node
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  integer ( kind = 4 ) xblk(nblks+1)
  integer ( kind = 4 ) xenv(n+1)

  env_size = 1
!
!  Loop through each block in the partitioning.
!
  do k = 1, nblks

    blkbeg = xblk(k)
    blkend = xblk(k+1) - 1
!
!  KFIRST stores the first node in the K-th block that has a neighbor
!  in the previous blocks.
!
    kfirst = blkend

    do i = blkbeg, blkend

      xenv(i) = env_size
      node = perm(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1
!
!  IFIRST stores the first nonzero in the I-th row within
!  the K-th block.
!
      ifirst = i

      do j = jstrt, jstop

        nbr = adj(j)
        nbr = perm_inv(nbr)

        if ( blkbeg <= nbr ) then
          ifirst = min ( ifirst, nbr )
        else
          ifirst = min ( ifirst, kfirst )
          kfirst = min ( kfirst, i )
        end if

      end do

      env_size = env_size + i - ifirst

    end do

  end do

  xenv(blkend+1) = env_size
  env_size = env_size - 1

  return
end
subroutine fn1wd ( root, adj_row, adj, mask, nsep, sep, nlvl, xls, ls, n )

!*****************************************************************************80
!
!! FN1WD finds one-way dissectors of a connected component for the 1WD method.
!
!  Discussion:
!
!    The connected component is specified by MASK and ROOT.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the node that defines (along with MASK) the
!    component to be processed.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) NSEP, the number of nodes in the one-way dissectors.
!
!    Output, integer ( kind = 4 ) SEP(*), the dissector nodes.
!
!    Input/output, integer ( kind = 4 ) MASK(N), nodes in the dissector have their
!    mask values set to zero.
!
!    Workspace, integer ( kind = 4 ) XLS(N+1), LS(N), the level structure used by
!    the routine ROOT_FIND.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  real ( kind = 8 ) deltp1
  real ( kind = 8 ) fnlvl
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kstop
  integer ( kind = 4 ) kstrt
  integer ( kind = 4 ) lp1beg
  integer ( kind = 4 ) lp1end
  integer ( kind = 4 ) ls(n)
  integer ( kind = 4 ) lvl
  integer ( kind = 4 ) lvlbeg
  integer ( kind = 4 ) lvlend
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nlvl
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nsep
  integer ( kind = 4 ) root
  integer ( kind = 4 ) sep(*)
  real ( kind = 8 ) width
  integer ( kind = 4 ) xls(n+1)

  call root_find ( root, adj_row, adj, mask, nlvl, xls, ls, n )

  fnlvl = real ( nlvl, kind = 8 )
  nsep = xls(nlvl+1) - 1
  width = real ( nsep, kind = 8 ) / fnlvl
  deltp1 = 1.0D+00 + sqrt ( ( 3.0D+00 * width + 13.0D+00 ) / 2.0D+00 )
!
!  The component is too small, or the level structure
!  is very long and narrow.  Return the whole component.
!
  if ( nsep < 50 .or. 0.5D+00 * fnlvl < deltp1 ) then

    do i = 1, nsep
      node = ls(i)
      sep(i) = node
      mask(node) = 0
    end do

    return

  end if
!
!  Find the parallel dissectors.
!
  nsep = 0
  i = 0

  do

    i = i + 1
    lvl = int ( real ( i, kind = 8 ) * deltp1 + 0.5D+00 )

    if ( nlvl <= lvl ) then
      exit
    end if

    lvlbeg = xls(lvl)
    lp1beg = xls(lvl+1)
    lvlend = lp1beg - 1
    lp1end = xls(lvl+2) - 1

    do j = lp1beg, lp1end
      node = ls(j)
      adj_row(node) = -adj_row(node)
    end do
!
!  Nodes in level LVL are chosen to form dissector.
!  Include only those with neighbors in LVL + 1 level.
!  ADJ_ROW is used temporarily to mark nodes in LVL + 1.
!
    do j = lvlbeg, lvlend

      node = ls(j)
      kstrt = adj_row(node)
      kstop = abs ( adj_row(node+1) ) - 1

      do k = kstrt, kstop

        nbr = adj(k)

        if ( adj_row(nbr) <= 0 ) then
          nsep = nsep + 1
          sep(nsep) = node
          mask(node) = 0
          exit
        end if

      end do

    end do

    do j = lp1beg, lp1end
      node = ls(j)
      adj_row(node) = -adj_row(node)
    end do

  end do

  return
end
subroutine gennd ( n, adj_row, adj, perm )

!*****************************************************************************80
!
!! GENND finds a nested dissection ordering for a general graph.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) PERM(N), the nested dissection ordering.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) MASK(N), masks off variables that have
!    been numbered during the orderng process.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) nsep
  integer ( kind = 4 ) num
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) root

  mask(1:n) = 1

  num = 0

  do i = 1, n
!
!  For each masked component, find a separator and number the nodes next.
!
    if ( mask(i) /= 0 ) then

      root = i

      call fndsep ( root, adj_row, adj, mask, nsep, perm(num+1), n )

      num = num + nsep

      if ( n <= num ) then
        call i4vec_reverse ( n, perm )
        return
      end if

    end if

  end do
!
!  Separators found first should be ordered last.
!
  call i4vec_reverse ( n, perm )

  return
end
subroutine genqmd ( n, adj_row, adj, perm, perm_inv, marker, rchset, nbrhd, &
  qsize, qlink, nofsub )

!*****************************************************************************80
!
!! GENQMD implements the quotient minimum degree algorithm.
!
!  Discussion:
!
!    The routine uses the implicit representation of the elimination
!    graphs by quotient graphs, and the notion of indistinguishable nodes.
!
!    The adjacency vector ADJ will be destroyed.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) PERM(N), the minimum degree ordering.
!
!    Output, integer ( kind = 4 ) PERM_INV(N), the inverse of PERM.
!
!    Workspace, integer ( kind = 4 ) MARKER(N), a marker vector, where MARKER(I) is
!    negative means node I has been merged with another node and thus
!    can be ignored.
!
!    Workspace, integer ( kind = 4 ) RCHSET(*), used for the reachable set.
!
!    Workspace, integer ( kind = 4 ) NBRHD(*), a vector used for the neighborhood set.
!
!    Workspace, integer ( kind = 4 ) QSIZE(N), the sizes of indistinguishable supernodes.
!
!    Workspace, integer ( kind = 4 ) QLINK(N), vector to store indistinguishable nodes,
!    I, QLINK(I), QLINK(QLINK(I)) are the members of the supernode
!    represented by I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) deg(n)
  integer ( kind = 4 ) inode
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) irch
  integer ( kind = 4 ) j
  integer ( kind = 4 ) marker(n)
  integer ( kind = 4 ) mindeg
  integer ( kind = 4 ) nbrhd(*)
  integer ( kind = 4 ) ndeg
  integer ( kind = 4 ) nhdsze
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nofsub
  integer ( kind = 4 ) np
  integer ( kind = 4 ) num
  integer ( kind = 4 ) nxnode
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  integer ( kind = 4 ) qlink(n)
  integer ( kind = 4 ) qsize(n)
  integer ( kind = 4 ) rchset(*)
  integer ( kind = 4 ) rchsze
  integer ( kind = 4 ) search
  integer ( kind = 4 ) thresh
!
!  Initialize degree vector and other working variables.
!
!  DEG(I) negative means node I has been numbered.
!
  mindeg = n
  nofsub = 0

  do node = 1, n
    perm(node) = node
  end do

  do node = 1, n
    perm_inv(node) = node
  end do

  marker(1:n) = 0
  qsize(1:n) = 1
  qlink(1:n) = 0

  do node = 1, n
    ndeg = adj_row(node+1) - adj_row(node)
    deg(node) = ndeg
    mindeg = min ( mindeg, ndeg )
  end do

  num = 0
!
!  Perform threshold search to get a node of minimum degree.
!  The variable SEARCH points to where the search should start.
!
  search = 1
  thresh = mindeg
  mindeg = n

  do

    do

      search = max ( search, num + 1 )

      do j = search, n

        node = perm(j)

        if ( 0 <= marker(node) ) then
          ndeg = deg(node)
          if ( ndeg <= thresh ) then
            go to 50
          end if
          mindeg = min ( mindeg, ndeg )
        end if

      end do

      search = 1
      thresh = mindeg
      mindeg = n

    end do
!
!  The node has minimum degree.
!  Find its reachable sets by calling QMDRCH.
!
50  continue

    search = j
    nofsub = nofsub + deg(node)
    marker(node) = 1

    call qmdrch ( node, adj_row, adj, deg, marker, rchsze, rchset, nhdsze, &
      nbrhd, n )
!
!  Eliminate all nodes indistinguishable from node.
!  They are given by node, qlink(node),
!
    nxnode = node

    do

      num = num + 1
      np = perm_inv(nxnode)
      ip = perm(num)
      perm(np) = ip
      perm_inv(ip) = np
      perm(num) = nxnode
      perm_inv(nxnode) = num
      deg(nxnode) = -1
      nxnode = qlink(nxnode)

      if ( nxnode <= 0 ) then
        exit
      end if

    end do

    if ( rchsze <= 0 ) then
      go to 80
    end if
!
!  Update the degrees of the nodes in the reachable
!  set and identify indistinguishable nodes.
!
    call qmdupd ( adj_row, adj, rchsze, rchset, deg, qsize, qlink, marker, &
      rchset(rchsze+1), nbrhd(nhdsze+1), n )
!
!  Reset marker value of nodes in reach set.
!  Update threshold value for cyclic search.
!  Call QMDQT to form new quotient graph.
!
    marker(node) = 0

    do irch = 1, rchsze

      inode = rchset(irch)

      if ( 0 <= marker(inode) ) then

        marker(inode) = 0
        ndeg = deg(inode)
        mindeg = min ( mindeg, ndeg )

        if ( ndeg <= thresh ) then
          mindeg = thresh
          thresh = ndeg
          search = perm_inv(inode)
        end if

      end if

    end do

    if ( 0 < nhdsze ) then
      call qmdqt ( node, adj_row, adj, marker, rchsze, rchset, nbrhd, n )
    end if

80  continue

    if ( n <= num ) then
      exit
    end if

  end do

  return
end
subroutine genrcm ( n, adj_row, adj, perm )

!*****************************************************************************80
!
!! GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
!
!  Discussion:
!
!    For each connected component in the graph, the routine obtains
!    an ordering by calling RCM.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) PERM(N), the RCM ordering.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) XLS(N+1), the index vector for a level structure.
!    The level structure is stored in the currently unused spaces in the
!    permutation vector PERM.
!
!    Local, integer ( kind = 4 ) MASK(N), marks variables that have been numbered.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ccsize
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) nlvl
  integer ( kind = 4 ) num
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) root
  integer ( kind = 4 ) xls(n+1)

  mask(1:n) = 1

  num = 1

  do i = 1, n
!
!  For each masked connected component...
!
    if ( mask(i) /= 0 ) then

      root = i
!
!  Find a pseudo-peripheral node ROOT.  The level structure found by
!  ROOT_FIND is stored starting at PERM(NUM).
!
      call root_find ( root, adj_row, adj, mask, nlvl, xls, perm(num), n )
!
!  RCM orders the component using ROOT as the starting node.
!
      call rcm ( root, adj_row, adj, mask, perm(num), ccsize, n )

      num = num + ccsize

      if ( n < num ) then
        return
      end if

    end if

  end do

  return
end
subroutine genrqt ( n, adj_row, adj, nblks, xblk, perm )

!*****************************************************************************80
!
!! GENRQT determines a partitioned ordering using refined quotient trees.
!
!  Discussion:
!
!    The graph being reordered may be disconnected.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) NBLKS, XBLK(*), the quotient tree partitioning.
!
!    Output, integer ( kind = 4 ) PERM(N), the permutation vector.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) XLS(N+1), LS(N), a level structure pair
!    used to find a pseudo-peripheral node.
!
!    Local, integer ( kind = 4 ) NODLVL(N), the level number of each node in a
!    level structure.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ixls
  integer ( kind = 4 ) leaf
  integer ( kind = 4 ) ls(n)
  integer ( kind = 4 ) nblks
  integer ( kind = 4 ) nlvl
  integer ( kind = 4 ) nodlvl(n)
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) root
  integer ( kind = 4 ) xblk(*)
  integer ( kind = 4 ) xls(n+1)

  nodlvl(1:n) = 1
  nblks = 0
  xblk(1) = 1

  do i = 1, n
!
!  For each connected component...
!
    if ( 0 < nodlvl(i) ) then
!
!  Find a rooted level structure.
!
      root = i
      call fnlvls ( root, adj_row, adj, nodlvl, nlvl, xls, ls, n )
      ixls = xls(nlvl)
      leaf = ls(ixls)
!
!  RQTREE gets the block order.
!
      call rqtree ( leaf, adj_row, adj, perm, nblks, xblk, nodlvl, xls, ls, n )

    end if

  end do

  return
end
subroutine gen1wd ( n, adj_row, adj, nblks, xblk, perm )

!*****************************************************************************80
!
!! GEN1WD partitions a general graph, for the 1WD method.
!
!  Discussion:
!
!    FN1WD is called once for each connected component.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) NBLKS, XBLK(*), the partitioning found.
!
!    Output, integer ( kind = 4 ) PERM(N), the one-way dissection ordering.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) XLS(N+1), LS(N), a level structure used
!    by LEVEL_SET.
!
!    Local, integer ( kind = 4 ) MASK(N), marks variables that have been numbered.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ccsize
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lnum
  integer ( kind = 4 ) ls(n)
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) nblks
  integer ( kind = 4 ) nlvl
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nsep
  integer ( kind = 4 ) num
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) root
  integer ( kind = 4 ) xblk(*)
  integer ( kind = 4 ) xls(n+1)

  mask(1:n) = 1

  nblks = 0
  num = 0

  do i = 1, n

    if ( mask(i) /= 0 ) then
!
!  Find a one-way dissector for each component.
!
      root = i

      call fn1wd ( root, adj_row, adj, mask, nsep, perm(num+1), nlvl, &
        xls, ls, n )

      num = num + nsep
      nblks = nblks + 1
      xblk(nblks) = n - num + 1
      ccsize = xls(nlvl+1) - 1
!
!  Number the remaining nodes in the component.
!  Each component in the remaining subgraph forms
!  a new block in the partitioning.
!
      do j = 1, ccsize

        node = ls(j)

        if ( mask(node) /= 0 ) then

          call level_set ( node, adj_row, adj, mask, nlvl, xls, perm(num+1), n )

          lnum = num + 1
          num = num + xls(nlvl+1) - 1
          nblks = nblks + 1
          xblk(nblks) = n - num + 1

          do k = lnum, num
            node = perm(k)
            mask(node) = 0
          end do

          if ( n < num ) then
            go to 50
          end if

        end if

      end do

    end if

  end do
!
!  Since dissectors found first should be ordered last,
!  REVRSE is called to adjust the ordering
!  vector and the block index vector.
!
   50 continue

  call i4vec_reverse ( n, perm )

  call i4vec_reverse ( nblks, xblk )

  xblk(nblks+1) = n + 1

  return
end
subroutine gs_factor ( n, ixlnz, lnz, xnzsub, nzsub, diag )

!*****************************************************************************80
!
!! GS_FACTOR: symmetric factorization for a general sparse system.
!
!  Discussion:
!
!    The matrix is stored in the compressed subscript data format.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IXLNZ(N+1), the index vector for LNZ.  IZLNZ(I)
!    points to the start of the nonzeros in column I of factor L.
!
!    Input/output, real ( kind = 8 ) LNZ(*), on input, the nonzeros of A,
!    and on output, the nonzeros of L.
!
!    Input, integer ( kind = 4 ) XNZSUB(*), NZSUB(*), the compressed subscript data
!    structure for the factor L.
!
!    Input/output, real ( kind = 8 ) DIAG(N), on input, the diagonal of A,
!    and on output, the diagonal of L.
!
!  Local parameters:
!
!    integer ( kind = 4 ) FIRST(N), temporary vector to point to the first
!    nonzero in each column that will be used
!    next for modification.
!
!    Workspace, integer ( kind = 4 ) LINK(N), at step J, the list in LINK(J),
!    LINK(LINK(J)), ... consists of those columns that will modify
!    the column L(*,J).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) diagj
  integer ( kind = 4 ) first(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) istrt
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) ixlnz(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kfirst
  integer ( kind = 4 ) link(n)
  real ( kind = 8 ) ljk
  real ( kind = 8 ) lnz(*)
  integer ( kind = 4 ) newk
  integer ( kind = 4 ) nzsub(*)
  real ( kind = 8 ) temp(n)
  integer ( kind = 4 ) xnzsub(*)
!
!  Initialize the working vectors.
!
  link(1:n) = 0
  temp(1:n) = 0.0D+00
!
!  Compute column L(*,j) for j = 1,....,N.
!
  do j = 1, n
!
!  For each column L(*,K) that affects L(*,J).
!
    diagj = 0.0D+00
    newk = link(j)

    do

      k = newk
      if ( k == 0 ) then
        exit
      end if

      newk = link(k)
!
!  Outer product modification of L(*,J) by L(*,K) starting at FIRST(K)
!  of L(*,K).
!
      kfirst = first(k)
      ljk = lnz(kfirst)
      diagj = diagj + ljk**2
      istrt = kfirst + 1
      istop = ixlnz(k+1) - 1

      if ( istop < istrt ) then
        cycle
      end if
!
!  Before modification, update vectors first,
!  and link for future modification steps.
!
      first(k) = istrt
      i = xnzsub(k) + ( kfirst - ixlnz(k) ) + 1
      isub = nzsub(i)
      link(k) = link(isub)
      link(isub) = k
!
!  The actual modification is saved in vector TEMP.
!
      do ii = istrt, istop
        isub = nzsub(i)
        temp(isub) = temp(isub) + lnz(ii) * ljk
        i = i + 1
      end do

    end do
!
!  Apply the modifications accumulated in TEMP to column L(*,J).
!
    diagj = diag(j) - diagj

    if ( diagj <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GS_FACTOR - Fatal error'
      write ( *, '(a)' ) '  Zero or negative diagonal entry!'
      write ( *, '(a,g14.6)' ) '  DIAG(J) = ', diagj
      write ( *, '(a,i8)' ) '  for diagonal J = ', j
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix is not positive definite.'
      stop
    end if

    diagj = sqrt ( diagj )
    diag(j) = diagj
    istrt = ixlnz(j)
    istop = ixlnz(j+1) - 1

    if ( istrt <= istop ) then

      first(j) = istrt
      i = xnzsub(j)
      isub = nzsub(i)
      link(j) = link(isub)
      link(isub) = j

      do ii = istrt, istop
        isub = nzsub(i)
        lnz(ii) = ( lnz(ii) - temp(isub) ) / diagj
        temp(isub) = 0.0D+00
        i = i + 1
      end do

    end if

  end do

  return
end
subroutine gs_solve ( n, ixlnz, lnz, xnzsub, nzsub, diag, rhs )

!*****************************************************************************80
!
!! GS_SOLVE solves a factored general sparse system.
!
!  Discussion:
!
!    The matrix is stored in compressed subscript format.
!
!    The matrix must have been factored by GS_FACTOR.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IXLNZ(N+1), real ( kind = 8 ) LNZ(*), the structure of
!    nonzeros in L.
!
!    Input, integer ( kind = 4 ) XNZSUB(*), NZSUB(*), the compressed subscript structure.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal components of L.
!
!    Input/output, real ( kind = 8 ) RHS(N), on input, the right hand side
!    of the linear system; on output, the solution vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) istrt
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) ixlnz(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  real ( kind = 8 ) lnz(*)
  integer ( kind = 4 ) nzsub(*)
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) rhsj
  real ( kind = 8 ) s
  integer ( kind = 4 ) xnzsub(*)
!
!  Forward substitution.
!
  do j = 1, n

    rhsj = rhs(j) / diag(j)
    rhs(j) = rhsj
    istrt = ixlnz(j)
    istop = ixlnz(j+1) - 1

    i = xnzsub(j)

    do ii = istrt, istop
      isub = nzsub(i)
      rhs(isub) = rhs(isub) - lnz(ii) * rhsj
      i = i + 1
    end do

  end do
!
!  Backward substitution.
!
  j = n

  do jj = 1, n

    s = rhs(j)
    istrt = ixlnz(j)
    istop = ixlnz(j+1) - 1

    i = xnzsub(j)

    do ii = istrt, istop
      isub = nzsub(i)
      s = s - lnz(ii) * rhs(isub)
      i = i + 1
    end do

    rhs(j) = s / diag(j)
    j = j - 1

  end do

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP switches two integer ( kind = 4 ) values.
!
!  Modified:
!
!    24 February 2007
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
subroutine i4vec_copy ( n, a, b )

!*****************************************************************************80
!
!! I4VEC_COPY copies one integer ( kind = 4 ) vector into another.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be copied.
!
!    Output, integer ( kind = 4 ) B(N), a copy of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)

  b(1:n) = a(1:n)

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an integer ( kind = 4 ) vector to the indicator vector.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! I4VEC_REVERSE reverses the elements of an integer ( kind = 4 ) vector.
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n/2
    call i4_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine i4vec_sort_insert_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_INSERT_A uses an ascending insertion sort on an integer ( kind = 4 ) vector.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, data to be sorted.
!    On output, the entries of A have been sorted in ascending order.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( a(j) <= x ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine level_set ( root, adj_row, adj, mask, nlvl, xls, ls, n )

!*****************************************************************************80
!
!! LEVEL_SET generates the connected level structure rooted at a given node.
!
!  Discussion:
!
!    Only nodes for which MASK is nonzero will be considered.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the node at which the level structure
!    is to be rooted.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input/output, integer ( kind = 4 ) MASK(N).  On input, only nodes with nonzero
!    MASK are to be processed.  On output, those nodes which were included
!    in the level set have MASK set to 1.
!
!    Output, integer ( kind = 4 ) NLVL, the number of levels in the level
!    structure.  ROOT is in level 1.  The neighbors of ROOT
!    are in level 2, and so on.
!
!    Output, integer ( kind = 4 ) XLS(N+1), LS(N), the rooted level structure.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ccsize
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) lbegin
  integer ( kind = 4 ) ls(n)
  integer ( kind = 4 ) lvlend
  integer ( kind = 4 ) lvsize
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nlvl
  integer ( kind = 4 ) node
  integer ( kind = 4 ) root
  integer ( kind = 4 ) xls(n+1)

  mask(root) = 0
  ls(1) = root
  nlvl = 0
  lvlend = 0
  ccsize = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
  do

    lbegin = lvlend + 1
    lvlend = ccsize
    nlvl = nlvl + 1
    xls(nlvl) = lbegin
!
!  Generate the next level by finding all the masked neighbors of nodes
!  in the current level.
!
    do i = lbegin, lvlend

      node = ls(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then
          ccsize = ccsize + 1
          ls(ccsize) = nbr
          mask(nbr) = 0
        end if

      end do

    end do
!
!  Compute the current level width.
!  If it is positive, generate the next level.
!
    lvsize = ccsize - lvlend

    if ( lvsize <= 0 ) then
      exit
    end if

  end do
!
!  Reset MASK to one for the nodes in the level structure.
!
  xls(nlvl+1) = lvlend + 1

  do i = 1, ccsize
    node = ls(i)
    mask(node) = 1
  end do

  return
end
subroutine perm_inverse ( n, perm, perm_inv )

!*****************************************************************************80
!
!! PERM_INVERSE produces the inverse of a given permutation.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items permuted.
!
!    Input, integer ( kind = 4 ) PERM(N), a permutation.
!
!    Output, integer ( kind = 4 ) PERM_INV(N), the inverse permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)

  do i = 1, n
    perm_inv(perm(i)) = i
  end do

  return
end
subroutine perm_rv ( n, rhs, perm )

!*****************************************************************************80
!
!! PERM_RV undoes the permutation of the right hand side.
!
!  Discussion:
!
!    This routine should be called once the linear system has been solved and
!    the solution returned in RHS.  The routine then undoes the permutation
!    of RHS, restoring the original ordering.  To do this, it needs the
!    PERM vector which defined the reordering used by the solver.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input/output, real ( kind = 8 ) RHS(N).
!    On input, the solution of the permuted linear system.
!    On output, the solution of the original linear system.
!
!    Input, integer ( kind = 4 ) PERM(N), the permutation information.
!    PERM(I) = K means that the K-th equation and variable in the
!    original ordering became the I-th equation and variable in the
!    reordering.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) perm(n)
  real ( kind = 8 ) pull
  real ( kind = 8 ) put
  real ( kind = 8 ) rhs(n)
!
!  Mark PERM with negative signs which will be removed
!  as each permuted element is restored to its rightful place
!
  perm(1:n) = -perm(1:n)
!
!  Search for the next element of PERM which is the first
!  element of a permutation cycle.
!
  istart = 0

20 continue

  do

    istart = istart + 1

    if ( n < istart ) then
      return
    end if

    if ( 0 < perm(istart) ) then
      cycle
    end if

    if ( abs ( perm(istart) ) /= istart ) then
      exit
    end if

    perm(istart) = abs ( perm(istart) )

  end do
!
!  Begin a cycle.
!
  perm(istart) = abs ( perm(istart) )
  iput = istart
  pull = rhs(iput)

  do

    iput = abs ( perm(iput) )
    put = rhs(iput)
    rhs(iput) = pull
    pull = put

    if ( 0 < perm(iput) ) then
      go to 20
    end if

    perm(iput) = abs ( perm(iput) )

  end do

end
subroutine qmdmrg ( adj_row, adj, deg, qsize, qlink, marker, deg0, nhdsze, &
  nbrhd, rchset, ovrlp, n )

!*****************************************************************************80
!
!! QMDMRG merges indistinguishable nodes for the QMD method.
!
!  Discussion:
!
!    QMDMRG merges indistinguishable nodes in the minimum degree ordering
!    algorithm, and computes the new degrees of these new supernodes.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) DEG0, the number of nodes in the given set.
!
!    Input, integer ( kind = 4 ) NHDSZE, NBHRD(*), the set of eliminated supernodes
!    adjacent to some nodes in the set.
!
!    Input/output, integer ( kind = 4 ) DEG(N), the degree vector.
!
!    Output, integer ( kind = 4 ) QSIZE(N), the size of indistinguishable nodes.
!
!    Output, integer ( kind = 4 ) QLINK(N), the linked list for indistinguishable nodes.
!
!    Output, integer ( kind = 4 ) MARKER(N), the given set is given by
!    those nodes with marker value set to 1.  Those nodes with degree
!    updated will have marker value set to 2.
!
!    Workspace, integer ( kind = 4 ) RCHSET(*), the reachable set.
!
!    Workspace, integer ( kind = 4 ) OVRLP(*), stores the intersection
!    of two reachable sets.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) deg(n)
  integer ( kind = 4 ) deg0
  integer ( kind = 4 ) deg1
  integer ( kind = 4 ) head
  integer ( kind = 4 ) inhd
  integer ( kind = 4 ) iov
  integer ( kind = 4 ) irch
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) link
  integer ( kind = 4 ) lnode
  integer ( kind = 4 ) mark
  integer ( kind = 4 ) marker(n)
  integer ( kind = 4 ) mrgsze
  integer ( kind = 4 ) nabor
  integer ( kind = 4 ) nbrhd(*)
  integer ( kind = 4 ) nhdsze
  integer ( kind = 4 ) node
  integer ( kind = 4 ) novrlp
  integer ( kind = 4 ) ovrlp(*)
  integer ( kind = 4 ) qlink(n)
  integer ( kind = 4 ) qsize(n)
  integer ( kind = 4 ) rchset(*)
  integer ( kind = 4 ) rchsze
  integer ( kind = 4 ) root

  if ( nhdsze <= 0 ) then
    return
  end if

  do inhd = 1, nhdsze
    root = nbrhd(inhd)
    marker(root) = 0
  end do
!
!  Loop through each eliminated supernode in the set (NHDSZE,NBRHD).
!
  do inhd = 1, nhdsze

    root = nbrhd(inhd)
    marker(root) = - 1
    rchsze = 0
    novrlp = 0
    deg1 = 0

   20   continue

    jstrt = adj_row(root)
    jstop = adj_row(root+1) - 1
!
!  Determine the reachable set and its intersection with the input
!  reachable set.
!
    do j = jstrt, jstop

      nabor = adj(j)
      root = -nabor

      if ( nabor < 0 ) then
        go to 20
      end if

      if ( nabor == 0 ) then
        exit
      end if

      mark = marker(nabor)

      if ( mark == 0 ) then

        rchsze = rchsze + 1
        rchset(rchsze) = nabor
        deg1 = deg1+qsize(nabor)
        marker(nabor) = 1

      else if ( mark == 1 ) then

        novrlp = novrlp + 1
        ovrlp(novrlp) = nabor
        marker(nabor) = 2

      end if

    end do
!
!  From the overlapped set, determine the nodes that can be merged.
!
    head = 0
    mrgsze = 0

    do iov = 1, novrlp

      node = ovrlp(iov)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1

      do j = jstrt, jstop

        nabor = adj(j)

        if ( marker(nabor) == 0 ) then
          marker(node) = 1
          go to 110
        end if

      end do
!
!  NODE belongs to the new merged supernode.
!  Update the vectors QLINK and QSIZE.
!
      mrgsze = mrgsze + qsize(node)
      marker(node) = -1
      lnode = node

      do

        link = qlink(lnode)

        if ( link <= 0 ) then
          exit
        end if

        lnode = link

      end do

      qlink(lnode) = head
      head = node

  110     continue

    end do

    if ( 0 < head ) then
      qsize(head) = mrgsze
      deg(head) = deg0 + deg1 - 1
      marker(head) = 2
    end if
!
!  Reset marker values.
!
    root = nbrhd(inhd)
    marker(root) = 0

    do irch = 1, rchsze
      node = rchset(irch)
      marker(node) = 0
    end do

  end do

  return
end
subroutine qmdqt ( root, adj_row, adj, marker, rchsze, rchset, nbrhd, n )

!*****************************************************************************80
!
!! QMDQT performs the quotient graph transformation.
!
!  Discussion:
!
!    QMDQT is called after a node has been eliminated.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the node just eliminated.  It becomes the
!    representative of the new supernode.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) RCHSZE, RCHSET(*), the reachable set of root in the
!    old quotient graph.
!
!    Input, integer ( kind = 4 ) NBRHD(*), the neighborhood set which will be merged
!    with root to form the new supernode.
!
!    Input, integer ( kind = 4 ) MARKER(N), the marker vector.
!
!    Input/output, integer ( kind = 4 ) ADJ(*), becomes the adjacency structure of
!    the quotient graph.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) inhd
  integer ( kind = 4 ) irch
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) link
  integer ( kind = 4 ) marker(n)
  integer ( kind = 4 ) nabor
  integer ( kind = 4 ) nbrhd(*)
  integer ( kind = 4 ) node
  integer ( kind = 4 ) rchset(*)
  integer ( kind = 4 ) rchsze
  integer ( kind = 4 ) root

  irch = 0
  inhd = 0
  node = root

  do

    jstrt = adj_row(node)
    jstop = adj_row(node+1) - 2
!
!  Place reach nodes into the adjacent list of node.
!
    do j = jstrt, jstop
      irch = irch + 1
      adj(j) = rchset(irch)
      if ( rchsze <= irch ) then
        go to 40
      end if
    end do
!
!  Link to other space provided by the nbrhd set.
!
    link = adj(jstop+1)
    node = -link

    if ( 0 <= link ) then
      inhd = inhd + 1
      node = nbrhd(inhd)
      adj(jstop+1) = -node
    end if

  end do
!
!  All reachable nodes have been saved.  End the adjacency list.
!  Add root to the neighbor list of each node in the reach set.
!
40 continue

  adj(j+1) = 0

  do irch = 1, rchsze

    node = rchset(irch)

    if ( 0 <= marker(node) ) then

      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1

      do j = jstrt, jstop

        nabor = adj(j)

        if ( marker(nabor) < 0 ) then
          adj(j) = root
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine qmdrch ( root, adj_row, adj, deg, marker, rchsze, rchset, nhdsze, &
  nbrhd, n )

!*****************************************************************************80
!
!! QMDRCH determines the reachable set of a node, for the QMD method.
!
!  Discussion:
!
!    The adjacency structure is assumed to be stored in a quotient graph format.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the given node, which is not in the subset.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) DEG(N), the degree vector.  DEG(I) < 0 means node I
!    belongs to the given subset.
!
!    Output, integer ( kind = 4 ) RCHSZE, RCHSET(*), the reachable set.
!
!    Output, integer ( kind = 4 ) NHDSZE, NBRHD(*), the neighborhood set.
!
!    Input/output, integer ( kind = 4 ) MARKER(N), the marker vector for the reachable
!    and neighborhood sets.
!    > 0 means the node is in reach set.
!    < 0 means the node has been merged with
!    others in the quotient or it is in nbrhd set.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) deg(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) istrt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) marker(n)
  integer ( kind = 4 ) nabor
  integer ( kind = 4 ) nbrhd(*)
  integer ( kind = 4 ) nhdsze
  integer ( kind = 4 ) node
  integer ( kind = 4 ) rchset(*)
  integer ( kind = 4 ) rchsze
  integer ( kind = 4 ) root
!
!  Loop through the neighbors of ROOT in the quotient graph.
!
  nhdsze = 0
  rchsze = 0
  istrt = adj_row(root)
  istop = adj_row(root+1) - 1

  do i = istrt, istop

    nabor = adj(i)
    if ( nabor == 0 ) then
      return
    end if

    if ( marker(nabor) /= 0 ) then
      cycle
    end if
!
!  Include NABOR in the reachable set.
!
    if ( 0 <= deg(nabor) ) then
      rchsze = rchsze + 1
      rchset(rchsze) = nabor
      marker(nabor) = 1
      cycle
    end if
!
!  NABOR has been eliminated.  Find nodes reachable from it.
!
    marker(nabor) = -1
    nhdsze = nhdsze + 1
    nbrhd(nhdsze) = nabor

   20   continue

    jstrt = adj_row(nabor)
    jstop = adj_row(nabor+1) - 1

    do j = jstrt, jstop

      node = adj(j)
      nabor = -node

      if ( node < 0 ) then
        go to 20
      end if

      if ( node == 0 ) then
        cycle
      end if

      if ( marker(node) == 0 ) then
        rchsze = rchsze + 1
        rchset(rchsze) = node
        marker(node) = 1
      end if

    end do

  end do

  return
end
subroutine qmdupd ( adj_row, adj, nlist, list, deg, qsize, qlink, marker, &
  rchset, nbrhd, n )

!*****************************************************************************80
!
!! QMDUPD updates the node degrees, for the QMD method.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) NLIST, LIST(*), the nodes whose degree has to be updated.
!
!    Input/output, integer ( kind = 4 ) DEG(N), the degree vector.
!
!    Input/output, integer ( kind = 4 ) QSIZE, the size of indistinguishable supernodes.
!
!    Input/output, integer ( kind = 4 ) QLINK(N), linked list for indistinguishable nodes.
!
!    Input/output, integer ( kind = 4 ) MARKER(N), marks nodes in reach
!    and neighborhood sets.
!
!    Workspace, integer ( kind = 4 ) RCHSET(*), the reachable set.
!
!    Workspace, integer ( kind = 4 ) NBRHD(*), the neighborhood set.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) deg(n)
  integer ( kind = 4 ) deg0
  integer ( kind = 4 ) deg1
  integer ( kind = 4 ) il
  integer ( kind = 4 ) inhd
  integer ( kind = 4 ) inode
  integer ( kind = 4 ) irch
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) mark
  integer ( kind = 4 ) marker(n)
  integer ( kind = 4 ) nabor
  integer ( kind = 4 ) nbrhd(*)
  integer ( kind = 4 ) nhdsze
  integer ( kind = 4 ) nlist
  integer ( kind = 4 ) node
  integer ( kind = 4 ) qlink(n)
  integer ( kind = 4 ) qsize(n)
  integer ( kind = 4 ) rchset(*)
  integer ( kind = 4 ) rchsze
!
!  Find all eliminated supernodes that are adjacent to some nodes in the
!  given list.  Put them into (nhdsze,nbrhd).  DEG0 is the number of
!  nodes in the list.
!
  if ( nlist <= 0 ) then
    return
  end if

  deg0 = 0
  nhdsze = 0

  do il = 1, nlist

    node = list(il)
    deg0 = deg0 + qsize(node)
    jstrt = adj_row(node)
    jstop = adj_row(node+1) - 1

    do j = jstrt, jstop

      nabor = adj(j)

      if ( marker(nabor) == 0 .and. deg(nabor) < 0 ) then
        marker(nabor) = -1
        nhdsze = nhdsze + 1
        nbrhd(nhdsze) = nabor
      end if

    end do

  end do
!
!  Merge indistinguishable nodes in the list.
!
  if ( 0 < nhdsze ) then
    call qmdmrg ( adj_row, adj, deg, qsize, qlink, marker, deg0, nhdsze, &
      nbrhd, rchset, nbrhd(nhdsze+1), n )
  end if
!
!  Find the new degrees of the nodes that have not been merged.
!
  do il = 1, nlist

    node = list(il)
    mark = marker(node)

    if ( mark == 0 .or. mark == 1 ) then

      marker(node) = 2

      call qmdrch ( node, adj_row, adj, deg, marker, rchsze, rchset, nhdsze, &
        nbrhd, n )

      deg1 = deg0

      do irch = 1, rchsze
        inode = rchset(irch)
        deg1 = deg1 + qsize(inode)
        marker(inode) = 0
      end do

      deg(node) = deg1 - 1

      do inhd = 1, nhdsze
        inode = nbrhd(inhd)
        marker(inode) = 0
      end do

    end if

  end do

  return
end
subroutine rcm ( root, adj_row, adj, mask, perm, ccsize, n )

!*****************************************************************************80
!
!! RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
!
!  Discussion:
!
!    The connected component is specified by a node ROOT and a mask.
!    The numbering starts at the root node.
!
!    An outline of the algorithm is as follows:
!
!    X(1) = ROOT.
!
!    for ( I = 1 to N-1)
!      Find all unlabeled neighbors of X(I),
!      assign them the next available labels, in order of increasing degree.
!
!    When done, reverse the ordering.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the node that defines the connected component.
!    It is used as the starting point for the RCM ordering.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input/output, integer ( kind = 4 ) MASK(N), a mask for the nodes.  Only those nodes
!    with nonzero input mask values are considered by the routine.  The
!    nodes numbered by RCM will have their mask values set to zero.
!
!    Output, integer ( kind = 4 ) PERM(N), the RCM ordering.
!
!    Output, integer ( kind = 4 ) CCSIZE, the size of the connected component
!    that has been numbered.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!  Local parameters:
!
!    Workspace, integer ( kind = 4 ) DEG(N), a temporary vector used to
!    hold the degree of the nodes in the section graph specified by MASK and ROOT.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) ccsize
  integer ( kind = 4 ) deg(n)
  integer ( kind = 4 ) fnbr
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lbegin
  integer ( kind = 4 ) lnbr
  integer ( kind = 4 ) lperm
  integer ( kind = 4 ) lvlend
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) node
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) root
!
!  Find the degrees of the nodes in the component specified by MASK and ROOT.
!
  call degree ( root, adj_row, adj, mask, deg, ccsize, perm, n )

  mask(root) = 0

  if ( ccsize <= 1 ) then
    return
  end if

  lvlend = 0
  lnbr = 1
!
!  LBEGIN and LVLEND point to the beginning and
!  the end of the current level respectively.
!
  do while ( lvlend < lnbr )

    lbegin = lvlend + 1
    lvlend = lnbr

    do i = lbegin, lvlend
!
!  For each node in the current level...
!
      node = perm(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1
!
!  Find the unnumbered neighbors of NODE.
!
!  FNBR and LNBR point to the first and last neighbors
!  of the current node in PERM.
!
      fnbr = lnbr + 1

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then
          lnbr = lnbr + 1
          mask(nbr) = 0
          perm(lnbr) = nbr
        end if

      end do
!
!  Sort the neighbors of NODE in increasing order by degree.
!  Linear insertion is used.
!
      if ( fnbr < lnbr ) then

        k = fnbr

        do while ( k < lnbr )

          l = k
          k = k + 1
          nbr = perm(k)

          do while ( fnbr < l )

            lperm = perm(l)

            if ( deg(lperm) <= deg(nbr) ) then
              exit
            end if

            perm(l+1) = lperm
            l = l - 1

          end do

          perm(l+1) = nbr

        end do

      end if

    end do

  end do
!
!  We now have the Cuthill-McKee ordering.  Reverse it.
!
  call i4vec_reverse ( ccsize, perm )

  return
end
subroutine rcm_sub ( adj_row, adj, nsubg, subg, perm, n, mask )

!*****************************************************************************80
!
!! RCM_SUB finds the reverse Cuthill McKee ordering for a given subgraph.
!
!  Discussion:
!
!    The subgraph may be disconnected.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) NSUBG, the size of the subgraph.
!
!    Input, integer ( kind = 4 ) SUBG(NSUBG), the nodes in the subgraph.
!
!    Output, integer ( kind = 4 ) PERM(N), the permutation vector.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) MASK(N), a mask vector used to specify nodes in the
!    subgraph.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) XLS(N+1), index to a level structure.  The level
!    structure is stored in part of PERM.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ccsize
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) nlvl
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nsubg
  integer ( kind = 4 ) num
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) subg(n)
  integer ( kind = 4 ) xls(n+1)

  do i = 1, nsubg
    mask(subg(i)) = 1
  end do
!
!  For each connected component in the subgraph, call ROOT_FIND and RCM
!  for the ordering.
!
  num = 0

  do i = 1, nsubg

    node = subg(i)

    if ( 0 < mask(node) ) then

      call root_find ( node, adj_row, adj, mask, nlvl, xls, perm(num+1), n )

      call rcm ( node, adj_row, adj, mask, perm(num+1), ccsize, n )

      num = num + ccsize

      if ( nsubg <= num ) then
        return
      end if

    end if

  end do

  return
end
subroutine reach ( root, adj_row, adj, mask, marker, rchsze, rchset, nhdsze, &
  nbrhd, n )

!*****************************************************************************80
!
!! REACH determines the reachable set of a node through a subset in a subgraph.
!
!  Discussion:
!
!    The routine returns the neighborhood set of node Y in subset S, that is,
!    NBRHD(Y,S), the set of nodes in S that can be reached from Y through
!    a subset of S.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the given node, which is not in the subset S.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) MASK(N), the mask vector for the set S.
!    = 0, if the node is not in S,
!    > 0, if the node is in S.
!
!    Input/output, integer ( kind = 4 ) MARKER(N), the marker vector used to define
!    the subgraph.  Nodes in the subgraph have marker value 0.
!    On return, the reachable and neighborhood node sets have their
!    marker values reset to ROOT.
!
!    Output, integer ( kind = 4 ) RCHSZE, RCHSET(*), the reachable set.
!
!    Output, integer ( kind = 4 ) NHDSZE, NBRHD(*), the neighborhood set.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) istrt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) marker(n)
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) nabor
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nbrhd(*)
  integer ( kind = 4 ) nhdptr
  integer ( kind = 4 ) nhdsze
  integer ( kind = 4 ) node
  integer ( kind = 4 ) rchset(*)
  integer ( kind = 4 ) rchsze
  integer ( kind = 4 ) root
!
!  Initialization.
!
  nhdsze = 0
  rchsze = 0

  if ( marker(root) <= 0 ) then
    rchsze = 1
    rchset(1) = root
    marker(root) = root
  end if

  istrt = adj_row(root)
  istop = adj_row(root+1) - 1
  if ( istop < istrt ) then
    return
  end if
!
!  Loop through the neighbors of ROOT.
!
  do i = istrt, istop

    nabor = adj(i)

    if ( marker(nabor) /= 0 ) then
      cycle
    end if
!
!  If NABOR is not in subset S, include it in the reach set.
!
    if ( mask(nabor) <=  0 ) then
      rchsze = rchsze + 1
      rchset(rchsze) = nabor
      marker(nabor) = root
      cycle
    end if
!
!  NABOR is in subset S, and has not been considered.
!  Include it into the neighborhood set and find the nodes
!  reachable from ROOT through NABOR.
!
    nhdsze = nhdsze + 1
    nbrhd(nhdsze) = nabor
    marker(nabor) = root
    nhdptr = nhdsze

    do

      node = nbrhd(nhdptr)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1

      do j = jstrt, jstop

        nbr = adj(j)

        if ( marker(nbr) == 0 ) then

          if ( mask(nbr) /= 0 ) then
            nhdsze = nhdsze + 1
            nbrhd(nhdsze) = nbr
            marker(nbr) = root
          else
            rchsze = rchsze + 1
            rchset(rchsze) = nbr
            marker(nbr) = root
          end if

        end if

      end do

      nhdptr = nhdptr + 1

      if ( nhdsze < nhdptr ) then
        exit
      end if

    end do

  end do

  return
end
subroutine root_find ( root, adj_row, adj, mask, nlvl, xls, ls, n )

!*****************************************************************************80
!
!! ROOT_FIND finds pseudo-peripheral nodes.
!
!  Discussion:
!
!    The diameter of a graph is the maximum distance (number of edges)
!    between any two nodes of the graph.
!
!    The eccentricity of a node is the maximum distance between that
!    node and any other node of the graph.
!
!    A peripheral node is a node whose eccentricity equals the
!    diameter of the graph.
!
!    A pseudo-peripheral node is an approximation to a peripheral node;
!    it may be a peripheral node, but all we know is that we tried our
!    best.
!
!    The routine is given a graph, and seeks pseudo-peripheral nodes,
!    using a modified version of the scheme of Gibbs, Poole and
!    Stockmeyer.  It determines such a node for the section subgraph
!    specified by MASK and ROOT.
!
!    The routine also determines the level structure associated with
!    the given pseudo-peripheral node; that is, how far each node
!    is from the pseudo-peripheral node.  The level structure is
!    returned as a list of nodes LS, and pointers to the beginning
!    of the list of nodes that are at a distance of 0, 1, 2, ...,
!    N-1 from the pseudo-peripheral node.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!    Norman Gibbs, William Poole, Paul Stockmeyer,
!    An Algorithm for Reducing the Bandwidth
!    and Profile of a Sparse Matrix,
!    SIAM Journal on Numerical Analysis,
!    Volume 13, Number 2, April 1976, pages 236-250.
!
!    Norman Gibbs,
!    Algorithm 509:
!    A Hybrid Profile Reduction Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 2, Number 4, December 1976, pages 378-387.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ROOT.  On input, ROOT is a node in the
!    the component of the graph for which a pseudo-peripheral node is
!    sought.  On output, ROOT is the pseudo-peripheral node obtained.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) MASK(N), specifies a section subgraph.  Nodes for which
!    MASK is zero are ignored by ROOT_FIND.
!
!    Output, integer ( kind = 4 ) NLVL, is the number of levels in the level structure
!    rooted at the node ROOT.
!
!    Output, integer ( kind = 4 ) XLS(N+1), LS(N), the level structure array
!    pair containing the level structure found.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) ccsize
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kstop
  integer ( kind = 4 ) kstrt
  integer ( kind = 4 ) ls(n)
  integer ( kind = 4 ) mask(n)
  integer ( kind = 4 ) mindeg
  integer ( kind = 4 ) nabor
  integer ( kind = 4 ) ndeg
  integer ( kind = 4 ) nlvl
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nunlvl
  integer ( kind = 4 ) root
  integer ( kind = 4 ) xls(n+1)
!
!  Determine the level structure rooted at ROOT.
!
  call level_set ( root, adj_row, adj, mask, nlvl, xls, ls, n )
!
!  Count the number of nodes in this level structure.
!
  ccsize = xls(nlvl+1) - 1
!
!  Extreme case:
!    A complete graph has a level set of only a single level.
!    Every node is equally good (or bad).
!
  if ( nlvl == 1 ) then
    return
  end if
!
!  Extreme case:
!    A "line graph" 0--0--0--0--0 has every node in its only level.
!    By chance, we've stumbled on the ideal root.
!
  if ( nlvl == ccsize ) then
    return
  end if
!
!  Pick any node from the last level that has minimum degree
!  as the starting point to generate a new level set.
!
  do

    mindeg = ccsize

    jstrt = xls(nlvl)
    root = ls(jstrt)

    if ( jstrt < ccsize ) then

      do j = jstrt, ccsize

        node = ls(j)
        ndeg = 0
        kstrt = adj_row(node)
        kstop = adj_row(node+1) - 1

        do k = kstrt, kstop
          nabor = adj(k)
          if ( 0 < mask(nabor) ) then
            ndeg = ndeg + 1
          end if
        end do

        if ( ndeg < mindeg ) then
          root = node
          mindeg = ndeg
        end if

      end do

    end if
!
!  Generate the rooted level structure associated with this node.
!
    call level_set ( root, adj_row, adj, mask, nunlvl, xls, ls, n )
!
!  If the number of levels did not increase, then accept
!  the new ROOT.
!
    if ( nunlvl <= nlvl ) then
      exit
    end if

    nlvl = nunlvl
!
!  In the unlikely case that ROOT is one endpoint of a line graph,
!  we can exit now.
!
    if ( ccsize <= nlvl ) then
      exit
    end if

  end do

  return
end
subroutine rqtree ( leaf, adj_row, adj, perm, nblks, xblk, nodlvl, adjs, &
  stack, n )

!*****************************************************************************80
!
!! RQTREE finds a quotient tree ordering for a component, in the RQT method.
!
!  Discussion:
!
!    The component is specified by leaf and nodlvl.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LEAF, the input node that defines the connected
!    component.  It is also a leaf node in the rooted level structure
!    passed to RQTREE, that is, it has no neighbor in the next level.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) PERM(N), the ordering.
!
!    Output, integer ( kind = 4 ) NBLKS, XBLK(*), the quotient tree partitioning.
!
!    Input/output, integer ( kind = 4 ) NODLVL(N), the node level number vector.
!    Nodes in the component have their entry of NODLVL set to zero as
!    they are numbered.
!
!    Workspace, integer ( kind = 4 ) ADJS(*), the adjacent set of nodes in a
!    particular level.
!
!    Workspace, integer ( kind = 4 ) STACK(*), temporary vector used to
!    maintain the stack of node subsets.  it is organized as:
!    ( subset nodes, subset size, subset level ).
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) nblks
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) adjs(*)
  integer ( kind = 4 ) blksze
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) leaf
  integer ( kind = 4 ) level
  integer ( kind = 4 ) nadjs
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nodlvl(n)
  integer ( kind = 4 ) npop
  integer ( kind = 4 ) nuleaf
  integer ( kind = 4 ) num
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) stack(*)
  integer ( kind = 4 ) toplvl
  integer ( kind = 4 ) topstk
  integer ( kind = 4 ) xblk(*)
!
!  Initialize the stack vector and its pointers.
!
  stack(1) = 0
  stack(2) = 0
  topstk = 2
  toplvl = 0
  num = xblk(nblks+1) - 1
!
!  Form a leaf block, that is, one with no neighbors
!  in its next higher level.
!
  do

    level = nodlvl(leaf)
    nodlvl(leaf) = 0
    perm(num+1) = leaf
    blksze = 1

    call fnspan ( adj_row, adj, nodlvl, blksze, perm(num+1), level, nadjs, &
      adjs, nuleaf, n )

    if ( 0 < nuleaf ) then

      jp = num

      do j = 1, blksze
        jp = jp + 1
        node = perm(jp)
        nodlvl(node) = level
      end do

      leaf = nuleaf
      cycle

    end if
!
!  A new block has been found.
!
    do

      nblks = nblks + 1
      xblk(nblks) = num + 1
      num = num+blksze
!
!  Find the next possible block by using the adjacent
!  set in the lower level and the top node subset (if
!  appropriate) in the stack.
!
      level = level - 1

      if ( level <= 0 ) then
        xblk(nblks+1) = num + 1
        return
      end if

      call i4vec_copy ( nadjs, adjs, perm(num+1) )
      blksze = nadjs
!
!  The level of the node subset at the top of the
!  stack is the same as that of the adjacent set.
!  Pop the node subset from the stack.
!
      if ( level == toplvl ) then
        npop = stack(topstk-1)
        topstk = topstk - npop - 2
        ip = num + blksze + 1
        call i4vec_copy ( npop, stack(topstk+1), perm(ip) )
        blksze = blksze + npop
        toplvl = stack(topstk)
      end if

      call fnspan ( adj_row, adj, nodlvl, blksze, perm(num+1), level, nadjs, &
        adjs, nuleaf, n )

      if ( 0 < nuleaf ) then
        exit
      end if

    end do
!
!  Push the current node set into the stack.
!
    call i4vec_copy ( blksze, perm(num+1), stack(topstk+1) )
    topstk = topstk + blksze + 2
    stack(topstk-1) = blksze
    stack(topstk) = level
    toplvl = level
    leaf = nuleaf

  end do

end
subroutine smb_factor ( n, adj_row, adj, perm, perm_inv, ixlnz, nofnz, xnzsub, &
  nzsub, maxsub, rchlnk, mrglnk )

!*****************************************************************************80
!
!! SMB_FACTOR performs symbolic factorization on a permuted linear system.
!
!  Discussion:
!
!    The routine also sets up the compressed data structure for the system.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(N+1).  Information about row I is
!    stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(*), the adjacency structure.
!    For each row, the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) PERM(N), the permutation vector.
!
!    Input, integer ( kind = 4 ) PERM_INV(N), the inverse permutation vector.
!
!    Input/output, integer ( kind = 4 ) MAXSUB, the size of the subscript array NZSUB.
!    On return, the actual number of subscripts used.
!
!    Output, integer ( kind = 4 ) IXLNZ(N+1), index into the nonzero storage vector LNZ.
!
!    Output, integer ( kind = 4 ) XNZSUB(N+1), NZSUB(*), the compressed subscript vectors.
!
!    Output, integer ( kind = 4 ) NOFNZ, the number of nonzeros found.
!
!    Workspace, integer ( kind = 4 ) MRGLNK(N).  At the Kth step,
!    MRGLNK(K), MRGLNK(MRGLNK(K)),
!    is a list of those columns L(*,J) with J less than K,
!    such that its first off-diagonal nonzero is L(K,J).  Thus, the
!    nonzero structure of column L(*,K) can be found
!    by merging that of such columns L(*,J) with
!    the structure of A(*,K).
!
!    Workspace, integer ( kind = 4 ) RCHLNK(N), used to accumulate the structure of each
!    column L(*,K).  At the end of the K-th step,
!    RCHLNK(K), RCHLNK(RCHLNK(K)),
!    is the list of positions of nonzeros in column K of the factor L.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(*)
  integer ( kind = 4 ) adj_row(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inz
  integer ( kind = 4 ) ixlnz(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) knz
  integer ( kind = 4 ) kxsub
  integer ( kind = 4 ) lmax
  integer ( kind = 4 ) m
  integer ( kind = 4 ) marker(n)
  integer ( kind = 4 ) maxsub
  integer ( kind = 4 ) mrgk
  integer ( kind = 4 ) mrglnk(n)
  integer ( kind = 4 ) mrkflg
  integer ( kind = 4 ) nabor
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nofnz
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) nzbeg
  integer ( kind = 4 ) nzend
  integer ( kind = 4 ) nzsub(*)
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  integer ( kind = 4 ) rchlnk(n)
  integer ( kind = 4 ) rchm
  integer ( kind = 4 ) xnzsub(*)
!
!  Initialization.
!
  nzbeg = 1
  nzend = 0
  ixlnz(1) = 1
  mrglnk(1:n) = 0
!
!  MARKER is used to test if mass symbolic elimination can be performed.
!  That is, it is used to check whether the structure of the current
!  column K being processed is completely determined by the single
!  column MRGLNK(K).
!
  marker = 0
!
!  For each column KNZ counts the number
!  of nonzeros in column K accumulated in RCHLNK.
!
  np1 = n + 1

  do k = 1, n

    knz = 0
    mrgk = mrglnk(k)
    mrkflg = 0
    marker(k) = k
    if ( mrgk /= 0 ) then
      marker(k) = marker(mrgk)
    end if
    xnzsub(k) = nzend
    node = perm(k)
    jstrt = adj_row(node)
    jstop = adj_row(node+1) - 1

    if ( jstop < jstrt ) then
      go to 160
    end if
!
!  Use RCHLNK to link through the structure of A(*,K) below diagonal.
!
    rchlnk(k) = np1

    do j = jstrt, jstop

      nabor = adj(j)
      nabor = perm_inv(nabor)

      if ( k < nabor ) then

        rchm = k

        do

          m = rchm
          rchm = rchlnk(m)

          if ( nabor < rchm ) then
            exit
          end if

        end do

        knz = knz + 1
        rchlnk(m) = nabor
        rchlnk(nabor) = rchm

        if ( marker(nabor) /= marker(k) ) then
          mrkflg = 1
        end if

      end if

    end do
!
!  Test for mass symbolic elimination.
!
    lmax = 0

    if ( mrkflg /= 0 .or. mrgk == 0 ) then
      go to 40
    end if

    if ( mrglnk(mrgk) /= 0 ) then
      go to 40
    end if

    xnzsub(k) = xnzsub(mrgk) + 1
    knz = ixlnz(mrgk+1) - ( ixlnz(mrgk) + 1 )
    go to 150
!
!  Link through each column I that affects L(*,K).
!
   40   continue

    i = k

   50   continue

    i = mrglnk(i)
    if ( i == 0 ) then
      go to 90
    end if

    inz = ixlnz(i+1) - ( ixlnz(i) + 1 )
    jstrt = xnzsub(i) + 1
    jstop = xnzsub(i) + inz

    if ( lmax < inz ) then
      lmax = inz
      xnzsub(k) = jstrt
    end if
!
!  Merge structure of L(*,I) in NZSUB into RCHLNK.
!
    rchm = k

    do j = jstrt, jstop

      nabor = nzsub(j)

      do

        m = rchm
        rchm = rchlnk(m)

        if ( nabor <= rchm ) then
          exit
        end if

      end do

      if ( rchm /= nabor ) then
        knz = knz + 1
        rchlnk(m) = nabor
        rchlnk(nabor) = rchm
        rchm = nabor
      end if

    end do

    go to 50
!
!  Check if subscripts duplicate those of another column...
!
90   continue

    if ( knz == lmax ) then
      go to 150
    end if
!
!  ...or if tail of column K-1 matches head of column K.
!
    if ( nzend < nzbeg ) then
      go to 130
    end if

    i = rchlnk(k)

    do jstrt = nzbeg, nzend

      if ( nzsub(jstrt) == i ) then
        go to 110
      end if

      if ( i < nzsub(jstrt) ) then
        go to 130
      end if

    end do

    go to 130

  110   continue

    xnzsub(k) = jstrt

    do j = jstrt, nzend
      if ( nzsub(j) /= i ) then
        go to 130
      end if
      i = rchlnk(i)
      if ( n < i ) then
        go to 150
      end if
    end do

    nzend = jstrt - 1
!
!  Copy the structure of L(*,K) from rchlnk
!  to the data structure (xnzsub,nzsub).
!
  130   continue

    nzbeg = nzend + 1
    nzend = nzend + knz

    if ( maxsub < nzend ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SMB_FACTOR - Fatal error!'
      write ( *, '(a)' ) '  Insufficient storage for nonzero entries.'
      write ( *, '(a,i8)' ) '  MAXSUB  = ', maxsub
      write ( *, '(a,i8)' ) '  Exceeded by NZEND = ', nzend
      stop
    end if

    i = k

    do j = nzbeg, nzend
      i = rchlnk(i)
      nzsub(j) = i
      marker(i) = k
    end do

    xnzsub(k) = nzbeg
    marker(k) = k
!
!  Update the vector mrglnk.  Note column L(*,K) just found
!  is required to determine column L(*,J), where
!  L(J,K) is the first nonzero in L(*,K) below diagonal.
!
  150   continue

    if ( 1 < knz ) then
      kxsub = xnzsub(k)
      i = nzsub(kxsub)
      mrglnk(k) = mrglnk(i)
      mrglnk(i) = k
    end if

  160   continue

    ixlnz(k+1) = ixlnz(k) + knz

  end do

  nofnz = ixlnz(n) - 1
  maxsub = xnzsub(n)
  xnzsub(n+1) = xnzsub(n)

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
subroutine ts_factor ( nblks, xblk, father, diag, xenv, env, xnonz, nonz, &
  nzsubs, n, ierror )

!*****************************************************************************80
!
!! TS_FACTOR performs the symmetric factorization of a tree-partitioned system.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBLKS, XBLK(NBLKS+1), FATHER(N),
!    the tree partitioning.
!
!    Input/output, real ( kind = 8 ) DIAG(N), integer ( kind = 4 ) XENV(N+1),
!    real ( kind = 8 ) ENV(*),
!    storage arrays for the diagonal blocks of the matrix.  On output, the
!    diagonal blocks of the factor.
!
!    Input, integer ( kind = 4 ) XNONZ(N+1), real ( kind = 8 ) NONZ(*),
!    integer ( kind = 4 ) NZSUBS(*),
!    the off-diagonal nonzeros in the original matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error, the factorization was carried out.
!    1, the matrix is not positive definite.
!
!  Local parameters:
!
!    integer ( kind = 4 ) FIRST(N), temporary vector used to facilitate the
!    indexing to the vector NONZ (or NZSUBS) for non-null subcolumns
!    in off-diagonal blocks.
!
  implicit none

  integer ( kind = 4 ) nblks
  integer ( kind = 4 ) n

  integer ( kind = 4 ) blksze
  integer ( kind = 4 ) col
  integer ( kind = 4 ) col1
  integer ( kind = 4 ) colbeg
  integer ( kind = 4 ) colend
  integer ( kind = 4 ) colsze
  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) env(*)
  integer ( kind = 4 ) father(n)
  integer ( kind = 4 ) first(n)
  integer ( kind = 4 ) fnz
  integer ( kind = 4 ) fnz1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) istrt
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kenv
  integer ( kind = 4 ) kenv0
  integer ( kind = 4 ) kfathr
  real ( kind = 8 ) nonz(*)
  integer ( kind = 4 ) nzsubs(*)
  integer ( kind = 4 ) row
  integer ( kind = 4 ) rowbeg
  integer ( kind = 4 ) rowend
  real ( kind = 8 ) s
  real ( kind = 8 ) temp(n)
  integer ( kind = 4 ) xblk(nblks+1)
  integer ( kind = 4 ) xenv(n+1)
  integer ( kind = 4 ) xnonz(n+1)

  ierror = 0

  temp(1:n) = 0.0D+00
  first(1:n) = xnonz(1:n)
!
!  Loop through the blocks.
!
  write (  *, * ) 'DEBUG: TS_FACTOR has NBLKS = ', nblks

  do k = 1, nblks

    write ( *, * ) 'DEBUG: TS_FACTOR has K = ', k

    rowbeg = xblk(k)
    rowend = xblk(k+1) - 1
    blksze = rowend - rowbeg + 1

    write ( *, * ) 'DEBUG: TS_FACTOR calls ES_FACTOR'

    call es_factor ( blksze, xenv(rowbeg), env, diag(rowbeg), ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TS_FACTOR - Fatal error!'
      write ( *, '(a,i8)' ) '  ES_FACTOR returns IERROR = ', ierror
      return
    end if
!
!  Modify the father diagonal block A(FATHER(K),FATHER(K)) from the
!  off-diagonal block A(K,FATHER(K)).
!
    write ( *, * ) 'DEBUG: TS_FACTOR modifies father diagonal.'

    kfathr = father(k)

    write ( *, * ) 'DEBUG: TS_FACTOR has KFATHR = ', kfathr

    if ( kfathr <= 0 ) then
      cycle
    end if

    colbeg = xblk(kfathr)
    colend = xblk(kfathr+1) - 1

    write ( *, * ) 'DEBUG: TSFACTOR has COLBEG = ', colbeg, ' COLEND = ', colend
!
!  Find the first and last non-null column in the off-diagonal block.
!  Reset COLBEG and COLEND.
!
    do col = colbeg, colend

      jstrt = first(col)
      jstop = xnonz(col+1) - 1

      if ( jstrt <= jstop .and. nzsubs(jstrt) <= rowend ) then
        exit
      end if

    end do

    colbeg = col
    col = colend

    do col1 = colbeg, colend

      jstrt = first(col)
      jstop = xnonz(col+1) - 1

      if ( jstrt <= jstop .and. nzsubs(jstrt) <= rowend ) then
        exit
      end if

      col = col - 1

    end do

    colend = col

    do col = colbeg, colend

      jstrt = first(col)
      jstop = xnonz(col+1) - 1
!
!  Test for null subcolumn.  FNZ stores the first nonzero subscript
!  in the block column.
!
      if ( jstop < jstrt ) then
        go to 130
      end if

      fnz = nzsubs(jstrt)
      if ( rowend < fnz ) then
        go to 130
      end if
!
!  Unpack a column in the off-diagonal block and perform upper and
!  lower solves on the unpacked column.
!
      do j = jstrt, jstop
        row = nzsubs(j)
        if ( rowend < row ) then
          exit
        end if
        temp(row) = nonz(j)
      end do

      colsze = rowend - fnz + 1

      call el_solve ( colsze, xenv(fnz), env, diag(fnz), temp(fnz) )

      call eu_solve ( colsze, xenv(fnz), env, diag(fnz), temp(fnz) )
!
!  Do the modification by looping through
!  the columns and forming inner products.
!
      kenv0 = xenv(col+1) - col

      do col1 = colbeg, colend

        istrt = first(col1)
        istop = xnonz(col1+1) - 1
!
!  Check to see if subcolumn is null.
!
        fnz1 = nzsubs(istrt)
        if ( istop < istrt .or. rowend < fnz1 ) then
          cycle
        end if
!
!  Check if inner product should be done.
!
        if ( fnz1 < fnz ) then
          cycle
        end if

        if ( fnz1 == fnz .and. col1 < col ) then
          cycle
        end if

        s = 0.0D+00
        do i = istrt, istop
          isub = nzsubs(i)
          if ( rowend < isub ) then
            exit
          end if
          s = s + temp(isub) * nonz(i)
        end do
!
!  Modify ENV or DIAG.
!
        if ( col1 /= col ) then

          kenv = kenv0 + col1
          if ( col < col1 ) then
            kenv = xenv(col1+1) - col1 + col
          end if
          env(kenv) = env(kenv) - s

        else

          diag(col1) = diag(col1) - s

        end if

      end do
!
!  Reset part of the TEMP vector to zero.
!
      temp(fnz:rowend) = 0.0D+00

130   continue

    end do

    write ( *, * ) 'DEBUG: TS_FACTOR got this far.'
!
!  Update the first vector for columns in FATHER(K) block, so that it
!  will index to the beginning of the next off-diagonal block to be
!  considered.
!
    do col = colbeg, colend

      jstrt = first(col)
      jstop = xnonz(col+1) - 1

      if ( jstrt <= jstop ) then

        do j = jstrt, jstop

          row = nzsubs(j)

          if ( rowend < row ) then
            first(col) = j
            go to 150
          end if

        end do

        first(col) = jstop + 1

      end if

150   continue

    end do

  end do

  return
end
subroutine ts_solve ( nblks, xblk, diag, xenv, env, xnonz, nonz, nzsubs, rhs, &
  n )

!*****************************************************************************80
!
!! TS_SOLVE solves a tree-partitioned factored system.
!
!  Discussion:
!
!    Implicit back substitution is used.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBLKS, XBLK(NBLKS+1), the partitioning.
!
!    Input, integer ( kind = 4 ) XENV(N+1), real ( kind = 8 ) ENV(*), the envelope of the
!    diagonal blocks.
!
!    Input, integer ( kind = 4 ) XNONZ(N+1), real ( kind = 8 ) NONZ(*),
!    integer ( kind = 4 ) NZSUBS(*),
!    the data structure for the off-block diagonal nonzeros.
!
!    Input/output, real ( kind = 8 ) RHS(N), on input the right hand vector.
!    on output, the solution vector.
!
  implicit none

  integer ( kind = 4 ) nblks
  integer ( kind = 4 ) n

  integer ( kind = 4 ) col
  integer ( kind = 4 ) col1
  integer ( kind = 4 ) col2
  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) env(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) last
  integer ( kind = 4 ) ncol
  real ( kind = 8 ) nonz(*)
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nzsubs(*)
  real ( kind = 8 ) rhs(n)
  integer ( kind = 4 ) row
  integer ( kind = 4 ) row1
  integer ( kind = 4 ) row2
  real ( kind = 8 ) s
  real ( kind = 8 ) temp(n)
  integer ( kind = 4 ) xblk(nblks+1)
  integer ( kind = 4 ) xenv(n+1)
  integer ( kind = 4 ) xnonz(n+1)
!
!  Forward substitution.
!
  do i = 1, nblks

    row1 = xblk(i)
    row2 = xblk(i+1) - 1
    last = xnonz(row2+1)
!
!  Modify the right hand side vector by the product of the
!  off-diagonal block with the corresponding part of the right hand side.
!
    if ( i /= 1 .and. last /= xnonz(row1) ) then

      do row = row1, row2

        jstrt = xnonz(row)

        if ( jstrt /= last ) then

          jstop = xnonz(row+1) - 1

          if ( jstrt <= jstop ) then

            s = 0.0D+00
            do j = jstrt, jstop
              col = nzsubs(j)
              s = s + rhs(col) * nonz(j)
            end do

            rhs(row) = rhs(row) - s

          end if

        end if

      end do

    end if

    nrow = row2 - row1 + 1

    call el_solve ( nrow, xenv(row1), env, diag(row1), rhs(row1) )

    call eu_solve ( nrow, xenv(row1), env, diag(row1), rhs(row1) )

  end do
!
!  Backward solution.
!
  if ( nblks == 1 ) then
    return
  end if

  last = xblk(nblks) - 1

  temp(1:last) = 0.0D+00

  i = nblks
  col1 = xblk(i)
  col2 = xblk(i+1) - 1

  do

    if ( i == 1 ) then
      exit
    end if

    last = xnonz(col2+1)
!
!  Multiply the off-diagonal block by the corresponding
!  part of the solution vector and store in TEMP.
!
    if ( last /= xnonz(col1) ) then

      do col = col1, col2

        s = rhs(col)

        if ( s /= 0.0D+00 ) then

          jstrt = xnonz(col)
          if ( jstrt == last ) then
            exit
          end if

          jstop = xnonz(col+1) - 1

          do j = jstrt, jstop
            row = nzsubs(j)
            temp(row) = temp(row) + s * nonz(j)
          end do

        end if

      end do

    end if

    i = i - 1
    col1 = xblk(i)
    col2 = xblk(i+1) - 1
    ncol = col2 - col1 + 1

    call el_solve ( ncol, xenv(col1), env, diag(col1), temp(col1) )

    call eu_solve ( ncol, xenv(col1), env, diag(col1), temp(col1) )

    do j = col1, col2
      rhs(j) = rhs(j) - temp(j)
    end do

  end do

  return
end
