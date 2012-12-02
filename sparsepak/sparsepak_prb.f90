program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSEPAK_PRB.
!
!  Discussion:
!
!    SPARSEPAK_PRB runs the SPARSEPAK tests.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEPAK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPARSEPAK library.'

  if ( .true. ) then
    write ( *, * ) ' '
    write ( *, * ) 'TEST01 is temporarily being skipped!'
  else
    call test01 ( )
  end if
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEPAK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests the 1WD method.
!
!  Discussion:
!
!    This example involves the following equations:
!
!    2*x1  - x10       =  0
!    2*x2  - x9  - x10 =  0
!    2*x3  - x8  - x9  =  0
!    2*x4  - x7  - x8  =  0
!    2*x5  - x6  - x7  =  0
!    2*x6  - x5        = 11
!    2*x7  - x4  - x5  =  0
!    2*x8  - x3  - x4  =  0
!    2*x9  - x2  - x3  =  0
!    2*x10 - x1  - x2  =  0
!
!    with solution (1,3,5,7,9,10,8,6,4,2).
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxadj = 300
  integer ( kind = 4 ), parameter :: maxblk = 10
  integer ( kind = 4 ), parameter :: maxenv = 300
  integer ( kind = 4 ), parameter :: maxnon = 300
  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) adj(maxadj)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) env(maxenv)
  integer ( kind = 4 ) father(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iband
  integer ( kind = 4 ) env_size
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) marker(n)
  integer ( kind = 4 ) nadj
  integer ( kind = 4 ) nblks
  real    ( kind = 8 ) nonz(maxnon)
  integer ( kind = 4 ) nzsub(maxnon)
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  integer ( kind = 4 ) rchset(n)
  real    ( kind = 8 ) rhs(n)
  integer ( kind = 4 ) xadj(n+1)
  integer ( kind = 4 ) xblk(maxblk+1)
  integer ( kind = 4 ) xenv(n+1)
  integer ( kind = 4 ) xnonz(n+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use the 1WD method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The matrix order is N = ', n
!
!  Initialize the permutation vectors.
!
  call i4vec_indicator ( n, perm )
  call i4vec_indicator ( n, perm_inv )
!
!  Store the adjacency information.
!
  call adj_set_1 ( adj, maxadj, nadj, n, xadj )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Number of adjacency entries NADJ = ', nadj
!
!  Display the adjacency information.
!
  call adj_print ( n, nadj, xadj, adj )
!
!  Determine the initial envelope size.
!
  call adj_env_size ( n, xadj, nadj, adj, perm, perm_inv, env_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial envelope size is ', env_size
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Generate the 1WD ordering.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GEN1WD generates the 1WD ordering.'

  call gen1wd ( n, xadj, adj, nblks, xblk, perm )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The number of blocks is ', nblks
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Compute the inverse ordering.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PERM_INVERSE computes the inverse ordering.'

  call perm_inverse ( n, perm, perm_inv )
!
!  Print orderings
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I    Perm(I)   InvPerm(I)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(3i6)' ) i, perm(i), perm_inv(i)
  end do
!
!  Determine the reordered envelope size.
!
  call adj_env_size ( n, xadj, nadj, adj, perm, perm_inv, env_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The reordered envelope size is ', env_size
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Determine the envelope.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FNBENV determines the envelope'

  call fnbenv ( xadj, adj, perm, perm_inv, nblks, xblk, xenv, &
    env_size, marker, rchset, n )
!
!  Set RHS, DIAG, ENV.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SETSY1 sets RHS, DIAG and ENV.'

  call setsy1 ( diag, env, adj, perm_inv, xadj, maxadj, maxenv, &
    n, nonz, nzsub, rhs, xenv, xnonz )
!
!  Factor the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TS_FACTOR factors the system.'

  call ts_factor ( nblks, xblk, father, diag, xenv, env, xnonz, &
    nonz, nzsub, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a,i6)' ) '  TS_FACTOR returns error flag IERROR = ', ierror
    return
  end if
!
!  Solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TS_SOLVE solves the system.'

  call ts_solve ( nblks, xblk, diag, xenv, env, xnonz, nonz, &
    nzsub, rhs, n )
!
!  Unpermute the solution.
!
  call perm_rv ( n, rhs, perm )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution:'
  write ( *, '(a)' ) ' '
  write ( *, '(5g14.6)' ) rhs(1:n)

  return
end
subroutine setsy1 ( diag, env, adj, perm_inv, xadj, maxadj, maxenv, n, nonz, &
  nzsub, rhs, xenv, xnonz )

!*****************************************************************************80
!
!! SETSY1 stores the numerical values defining problem 1.
!
!  Discussion:
!
!    There is only one nonzero right hand side entry.
!    The matrix diagonal entries are all 2.
!    The nonzero offdiagonal entries are all -1.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) maxadj
  integer ( kind = 4 ) maxenv
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(maxadj)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) env(maxenv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) j
  integer ( kind = 4 ) xadj(n+1)
  integer ( kind = 4 ) jsub
  real    ( kind = 8 ) nonz(*)
  integer ( kind = 4 ) nzsub(*)
  integer ( kind = 4 ) perm_inv(n)
  real    ( kind = 8 ) rhs(n)
  real    ( kind = 8 ) value
  integer ( kind = 4 ) xenv(n+1)
  integer ( kind = 4 ) xnonz(n+1)
!
!  Zero out storage.
!
  rhs(1:n) = 0.0D+00
  diag(1:n) = 0.0D+00
  env(1:maxenv) = 0.0D+00
!
!  Set the nonzero elements of the right hand side vector.
!
  isub = 6
  value = 11.0D+00

  call addrhs ( perm_inv, isub, n, rhs, value )
!
!  Set the diagonal entries of the matrix.
!
  diag(1:n) = 2.0D+00
!
!  Set the off diagonal terms of the matrix.
!
  do i = 1, n

    isub = i

    do j = xadj(i), xadj(i+1) - 1

      jsub = adj(j)
      value = -1.0D+00
      call addrqt ( isub, jsub, value, perm_inv, diag, xenv, env, &
        xnonz, nonz, nzsub, n )

    end do

  end do

  return
end
subroutine adj_set_1 ( adj, maxadj, nadj, n, xadj )

!*****************************************************************************80
!
!! ADJ_SET_1 sets up the adjacency structure for problem 1.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nonz = 19

  integer ( kind = 4 ) maxadj
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(maxadj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter, dimension ( nonz ) :: ilist = (/ &
    -1,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10/)
  integer ( kind = 4 ), parameter, dimension ( nonz) :: jlist = (/ &
    -1,10,10,9,9,8,7,8,7,6,5,4,5,3,4,2,3,1,2/)
  integer ( kind = 4 ) nadj
  integer ( kind = 4 ) xadj(n+1)

  do i = 1, nonz
    call adj_set ( adj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests the ND method.
!
!  Discussion:
!
!    This example involves the following equations:
!
!      2*x1  - x10       = 0
!      2*x2  - x9  - x10 = 0
!      2*x3  - x8  - x9  = 0
!      2*x4  - x7  - x8  = 0
!      2*x5  - x6  - x7  = 0
!      2*x6  - x5        = 11
!      2*x7  - x4  - x5  = 0
!      2*x8  - x3  - x4  = 0
!      2*x9  - x2  - x3  = 0
!      2*x10 - x1  - x2  = 0
!
!    with solution
!
!      x = (1,3,5,7,9,10,8,6,4,2).
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxadj = 300
  integer ( kind = 4 ), parameter :: maxenv = 300
  integer ( kind = 4 ), parameter :: maxnzsub = 300
  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) adj(maxadj)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) env(maxenv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iband
  integer ( kind = 4 ) env_size
  integer ( kind = 4 ) ixlnz(n+1)
  integer ( kind = 4 ) maxsub
  integer ( kind = 4 ) mrglnk(n)
  integer ( kind = 4 ) nadj
  integer ( kind = 4 ) nofnz
  integer ( kind = 4 ) nzsub(maxnzsub)
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  integer ( kind = 4 ) rchlnk(n)
  real    ( kind = 8 ) rhs(n)
  integer ( kind = 4 ) xadj(n+1)
  real    ( kind = 8 ) xlnz(99)
  integer ( kind = 4 ) xnzsub(99)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Use the ND method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The matrix order is N = ', n
!
!  Initialize the permutation vectors.
!
  call i4vec_indicator ( n, perm )
  call i4vec_indicator ( n, perm_inv )
!
!  Store the adjacency information.
!
  call adj_set_2 ( adj, maxadj, nadj, n, xadj )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Number of adjacency entries NADJ = ', nadj
!
!  Display adjacency information.
!
  call adj_print ( n, nadj, xadj, adj )
!
!  Determine the initial envelope size.
!
  call adj_env_size ( n, xadj, nadj, adj, perm, perm_inv, env_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial envelope size is ', env_size
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Generate the 1WD ordering.
!
  call gennd ( n, xadj, adj, perm )
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Compute the inverse ordering.
!
  call perm_inverse ( n, perm, perm_inv )
!
!  Print orderings
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I    Perm(I)   InvPerm(I)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(3i6)' ) i, perm(i), perm_inv(i)
  end do
!
!  Determine the reordered envelope size.
!
  call adj_env_size ( n, xadj, nadj, adj, perm, perm_inv, env_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The reordered envelope size is ', env_size
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Symbolic factorization.
!
  maxsub = maxnzsub

  call smb_factor ( n, xadj, adj, perm, perm_inv, ixlnz, nofnz, &
    xnzsub, nzsub, maxsub, rchlnk, mrglnk )
!
!  Set RHS, DIAG, ENV.
!
  xlnz(1:99) = 0.0D+00

  call setsy2 ( diag, env, adj, perm_inv, xadj, maxadj, maxenv, &
    n, rhs, nzsub, xnzsub, xlnz, ixlnz )
!
!  Factor the system.
!
  call gs_factor ( n, ixlnz, xlnz, xnzsub, nzsub, diag )
!
!  Solve the system.
!
  call gs_solve ( n, ixlnz, xlnz, xnzsub, nzsub, diag, rhs )
!
!  Unpermute the solution.
!
  call perm_rv ( n, rhs, perm )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution:'
  write ( *, '(a)' ) ' '
  write ( *, '(5g14.6)' ) rhs(1:n)

  return
end
subroutine setsy2 ( diag, env, adj, perm_inv, xadj, maxadj, maxenv, n, rhs, &
  nzsub, xnzsub, xlnz, ixlnz )

!*****************************************************************************80
!
!! SETSY2 stores the numerical values defining problem 2.
!
!  Discussion:
!
!    There is only one nonzero right hand side entry.
!    The matrix diagonal entries are all 2.
!    The nonzero offdiagonal entries are all -1.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) maxadj
  integer ( kind = 4 ) maxenv
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(maxadj)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) env(maxenv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jsub
  integer ( kind = 4 ) ixlnz(n+1)
  integer ( kind = 4 ) nzsub(*)
  integer ( kind = 4 ) perm_inv(n)
  real    ( kind = 8 ) rhs(n)
  real    ( kind = 8 ) value
  integer ( kind = 4 ) xadj(n+1)
  real    ( kind = 8 ) xlnz(*)
  integer ( kind = 4 ) xnzsub(n)
!
!  Zero out storage.
!
  rhs(1:n) = 0.0D+00
  diag(1:n) = 0.0D+00
  env(1:maxenv) = 0.0D+00
!
!  Set the nonzero elements of the right hand side vector.
!
  isub = 6
  value = 11.0D+00

  call addrhs ( perm_inv, isub, n, rhs, value )
!
!  Set the diagonal entries of the matrix.
!
  diag(1:n) = 2.0D+00
!
!  Set the off diagonal terms of the matrix.
!
  do i = 1, n

    isub = i

    do j = xadj(i), xadj(i+1) - 1

      jsub = adj(j)
      value = -1.0D+00

      call addcom ( isub, jsub, value, perm_inv, diag, xlnz, ixlnz, nzsub, &
        xnzsub, n )

    end do

  end do

  return
end
subroutine adj_set_2 ( adj, maxadj, nadj, n, xadj )

!*****************************************************************************80
!
!! ADJ_SET_2 sets up the adjacency structure for problem 2.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nonz = 19

  integer ( kind = 4 ) maxadj
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(maxadj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter, dimension ( nonz ) :: ilist = (/ &
    -1,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10/)
  integer ( kind = 4 ), parameter, dimension ( nonz ) :: jlist = (/ &
    -1,10,10,9,9,8,7,8,7,6,5,4,5,3,4,2,3,1,2/)
  integer ( kind = 4 ) nadj
  integer ( kind = 4 ) xadj(n+1)

  do i = 1, nonz
    call adj_set ( adj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests the QMD method.
!
!  Discussion:
!
!    This example involves the following equations:
!
!      2*x1  - x10       = 0
!      2*x2  - x9  - x10 = 0
!      2*x3  - x8  - x9  = 0
!      2*x4  - x7  - x8  = 0
!      2*x5  - x6  - x7  = 0
!      2*x6  - x5        = 11
!      2*x7  - x4  - x5  = 0
!      2*x8  - x3  - x4  = 0
!      2*x9  - x2  - x3  = 0
!      2*x10 - x1  - x2  = 0
!
!    with solution
!
!      x = (1,3,5,7,9,10,8,6,4,2).
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxadj = 300
  integer ( kind = 4 ), parameter :: maxenv = 300
  integer ( kind = 4 ), parameter :: maxnzsub = 300
  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) adj(maxadj)
  integer ( kind = 4 ) adj2(maxadj)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) env(maxenv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iband
  integer ( kind = 4 ) env_size
  integer ( kind = 4 ) ixlnz(n+1)
  integer ( kind = 4 ) marker(n)
  integer ( kind = 4 ) maxsub
  integer ( kind = 4 ) mrglnk(n)
  integer ( kind = 4 ) nadj
  integer ( kind = 4 ) nbrhd(n)
  integer ( kind = 4 ) nofnz
  integer ( kind = 4 ) nofsub
  integer ( kind = 4 ) nzsub(maxnzsub)
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  integer ( kind = 4 ) qlink(n)
  integer ( kind = 4 ) qsize(n)
  integer ( kind = 4 ) rchlnk(n)
  integer ( kind = 4 ) rchset(n)
  real    ( kind = 8 ) rhs(n)
  integer ( kind = 4 ) xadj(n+1)
  real    ( kind = 8 ) xlnz(99)
  integer ( kind = 4 ) xnzsub(n+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Use the QMD method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The matrix order is N = ', n
!
!  Initialize the permutation vectors.
!
  call i4vec_indicator ( n, perm )
  call i4vec_indicator ( n, perm_inv )
!
!  Store the adjacency information.
!
  call adj_set_3 ( adj, maxadj, nadj, n, xadj )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The number of adjacency entries NADJ = ', nadj
!
!  Display adjacency information.
!
  call adj_print ( n, nadj, xadj, adj )
!
!  Determine the initial envelope size.
!
  call adj_env_size ( n, xadj, nadj, adj, perm, perm_inv, env_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial envelope size is ', env_size
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Generate the QMD ordering.
!
  adj2(1:nadj) = adj(1:nadj)

  call genqmd ( n, xadj, adj2, perm, perm_inv, marker, rchset, nbrhd, qsize, &
    qlink, nofsub )
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Compute the inverse ordering.
!
  call perm_inverse ( n, perm, perm_inv )
!
!  Print orderings
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I    Perm(I)   InvPerm(I)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(3i6)' ) i, perm(i), perm_inv(i)
  end do
!
!  Determine the reordered envelope size.
!
  call adj_env_size ( n, xadj, nadj, adj, perm, perm_inv, env_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The reordered envelope size is ', env_size
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Symbolic factorization.
!
  maxsub = maxnzsub

  call smb_factor ( n, xadj, adj, perm, perm_inv, ixlnz, nofnz, xnzsub, &
    nzsub, maxsub, rchlnk, mrglnk )
!
!  Set RHS, DIAG, ENV.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IXLNZ array:'
  write ( *, '(a)' ) ' '
  do i = 1, n + 1
    write ( *, '(2x,i8,2x,i8)' ) i, ixlnz(i)
  end do

  xlnz(1:99) = 0.0D+00

  call setsy3 ( diag, env, adj, perm_inv, xadj, maxadj, maxenv, n, rhs, &
    nzsub, xnzsub, xlnz, ixlnz )
!
!  Factor the system.
!
  call gs_factor ( n, ixlnz, xlnz, xnzsub, nzsub, diag )
!
!  Solve the system.
!
  call gs_solve ( n, ixlnz, xlnz, xnzsub, nzsub, diag, rhs )
!
!  Unpermute the solution.
!
  call perm_rv ( n, rhs, perm )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution:'
  write ( *, '(a)' ) ' '
  write ( *, '(5g14.6)' ) rhs(1:n)

  return
end
subroutine setsy3 ( diag, env, adj, perm_inv, xadj, maxadj, maxenv, n, rhs, &
  nzsub, xnzsub, xlnz, ixlnz )

!*****************************************************************************80
!
!! SETSY3 stores the numerical values defining problem 3.
!
!  Discussion:
!
!    There is only one nonzero right hand side entry.
!    The matrix diagonal entries are all 2.
!    The nonzero offdiagonal entries are all -1.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) maxadj
  integer ( kind = 4 ) maxenv
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(maxadj)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) env(maxenv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) ixlnz(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jsub
  integer ( kind = 4 ) nzsub(*)
  integer ( kind = 4 ) perm_inv(n)
  real    ( kind = 8 ) rhs(n)
  real    ( kind = 8 ) value
  integer ( kind = 4 ) xadj(n+1)
  real    ( kind = 8 ) xlnz(*)
  integer ( kind = 4 ) xnzsub(n)
!
!  Zero out storage.
!
  rhs(1:n) = 0.0D+00
  diag(1:n) = 0.0D+00
  env(1:maxenv) = 0.0D+00
!
!  Set the nonzero elements of the right hand side vector.
!
  isub = 6
  value = 11.0D+00

  call addrhs ( perm_inv, isub, n, rhs, value )
!
!  Set the diagonal entries of the matrix.
!
  diag(1:n) = 2.0D+00
!
!  Set the off diagonal terms of the matrix.
!
  do i = 1, n

    isub = i

    do j = xadj(i), xadj(i+1) - 1

      jsub = adj(j)
      value = -1.0D+00

      call addcom ( isub, jsub, value, perm_inv, diag, xlnz, ixlnz, nzsub, &
        xnzsub, n )

    end do

  end do

  return
end
subroutine adj_set_3 ( adj, maxadj, nadj, n, xadj )

!*****************************************************************************80
!
!! ADJ_SET_3 sets up the adjacency structure for problem 3.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nonz = 19

  integer ( kind = 4 ) maxadj
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(maxadj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save, dimension ( nonz ) :: ilist = (/ &
    -1,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10/)
  integer ( kind = 4 ), save, dimension ( nonz ) :: jlist = (/ &
    -1,10,10,9,9,8,7,8,7,6,5,4,5,3,4,2,3,1,2/)
  integer ( kind = 4 ) nadj
  integer ( kind = 4 ) xadj(n+1)

  do i = 1, nonz
    call adj_set ( adj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests the RCM method.
!
!  Discussion:
!
!    The example involves the following equations:
!
!      2*x1  - x10       = 0
!      2*x2  - x9  - x10 = 0
!      2*x3  - x8  - x9  = 0
!      2*x4  - x7  - x8  = 0
!      2*x5  - x6  - x7  = 0
!      2*x6  - x5        = 11
!      2*x7  - x4  - x5  = 0
!      2*x8  - x3  - x4  = 0
!      2*x9  - x2  - x3  = 0
!      2*x10 - x1  - x2  = 0
!
!    with solution
!
!      x = (1,3,5,7,9,10,8,6,4,2).
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxadj = 3000
  integer ( kind = 4 ), parameter :: maxenv = 3000
  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) adj(maxadj)
  integer ( kind = 4 ) bandw
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) env(maxenv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iband
  integer ( kind = 4 ) env_size
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) nadj
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  real    ( kind = 8 ) rhs(n)
  integer ( kind = 4 ) xadj(n+1)
  integer ( kind = 4 ) xenv(n+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Use the RCM method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The matrix order is N = ', n
!
!  Initialize the permutation vectors.
!
  call i4vec_indicator ( n, perm )
  call i4vec_indicator ( n, perm_inv )
!
!  Store the adjacency information.
!
  call adj_set_4 ( adj, maxadj, nadj, n, xadj )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Number of adjacency entries NADJ = ', nadj
!
!  Display adjacency information
!
  call adj_print ( n, nadj, xadj, adj )
!
!  Determine the initial envelope size.
!
  call adj_env_size ( n, xadj, nadj, adj, perm, perm_inv, env_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial envelope size is ', env_size
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Generate the RCM ordering.
!
  call genrcm ( n, xadj, adj, perm )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The starting node was originally labeled ', perm(n)
!
!  Get inverse ordering
!
  call perm_inverse ( n, perm, perm_inv )
!
!  Print orderings
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I    Perm(I)   InvPerm(I)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(3i6)' ) i, perm(i), perm_inv(i)
  end do
!
!  Determine the reordered envelope size.
!
  call adj_env_size ( n, xadj, nadj, adj, perm, perm_inv, env_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The reordered envelope size is ', env_size
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Compute the envelope.
!
  call fnenv ( n, xadj, adj, perm, perm_inv, xenv, env_size, bandw )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The envelope size is ', env_size
  write ( *, '(a,i12)' ) '  The bandwidth is ', bandw
!
!  Set RHS, DIAG, ENV.
!
  call setsy4 ( diag, env, adj, perm_inv, xadj, maxadj, maxenv, n, rhs, xenv )
!
!  Factor the matrix.
!
  call es_factor ( n, xenv, env, diag, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04 - Fatal error!'
    write ( *, '(a)' ) '  The matrix is not positive definite.'
    return
  end if
!
!  Solve the system.
!
  call el_solve ( n, xenv, env, diag, rhs )

  call eu_solve ( n, xenv, env, diag, rhs )
!
!  Unpermute the solution.
!
  call perm_rv ( n, rhs, perm )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution:'
  write ( *, '(a)' ) ' '
  write ( *, '(5g14.6)' ) rhs(1:n)

  return
end
subroutine setsy4 ( diag, env, adj, perm_inv, xadj, maxadj, maxenv, n, rhs, &
  xenv )

!*****************************************************************************80
!
!! SETSY4 stores the numerical values defining problem 4.
!
!  Discussion:
!
!    There is only one nonzero right hand side entry.
!    The matrix diagonal entries are all 2.
!    The nonzero offdiagonal entries are all -1.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) maxadj
  integer ( kind = 4 ) maxenv
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(maxadj)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) env(maxenv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jsub
  integer ( kind = 4 ) perm_inv(n)
  real    ( kind = 8 ) rhs(n)
  real    ( kind = 8 ) value
  integer ( kind = 4 ) xadj(n+1)
  integer ( kind = 4 ) xenv(n+1)
!
!  Zero out storage.
!
  rhs(1:n) = 0.0D+00
  diag(1:n) = 0.0D+00
  env(1:maxenv) = 0.0D+00
!
!  Set the nonzero elements of the right hand side vector.
!
  isub = 6
  value = 11.0D+00

  call addrhs ( perm_inv, isub, n, rhs, value )
!
!  Set the diagonal entries of the matrix.
!
  diag(1:n) = 2.0D+00
!
!  Set the off diagonal terms of the matrix.
!
  do i = 1, n

    isub = i

    do j = xadj(i), xadj(i+1) - 1

      jsub = adj(j)
      value = -1.0D+00

      call addrcm ( isub, jsub, value, perm_inv, diag, xenv, env, n )

    end do

  end do

  return
end
subroutine adj_set_4 ( adj, maxadj, nadj, n, xadj )

!*****************************************************************************80
!
!! ADJ_SET_4 sets up the adjacency structure for problem 4.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nonz = 19

  integer ( kind = 4 ) maxadj
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(maxadj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save, dimension ( nonz ) :: ilist = (/ &
    -1,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10 /)
  integer ( kind = 4 ), save, dimension ( nonz ) :: jlist = (/ &
    -1,10,10,9,9,8,7,8,7,6,5,4,5,3,4,2,3,1,2 /)
  integer ( kind = 4 ) nadj
  integer ( kind = 4 ) xadj(n+1)

  do i = 1, nonz
    call adj_set ( adj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests the RQT method.
!
!  Discussion:
!
!    This example involves the following equations:
!
!      2*x1  - x10       = 0
!      2*x2  - x9  - x10 = 0
!      2*x3  - x8  - x9  = 0
!      2*x4  - x7  - x8  = 0
!      2*x5  - x6  - x7  = 0
!      2*x6  - x5        = 11
!      2*x7  - x4  - x5  = 0
!      2*x8  - x3  - x4  = 0
!      2*x9  - x2  - x3  = 0
!      2*x10 - x1  - x2  = 0
!
!    with solution
!
!      x = (1,3,5,7,9,10,8,6,4,2).
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxadj = 300
  integer ( kind = 4 ), parameter :: maxblk = 10
  integer ( kind = 4 ), parameter :: maxenv = 300
  integer ( kind = 4 ), parameter :: maxnon = 300
  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) adj(maxadj)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) env(maxenv)
  integer ( kind = 4 ) father(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iband
  integer ( kind = 4 ) env_size
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) nadj
  integer ( kind = 4 ) nblks
  integer ( kind = 4 ) nofnz
  real    ( kind = 8 ) nonz(maxnon)
  integer ( kind = 4 ) nzsub(maxnon)
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  real    ( kind = 8 ) rhs(n)
  integer ( kind = 4 ) xadj(n+1)
  integer ( kind = 4 ) xblk(maxblk+1)
  integer ( kind = 4 ) xenv(n+1)
  integer ( kind = 4 ) xnonz(n+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Use the RQT method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The matrix order is N = ', n
!
!  Initialize the permutation vectors.
!
  call i4vec_indicator ( n, perm )
  call i4vec_indicator ( n, perm_inv )
!
!  Store the adjacency information.
!
  call adj_set_5 ( adj, maxadj, nadj, n, xadj )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Number of adjacency entries NADJ = ', nadj
!
!  Display the adjacency information.
!
  call adj_print ( n, nadj, xadj, adj )
!
!  Determine the initial envelope size.
!
  call adj_env_size ( n, xadj, nadj, adj, perm, perm_inv, env_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial envelope size is ', env_size
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Generate the RQT ordering.
!
  call genrqt ( n, xadj, adj, nblks, xblk, perm )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  After GENRQT, the  number of blocks is ', nblks
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The envelope size is ', env_size

  call block_shuffle ( xadj, adj, perm, nblks, xblk, n )
!
!  Compute the inverse ordering.
!
  call perm_inverse ( n, perm, perm_inv )
!
!  Print orderings
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I    Perm(I)   InvPerm(I)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(3i6)' ) i, perm(i), perm_inv(i)
  end do
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Determine the quotient tree adjacency structure.
!
  call fntadj ( xadj, adj, perm, nblks, xblk, father, n )
!
!  Determine the envelope index vector.
!
  call fntenv ( xadj, adj, perm, perm_inv, nblks, xblk, xenv, env_size, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The reordered envelope size is ', env_size

  nofnz = maxnon
  call fnofnz ( xadj, adj, perm, perm_inv, nblks, xblk, xnonz, nzsub, nofnz, n )
!
!  Set RHS, DIAG, ENV.
!
  call setsy5 ( diag, env, adj, perm_inv, xadj, maxadj, maxenv, n, nonz, &
    nzsub, rhs, xenv, xnonz )
!
!  Factor the system.
!
  call ts_factor ( nblks, xblk, father, diag, xenv, env, xnonz, nonz, nzsub, &
    n, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Fatal error!'
    write ( *, '(a)' ) '  The matrix is not positive definite.'
    return
  end if
!
!  Solve the system.
!
  call ts_solve ( nblks, xblk, diag, xenv, env, xnonz, nonz, nzsub, rhs, n )
!
!  Unpermute the solution.
!
  call perm_rv ( n, rhs, perm )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution:'
  write ( *, '(a)' ) ' '
  write ( *, '(5g14.6)' ) rhs(1:n)

  return
end
subroutine setsy5 ( diag, env, adj, perm_inv, xadj, maxadj, maxenv, n, nonz, &
  nzsub, rhs, xenv, xnonz )

!*****************************************************************************80
!
!! SETSY5 stores the numerical values defining problem 5.
!
!  Discussion:
!
!    There is only one nonzero right hand side entry.
!    The matrix diagonal entries are all 2.
!    The nonzero offdiagonal entries are all -1.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) maxadj
  integer ( kind = 4 ) maxenv
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(maxadj)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) env(maxenv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) j
  integer ( kind = 4 ) xadj(n+1)
  integer ( kind = 4 ) jsub
  real    ( kind = 8 ) nonz(*)
  integer ( kind = 4 ) nzsub(*)
  integer ( kind = 4 ) perm_inv(n)
  real    ( kind = 8 ) rhs(n)
  real    ( kind = 8 ) value
  integer ( kind = 4 ) xenv(n+1)
  integer ( kind = 4 ) xnonz(n+1)
!
!  Zero out storage.
!
  rhs(1:n) = 0.0D+00
  diag(1:n) = 0.0D+00
  env(1:maxenv) = 0.0D+00
!
!  Set the nonzero elements of the right hand side vector.
!
  isub = 6
  value = 11.0D+00

  call addrhs ( perm_inv, isub, n, rhs, value )
!
!  Set the diagonal entries of the matrix.
!
  diag(1:n) = 2.0D+00
!
!  Set the off diagonal terms of the matrix.
!
  do i = 1, n

    isub = i

    do j = xadj(i), xadj(i+1) - 1

      jsub = adj(j)
      value = -1.0D+00

      call addrqt ( isub, jsub, value, perm_inv, diag, xenv, env, xnonz, nonz, &
        nzsub, n )

    end do

  end do

  return
end
subroutine adj_set_5 ( adj, maxadj, nadj, n, xadj )

!*****************************************************************************80
!
!! ADJ_SET_5 sets up the adjacency structure for problem 5.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nonz = 19

  integer ( kind = 4 ) maxadj
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(maxadj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save, dimension ( nonz ) :: ilist = (/ &
    -1,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10 /)
  integer ( kind = 4 ), save, dimension ( nonz ) :: jlist = (/ &
    -1,10,10,9,9,8,7,8,7,6,5,4,5,3,4,2,3,1,2 /)
  integer ( kind = 4 ) nadj
  integer ( kind = 4 ) xadj(n+1)

  do i = 1, nonz
    call adj_set ( adj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests the RCM method.
!
!  Discussion:
!
!    This example corresponds to an implementation of the
!    5 point finite difference operator that approximates
!    the second derivative at a node by a centered difference.
!
!    The 8 by 8 grid is poorly numbered, to demonstrate how
!    SPARSEPAK is able to improve a terrible bandwidth.
!
!    Here is the grid as defined by the user:
!
!      18 02 57 41 11 29 48 21
!      49 38 12 25 45 63 37 01
!      30 62 53 07 60 52 15 56
!      08 46 64 16 34 28 24 40
!      42 22 26 33 17 06 44 10
!      58 13 05 61 27 51 59 32
!      03 35 54 43 23 14 36 47
!      19 50 31 09 39 55 04 20
!
!    The matrix consists of 4's on the diagonal, -1's on the
!    off-diagonals corresponding to the 4 immediate neighbors.
!    the right hand side is 2 on the corners, 1 on the sides,
!    and 0 in the interior.  The solution is (1,1,1,...,1,1,1).
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxadj = 3000
  integer ( kind = 4 ), parameter :: maxenv = 3000
  integer ( kind = 4 ), parameter :: n = 64

  integer ( kind = 4 ) adj(maxadj)
  integer ( kind = 4 ) bandw
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) env(maxenv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iband
  integer ( kind = 4 ) env_size
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) nadj
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  real    ( kind = 8 ) rhs(n)
  integer ( kind = 4 ) xadj(n+1)
  integer ( kind = 4 ) xenv(n+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Use the RCM method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The matrix order is N = ', n
!
!  Initialize the permutation vectors.
!
  call i4vec_indicator ( n, perm )
  call i4vec_indicator ( n, perm_inv )
!
!  Store the adjacency information.
!
  call adj_set_6 ( adj, maxadj, nadj, n, xadj )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Number of adjacency entries NADJ = ', nadj
!
!  Display the adjacency information
!
  call adj_print ( n, nadj, xadj, adj )
!
!  Determine the initial envelope size.
!
  call adj_env_size ( n, xadj, nadj, adj, perm, perm_inv, env_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial envelope size is ', env_size
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Generate the RCM ordering.
!
  call genrcm ( n, xadj, adj, perm )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The starting node was originally labeled ', perm(n)
!
!  Get the inverse ordering.
!
  call perm_inverse ( n, perm, perm_inv )
!
!  Print the orderings.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I    Perm(I)   Perm_Inv(I)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(3i6)' ) i, perm(i), perm_inv(i)
  end do
!
!  Determine the reordered envelope size.
!
  call adj_env_size ( n, xadj, nadj, adj, perm, perm_inv, env_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The reordered envelope size is ', env_size
!
!  Get a picture of the matrix.
!
  call adj_show ( adj, iband, perm_inv, nadj, n, perm, xadj )
!
!  Compute the envelope.
!
  call fnenv ( n, xadj, adj, perm, perm_inv, xenv, env_size, bandw )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The envelope size is ', env_size
  write ( *, '(a,i12)' ) '  The bandwidth is ', bandw
!
!  Set the right hand side and the matrix.
!
  call setsy6 ( diag, env, adj, perm_inv, xadj, maxadj, maxenv, n, rhs, xenv )
!
!  Factor the matrix.
!
  call es_factor ( n, xenv, env, diag, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06 - Fatal error!'
    write ( *, '(a)' ) '  The matrix is not positive definite.'
    return
  end if
!
!  Solve the system.
!
  call el_solve ( n, xenv, env, diag, rhs )

  call eu_solve ( n, xenv, env, diag, rhs )
!
!  Unpermute the solution.
!
  call perm_rv ( n, rhs, perm )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution:'
  write ( *, '(a)' ) ' '
  write ( *, '(5g14.6)' ) rhs(1:n)

  return
end
subroutine setsy6 ( diag, env, adj, perm_inv, xadj, maxadj, maxenv, n, rhs, &
  xenv )

!*****************************************************************************80
!
!! SETSY6 stores the numerical values defining problem 6.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) maxadj
  integer ( kind = 4 ) maxenv
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(maxadj)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) env(maxenv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter, dimension ( 28 ) :: ilist = (/ &
    1, 2, 3, 4, 8, 9,10,11,18,19,20,21,29,30,31,32,39, &
   40,41,42,47,48,49,50,55,56,57,58 /)
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) j
  integer ( kind = 4 ) xadj(n+1)
  integer ( kind = 4 ) jsub
  integer ( kind = 4 ) perm_inv(n)
  real    ( kind = 8 ) rhs(n)
  real    ( kind = 8 ), parameter, dimension ( 28 ) :: rlist = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 2.0D+00, 2.0D+00, &
    2.0D+00, 2.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00 /)
  real    ( kind = 8 ) value
  integer ( kind = 4 ) xenv(n+1)
!
!  Zero out the arrays.
!
  rhs(1:n) = 0.0D+00
  diag(1:n) = 4.0D+00
  env(1:maxenv) = 0.0D+00
!
!  Set the nonzero elements of the right hand side.
!
  do i = 1, 28
    call addrhs ( perm_inv, ilist(i), n, rhs, rlist(i) )
  end do
!
!  Set the off diagonal terms of the matrix.
!
  do i = 1, n

    isub = i

    do j = xadj(i), xadj(i+1) - 1

      jsub = adj(j)
      value = -1.0D+00

      call addrcm ( isub, jsub, value, perm_inv, diag, xenv, env, n )

    end do

  end do

  return
end
subroutine adj_set_6 ( adj, maxadj, nadj, n, xadj )

!*****************************************************************************80
!
!! ADJ_SET_6 sets up the adjacency structure for problem 6.
!
!  Modified:
!
!    24 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) maxadj
  integer ( kind = 4 ) n

  integer ( kind = 4 ) adj(maxadj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter, dimension ( 113 ) :: ilist = (/ &
    -1,13,16,17,18,19,20,21,22,23,24,25,25,26,26,27,27,28,28,29, &
    30,31,32,33,33,33,34,34,34,35,35,36,36,37,37,38,38,39,39,40, &
    40,41,41,42,42,43,43,44,44,44,45,45,46,46,47,47,47,48,48,48, &
    49,49,49,50,50,50,51,51,51,52,52,53,53,54,54,54,54,55,55,55, &
    56,56,56,57,57,57,58,58,58,59,59,59,59,60,60,60,60,61,61,61, &
    61,62,62,62,62,63,63,63,63,64,64,64,64 /)
  integer ( kind = 4 ) xadj(n+1)
  integer ( kind = 4 ), parameter, dimension ( 113 ) :: jlist = (/ &
    -1, 5, 7, 6, 2, 3, 4, 1,13,14,15, 7,12, 5,22,17,23, 6,24,11, &
     8, 9,10,16,17,26,16,17,28, 3,13, 4,14, 1,15, 2,12, 9,23,10, &
    24,11,25, 8,22, 9,23, 6,10,24,11,25, 8,22,20,32,36,21,29,37, &
    18,30,38,19,31,35, 6,14,27,15,28, 7,12, 5,31,35,43, 4,14,39, &
     1,15,40, 2,12,41, 3,13,42,32,36,44,51, 7,34,45,52, 5,27,33, &
    43,30,38,46,53,29,37,45,52,16,26,46,53 /)
  integer ( kind = 4 ) nadj

  do i = 1, 113
    call adj_set ( adj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
  end do

  return
end
