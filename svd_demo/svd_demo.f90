program main

!*****************************************************************************80
!
!! MAIN is the main program for SVD_DEMO.
!
!  Discussion:
!
!    SVD_DEMO demonstrates the SVD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modififed:
!
!    19 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Usage:
!
!    svd_demo m n [seed]
!
!  Command Parameters:
!
!    Command parameter, integer ( kind = 4 ) M, N, the number of rows and 
!    columns of the matrix.  If M or N is not supplied on the command line,
!    the user is prompted to supply them.
!
!    Command parameter, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.  If SEED is not supplied on the command line, a value
!    is generated internally.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Local, real ( kind = 8 ) S(M,N), the diagonal factor
!    in the singular value decomposition of A.
!
!    Output, real ( kind = 8 ) U(M,M), the first orthogonal factor
!    in the singular value decomposition of A.
!
!    Output, real ( kind = 8 ) V(N,N), the second orthogonal factor
!    in the singular value decomposition of A.
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a_pseudo
  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: s
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: s2
  integer ( kind = 4 ) seed
  character ( len = 80 ) string
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: u
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: u2
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v2

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SVD_DEMO'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the singular value decomposition (SVD)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A real MxN matrix A can be factored as:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    A = U * S * V'''
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    U = MxM orthogonal,'
  write ( *, '(a)' ) '    S = MxN zero except for diagonal,'
  write ( *, '(a)' ) '    V = NxN orthogonal.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The diagonal of S contains only nonnegative numbers'
  write ( *, '(a)' ) '  and these are arranged in descending order.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  If at least one command line argument, it's M.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, m, ierror, length )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_DEMO:'
    write ( *, '(a)' ) '  Please enter the value of M.'
    write ( *, '(a)' ) '  (Number of rows in matrix A)'
    write ( *, '(a)' ) '  (We prefer M <= 10!).'

    read ( *, * ) m

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix row order    M = ', m
!
!  If at least two command line argument, it's N.
!
  if ( 2 <= arg_num ) then

    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, n, ierror, length )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_DEMO:'
    write ( *, '(a)' ) '  Please enter the value of N.'
    write ( *, '(a)' ) '  (Number of columns in matrix A)'
    write ( *, '(a)' ) '  (We prefer N <= 10!).'

    read ( *, * ) n

  end if

  write ( *, '(a,i8)' ) '  Matrix column order N = ', n
!
!  If a third command line argument, it's SEED.
!
  if ( 3 <= arg_num ) then

    iarg = 3
    call getarg ( iarg, string )
    call s_to_i4 ( string, seed, ierror, length )
    write ( *, '(a,i12)' ) '  Random number SEED    = ', seed
    write ( *, '(a)' ) '  (Chosen by user.)'

  else

    call get_seed ( seed )
    write ( *, '(a,i12)' ) '  Random number SEED    = ', seed
    write ( *, '(a)' ) '  (Chosen by program.)'

  end if
!
!  Set aside space for the arrays.
!
  allocate ( a(1:m,1:n) )
  allocate ( a_pseudo(1:n,1:m) )
  allocate ( u(1:m,1:m) )
  allocate ( u2(1:m,1:m) )
  allocate ( s(1:m,1:n) )
  allocate ( s2(1:m,1:n) )
  allocate ( v(1:n,1:n) )
  allocate ( v2(1:n,1:n) )
!
!  Generate the matrix A.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We choose a "random" matrix A, with integral'
  write ( *, '(a)' ) '  values between 0 and 10.'

  call r8mat_uniform_01 ( m, n, seed, a )

  a(1:m,1:n) = anint ( 10.0D+00 * a(1:m,1:n) )

  call r8mat_print ( m, n, a, '  The matrix A:' )
!
!  Get the SVD from LAPACK.
!
  call get_svd_lapack ( m, n, a, u, s, v )
!
!  Print the SVD.
!
  call r8mat_print ( m, m, u, '  The orthogonal factor U:' )

  call r8mat_print ( m, n, s, '  The diagonal factor S:' )

  call r8mat_print ( n, n, v, '  The orthogonal factor V:' )
!
!  Check that A = U * S * V'.
!
  call svd_product_test ( m, n, a, u, s, v )
!
!  Compute the norm of the difference between A and the successive
!  sums of rank one approximants.
!
  call rank_one_test ( m, n, a, u, s, v )
!
!  Actually print the sums of rank one approximants.
!
  call rank_one_print_test ( m, n, a, u, s, v )
!
!  Compute the pseudoinverse.
!
  call pseudo_inverse ( m, n, u, s, v, a_pseudo )

  call r8mat_print ( n, m, a_pseudo, '  The pseudoinverse of A:' )
!
!  Test A*A+ = I+, A+*A = I+
!
  call pseudo_product_test ( m, n, a, a_pseudo )
!
!  Demonstrate the use of the pseudoinverse for linear systems.
!
  call pseudo_linear_solve_test ( m, n, a, a_pseudo, seed )
!
!  Get the SVD from LINPACK.
!
  call get_svd_linpack ( m, n, a, u2, s2, v2 )

  call compare_linpack_lapack ( m, n, u, s, v, u2, s2, v2 )
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( a_pseudo )
  deallocate ( s )
  deallocate ( s2 )
  deallocate ( u )
  deallocate ( u2 )
  deallocate ( v )
  deallocate ( v2 )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SVD_DEMO:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine compare_linpack_lapack ( m, n, u, s, v, u2, s2, v2 )

!*****************************************************************************80
!
!! COMPARE_LINPACK_LAPACK compares the SVD's from LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A, computed by LAPACK.
!
!    Input, real ( kind = 8 ) U2(M,M), S2(M,N), V2(N,N), the factors
!    that form the singular value decomposition of A, computed by LINPACK.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) r8mat_dif_fro
  real ( kind = 8 ) s(m,n)
  real ( kind = 8 ) s2(m,n)
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) u2(m,m)
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) v2(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMPARE_LINPACK_LAPACK:'
  write ( *, '(a)' ) '  While the singular values should be identical,'
  write ( *, '(a)' ) '  the orthogonal factors may have some differences.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Maximum differences:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  U(LAPACK) - U(LINPACK): ', &
    r8mat_dif_fro ( m, m, u, u2 )
  write ( *, '(a,g14.6)' ) '  S(LAPACK) - S(LINPACK): ', &
    r8mat_dif_fro ( m, n, s, s2 )
  write ( *, '(a,g14.6)' ) '  V(LAPACK) - V(LINPACK): ', &
    r8mat_dif_fro ( n, n, v, v2 )

  return
end
subroutine get_seed ( seed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
!
  implicit none

  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 ) today
  integer ( kind = 4 ) values(8)
  character ( len = 5 ) zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) /  11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) /  30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) /  23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp                                    /   6.0D+00

  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( huge ( 1 ), kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == huge ( 1 ) ) then
    seed = seed - 1
  end if

  return
end
subroutine get_svd_lapack ( m, n, a, u, s, v )

!*****************************************************************************80
!
!! GET_SVD_LAPACK gets the SVD of a matrix using a call to LAPACK.
!
!  Discussion:
!
!    The singular value decomposition of a real MxN matrix A has the form:
!
!      A = U * S * V'
!
!    where
!
!      U is MxM orthogonal,
!      S is MxN, and entirely zero except for the diagonal;
!      V is NxN orthogonal.
!
!    Moreover, the nonzero entries of S are positive, and appear
!    in order, from largest magnitude to smallest.
!
!    This routine calls the LAPACK routine DGESVD to compute the
!    factorization.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Output, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_copy(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  character jobu
  character jobv
  integer ( kind = 4 ) lwork
  real ( kind = 8 ) sdiag(min(m,n))
  real ( kind = 8 ) s(m,n)
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ), allocatable, dimension ( : ) :: work

  lwork = max ( 3 * min ( m, n ) + max ( m, n ), 5 * min ( m, n ) )

  allocate ( work(1:lwork) )
!
!  Compute the eigenvalues and eigenvectors.
!
  jobu = 'A'
  jobv = 'A'
  lda = m
  ldu = m
  ldv = n
!
!  The input matrix is destroyed by the routine.  Since we need to keep
!  it around, we only pass a copy to the routine.
!
  a_copy(1:m,1:n) = a(1:m,1:n)

  call dgesvd ( jobu, jobv, m, n, a_copy, lda, sdiag, u, ldu, v, ldv, work, &
    lwork, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GET_SVD_LAPACK - Failure!'
    write ( *, '(a)' ) '  The SVD could not be calculated.'
    write ( *, '(a)' ) '  LAPACK routine DGESVD returned a nonzero'
    write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
    return
  end if
!
!  Make the MxN matrix S from the diagonal values in SDIAG.
!
  s(1:m,1:n) = 0.0D+00
  do i = 1, min ( m, n )
    s(i,i) = sdiag(i)
  end do
!
!  Transpose V.
!
  v = transpose ( v )

  deallocate ( work )

  return
end
subroutine get_svd_linpack ( m, n, a, u, s, v )

!*****************************************************************************80
!
!! GET_SVD_LINPACK gets the SVD of a matrix using a call to LINPACK.
!
!  Discussion:
!
!    The singular value decomposition of a real MxN matrix A has the form:
!
!      A = U * S * V'
!
!    where
!
!      U is MxM orthogonal,
!      S is MxN, and entirely zero except for the diagonal;
!      V is NxN orthogonal.
!
!    Moreover, the nonzero entries of S are positive, and appear
!    in order, from largest magnitude to smallest.
!
!    This routine calls the LINPACK routine DSVDC to compute the
!    factorization.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Output, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_copy(m,n)
  real ( kind = 8 ) e(max(m+1,n))
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lwork
  real ( kind = 8 ) s(m,n)
  real ( kind = 8 ) sdiag(max(m+1,n))
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ), allocatable, dimension ( : ) :: work

  allocate ( work(1:m) )
!
!  Compute the eigenvalues and eigenvectors.
!
  job = 11
  lda = m
  ldu = m
  ldv = n
!
!  The input matrix is destroyed by the routine.  Since we need to keep
!  it around, we only pass a copy to the routine.
!
  a_copy(1:m,1:n) = a(1:m,1:n)

  call dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GET_SVD_LINPACK - Failure!'
    write ( *, '(a)' ) '  The SVD could not be calculated.'
    write ( *, '(a)' ) '  LINPACK routine DSVDC returned a nonzero'
    write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
    return
  end if
!
!  Make the MxN matrix S from the diagonal values in SDIAG.
!
  s(1:m,1:n) = 0.0D+00
  do i = 1, min ( m, n )
    s(i,i) = sdiag(i)
  end do
!
!  Note that we do NOT need to transpose the V that comes out of LINPACK!
!
  deallocate ( work )

  return
end
subroutine pseudo_inverse ( m, n, u, s, v, a_pseudo )

!*****************************************************************************80
!
!! PSEUDO_INVERSE computes the pseudoinverse.
!
!  Discussion:
!
!    Given the singular value decomposition of a real MxN matrix A:
!
!      A = U * S * V'
!
!    where
!
!      U is MxM orthogonal,
!      S is MxN, and entirely zero except for the diagonal;
!      V is NxN orthogonal.
!
!    the pseudo inverse is the NxM matrix A+ with the form
!
!      A+ = V * S+ * U'
!
!    where
!
!      S+ is the NxM matrix whose nonzero diagonal elements are
!      the inverses of the corresponding diagonal elements of S.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Input, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A.
!
!    Output, real ( kind = 8 ) A_PSEUDO(N,M), the pseudo_inverse of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_pseudo(n,m)
  integer ( kind = 4 ) i
  real ( kind = 8 ) s(m,n)
  real ( kind = 8 ) sp(n,m)
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) v(n,n)

  sp(1:n,1:m) = 0.0D+00
  do i = 1, min ( m, n )
    if ( s(i,i) /= 0.0D+00 ) then
      sp(i,i) = 1.0D+00 / s(i,i)
    end if
  end do

  a_pseudo(1:n,1:m) = matmul ( v(1:n,1:n), &
    matmul ( sp(1:n,1:m), transpose ( u(1:m,1:m) ) ) )

  return
end
subroutine pseudo_linear_solve_test ( m, n, a, a_pseudo, seed )

!*****************************************************************************80
!
!! PSEUDO_LINEAR_SOLVE_TEST uses the pseudoinverse for linear systems.
!
!  Discussion:
!
!    Given an MxN matrix A, and its pseudoinverse A+:
!
!      "Solve" A  * x = b by x = A+  * b.
!
!      "Solve" A' * x = b by x = A+' * b.
!
!    When the system is overdetermined, the solution minimizes the
!    L2 norm of the residual.
!
!    When the system is underdetermined, the solution
!    is the solution of minimum L2 norm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Input, real ( kind = 8 ) A_PSEUDO(N,M), the pseudo_inverse of A.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_pseudo(n,m)
  real ( kind = 8 ) bm(m)
  real ( kind = 8 ) bn(n)
  real ( kind = 8 ) r8vec_norm_l2
  real ( kind = 8 ) rm(m)
  real ( kind = 8 ) rn(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) xm1(m)
  real ( kind = 8 ) xm2(m)
  real ( kind = 8 ) xn1(n)
  real ( kind = 8 ) xn2(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PSEUDO_LINEAR_SOLVE_TEST'
!
!  A * x = b, b in range of A.
!
  call r8vec_uniform_01 ( n, seed, xn1 )
  xn1(1:n) = anint ( 10.0D+00 * xn1(1:n) )
  bm(1:m) = matmul ( a(1:m,1:n), xn1(1:n) )
  xn2(1:n) = matmul ( a_pseudo(1:n,1:m), bm(1:m) )
  rm(1:m) = matmul ( a(1:m,1:n), xn2(1:n) ) - bm(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Given:'
  write ( *, '(a)' ) '    b = A * x1'
  write ( *, '(a)' ) '  so that b is in the range of A, solve'
  write ( *, '(a)' ) '    A * x = b'
  write ( *, '(a)' ) '  using the pseudoinverse:'
  write ( *, '(a)' ) '    x2 = A+ * b.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Norm of x1 = ', r8vec_norm_l2 ( n, xn1 )
  write ( *, '(a,g14.6)' ) '  Norm of x2 = ', r8vec_norm_l2 ( n, xn2 )
  write ( *, '(a,g14.6)' ) '  Norm of residual = ', r8vec_norm_l2 ( m, rm )
!
!  A * x = b, b not in range of A.
!
  if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  For N < M, most systems A*x=b will not be'
    write ( *, '(a)' ) '  exactly and uniquely solvable, except in the '
    write ( *, '(a)' ) '  least squares sense.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Here is an example:'

    call r8vec_uniform_01 ( m, seed, bm )
    xn2(1:n) = matmul ( a_pseudo(1:n,1:m), bm(1:m) )
    rm(1:m) = matmul ( a(1:m,1:n), xn2(1:n) ) - bm(1:m)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Given b is NOT in the range of A, solve'
    write ( *, '(a)' ) '    A * x = b'
    write ( *, '(a)' ) '  using the pseudoinverse:'
    write ( *, '(a)' ) '    x2 = A+ * b.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Norm of x2 = ', r8vec_norm_l2 ( n, xn2 )
    write ( *, '(a,g14.6)' ) '  Norm of residual = ', r8vec_norm_l2 ( m, rm )
  end if
!
!  A' * x = b, b is in the range of A'.
!
  call r8vec_uniform_01 ( m, seed, xm1 )
  xm1(1:m) = anint ( 10.0D+00 * xm1(1:m) )
  bn(1:n) = matmul ( transpose ( a(1:m,1:n) ), xm1(1:m) )
  xm2(1:m) = matmul ( transpose ( a_pseudo(1:n,1:m) ), bn(1:n) )
  rn(1:n) = matmul ( transpose ( a(1:m,1:n) ), xm2(1:m) ) - bn(1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Given:'
  write ( *, '(a)' ) '    b = A'' * x1'
  write ( *, '(a)' ) '  so that b is in the range of A'', solve'
  write ( *, '(a)' ) '    A'' * x = b'
  write ( *, '(a)' ) '  using the pseudoinverse:'
  write ( *, '(a)' ) '    x2 = A+'' * b.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Norm of x1 = ', r8vec_norm_l2 ( m, xm1 )
  write ( *, '(a,g14.6)' ) '  Norm of x2 = ', r8vec_norm_l2 ( m, xm2 )
  write ( *, '(a,g14.6)' ) '  Norm of residual = ', r8vec_norm_l2 ( n, rn )
!
!  A' * x = b, b is not in the range of A'.

  if ( m < n ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  For M < N, most systems A''*x=b will not be'
    write ( *, '(a)' ) '  exactly and uniquely solvable, except in the'
    write ( *, '(a)' ) '  least squares sense.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Here is an example:'

    call r8vec_uniform_01 ( n, seed, bn )
    xm2(1:m) = matmul ( transpose ( a_pseudo(1:n,1:m) ), bn(1:n) )
    rn(1:n) = matmul ( transpose ( a(1:m,1:n) ), xm2(1:m) ) - bn(1:n)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Given b is NOT in the range of A'', solve'
    write ( *, '(a)' ) '    A'' * x = b'
    write ( *, '(a)' ) '  using the pseudoinverse:'
    write ( *, '(a)' ) '    x2 = A+'' * b.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Norm of x2 = ', r8vec_norm_l2 ( m, xm2 )
    write ( *, '(a,g14.6)' ) '  Norm of residual = ', r8vec_norm_l2 ( n, rn )

  end if

  return
end
subroutine pseudo_product_test ( m, n, a, a_pseudo )

!*****************************************************************************80
!
!! PSEUDO_PRODUCT_TEST examines pseudoinverse products.
!
!  Discussion:
!
!    Given an MxN matrix A, and its pseudoinverse A+, we must have
!
!      A+ * A * A+ = A+
!      A * A+ * A = A
!      ( A * A+ )' = A * A+ (MxM symmetry)
!      ( A+ * A )' = A+ * A (NxN symmetry)
!
!    If M <= N, A * A+ may be "interesting" (equal to or "like" the identity),
!    if N <= M, A+ * A may be "interesting" (equal to or "like" the identity).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Input, real ( kind = 8 ) A_PSEUDO(N,M), the pseudo_inverse of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_pseudo(n,m)
  real ( kind = 8 ) bmm(m,m)
  real ( kind = 8 ) bmn(m,n)
  real ( kind = 8 ) bnm(n,m)
  real ( kind = 8 ) bnn(n,n)
  real ( kind = 8 ) dif1
  real ( kind = 8 ) dif2
  real ( kind = 8 ) dif3
  real ( kind = 8 ) dif4
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8mat_dif_fro

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PSEUDO_PRODUCT_TEST'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The following relations MUST hold:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   A  * A+ * A  = A'
  write ( *, '(a)' ) '   A+ * A  * A+ = A+'
  write ( *, '(a)' ) ' ( A  * A+ ) is MxM symmetric;'
  write ( *, '(a)' ) ' ( A+ * A  ) is NxN symmetric'

  bmn(1:m,1:n) = matmul ( a(1:m,1:n), &
                 matmul ( a_pseudo(1:n,1:m), a(1:m,1:n) ) )

  dif1 = r8mat_dif_fro ( m, n, a, bmn )

  bnm(1:n,1:m) = matmul ( a_pseudo(1:n,1:m), &
                 matmul ( a(1:m,1:n), a_pseudo(1:n,1:m) ) )

  dif2 = r8mat_dif_fro ( n, m, a_pseudo, bnm )

  bmm(1:m,1:m) = matmul ( a(1:m,1:n), a_pseudo(1:n,1:m) )

  dif3 = r8mat_dif_fro ( m, m, bmm, transpose ( bmm ) )

  bnn(1:n,1:n) = matmul ( a_pseudo(1:n,1:m), a(1:m,1:n) )

  dif4 = r8mat_dif_fro ( n, n, bnn, transpose ( bnn ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here are the Frobenius norms of the errors'
  write ( *, '(a)' ) '  in these relationships:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '   A  * A+ * A  = A            ', dif1
  write ( *, '(a,g14.6)' ) '   A+ * A  * A+ = A+           ', dif2
  write ( *, '(a,g14.6)' ) ' ( A  * A+ ) is MxM symmetric; ', dif3
  write ( *, '(a,g14.6)' ) ' ( A+ * A  ) is NxN symmetric; ', dif4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In some cases, the matrix A * A+'
  write ( *, '(a)' ) '  may be interesting (if M <= N, then'
  write ( *, '(a)' ) '  it MIGHT look like the identity.)'
  write ( *, '(a)' ) ' '

  call r8mat_print ( m, m, bmm, '  A * A+:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In some cases, the matrix A+ * A'
  write ( *, '(a)' ) '  may be interesting (if N <= M, then'
  write ( *, '(a)' ) '  it MIGHT look like the identity.)'
  write ( *, '(a)' ) ' '

  call r8mat_print ( n, n, bnn, '  A+ * A:' )

  return
end
function r8mat_dif_fro ( m, n, a, b )

!*****************************************************************************80
!
!! R8MAT_DIF_FRO returns the Frobenius norm of the difference of two R8MAT's.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = 8 ) values.
!
!    The Frobenius norm is defined as
!
!      R8MAT_DIF_FRO = sqrt (
!        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) ( A(I,J) - B(I,J) )^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A and B.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A and B.
!
!    Input, real ( kind = 8 ) A(M,N), B(M,N), the matrices
!    for which we want the Frobenius norm of the difference.
!
!    Output, real ( kind = 8 ) R8MAT_DIF_FRO, the Frobenius norm of
!    the difference of A and B.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m,n)
  real ( kind = 8 ) r8mat_dif_fro

  r8mat_dif_fro = sqrt ( sum ( ( a(1:m,1:n) - b(1:m,1:n) )**2 ) )

  return
end
function r8mat_norm_fro ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = 8 ) values.
!
!    The Frobenius norm is defined as
!
!      R8MAT_NORM_FRO = sqrt (
!        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
!
!    The matrix Frobenius norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2000
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
!    Input, real ( kind = 8 ) A(M,N), the matrix whose Frobenius
!    norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_FRO, the Frobenius norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8mat_norm_fro

  r8mat_norm_fro = sqrt ( sum ( a(1:m,1:n)**2 ) )

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = 8 ) values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
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
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = 8 ) values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2005
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

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = 8 ) values.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
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
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
function r8vec_norm_l2 ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_L2, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm_l2

  r8vec_norm_l2 = sqrt ( sum ( a(1:n)**2 ) )

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
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
!    Input, integer ( kind = 4 ) N, the number of entries
!    in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value,
!    which should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine rank_one_print_test ( m, n, a, u, s, v )

!*****************************************************************************80
!
!! RANK_ONE_PRINT_TEST prints the sums of rank one matrices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Input, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_norm
  real ( kind = 8 ) dif_norm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8mat_dif_fro
  real ( kind = 8 ) r8mat_norm_fro
  real ( kind = 8 ) s(m,n)
  character ( len = 80 ) title
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) usv(m,n)
  real ( kind = 8 ) v(n,n)

  a_norm = r8mat_norm_fro ( m, n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANK_ONE_PRINT_TEST:'
  write ( *, '(a)' ) '  Print the sums of R rank one matrices.'

  do r = 0, min ( m, n )

    usv(1:m,1:n) = matmul ( u(1:m,1:r), &
                   matmul ( s(1:r,1:r), transpose ( v(1:n,1:r) ) ) )

    write ( title, '(a,i8)' ) '  Rank R = ', r
    call r8mat_print ( m, n, usv, title )

  end do

  call r8mat_print ( m, n, a, '  Original matrix A:' )

  return
end
subroutine rank_one_test ( m, n, a, u, s, v )

!*****************************************************************************80
!
!! RANK_ONE_TEST compares A to the sum of rank one matrices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Input, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_norm
  real ( kind = 8 ) dif_norm
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8mat_dif_fro
  real ( kind = 8 ) r8mat_norm_fro
  real ( kind = 8 ) s(m,n)
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) usv(m,n)
  real ( kind = 8 ) v(n,n)

  a_norm = r8mat_norm_fro ( m, n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANK_ONE_TEST:'
  write ( *, '(a)' ) '  Compare A to the sum of R rank one matrices.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         R    Absolute      Relative'
  write ( *, '(a)' ) '              Error         Error'
  write ( *, '(a)' ) ' '

  do r = 0, min ( m, n )

    usv(1:m,1:n) = matmul ( u(1:m,1:r), &
                   matmul ( s(1:r,1:r), transpose ( v(1:n,1:r) ) ) )

    dif_norm = r8mat_dif_fro ( m, n, a, usv )

    write ( *, '(2x,i8,2x,g14.6,2x, g14.6)' ) r, dif_norm, dif_norm / a_norm

  end do

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
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2006
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
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S
!    used to make the integer.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * ) s
  integer ( kind = 4 ) state
  integer ( kind = 4 ) value

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
subroutine svd_product_test ( m, n, a, u, s, v )

!*****************************************************************************80
!
!! SVD_PRODUCT_TEST tests that A = U * S * V'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Input, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_norm
  real ( kind = 8 ) dif_norm
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8mat_dif_fro
  real ( kind = 8 ) r8mat_norm_fro
  real ( kind = 8 ) s(m,n)
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) usv(m,n)
  real ( kind = 8 ) v(n,n)

  a_norm = r8mat_norm_fro ( m, n, a )

  usv(1:m,1:n) = matmul ( u(1:m,1:m), &
                 matmul ( s(1:m,1:n), transpose ( v(1:n,1:n) ) ) )

  call r8mat_print ( m, n, usv, '  The product U * S * V'':' )

  dif_norm = r8mat_dif_fro ( m, n, a, usv )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Frobenius Norm of A, A_NORM = ', a_norm
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ABSOLUTE ERROR for A = U*S*V'':'
  write ( *, '(a,g14.6)' ) '  Frobenius norm of difference A-U*S*V'' = ', &
    dif_norm
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RELATIVE ERROR for A = U*S*V'':'
  write ( *, '(a,g14.6)' ) '  Ratio of DIF_NORM / A_NORM = ', &
    dif_norm / a_norm

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
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
