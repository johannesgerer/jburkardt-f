program main

!*****************************************************************************80
!
!! MAIN is the main program for the matrix multiplication example.
!
!  Discussion:
!
!    This program is part of a demonstration of how a FORTRAN90
!    program can interact with MATLAB.
!
!    In this example, the FORTRAN90 program generates matrices A
!    and B, writes them to a file, asks MATLAB to multiply them,
!    and reads the result back from another file.
!
!    The interaction with MATLAB is done using the (nonstandard)
!    SYSTEM call, which is available with certain FORTRAN compilers,
!    including GNU G95 and IBM XLF.
!
!  Modified:
!
!    02 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  character ( len = 80 ) a_file_name
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: b
  character ( len = 80 ) b_file_name
  real    ( kind = 8 ) beta
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: c
  character ( len = 80 ) c_file_name
  character ( len = 120 ) command
  real    ( kind = 8 ) error_frobenius
  integer ( kind = 4 ) n
  integer ( kind = 4 ) result

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATRIX_MULTIPLY:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program:'
  write ( *, '(a)' ) '  * generates matrices A and B;'
  write ( *, '(a)' ) '  * writes them to files;'
  write ( *, '(a)' ) '  * calls MATLAB (via the SYSTEM routine) which'
  write ( *, '(a)' ) '    ** reads A and B from files;'
  write ( *, '(a)' ) '    ** computes C = A * B;'
  write ( *, '(a)' ) '    ** writes C to a file;'
  write ( *, '(a)' ) '  * reads C from the file;'
  write ( *, '(a)' ) '  * reports the success of the computation.'
!
!  Matrix generation.
!
  n = 10
  alpha = 2.0D+00
  beta = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,g14.6)' ) '  Matrix parameter ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  Matrix parameter BETA =  ', beta

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )

  call combin ( alpha, beta, n, a )
  call combin_inverse ( alpha, beta, n, b )
!
!  Write the matrices to files.
!
  a_file_name = 'a.txt'
  call r8mat_write ( a_file_name, n, n, a )

  b_file_name = 'b.txt'
  call r8mat_write ( b_file_name, n, n, b )

  deallocate ( a )
  deallocate ( b )
!
!  Call MATLAB to multiply the matrices.
!
  command = '/usr/local/bin/matlab -nosplash -nodisplay ' // &
    '< matrix_multiply.m > matrix_multiply_matlab_output.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Issued call to MATLAB.'

  call system ( command, result )

  if ( result /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MATRIX_MULTIPLY:'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Failure of command issued through SYSTEM.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The command was:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "' // trim ( command ) // '".'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error return from command = ', result

    write ( *, '(a)' ) ' '
    call timestamp ( )
    stop

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  MATLAB sesion terminated normally.'
!
!  Read the result matrix from the file.
!
  allocate ( c(1:n,1:n) )
  c_file_name = 'c.txt'

  call r8mat_read ( c_file_name, n, n, c )
!
!  Determine ||A*B-I|| which should be zero.
!
  call r8mat_is_identity ( n, c, error_frobenius )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Frobenius norm of A * B - I is ', error_frobenius

  deallocate ( c )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATRIX_MULTIPLY:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine combin ( alpha, beta, n, a )

!*****************************************************************************80
!
!! COMBIN returns the COMBIN matrix.
!
!  Discussion:
!
!    This matrix is known as the combinatorial matrix.
!
!  Formula:
!
!    If ( I = J ) then
!      A(I,J) = ALPHA + BETA
!    else
!      A(I,J) = BETA
!
!  Example:
!
!    N = 5, ALPHA = 2, BETA = 3
!
!    5 3 3 3 3
!    3 5 3 3 3
!    3 3 5 3 3
!    3 3 3 5 3
!    3 3 3 3 5
!
!  Properties:
!
!    A is symmetric: A' = A.
!
!    Because A is symmetric, it is normal.
!
!    Because A is normal, it is diagonalizable.
!
!    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
!
!    A is a circulant matrix: each row is shifted once to get the next row.
!
!    det ( A ) = ALPHA**(N-1) * ( ALPHA + N * BETA ).
!
!    A has constant row sums.
!
!    Because A has constant row sums,
!    it has an eigenvalue with this value,
!    and a (right) eigenvector of ( 1, 1, 1, ..., 1 ).
!
!    A has constant column sums.
!
!    Because A has constant column sums,
!    it has an eigenvalue with this value,
!    and a (left) eigenvector of ( 1, 1, 1, ..., 1 ).
!
!    LAMBDA(1:N-1) = ALPHA,
!    LAMBDA(N) = ALPHA + N * BETA.
!
!    The eigenvector associated with LAMBDA(N) is (1,1,1,...,1)/sqrt(N).
!
!    The other N-1 eigenvectors are simply any (orthonormal) basis
!    for the space perpendicular to (1,1,1,...,1).
!
!    A is nonsingular if ALPHA /= 0 and ALPHA + N * BETA /= 0.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Gregory, David Karney,
!    Example 3.25,
!    A Collection of Matrices for Testing Computational Algorithms,
!    Wiley, 1969, page 53, QA263 G862.
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1, Fundamental Algorithms, Second Edition,
!    Addison-Wesley, Reading, Massachusetts, 1973, page 36.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, BETA, scalars that define A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) i

  a(1:n,1:n) = beta

  do i = 1, n
    a(i,i) = a(i,i) + alpha;
  end do

  return
end
subroutine combin_inverse ( alpha, beta, n, a )

!*****************************************************************************80
!
!! COMBIN_INVERSE returns the inverse of the COMBIN matrix.
!
!  Formula:
!
!    if ( I = J )
!      A(I,J) = (ALPHA+(N-1)*BETA) / (ALPHA*(ALPHA+N*BETA))
!    else
!      A(I,J) =             - BETA / (ALPHA*(ALPHA+N*BETA))
!
!  Example:
!
!    N = 5, ALPHA = 2, BETA = 3
!
!           14 -3 -3 -3 -3
!           -3 14 -3 -3 -3
!   1/34 *  -3 -3 14 -3 -3
!           -3 -3 -3 14 -3
!           -3 -3 -3 -3 14
!
!  Properties:
!
!    A is symmetric: A' = A.
!
!    Because A is symmetric, it is normal.
!
!    Because A is normal, it is diagonalizable.
!
!    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
!
!    A is a circulant matrix: each row is shifted once to get the next row.
!
!    A is Toeplitz: constant along diagonals.
!
!    det ( A ) = 1 / (ALPHA**(N-1) * (ALPHA+N*BETA)).
!
!    A is well defined if ALPHA /= 0D+00 and ALPHA+N*BETA /= 0.
!
!    A is also a combinatorial matrix.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1, Fundamental Algorithms, Second Edition,
!    Addison-Wesley, Reading, Massachusetts, 1973, page 36.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, BETA, scalars that define the matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  real    ( kind = 8 ) bot
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( alpha == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMBIN_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  The entries of the matrix are undefined '
    write ( *, '(a)' ) '  because ALPHA = 0.'
    stop
  else if ( alpha + n * beta == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMBIN_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  The entries of the matrix are undefined '
    write ( *, '(a)' ) '  because ALPHA+N*BETA is zero.'
    stop
  end if

  bot = alpha * ( alpha + real ( n, kind = 8 ) * beta )

  do i = 1, n
    do j = 1, n

      if ( i == j ) then
        a(i,j) = ( alpha + real ( n - 1, kind = 8 ) * beta ) / bot
      else
        a(i,j) = - beta / bot
      end if

    end do
  end do

  return
end
subroutine r8mat_is_identity ( n, a, error_frobenius )

!*****************************************************************************80
!
!! R8MAT_IS_IDENTTITY determines if a matrix is the identity.
!
!  Discussion:
!
!    The routine returns the Frobenius norm of A - I.
!
!  Modified:
!
!    02 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) ERROR_FROBENIUS, the Frobenius norm
!    of the difference matrix A - I, which would be exactly zero
!    if A were the identity matrix.
!
  implicit none

  integer n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) error_frobenius
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) value

  error_frobenius = 0.0D+00

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        error_frobenius = error_frobenius + ( a(i,j) - 1.0D+0 )**2
      else
        error_frobenius = error_frobenius + a(i,j)**2
      end if
    end do
  end do

  error_frobenius = sqrt ( error_frobenius )

  return
end
subroutine r8mat_read ( file_name, m, n, a )

!*****************************************************************************80
!
!! R8MAT_READ reads an R8MAT from a file.
!
!  Modified:
!
!    02 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  character ( len = * ) file_name
  integer ( kind = 4 ), parameter :: file_unit = 1
  integer ( kind = 4 ) i

  open ( unit = file_unit, file = file_name, status = 'old' )

  do i = 1, m
    read ( file_unit, * ) a(i,1:n)
  end do

  close ( unit = file_unit )

  return
end
subroutine r8mat_write ( file_name, m, n, a )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT to a file.
!
!  Modified:
!
!    02 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  character ( len = * ) file_name
  integer ( kind = 4 ), parameter :: file_unit = 1
  integer ( kind = 4 ) i
  character ( len = 40 ) string

  write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', n, 'g', 14, '.', 6, ')'

  open ( unit = file_unit, file = file_name, status = 'replace' )

  do i = 1, m
    write ( file_unit, string ) a(i,1:n)
  end do

  close ( unit = file_unit )

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
