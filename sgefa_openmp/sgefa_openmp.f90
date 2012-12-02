program main

!*****************************************************************************80
!
!! MAIN is the main program for the SGEFA_OPENMP test program.
!
!  Discussion:
!
!    We want to compare methods of solving the linear system A*x=b.
!
!    The first way uses the standard sequential algorithm "SGEFA".
!
!    The second way uses a variant of SGEFA that has been modified to
!    take advantage of OpenMP.
!
!    The third way reruns the variant code, but with OpenMP turned off.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) proc_num
  integer ( kind = 4 ) thread_num

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGEFA_OPENMP'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'

  proc_num = omp_get_num_procs ( )
  thread_num = omp_get_max_threads ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of processors available = ', proc_num
  write ( *, '(a,i8)' ) '  The number of threads available =    ', thread_num
  

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Algorithm        Mode          N    Error       Time'

  write ( *, '(a)' ) ' '
  n = 10
  call test_sgefa ( n )
  call test_sgefa_c ( n )
  call test_sgefa_c_omp ( n )
  call test_sgefa_r ( n )
  call test_sgefa_r_omp ( n )

  write ( *, '(a)' ) ' '
  n = 100
  call test_sgefa ( n )
  call test_sgefa_c ( n )
  call test_sgefa_c_omp ( n )
  call test_sgefa_r ( n )
  call test_sgefa_r_omp ( n )

  write ( *, '(a)' ) ' '
  n = 1000
  call test_sgefa ( n )
  call test_sgefa_c ( n )
  call test_sgefa_c_omp ( n )
  call test_sgefa_r ( n )
  call test_sgefa_r_omp ( n )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGEFA_OPENMP'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test_sgefa ( n )

!*****************************************************************************80
!
!! TEST_SGEFA tests the original version of SGEFA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  real ( kind = 4 ) err
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ), allocatable, dimension ( : ) :: ipvt
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
  real ( kind = 8 ) wtime
  real ( kind = 4 ), allocatable, dimension ( : ) :: x
!
!  Generate the linear system A * x = b.
!
  lda = n

  allocate ( a(1:lda,1:n) )
  allocate ( b(1:n) )
  allocate ( ipvt(1:n) )
  allocate ( x(1:n) )

  call matgen ( lda, n, a, x, b )
!
!  Factor the linear system.
!
  wtime = omp_get_wtime ( )
  call sgefa ( a, lda, n, ipvt, info )
  wtime = omp_get_wtime ( ) - wtime

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST_SGEFA - Fatal error!'
    write ( *, '(a)' ) '  SGEFA reports the matrix is singular.'
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call sgesl ( a, lda, n, ipvt, b, job )

  err = sum ( abs ( x(1:n) - b(1:n) ) )

  write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) &
    '  SGEFA     Sequential   ', n, err, wtime

  deallocate ( a )
  deallocate ( b )
  deallocate ( ipvt )
  deallocate ( x )

  return
end
subroutine test_sgefa_c ( n )

!*****************************************************************************80
!
!! TEST_SGEFA_C tests a column oriented version of SGEFA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  real ( kind = 4 ) err
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ), allocatable, dimension ( : ) :: ipvt
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
  real ( kind = 8 ) wtime
  real ( kind = 4 ), allocatable, dimension ( : ) :: x
!
!  Generate the linear system A * x = b.
!
  lda = n

  allocate ( a(1:lda,1:n) )
  allocate ( b(1:n) )
  allocate ( ipvt(1:n) )
  allocate ( x(1:n) )

  call matgen ( lda, n, a, x, b )
!
!  Factor the linear system.
!
  wtime = omp_get_wtime ( )
  call sgefa_c ( a, lda, n, ipvt, info )
  wtime = omp_get_wtime ( ) - wtime

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST_SGEFA_C - Fatal error!'
    write ( *, '(a)' ) '  SGEFA_C reports the matrix is singular.'
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call sgesl ( a, lda, n, ipvt, b, job )

  err = sum ( abs ( x(1:n) - b(1:n) ) )

  write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) &
    '  SGEFA_C   Sequential   ', n, err, wtime

  deallocate ( a )
  deallocate ( b )
  deallocate ( ipvt )
  deallocate ( x )

  return
end
subroutine test_sgefa_c_omp ( n )

!*****************************************************************************80
!
!! TEST_SGEFA_C_OMP tests a column oriented version of SGEFA + OpenMP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  real ( kind = 4 ) err
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer, allocatable, dimension ( : ) :: ipvt
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
  real ( kind = 8 ) wtime
  real ( kind = 4 ), allocatable, dimension ( : ) :: x
!
!  Generate the linear system A * x = b.
!
  lda = n

  allocate ( a(1:lda,1:n) )
  allocate ( b(1:n) )
  allocate ( ipvt(1:n) )
  allocate ( x(1:n) )

  call matgen ( lda, n, a, x, b )
!
!  Factor the linear system.
!
  wtime = omp_get_wtime ( )
  call sgefa_c_omp ( a, lda, n, ipvt, info )
  wtime = omp_get_wtime ( ) - wtime

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST_SGEFA_C_OMP - Fatal error!'
    write ( *, '(a)' ) '  SGEFA_C_OMP reports the matrix is singular.'
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call sgesl ( a, lda, n, ipvt, b, job )

  err = sum ( abs ( x(1:n) - b(1:n) ) )

  write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) &
    '  SGEFA_C     Parallel   ', n, err, wtime

  deallocate ( a )
  deallocate ( b )
  deallocate ( ipvt )
  deallocate ( x )

  return
end
subroutine test_sgefa_r ( n )

!*****************************************************************************80
!
!! TEST_SGEFA_R tests a row oriented version of SGEFA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  real ( kind = 4 ) err
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ), allocatable, dimension ( : ) :: ipvt
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
  real ( kind = 8 ) wtime
  real ( kind = 4 ), allocatable, dimension ( : ) :: x
!
!  Generate the linear system A * x = b.
!
  lda = n

  allocate ( a(1:lda,1:n) )
  allocate ( b(1:n) )
  allocate ( ipvt(1:n) )
  allocate ( x(1:n) )

  call matgen ( lda, n, a, x, b )
!
!  Factor the linear system.
!
  wtime = omp_get_wtime ( )
  call sgefa_r ( a, lda, n, ipvt, info )
  wtime = omp_get_wtime ( ) - wtime

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST_SGEFA_R - Fatal error!'
    write ( *, '(a)' ) '  SGEFA_R reports the matrix is singular.'
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call sgesl ( a, lda, n, ipvt, b, job )

  err = sum ( abs ( x(1:n) - b(1:n) ) )

  write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) &
    '  SGEFA_R   Sequential   ', n, err, wtime

  deallocate ( a )
  deallocate ( b )
  deallocate ( ipvt )
  deallocate ( x )

  return
end
subroutine test_sgefa_r_omp ( n )

!*****************************************************************************80
!
!! TEST_SGEFA_R_OMP tests a row-oriented version of SGEFA + OpenMP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  real ( kind = 4 ) err
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ), allocatable, dimension ( : ) :: ipvt
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
  real ( kind = 8 ) wtime
  real ( kind = 4 ), allocatable, dimension ( : ) :: x
!
!  Generate the linear system A * x = b.
!
  lda = n

  allocate ( a(1:lda,1:n) )
  allocate ( b(1:n) )
  allocate ( ipvt(1:n) )
  allocate ( x(1:n) )

  call matgen ( lda, n, a, x, b )
!
!  Factor the linear system.
!
  wtime = omp_get_wtime ( )
  call sgefa_r_omp ( a, lda, n, ipvt, info )
  wtime = omp_get_wtime ( ) - wtime

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST_SGEFA_R_OMP - Fatal error!'
    write ( *, '(a)' ) '  SGEFA_R_OMP reports the matrix is singular.'
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call sgesl ( a, lda, n, ipvt, b, job )

  err = sum ( abs ( x(1:n) - b(1:n) ) )

  write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) &
    '  SGEFA_R     Parallel   ', n, err, wtime

  deallocate ( a )
  deallocate ( b )
  deallocate ( ipvt )
  deallocate ( x )

  return
end
subroutine matgen ( lda, n, a, x, b )

!*****************************************************************************80
! 
!! MATGEN generates a "random" matrix for testing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix, and the length of the vector.
!
!    Output, real ( kind = 4 ) A(LDA,N), the matrix.
!
!    Output, real ( kind = 4 ) X(N), the solution vector.
!
!    Output, real ( kind = 4 ) B(N), the right hand side vector.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 4 ) value
  real ( kind = 4 ) x(n)

  seed = 1325
!
!  Set the matrix A.
!
  do j = 1, n
    do i = 1, n
      seed = mod ( 3125 * seed, 65536 )
      value = ( real ( seed, kind = 4 ) - 32768.0 ) / 16384.0
      a(i,j) = value
    end do
  end do
!
!  Set x.
!
  do i = 1, n
    x(i) = real ( i, kind = 4 ) / real ( n, kind = 4 )
  end do
!
!  Set b = A * x.
!
  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

  return
end
subroutine sgefa ( a, lda, n, ipvt, info )

!*********************************************************************72
!
!! SGEFA factors a real matrix by gaussian elimination.
!
!  Discussion:
!
!    This is a version of the LINPACK routine SGEFA which has been
!    simplified by replacing all the calls to BLAS routines by 
!    equivalent code.
!
!     on entry
!
!        a       real ( kind = 4 ) (lda, n)
!                the matrix to be factored.
!
!        lda     integer ( kind = 4 )
!                the leading dimension of the array  a .
!
!        n       integer ( kind = 4 )
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer ( kind = 4 ) (n)
!                an integer vector of pivot indices.
!
!        info    integer ( kind = 4 )
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that sgesl or sgedi will divide by zero
!                     if called.  use  rcond  in sgeco for a reliable
!                     indication of singularity.
!
  implicit none

  integer ( kind = 4 ) lda
  integer  ( kind = 4 )n

  real ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 4 ) t

  info = 0
  do k = 1, n - 1
!
!  Find the pivot index L.
!
    l = k
    t = abs ( a(k,k) )
    do i = k + 1, n
      if ( t < abs ( a(i,k) ) ) then
        l = i
        t = abs ( a(i,k) )
      end if
    end do
    ipvt(k) = l

    if ( a(l,k) == 0.0E+00 ) then
      info = k
      return
    end if
!
!  Interchange rows K and L.
!
    if ( l /= k ) then
      do j = k, n
        t      = a(l,j)
        a(l,j) = a(k,j)
        a(k,j) = t
      end do
    end if
!
!  Compute column K of the lower triangular factor.
!
    do i = k + 1, n
      a(i,k) = - a(i,k) / a(k,k)
    end do
!
!  Add multiples of the pivot row to the remaining rows.
!
    do j = k + 1, n
      do i = k + 1, n
        a(i,j) = a(i,j) + a(i,k) * a(k,j)
      end do
    end do

  end do

  ipvt(n) = n
  if ( a(n,n) == 0.0E+00 ) then
    info = n
  end if

  return
end
subroutine sgefa_c ( a, lda, n, ipvt, info )

!*****************************************************************************80
! 
!! SGEFA_C factors a matrix by gaussian elimination.
!
!  Discussion:
!
!    The step in which multiples of the pivot row are added to individual
!    rows has been replaced by a single call which updates the entire
!    matrix sub-block.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(LDA,N).  On input, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers which were 
!    used to obtain it.  The factorization can be written A = L * U where
!    L is a product of permutation and unit lower triangular matrices and
!    U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, indicates singularity.
!    If 0, this is the normal value, and the algorithm succeeded.
!    If K, then on the K-th elimination step, a zero pivot was encountered.
!    The matrix is numerically not invertible.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 4 ) t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1
!
!  Find the pivot index L.
!
    l = k
    t = abs ( a(k,k) )
    do i = k + 1, n
      if ( t < abs ( a(i,k) ) ) then
        l = i
        t = abs ( a(i,k) )
      end if
    end do
    ipvt(k) = l

    if ( a(l,k) == 0.0E+00 ) then
      info = k
      return
    end if
!
!  Interchange rows K and L.
!
    if ( l /= k ) then
      do j = k, n
        t      = a(l,j)
        a(l,j) = a(k,j)
        a(k,j) = t
      end do
    end if
!
!  Compute column K of the lower triangular factor.
!
    a(k+1:n,k) = - a(k+1:n,k) / a(k,k)
!
!  Add multiples of the pivot row to the remaining rows.
!
    do j = k+1, n
      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)
    end do

  end do

  ipvt(n) = n

  if ( a(n,n) == 0.0 ) then
    info = n
  end if

  return
end
subroutine sgefa_c_omp ( a, lda, n, ipvt, info )

!*****************************************************************************80
! 
!! SGEFA_C_OMP factors a matrix by gaussian elimination.
!
!  Discussion:
!
!    This variant of SGEFA uses OpenMP for improved parallel execution.
!
!    The step in which multiples of the pivot row are added to individual
!    rows has been replaced by a single call which updates the entire
!    matrix sub-block.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(LDA,N).  On input, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers which were 
!    used to obtain it.  The factorization can be written A = L * U where
!    L is a product of permutation and unit lower triangular matrices and
!    U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, indicates singularity.
!    If 0, this is the normal value, and the algorithm succeeded.
!    If K, then on the K-th elimination step, a zero pivot was encountered.
!    The matrix is numerically not invertible.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 4 ) t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1
!
!  Find the pivot index L.
!
    l = k
    t = abs ( a(k,k) )
    do i = k + 1, n
      if ( t < abs ( a(i,k) ) ) then
        l = i
        t = abs ( a(i,k) )
      end if
    end do
    ipvt(k) = l

    if ( a(l,k) == 0.0E+00 ) then
      info = k
      return
    end if
!
!  Interchange rows K and L.
!
    if ( l /= k ) then
!$omp parallel &
!$omp shared ( a, k, l, n ) &
!$omp private ( j, t )

!$omp do
      do j = k, n
        t      = a(l,j)
        a(l,j) = a(k,j)
        a(k,j) = t
      end do
!$omp end do

!$omp end parallel
    end if
!
!  Compute column K of the lower triangular factor.
!
!$omp parallel &
!$omp shared ( a, k, n ) &
!$omp private ( j )

!$omp do
    do j = k + 1, n
      a(j,k) = - a(j,k) / a(k,k)
    end do
!$omp end do

!$omp end parallel
!
!  Add multiples of the pivot row to the remaining rows.
!
!$omp parallel &
!$omp shared ( a, k, n ) &
!$omp private ( j )

!$omp do
    do j = k + 1, n
      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)
    end do
!$omp end do

!$omp end parallel

  end do

  ipvt(n) = n

  if ( a(n,n) == 0.0 ) then
    info = n
  end if

  return
end
subroutine sgefa_r ( a, lda, n, ipvt, info )

!*****************************************************************************80
! 
!! SGEFA_R factors a matrix by gaussian elimination.
!
!  Discussion:
!
!    The step in which multiples of the pivot row are added to individual
!    rows has been replaced by a single call which updates the entire
!    matrix sub-block.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(LDA,N).  On input, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers which were 
!    used to obtain it.  The factorization can be written A = L * U where
!    L is a product of permutation and unit lower triangular matrices and
!    U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, indicates singularity.
!    If 0, this is the normal value, and the algorithm succeeded.
!    If K, then on the K-th elimination step, a zero pivot was encountered.
!    The matrix is numerically not invertible.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 4 ) t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1
!
!  Find the pivot index L.
!
    l = k
    t = abs ( a(k,k) )
    do i = k + 1, n
      if ( t < abs ( a(i,k) ) ) then
        l = i
        t = abs ( a(i,k) )
      end if
    end do
    ipvt(k) = l

    if ( a(l,k) == 0.0E+00 ) then
      info = k
      return
    end if
!
!  Interchange rows K and L.
!
    if ( l /= k ) then
      do j = k, n
        t      = a(l,j)
        a(l,j) = a(k,j)
        a(k,j) = t
      end do
    end if
!
!  Compute column K of the lower triangular factor.
!
    a(k+1:n,k) = - a(k+1:n,k) / a(k,k)
!
!  Add multiples of the pivot row to the remaining rows.
!
   do i = k+1, n
     a(i,k+1:n) = a(i,k+1:n) + a(i,k) * a(k,k+1:n)
   end do

  end do

  ipvt(n) = n

  if ( a(n,n) == 0.0 ) then
    info = n
  end if

  return
end
subroutine sgefa_r_omp ( a, lda, n, ipvt, info )

!*****************************************************************************80
! 
!! SGEFA_R_OMP factors a matrix by gaussian elimination.
!
!  Discussion:
!
!    This variant of SGEFA uses OpenMP for improved parallel execution.
!
!    The step in which multiples of the pivot row are added to individual
!    rows has been replaced by a single call which updates the entire
!    matrix sub-block.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(LDA,N).  On input, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers which were 
!    used to obtain it.  The factorization can be written A = L * U where
!    L is a product of permutation and unit lower triangular matrices and
!    U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, indicates singularity.
!    If 0, this is the normal value, and the algorithm succeeded.
!    If K, then on the K-th elimination step, a zero pivot was encountered.
!    The matrix is numerically not invertible.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 4 ) t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1
!
!  Find the pivot index L.
!
    l = k
    t = abs ( a(k,k) )
    do i = k + 1, n
      if ( t < abs ( a(i,k) ) ) then
        l = i
        t = abs ( a(i,k) )
      end if
    end do
    ipvt(k) = l

    if ( a(l,k) == 0.0E+00 ) then
      info = k
      return
    end if
!
!  Interchange rows K and L.
!
    if ( l /= k ) then
      do j = k, n
        t      = a(l,j)
        a(l,j) = a(k,j)
        a(k,j) = t
      end do
    end if
!
!  Compute column K of the lower triangular factor.
!
!$omp parallel

!$omp workshare
    a(k+1:n,k) = - a(k+1:n,k) / a(k,k)
!$omp end workshare

!$omp end parallel
!
!  Add multiples of the pivot row to the remaining rows.
!
!$omp parallel &
!$omp shared ( a, k, n ) &
!$omp private ( i )

!$omp do
   do i = k+1, n
     a(i,k+1:n) = a(i,k+1:n) + a(i,k) * a(k,k+1:n)
   end do
!$omp end do

!$omp end parallel

  end do

  ipvt(n) = n

  if ( a(n,n) == 0.0 ) then
    info = n
  end if

  return
end
subroutine sgesl ( a, lda, n, ipvt, b, job )

!*****************************************************************************80
!
!! SGESL solves a real general linear system A * X = B.
!
!  Discussion:
!
!    SGESL can solve either of the systems A * X = B or A' * X = B.
!
!    The system matrix must have been factored by SGECO or SGEFA.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if SGECO has set 0.0 < RCOND
!    or SGEFA has set INFO == 0.
!
!  Modified:
!
!    13 May 2007
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(LDA,N), the output from SGECO or SGEFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from SGECO or SGEFA.
!
!    Input/output, real ( kind = 4 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB.
!    0, solve A * X = B;
!    nonzero, solve A' * X = B.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 4 ) t
!
!  Solve A * X = B.
!
  if ( job == 0 ) then

    do k = 1, n - 1

      l = ipvt(k)

      t = b(l)
      if ( l /= k ) then
        b(l) = b(k)
        b(k) = t
      end if

      b(k+1:n) = b(k+1:n) + t * a(k+1:n,k)

    end do

    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      t = - b(k)
      b(1:k-1) = b(1:k-1) + t * a(1:k-1,k)
    end do

  end if

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
