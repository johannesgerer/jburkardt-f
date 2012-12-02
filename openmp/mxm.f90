program main

!****************************************************************************80
!
!! MAIN is the main program for MXM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer ( kind = 4 ) id
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) thread_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MXM'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix multiplication tests.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of processors available = ', omp_get_num_procs ( )
  write ( *, '(a,i8)' ) '  The number of threads available    = ', omp_get_max_threads ( )

  l = 500
  m = 500
  n = 500

  call r8_mxm ( l, m, n )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MXM:'
  write ( *, '(a)' ) '  Normal end of execution.'

  return
end
subroutine r8_mxm ( l, m, n )

!****************************************************************************80
!
!  Purpose:
!
!    R8_MXM carries out a matrix-matrix multiplication in R8 arithmetic.
!
!  Discussion:
!
!    A(LxN) = B(LxM) * C(MxN).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer L, M, N, the dimensions that specify the sizes of the
!    A, B, and C matrices.
!
  use omp_lib

  real    ( kind = 8 ) a(l,n)
  real    ( kind = 8 ) b(l,m)
  real    ( kind = 8 ) c(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ops
  real    ( kind = 8 ) r8_uniform_01
  real    ( kind = 8 ) rate
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) time_begin
  real    ( kind = 8 ) time_elapsed
  real    ( kind = 8 ) time_stop
!
!  Assign values to the B and C matrices.
!
  seed = 123456789

  do j = 1, m
    do i = 1, l
      b(i,j) = r8_uniform_01 ( seed )
    end do
  end do

  do j = 1, n
    do i = 1, m
      c(i,j) = r8_uniform_01 ( seed )
    end do
  end do
!
!  Compute A = B * C.
!
  time_begin = omp_get_wtime ( )

!$omp parallel &
!$omp shared ( a, b, c, l, m, n ) &
!$omp private ( i, j, k )

!$omp do
  do j = 1, n
    do i = 1, l
      a(i,j) = 0.0D+00
      do k = 1, m
        a(i,j) = a(i,j) + b(i,k) * c(k,j)
      end do
    end do
  end do
!$omp end do

!$omp end parallel

  time_stop = omp_get_wtime ( )
!
!  Report.
!
  ops = l * n * ( 2 * m )
  time_elapsed = time_stop - time_begin
  rate = dble ( ops ) / time_elapsed / 1000000.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_MXM matrix multiplication timing.'
  write ( *, '(a)' ) '  A(LxN) = B(LxM) * C(MxN).'
  write ( *, '(a,i8)' ) '  L = ', l
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a,i12)' ) '  Floating point OPS roughly ', ops
  write ( *, '(a,g14.6)' ) '  Elapsed time dT = ', time_elapsed
  write ( *, '(a,g14.6)' ) '  Rate = MegaOPS/dT = ', rate

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
  real    ( kind = 8 ) r8_uniform_01
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

