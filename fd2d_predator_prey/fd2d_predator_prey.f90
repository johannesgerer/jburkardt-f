program main

!*****************************************************************************80
!
!! MAIN is the main program for FD2D_PREDATOR_PREY.
!
!  Discussion:
!
!    The program simulates the dynamics of a 2D predator-prey system.
!
!  Modified:
!
!    29 November 2004
!
!  Author:
!
!    Marcus Garvie (MATLAB version)
!    John Burkardt (FORTRAN90 translation)
!
!  Reference:
!
!    M R Garvie,
!    Computational Algorithms for Spatially Extended Predator-Prey
!      Systems with a Holling Type II Functional Response,
!    To appear.
!
!  Parameters:
!
!    User input, real ( kind = 8 ) A, the left endpoints of the interval.
!
!    User input, real ( kind = 8 ) ALPHA, the value of the parameter ALPHA.
!
!    User input, real ( kind = 8 ) B, the right endpoints of the interval.
!
!    User input, real ( kind = 8 ) BETA, the value of the parameter BETA.
!
!    Local parameter, integer BIGJ, the number of intervals in either
!    spatial coordinate.
!
!    User input, real ( kind = 8 ) DELT, the time step to use.
!
!    User input, real ( kind = 8 ) DELTA, the value of the parameter DELTA.
!
!    Local parameter, integer DIMJ, the number of nodes in either
!
!    User input, real ( kind = 8 ) GAMMA, the value of the parameter GAMMA.
!
!    User input, real ( kind = 8 ) H, the "space step", the desired spacing
!    between nodes in [A,B].
!
!    Local parameter, integer N, the total number of spatial nodes.
!
!    Local parameter, integer SOLVE, 0 for GMRES routine, 1 for Jacobi.
!
!    User input, real ( kind = 8 ) T, the final time.  The problem is to be
!    integrated from 0 to T.
!
  implicit none

  integer, parameter :: maxl = 10

  integer, parameter :: ligw = 20

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  integer, allocatable, dimension ( : ) :: b_i
  integer, allocatable, dimension ( : ) :: b_j
  integer, allocatable, dimension ( : ) :: b_ppm
  real ( kind = 8 ), allocatable, dimension ( : ) :: b_u
  real ( kind = 8 ), allocatable, dimension ( : ) :: b_v
  real ( kind = 8 ) beta
  integer bigj
  real ( kind = 8 ) delt
  real ( kind = 8 ) delta
  real ( kind = 8 ) diff
  integer dimj
  real ( kind = 8 ) dummy(1)
  real ( kind = 8 ) err
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  real ( kind = 8 ), allocatable, dimension ( : ) :: g
  integer, allocatable, dimension ( : ) :: g_ppm
  real ( kind = 8 ) gamma
  real ( kind = 8 ) h
  real ( kind = 8 ), allocatable, dimension ( : ) :: hhat
  integer i
  integer i2
  integer idummy(1)
  integer ierr
  integer ierror
  integer igwk(ligw)
  integer isym
  integer it
  integer it_max
  integer itol
  integer iunit
  integer j
  integer job
  integer k
  integer lrgw
  external matvec_triad
  integer mode
  external msolve_identity
  real ( kind = 8 ) mu
  integer n
  integer nz
  integer nz_num
  integer, allocatable, dimension ( : ) :: r_ppm
  real ( kind = 8 ), allocatable, dimension ( : ) :: rgwk
  integer solve
  real ( kind = 8 ) t
  integer time_step
  integer time_steps
  real ( kind = 8 ) tol
  real ( kind = 8 ), allocatable, dimension ( : ) :: u
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: u_grid
  real ( kind = 8 ) u_max
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: u0
  real ( kind = 8 ), allocatable, dimension ( : ) :: v
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v_grid
  real ( kind = 8 ) v_max
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v0
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x_grid
  real ( kind = 8 ), allocatable, dimension ( : ) :: y
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: y_grid
  real ( kind = 8 ), allocatable, dimension ( : ) :: y1
  real ( kind = 8 ), allocatable, dimension ( : ) :: y2

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD2D_PREDATOR_PREY'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A two dimensional finite difference algorithm'
  write ( *, '(a)' ) '  for a predator-prey system.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter parameter alpha:'
  read ( *, * ) alpha
  write ( *, '(a)' ) '  Enter parameter beta:'
  read ( *, * ) beta
  write ( *, '(a)' ) '  Enter parameter gamma:'
  read ( *, * ) gamma
  write ( *, '(a)' ) '  Enter parameter delta:'
  read ( *, * ) delta

  write ( *, '(a)' ) '  Enter a in [a,b]:'
  read ( *, * ) a
  write ( *, '(a)' ) '  Enter b in [a,b]:'
  read ( *, * ) b

  write ( *, '(a)' ) '  Enter space-step h:'
  read ( *, * ) h
  write ( *, '(a)' ) '  Enter maximum time t:'
  read ( *, * ) t
  write ( *, '(a)' ) '  Enter time-step Delta t:'
  read ( *, * ) delt
  write ( *, '(a)' ) '  Enter SOLVE ( 0 = GMRES, 1 = Jacobi ):'
  read ( *, * ) solve

  if ( solve /= 0 ) then
    solve = 1
  end if

  time_steps = nint ( t / delt )

  mu = delt / ( h * h )
  bigj = ( b - a ) / h
  dimj = 1 + int ( ( b - a ) / h )
  n = dimj * dimj

  write ( *, '(a)'       ) ' '
  write ( *, '(a)'       ) '  Parameters:'
  write ( *, '(a)'       ) ' '
  write ( *, '(a,g14.6)' ) '    ALPHA =      ', alpha
  write ( *, '(a,g14.6)' ) '    BETA =       ', beta
  write ( *, '(a,g14.6)' ) '    GAMMA =      ', gamma
  write ( *, '(a,g14.6)' ) '    DELTA =      ', delta
  write ( *, '(a)'       ) ' '
  write ( *, '(a,g14.6)' ) '    A =          ', a
  write ( *, '(a,g14.6)' ) '    B =          ', b
  write ( *, '(a,g14.6)' ) '    H =          ', h
  write ( *, '(a,i6)'    ) '    BIGJ =       ', bigj
  write ( *, '(a)'       ) '    (number of intervals on one side)'
  write ( *, '(a,i6)'    ) '    DIMJ =       ', dimj
  write ( *, '(a)'       ) '    (number of nodes on one side)'
  write ( *, '(a,i12)'   ) '    N =          ', n
  write ( *, '(a)'       ) '    (total number of nodes)'
  write ( *, '(a)'       ) ' '
  write ( *, '(a,g14.6)' ) '    T =          ', t
  write ( *, '(a)'       ) '    (final time)'
  write ( *, '(a,g14.6)' ) '    DELT =       ', delt
  write ( *, '(a)'       ) '    (time step)'
  write ( *, '(a,i12)'   ) '    TIME_STEPS = ', time_steps
  write ( *, '(a)'       ) ' '
  write ( *, '(a,g14.6)' ) '    MU =         ', mu
  write ( *, '(a,g14.6)' ) '    SOLVE =      ', mu
  if ( solve == 0 ) then
    write ( *, '(a)'       ) '          = GMRES iteration.'
  else
    write ( *, '(a)'       ) '          = Jacobi iteration.'
  end if
  allocate ( f(1:n) )
  allocate ( g(1:n) )
  allocate ( hhat(1:n) )
  allocate ( u(1:n) )
  allocate ( u_grid(1:dimj,1:dimj) )
  allocate ( u0(1:dimj,1:dimj) )
  allocate ( v(1:n) )
  allocate ( v_grid(1:dimj,1:dimj) )
  allocate ( v0(1:dimj,1:dimj) )
  allocate ( x(1:dimj) )
  allocate ( x_grid(1:dimj,1:dimj) )
  allocate ( y(1:dimj) )
  allocate ( y_grid(1:dimj,1:dimj) )
  allocate ( y1(1:n) )
  allocate ( y2(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Arrays allocated.'

  do i = 1, dimj
    x(i) = ( i - 1 ) * h
  end do

  do i = 1, dimj
    y(i) = ( i - 1 ) * h
  end do

  do j = 1, dimj
    x_grid(1:dimj,j) = x(j)
  end do

  do i = 1, dimj
    y_grid(i,1:dimj) = y(i)
  end do

  call u_init ( alpha, beta, gamma, delta, dimj, dimj, x_grid, y_grid, u0 )
  call v_init ( alpha, beta, gamma, delta, dimj, dimj, x_grid, y_grid, v0 )

  k = 0
  do j = 1, dimj
    do i = 1, dimj
      k = k + 1
      u(k) = u0(i,j)
      v(k) = v0(i,j)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Arrays initialized.'
!
!  MODE = 0, initialize SPARSE.
!
  allocate ( b_i(1) )
  allocate ( b_j(1) )
  allocate ( b_v(1) )

  mode = 0
  call sparse ( 0, 0, 0.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
!
!  MODE = 1, count entries.
!  MODE = 2, store entries into B_V.
!
  do mode = 1, 2

    call sparse ( 1, 1, 3.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    call sparse ( 1, 2, -1.5D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    call sparse ( bigj+1, bigj+1, 6.0D+00, n, n, mode, nz_num, nz, &
      b_i, b_j, b_v )
    call sparse ( bigj+1, bigj, -3.0D+00, n, n, mode, nz_num, nz, &
      b_i, b_j, b_v )

    do i = 2, bigj
      j = i+1
      call sparse ( i, j, -1.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
    do i = 2, bigj
      j = i
      call sparse ( i, j, 4.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
    do i = 2, bigj
      j = i-1
      call sparse ( i, j, -1.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
!
!  Construct block T.
!
    i = 1
    j = bigj + 2
    call sparse ( i, j, -1.5D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    i = bigj + 1
    j = 2 * bigj + 2
    call sparse ( i, j, -3.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    do i = 2, bigj
      j = bigj + 1 + i
      call sparse ( i, j, -2.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
!
!  Construct block Z.
!
    i = n - bigj
    j = n - bigj
    call sparse ( i, j, 6.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    i = n - bigj
    j = n - bigj + 1
    call sparse ( i, j, -3.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    i = n
    j = n
    call sparse ( i, j, 3.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    i = n
    j = n - 1
    call sparse ( i, j, -1.5D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )

    do i = n-bigj+1, n-1
      j = i + 1
      call sparse ( i, j, -1.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
    do i = n-bigj+1, n-1
      j = i
      call sparse ( i, j,  4.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
    do i = n-bigj+1, n-1
      j = i - 1
      call sparse ( i, j, -1.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
!
!  Construct block Y.
!
    i = n - bigj
    j = n - ( 2 * bigj + 1 )
    call sparse ( i, j, -3.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    i = n
    j = n - ( bigj + 1 )
    call sparse ( i, j, -1.5D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    do i = n-bigj+1, n-1
      j = i - bigj - 1
      call sparse ( i, j, -2.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
!
!  Upper W blocks.
!
    do i = bigj+2, n-bigj-1
      j = i + bigj + 1
      call sparse ( i, j, -1.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
!
!  Lower W blocks.
!
    do i = bigj+2, n-bigj-1
      j = i - bigj - 1
      call sparse ( i, j, -1.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do

    do i = bigj+2, n-bigj-1
      j = i
      call sparse ( i, j,  4.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
!
!  Upper diagonals of X blocks.
!
    do i = bigj+2, n-bigj-2
      j = i + 1
      call sparse ( i, j, -1.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
    do i = bigj+2, n-2*bigj-1, bigj+1
      j = i + 1
      call sparse ( i, j, -1.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
    do i = 2*bigj+2, n-2*bigj-2, bigj+1
      j = i + 1
      call sparse ( i, j, +1.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
!
!  Lower diagonals of X blocks.
!
    do i = bigj+3, n-bigj-1
      j = i - 1
      call sparse ( i, j, -1.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
    do i = 2*bigj+2, n-bigj-1, bigj+1
      j = i - 1
      call sparse ( i, j, -1.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
    do i = 2*bigj+3, n-2*bigj-1, bigj+1
      j = i - 1
      call sparse ( i, j, +1.0D+00, n, n, mode, nz_num, nz, b_i, b_j, b_v )
    end do
!
!  Once the entries are counted, you can allocate space.
!
    if ( mode == 1 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD2D_PREDATOR_PREY:'
      write ( *, '(a,i12)' ) '  Maximum number of nonzeros = ', nz_num
      deallocate ( b_i )
      deallocate ( b_j )
      deallocate ( b_v )
      allocate ( b_i(1:nz_num) )
      allocate ( b_j(1:nz_num) )
      allocate ( b_u(1:nz_num) )
      allocate ( b_v(1:nz_num) )

    end if

  end do
!
!  Construct B_U = I + mu * L.
!
  do k = 1, nz_num

    b_u(k) = mu * b_v(k)

    i = b_i(k)
    j = b_j(k)
    if ( i == j ) then
      b_u(k) = b_u(k) + 1.0D+00
    end if

  end do
!
!  Construct B_V = I + delta * mu * L.
!
  do k = 1, nz_num

    b_v(k) = delta * mu * b_v(k)

    i = b_i(k)
    j = b_j(k)
    if ( i == j ) then
      b_v(k) = b_v(k) + 1.0D+00
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sparse matrix values assigned.'
!
!  Now sort the (I,J) data, and adjust the B_U and B_V data accordingly.
!
  write ( *, * ) ' '
  write ( *, '(a)' ) '  Begin sorting.'
  call timestamp ( )

  call i4vec2_sort_a_plus2 ( nz_num, b_i, b_j, b_u, b_v )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sorting done.'

  if ( .false. ) then
    do i = 1, nz_num
      write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) b_i(i), b_j(i), b_v(i)
    end do
  end if
!
!  Merge any multiple occurrences of (I,J) information.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compress the matrix data.'
  write ( *, '(a,i12)' ) '  Initial number of matrix entries, NZ_NUM = ', nz_num

  i2 = 1
  do i = 2, nz_num
    if ( b_i(i2) == b_i(i) .and. b_j(i2) == b_j(i) ) then
      b_u(i2) = b_u(i2) + b_u(i)
      b_v(i2) = b_v(i2) + b_v(i)
    else
      i2 = i2 + 1
      b_i(i2) = b_i(i)
      b_j(i2) = b_j(i)
      b_u(i2) = b_u(i)
      b_v(i2) = b_v(i)
    end if
  end do

  nz_num = i2

  write ( *, '(a,i12)' ) '  After compression,                NZ_NUM = ', nz_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Begin time loop.'
!
!  Time loop.
!
  if ( solve == 0 ) then

    isym = 0
    itol = 0
    tol = 0.00001D+00
    it_max = 10
    iunit = 0
    lrgw = 1 + n * ( maxl + 6 ) + maxl * ( maxl + 3 )
    allocate ( rgwk(1:lrgw) )
    igwk(1) = maxl
    igwk(2) = maxl
    igwk(3) = 0
    igwk(4) = 0
    igwk(5) = 10

  else

    tol = 0.00001D+00
    isym = 0
    it_max = 10
    job = 0

    call ds3_diagonal2 ( n, nz_num, isym, b_i, b_j, b_u, b_v )

  end if

  do time_step = 1, time_steps
!
!  About 10 times during the iteration, print out the current time step.
!
    if ( time_steps * ( ( 10 * time_step ) / time_steps ) &
      == 10 * time_step ) then
      write ( *, * ) '  TIME_STEP = ', time_step
      call timestamp ( )
    end if

    hhat(1:n) = u(1:n) / ( alpha + abs ( u(1:n) ) )
    f(1:n) = u(1:n) - u(1:n) * abs ( u(1:n) ) - v(1:n) * hhat(1:n)
    g(1:n) = beta * v(1:n) * hhat(1:n) - gamma * v(1:n)
    y1(1:n) = u(1:n) + delt * f(1:n)

    if ( .false. ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Right hand side Y1'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, * ) i, y1(i)
      end do
    end if

    y2(1:n) = v(1:n) + delt * g(1:n)
!
!  Solve BU * u_new = y1
!        BV * v_new = y2
!
    if ( solve == 0 ) then

      call dgmres ( n, y1, u, nz_num, b_i, b_j, b_u, isym, matvec_triad, &
        msolve_identity, itol, tol, it_max, it, err, ierr, iunit, dummy, &
        dummy, rgwk, lrgw, igwk, ligw, dummy, idummy )

      if ( time_steps * ( ( 10 * time_step ) / time_steps ) &
        == 10 * time_step ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  Number of iterations:  ', it
        write ( *, '(a,g14.6)' ) '  Convergence measure is ', rgwk(1)
        write ( *, '(a,g14.6)' ) '  Error estimate ', err
        write ( *, '(a,i6)' ) '  Error code is ', ierr
      end if

    else

      call ds3_jac_sl ( n, nz_num, isym, b_i, b_j, b_u, y1, u, tol, it_max, &
        job, it, diff )

      if ( time_steps * ( ( 10 * time_step ) / time_steps ) &
        == 10 * time_step ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  Number of iterations:  ', it
        write ( *, '(a,g14.6)' ) '  Last change in solution ', diff
      end if

    end if

    if ( .false. ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution U'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, * ) i, u(i)
      end do
    end if

    if ( solve == 0 ) then

      call dgmres ( n, y2, v, nz_num, b_i, b_j, b_v, isym, matvec_triad, &
        msolve_identity, itol, tol, it_max, it, err, ierr, iunit, dummy, &
        dummy, rgwk, lrgw, igwk, ligw, dummy, idummy )

      if ( .false. ) then
        write ( *, '(a)'       ) ' '
        write ( *, '(a,i6)'    ) '  Number of iterations:  ', it
        write ( *, '(a,g14.6)' ) '  Convergence measure is ', rgwk(1)
        write ( *, '(a,g14.6)' ) '  Error estimate ', err
        write ( *, '(a,i6)'    ) '  Error code is ', ierr
      end if

    else

      call ds3_jac_sl ( n, nz_num, isym, b_i, b_j, b_v, y2, v, tol, it_max, &
        job, it, diff )

      if ( time_steps * ( ( 10 * time_step ) / time_steps ) &
        == 10 * time_step ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  Number of iterations:  ', it
        write ( *, '(a,g14.6)' ) '  Last change in solution ', diff
      end if

    end if

    if ( .false. ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution V'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, * ) i, v(i)
      end do
    end if

  end do
!
!  Make a PPM plot of the solution components U and V at the final time.
!
  allocate ( b_ppm(1:n) )
  allocate ( g_ppm(1:n) )
  allocate ( r_ppm(1:n) )

  u_max = maxval ( u(1:n) )
  r_ppm(1:n) = nint ( 255 *           u(1:n)   / u_max )
  g_ppm(1:n) = 0
  b_ppm(1:n) = nint ( 255 * ( u_max - u(1:n) ) / u_max )

  call ppma_write ( 'u2d.ppma', dimj, dimj, r_ppm, g_ppm, b_ppm, ierror )

  v_max = maxval ( v(1:n) )
  r_ppm(1:n) = nint ( 255 *           v(1:n)   / v_max )
  g_ppm(1:n) = 0
  b_ppm(1:n) = nint ( 255 * ( v_max - v(1:n) ) / v_max )

  call ppma_write ( 'v2d.ppma', dimj, dimj, r_ppm, g_ppm, b_ppm, ierror )
!
!  Free up the allocated memory.
!
  deallocate ( b_i )
  deallocate ( b_j )
  deallocate ( b_ppm )
  deallocate ( b_u )
  deallocate ( b_v )
  deallocate ( f )
  deallocate ( g )
  deallocate ( g_ppm )
  deallocate ( hhat )
  deallocate ( r_ppm )
  if ( solve == 0 ) then
    deallocate ( rgwk )
  end if
  deallocate ( u )
  deallocate ( u_grid )
  deallocate ( u0 )
  deallocate ( v )
  deallocate ( v_grid )
  deallocate ( v0 )
  deallocate ( x )
  deallocate ( x_grid )
  deallocate ( y )
  deallocate ( y_grid )
  deallocate ( y1 )
  deallocate ( y2 )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD2D_PREDATOR_PREY:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ds3_diagonal2 ( n, nz_num, isym, row, col, a, b )

!*****************************************************************************80
!
!! DS3_DIAGONAL2 reorders two square DS3 matrices so diagonal entries are first.
!
!  Discussion:
!
!    The DS3 storage format corresponds to the SLAP Triad format.
!
!    The DS3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!    This routine reorders the entries of A so that the first N entries
!    are exactly the diagonal entries of the matrix, in order.
!
!  Modified:
!
!    24 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in the matrix.
!
!    Input, integer ISYM, is 0 if the matrix is not symmetric, and 1
!    if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column indices
!    of the nonzero elements.
!
!    Input/output, real ( kind = 8 ) A(NZ_NUM), the nonzero elements
!    of the matrix.
!
!    Input/output, real ( kind = 8 ) B(NZ_NUM), the nonzero elements
!    of the second matrix.
!
  implicit none

  integer n
  integer nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(nz_num)
  integer col(nz_num)
  integer found
  integer i
  integer isym
  integer k
  integer row(nz_num)

  found = 0

  do k = 1, nz_num

    do while ( row(k) == col(k) )

      if ( row(k) == k ) then
        found = found + 1
        exit
      end if

      i = row(k)

      call i4_swap ( row(i), row(k) )
      call i4_swap ( col(i), col(k) )
      call r8_swap ( a(i), a(k) )
      call r8_swap ( b(i), b(k) )

      found = found + 1

      if ( n <= found ) then
        exit
      end if

    end do

    if ( n <= found ) then
      exit
    end if

  end do

  if ( found < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DS3_DIAGONAL - Warning!'
    write ( *, '(a,i6)' ) '  Number of diagonal entries expected was ', n
    write ( *, '(a,i6)' ) '  Number found was ', found
    stop
  end if

  return
end
subroutine ds3_jac_sl ( n, nz_num, isym, row, col, a, b, x, tol, it_max, &
  job, it, diff )

!*****************************************************************************80
!
!! DS3_JAC_SL solves a DS3 system using Jacobi iteration.
!
!  Discussion:
!
!    The DS3 storage format corresponds to the SLAP Triad format.
!
!    The DS3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!    This routine REQUIRES that the matrix be square, that the matrix
!    have nonzero diagonal entries, and that the first N entries of
!    the array A be exactly the diagonal entries of the matrix, in order.
!
!  Modified:
!
!    27 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in the matrix.
!
!    Input, integer ISYM, is 0 if the matrix is not symmetric, and 1
!    if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column indices
!    of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = 8 ) X(N), an approximate solution
!    to the system.
!
!    Input, real ( kind = 8 ) TOL, a tolerance.  If the maximum change in
!    the solution is less than TOL, the iteration is terminated early.
!
!    Input, integer IT_MAX, the maximum number of iterations.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
!    Output, integer IT, the number of iterations taken.
!
!    Output, real ( kind = 8 ) DIFF, the maximum change in the solution
!    on the last iteration.
!
  implicit none

  integer n
  integer nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(n)
  integer col(nz_num)
  real ( kind = 8 ) diff
  integer i
  integer isym
  integer it
  integer it_max
  integer it_num
  integer j
  integer job
  integer k
  integer row(nz_num)
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xnew(n)
  real ( kind = 8 ) x_norm

  if ( job == 0 ) then

    do it_num = 1, it_max

      it = it_num
!
!  Initialize to right hand side.
!
      xnew(1:n) = b(1:n)
!
!  Subtract off-diagonal terms.
!
      do k = n+1, nz_num
        i = row(k)
        j = col(k)
        xnew(i) = xnew(i) - a(k) * x(j)
        if ( isym == 1 ) then
          xnew(j) = xnew(j) - a(k) * x(i)
        end if
      end do
!
!  Divide by diagonal terms.
!
      xnew(1:n) = xnew(1:n) / a(1:n)
!
!  Measure change:
!
      x_norm = maxval ( abs ( x(1:n) ) )
      diff = maxval ( abs ( xnew(1:n) - x(1:n) ) )
!
!  Update.
!
      x(1:n) = xnew(1:n)
!
!  Test for early termination.
!
      if ( diff <= tol  * ( x_norm + 1.0D+00 ) ) then
        exit
      end if

    end do

  else

    do it_num = 1, it_max

      it = it_num
!
!  Initialize to right hand side.
!
      xnew(1:n) = b(1:n)
!
!  Subtract off-diagonal terms.
!
      do k = n + 1, nz_num
        i = col(k)
        j = row(k)
        xnew(i) = xnew(i) - a(k) * x(j)
        if ( isym == 1 ) then
          xnew(j) = xnew(j) - a(k) * x(i)
        end if
      end do
!
!  Divide by diagonal terms.
!
      xnew(1:n) = xnew(1:n) / a(1:n)
!
!  Measure change:
!
      x_norm = maxval ( abs ( x(1:n) ) )
      diff = maxval ( abs ( xnew(1:n) - x(1:n) ) )
!
!  Update.
!
      x(1:n) = xnew(1:n)
!
!  Test for early termination.
!
      if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
        exit
      end if

    end do

  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical              lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
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
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer i
  integer j
  integer k

  k = i
  i = j
  j = k

  return
end
subroutine i4vec2_sort_a_plus2 ( n, a1, a2, a3, a4 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A_PLUS2 ascending sorts integer pairs, and adjusts real vectors.
!
!  Discussion:
!
!    The data items have the form (I,J,K,L), where I and J are integers
!    and K and L are real values.  The data is to be ascending sorted on the
!    I and J fields.
!
!  Modified:
!
!    05 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items of data.
!
!    Input/output, integer A1(N), integer A2(N), the data to be sorted.
!
!    Input/output, real ( kind = 8 ) A3(N), A4(N), the data to be rearranged
!    in accordance with the sorting of A1 and A2.
!
  implicit none

  integer n

  integer a1(n)
  integer a2(n)
  real ( kind = 8 ) a3(n)
  real ( kind = 8 ) a4(n)
  integer i
  integer indx
  integer isgn
  integer j
  integer temp
  real ( kind = 8 ) rtemp

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects and "carry" K.
!
    if ( 0 < indx ) then

      temp = a1(i)
      a1(i) = a1(j)
      a1(j) = temp

      temp = a2(i)
      a2(i) = a2(j)
      a2(j) = temp

      rtemp = a3(i)
      a3(i) = a3(j)
      a3(j) = rtemp

      rtemp = a4(i)
      a4(i) = a4(j)
      a4(j) = rtemp
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      if ( a1(i) < a1(j) ) then
        isgn = -1
      else if ( a1(i) == a1(j) .and. a2(i) <= a2(j) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine matvec_triad ( n, x, y, nelt, ia, ja, a, isym )

!*****************************************************************************80
!
!! MATVEC_TRIAD computes A*X for a matrix A stored in SLAP Triad form.
!
!  Modified:
!
!    21 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vectors.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) Y(N), the product A * X.
!
!    Input, integer NELT, the number of nonzero entries in A.
!
!    Input, integer IA(NELT), JA(NELT), real ( kind = 8 ) A(NELT), the data
!    structure storing the sparse matrix.
!
!    Input, integer ISYM, is 0 if all nonzero entries of the matrix
!    are stored, and 1 if only the diagonal and upper or lower triangle
!    are stored.
!
  implicit none

  integer n
  integer nelt

  real ( kind = 8 ) a(nelt)
  integer ia(nelt)
  integer isym
  integer ja(nelt)
  integer k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  y(1:n) = 0.0D+00

  do k = 1, nelt

    if ( ia(k) < 1 .or. n < ia(k) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MATVEC_TRIAD - Fatal error!'
      write ( *, '(a,i6,a,i6)' ) '  IA(',k,') = ', ia(k)
      stop
    end if

    if ( ja(k) < 1 .or. n < ja(k) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MATVEC_TRIAD - Fatal error!'
      write ( *, '(a,i6,a,i6)' ) '  JA(',k,') = ', ja(k)
      stop
    end if

    y(ia(k)) = y(ia(k)) + a(k) * x(ja(k))

  end do

  return
end
subroutine msolve_identity ( n, r, z, nelt, ia, ja, a, isym, rwork, iwork )

!*****************************************************************************80
!
!! MSOLVE_IDENTITY applies the identity matrix preconditioner.
!
!  Discussion:
!
!    Most SLAP solver routines require a preconditioner routine
!    that can solve M * Z = R.  If no preconditioning is required,
!    then you can simply behave as though the preconditioning matrix
!    M was the identity matrix.
!
!  Modified:
!
!    21 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vectors.
!
!    Input, real ( kind = 8 ) R(N), the right hand side.
!
!    Output, real ( kind = 8 ) Z(N), the solution of M * Z = R.
!
!    Input, integer NELT, the number of nonzero entries in A.
!
!    Input, integer IA(NELT), JA(NELT), real ( kind = 8 ) A(NELT), the data
!    structure storing the sparse matrix.
!
!    Input, integer ISYM, is 0 if all nonzero entries of the matrix
!    are stored, and 1 if only the diagonal and upper or lower triangle
!    are stored.
!
!    Input, real ( kind = 8 ) RWORK(*), integer IWORK(*), arrays that can be
!    used to pass information to the preconditioner.
!
  implicit none

  integer n
  integer nelt

  real ( kind = 8 ) a(nelt)
  integer ia(nelt)
  integer isym
  integer iwork(*)
  integer ja(nelt)
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rwork(*)
  real ( kind = 8 ) z(n)

  z(1:n) = r(1:n)

  return
end
subroutine ppma_write ( file_out_name, row_num, col_num, r, g, b, ierror )

!*****************************************************************************80
!
!! PPMA_WRITE writes an ASCII portable pixel map file.
!
!  Discussion:
!
!    PPM files can be viewed by XV.
!
!    Programs to convert files to this format include:
!
!      GIFTOPPM  - GIF file
!      PGMTOPPM  - Portable Gray Map file
!      PICTTOPPM - Macintosh PICT file
!      XPMTOPPM  - X11 pixmap file
!
!    Various programs can convert other formats to PPM format, including:
!
!      BMPTOPPM - Microsoft Windows BMP file.
!
!    A PPM file can also be converted to other formats, by programs:
!
!      PPMTOACAD - AutoCAD file
!      PPMTOGIF  - GIF file
!      PPMTOPGM  - Portable Gray Map file
!      PPMTOPICT - Macintosh PICT file
!      PPMTOPUZZ - X11 puzzle file
!      PPMTORGB3 - 3 Portable Gray Map files
!      PPMTOXPM  - X11 pixmap file
!      PPMTOYUV  - Abekas YUV file
!
!  Example:
!
!    P3
!    # feep.ppma created by PBMLIB(PPMA_WRITE).
!    4 4
!    15
!     0  0  0    0  0  0    0  0  0   15  0 15
!     0  0  0    0 15  7    0  0  0    0  0  0
!     0  0  0    0  0  0    0 15  7    0  0  0
!    15  0 15    0  0  0    0  0  0    0  0  0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the file
!    to which the data should be written.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows
!    and columns of data.
!
!    Input, integer ( kind = 4 ) R(ROW_NUM,COL_NUM), G(ROW_NUM,COL_NUM),
!    B(ROW_NUM,COL_NUM), the red, green and blue values of each pixel.  These
!    should be positive.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  integer   ( kind = 4 ) col_num
  integer   ( kind = 4 ) row_num

  integer   ( kind = 4 ) b(row_num,col_num)
  logical, parameter :: debug = .false.
  character ( len = * ) file_out_name
  integer   ( kind = 4 ) file_out_unit
  integer   ( kind = 4 ) g(row_num,col_num)
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ios
  integer   ( kind = 4 ) r(row_num,col_num)
  integer   ( kind = 4 ) rgb_max

  ierror = 0
!
!  Compute the maximum color value.
!
  rgb_max = max ( &
    maxval ( r(1:row_num,1:col_num) ), &
    maxval ( g(1:row_num,1:col_num) ), &
    maxval ( b(1:row_num,1:col_num) ) )
!
!  Open the file.
!
  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    ierror = 2
    return
  end if
!
!  Write the header.
!
  call ppma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
    rgb_max, ierror )
!
!  Write the data.
!
  call ppma_write_data ( file_out_unit, row_num, col_num, r, g, b, ierror )
!
!  Close the file.
!
  close ( unit = file_out_unit )
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE - Note:'
    write ( *, '(a)' ) '  The data was checked and written.'
    write ( *, '(a,i8)' ) '  Number of data rows ROW_NUM =    ', row_num
    write ( *, '(a,i8)' ) '  Number of data columns COL_NUM = ', col_num
    write ( *, '(a,i8)' ) '  Maximum RGB value RGB_MAX =      ', rgb_max
  end if

  return
end
subroutine ppma_write_data ( file_out_unit, row_num, col_num, r, g, b, ierror )

!*****************************************************************************80
!
!! PPMA_WRITE_DATA writes the data of a PPMA file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows
!    and columns of data.
!
!    Input, integer ( kind = 4 ) R(ROW_NUM,COL_NUM), G(ROW_NUM,COL_NUM),
!    B(ROW_NUM,COL_NUM), the red, green and blue values of each pixel.  These
!    should be positive.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  integer ( kind = 4 ) b(row_num,col_num)
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) r(row_num,col_num)

  ierror = 0
!
!  Write the header.
!
  do i = 1, row_num
    do jlo = 1, col_num, 4
      jhi = min ( jlo + 3, col_num )
      write ( file_out_unit, '(12i5)' ) ( r(i,j), g(i,j), b(i,j), j = jlo, jhi )
    end do
  end do

  return
end
subroutine ppma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
  rgb_max, ierror )

!*****************************************************************************80
!
!! PPMA_WRITE_HEADER writes the header of a PPMA file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and
!    columns of data.
!
!    Input, integer ( kind = 4 ) RGB_MAX, the maximum value of any
!    data component.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  character ( len = * )  file_out_name
  integer   ( kind = 4 ) file_out_unit
  integer   ( kind = 4 ) ierror
  character ( len = 2 ) :: magic = 'P3'
  integer   ( kind = 4 ) col_num
  integer   ( kind = 4 ) row_num
  integer   ( kind = 4 ) rgb_max

  ierror = 0
!
!  Write the header.
!
  write ( file_out_unit, '(a2)' ) magic
  write ( file_out_unit, '(a)' ) '# ' // trim ( file_out_name ) &
    // ' created by PPMA_IO::PPMA_WRITE.F90.'
  write ( file_out_unit, '(i5,2x,i5)' ) col_num, row_num
  write ( file_out_unit, '(i5)' ) rgb_max

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine sparse ( i, j, v, m, n, mode, nz_num, nz, b_i, b_j, b_v )

!*****************************************************************************80
!
!! SPARSE manages the storage of sparse matrix information.
!
!  Discussion:
!
!    To use this routine, first call it with MODE = 0 to initialize.
!
!    Then, for every nonzero entry of the matrix, call with MODE = 1,
!    just to count the entries.  The value of NZ_NUM at the end of
!    this step can be used to allocate storage for B_I, B_J, and B_V.
!
!    Then, for every nonzero entry, call again with MODE = 2, and the
!    I, J, and V values set.  These will be transferred into the
!    B_I, B_J and B_V arrays.
!
!    At the moment, the data still needs to be sorted and merged, but
!    that must be done externally to this routine.
!
!  Modified:
!
!    30 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the row and column index of the entry.
!
!    Input, real ( kind = 8 ) V, the value.
!
!    Input, integer M, N, the number of rows and columns in the matrix.
!
!    Input, integer MODE
!    0, initialize NZ and NZ_NUM.
!    1, "report" each entry, return a count in NZ_NUM.
!    2, store each entry, and keep a running count in NZ.
!    3, sort all the entries (not set up yet).
!    4, compress a sorted set of entries (not set up yet).
!
!    Input/output, integer NZ_NUM, is used in MODE = 0 and 1 only.
!    Mode 0 initializes NZ_NUM to 0.  Then, every call with MODE = 1
!    increments NZ_NUM by 1.  This results in a crude count of the maximum
!    number of storage entries needed.
!
!    Input/output, integer NZ, is used only in MODE = 0 and 2.
!    A call with MODE = 0 initializes NZ to 0.  Then, each call with
!    MODE = 2 presumably is providing another set of (I,J,V) data
!    to be stored in the sparse matrix array.  NZ is incremented
!    by 1, and the data is stored at that index.
!
!    Input/output, integer B_I(*), B_J(*), real ( kind = 8 ) B_V(*).  These
!    vectors are used to store the sparse matrix data.  Presumably, a dimension
!    of NZ_NUM (determined with MODE = 1) would be enough to store
!    the entries (assigned with MODE = 2).  Initially, the data is
!    stored as the user reports it, which means that the data is unsorted,
!    and that the same row and column indices may appear multiple times.
!
  implicit none

  integer b_i(*)
  integer b_j(*)
  real ( kind = 8 ) b_v(*)
  integer i
  integer j
  integer m
  integer mode
  integer n
  integer nz
  integer nz_num
  real ( kind = 8 ) v
!
!  Initialize.
!
  if ( mode == 0 ) then

    nz_num = 0
    nz = 0
!
!  Just count entries one by one.
!
  else if ( mode == 1 ) then

    nz_num = nz_num + 1
!
!  Store.
!
  else if ( mode == 2 ) then

    nz = nz + 1
    b_i(nz) = i
    b_j(nz) = j
    b_v(nz) = v

  else if ( mode == 3 ) then

  else if ( mode == 4 ) then

  end if

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Modified:
!
!    05 February 2004
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer i
  integer, save :: i_save = 0
  integer indx
  integer isgn
  integer j
  integer, save :: j_save = 0
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
