program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_ZERO_PRB.
!
!  Discussion:
!
!    TEST_ZERO_PRB demonstrates the use of the TEST_ZERO scalar test functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_root = 4
  integer ( kind = 4 ), parameter :: max_start = 4

  real ( kind = 8 ), parameter :: fatol = 1.0D-06
  real ( kind = 8 ) fx
  real ( kind = 8 ) fxa
  real ( kind = 8 ) fxb
  real ( kind = 8 ) fxc
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: max_step = 25
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) r8_sign
  real ( kind = 8 ) range(2)
  integer ( kind = 4 ) root_num
  integer ( kind = 4 ) start_num
  character ( len = 80 ) title
  real ( kind = 8 ) x
  real ( kind = 8 ) xa
  real ( kind = 8 ), parameter :: xatol = 1.0D-06
  real ( kind = 8 ) xb
  real ( kind = 8 ) xc
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xroot(max_root)
  real ( kind = 8 ), parameter :: xrtol = 1.0D-06
  real ( kind = 8 ) xstart(max_start)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ZERO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_ZERO library.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Function value tolerance = ', fatol
  write ( *, '(a,g14.6)' ) '  Root absolute tolerance =  ', xatol
  write ( *, '(a,g14.6)' ) '  Root relative tolerance =  ', xrtol
  write ( *, '(a,i4)' ) '  Maximum number of steps =  ', max_step
!
!  Find out how many problems there are
!
  call p00_prob_num ( prob_num )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of problems available is ', prob_num

  do prob = 1, prob_num
!
!  Get the problem title.
!
    call p00_title ( prob, title )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Problem number ', prob
    write ( *, '(2x,a)' ) trim ( title )
!
!  Get the problem interval.
!
    call p00_range ( prob, range )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  We seek roots between'
    write ( *, '(2x,g14.6)' ) range(1)
    write ( *, '(a)' ) '  and'
    write ( *, '(2x,g14.6)' ) range(2)
!
!  Get the number of roots.
!
    call p00_root_num ( prob, root_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Number of known roots = ', root_num 
!
!  Get the roots.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I         X            F(X)'
    write ( *, '(a)' ) ' '
    do i = 1, root_num
      call p00_root ( prob, i, x )
      call p00_fx ( prob, x, fx )
      write ( *, '(2x,i4,2g16.8)' ) i, x, fx
    end do
!
!  Get the number of starting points.
!
    call p00_start_num ( prob, start_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Number of starting points = ', start_num 
!
!  Get the starting points.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  I  XSTART(I), F(XSTART(I))'
    write ( *, '(a)' ) ' '
    do i = 1, start_num
      call p00_start ( prob, i, xstart(i) )
      call p00_fx ( prob, xstart(i), fx )
      write ( *, '(2x,i2,2x,2g16.8)' ) i, xstart(i), fx
    end do
!
!  Bisection.
!
    call p00_start ( prob, 1, xa )
    call p00_fx ( prob, xa, fxa )

    do i = 2, start_num
      call p00_start ( prob, i, xb )
      call p00_fx ( prob, xb, fxb )
      if ( r8_sign ( fxa ) /= r8_sign ( fxb ) ) then
        call bisection ( fatol, max_step, prob, xatol, xa, xb, fxa, fxb )
        exit
      end if

    end do
!
!  Brent's method.
!
    call p00_start ( prob, 1, xa )
    call p00_fx ( prob, xa, fxa )

    do i = 2, start_num
      call p00_start ( prob, i, xb )
      call p00_fx ( prob, xb, fxb )
      if ( r8_sign ( fxa ) /= r8_sign ( fxb ) ) then
        call brent ( fatol, max_step, prob, xatol, xrtol, xa, xb, fxa, fxb )
        exit
      end if

    end do
!
!  Muller's method.
!
    if ( 3 <= start_num ) then

      call p00_start ( prob, 1, xa )
      call p00_fx ( prob, xa, fxa )
      call p00_start ( prob, 2, xb )
      call p00_fx ( prob, xb, fxb )
      call p00_start ( prob, 3, xc )
      call p00_fx ( prob, xc, fxc )

      call muller ( fatol, max_step, prob, xatol, xrtol, xa, xb, xc, &
        fxa, fxb, fxc )

    end if
!
!  Newton.
!
    do i = 1, start_num

      call p00_start ( prob, i, xa )
      call p00_fx ( prob, xa, fxa )

      call newton ( fatol, max_step, prob, xatol, xmin, xmax, xa, fxa )

    end do
!
!  Regula Falsi.
!
    call p00_start ( prob, 1, xa )
    call p00_fx ( prob, xa, fxa )

    do i = 2, start_num

      call p00_start ( prob, i, xb )
      call p00_fx ( prob, xb, fxb )
      if ( r8_sign ( fxa ) /= r8_sign ( fxb ) ) then
        call regula_falsi ( fatol, max_step, prob, xatol, xa, xb, fxa, fxb )
        exit
      end if

    end do
!
!  Secant.
!
    do i = 1, start_num - 1

      call p00_start ( prob, i, xa )
      call p00_fx ( prob, xa, fxa )

      call p00_start ( prob, i + 1, xb )
      call p00_fx ( prob, xb, fxb )

      call secant ( fatol, max_step, prob, xatol, xmin, xmax, xa, xb, fxa, fxb )

    end do

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ZERO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end

