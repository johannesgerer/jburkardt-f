program main

!*****************************************************************************80
!
!! TOMS453_PRB tests BROMIN, ACM TOMS algorithm 453.
!
!  Modified:
!
!    12 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS453_PRB'
  write ( *, '(a)' ) '  FORTRAN80 version'
  write ( *, '(a)' ) '  Test the TOMS453 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS453_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BROMIN, ACM TOMS algorithm 453.
!
!  Modified:
!
!    11 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) nhalf_max
  integer ( kind = 4 ) n_num
  integer ( kind = 4 ) s_num

  parameter ( nhalf_max = 6 )
  parameter ( n_num = 3 )
  parameter ( s_num = 4 )

  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_half
  integer ( kind = 4 ) n_vec(n_num)
  real ( kind = 8 ) s
  real ( kind = 8 ) s_vec(s_num)
  real ( kind = 8 ) tol
  real ( kind = 8 ) total
  real ( kind = 8 ) wi(nhalf_max)
  real ( kind = 8 ) wr(nhalf_max)
  real ( kind = 8 ) xi(nhalf_max)
  real ( kind = 8 ) xr(nhalf_max)

  save n_vec
  save s_vec

  data n_vec / 6, 9, 12 /
  data s_vec / 0.0D+00, 0.1D+00, 1.0D+00, 4.0D+00 /

  tol = 0.1D-8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Determine abscissas and weights for'
  write ( *, '(a)' ) '  a variety of values of S and N.'

  do i = 1, n_num

    n = n_vec(i)
    n_half = ( n + 1 ) / 2

    do j = 1, s_num

      s = s_vec(j)

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  N = ', n
      write ( *, '(a,g14.6)' ) '  S = ', s

      call bromin ( n, s, tol, xr, xi, wr, wi, eps, ier )

      if ( 0 < ier ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) 'BROMIN returned IER = ', ier

      else

        if ( ier == -1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Note that the requested accuracy'
          write ( *, '(a)' ) '  was not achieved.'
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a,a)' ) '                           ', &
          'XR              XI              WR              WI'
        write ( *, '(a)' ) ' '
        total = 0.0D+00
        do kk = 1, n
          if ( kk .le. ( n - n_half ) ) then
            k = n_half + 1 - kk
            write ( *, '(2x,i8,2x,i8,4(2x,g14.6))' ) &
              kk, k, xr(k), - xi(k), wr(k), - wi(k)
            total = total + wr(k)
          else
            k = kk - ( n - n_half )
            write ( *, '(2x,i8,2x,i8,4(2x,g14.6))' ) &
              kk, k, xr(k),   xi(k), wr(k),   wi(k)
            total = total + wr(k)
          end if
        end do
        write ( *, '(2x,a8,2x,8x,16x,16x,2x,g14.6)' ) 'WR total', total

      end if

    end do

  end do

  return
end
