program main

!*****************************************************************************80
!
!! MAIN is the main program for GETWGT_PRB.
!
!  Discussion:
!
!    GETWGT_PRB tests the GETWGT routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: acid_num = 20

  integer ( kind = 4 ) aacnts(acid_num)
  integer ( kind = 4 ) acid_i
  integer ( kind = 4 ) event_i
  integer ( kind = 4 ) num
  real ( kind = 4 ) pseudocount(acid_num)

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GETWGT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the GETWGT library.'
!
!  Now observe a number of events, reporting each and updating
!  the pseudocount vector.
!
  do num = 1, 10

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) num, ' #8'
    write ( *, '(a)' ) ' '

    aacnts(1:acid_num) = 0
    aacnts(8) = num

    call getwgt ( aacnts, pseudocount )

    call r4vec_print ( acid_num, pseudocount, 'Pseudocounts:' )

  end do
!
!  Now observe a number of events, reporting each and updating
!  the pseudocount vector.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  50 #18'
  write ( *, '(a)' ) ' '

  aacnts(1:acid_num) = 0
  aacnts(18) = 50

  call getwgt ( aacnts, pseudocount )

  call r4vec_print ( acid_num, pseudocount, 'Pseudocounts:' )
!
!  Now observe a number of events, reporting each and updating
!  the pseudocount vector.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  50 #8'
  write ( *, '(a)' ) ' '

  aacnts(1:acid_num) = 0
  aacnts(8) = 50

  call getwgt ( aacnts, pseudocount )

  call r4vec_print ( acid_num, pseudocount, 'Pseudocounts:' )
!
!  Now observe a number of events, reporting each and updating
!  the pseudocount vector.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  50 #18'
  write ( *, '(a)' ) ' '

  aacnts(1:acid_num) = 0
  aacnts(18) = 50

  call getwgt ( aacnts, pseudocount )

  call r4vec_print ( acid_num, pseudocount, 'Pseudocounts:' )
!
!  Now observe a number of events, reporting each and updating
!  the pseudocount vector.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  500 #18'
  write ( *, '(a)' ) ' '

  aacnts(1:acid_num) = 0
  aacnts(18) = 500

  call getwgt ( aacnts, pseudocount )

  call r4vec_print ( acid_num, pseudocount, 'Pseudocounts:' )
!
!  Now observe a number of events, reporting each and updating
!  the pseudocount vector.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  49 #20'
  write ( *, '(a)' ) ' '

  aacnts(1:acid_num) = 0
  aacnts(20) = 49

  call getwgt ( aacnts, pseudocount )

  call r4vec_print ( acid_num, pseudocount, 'Pseudocounts:' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GETWGT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
