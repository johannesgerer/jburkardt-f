program main

!*****************************************************************************80
!
!! MAIN is the main program for GRAPH_PAPER_PRB.
!
!  Discussion:
!
!    GRAPH_PAPER_PRB is a sample calling program for GRAPH_PAPER routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 )  angle_num
  integer   ( kind = 4 )  circle_num
  character ( len = 255 ) file_name
  integer   ( kind = 4 )  n

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRAPH_PAPER_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the GRAPH_PAPER library.'
!
!  Make hexagonal graph paper.
!
  n = 10
  file_name = 'hexagonal_1.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling HEXAGONAL_1'
  write ( *, '(a,i8)' ) '    N = ', n
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call hexagonal_1 ( n, file_name )
!
!  Make hexagonal graph paper.
!
  n = 21
  file_name = 'hexagonal_2.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling HEXAGONAL_2'
  write ( *, '(a,i8)' ) '    N = ', n
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call hexagonal_2 ( n, file_name )
!
!  Make hexagonal graph paper.
!
  n = 21
  file_name = 'hexagonal_3.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling HEXAGONAL_3'
  write ( *, '(a,i8)' ) '    N = ', n
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call hexagonal_3 ( n, file_name )
!
!  Make hexagonal graph paper.
!
  n = 21
  file_name = 'hexagonal_4.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling HEXAGONAL_4'
  write ( *, '(a,i8)' ) '    N = ', n
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call hexagonal_4 ( n, file_name )
!
!  Make hexagonal graph paper for the game HEX.
!
  n = 14
  file_name = 'hexagonal_5.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling HEXAGONAL_5'
  write ( *, '(a,i8)' ) '    N = ', n
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call hexagonal_5 ( n, file_name )
!
!  Make polar graph paper.
!
  angle_num = 12
  circle_num = 10
  file_name = 'polar_1.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling POLAR_1'
  write ( *, '(a,i8)' ) '    ANGLE_NUM =  ', angle_num
  write ( *, '(a,i8)' ) '    CIRCLE_NUM = ', circle_num
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call polar_1 ( angle_num, circle_num, file_name )
!
!  Make staggered dot graph paper.
!
  n = 21
  file_name = 'staggered_2.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling STAGGERED_2'
  write ( *, '(a,i8)' ) '    N =  ', n
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call staggered_2 ( n, file_name )
!
!  Make a blank sudoku sheet.
!
  file_name = 'sudoku_blank.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling SUDOKU_SHEET_BLANK'
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call sudoku_sheet_blank ( file_name )
!
!  Make a filled sudokuo sheet.
!
  file_name = 'sudoku_filled.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling SUDOKU_SHEET_FILLED'
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call sudoku_sheet_filled ( file_name )
!
!  Make trangular graph paper.
!
  n = 21
  file_name = 'triangular_1.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling TRIANGULAR_1'
  write ( *, '(a,i8)' ) '    N =  ', n
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call triangular_1 ( n, file_name )
!
!  Make trangular graph paper.
! 
  n = 21
  file_name = 'triangular_2.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling TRIANGULAR_2'
  write ( *, '(a,i8)' ) '    N =  ', n
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call triangular_2 ( n, file_name )
!
!  Make uniform graph paper.
!
  file_name = 'uniform_1.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling UNIFORM_1'
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call uniform_1 ( file_name )
!
!  Make uniform graph paper.
!
  n = 20
  file_name = 'uniform_2.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling UNIFORM_2'
  write ( *, '(a,i8)' ) '    N =  ', n
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call uniform_2 ( n, file_name )
!
!  Make uniform graph paper.
!
  file_name = 'uniform_3.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling UNIFORM_3'
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call uniform_3 ( file_name )
!
!  Make uniform graph paper.
!
  n = 12
  file_name = 'uniform_4.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling UNIFORM_4'
  write ( *, '(a,i8)' ) '    N =  ', n
  write ( *, '(a)' ) '    Filename = "' // trim ( file_name ) // '".'

  call uniform_4 ( n, file_name )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRAPH_PAPER_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
 
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
