program main

!*****************************************************************************80
!
!! MAIN is the main program for GRID_DATASET.
!
!  Discussion:
!
!    GRID_DATASET generates a grid dataset and writes it to a file.
!
!    Interesting features of this problem are the determination
!    of the side of a grid that will generate "about" N points,
!    the method of dropping the extra points at random, and the
!    ability to center the grid inside the unit hypercube in a
!    number of ways.
!
!  Usage:
!
!    grid_dataset ( m, n, seed, center )
!
!    where
!
!    * M, the spatial dimension,
!    * N, the number of points to generate,
!    * SEED, the seed, a positive integer.
!    * CENTER, the grid centering option.
!      1: 0/(  N-1) ... (  N-1)/(  N-1)
!      2: 1/(  N+1) ...    N   /(  N+1)
!      3: 0/   N    ... (  N-1)/   N
!      4: 1/   N    ...    N   /   N
!      5: 1/(2*N)   ... (2*N-1)/(2*N  )
!
!    The program generates the data, writes it to the file
!
!      grid_M_N_CENTER.txt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) center
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) last
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = 255 ) output_filename
  real ( kind = 8 ), allocatable :: r(:,:)
  integer ( kind = 4 ) seed
  character ( len = 255 ) string

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRID_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Generate a grid dataset.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the spatial dimension M.
!
  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, m, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the spatial dimension M (1 or greater)'
    read ( *, * ) m
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
!
!  Get the number of points N.
!
  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, n, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the number of points N (1 or greater)'
    read ( *, * ) n
  end if

  write ( *, '(a,i8)' ) '  Number of points N = ', n
!
!  Get the seed, SEED
!
  if ( 3 <= arg_num ) then
    iarg = 3
    call getarg ( iarg, string )
    call s_to_i4 ( string, seed, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the seed SEED (1 or greater)'
    read ( *, * ) seed
  end if

  write ( *, '(a,i12)' ) '  SEED = ', seed

  if ( seed == 0 ) then
    call get_seed ( seed )
    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  Chosen value of SEED = ', seed
  end if
!
!  Get CENTER.
!
  if ( 4 <= arg_num ) then
    iarg = 4
    call getarg ( iarg, string )
    call s_to_i4 ( string, center, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter CENTER, the grid centering option.'
    write ( *, '(a)' ) '  Normal values are between 1 and 5:'
    write ( *, '(a)' ) '  1: 0/(  N-1) ... (  N-1)/(  N-1)'
    write ( *, '(a)' ) '  2: 1/(  N+1) ...    N   /(  N+1)'
    write ( *, '(a)' ) '  3: 0/   N    ... (  N-1)/   N'
    write ( *, '(a)' ) '  4: 1/   N    ...    N   /   N'
    write ( *, '(a)' ) '  5: 1/(2*N)   ... (2*N-1)/(2*N  )'
    read ( *, *, iostat = ios ) center
  end if

  write ( *, '(a,i8)' ) '  CENTER = ', center
!
!  Compute the data.
!
  allocate ( r(1:m,1:n) )

  call grid_generate ( m, n, center, seed, r )
!
!  Write the data to a file.
!
  write ( output_filename, '(a,i2.2,a,i5.5,a,i1,a)' ) &
    'grid_', m, '_', n, '_', center, '.txt'

  call r8mat_write ( output_filename, m, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The grid data was written to the file "' &
     // trim ( output_filename ) // '".'
!
!  Free memory.
!
  deallocate ( r )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRID_DATASET'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
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
  logical lopen

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
subroutine grid_generate ( dim_num, n, center, seed, r )

!*****************************************************************************80
!
!! GRID_GENERATE generates a grid dataset.
!
!  Discussion:
!
!    N points are needed in a DIM_NUM-dimensional space.
!
!    The points are to lie on a uniform grid of side N_SIDE.
!
!    Unless the N = N_SIDE**DIM_NUM for some N_SIDE, we can't use all the
!    points on a grid.  What we do is find the smallest N_SIDE
!    that's big enough, and randomly omit some points.
!
!    If N_SIDE is 4, then the choices in 1D are:
!
!    A: 0,   1/3, 2/3, 1
!    B: 1/5, 2/5, 3/5, 4/5
!    C: 0,   1/4, 2/4, 3/4
!    D: 1/4, 2/4, 3/4, 1
!    E: 1/8, 3/8, 5/8, 7/8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points to generate.
!
!    Input, integer ( kind = 4 ) CENTER, specifies the 1D grid centering:
!    1: first point is 0.0, last point is 1.0;
!    2: first point is 1/(N+1), last point is N/(N+1);
!    3: first point is 0, last point is (N-1)/N;
!    4: first point is 1/N, last point is 1;
!    5: first point is 1/(2*N), last point is (2*N-1)/(2*N);
!
!    Input/output, integer ( kind = 4 ) SEED, the random number seed.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) center
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n_grid
  integer ( kind = 4 ) n_side
  real ( kind = 8 ), dimension ( dim_num, n ) :: r
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_list(n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) tuple(dim_num)
!
!  Find the dimension of the smallest grid with N points.
!
  call grid_side ( dim_num, n, n_side )
!
!  We need to select N points out of N_SIDE**DIM_NUM set.
!
  n_grid = n_side**dim_num
!
!  Generate a random subset of N items from a set of size N_GRID.
!
  call ksub_random2 ( n_grid, n, seed, rank_list )
!
!  Must make one dummy call to TUPLE_NEXT_FAST with RANK = 0.
!
  rank = 0
  call tuple_next_fast ( n_side, dim_num, rank, tuple )
!
!  Now generate the appropriate indices, and "center" them.
!
  do j = 1, n

    rank = rank_list(j) - 1

    call tuple_next_fast ( n_side, dim_num, rank, tuple )

    if ( center == 1 ) then
      r(1:dim_num,j) = real (     tuple(1:dim_num) - 1, kind = 8 ) & 
                     / real (     n_side - 1,     kind = 8 )
    else if ( center == 2 ) then
      r(1:dim_num,j) = real (     tuple(1:dim_num),     kind = 8 ) &
                     / real (     n_side + 1,     kind = 8 )
    else if ( center == 3 ) then
      r(1:dim_num,j) = real (     tuple(1:dim_num) - 1, kind = 8 ) &
                     / real (     n_side,         kind = 8 )
    else if ( center == 4 ) then
      r(1:dim_num,j) = real (     tuple(1:dim_num),     kind = 8 ) &
                     / real (     n_side,         kind = 8 )
    else if ( center == 5 ) then
      r(1:dim_num,j) = real ( 2 * tuple(1:dim_num) - 1, kind = 8 ) &
                     / real ( 2 * n_side,         kind = 8 )
    end if

  end do

  return
end
subroutine grid_side ( dim_num, n, n_side )

!*****************************************************************************80
!
!! GRID_SIDE finds the smallest grid containing at least N points.
!
!  Discussion:
!
!    Each coordinate of the grid will have N_SIDE distinct values.
!    Thus the total number of points in the grid is N_SIDE**DIM_NUM.
!    This routine seeks the smallest N_SIDE such that N <= N_SIDE**DIM_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points to generate.
!
!    Output, integer ( kind = 4 ) N_SIDE, the length of one side of the smallest 
!    grid in DIM_NUM dimensions that contains at least N points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) exponent
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_side

  if ( n <= 0 ) then
    n_side = 0
    return
  end if

  if ( dim_num <= 0 ) then
    n_side = -1
    return
  end if

  exponent = 1.0D+00 / real ( dim_num, kind = 8 )

  n_side = int ( ( real ( n, kind = 8 ) )**exponent )

  if ( n_side**dim_num < n ) then
    n_side = n_side + 1
  end if

  return
end
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
!    use I4_HUGE() and HUGE interchangeably.
!
!    However, when using the G95, the values returned by HUGE were
!    not equal to 2147483647, apparently, and were causing severe
!    and obscure errors in my random number generator, which needs to
!    add I4_HUGE to the seed whenever the seed is negative.  So I
!    am backing away from invoking HUGE, whereas I4_HUGE is under
!    my control.
!
!    Explanation: because under G95 the default integer type is 64 bits!
!    So HUGE ( 1 ) = a very very huge integer indeed, whereas
!    I4_HUGE ( ) = the same old 32 bit big value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end
subroutine ksub_random2 ( n, k, seed, a )

!*****************************************************************************80
!
!! KSUB_RANDOM2 selects a random subset of size K from a set of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 April 2003
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
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
!    Input, integer ( kind = 4 ) N, the size of the set.
!
!    Input, integer ( kind = 4 ) K, the size of the subset, between 0 and N.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) A(K), the indices of the selected elements.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) available
  integer ( kind = 4 ) candidate
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) have
  integer ( kind = 4 ) n
  integer ( kind = 4 ) need
  real    ( kind = 8 ) r
  integer ( kind = 4 ) seed

  if ( k < 0 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM2 - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop
  end if

  if ( k == 0 ) then
    return
  end if

  need = k
  have = 0

  available = n
  candidate = 0

  do

    candidate = candidate + 1

    r = r8_uniform_01 ( seed )

    if ( real ( available, kind = 8 ) * r <= real ( need, kind = 8 ) ) then

      need = need - 1
      have = have + 1
      a(have) = candidate

      if ( need <= 0 ) then
        exit
      end if

    end if

    available = available - 1

  end do

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
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
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

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
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
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * )  output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

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
!    12 January 2009
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
!    Output, integer ( kind = 4 ) LENGTH, the number of characters 
!    of S used to make the integer.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * )  s
  integer ( kind = 4 ) state
  character :: TAB = achar ( 9 )
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

      if ( c == ' ' .or. c == TAB ) then

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

      if ( c == ' ' .or. c == TAB ) then

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
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
subroutine tuple_next_fast ( m, n, rank, x )

!*****************************************************************************80
!
!! TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
!
!  Discussion:
!
!    The elements are N vectors.  Each entry is constrained to lie
!    between 1 and M.  The elements are produced one at a time.
!    The first element is
!      (1,1,...,1)
!    and the last element is
!      (M,M,...,M)
!    Intermediate elements are produced in lexicographic order.
!
!    This code was written as a possibly faster version of TUPLE_NEXT.
!
!  Example:
!
!    N = 2,
!    M = 3
!
!    INPUT        OUTPUT
!    -------      -------
!    Rank  X      X
!    ----  ---    ---
!    0     * *    1 1
!    1     1 1    1 2
!    2     1 2    1 3
!    3     1 3    2 1
!    4     2 1    2 2
!    5     2 2    2 3
!    6     2 3    3 1
!    7     3 1    3 2
!    8     3 2    3 3
!    9     3 3    1 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 April 2003
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the maximum entry.
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, integer ( kind = 4 ) RANK, indicates the rank of the tuples.
!    On the very first call only, it is necessary that
!    the user set RANK = 0.
!
!    Input/output, integer ( kind = 4 ) X(N), on input the previous tuple.
!    On output, the next tuple.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), save, allocatable, dimension ( : ) :: base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  if ( rank == 0 ) then

    if ( allocated ( base ) ) then
      deallocate ( base )
    end if
    allocate ( base(1:n) )
    base(n) = 1
    do i = n-1, 1, -1
      base(i) = base(i+1) * m
    end do

  end if

  x(1:n) = mod ( rank / base(1:n), m ) + 1

  return
end
