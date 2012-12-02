program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSE_GRID_MIXED_DATASET.
!
!  Discussion:
!
!    This program computes a sparse grid quadrature rule based on a mixture
!    of 1D rules, and writes it to a file.
!
!    The user specifies:
!    * M, the spatial dimension of the quadrature region,
!    * L, the level that defines the Smolyak grid.
!
!    Then the user specifies rules for each of the M dimensions.
!    A rule, when specified, may be used for one, or for multiple consecutive
!    dimensions.
!
!    * RULE identifies the 1D rule.
!      "CC", "F2", "GP", "GL", "GH", "GGH", "LG", "GLG", "GJ", "GW",
!      "CCS", "F2S", "GPS".
!    * the number of times the rule is to be used.
!    * ALPHA parameter for that rule;
!    * BETA parameter for that rule.
!
!  Licensing:
!
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    24 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( : ) :: alpha
  real ( kind = 8 )  alpha_1d
  integer   ( kind = 4 )  arg_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: beta
  real ( kind = 8 )  beta_1d
  integer   ( kind = 4 )  dim_inc
  integer   ( kind = 4 )  dim_index
  integer   ( kind = 4 )  dim_num
  character ( len = 255 ) file_name
  integer   ( kind = 4 )  iarg
  integer   ( kind = 4 )  iargc
  integer   ( kind = 4 )  ierror
  integer   ( kind = 4 )  last
  integer   ( kind = 4 )  level_max
  integer   ( kind = 4 ), allocatable, dimension ( : ) :: rule
  integer   ( kind = 4 )  rule_1d
  character ( len = 10 )  rule_string
  character ( len = 255 ) string
  real ( kind = 8 )  tol

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compute the abscissas and weights of a quadrature rule'
  write ( *, '(a)' ) '  associated with a sparse grid derived from a Smolyak'
  write ( *, '(a)' ) '  construction based on a mixture of 1D rules.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Inputs to the program include:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    DIM_NUM, the spatial dimension.'
  write ( *, '(a)' ) '    (typically in the range of 2 to 10)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    LEVEL_MAX, the "level" of the sparse grid.'
  write ( *, '(a)' ) '    (typically in the range of 0, 1, 2, 3, ...'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Then the user must define 1D quadrature rules to be used.'
  write ( *, '(a)' ) '   Each rule is used for at least the "next" dimension, but can be'
  write ( *, '(a)' ) '   used for several or all the remaining consecutive dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rule definition requires:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  * Rule identifier:'
  write ( *, '(a)' ) '    CC, F2, GP, GL, GH, GGH, LG, GLG, GJ, GW, CCS, F2S, GPS.'
  write ( *, '(a)' ) '  * Repetition factor (consecutive dimensions with same rule):'
  write ( *, '(a)' ) '  * ALPHA, (only for GGH, GLG, GJ rules)'
  write ( *, '(a)' ) '  * BETA, (only for GJ rule.)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output from the program includes:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    * Files that define the alphas, betas, ranges, weights, abscissas.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the spatial dimension.
!
  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, dim_num, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the value of DIM_NUM (1 or greater)'
    read ( *, * ) dim_num
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension requested is = ', dim_num
!
!  Get the level.
!
  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, level_max, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the value of LEVEL_MAX (0 or greater).'
    read ( *, * ) level_max
  end if
!
!  Now get the rules.
!
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )

  dim_index = 0

  do while ( dim_index < dim_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Rule identifiers include:'
    write ( *, '(a)' ) '  CC, F2, GP, GL, GH, GGH, LG, GLG, GJ, GW, CCS, F2S, GPS'
    write ( *, '(a,i8)' ) '  Enter the rule identifier for dimension ', dim_index + 1
    read ( *, * ) rule_string
    call rule_string_to_index ( rule_string, rule_1d )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  How many consecutive dimensions will this same rule be used?'
    read ( *, * ) dim_inc

    if ( dim_num < dim_index + dim_inc ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSE_GRID_MIXED_DATASET - Fatal error!'
      write ( *, '(a)' ) '  Dimension count exceeds limit.'
      stop
    end if

    if ( rule_1d == 6 .or. rule_1d == 8 .or. rule_1d == 9 ) then
      write ( *, * ) '  Enter the parameter ALPHA for this rule:'
      read ( *, * ) alpha_1d
    else
      alpha_1d = 0.0D+00
    end if

    if ( rule_1d == 9 ) then
      write ( *, * ) '  Enter the parameter BETA for this rule:'
      read ( *, * ) beta_1d
    else
      beta_1d = 0.0D+00
    end if

    rule(dim_index+1:dim_index+dim_inc) = rule_1d
    alpha(dim_index+1:dim_index+dim_inc) = alpha_1d
    beta(dim_index+1:dim_index+dim_inc) = beta_1d

    dim_index = dim_index + dim_inc

  end do
!
!  Get the filename.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter an identifier to use for the filenames:'
  read ( *, * ) file_name
!
!  Create the dataset.
!
  tol = sqrt ( epsilon ( tol ) )

  call sparse_grid_mixed_dataset_handle ( dim_num, level_max, rule, &
    alpha, beta, tol, file_name )
!
!  Free memory.
!
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_DATASET:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
!    which guarantee the ASCII collating sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  character ch
  integer   ( kind = 4 ) itemp

  itemp = iachar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
  end if

  return
end
subroutine rule_string_to_index ( rule_string, rule_1d )

!****************************************************************************80
!
!! RULE_STRING_TO_INDEX converts a string identifying a rule to an index.
!
!  Discussion:
!
!     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
!     2, "F2",  Fejer Type 2, Open Fully Nested rule.
!     3, "GP",  Gauss Patterson, Open Fully Nested rule.
!     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
!     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
!     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
!     7, "LG",  Gauss Laguerre, Open Non Nested rule.
!     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
!     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
!    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
!    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
!    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
!    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) RULE_STRING, a string.
!
!    Output, integer ( kind = 4 ) RULE_1D, the rule index.
!
  implicit none

  integer   ( kind = 4 ) rule_1d
  character ( len = * ) rule_string
  logical s_eqi

  if ( s_eqi ( rule_string, 'CC' ) ) then
    rule_1d = 1
  else if ( s_eqi ( rule_string, 'F2' ) ) then
    rule_1d = 2
  else if ( s_eqi ( rule_string, 'GP' ) ) then
    rule_1d = 3
  else if ( s_eqi ( rule_string, 'GL' ) ) then
    rule_1d = 4
  else if ( s_eqi ( rule_string, 'GH' ) ) then
    rule_1d = 5
  else if ( s_eqi ( rule_string, 'GGH' ) ) then
    rule_1d = 6
  else if ( s_eqi ( rule_string, 'LG' ) ) then
    rule_1d = 7
  else if ( s_eqi ( rule_string, 'GLG' ) ) then
    rule_1d = 8
  else if ( s_eqi ( rule_string, 'GJ' ) ) then
    rule_1d = 9
  else if ( s_eqi ( rule_string, 'GW' ) ) then
    rule_1d = 10
  else if ( s_eqi ( rule_string, 'CCS' ) ) then
    rule_1d = 11
  else if ( s_eqi ( rule_string, 'F2S' ) ) then
    rule_1d = 12
  else if ( s_eqi ( rule_string, 'GPS' ) ) then
    rule_1d = 13
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RULE_STRING_TO_INDEX - Fatal error!'
    write ( *, '(a)' ) '  Unexepected string.'
    stop
  end if

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Discussion:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character              c1
  character              c2
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lenc
  logical s_eqi
  character ( len = *  ) s1
  integer   ( kind = 4 ) s1_length
  character ( len = *  ) s2
  integer   ( kind = 4 ) s2_length

  s1_length = len ( s1 )
  s2_length = len ( s2 )
  lenc = min ( s1_length, s2_length )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, s1_length
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, s2_length
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

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

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) length
  character ( len = * )  s
  integer   ( kind = 4 ) state
  character :: TAB = achar ( 9 )
  integer   ( kind = 4 ) value

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
subroutine sparse_grid_mixed_dataset_handle ( dim_num, level_max, rule, &
  alpha, beta, tol, file_name )

!****************************************************************************80
!
!! SPARSE_GRID_MIXED_DATASET_HANDLE handles the creation of the dataset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the level that defines the grid.
!
!    Input, integer ( kind = 4 ) RULE(DIM_NUM), the rule in each dimension.
!     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
!     2, "F2",  Fejer Type 2, Open Fully Nested rule.
!     3, "GP",  Gauss Patterson, Open Fully Nested rule.
!     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
!     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
!     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
!     7, "LG",  Gauss Laguerre, Open Non Nested rule.
!     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
!     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
!    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
!    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
!    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
!    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.
!
!    Input, real ( kind = 8 ) ALPHA(DIM_NUM), BETA(DIM_NUM), parameters used for
!    Generalized Gauss Hermite, Generalized Gauss Laguerre, and Gauss Jacobi rules.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for point equality.
!
!    Input, character ( len = * ) FILE_NAME, the main name of the
!    output files.
!
  implicit none

  integer   ( kind = 4 ) dim_num

  real ( kind = 8 ) alpha(dim_num)
  real ( kind = 8 ) beta(dim_num)
  character ( len = *  ) file_name
  integer   ( kind = 4 ) level_max
  integer   ( kind = 4 ) point_num
  integer   ( kind = 4 ) point_total_num
  integer   ( kind = 4 ) rule(dim_num)
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_index
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_order
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sparse_point
  integer   ( kind = 4 ), allocatable, dimension(:) :: sparse_unique_index
  real ( kind = 8 ), allocatable, dimension ( : ) :: sparse_weight
  real ( kind = 8 ) tol
!
!  Compute necessary data.
!
  call sparse_grid_mixed_size_total ( dim_num, level_max, rule, &
    point_total_num )

  call sparse_grid_mixed_size ( dim_num, level_max, rule, alpha, beta, &
    tol, point_num )

  allocate ( sparse_unique_index(1:point_total_num) )

  call sparse_grid_mixed_unique_index ( dim_num, level_max, rule, alpha, &
    beta, tol, point_num, point_total_num, sparse_unique_index )

  allocate ( sparse_order(1:dim_num,1:point_num) )
  allocate ( sparse_index(1:dim_num,1:point_num) )

  call sparse_grid_mixed_index ( dim_num, level_max, rule, point_num, &
    point_total_num, sparse_unique_index, sparse_order, sparse_index )
!
!  Compute points and weights.
!
  allocate ( sparse_point(1:dim_num,1:point_num) )

  call sparse_grid_mixed_point ( dim_num, level_max, rule, alpha, beta, &
    point_num, sparse_order, sparse_index, sparse_point )

  allocate ( sparse_weight(1:point_num) )

  call sparse_grid_mixed_weight ( dim_num, level_max, rule, alpha, beta, &
    point_num, point_total_num, sparse_unique_index, sparse_weight )
!
!  Write points and weights to files.
!
  call sparse_grid_mixed_write ( dim_num, rule, alpha, beta, point_num, &
    sparse_weight, sparse_point, file_name )

  deallocate ( sparse_index )
  deallocate ( sparse_order )
  deallocate ( sparse_point )
  deallocate ( sparse_unique_index )
  deallocate ( sparse_weight )

  return
end
