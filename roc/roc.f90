program main

!*****************************************************************************80
!
!! MAIN is the main program for ROC.
!
!  Discussion:
!
!    ROC computes the Receiver Operator Characteristic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_score = 600

  logical debug
  character ( len = 256 ) debug_file_name
  integer ( kind = 4 ) debug_unit
  character ( len = 35 ) good(max_score)
  character ( len = 256 ) good_file_name
  character ( len = 256 ) graph_file_name
  integer ( kind = 4 ) graph_file_unit
  character ( len = 35 ) identifier(max_score)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(max_score)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) max_score_neg_select
  integer ( kind = 4 ) max_score_neg_select_hi
  integer ( kind = 4 ) max_score_neg_select_inc
  integer ( kind = 4 ) max_score_neg_select_lo
  integer ( kind = 4 ) max_score_select
  integer ( kind = 4 ) max_score_select_hi
  integer ( kind = 4 ) max_score_select_inc
  integer ( kind = 4 ) max_score_select_lo
  integer ( kind = 4 ) num_good
  integer ( kind = 4 ) num_score
  integer ( kind = 4 ) num_score_neg
  integer ( kind = 4 ) num_score_neg_select
  integer ( kind = 4 ) num_score_pos
  integer ( kind = 4 ) num_score_pos_select
  integer ( kind = 4 ) num_score_select
  integer ( kind = 4 ) num_score_skip
  logical plot
  integer ( kind = 4 ) result(max_score)
  real ( kind = 8 ) roc_val
  real ( kind = 8 ) score(max_score)
  real ( kind = 8 ) score_max
  real ( kind = 8 ) score_max_neg
  real ( kind = 8 ) score_max_pos
  real ( kind = 8 ) score_max_selected
  real ( kind = 8 ) score_min
  real ( kind = 8 ) score_min_neg
  real ( kind = 8 ) score_min_pos
  real ( kind = 8 ) score_min_selected
  real ( kind = 8 ) score_cutoff
  character ( len = 256 ) search_file_name
  character ( len = 10 ) search_type
  character ( len = 35 ) selection

  score_max_selected = 0.0D+00
  score_min_selected = 0.0D+00

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ROC'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Receiver-Operator-Characteristic (ROC) area calculator.'
  write ( *, '(a)' ) ' '
!
!  Set the ROC analysis parameters to default values.
!
  call params_default ( debug, debug_file_name, good_file_name, &
    graph_file_name, max_score_neg_select_hi, max_score_neg_select_inc, &
    max_score_neg_select_lo, max_score_select_hi, max_score_select_inc, &
    max_score_select_lo, plot, score_cutoff, search_file_name, search_type, &
    selection )
!
!  If DEBUG, open the debug file.
!
  if ( debug ) then

    call get_unit ( debug_unit )

    open ( unit = debug_unit, file = debug_file_name, status = 'replace', &
      iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ROC - Fatal error!'
      write ( *, '(a)' ) '  Could not open the debug file.'
      stop
    end if

  end if
!
!  Get the good file information.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter the good identifier file name:'
  read ( *, '(a)', iostat = ios ) good_file_name

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROC - Fatal error!'
    write ( *, '(a)' ) '  Unexpected end of input.'
    stop
  end if

  call good_get ( good, good_file_name, ierror, max_score, num_good )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROC - Fatal error!'
    write ( *, '(a)' ) '  Error reading the good identifiers.'
    stop
  end if

  if ( debug ) then
    call good_print ( debug_unit, good, good_file_name, max_score, num_good )
  end if
!
!  Get the data.
!
  call data_params_get ( ierror, search_file_name, search_type )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROC - Fatal error!'
    write ( *, '(a)' ) '  Error reading the data parameters.'
    stop
  end if

  if ( debug ) then
    call data_params_print ( debug_unit, search_file_name, search_type )
  end if

  call data_get ( debug, debug_unit, identifier, ierror, indx, max_score, &
    num_score, score, search_file_name, search_type )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROC - Fatal error!'
    write ( *, '(a)' ) '  Error reading the score data.'
    stop
  end if

  if ( debug ) then
    call data_print ( debug_unit, identifier, max_score, num_score, score )
  end if

  call data_count ( good, identifier, max_score, num_good, num_score, &
    num_score_neg, num_score_pos, score, score_max, score_max_neg, &
    score_max_pos, score_min, score_min_neg, score_min_pos )

  if ( debug ) then
    call data_count_print ( debug_unit, num_score, num_score_neg, &
      num_score_pos, score_max, score_max_neg, score_max_pos, score_min, &
      score_min_neg, score_min_pos )
  end if
!
!  Get the ROC analysis parameters.
!
  call analysis_params_get ( ierror, max_score_neg_select_hi, &
    max_score_neg_select_inc, max_score_neg_select_lo, max_score_select_hi, &
    max_score_select_inc, max_score_select_lo, score_cutoff, score_max_neg, &
    score_min, selection )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROC - Fatal error!'
    write ( *, '(a)' ) '  Error reading the parameters.'
    stop
  end if

  if ( debug ) then
    call analysis_params_print ( debug_unit, max_score_neg_select_hi, &
      max_score_neg_select_inc, max_score_neg_select_lo, max_score_select_hi, &
      max_score_select_inc, max_score_select_lo, score_cutoff, selection )
  end if
!
!  Analyze the data.
!
  call analysis_pre ( graph_file_name, graph_file_unit, ierror, plot )
!
!  Prepare for DO loop.
!
  if ( selection == 'NEGATIVE_RESTRICT' ) then
    max_score_neg_select = max_score_neg_select_lo
  else
    max_score_neg_select = 0
  end if

  if ( selection == 'TOTAL_RESTRICT' ) then
    max_score_select = max_score_select_lo
  else
    max_score_select = 0
  end if
!
!  One iteration of DO loop.
!
  do

    call analysis ( debug, debug_unit, good, good_file_name, graph_file_unit, &
      identifier, max_score, max_score_neg_select, max_score_select, num_good, &
      num_score, num_score_neg_select, num_score_pos, num_score_pos_select, &
      num_score_select, num_score_skip, plot, result, roc_val, score, &
      score_max_selected, score_min_selected, score_cutoff, search_file_name, &
      selection )
!
!  Print the analysis data.
!
    if ( debug ) then

      call analysis_print ( debug_unit, identifier, max_score, num_score, &
        num_score_neg, num_score_neg_select, num_score_pos, &
        num_score_pos_select, num_score_select, result, roc_val, score, &
        score_max_selected, score_min_selected )

    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  The ROC area integral is ', roc_val
    write ( *, '(a,g14.6,a,g14.6)' ) '  The range of the selected data was ', &
      score_min_selected, ' to ' , score_max_selected
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a,i8)' ) '  Selected ', num_score_pos_select, &
      ' positive items out of ', num_score_pos
    write ( *, '(a,i8,a,i8)' ) '  Selected ', num_score_neg_select, &
      ' negative items out of ', num_score_neg
    write ( *, '(a,i8,a,i8)' ) '  Selected ', num_score_select, &
      ' items out of          ', num_score

    if ( selection == 'NEGATIVE_RESTRICT' ) then
      if ( max_score_neg_select < max_score_neg_select_hi ) then
        max_score_neg_select = max_score_neg_select + max_score_neg_select_inc
        cycle
      end if
    end if

    if ( selection == 'TOTAL_RESTRICT' ) then
      if ( max_score_select < max_score_select_hi ) then
        max_score_select = max_score_select + max_score_select_inc
        cycle
      end if
    end if

    exit

  end do

  call analysis_post ( graph_file_unit, plot )
!
!  If DEBUG, close the debug file.
!
  if ( debug ) then
    close ( unit = debug_unit )
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ROC'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine analysis ( debug, debug_unit, good, good_file_name, &
  graph_file_unit, identifier, max_score, max_score_neg_select, &
  max_score_select, num_good, num_score, num_score_neg_select, num_score_pos, &
  num_score_pos_select, num_score_select, num_score_skip, plot, result, &
  roc_val, score, score_max_selected, score_min_selected, score_cutoff, &
  search_file_name, selection )

!*****************************************************************************80
!
!! ANALYSIS carries out the ROC analysis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical DEBUG, is TRUE if debug output is desired.
!
!    Input, integer ( kind = 4 ) DEBUG_UNIT, the unit number of the debug file.
!
!    Input, character ( len = 35 ) GOOD(MAX_SCORE), a list of identifiers
!    of all the objects which should receive a positive score.
!
!    Input, character ( len = * ) GOOD_FILE_NAME, the name of the good ID
!    file to be read.
!
!    Input, integer ( kind = 4 ) GRAPH_FILE_UNIT, the unit number associated with the
!    graph file.
!
!    Input, character ( len = 35 ) IDENTIFIER(MAX_SCORE), a list of the
!    (unique) identifiers for which a score was computed.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Input, integer ( kind = 4 ) MAX_SCORE_NEG_SELECT, the maximum number of scores with
!    negative results that will be considered.
!
!    Input, integer ( kind = 4 ) MAX_SCORE_SELECT, the maximum number of scores that
!    will be considered.
!
!    Input, integer ( kind = 4 ) NUM_GOOD, the number of "good" identifiers.
!
!    Input, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_NEG_SELECT, the number of selected scores
!    associated with negative results.
!
!    Input, integer ( kind = 4 ) NUM_SCORE_POS, the number of scores associated
!    with positive results.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_POS_SELECT, the number of selected scores
!    associated with positive results.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_SKIP, the number of unselected scores.
!
!    Input, logical PLOT, is TRUE if plots are desired.
!
!    Output, integer ( kind = 4 ) RESULT(MAX_SCORE), the result code.  RESULT(I) is:
!      J, item I is selected as the J-th POSITIVE representative;
!      0, item I is omitted;
!    - J, item I is selected as the J-th NEGATIVE representative.
!
!    Input, real ( kind = 8 ) ROC_VAL, the value of the ROC integral.
!
!    Input, real ( kind = 8 ) SCORE(MAX_SCORE), the scores.
!
!    Input, real ( kind = 8 ) SCORE_MAX_SELECTED, the maximum score of the selected data.
!
!    Input, real ( kind = 8 ) SCORE_MIN_SELECTED, the minimum score of the selected data.
!
!    Input, real ( kind = 8 ) SCORE_CUTOFF, the minimum score to be considered.
!
!    Input, character ( len = * ) SEARCH_FILE_NAME, the name of the search file.
!
!    Input, character ( len = * ) SELECTION, the selection criterion:
!    'ALL',               select all data;
!    'SCORE_MIN'          select all data with scores no less than
!                         SCORE_CUTOFF;
!    'POSITIVE_ALL'       select all data with scores no less than the lowest
!                         score of a positive data item;
!    'NEGATIVE_RESTRICT'  select all data with scores no less than that of
!                         the MAX_SCORE_NEG_SELECT-th negative item;
!    'TOTAL_RESTRICT'     select the MAX_SCORE_SELECT highest ranking items.
!
  implicit none

  integer ( kind = 4 ) max_score
  integer ( kind = 4 ), parameter :: num_plot = 101

  logical debug
  integer ( kind = 4 ) debug_unit
  character ( len = 35 ) good(max_score)
  character ( len = * ) good_file_name
  integer ( kind = 4 ) graph_file_unit
  integer ( kind = 4 ) i
  character ( len = 35 ) identifier(max_score)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) max_score_neg_select
  integer ( kind = 4 ) max_score_select
  real ( kind = 8 ) neg_plot(num_plot)
  integer ( kind = 4 ) num_good
  integer ( kind = 4 ) num_score
  integer ( kind = 4 ) num_score_neg_select
  integer ( kind = 4 ) num_score_pos
  integer ( kind = 4 ) num_score_pos_select
  integer ( kind = 4 ) num_score_select
  integer ( kind = 4 ) num_score_skip
  logical plot
  real ( kind = 8 ) pn_max_diff
  real ( kind = 8 ) pos_plot(num_plot)
  integer ( kind = 4 ) result(max_score)
  real ( kind = 8 ) roc_val
  real ( kind = 8 ) s_max_diff
  real ( kind = 8 ) score(max_score)
  real ( kind = 8 ) score_max_selected
  real ( kind = 8 ) score_min_selected
  real ( kind = 8 ) score_cutoff
  real ( kind = 8 ) score_plot(num_plot)
  character ( len = * ) search_file_name
  character ( len = 35 ) selection
!
!  Select the data to be analyzed.
!
  call data_select ( good, identifier, max_score, &
    max_score_neg_select, max_score_select, num_good, num_score, &
    num_score_neg_select, num_score_pos, &
    num_score_pos_select, num_score_select, num_score_skip, &
    result, score, score_cutoff, selection )
!
!  Compute the score range of the selected data.
!
  score_max_selected = - huge ( 1.0D+00 )
  score_min_selected = huge ( 1.0D+00 )
  do i = 1, num_score
    if ( result(i) /= 0 ) then
      score_max_selected = max ( score_max_selected, score(i) )
      score_min_selected = min ( score_min_selected, score(i) )
    end if
  end do
!
!  Compute the ROC integral.
!
  call pn_roc_int ( max_score, num_score, num_score_neg_select, &
    num_score_pos_select, result, roc_val )

  if ( .not. plot ) then
    return
  end if
!
!  Make the PN plot.
!  Here, we want the data to be DESCENDING with respect to SCORE.
!
  call pn_graph_file_write ( graph_file_unit, good_file_name, &
    ierror, max_score, num_good, num_score, num_score_neg_select, &
    num_score_pos_select, result, roc_val, search_file_name )
!
!  Make the PNS plot.
!  Here, we want the data to be ASCENDING with respect to SCORE.
!  Then adjust the SCORE limits.
!
  call r8vec_reverse ( num_score, score )

  call i4vec_reverse ( num_score, result )

  call svec_reverse ( num_score, identifier )

  if ( score_min_selected == score_max_selected ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ANALYSIS - Fatal error!'
    write ( *, '(a)' ) '  Minimum score = maximum score.'
    write ( *, '(a)' ) '  No plot is possible.'
    return
  end if

  call pns_data_set ( debug, debug_unit, max_score, neg_plot, &
    num_plot, num_score, pn_max_diff, pos_plot, result, &
    s_max_diff, score, score_max_selected, score_min_selected, score_plot )

  call pns_graph_file_write ( good_file_name, graph_file_unit, ierror, &
    neg_plot, num_good, num_plot, num_score_neg_select, &
    num_score_pos_select, pn_max_diff, pos_plot, s_max_diff, &
    search_file_name, score_plot )
!
!  In case we're in a loop, restore the DESCENDING order.
!
  call r8vec_reverse ( num_score, score )

  call i4vec_reverse ( num_score, result )

  call svec_reverse ( num_score, identifier )

  return
end
subroutine analysis_params_get ( ierror, max_score_neg_select_hi, &
  max_score_neg_select_inc, max_score_neg_select_lo, max_score_select_hi, &
  max_score_select_inc, max_score_select_lo, score_cutoff, score_max_neg, &
  score_min, selection )

!*****************************************************************************80
!
!! ANALYSIS_PARAMS_GET gets the ROC analysis parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Output, integer ( kind = 4 ) MAX_SCORE_NEG_SELECT_HI, MAX_SCORE_NEG_SELECT_INC,
!    MAX_SCORE_NEG_SELECT_LO, the maximum, increment, and minimum for the
!    DO loop that determines the number of scores with negative results
!    that will be considered.
!
!    Output, integer ( kind = 4 ) MAX_SCORE_SELECT_HI, MAX_SCORE_SELECT_INC,
!    MAX_SCORE_SELECT_LO, the maximum, increment, and minimum for the DO
!    loop that determines the number of scores that will be considered.
!
!    Output, real ( kind = 8 ) SCORE_CUTOFF, the minimum score to be considered.
!
!    Output, character ( len = * ) SELECTION, the selection criterion:
!    'ALL',               select all data;
!    'SCORE_MIN'          select all data with scores no less than
!                         SCORE_CUTOFF;
!    'POSITIVE_ALL'       select all data with scores no less than the lowest
!                         score of a positive data item;
!    'NEGATIVE_RESTRICT'  select all data with scores no less than that of
!                         the MAX_SCORE_NEG_SELECT-th negative item;
!    'TOTAL_RESTRICT'     select the MAX_SCORE_SELECT highest ranking items.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) max_score_neg_select_hi
  integer ( kind = 4 ) max_score_neg_select_inc
  integer ( kind = 4 ) max_score_neg_select_lo
  integer ( kind = 4 ) max_score_select_hi
  integer ( kind = 4 ) max_score_select_inc
  integer ( kind = 4 ) max_score_select_lo
  integer ( kind = 4 ) nchar
  real ( kind = 8 ) score_cutoff
  real ( kind = 8 ) score_max_neg
  real ( kind = 8 ) score_min
  character ( len = 35 ) selection
  character ( len = 100 ) string

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter the number or initial of the data selection:'
  write ( *, '(a)' ) '  1 ALL                 all data;'
  write ( *, '(a)' ) '  2 SCORE_MIN           all above min score;'
  write ( *, '(a)' ) '  3 POSITIVE_ALL        down to lowest positive;'
  write ( *, '(a)' ) '  4 NEGATIVE_RESTRICT   down to N-th negative;'
  write ( *, '(a)' ) '  5 TOTAL_RESTRICT      N data items;'
  write ( *, '(a)' ) '  0 QUIT'

  read ( *, '(a)', iostat = ios ) selection

  if ( ios /= 0 ) then
    ierror = 1
    return
  end if

  call s_cap ( selection )

  if ( selection == '0' ) then
    selection = 'QUIT'
  else if ( selection == '1' ) then
    selection = 'ALL'
  else if ( selection == '2' ) then
    selection = 'SCORE_MIN'
  else if ( selection == '3' ) then
    selection = 'POSITIVE_ALL'
  else if ( selection == '4' ) then
    selection = 'NEGATIVE_RESTRICT'
  else if ( selection == '5' ) then
    selection = 'TOTAL_RESTRICT'
  end if

  if ( selection(1:1) == 'Q' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'User-requested stop.'
    stop
  end if

  if ( selection(1:1) == 'A' ) then

  else if ( selection(1:1) == 'S' ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) 'The maximum score of a negative item is ', &
      score_max_neg
    write ( *, '(a)' ) 'Your cutoff score should be SMALLER than this.'
    write ( *, '(a,g14.6)' ) 'The minimum score of all items is ', score_min
    write ( *, '(a)' ) 'Your cutoff score should be BIGGER than this.'
    write ( *, '(a)' ) 'Enter the cutoff score:'

    read ( *, *, iostat = ios ) score_cutoff

    if ( ios /= 0 ) then
      ierror = 1
      return
    end if

  else if ( selection(1:1) == 'P' ) then

  else if ( selection(1:1) == 'N' ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'You may control the number of negative items.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter:'
    write ( *, '(a)' ) '      RETURN for no control;'
    write ( *, '(a)' ) '  or  one value, for a single limit;'
    write ( *, '(a)' ) '  or  minimum, maximum, increment, for a DO loop.'

    read ( *, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      ierror = 1
      return
    end if

    if ( string == ' ' ) then

      max_score_neg_select_lo = huge ( 1 )
      max_score_neg_select_hi = huge ( 1 )
      max_score_neg_select_inc = 1

    else

      call s_to_i4 ( string, intval, ierror, nchar )

      if ( ierror /= 0 ) then
        ierror = 0
        return
      end if

      max_score_neg_select_lo = intval
      max_score_neg_select_hi = intval
      max_score_neg_select_inc = 1

      string(1:nchar) = ' '

      call s_to_i4 ( string, intval, ierror, nchar )

      if ( ierror /= 0 ) then
        ierror = 0
        return
      end if

      max_score_neg_select_hi = intval

      string(1:nchar) = ' '
      call s_to_i4 ( string, intval, ierror, nchar )

      if ( ierror /= 0 ) then
        ierror = 0
        return
      end if

      max_score_neg_select_inc = intval

    end if

  else if ( selection(1:1) == 'T' ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'You may control the total number of items.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter:'
    write ( *, '(a)' ) '      RETURN for no control;'
    write ( *, '(a)' ) '  or  one value, for a single limit;'
    write ( *, '(a)' ) '  or  minimum, maximum, increment, for a DO loop.'

    read ( *, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      ierror = 1
      return
    end if

    if ( string == ' ' ) then

      max_score_select_lo = huge ( 1 )
      max_score_select_hi = huge ( 1 )
      max_score_select_inc = 1

    else

      call s_to_i4 ( string, intval, ierror, nchar )

      if ( ierror /= 0 ) then
        ierror = 0
        return
      end if

      max_score_select_lo = intval
      max_score_select_hi = intval
      max_score_select_inc = 1

      string(1:nchar) = ' '

      call s_to_i4 ( string, intval, ierror, nchar )

      if ( ierror /= 0 ) then
        ierror = 0
        return
      end if

      max_score_select_hi = intval

      string(1:nchar) = ' '
      call s_to_i4 ( string, intval, ierror, nchar )

      if ( ierror /= 0 ) then
        ierror = 0
        return
      end if

      max_score_select_inc = intval

    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ANALYSIS_PARAMS_GET - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized choice!'
    ierror = 1
    return

  end if

  return
end
subroutine analysis_params_print ( debug_unit, max_score_neg_select_hi, &
  max_score_neg_select_inc, max_score_neg_select_lo, max_score_select_hi, &
  max_score_select_inc, max_score_select_lo, score_cutoff, selection )

!*****************************************************************************80
!
!! ANALYSIS_PARAMS_PRINT prints the analysis parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DEBUG_UNIT, the unit number of the debug file.
!
!    Input, integer ( kind = 4 ) MAX_SCORE_NEG_SELECT_HI, MAX_SCORE_NEG_SELECT_INC,
!    MAX_SCORE_NEG_SELECT_LO, the maximum, increment, and minimum for the
!    DO loop that determines the number of scores with negative results
!    that will be considered.
!
!    Input, integer ( kind = 4 ) MAX_SCORE_SELECT_HI, MAX_SCORE_SELECT_INC,
!    MAX_SCORE_SELECT_LO, the maximum, increment, and minimum for the DO
!    loop that determines the number of scores that will be considered.
!
!    Input, real ( kind = 8 ) SCORE_CUTOFF, the minimum score to be considered.
!
!    Input, character ( len = * ) SELECTION, the selection criterion:
!    'ALL',               select all data;
!    'SCORE_MIN'          select all data with scores no less than
!                         SCORE_CUTOFF;
!    'POSITIVE_ALL'       select all data with scores no less than the lowest
!                         score of a positive data item;
!    'NEGATIVE_RESTRICT'  select all data with scores no less than that of
!                         the MAX_SCORE_NEG_SELECT-th negative item;
!    'TOTAL_RESTRICT'     select the MAX_SCORE_SELECT highest ranking items.
!
  implicit none

  integer ( kind = 4 ) debug_unit
  integer ( kind = 4 ) max_score_neg_select_hi
  integer ( kind = 4 ) max_score_neg_select_inc
  integer ( kind = 4 ) max_score_neg_select_lo
  integer ( kind = 4 ) max_score_select_hi
  integer ( kind = 4 ) max_score_select_inc
  integer ( kind = 4 ) max_score_select_lo
  real ( kind = 8 ) score_cutoff
  character ( len = 35 ) selection

  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) 'ANALYSIS_PARAMS_PRINT:'

  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) '  Data selection criterion:'
  write ( debug_unit, * ) ' '

  if ( selection == 'ALL' ) then

    write ( debug_unit, * ) '    Consider all scores.'

  else if ( selection == 'SCORE_MIN' ) then

    write ( debug_unit, * ) '    Only consider data with scores of at least ', &
      score_cutoff

  else if ( selection == 'POSITIVE_ALL' ) then

    write ( debug_unit, * ) '    Only consider data down to the lowest ' // &
      'scoring positive data item.'

  else if ( selection == 'NEGATIVE_RESTRICT' ) then

    if ( max_score_neg_select_lo == max_score_neg_select_hi ) then
      write ( debug_unit, * ) '    Only consider data down to a total number ' &
        // 'of negative data items of ', max_score_neg_select_lo
    else
      write ( debug_unit, * ) '    Only consider data down to a total number ' &
        // 'of negative data items of MAX_SCORE_NEG_SELECT,'
      write ( debug_unit, * ) '  where MAXSCORE_NEG_SELECT is controlled by'
      write ( debug_unit, * ) '  a DO loop from ', max_score_neg_select_lo, &
        ' to ', max_score_neg_select_hi, ' by ', max_score_neg_select_inc
    end if

  else if ( selection == 'TOTAL_RESTRICT' ) then

    if ( max_score_select_lo == max_score_select_hi ) then
      write ( debug_unit, * ) '  Only consider data down to a total number ' &
        // 'of data items of ', max_score_select_lo
    else
      write ( debug_unit, * ) '  Only consider data down to a total number ' &
        // 'of data items of MAX_SCORE_SELECT,'
      write ( debug_unit, * ) '  where MAXSCORE_SELECT is controlled by'
      write ( debug_unit, * ) '  a DO loop from ', max_score_select_lo, &
        ' to ', max_score_select_hi, ' by ', max_score_select_inc
    end if

  end if

  return
end
subroutine analysis_post ( graph_file_unit, plot )

!*****************************************************************************80
!
!! ANALYSIS_POST does some post-analysis tasks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GRAPH_FILE_UNIT, the unit number associated
!    with the graph file.
!
!    Input, logical PLOT, is TRUE if plots are desired.
!
  implicit none

  integer ( kind = 4 ) graph_file_unit
  logical              plot

  if ( .not. plot ) then
    return
  end if

  write ( graph_file_unit, '(a)' ) 'endfile'

  close ( unit = graph_file_unit )

  return
end
subroutine analysis_pre ( graph_file_name, graph_file_unit, ierror, plot )

!*****************************************************************************80
!
!! ANALYSIS_PRE carries out some pre-analysis tasks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) GRAPH_FILE_NAME, the name of the graph file
!    to be created.
!
!    Output, integer ( kind = 4 ) GRAPH_FILE_UNIT, the unit number associated with the
!    graph file.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    nonzero, an error occurred.
!
!    Input, logical PLOT, is TRUE if plots are desired.
!
  implicit none

  character ( len = * ) graph_file_name
  integer ( kind = 4 ) graph_file_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  logical plot

  ierror = 0

  if ( .not. plot ) then
    return
  end if
!
!  Get a free unit number.
!
  call get_unit ( graph_file_unit )

  if ( graph_file_unit == 0 ) then
    ierror = 99
    return
  end if
!
!  Open the new graphics file.
!
  open ( unit = graph_file_unit, file = graph_file_name, status = 'replace', &
    access = 'sequential', form = 'formatted', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PN_GRAPH_FILE_WRITE - Warning!'
    write ( *, '(a)' ) '  Could not open the graph file ' // &
      trim ( graph_file_name )
    write ( *, '(a)' ) '  No graph file will be created.'
    return
  end if

  write ( graph_file_unit, '(a)' ) '# ' // trim ( graph_file_name ) // &
    ' created by ROC.'
  write ( graph_file_unit, '(a)' ) '#'
  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) 'file ' // trim ( graph_file_name )
  write ( graph_file_unit, '(a)' ) '  space 0.0 0.0 8.5 11.0'

  return
end
subroutine analysis_print ( debug_unit, identifier, max_score, num_score, &
  num_score_neg, num_score_neg_select, num_score_pos, num_score_pos_select, &
  num_score_select, result, roc_val, score, score_max_selected, &
  score_min_selected )

!*****************************************************************************80
!
!! ANALYSIS_PRINT prints the analysis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DEBUG_UNIT, the unit number of the debug file.
!
!    Input, character ( len = 35 ) IDENTIFIER(MAX_SCORE), a list of the
!    (unique) identifiers for which a score was computed.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Input, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Input, integer ( kind = 4 ) NUM_SCORE_NEG, the number of scores associated
!    with negative results.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_NEG_SELECT, the number of selected scores
!    associated with negative results.
!
!    Input, integer ( kind = 4 ) NUM_SCORE_POS, the number of scores associated
!    with positive results.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_POS_SELECT, the number of selected scores
!    associated with positive results.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_SELECT, the number of selected scores.
!
!    Input, integer ( kind = 4 ) RESULT(MAX_SCORE), the result code.  RESULT(I) is:
!      J, item I is selected as the J-th POSITIVE representative;
!      0, item I is omitted;
!    - J, item I is selected as the J-th NEGATIVE representative.
!
!    Input, real ( kind = 8 ) ROC_VAL, the value of the ROC integral.
!
!    Input, real ( kind = 8 ) SCORE(MAX_SCORE), the scores.
!
!    Input, real ( kind = 8 ) SCORE_MAX_SELECTED, SCORE_MIN_SELECTED, the maximum and
!    minimum scores to be plotted.
!
  implicit none

  integer ( kind = 4 ) max_score

  integer ( kind = 4 ) debug_unit
  integer ( kind = 4 ) i
  character ( len = 35 ) identifier(max_score)
  integer ( kind = 4 ) num_score
  integer ( kind = 4 ) num_score_neg
  integer ( kind = 4 ) num_score_neg_select
  integer ( kind = 4 ) num_score_pos
  integer ( kind = 4 ) num_score_pos_select
  integer ( kind = 4 ) num_score_select
  integer ( kind = 4 ) result(max_score)
  real ( kind = 8 ) roc_val
  real ( kind = 8 ) score(max_score)
  real ( kind = 8 ) score_max_selected
  real ( kind = 8 ) score_min_selected

  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) 'ANALYSIS_PRINT:'

  write ( debug_unit, * ) ' '
  write ( debug_unit, '(a)' ) '  Result  Identifier          Score'
  write ( debug_unit, * ) ' '
  do i = 1, num_score
    write ( debug_unit, '(i8,2x,a30,2x,g14.6)' ) result(i), identifier(i), &
      score(i)
  end do

  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) '  The ROC area integral is ', roc_val
  write ( debug_unit, * ) '  The range of the selected data was ', &
    score_min_selected, ' to ' , score_max_selected
  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) '  Selected ', num_score_pos_select, &
    ' positive items out of ', num_score_pos
  write ( debug_unit, * ) '  Selected ', num_score_neg_select, &
    ' negative items out of ', num_score_neg
  write ( debug_unit, * ) '  Selected ', num_score_select, &
    ' items out of          ', num_score

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
  integer ( kind = 4 ) itemp

  itemp = iachar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!     CH  DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character CH, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  If CH was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then

    digit = iachar ( ch ) - 48

  else if ( ch == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine chvec_permute ( n, a, p )

!*****************************************************************************80
!
!! CHVEC_PERMUTE permutes a character vector in place.
!
!  Discussion:
!
!    This routine permutes an array of character "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,   4,   5,   1,   3 )
!      A = (  'B', 'D', 'E', 'A', 'C' )
!
!    Output:
!
!      A    = ( 'A', 'B', 'C', 'D', 'E' ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, character A(N), the array to be permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  implicit none

  integer ( kind = 4 ) n

  character a(n)
  character a_temp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) get
  integer ( kind = 4 ) put
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHVEC_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i8)' ) '  array is missing the value ', ierror
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = -p(istart)
      cycle

    else

      a_temp = a(istart)
      get = istart
!
!  Copy the new value into the vacated entry.
!
      do

        put = get
        get = p(get)

        p(put) = -p(put)

        if ( get < 1 .or. n < get ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CHVEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  "get" character is out of bounds.'
          stop
        end if

        if ( get == istart ) then
          a(put) = a_temp
          exit
        end if

        a(put) = a(get)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine chvec_reverse ( n, x )

!*****************************************************************************80
!
!! CHVEC_REVERSE reverses the elements of a character vector.
!
!  Example:
!
!    Input:
!
!      N = 4, X = ( 'L', 'I', 'V', 'E' ).
!
!    Output:
!
!      X = ( 'E', 'V', 'I', 'L' ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, character X(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  character cval
  integer ( kind = 4 ) i
  character x(n)

  do i = 1, n/2
    cval = x(i)
    x(i) = x(n+1-i)
    x(n+1-i) = cval
  end do

  return
end
subroutine data_count ( good, identifier, max_score, num_good, num_score, &
  num_score_neg, num_score_pos, score, score_max, score_max_neg, &
  score_max_pos, score_min, score_min_neg, score_min_pos )

!*****************************************************************************80
!
!! DATA_COUNT counts the positive and negative data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 35 ) GOOD(MAX_SCORE), a list of the NUM_GOOD
!    "good" identifiers.
!
!    Input, character ( len = 35 ) IDENTIFIER(MAX_SCORE), a list of the
!    NUM_SCORE identifiers for which a score was computed.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scored items.
!
!    Input, integer ( kind = 4 ) NUM_GOOD, the number of good identifiers.
!
!    Input, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_NEG, the number of scores associated with
!    negative results.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_POS, the number of scores associated with
!    positive results.
!
  implicit none

  integer ( kind = 4 ) max_score

  character ( len = 35 ) good(max_score)
  integer ( kind = 4 ) i
  character ( len = 35 ) identifier(max_score)
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) num_good
  integer ( kind = 4 ) num_score
  integer ( kind = 4 ) num_score_neg
  integer ( kind = 4 ) num_score_pos
  real ( kind = 8 ) score(max_score)
  real ( kind = 8 ) score_max
  real ( kind = 8 ) score_max_neg
  real ( kind = 8 ) score_max_pos
  real ( kind = 8 ) score_min
  real ( kind = 8 ) score_min_neg
  real ( kind = 8 ) score_min_pos

  num_score_neg = 0
  num_score_pos = 0

  score_max = 0.0D+00
  score_max_neg = 0.0D+00
  score_max_pos = 0.0D+00
  score_min = 0.0D+00
  score_min_neg = 0.0D+00
  score_min_pos = 0.0D+00

  do i = 1, num_score

    call sveci_search_binary_a ( num_good, good, identifier(i), indx )

    if ( i == 1 ) then
      score_max = score(i)
      score_min = score(i)
    else
      score_max = max ( score_max, score(i) )
      score_min = min ( score_min, score(i) )
    end if

    if ( indx == 0 ) then

      num_score_neg = num_score_neg + 1

      if ( num_score_neg == 1 ) then
        score_max_neg = score(i)
        score_min_neg = score(i)
      else
        score_max_neg = max ( score_max_neg, score(i) )
        score_min_neg = min ( score_min_neg, score(i) )
      end if

    else

      num_score_pos = num_score_pos + 1

      if ( num_score_pos == 1 ) then
        score_max_pos = score(i)
        score_min_pos = score(i)
      else
        score_max_pos = max ( score_max_pos, score(i) )
        score_min_pos = min ( score_min_pos, score(i) )
      end if

    end if

  end do

  return
end
subroutine data_count_print ( debug_unit, num_score, num_score_neg, &
  num_score_pos, score_max, score_max_neg, score_max_pos, score_min, &
  score_min_neg, score_min_pos )

!*****************************************************************************80
!
!! DATA_COUNT_PRINT prints the data counts.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_NEG, the number of scores associated with
!    negative results.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_POS, the number of scores associated with
!    positive results.
!
  implicit none

  integer ( kind = 4 ) debug_unit
  integer ( kind = 4 ) num_score
  integer ( kind = 4 ) num_score_neg
  integer ( kind = 4 ) num_score_pos
  real ( kind = 8 ) score_max
  real ( kind = 8 ) score_max_neg
  real ( kind = 8 ) score_max_pos
  real ( kind = 8 ) score_min
  real ( kind = 8 ) score_min_neg
  real ( kind = 8 ) score_min_pos

  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) 'DATA_COUNT_PRINT'
  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) '  Positive data'
  write ( debug_unit, * ) '    Number        ', num_score_pos
  write ( debug_unit, * ) '    Minimum score ', score_min_pos
  write ( debug_unit, * ) '    Maximum score ', score_max_pos
  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) '  Negative data'
  write ( debug_unit, * ) '    Number        ', num_score_neg
  write ( debug_unit, * ) '    Minimum score ', score_min_neg
  write ( debug_unit, * ) '    Maximum score ', score_max_neg
  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) '  Total'
  write ( debug_unit, * ) '    Number        ', num_score
  write ( debug_unit, * ) '    Minimum score ', score_min
  write ( debug_unit, * ) '    Maximum score ', score_max

  return
end
subroutine data_get ( debug, debug_unit, identifier, ierror, indx, max_score, &
  num_score, score, search_file_name, search_type )

!*****************************************************************************80
!
!! DATA_GET retrieves the identifiers, scores, and good identifiers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical DEBUG, is TRUE if debug output is desired.
!
!    Input, integer ( kind = 4 ) DEBUG_UNIT, the unit number of the debug file.
!
!    Output, character ( len = 35 ) IDENTIFIER(MAX_SCORE), a list of the
!    (unique) identifiers for which a score was computed.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    nonzero, an error occurred.
!
!    Workspace, integer INDX(MAX_SCORE).
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Input, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Output, real ( kind = 8 ) SCORE(MAX_SCORE), the scores.
!
!    Input, character ( len = * ) SEARCH_FILE_NAME, the name of the search file.
!
!    Input, character ( len = 10 ) SEARCH_TYPE, identifies the search used to
!    create the search file:
!    'BLAST', 'FASTA', 'GENERIC', 'MAXSEGS', 'PEARSON'.
!
  implicit none

  integer ( kind = 4 ) max_score

  logical debug
  integer ( kind = 4 ) debug_unit
  character ( len = 35 ) identifier(max_score)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(max_score)
  integer ( kind = 4 ) num_score
  real ( kind = 8 ) score(max_score)
  integer ( kind = 4 ) score_order
  character ( len = * ) search_file_name
  character ( len = 10 ) search_type
!
!  Read the data.
!
  if ( search_type == 'BLAST' ) then

    call score_blast_file_read ( identifier, ierror, max_score, num_score, &
      score, search_file_name )

  else if ( search_type == 'FASTA' ) then

    call score_fasta_file_read ( identifier, ierror, max_score, num_score, &
      score, search_file_name )

  else if ( search_type == 'GENERIC' ) then

    call score_generic_file_read ( identifier, ierror, max_score, num_score, &
      score, search_file_name )

  else if ( search_type == 'MAXSEGS' ) then

    call score_maxsegs_file_read ( identifier, ierror, max_score, num_score, &
      score, search_file_name )

  else if ( search_type == 'PEARSON' ) then

    call score_pearson_file_read ( identifier, ierror, max_score, num_score, &
      score, search_file_name )

  end if
!
!  Check before proceeding.
!
  if ( ierror /= 0 ) then
    return
  end if

  if ( num_score <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_GET - Fatal error!'
    write ( *, '(a)' ) '  There are no results in the file.'
    write ( *, '(a)' ) '  The ROC cannot be computed for this data.'
    ierror = 1
    return
  end if

  if ( max_score < num_score ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_GET - Warning!'
    write ( *, '(a)' ) '  Too much data in the file!'
    write ( *, '(a,i8)' ) '  NUM_SCORE = ', num_score
    write ( *, '(a,i8)' ) '  MAX_SCORE = ', max_score
    ierror = 1
    return
  end if
!
!  Determine if the scores are sorted ascending, sorted descending, or other.
!
  call r8vec_order_type ( num_score, score, score_order )
!
!  If descending, all is well...
!
  if ( score_order == 0 .or. score_order == 3 .or. score_order == 4 ) then

    if ( debug ) then
      write ( debug_unit, * ) ' '
      write ( debug_unit, * ) 'DATA_GET:'
      write ( debug_unit, * ) '  The scores are already in descending order.'
    end if
!
!  ...else, if ascending, make them descending by reversal...
!
  else if ( score_order == 1 .or. score_order == 2 ) then

    if ( debug ) then
      write ( debug_unit, * ) ' '
      write ( debug_unit, * ) 'DATA_GET:'
      write ( debug_unit, * ) '  The scores are in ascending order.'
      write ( debug_unit, * ) '  A reversal will be applied.'
    end if

    call r8vec_reverse ( num_score, score )

    call svec_reverse ( num_score, identifier )
!
!  ...else, if not sorted at all, sort them into descending order.
!
  else

    if ( debug ) then
      write ( debug_unit, * ) ' '
      write ( debug_unit, * ) 'DATA_GET:'
      write ( debug_unit, * ) '  The scores are not in order.'
      write ( debug_unit, * ) '  The scores will be ascending'
      write ( debug_unit, * ) '  sorted, then reversed.'
    end if

    call r8vec_sort_heap_index_a ( num_score, score, indx )

    call r8vec_permute ( num_score, score, indx )

    call r8vec_reverse ( num_score, score )

    call svec_permute ( num_score, identifier, indx )

    call svec_reverse ( num_score, identifier )

  end if

  return
end
subroutine data_params_get ( ierror, search_file_name, search_type )

!*****************************************************************************80
!
!! DATA_PARAMS_GET gets the data parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Output, character ( len = * ) SEARCH_FILE_NAME, the name of the search
!    file.
!
!    Output, character ( len = 10 ) SEARCH_TYPE, identifies the search used to
!    create the search file:
!    'BLAST', 'FASTA', 'GENERIC', 'MAXSEGS', 'PEARSON'.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 256 ) search_file_name
  character ( len = 10 ) search_type

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter the name of the search file:'
  read ( *, '(a)', iostat = ios ) search_file_name

  if ( ios /= 0 ) then
    ierror = 1
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter the number or initial of the search type:'
  write ( *, '(a)' ) '  1 BLAST'
  write ( *, '(a)' ) '  2 FASTA'
  write ( *, '(a)' ) '  3 GENERIC'
  write ( *, '(a)' ) '  4 MAXSEGS'
  write ( *, '(a)' ) '  5 PEARSON'
  write ( *, '(a)' ) '  0 QUIT'

  read ( *, '(a)', iostat = ios ) search_type

  if ( ios /= 0 ) then
    ierror = 1
    return
  end if

  call s_cap ( search_type )

  if ( search_type == '0' .or. search_type(1:1) == 'Q' ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_PARAMS_GET:'
    write ( *, '(a)' ) '  User-requested stop.'
    stop

  else if ( search_type == '1' .or. search_type(1:1) == 'B' ) then

    search_type = 'BLAST'

  else if ( search_type == '2' .or. search_type(1:1) == 'F' ) then

    search_type = 'FASTA'

  else if ( search_type == '3' .or. search_type(1:1) == 'G' ) then

    search_type = 'GENERIC'

  else if ( search_type == '4' .or. search_type(1:1) == 'M' ) then

    search_type = 'MAXSEGS'

  else if ( search_type == '5' .or. search_type(1:1) == 'P' ) then

    search_type = 'PEARSON'

  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARAMS_GET - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized search type!'
    ierror = 1
    return
  end if

  return
end
subroutine data_params_print ( debug_unit, search_file_name, search_type )

!*****************************************************************************80
!
!! DATA_PARAMS_PRINT prints the data parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DEBUG_UNIT, the unit number of the debug file.
!
!    Input, character ( len = * ) SEARCH_FILE_NAME, the name of the search file.
!
!    Input, character ( len = 10 ) SEARCH_TYPE, identifies the search used to
!    create the search file:
!    'BLAST', 'FASTA', 'GENERIC', 'MAXSEGS', 'PEARSON'.
!
  implicit none

  integer ( kind = 4 ) debug_unit
  character ( len = 256 ) search_file_name
  character ( len = 10 ) search_type

  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) 'DATA_PARAMS_PRINT:'
  write ( debug_unit, * ) '  Search file = ' // trim ( search_file_name )
  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) '  The search type was: ' // trim ( search_type )
  write ( debug_unit, * ) ' '

  return
end
subroutine data_print ( debug_unit, identifier, max_score, num_score, score )

!*****************************************************************************80
!
!! DATA_PRINT prints the data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DEBUG_UNIT, the unit number of the debug file.
!
!    Input, character ( len = 35 ) IDENTIFIER(MAX_SCORE), a list of the
!    (unique) identifiers for which a score was computed.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Input, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Input, real ( kind = 8 ) SCORE(MAX_SCORE), the scores.
!
  implicit none

  integer ( kind = 4 ) max_score

  integer ( kind = 4 ) debug_unit
  integer ( kind = 4 ) i
  character ( len = 35 ) identifier(max_score)
  integer ( kind = 4 ) num_score
  real ( kind = 8 ) score(max_score)

  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) 'DATA_PRINT:'
  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) '  Number of scored data items:  ',  num_score
  write ( debug_unit, * ) ' '
  write ( debug_unit, '(a)' ) ' Index          Identifier          Score'
  write ( debug_unit, * ) ' '
  do i = 1, num_score
    write ( debug_unit, '(i6,2x,a35,2x,g14.6)' ) i, identifier(i), score(i)
  end do

  return
end
subroutine data_select ( good, identifier, max_score, max_score_neg_select, &
  max_score_select, num_good, num_score, num_score_neg_select, num_score_pos, &
  num_score_pos_select, num_score_select, num_score_skip, result, score, &
  score_cutoff, selection )

!*****************************************************************************80
!
!! DATA_SELECT selects the data to be analyzed.
!
!  Discussion:
!
!    Data item I has an IDENTIFIER and a SCORE.
!    Data items are listed in descending order of SCORE.
!
!    Starting with the highest scoring item I = 1, the scores
!    are considered for selection, based on one of several
!    selection criteria.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 35 ) GOOD(MAX_SCORE), a list of the NUM_GOOD
!    "good" identifiers.
!
!    Input, character ( len = 35 ) IDENTIFIER(MAX_SCORE), a list of the
!    NUM_SCORE identifiers for which a score was computed.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scored items.
!
!    Input, integer ( kind = 4 ) MAX_SCORE_NEG_SELECT, the maximum number of scores with
!    negative results that will be considered.
!
!    Input, integer ( kind = 4 ) MAX_SCORE_SELECT, the maximum number of scores that
!    will be considered.
!
!    Input, integer ( kind = 4 ) NUM_GOOD, the number of good identifiers.
!
!    Input, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_NEG_SELECT, the number of selected scores
!    associated with negative results.
!
!    Input, integer ( kind = 4 ) NUM_SCORE_POS, the number of scores associated with
!    positive results.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_POS_SELECT, the number of selected scores
!    associated with positive results.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_SELECT, the number of selected scores.
!
!    Output, integer ( kind = 4 ) NUM_SCORE_SKIP, the number of unselected scores.
!
!    Output, integer ( kind = 4 ) RESULT(MAX_SCORE), the result code.  RESULT(I) is:
!      J, item I is selected as the J-th POSITIVE representative;
!      0, item I is omitted;
!    - J, item I is selected as the J-th NEGATIVE representative.
!
!    Input, real ( kind = 8 ) SCORE(MAX_SCORE), the scores.
!
!    Input, real ( kind = 8 ) SCORE_CUTOFF, the minimum score to be considered.
!
!    Input, character ( len = * ) SELECTION, the selection criterion:
!    'ALL',               select all data;
!    'SCORE_MIN'          select all data with scores no less than
!                         SCORE_CUTOFF;
!    'POSITIVE_ALL'       select all data with scores no less than the lowest
!                         score of a positive data item;
!    'NEGATIVE_RESTRICT'  select all data with scores no less than that of
!                         the MAX_SCORE_NEG_SELECT-th negative item;
!    'TOTAL_RESTRICT'     select the MAX_SCORE_SELECT highest ranking items.
!
  implicit none

  integer ( kind = 4 ) max_score

  character ( len = 35 ) good(max_score)
  integer ( kind = 4 ) i
  character ( len = 35 ) identifier(max_score)
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) max_score_neg_select
  integer ( kind = 4 ) max_score_select
  integer ( kind = 4 ) num_good
  integer ( kind = 4 ) num_score
  integer ( kind = 4 ) num_score_neg_select
  integer ( kind = 4 ) num_score_pos
  integer ( kind = 4 ) num_score_pos_select
  integer ( kind = 4 ) num_score_select
  integer ( kind = 4 ) num_score_skip
  integer ( kind = 4 ) result(max_score)
  real ( kind = 8 ) score(max_score)
  real ( kind = 8 ) score_cutoff
  character ( len = 35 ) selection
  logical selected

  num_score_neg_select = 0
  num_score_pos_select = 0
  num_score_select = 0
  num_score_skip = 0

  do i = 1, num_score

    selected = .true.

    if ( selection == 'ALL' ) then

    else if ( selection == 'SCORE_MIN' ) then

      if ( score(i) < score_cutoff ) then
        selected = .false.
      end if

    else if ( selection == 'POSITIVE_ALL' ) then

      if ( num_score_pos <= num_score_pos_select ) then
        selected = .false.
      end if

    else if ( selection == 'NEGATIVE_RESTRICT' ) then

      if ( max_score_neg_select <= num_score_neg_select ) then
        selected = .false.
      end if

    else if ( selection == 'TOTAL_RESTRICT' ) then

      if ( max_score_select <= num_score_select ) then
        selected = .false.
      end if

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DATA_SELECT - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of SELECTION.'
      stop

    end if

    if ( selected ) then

      num_score_select = num_score_select + 1

      call sveci_search_binary_a ( num_good, good, identifier(i), indx )

      if ( indx == 0 ) then
        num_score_neg_select = num_score_neg_select + 1
        result(i) = - num_score_neg_select
      else
        num_score_pos_select = num_score_pos_select + 1
        result(i) = num_score_pos_select
      end if

    else

      result(i) = 0
      num_score_skip = num_score_skip + 1

    end if

  end do

  return
end
subroutine digit_to_ch ( digit, ch )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Discussion:
!
!    Instead of CHAR, we now use the ACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!    DIGIT   CH
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character CH, the corresponding character.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    ch = achar ( digit + 48 )

  else

    ch = '*'

  end if

  return
end
subroutine file_advance_to_string ( iunit, string, line, ierror )

!*****************************************************************************80
!
!! FILE_ADVANCE_TO_STRING searches ahead in a text file for a string.
!
!  Discussion:
!
!    The file should already have been opened.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the unit number associated with the open file.
!
!    Input, character ( len = * ) STRING, a string to search for.
!
!    Output, character ( len = * ) LINE:
!    If IERROR = 0, the line of the file that was just read,
!      and which contains the string.
!    If IERROR = 1, then LINE is blank.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error, the string was found.
!    1, error, the end of the file was reached.
!
  implicit none

  logical, parameter :: debug = .false.
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) lens
  character ( len = * ) line
  integer ( kind = 4 ) num_read
  character ( len = * ) string

  ierror = 0
  num_read = 0
  lens = len_trim ( string )

  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    num_read = num_read + 1

    if ( index ( line, string(1:lens) ) /= 0 ) then
      return
    end if

  end do

  line = ' '
  ierror = 1

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ADVANCE_TO_STRING - Warning!'
    write ( *, '(a)' ) '  Did not find string:'
    write ( *, '(a)' ) trim ( string )
    write ( *, '(a,i8)' ) '  Number of lines read was ', num_read
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

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
subroutine good_file_read ( good, good_file_name, ierror, max_score, num_good )

!*****************************************************************************80
!
!! GOOD_FILE_READ reads a  "good" ID file.
!
!  Discussion:
!
!    The "good" ID file contains a list of all the identifiers for which
!    a positive score should be recorded.  This file can be used in
!    conjunction with an "IS" or "identifier score" file, to compute the
!    result, which should be positive if the identifier is on the good ID list,
!    negative if it is not, and 0 if it was omitted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 35 ) GOOD(MAX_SCORE), a list of identifiers
!    of all the objects which should receive a positive score.
!
!    Input, character ( len = * ) GOOD_FILE_NAME, the name of the good ID
!    file to be read.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Output, integer ( kind = 4 ) NUM_GOOD, the number of "good" identifiers.
!
  implicit none

  integer ( kind = 4 ) max_score

  character ( len = 35 ) good(max_score)
  character ( len = * ) good_file_name
  integer ( kind = 4 ) good_file_unit
  integer ( kind = 4 ), parameter :: id_column = 1
  character ( len = 35 ) id_temp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 80 ) line
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) num_good
  integer ( kind = 4 ) num_lines

  ierror = 0
  num_lines = 0
  num_good = 0

  call get_unit ( good_file_unit )

  open ( unit = good_file_unit, file = good_file_name, status = 'old', &
    access = 'sequential', form = 'formatted', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GOOD_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the "good ID" file:'
    write ( *, '(a)' ) trim ( good_file_name )
    return
  end if

  do

    read ( good_file_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    num_lines = num_lines + 1
!
!  Comment lines begin with '#' and are ignored.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if
!
!  Blank lines are ignored.
!
    if ( line == ' ' ) then
      cycle
    end if
!
!  Read the identifier.
!
    call s_word_find ( line, id_column, id_temp, nchar )

    num_good = num_good + 1

    if ( num_good <= max_score ) then
      good(num_good) = id_temp
    end if

  end do

  close ( unit = good_file_unit )

  if ( num_good <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GOOD_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  There are no good ID''s in the file.'
    ierror = 4
  end if

  if ( max_score < num_good ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GOOD_FILE_READ - Warning!'
    write ( *, '(a)' ) '  Too much data in the file!'
    write ( *, '(a,i8)' ) '  NUM_GOOD = ', num_good
    write ( *, '(a,i8)' ) '  MAX_SCORE = ', max_score
    ierror = 5
  end if

  return
end
subroutine good_get ( good, good_file_name, ierror, max_score, num_good )

!*****************************************************************************80
!
!! GOOD_GET retrieves the good identifiers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 35 ) GOOD(MAX_SCORE), a list of identifiers
!    of all the objects which should receive a positive score.
!
!    Input, character ( len = * ) GOOD_FILE_NAME, the name of the good ID
!    file to be read.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    nonzero, an error occurred.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Input, integer ( kind = 4 ) NUM_GOOD, the number of "good" identifiers.
!
  implicit none

  integer ( kind = 4 ) max_score

  character ( len = 35 ) good(max_score)
  character ( len = * ) good_file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) num_good

  call good_file_read ( good, good_file_name, ierror, max_score, num_good )

  call sveci_sort_heap_a ( num_good, good )

  return
end
subroutine good_print ( debug_unit, good, good_file_name, max_score, num_good )

!*****************************************************************************80
!
!! GOOD_PRINT prints the good identifiers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DEBUG_UNIT, the unit number of the debug file.
!
!    Input, character ( len = 35 ) GOOD(MAX_SCORE), a list of identifiers
!    of all the objects which should receive a positive score.
!
!    Input, character ( len = * ) GOOD_FILE_NAME, the name of the good ID
!    file to be read.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Input, integer ( kind = 4 ) NUM_GOOD, the number of "good" identifiers.
!
  implicit none

  integer ( kind = 4 ) max_score

  integer ( kind = 4 ) debug_unit
  integer ( kind = 4 ) i
  character ( len = 35 ) good(max_score)
  character ( len = * ) good_file_name
  integer ( kind = 4 ) num_good

  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) 'GOOD_PRINT:'
  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) '  "Good" identifier file = ' // &
    trim ( good_file_name )
  write ( debug_unit, * ) ' '
  write ( debug_unit, * ) '  Number of "good" identifiers: ', num_good
  write ( debug_unit, * ) ' '
  write ( debug_unit, '(a)' ) '  Good identifiers'
  write ( debug_unit, * ) ' '
  do i = 1, num_good
    write ( debug_unit, '(i6,2x,a35)' ) i, good(i)
  end do

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine i4_to_s_left ( intval, s )

!*****************************************************************************80
!
!! I4_TO_S_LEFT converts an I4 to a left-justified string.
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  1
!        -1  -1
!         0  0
!      1952  1952
!    123456  123456
!   1234567  ******  <-- Not enough room!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be left-justified.  If there is not enough space,
!    the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * ) s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = -ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Find the last digit of IVAL, strip it off, and stick it into the string.
!
  do

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

    if ( ival == 0 ) then
      exit
    end if

  end do
!
!  Shift the string to the left.
!
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '

  return
end
subroutine i4vec_cum ( n, a, a_cum )

!*****************************************************************************80
!
!! I4VEC_CUM computes the cumulutive sum of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Example:
!
!    Input:
!
!      A = (/ 1, 2, 3, 4 /)
!
!    Output:
!
!      A_CUM = (/ 0, 1, 3, 6, 10 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be summed.
!
!    Output, integer ( kind = 4 ) A_CUM(0:N), the cumulative sum of the entries of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_cum(0:n)
  integer ( kind = 4 ) i

  a_cum(0) = 0

  do i = 1, n
    a_cum(i) = a_cum(i-1) + a(i)
  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_permute ( n, a, p )

!*****************************************************************************80
!
!! I4VEC_PERMUTE permutes an I4VEC in place.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    This routine permutes an array of integer "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,   4,   5,   1,   3 )
!      A = (   1,   2,   3,   4,   5 )
!
!    Output:
!
!      A    = (   2,   4,   5,   1,   3 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  Pmust be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_temp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = -p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = -p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! I4VEC_REVERSE reverses the elements of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    In FORTRAN90, call I4VEC_REVERSE is equivalent to:
!
!      A(1:N) = A(N:1:-1)
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n/2
    call i4_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine i4vec_negative_index ( n, x, max_index, n_index, value_index )

!*****************************************************************************80
!
!! I4VEC_NEGATIVE_INDEX indexes negative integer vector entries.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, real ( kind = 8 ) X(N), the array to be indexed.
!
!    Input, integer ( kind = 4 ) MAX_INDEX, the maximum number of indices to find.
!
!    Output, integer ( kind = 4 ) N_INDEX, the number of negative entries of X.
!
!    Output, integer ( kind = 4 ) VALUE_INDEX(MAX_INDEX), the indices in X of negative
!    entries.
!
  implicit none

  integer ( kind = 4 ) max_index
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n_index
  integer ( kind = 4 ) value_index(max_index)
  integer ( kind = 4 ) x(n)

  n_index = 0

  do i = 1, n

    if ( x(i) < 0 ) then

      if ( max_index <= n_index ) then
        return
      end if

      n_index = n_index + 1
      value_index(n_index) = i

    end if

  end do

  return
end
subroutine i4vec_positive_index ( n, x, max_index, n_index, value_index )

!*****************************************************************************80
!
!! I4VEC_POSITIVE_INDEX indexes positive integer vector entries.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, real ( kind = 8 ) X(N), the array to be indexed.
!
!    Input, integer ( kind = 4 ) MAX_INDEX, the maximum number of indices to find.
!
!    Output, integer ( kind = 4 ) N_INDEX, the number of positive entries of X.
!
!    Output, integer ( kind = 4 ) VALUE_INDEX(MAX_INDEX), the indices in X of positive
!    entries.
!
  implicit none

  integer ( kind = 4 ) max_index
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n_index
  integer ( kind = 4 ) value_index(max_index)
  integer ( kind = 4 ) x(n)

  n_index = 0

  do i = 1, n

    if ( 0 < x(i)) then

      if ( max_index <= n_index ) then
        return
      end if

      n_index = n_index + 1
      value_index(n_index) = i

    end if

  end do

  return
end
function lgti ( strng1, strng2 )

!*****************************************************************************80
!
!! LGTI = STRNG1 is lexically greater than STRNG2.
!
!  Discussion:
!
!    The comparison is done in a case-insensitive way.
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
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical LGTI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical lgti
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( lgt ( s1, s2 ) ) then
      lgti = .true.
      return
    else if ( llt ( s1, s2 ) ) then
      lgti = .false.
      return
    end if

  end do

  if ( len1 <= len2 ) then
    lgti = .false.
  else
    lgti = .true.
  end if

  return
end
function llei ( strng1, strng2 )

!*****************************************************************************80
!
!! LLEI = STRNG1 is lexically less than or equal to STRNG2.
!
!  Discussion:
!
!    The comparison is done in a case-insensitive way.
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
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical LLEI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical llei
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( llt ( s1, s2 ) ) then
      llei = .true.
      return
    else if ( lgt ( s1, s2 ) ) then
      llei = .false.
      return
    end if

  end do

  if ( len1 <= len2 ) then
    llei = .true.
  else
    llei = .false.
  end if

  return
end
function llti ( strng1, strng2 )

!*****************************************************************************80
!
!! LLTI = STRNG1 is lexically less than STRNG2.
!
!  Discussion:
!
!    The comparison is done in a case-insensitive way.
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
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical LLTI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical llti
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( llt ( s1, s2 ) ) then
      llti = .true.
      return
    else if ( lgt ( s1, s2 ) ) then
      llti = .false.
      return
    end if

  end do

  if ( len1 < len2 ) then
    llti = .true.
  else
    llti = .false.
  end if

  return
end
subroutine params_default ( debug, debug_file_name, good_file_name, &
  graph_file_name, max_score_neg_select_hi, max_score_neg_select_inc, &
  max_score_neg_select_lo, max_score_select_hi, max_score_select_inc, &
  max_score_select_lo, plot, score_cutoff, search_file_name, search_type, &
  selection )

!*****************************************************************************80
!
!! PARAMS_DEFAULT sets the parameters to default values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical DEBUG, is TRUE if debugging output is desired.
!
!    Output, character ( len = * ) GOOD_FILE_NAME, the name of the good ID
!    file to be read.
!
!    Output, character ( len = * ) GRAPH_FILE_NAME, the name of the graph file
!    to be created.
!
!    Output, integer ( kind = 4 ) MAX_SCORE_NEG_SELECT_HI, MAX_SCORE_NEG_SELECT_INC,
!    MAX_SCORE_NEG_SELECT_LO, the maximum, increment, and minimum for the
!    DO loop that determines the number of scores with negative results
!    that will be considered.
!
!    Output, integer ( kind = 4 ) MAX_SCORE_SELECT_HI, MAX_SCORE_SELECT_INC,
!    MAX_SCORE_SELECT_LO, the maximum, increment, and minimum for the DO
!    loop that determines the number of scores that will be considered.
!
!    Output, logical PLOT, is TRUE if plots are desired.
!
!    Output, real ( kind = 8 ) SCORE_CUTOFF, the minimum score to be considered.
!
!    Output, character ( len = * ) SEARCH_FILE_NAME, the name of the search
!    file.
!
!    Output, character ( len = 10 ) SEARCH_TYPE, identifies the search used to
!    create the search file:
!    'BLAST', 'FASTA', 'GENERIC', 'MAXSEGS', 'PEARSON'.
!
!    Output, character ( len = * ) SELECTION, the selection criterion:
!    'ALL',               select all data;
!    'SCORE_MIN'          select all data with scores no less than
!                         SCORE_CUTOFF;
!    'POSITIVE_ALL'       select all data with scores no less than the lowest
!                         score of a positive data item;
!    'NEGATIVE_RESTRICT'  select all data with scores no less than that of
!                         the MAX_SCORE_NEG_SELECT-th negative item;
!    'TOTAL_RESTRICT'     select the MAX_SCORE_SELECT highest ranking items.
!
  implicit none

  logical debug
  character ( len = 256 ) debug_file_name
  character ( len = 256 ) good_file_name
  character ( len = 256 ) graph_file_name
  integer ( kind = 4 ) max_score_neg_select_hi
  integer ( kind = 4 ) max_score_neg_select_inc
  integer ( kind = 4 ) max_score_neg_select_lo
  integer ( kind = 4 ) max_score_select_hi
  integer ( kind = 4 ) max_score_select_inc
  integer ( kind = 4 ) max_score_select_lo
  logical plot
  real ( kind = 8 ) score_cutoff
  character ( len = 256 ) search_file_name
  character ( len = 10 ) search_type
  character ( len = 35 ) selection

  debug = .true.
  debug_file_name = 'roc.debug'
  good_file_name = ' '
  graph_file_name = 'roc.plot'
  max_score_neg_select_hi = huge ( 1 )
  max_score_neg_select_inc = 1
  max_score_neg_select_lo = huge ( 1 )
  max_score_select_hi = huge ( 1 )
  max_score_select_inc = 1
  max_score_select_lo = huge ( 1 )
  plot = .true.
  score_cutoff = - huge ( 1.0D+00 )
  search_file_name = ' '
  search_type = ' '
  selection = 'ALL'

  return
end
subroutine perm_check ( n, p, ierror )

!*****************************************************************************80
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from 1
!    to N occurs among the N entries of the permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) P(N), the array to check.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifind
  integer ( kind = 4 ) iseek
  integer ( kind = 4 ) p(n)

  ierror = 0

  do iseek = 1, n

    ierror = iseek

    do ifind = 1, n
      if ( p(ifind) == iseek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      return
    end if

  end do

  return
end
subroutine pn_graph_file_write ( graph_file_unit, good_file_name, ierror, &
  max_score, num_good, num_score, num_score_neg_select, num_score_pos_select, &
  result, roc_val, search_file_name )

!*****************************************************************************80
!
!! PN_GRAPH_FILE_WRITE creates a positive/negative graphics file of the results.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) GOOD_FILE_NAME, the name of the good ID
!    file to be read.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Input, integer ( kind = 4 ) NUM_GOOD, the number of "good" identifiers.
!
!    Input, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Input, integer ( kind = 4 ) NUM_SCORE_NEG_SELECT, the number of selected scores
!    associated with negative results.
!
!    Input, integer ( kind = 4 ) NUM_SCORE_POS_SELECT, the number of selected scores
!    associated with positive results.
!
!    Input, integer ( kind = 4 ) RESULT(MAX_SCORE), the result code.  RESULT(I) is:
!      J, item I is selected as the J-th POSITIVE representative;
!      0, item I is omitted;
!    - J, item I is selected as the J-th NEGATIVE representative.
!
!    Input, real ( kind = 8 ) ROC_VAL, the value of the receiver operator characteristic.
!
!    Input, character ( len = * ) SEARCH_FILE_NAME, the name of the search file.
!
  implicit none

  integer ( kind = 4 ) max_score

  character ( len = * ) good_file_name
  integer ( kind = 4 ) graph_file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_neg
  integer ( kind = 4 ) i_pos
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) num_good
  integer ( kind = 4 ) num_score
  integer ( kind = 4 ) num_score_neg_select
  integer ( kind = 4 ) num_score_pos_select
  real ( kind = 8 ) px
  real ( kind = 8 ) pxmax
  real ( kind = 8 ) pxmin
  real ( kind = 8 ) py
  real ( kind = 8 ) pymax
  real ( kind = 8 ) pymin
  integer ( kind = 4 ) result(max_score)
  real ( kind = 8 ) roc_val
  character ( len = * ) search_file_name
  character ( len = 14 ) string
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  ierror = 0

  i_pos = 0
  i_neg = 0

  pxmin = 1.5D+00
  pxmax = 7.5D+00
  pymin = 2.0D+00
  pymax = 8.0D+00

  x = 0.0D+00
  y = 0.0D+00

  write ( graph_file_unit, '(a)' ) '  page'
  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    line_width 1'
  write ( graph_file_unit, '(a)' ) '    line_rgb 0.5 0.5 0.5'
  write ( graph_file_unit, '(a)' ) '    grid 1.5 2.0 7.5 8.0 11 11'
  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    line_width 2'
  write ( graph_file_unit, '(a)' ) '    line_rgb 0.0 0.0 1.0'

  px = pxmin + x * ( pxmax - pxmin )
  py = pymin + y * ( pymax - pymin )

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a,2g14.6)' ) '    moveto ', px, py

  do i = 1, num_score

    if ( 0  < result(i) ) then

      i_pos = i_pos + 1
      y = real ( i_pos, kind = 8 ) &
        / real ( num_score_pos_select, kind = 8 )
      py = pymin + y * ( pymax - pymin )
      write ( graph_file_unit, '(a,2g14.6)' ) '    lineto ', px, py

    else if ( result(i) < 0 ) then

      i_neg = i_neg + 1
      x = real ( i_neg, kind = 8 ) &
        / real ( num_score_neg_select, kind = 8 )
      px = pxmin + x * ( pxmax - pxmin )
      write ( graph_file_unit, '(a,2g14.6)' ) '    lineto ', px, py

    end if

  end do

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    line_width 2'
  write ( graph_file_unit, '(a)' ) '    line_rgb 0.0 0.0 0.0'

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    font_size 0.20'

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    moveto 1.1 2.0'
  write ( graph_file_unit, '(a)' ) '    label 0.0'
  write ( graph_file_unit, '(a)' ) '    moveto 1.1 3.2'
  write ( graph_file_unit, '(a)' ) '    label 0.2'
  write ( graph_file_unit, '(a)' ) '    moveto 1.1 4.4'
  write ( graph_file_unit, '(a)' ) '    label 0.4'
  write ( graph_file_unit, '(a)' ) '    moveto 1.1 5.6'
  write ( graph_file_unit, '(a)' ) '    label 0.6'
  write ( graph_file_unit, '(a)' ) '    moveto 1.1 6.8'
  write ( graph_file_unit, '(a)' ) '    label 0.8'
  write ( graph_file_unit, '(a)' ) '    moveto 1.1 8.0'
  write ( graph_file_unit, '(a)' ) '    label 1.0'

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    moveto 0.8 4.0'
  write ( graph_file_unit, '(a)' ) '    label_slant 90.0 Positive percentage'

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    moveto 1.4 1.75'
  write ( graph_file_unit, '(a)' ) '    label 0.0'
  write ( graph_file_unit, '(a)' ) '    moveto 2.6 1.75'
  write ( graph_file_unit, '(a)' ) '    label 0.2'
  write ( graph_file_unit, '(a)' ) '    moveto 3.8 1.75'
  write ( graph_file_unit, '(a)' ) '    label 0.4'
  write ( graph_file_unit, '(a)' ) '    moveto 5.0 1.75'
  write ( graph_file_unit, '(a)' ) '    label 0.6'
  write ( graph_file_unit, '(a)' ) '    moveto 6.2 1.75'
  write ( graph_file_unit, '(a)' ) '    label 0.8'
  write ( graph_file_unit, '(a)' ) '    moveto 7.4 1.75'
  write ( graph_file_unit, '(a)' ) '    label 1.0'

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    moveto 3.5 1.5'
  write ( graph_file_unit, '(a)' ) '    label Negative percentage'

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    font_size 0.20'
  write ( graph_file_unit, '(a)' ) '    moveto 2.0 9.0'
  write ( graph_file_unit, '(a)' ) '    label Receiver Operator Characteristic'
  write ( graph_file_unit, '(a)' ) '    moveto 2.0 8.5'
  write ( graph_file_unit, '(a)' ) &
    '    label Positive versus Negative Percentages'

  write ( graph_file_unit, '(a)' ) '    moveto 2.0 1.0'
  call r8_to_s_left ( roc_val, string )
  call s_blanks_delete ( string )
  write ( graph_file_unit, '(a)' ) '    label ROC integral = ' // &
    trim ( string )

  write ( graph_file_unit, '(a)' ) '    moveto 2.0 0.8'
  write ( graph_file_unit, '(a)' ) '    label Score file = ' // &
    trim ( search_file_name )

  write ( graph_file_unit, '(a)' ) '    moveto 2.0 0.6'
  call i4_to_s_left ( num_good, string )
  call s_blanks_delete ( string )
  write ( graph_file_unit, '(a)' ) '    label Good file = ' // &
    trim ( good_file_name ) // ' contains ' // trim ( string ) &
    // ' identifiers.'

  write ( graph_file_unit, '(a)' ) '    moveto 2.0 0.4'
  call i4_to_s_left ( num_score_pos_select, string )
  call s_blanks_delete ( string )
  write ( graph_file_unit, '(a)' ) '    label ' // trim ( string ) &
    // ' positive identifiers were plotted.'

  write ( graph_file_unit, '(a)' ) '    moveto 2.0 0.2'
  call i4_to_s_left ( num_score_neg_select, string )
  call s_blanks_delete ( string )
  write ( graph_file_unit, '(a)' ) '    label ' // trim ( string ) &
    // ' negative identifiers were plotted.'

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) 'endpage'

  return
end
subroutine pn_roc_int ( max_score, num_score, num_score_neg_select, &
  num_score_pos_select, result, roc_val )

!*****************************************************************************80
!
!! PN_ROC_INT computes the ROC integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Input, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Input, integer ( kind = 4 ) NUM_SCORE_NEG_SELECT, the number of selected scores
!    associated with negative results.
!
!    Input, integer ( kind = 4 ) NUM_SCORE_POS_SELECT, the number of selected scores
!    associated with positive results.
!
!    Input, integer ( kind = 4 ) RESULT(MAX_SCORE), the result code.  RESULT(I) is:
!      J, item I is selected as the J-th POSITIVE representative;
!      0, item I is omitted;
!    - J, item I is selected as the J-th NEGATIVE representative.
!
!    Output, real ( kind = 8 ) ROC_VAL, the value of the ROC.
!
  implicit none

  integer ( kind = 4 ) max_score

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_neg
  integer ( kind = 4 ) i_pos
  integer ( kind = 4 ) num_score
  integer ( kind = 4 ) num_score_neg_select
  integer ( kind = 4 ) num_score_pos_select
  integer ( kind = 4 ) result(max_score)
  real ( kind = 8 ) roc_val
  real ( kind = 8 ) x
  real ( kind = 8 ) xold
  real ( kind = 8 ) y

  i_pos = 0
  i_neg = 0

  roc_val = 0.0D+00
  x = 0.0D+00
  y = 0.0D+00
  xold = 0.0D+00

  do i = 1, num_score

    xold = x

    if ( 0 < result(i) ) then

      i_pos = i_pos + 1
      y = real ( i_pos, kind = 8 ) &
        / real ( num_score_pos_select, kind = 8 )

    else if ( result(i) < 0 ) then

      i_neg = i_neg + 1
      x = real ( i_neg, kind = 8 ) &
        / real ( num_score_neg_select, kind = 8 )
      roc_val = roc_val + y * ( x - xold )

    end if

  end do

  return
end
subroutine pns_data_set ( debug, debug_unit, max_score, neg_plot, &
  num_plot, num_score, pn_max_diff, pos_plot, result, s_max_diff, &
  score, score_max_selected, score_min_selected, score_plot )

!*****************************************************************************80
!
!! PNS_DATA_SET sets the data for a positive/negative/score graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical DEBUG, is TRUE if debug output is desired.
!
!    Input, integer ( kind = 4 ) DEBUG_UNIT, the unit number of the debug file.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Output, real ( kind = 8 ) NEG_PLOT(NUM_PLOT), the ordinate of the negative
!    percentage verse score plot.
!
!    Input, integer ( kind = 4 ) NUM_PLOT, the number of points to draw on each curve.
!
!    Input, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Output, real ( kind = 8 ) PN_MAX_DIFF, the maximum difference between the
!    positive and negative curves.
!
!    Output, real ( kind = 8 ) POS_PLOT(NUM_PLOT), the ordinate of the positive
!    percentage verse score plot.
!
!    Input, integer ( kind = 4 ) RESULT(MAX_SCORE), the result code.  RESULT(I) is:
!      J, item I is selected as the J-th POSITIVE representative;
!      0, item I is omitted;
!    - J, item I is selected as the J-th NEGATIVE representative.
!
!    Output, real ( kind = 8 ) S_MAX_DIFF, the score at which the maximum difference
!    between the positive and negative curves occurs.
!
!    Input, real ( kind = 8 ) SCORE(MAX_SCORE), the scores.
!
!    Input, real ( kind = 8 ) SCORE_MAX_SELECTED, SCORE_MIN_SELECTED, the maximum and
!    minimum scores to be plotted.
!
!    Output, real ( kind = 8 ) SCORE_PLOT(NUM_PLOT), the abscissa of both the
!    positive and negative percentage versus score plots.
!
  implicit none

  integer ( kind = 4 ), parameter :: local_max_score = 301
  integer ( kind = 4 ), parameter :: local_num_plot = 101

  integer ( kind = 4 ) max_score
  integer ( kind = 4 ) num_plot

  integer ( kind = 4 ) bin(0:local_num_plot)
  integer ( kind = 4 ) bin_cum(0:local_num_plot+1)
  logical debug
  integer ( kind = 4 ) debug_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin
  integer ( kind = 4 ) neg_index(local_max_score)
  real ( kind = 8 ) neg_plot(num_plot)
  integer ( kind = 4 ) num_score
  integer ( kind = 4 ) num_neg
  integer ( kind = 4 ) num_pos
  real ( kind = 8 ) pn_max_diff
  integer ( kind = 4 ) pos_index(local_max_score)
  real ( kind = 8 ) pos_plot(num_plot)
  integer ( kind = 4 ) result(max_score)
  real ( kind = 8 ) s_max_diff
  real ( kind = 8 ) score(max_score)
  real ( kind = 8 ) score_max_selected
  real ( kind = 8 ) score_min_selected
  real ( kind = 8 ) score_neg(local_max_score)
  real ( kind = 8 ) score_plot(num_plot)
  real ( kind = 8 ) score_pos(local_max_score)

  if ( local_max_score < num_score ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PNS_DATA_SET - Fatal error!'
    write ( *, '(a)' ) '  LOCAL_MAX_SCORE < NUM_SCORE.'
    stop
  end if
!
!  POS_INDEX will be a list of the indices in RESULT which have positive value.
!
  call i4vec_positive_index ( num_score, result, max_score, num_pos, pos_index )

  if ( 0 < num_pos ) then
!
!  Extract the positive scores.
!
    do i = 1, num_pos
      score_pos(i) = score(pos_index(i))
    end do
!
!  Bin the positive scores.
!
    nbin = num_plot - 1

    call r8vec_bin ( num_pos, score_pos, nbin, score_min_selected, &
      score_max_selected, bin, score_plot )

    call i4vec_cum ( num_plot + 1, bin, bin_cum )
!
!  Make the ordinate to be plotted.
!
    do i = 1, num_plot
      pos_plot(i) = real ( bin_cum(i), kind = 8 ) &
                  / real ( num_pos, kind = 8 )
    end do

  end if
!
!  Index negative RESULT values by NEG_INDEX.
!
  call i4vec_negative_index ( num_score, result, max_score, num_neg, neg_index )

  if ( 0 < num_neg ) then
!
!  Extract the negative SCORE's.
!
    do i = 1, num_neg
      score_neg(i) = score(neg_index(i))
    end do
!
!  Bin the negative scores.
!
    nbin = num_plot - 1

    call r8vec_bin ( num_neg, score_neg, nbin, score_min_selected, &
      score_max_selected, bin, score_plot )

    call i4vec_cum ( num_plot + 1, bin, bin_cum )
!
!  Make the ordinate to be plotted.
!
    do i = 1, num_plot
      neg_plot(i) = real ( bin_cum(i), kind = 8 ) &
                  / real ( num_neg, kind = 8 )
    end do

  end if
!
!  Crude estimate of max difference.
!
  s_max_diff = score_plot(1)
  pn_max_diff = neg_plot(1) - pos_plot(1)
  do i = 2, num_plot
    if ( pn_max_diff < neg_plot(i) - pos_plot(i) ) then
      s_max_diff = score_plot(i)
      pn_max_diff = neg_plot(i) - pos_plot(i)
    end if
  end do

  if ( debug ) then
    write ( debug_unit, * ) ' '
    write ( debug_unit, * ) 'PNS_DATA_SET:'
    write ( debug_unit, * ) '  S_MAX_DIFF = ', s_max_diff
    write ( debug_unit, * ) '  PN_MAX_DIFF = ', pn_max_diff
  end if

  if ( debug ) then
    write ( debug_unit, * ) ' '
    write ( debug_unit, * ) 'SCORE_PLOT, NEG_PLOT, POS_PLOT'
    write ( debug_unit, * ) ' '
    do i = 1, num_plot
      write ( debug_unit, * ) score_plot(i), neg_plot(i), pos_plot(i)
    end do
  end if

  return
end
subroutine pns_graph_file_write ( good_file_name, graph_file_unit, ierror, &
  neg_plot, num_good, num_plot, num_score_neg_select, num_score_pos_select, &
  pn_max_diff, pos_plot, s_max_diff, search_file_name, score_plot )

!*****************************************************************************80
!
!! PNS_GRAPH_FILE_WRITE creates a positive/negative/score graphics file of the r
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) GOOD_FILE_NAME, the name of the good ID
!    file to be read.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, real ( kind = 8 ) NEG_PLOT(NUM_PLOT), the ordinate of the negative
!    percentage verse score plot.
!
!    Input, integer ( kind = 4 ) NUM_GOOD, the number of "good" identifiers.
!
!    Input, integer ( kind = 4 ) NUM_PLOT, the number of points to draw on each curve.
!
!    Input, integer ( kind = 4 ) NUM_SCORE_NEG_SELECT, the number of selected scores
!    associated with negative results.
!
!    Input, integer ( kind = 4 ) NUM_SCORE_POS_SELECT, the number of selected scores
!    associated with positive results.
!
!    Output, real ( kind = 8 ) PN_MAX_DIFF, the maximum difference between the
!    positive and negative curves.
!
!    Input, character ( len = * ) PNS_GRAPH_FILE_NAME, the name of the graph
!    file to be created.
!
!    Input, real ( kind = 8 ) POS_PLOT(NUM_PLOT), the ordinate of the positive
!    percentage verse score plot.
!
!    Input, character ( len = * ) SEARCH_FILE_NAME, the name of the search file.
!
!    Input, real ( kind = 8 ) SCORE_PLOT(NUM_PLOT), the abscissa of both the
!    positive and negative percentage versus score plots.
!
  implicit none

  integer ( kind = 4 ) num_plot

  character ( len = * ) good_file_name
  integer ( kind = 4 ) graph_file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) neg_plot(num_plot)
  integer ( kind = 4 ) nlabel
  integer ( kind = 4 ) num_good
  integer ( kind = 4 ) num_score_neg_select
  integer ( kind = 4 ) num_score_pos_select
  real ( kind = 8 ) pn_max_diff
  real ( kind = 8 ) pos_plot(num_plot)
  real ( kind = 8 ) px
  real ( kind = 8 ) pxmax
  real ( kind = 8 ) pxmin
  real ( kind = 8 ) py
  real ( kind = 8 ) pymax
  real ( kind = 8 ) pymin
  real ( kind = 8 ) s
  real ( kind = 8 ) s_max_diff
  character ( len = * ) search_file_name
  real ( kind = 8 ) score_plot(num_plot)
  character ( len = 14 ) string
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yold

  ierror = 0

  x = 0.0D+00
  y = 0.0D+00

  write ( graph_file_unit, '(a)' ) '  page'
  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    line_width 1'
  write ( graph_file_unit, '(a)' ) '    line_rgb 0.5 0.5 0.5'
  write ( graph_file_unit, '(a)' ) '    grid 1.0 2.0 7.5 5.0 11 11'
  write ( graph_file_unit, '(a)' ) '    grid 1.0 6.0 7.5 9.0 11 11'
  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    line_width 2'
  write ( graph_file_unit, '(a)' ) ' '

  xmin = score_plot(1)
  xmax = score_plot(num_plot)
  ymin = 0.0D+00
  ymax = 1.0D+00
!
!  P% versus S plot.
!
  pxmin = 1.0D+00
  pxmax = 7.5D+00
  pymin = 6.0D+00
  pymax = 9.0D+00

  write ( graph_file_unit, '(a)' ) '    line_rgb 0.0 1.0 0.0'

  do i = 2, num_plot

    x = score_plot(i-1)
    y = pos_plot(i-1)
    px = pxmin + ( pxmax - pxmin ) * ( x - xmin ) / ( xmax - xmin )
    py = pymin + ( pymax - pymin ) * ( y - ymin ) / ( ymax - ymin )

    if ( i == 2 ) then
      write ( graph_file_unit, '(a,2g14.6)' ) '    moveto ', px, py
!       else
!         write ( graph_file_unit, '(a,2g14.6)' ) '    moveto ', px, py
    end if

    x = score_plot(i)
    y = pos_plot(i-1)
    px = pxmin + ( pxmax - pxmin ) * ( x - xmin ) / ( xmax - xmin )
    py = pymin + ( pymax - pymin ) * ( y - ymin ) / ( ymax - ymin )

    write ( graph_file_unit, '(a,2g14.6)' ) '    lineto ', px, py

    yold = y

    x = score_plot(i)
    y = pos_plot(i)
    px = pxmin + ( pxmax - pxmin ) * ( x - xmin ) / ( xmax - xmin )
    py = pymin + ( pymax - pymin ) * ( y - ymin ) / ( ymax - ymin )

    if ( yold /= y ) then
      write ( graph_file_unit, '(a,2g14.6)' ) '    lineto ', px, py
    end if

  end do
!
!  N% versus S plot.
!
  write ( graph_file_unit, '(a)' ) '    line_rgb 1.0 0.0 0.0'

  do i = 2, num_plot

    x = score_plot(i-1)
    y = neg_plot(i-1)
    px = pxmin + ( pxmax - pxmin ) * ( x - xmin ) / ( xmax - xmin )
    py = pymin + ( pymax - pymin ) * ( y - ymin ) / ( ymax - ymin )

    if ( i == 2 ) then
      write ( graph_file_unit, '(a,2g14.6)' ) '    moveto ', px, py
    else
      write ( graph_file_unit, '(a,2g14.6)' ) '    moveto ', px, py
    end if

    x = score_plot(i)
    y = neg_plot(i-1)
    px = pxmin + ( pxmax - pxmin ) * ( x - xmin ) / ( xmax - xmin )
    py = pymin + ( pymax - pymin ) * ( y - ymin ) / ( ymax - ymin )

    write ( graph_file_unit, '(a,2g14.6)' ) '    lineto ', px, py

    yold = y

    x = score_plot(i)
    y = neg_plot(i)
    px = pxmin + ( pxmax - pxmin ) * ( x - xmin ) / ( xmax - xmin )
    py = pymin + ( pymax - pymin ) * ( y - ymin ) / ( ymax - ymin )

    if ( yold /= y ) then
      write ( graph_file_unit, '(a,2g14.6)' ) '    lineto ', px, py
    end if

  end do
!
!  N% - P% difference plot.
!
  pxmin = 1.0D+00
  pxmax = 7.5D+00
  pymin = 2.0D+00
  pymax = 5.0D+00

  write ( graph_file_unit, '(a)' ) '    line_rgb 0.0 0.0 1.0'

  do i = 2, num_plot

    x = score_plot(i-1)
    y = neg_plot(i-1) - pos_plot(i-1)

    px = pxmin + ( pxmax - pxmin ) * ( x - xmin ) / ( xmax - xmin )
    py = pymin + ( pymax - pymin ) * ( y - ymin ) / ( ymax - ymin )

    if ( i == 2 ) then
      write ( graph_file_unit, '(a,2g14.6)' ) '    moveto ', px, py
    else
      write ( graph_file_unit, '(a,2g14.6)' ) '    lineto ', px, py
    end if

    x = score_plot(i)
    y = neg_plot(i-1) - pos_plot(i-1)

    px = pxmin + ( pxmax - pxmin ) * ( x - xmin ) / ( xmax - xmin )
    py = pymin + ( pymax - pymin ) * ( y - ymin ) / ( ymax - ymin )

    write ( graph_file_unit, '(a,2g14.6)' ) '    lineto ', px, py

    yold = y

    x = score_plot(i)
    y = neg_plot(i) - pos_plot(i)

    px = pxmin + ( pxmax - pxmin ) * ( x - xmin ) / ( xmax - xmin )
    py = pymin + ( pymax - pymin ) * ( y - ymin ) / ( ymax - ymin )

    if ( y /= yold ) then
      write ( graph_file_unit, '(a,2g14.6)' ) '    lineto ', px, py
    end if

  end do
!
!  Labels.
!
!  Labels at the top.
!
  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    font_size 0.20'
  write ( graph_file_unit, '(a)' ) '    line_rgb 0.0 0.0 0.0'
  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    moveto 2.0 10.0'
  write ( graph_file_unit, '(a)' )  '    label Receiver Operator Characteristic'
  write ( graph_file_unit, '(a)' ) '    moveto 2.0 9.5'
  write ( graph_file_unit, '(a)' ) '    label Percentages less than given Score'
!
!  Labels on the left side of the graphs.
!
  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    moveto 0.65 2.0'
  write ( graph_file_unit, '(a)' ) '    label 0.0'
  write ( graph_file_unit, '(a)' ) '    moveto 0.65 2.6'
  write ( graph_file_unit, '(a)' ) '    label 0.2'
  write ( graph_file_unit, '(a)' ) '    moveto 0.65 3.2'
  write ( graph_file_unit, '(a)' ) '    label 0.4'
  write ( graph_file_unit, '(a)' ) '    moveto 0.65 3.8'
  write ( graph_file_unit, '(a)' ) '    label 0.6'
  write ( graph_file_unit, '(a)' ) '    moveto 0.65 4.4'
  write ( graph_file_unit, '(a)' ) '    label 0.8'
  write ( graph_file_unit, '(a)' ) '    moveto 0.65 5.0'
  write ( graph_file_unit, '(a)' ) '    label 1.0'

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    moveto 0.4 2.5'
  write ( graph_file_unit, '(a)' ) '    label_slant 90.0 Percentage difference'

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    moveto 0.65 6.0'
  write ( graph_file_unit, '(a)' ) '    label 0.0'
  write ( graph_file_unit, '(a)' ) '    moveto 0.65 6.6'
  write ( graph_file_unit, '(a)' ) '    label 0.2'
  write ( graph_file_unit, '(a)' ) '    moveto 0.65 7.2'
  write ( graph_file_unit, '(a)' ) '    label 0.4'
  write ( graph_file_unit, '(a)' ) '    moveto 0.65 7.8'
  write ( graph_file_unit, '(a)' ) '    label 0.6'
  write ( graph_file_unit, '(a)' ) '    moveto 0.65 8.4'
  write ( graph_file_unit, '(a)' ) '    label 0.8'
  write ( graph_file_unit, '(a)' ) '    moveto 0.65 9.0'
  write ( graph_file_unit, '(a)' ) '    label 1.0'

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    moveto 0.4 7.0'
  write ( graph_file_unit, '(a)' ) '    label_slant 90.0 Percentage'
!
!  Labels under the graphs.
!
  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    moveto 4.0 1.5'
  write ( graph_file_unit, '(a)' ) '    label Score'

  nlabel = 5

  do i = 0, nlabel
    s = ( real ( nlabel - i, kind = 8 ) * score_plot(1) &
        + real (          i, kind = 8 ) * score_plot(num_plot) ) &
        / real ( nlabel,     kind = 8 )

    x = ( ( real ( nlabel - i, kind = 8 ) + 0.25 ) * pxmin &
        + ( real (          i, kind = 8 ) - 0.25 ) * pxmax ) &
         / real ( nlabel,      kind = 8 )

    y = 1.75D+00
    write ( graph_file_unit, '(a,2g14.6)' ) '    moveto ', x, y
    call r8_to_s_left ( s, string )
    write ( graph_file_unit, '(a)' ) '    label ' // string
  end do

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) '    moveto 4.0 5.5'
  write ( graph_file_unit, '(a)' ) '    label Score'

  nlabel = 5

  do i = 0, nlabel
    s = (   real ( nlabel - i, kind = 8 ) * score_plot(1) &
          + real (          i, kind = 8 ) * score_plot(num_plot) ) &
          / real ( nlabel,     kind = 8 )

    x = ( ( real ( nlabel - i, kind = 8 ) + 0.25 ) * pxmin &
        + ( real (          i, kind = 8 ) - 0.25 ) * pxmax ) &
          / real ( nlabel,     kind = 8 )

    y = 5.75D+00
    write ( graph_file_unit, '(a,2g14.6)' ) '    moveto ', x, y
    call r8_to_s_left ( s, string )
    write ( graph_file_unit, '(a)' ) '    label ' // string
  end do
!
!  Labels at the bottom.
!
  write ( graph_file_unit, '(a)' ) '    moveto 2.0 1.05'
  call r8_to_s_left ( pn_max_diff, string )
  call s_blanks_delete ( string )
  write ( graph_file_unit, '(a)' ) '    label Max difference N%-P% = ' // &
    trim ( string )

  write ( graph_file_unit, '(a)' ) '    moveto 2.0 0.85'
  call r8_to_s_left ( s_max_diff, string )
  call s_blanks_delete ( string )
  write ( graph_file_unit, '(a)' ) '    label Score at max difference = ' // &
    trim ( string )

  write ( graph_file_unit, '(a)' ) '    moveto 2.0 0.65'
  write ( graph_file_unit, '(a)' ) '    label Score file = ' // &
    trim ( search_file_name )

  write ( graph_file_unit, '(a)' ) '    moveto 2.0 0.45'
  call i4_to_s_left ( num_good, string )
  call s_blanks_delete ( string )
  write ( graph_file_unit, '(a)' ) '    label Good file = ' // &
    trim ( good_file_name ) // ' contains ' // trim ( string ) &
    // ' identifiers.'

  write ( graph_file_unit, '(a)' ) '    moveto 2.0 0.25'
  call i4_to_s_left ( num_score_pos_select, string )
  call s_blanks_delete ( string )
  write ( graph_file_unit, '(a)' ) '    label ' // trim ( string ) &
    // ' positive identifiers were plotted.'

  write ( graph_file_unit, '(a)' ) '    moveto 2.0 0.05'
  call i4_to_s_left ( num_score_neg_select, string )
  call s_blanks_delete ( string )
  write ( graph_file_unit, '(a)' ) '    label ' // trim ( string ) &
    // ' negative identifiers were plotted.'

  write ( graph_file_unit, '(a)' ) ' '
  write ( graph_file_unit, '(a)' ) 'endpage'

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 2000
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
subroutine r8_to_s_left ( d, s )

!*****************************************************************************80
!
!! R8_TO_S_LEFT writes an R8 into a left justified string.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    A 'G14.6' format is used with a WRITE statement.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) D, the number to be written into the string.
!
!    Output, character ( len = * ) S, the string into which
!    the real number is to be written.  If the string is less than 14
!    characters long, it will will be returned as a series of asterisks.
!
  implicit none

  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  character ( len = * ) s
  integer ( kind = 4 ) s_length
  character ( len = 14 ) s2

  s_length = len ( s )

  if ( s_length < 14 ) then

    do i = 1, s_length
      s(i:i) = '*'
    end do

  else if ( d == 0.0D+00 ) then
    s(1:14) = '     0.0      '
  else
    write ( s2, '(g14.6)' ) d
    s(1:14) = s2
  end if
!
!  Shift the string left.
!
  s = adjustl ( s )

  return
end
subroutine r8vec_bin ( n, x, nbin, bin_min, bin_max, bin, bin_limit )

!*****************************************************************************80
!
!! R8VEC_BIN computes bins based on a given real vector.
!
!  Discussion:
!
!    The user specifies minimum and maximum bin values, BIN_MIN and
!    BIN_MAX, and the number of bins, NBIN.  This determines a
!    "bin width":
!
!      H = ( BIN_MAX - BIN_MIN ) / NBIN
!
!    so that bin I will count all entries X(J) such that
!
!      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
!
!    The array X does NOT have to be sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of X.
!
!    Input, real ( kind = 8 ) X(N), an (unsorted) array to be binned.
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.  Two extra bins, #0 and
!    #NBIN+1, count extreme values.
!
!    Input, real ( kind = 8 ) BIN_MIN, BIN_MAX, define the range and size of the bins.
!    BIN_MIN and BIN_MAX must be distinct.
!    Normally, BIN_MIN < BIN_MAX, and the documentation will assume
!    this, but proper results will be computed if BIN_MAX < BIN_MIN.
!
!    Output, integer ( kind = 4 ) BIN(0:NBIN+1).
!    BIN(0) counts entries of X less than BIN_MIN.
!    BIN(NBIN+1) counts entries greater than or equal to BIN_MAX.
!    For 1 <= I <= NBIN, BIN(I) counts the entries X(J) such that
!      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
!    where H is the bin spacing.
!
!    Output, real ( kind = 8 ) BIN_LIMIT(0:NBIN), the "limits" of the bins.
!    BIN(I) counts the number of entries X(J) such that
!      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin

  integer ( kind = 4 ) bin(0:nbin+1)
  real ( kind = 8 ) bin_limit(0:nbin)
  real ( kind = 8 ) bin_max
  real ( kind = 8 ) bin_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)

  if ( bin_max == bin_min ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_BIN - Fatal error!'
    write ( *, '(a)' ) '  BIN_MIN = BIN_MAX.'
    stop
  end if

  bin(0:nbin+1) = 0

  do i = 1, n

    t = ( x(i) - bin_min ) / ( bin_max - bin_min )

    if ( t < 0.0D+00 ) then
      j = 0
    else if ( 1.0D+00 <= t ) then
      j = nbin + 1
    else
      j = 1 + int ( real ( nbin, kind = 8 ) * t )
    end if

    bin(j) = bin(j) + 1

  end do
!
!  Compute the bin limits.
!
  do i = 0, nbin
    bin_limit(i) = ( real ( nbin - i, kind = 8 ) * bin_min &
                   + real (        i, kind = 8 ) * bin_max ) &
                   / real ( nbin,     kind = 8 )
  end do

  return
end
subroutine r8vec_minmax ( n, x, xmin, xmax )

!*****************************************************************************80
!
!! R8VEC_MINMAX returns the minimum and maximum values in a real vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) X(N), the array.
!
!    Output, real ( kind = 8 ) XMIN, XMAX, the smallest and largest values in X.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  if ( n <= 0 ) then

    xmin = 0.0D+00
    xmax = 0.0D+00

  else

    xmin = x(1)
    xmax = x(1)

    do i = 2, n
      xmin = min ( xmin, x(i) )
      xmax = max ( xmax, x(i) )
    end do

  end if

  return
end
subroutine r8vec_order_type ( n, a, order )

!*****************************************************************************80
!
!! R8VEC_ORDER_TYPE determines if an R8VEC is (non)strictly ascending/descending.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the array.
!
!    Input, real ( kind = 8 ) A(N), the array to be checked.
!
!    Output, integer ( kind = 4 ) ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
!
!  Search for the first value not equal to A(1).
!
  i = 1

  do

    i = i + 1

    if ( n < i ) then
      order = 0
      return
    end if

    if ( a(1) < a(i) ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( a(i) < a(1) ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do while ( i < n )

    i = i + 1

    if ( order == 1 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 3
      end if

    end if

  end do

  return
end
subroutine r8vec_permute ( n, a, p )

!*****************************************************************************80
!
!! R8VEC_PERMUTE permutes an R8VEC in place.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    This routine permutes an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!    P(I) = J means that the I-th element of the output array should be
!    the J-th element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,   4,   5,   1,   3 )
!      A = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
!
!    Output:
!
!      A    = ( 2.0, 4.0, 5.0, 1.0, 3.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_temp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = -p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = -p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine r8vec_reverse ( n, a )

!*****************************************************************************80
!
!! R8VEC_REVERSE reverses the elements of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    In FORTRAN90, calling R8VEC_REVERSE is equivalent to
!
!      A(1:N) = A(N:1:-1)
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11.0, 12.0, 13.0, 14.0, 15.0 ).
!
!    Output:
!
!      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n/2
    call r8_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine r8vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I:N)) is sorted,
!
!    or explicitly, by the call
!
!      call r8vec_permute ( n, a, indx )
!
!    after which A(1:N) is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine s_blanks_delete ( s )

!*****************************************************************************80
!
!! S_BLANKS_DELETE replaces consecutive blanks by one blank.
!
!  Discussion:
!
!    Thanks to Bill Richmond for pointing out a programming flaw which
!    meant that, as characters were slid to the left through multiple
!    blanks, their original images were not blanked out.  This problem
!    is easiest resolved by using a copy of the string.
!
!    The remaining characters are left justified and right padded with blanks.
!    TAB characters are converted to spaces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character newchr
  character oldchr
  character ( len = * ) s
  character ( len = len ( s ) ) s_copy
  integer ( kind = 4 ) s_length
  character, parameter :: TAB = achar ( 9 )

  s_length = len ( s )

  j = 0
  s_copy(1:s_length) = s(1:s_length)
  s(1:s_length) = ' '

  newchr = ' '

  do i = 1, s_length

    oldchr = newchr
    newchr = s_copy(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      s(j:j) = newchr
    end if

  end do

  return
end
subroutine s_cap ( s )

!*****************************************************************************80
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character ch
  integer ( kind = 4 ) i
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len_trim ( s )

  do i = 1, s_length

    ch = s(i:i)
    call ch_cap ( ch )
    s(i:i) = ch

  end do

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
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

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  integer ( kind = 4 ) s1_length
  character ( len = * ) s2
  integer ( kind = 4 ) s2_length

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
function s_gei ( s1, s2 )

!*****************************************************************************80
!
!! S_GEI = ( S1 is lexically greater than or equal to S2 ).
!
!  Discussion:
!
!    The comparison is done in a case-insensitive way.
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
!    Output, logical S_GEI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lenc
  logical s_gei
  character ( len = * ) s1
  integer ( kind = 4 ) s1_length
  character ( len = * ) s2
  integer ( kind = 4 ) s2_length

  s1_length = len_trim ( s1 )
  s2_length = len_trim ( s2 )
  lenc = min ( s1_length, s2_length )

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( lgt ( c1, c2 ) ) then
      s_gei = .true.
      return
    else if ( llt ( c1, c2 ) ) then
      s_gei = .false.
      return
    end if

  end do

  if ( s1_length < s2_length ) then
    s_gei = .false.
  else
    s_gei = .true.
  end if

  return
end
function s_gti ( s1, s2 )

!*****************************************************************************80
!
!! S_GTI = S1 is lexically greater than S2.
!
!  Discussion:
!
!    The comparison is done in a case-insensitive way.
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
!    Output, logical S_GTI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lenc
  logical s_gti
  character ( len = * ) s1
  integer ( kind = 4 ) s1_length
  character ( len = * ) s2
  integer ( kind = 4 ) s2_length

  s1_length = len ( s1 )
  s2_length = len ( s2 )
  lenc = min ( s1_length, s2_length )

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( lgt ( c1, c2 ) ) then
      s_gti = .true.
      return
    else if ( llt ( s1, s2 ) ) then
      s_gti = .false.
      return
    end if

  end do

  if ( s1_length <= s2_length ) then
    s_gti = .false.
  else
    s_gti = .true.
  end if

  return
end
function s_lti ( s1, s2 )

!*****************************************************************************80
!
!! S_LTI = ( S1 is lexically less than S2 ).
!
!  Discussion:
!
!    The comparison is done in a case-insensitive way.
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
!    Output, logical S_LTI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_lti
  character ( len = * ) s1
  integer ( kind = 4 ) s1_length
  character ( len = * ) s2

  s1_length = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( s1_length, len2 )

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( llt ( c1, c2 ) ) then
      s_lti = .true.
      return
    else if ( lgt ( c1, c2 ) ) then
      s_lti = .false.
      return
    end if

  end do

  if ( s1_length < len2 ) then
    s_lti = .true.
  else
    s_lti = .false.
  end if

  return
end
subroutine s_swap ( s1, s2 )

!*****************************************************************************80
!
!! S_SWAP swaps two strings.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S1, S2.  On output, the values of S1
!    and S2 have been interchanged.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = 256 ) s3

  s3 = s1
  s1 = s2
  s2 = s3

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
!    13 January 2006
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
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S used to make
!    the integer.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * ) s
  integer ( kind = 4 ) state
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

      if ( c == ' ' ) then

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

      if ( c == ' ' ) then

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
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 value from a string.
!
!  Discussion:
!
!    An "R8" value is simply a real number to be stored as a
!    variable of type "real ( kind = 8 )".
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( s_length < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
  if ( iterm /= 1 .and. length+1 == s_length ) then
    length = s_length
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

  return
end
subroutine s_word_find ( s, iword, word, nchar )

!*****************************************************************************80
!
!! S_WORD_FIND finds the word of a given index in a string.
!
!  Discussion:
!
!    A "word" is any string of nonblank characters, separated from other
!    words by one or more blanks or TABS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, integer ( kind = 4 ) IWORD, the index of the word to be
!    searched for.  If IWORD is positive, then the IWORD-th
!    word is sought.  If IWORD is zero or negative, then
!    assuming that the string has N words in it, the
!    N+IWORD-th word will be sought.
!
!    Output, character ( len = * ) WORD, the IWORD-th word of the
!    string, or ' ' if the WORD could not be found.
!
!    Output, integer ( kind = 4 ) NCHAR, the number of characters in WORD,
!    or 0 if the word could not be found.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iblank
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) iword
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jword
  integer ( kind = 4 ) kword
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  integer ( kind = 4 ) s_len
  character, parameter :: TAB = achar ( 9 )
  character ( len = * ) word

  ilo = 0
  ihi = 0
  s_len = len_trim ( s )

  if ( s_len <= 0 ) then
    return
  end if

  if ( 0 < iword ) then

    if ( s(1:1) == ' ' .or. s(1:1) == TAB ) then
      iblank = 1
      jword = 0
      jlo = 0
      jhi = 0
    else
      iblank = 0
      jword = 1
      jlo = 1
      jhi = 1
    end if

    i = 1

    do

      i = i + 1

      if ( s_len < i ) then

        if ( jword == iword ) then
          ilo = jlo
          ihi = s_len
          nchar = s_len + 1 - jlo
          word = s(ilo:ihi)
        else
          ilo = 0
          ihi = 0
          nchar = 0
          word = ' '
        end if

        return

      end if

      if ( ( s(i:i) == ' ' .or. s(i:i) == TAB ) .and. iblank == 0 ) then

        jhi = i - 1
        iblank = 1
        if ( jword == iword ) then
          ilo = jlo
          ihi = jhi
          nchar = jhi + 1 - jlo
          word = s(ilo:ihi)
          return
        end if

      else if ( s(i:i) /= ' ' .and. s(i:i) /= TAB .and. iblank == 1 ) then

        jlo = i
        jword = jword + 1
        iblank = 0

      end if

    end do

  else

    iblank = 0
    kword = 1 - iword
    jword = 1
    jlo = s_len
    jhi = s_len
    i = s_len

    do

      i = i - 1

      if ( i <= 0 ) then

        if ( jword == kword ) then
          ilo = 1
          ihi = jhi
          nchar = jhi
          word = s(ilo:ihi)
        else
          ilo = 0
          ihi = 0
          nchar = 0
          word = ' '
        end if

        return

      end if

      if ( ( s(i:i) == ' ' .or. s == TAB ) .and. iblank == 0 ) then

        jlo = i + 1
        iblank = 1

        if ( jword == kword ) then
          ilo = jlo
          ihi = jhi
          nchar = jhi + 1 - jlo
          word = s(ilo:ihi)
          return
        end if

      else if ( s(i:i) /= ' ' .and. s(i:i) /= TAB .and. iblank == 1 ) then

        jhi = i
        jword = jword + 1
        iblank = 0

      end if

    end do

  end if

  return
end
subroutine score_blast_file_read ( identifier, ierror, max_score, num_score, &
  score, search_file_name )

!*****************************************************************************80
!
!! SCORE_BLAST_FILE_READ extracts information from a BLAST search file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 35 ) IDENTIFIER(MAX_SCORE), a list of the
!    (unique) identifiers for which a score was computed.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Output, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Output, real ( kind = 8 ) SCORE(MAX_SCORE), the scores.
!
!    Input, character ( len = * ) SEARCH_FILE_NAME, the name of the search file.
!
  implicit none

  integer ( kind = 4 ) max_score

  integer ( kind = 4 ) i_temp
  integer ( kind = 4 ) i1_temp
  integer ( kind = 4 ) i2_temp
  character ( len = 35 ) identifier(max_score)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) lchar
  character ( len = 100 ) line
  integer ( kind = 4 ) num_score
  real ( kind = 8 ) score(max_score)
  character ( len = * ) search_file_name
  integer ( kind = 4 ) score_unit
  character ( len = 20 ) string

  ierror = 0
  num_score = 0

  call get_unit ( score_unit )

  open ( unit = score_unit, file = search_file_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SCORE_BLAST_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the BLAST file:'
    write ( *, '(a)' ) trim ( search_file_name )
    return
  end if
!
!  Advance through the BLAST file to the line "Sequences producing",
!  read the next (almost blank) line, and then expect to read line of data.
!
  string = 'Sequences producing'

  call file_advance_to_string ( score_unit, string, line, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SCORE_BLAST_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  End-of-file, no marker string.'
    close ( unit = score_unit )
    return
  end if

  read ( score_unit, '(a)', iostat = ios ) line

  if ( ios == 0 ) then
!
!  The data lines have the format:
!
!    USELESS:Identifier xxxxxxxx score xxxxx
!
    do

      read ( score_unit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      if ( line == ' ' ) then
        exit
      end if

      if ( line(3:13) == 'End of List' ) then
        exit
      end if

      i1_temp = index ( line, ':' )
      i2_temp = index ( line, '!' )

      num_score = num_score + 1
      if ( num_score <= max_score ) then
        identifier(num_score) = line(i1_temp+1:i2_temp-2)
      end if

      call s_to_i4 ( line(63:66), i_temp, ierror, lchar )

      if ( ierror /= 0 ) then
        ierror = 3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SCORE_BLAST_FILE_READ - Fatal error!'
        write ( *, '(a)' ) '  Bad score data.'
        exit
      end if

      if ( num_score <= max_score ) then
        score(num_score) = real ( i_temp, kind = 8 )
      end if

    end do

  end if

  close ( unit = score_unit )

  if ( max_score < num_score ) then
    num_score = max_score
  end if

  return
end
subroutine score_fasta_file_read ( identifier, ierror, max_score, num_score, &
  score, search_file_name )

!*****************************************************************************80
!
!! SCORE_FASTA_FILE_READ extracts information from a FASTA search file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 35 ) IDENTIFIER(MAX_SCORE), a list of the
!    (unique) identifiers for which a score was computed.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Output, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Output, real ( kind = 8 ) SCORE(MAX_SCORE), the scores.
!
!    Input, character ( len = * ) SEARCH_FILE_NAME, the name of the search file.
!
  integer ( kind = 4 ) max_score

  integer ( kind = 4 ) i_temp
  character ( len = 35 ) identifier(max_score)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) lchar
  character ( len = 100 ) line
  integer ( kind = 4 ) num_score
  real ( kind = 8 ) score(max_score)
  character ( len = * ) search_file_name
  integer ( kind = 4 ) score_unit
  character ( len = 20 ) string

  ierror = 0
  num_score = 0

  call get_unit ( score_unit )

  open ( unit = score_unit, file = search_file_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SCORE_FASTA_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the FASTA file:'
    write ( *, '(a)' ) trim ( search_file_name )
    return
  end if
!
!  Advance through the FASTA file to the line "The best scores are:",
!  read the next (blank) line, and then expect to read pairs of lines.
!
  string = 'The best scores are:'

  call file_advance_to_string ( score_unit, string, line, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SCORE_FASTA_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  End-of-file, no marker string.'
    close ( unit = score_unit )
    return
  end if

  read ( score_unit, '(a)', iostat = ios ) line

  if ( ios == 0 ) then
!
!  The pairs of lines have the format:
!
!    USELESS:Identifier
!    xxxxxxxx score xxxxx
!
    do

      read ( score_unit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      if ( line == ' ' ) then
        exit
      end if

      if ( line(3:13) == 'End of List' ) then
        exit
      end if

      i_temp = index ( line, ':' )
      num_score = num_score + 1

      if ( num_score <= max_score ) then
        identifier(num_score) = line(i_temp+1:)
      end if

      read ( score_unit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      call s_to_i4 ( line(54:57), i_temp, ierror, lchar )

      if ( ierror /= 0 ) then
        ierror = 3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SCORE_FASTA_FILE_READ - Fatal error!'
        write ( *, '(a)' ) '  Bad score data.'
        exit
      end if

      if ( num_score <= max_score ) then
        score(num_score) = real ( i_temp, kind = 8 )
      end if

    end do

  end if

  close ( unit = score_unit )

  if ( max_score < num_score ) then
    num_score = max_score
  end if

  return
end
subroutine score_generic_file_read ( identifier, ierror, max_score, &
  num_score, score, search_file_name )

!*****************************************************************************80
!
!! SCORE_GENERIC_FILE_READ reads a "name, score" file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 35 ) IDENTIFIER(MAX_SCORE), a list of the
!    (unique) identifiers for which a score was computed.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Output, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Output, real ( kind = 8 ) SCORE(MAX_SCORE), the scores.
!
!    Input, character ( len = * ) SEARCH_FILE_NAME, the name of the search file.
!
  implicit none

  integer ( kind = 4 ) max_score

  integer ( kind = 4 ), parameter :: id_column = 1
  character ( len = 35 ) id_temp
  character ( len = 35 ) identifier(max_score)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 80 ) line
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) num_score
  real ( kind = 8 ) score(max_score)
  integer ( kind = 4 ), parameter :: score_column = 2
  character ( len = * ) search_file_name
  integer ( kind = 4 ) score_unit
  real ( kind = 8 ) value
  character ( len = 80 ) word

  ierror = 0
  num_score = 0

  call get_unit ( score_unit )

  open ( unit = score_unit, file = search_file_name, status = 'old', &
    access = 'sequential', form = 'formatted', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SCORE_GENERIC_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) trim ( search_file_name )
    return
  end if

  do

    read ( score_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Comment lines begin with '#' and are ignored.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if
!
!  Blank lines are ignored.
!
    if ( line == ' ' ) then
      cycle
    end if
!
!  Read the identifier:
!
    call s_word_find ( line, id_column, id_temp, nchar )

    if ( nchar <= 0 ) then
      ierror = 3
      exit
    end if
!
!  Read the score.
!
    call s_word_find ( line, score_column, word, nchar )

    if ( nchar <= 0 ) then
      ierror = 3
      exit
    end if

    call s_to_r8 ( word, value, ierror, nchar )

    if ( ierror /= 0 ) then
      ierror = 4
      exit
    end if
!
!  The line is "readable".
!
    num_score = num_score + 1

    if ( num_score <= max_score ) then
      identifier(num_score) = id_temp
      score(num_score) = value
    end if

  end do

  close ( unit = score_unit )

  return
end
subroutine score_maxsegs_file_read ( identifier, ierror, max_score, &
  num_score, score, search_file_name )

!*****************************************************************************80
!
!! SCORE_MAXSEGS_FILE_READ extracts information from a MAXSEGS search file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 35 ) IDENTIFIER(MAX_SCORE), a list of the
!    (unique) identifiers for which a score was computed.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Output, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Output, real ( kind = 8 ) SCORE(MAX_SCORE), the scores.
!
!    Input, character ( len = * ) SEARCH_FILE_NAME, the name of the search file.
!
  implicit none

  integer ( kind = 4 ) max_score

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_temp
  character ( len = 35 ) identifier(max_score)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) lchar
  character ( len = 100 ) line
  integer ( kind = 4 ) num_score
  real ( kind = 8 ) score(max_score)
  character ( len = * ) search_file_name
  integer ( kind = 4 ) score_unit
  character ( len = 20 ) string

  ierror = 0
  num_score = 0

  call get_unit ( score_unit )

  open ( unit = score_unit, file = search_file_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SCORE_MAXSEGS_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the MAXSEGS file:'
    write ( *, '(a)' ) trim ( search_file_name )
    return
  end if
!
!  The data lines have the format:
!
!     *****   Library=IDENTIFIER xxxxx scored SCORE.
!
  do

    string = 'Library='

    call file_advance_to_string ( score_unit, string, line, ierror )

    if ( ierror /= 0 ) then
      ierror = 0
      exit
    end if
!
!  Extract the identifier and score from the line.
!
    ilo = 19
    ihi = ilo

    do while ( line(ihi+1:ihi+1) /= ' ' )
      ihi = ihi + 1
    end do

    num_score = num_score + 1
    if ( num_score <= max_score ) then
      identifier(num_score) = line(ilo:ihi)
    end if

    i = index ( line, ' 1 scored' )
    call s_to_i4 ( line(i+9:), i_temp, ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = 3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCORE_MAXSEGS_FILE_READ - Fatal error!'
      write ( *, '(a)' ) '  Bad score data.'
      exit
    end if

    if ( num_score <= max_score ) then
      score(num_score) = real ( i_temp, kind = 8 )
    end if

  end do

  close ( unit = score_unit )

  if ( max_score < num_score ) then
    num_score = max_score
  end if

  return
end
subroutine score_pearson_file_read ( identifier, ierror, max_score, &
  num_score, score, search_file_name )

!*****************************************************************************80
!
!! SCORE_PEARSON_FILE_READ extracts information from a PEARSON search file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 35 ) IDENTIFIER(MAX_SCORE), a list of the
!    (unique) identifiers for which a score was computed.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer ( kind = 4 ) MAX_SCORE, the maximum number of scores.
!
!    Output, integer ( kind = 4 ) NUM_SCORE, the number of scores.
!
!    Output, real ( kind = 8 ) SCORE(MAX_SCORE), the scores.
!
!    Input, character ( len = * ) SEARCH_FILE_NAME, the name of the search file.
!
  implicit none

  integer ( kind = 4 ) max_score

  integer ( kind = 4 ) i_temp
  character ( len = 35 ) identifier(max_score)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) lchar
  character ( len = 100 ) line
  integer ( kind = 4 ) num_score
  real ( kind = 8 ) score(max_score)
  character ( len = * ) search_file_name
  integer ( kind = 4 ) score_unit
  character ( len = 20 ) string

  ierror = 0
  num_score = 0

  call get_unit ( score_unit )

  open ( unit = score_unit, file = search_file_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SCORE_PEARSON_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the PEARSON file:'
    write ( *, '(a)' ) trim ( search_file_name )
    return
  end if
!
!  Advance through the PEARSON file to the line "The best scores are:",
!  and then expect to read line of data.
!
  string = 'The best scores are:'

  call file_advance_to_string ( score_unit, string, line, ierror )

  if ( ierror /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SCORE_PEARSON_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  End-of-file, no marker string.'
    close ( unit = score_unit )
    return
  end if

  read ( score_unit, '(a)', iostat = ios ) line

  if ( ios == 0 ) then
!
!  The data lines have the format:
!
!    Identifier xxxxxxxx score xxxxx
!
    do

      read ( score_unit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      if ( line == ' ' ) then
        exit
      end if

      i_temp = index ( line, ' ' )

      num_score = num_score + 1
      if ( num_score <= max_score ) then
        identifier(num_score) = line(1:i_temp-1)
      end if

      call s_to_i4 ( line(60:63), i_temp, ierror, lchar )

      if ( ierror /= 0 ) then
        ierror = 3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SCORE_PEARSON_FILE_READ - Fatal error!'
        write ( *, '(a)' ) '  Bad score data.'
        exit
      end if

      if ( num_score <= max_score ) then
        score(num_score) = real ( i_temp, kind = 8 )
      end if

    end do

  end if

  close ( unit = score_unit )

  if ( max_score < num_score ) then
    num_score = max_score
  end if

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    Original FORTRAN77 version by A Nijenhuis, H Wilf.
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
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
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
!      * set ISGN = -1 if I precedes J, ISGN = +1 if J precedes I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I precedes J;
!    ISGN => 0 means J precedes I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0

  if ( n <= 1 ) then
    indx = 0
    return
  end if
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    n1 = n
    k = n / 2
    k1 = k
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i = i + 1
      end if

      j = k1
      k1 = i
      indx = - 1
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        indx = 0
      else
        i = n1
        n1 = n1 - 1
        j = 1
        indx = 1
      end if

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

    i = 2 * k1

    if ( i == n1 ) then
      j = k1
      k1 = i
      indx = - 1
      return
    else if ( i <= n1 ) then
      j = i + 1
      indx = - 2
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    indx = 0
  else
    i = n1
    n1 = n1 - 1
    j = 1
    indx = 1
  end if

  return
end
subroutine svec_permute ( n, a, p )

!*****************************************************************************80
!
!! SVEC_PERMUTE permutes a string vector in place.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (  3,     2,     4,       2,      1 )
!      A = ( 'ONE', 'TWO', 'THREE', 'FOUR', 'FIVE' )
!
!    Output:
!
!      A    = ( 'FIVE', 'FOUR', 'ONE', 'THREE', 'TWO' ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, character ( len = * ) A(N), the array to be permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = * ) a(n)
  character ( len = 256 ) a_temp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVEC_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i8)' ) '  array is missing the value ', ierror
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = -p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = -p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SVEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine svec_reverse ( n, a )

!*****************************************************************************80
!
!! SVEC_REVERSE reverses the elements of a string vector.
!
!  Example:
!
!    Input:
!
!      N = 4,
!      A = ( 'Bob', 'Carol', 'Ted', 'Alice' ).
!
!    Output:
!
!      A = ( 'Alice', 'Ted', 'Carol', 'Bob' ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, character ( len = * ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = * ) a(n)
  character ( len = 256 ) a_temp
  integer ( kind = 4 ) i

  do i = 1, n/2
    a_temp = a(i)
    a(i) = a(n+1-i)
    a(n+1-i) = a_temp
  end do

  return
end
subroutine sveci_search_binary_a ( n, a, b, indx )

!*****************************************************************************80
!
!! SVECI_SEARCH_BINARY_A searches an ascending sorted vector of implicitly capitalized strings.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.9,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the vector.
!
!    Input, character ( len = * ) A(N), the array to be searched.  A must
!    be sorted in increasing order.
!
!    Input, character ( len = * ) B, the value to be searched for.
!
!    Output, integer ( kind = 4 ) INDX, the result of the search.
!    0, B does not occur in A.
!    I, A(I) = B, ignoring capitalization.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = * ) a(n)
  character ( len = * ) b
  integer ( kind = 4 ) high
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid
  logical s_eqi
  logical s_gti
  logical s_lti

  indx = 0

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( s_eqi ( a(mid), b ) ) then
      indx = mid
      exit
    else if ( s_lti ( a(mid), b ) ) then
      low = mid + 1
    else if ( s_gti ( a(mid), b ) ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine sveci_sort_heap_a ( n, sarray )

!*****************************************************************************80
!
!! SVECI_SORT_HEAP_A heap sorts a vector of implicitly capitalized strings.
!
!  Discussion:
!
!    The characters in an implicitly capitalized string are treated as
!    though they had been capitalized.  Thus, the letters 'a' and 'A'
!    are considered equal, both 'a' and 'A' precede 'B', and
!    'Fox' and 'fOx' are considered equal.
!
!    The ASCII collating sequence is used, except that all
!    alphabetic characters are treated as though they were uppercase.
!
!    This means
!
!      A = a < B = b < C = c < .... < Y = y < Z = z.
!
!    Numbers and other symbols may also occur, and will be sorted
!    according to the ASCII ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in SARRAY.
!
!    Input/output, character ( len = * ) SARRAY(N), the array to be sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: MAX_CHAR = 255
  integer ( kind = 4 ) n

  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  logical s_gei
  logical s_lti
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1
  character ( len = * ) sarray(n)
  character ( len = MAX_CHAR ) s

  n1 = n
  l = n / 2
  s = sarray(l)
  l1 = l

  do

    m = 2 * l1

    if ( m <= n1 ) then

      if ( m < n1 ) then
        if ( s_gei ( sarray(m+1), sarray(m) ) ) then
          m = m + 1
        end if
      end if

      if ( s_lti ( s, sarray(m) ) ) then
        sarray(l1) = sarray(m)
        l1 = m
        cycle
      end if

    end if

    sarray(l1) = s

    if ( 1 < l ) then
      l = l - 1
      s = sarray(l)
      l1 = l
      cycle
    end if

    if ( n1 < 2 ) then
      exit
    end if

    s = sarray(n1)
    sarray(n1) = sarray(1)

    n1 = n1 - 1
    l1 = l

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
