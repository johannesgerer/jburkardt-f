program main

!*****************************************************************************80
!
!! MAIN is the main program for FOOTBALL.
!
!  Discussion:
!
!    FOOTBALL ranks college football teams.
!
!    Only some of the options have been implemented.
!
!    I could use a SMALL database for testing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Keener,
!    The Perron-Frobenius Theorem and the Ranking of Football Teams,
!    SIAM Review, Volume 35, Number 1, pages 80-93, March 1993.
!
!    Donald Knuth,
!    The Stanford Graph Base,
!    ACM Press, 1993.
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: game_max = 700
  integer ( kind = 4 ), parameter :: team_num = 120

  integer ( kind = 4 ) game_num
  integer ( kind = 4 ) game_score1(game_max)
  integer ( kind = 4 ) game_score2(game_max)
  integer ( kind = 4 ) game_team1(game_max)
  integer ( kind = 4 ) game_team2(game_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) rank(team_num)
  real ( kind = 8 ) rank_factor
  real ( kind = 8 ) rank_vector(team_num)
  integer ( kind = 4 ) rank_option
  real ( kind = 8 ), dimension ( team_num, team_num ) :: score_matrix
  integer ( kind = 4 ) :: score_option = 1
  integer ( kind = 4 ) team
  character ( len = 5 ) team_abbrev(team_num)
  character ( len = 30 ) team_name(team_num)
  integer ( kind = 4 ) team_games(team_num)
  integer ( kind = 4 ) team_loss(team_num)
  integer ( kind = 4 ) team_tie(team_num)
  integer ( kind = 4 ) team_won(team_num)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FOOTBALL'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read a set of football scores and rank the teams.'
!
!  Read the game score data file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = 'football.dat', status = 'old' )

  call games_read ( game_max, game_num, game_score1, game_score2, game_team1, &
    game_team2, iunit, team_abbrev, team_name, team_num )
  close ( unit = iunit )
!
!  Compute win/loss/tie
!
  team_games(1:team_num) = 0
  team_loss(1:team_num) = 0
  team_tie(1:team_num) = 0
  team_won(1:team_num) = 0

  do i = 1, game_num
    i1 = game_team1(i)
    i2 = game_team2(i)
    team_games(i1) = team_games(i1) + 1
    team_games(i2) = team_games(i2) + 1
    if ( game_score2(i) < game_score1(i) ) then
      team_won(i1) = team_won(i1) + 1
      team_loss(i2) = team_loss(i2) + 1
    else if ( game_score1(i) < game_score2(i) ) then
      team_loss(i1) = team_loss(i1) + 1
      team_won(i2) = team_won(i2) + 1
    else if ( game_score1(i) == game_score2(i) ) then
      team_tie(i1) = team_tie(i1) + 1
      team_tie(i2) = team_tie(i2) + 1
    end if
  end do
!
!  Print information from the file.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) 'Number of scores in the data file is ', game_num

  game_num = min ( game_num, game_max )

  call team_print ( team_abbrev, team_name, team_num )

  call game_print ( game_num, game_score1, game_score2, game_team1, &
    game_team2, team_abbrev, team_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Choose the ranking option:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  1: compute a scoring matrix, find its eigenvector;'
  write ( *, '(a)' ) '  2: find a fixed point of a nonlinear mapping.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter option, 1 or 2:'
  read ( *, * ) rank_option
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Using RANK_OPTION = ', rank_option

  if ( rank_option == 1 ) then
!
!  Compute the scoring matrix.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Prepare the scoring matrix.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  If team I scores Sij and team J scores Sji,'
    write ( *, '(a)' ) '    then the formula for the score matrix entry S(I,J)'
    write ( *, '(a)' ) '    can be:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  1:  0, 1/2, or 1 for loss, tie or win;'
    write ( *, '(a)' ) '  2:  Sij;'
    write ( *, '(a)' ) '  3:  Sij / (Sij+Sji);'
    write ( *, '(a)' ) '  4:  (Sij+1) / (Sij+Sji+2);'
    write ( *, '(a)' ) '  5:  H(Sij,Sji).'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the scoring option you want, from 1 to 5:'

    read ( *, * ) score_option

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Using SCORE_OPTION = ', score_option

    call score_matrix_set ( game_num, game_score1, game_score2, game_team1, &
      game_team2, score_option, team_num, score_matrix )
!
!  Determine the rankings.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Determine the rankings.'

    call r8mat_power_method ( team_num, team_num, score_matrix, &
      rank_factor, rank_vector )
!
!  METHOD 2
!
  else if ( rank_option == 2 ) then

    rank_vector(1:team_num) = 1.0D+00

    call score_matrix_set2 ( game_num, game_score1, game_score2, game_team1, &
      game_team2, team_num, score_matrix )

    call picard ( team_num, score_matrix, rank_vector, team_games )

  end if
!
!  Compute RANK from RANK_VECTOR.
!
  call r8vec_sort_insert_index_d ( team_num, rank_vector, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  #  Team  Games/Won/Loss/Tie  Strength'
  write ( *, '(a)' ) ' '

  do i = 1, team_num
    team = rank(i)
    write ( *, '(i3,2x,a5,2x,4i3,f10.4)' ) i, team_abbrev(team), &
      team_games(team), team_won(team), team_loss(team), team_tie(team), &
      rank_vector(team)
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FOOTBALL'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
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
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine ch_extract ( string, c )

!*****************************************************************************80
!
!! CH_EXTRACT extracts the next nonblank character from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) STRING, the string.  On output, the
!    first nonblank character of STRING has been removed, and STRING
!    has been shifted left.
!
!    Output, character C, the leading character of STRING.
!
  implicit none

  character c
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) lchar
  character ( len = * ) string

  c = ' '

  lchar = len_trim ( string )
!
!  Find the first nonblank.
!
  iget = 0

  do

    iget = iget + 1

    if ( lchar < iget ) then
      return
    end if

    if ( string(iget:iget) /= ' ' ) then
      exit
    end if

  end do
!
!  Copy the nonblank character.
!
  c = string(iget:iget)
!
!  Shift the string.
!
  string(1:iget) = ' '
  string = adjustl ( string(iget+1:) )

  return
end
subroutine game_print ( game_num, game_score1, game_score2, game_team1, &
  game_team2, team_abbrev, team_num )

!*****************************************************************************80
!
!! GAME_PRINT prints the game information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GAME_NUM, the number of games.
!
!    Input, integer ( kind = 4 ) GAME_SCORE1(GAME_NUM), GAME_SCORE2(GAME_NUM),
!    the scores for each game.
!
!    Input, integer ( kind = 4 ) GAME_TEAM1(GAME_NUM), GAME_TEAM2(GAME_NUM),
!    the teams for each game.
!
!    Input, character ( len = 5 ) TEAM_ABBREV(TEAM_NUM), the abbreviations for
!    each team.
!
!    Input, integer ( kind = 4 ) TEAM_NUM, the number of teams.
!
  implicit none

  integer ( kind = 4 ) game_num
  integer ( kind = 4 ) team_num

  integer ( kind = 4 ) game_i
  integer ( kind = 4 ) game_score1(game_num)
  integer ( kind = 4 ) game_score2(game_num)
  integer ( kind = 4 ) game_team1(game_num)
  integer ( kind = 4 ) game_team2(game_num)
  character ( len = 5 ) team_abbrev(team_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Game scores:'
  write ( *, '(a)' ) ' '
  do game_i = 1, game_num

    write ( *, '(i3,2x,a,2x,i3,4x,a,2x,i3)' ) game_i, &
      team_abbrev(game_team1(game_i)), game_score1(game_i), &
      team_abbrev(game_team2(game_i)), game_score2(game_i)

  end do

  return
end
subroutine games_read ( game_max, game_num, game_score1, game_score2, &
  game_team1, game_team2, iunit, team_abbrev, team_name, team_num )

!*****************************************************************************80
!
!! GAMES_READ reads the game score data file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GAME_MAX, the maximum number of games.
!
!    Output, integer ( kind = 4 ) GAME_NUM, the number of games.
!
!    Output, integer ( kind = 4 ) GAME_SCORE1(GAME_NUM), GAME_SCORE2(GAME_NUM),
!    the scores for each game.
!
!    Output, integer ( kind = 4 ) GAME_TEAM1(GAME_NUM), GAME_TEAM2(GAME_NUM),
!    the teams for each game.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit number of the file being read.
!
!    Input, character ( len = 5 ) TEAM_ABBREV(TEAM_NUM), the abbreviations for
!    each team.
!
!    Output, character ( len = 30 ) TEAM_NAME(TEAM_NUM), the full name of
!    each team.
!
!    Input, integer ( kind = 4 ) TEAM_NUM, the number of teams.
!
  implicit none

  integer ( kind = 4 ) game_max
  integer ( kind = 4 ) team_num

  character c
  integer ( kind = 4 ) game_num
  integer ( kind = 4 ) game_score1(game_max)
  integer ( kind = 4 ) game_score2(game_max)
  integer ( kind = 4 ) game_team1(game_max)
  integer ( kind = 4 ) game_team2(game_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ival1
  integer ( kind = 4 ) ival2
  character ( len = 100 ) line
  integer ( kind = 4 ) match1
  integer ( kind = 4 ) match2
  character ( len = 5 ) team_abbrev(team_num)
  character ( len = 30 ) team_name(team_num)
!
!  Phase 0: skip header
!
  do i = 1, 4
    read ( iunit, '(a)' ) line
  end do
!
!  Phase 1: Read school abbreviations and names
!
  call team_read ( iunit, team_abbrev, team_name, team_num )
!
!  Phase 2: Read scores
!
  game_num = 0

  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Phase 2.1: Is this a comment?
!
    if ( line(1:1) == '*' ) then
      cycle
    end if
!
!  Phase 2.2: Is this a date?
!
    if ( line(1:1) == '>' ) then
      cycle
    end if
!
!  Phase 2.3: This is a game.
!
    call token_extract ( line, team_num, team_abbrev, match1 )
    call i4_extract ( line, ival1, ierror )
    call ch_extract ( line, c )
    call token_extract ( line, team_num, team_abbrev, match2 )
    call i4_extract ( line, ival2, ierror )

    game_num = game_num + 1

    if ( game_num <= game_max ) then
      game_team1(game_num) = match1
      game_score1(game_num) = ival1
      game_team2(game_num) = match2
      game_score2(game_num) = ival2
    end if

  end do

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer ( kind = 4 ) between 1 and 99 which
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
!    Otherwise, IUNIT is an integer ( kind = 4 ) between 1 and 99, representing a
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
subroutine i4_extract ( string, ival, ierror )

!*****************************************************************************80
!
!! I4_EXTRACT "extracts" an I4 from the beginning of a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) STRING; on input, a string from
!    whose beginning an integer ( kind = 4 ) is to be extracted.  On output,
!    the integer ( kind = 4 ), if found, has been removed.
!
!    Output, integer ( kind = 4 ) IVAL.  If IERROR is 0, then IVAL contains the
!    "next" integer ( kind = 4 ) read from LINE; otherwise IVAL is 0.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, no error.
!    nonzero, an integer ( kind = 4 ) could not be extracted from the beginning of the
!    string.  IVAL is 0 and STRING is unchanged.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) lchar
  character ( len = * ) string

  ival = 0

  call s_to_i4 ( string, ival, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    ierror = 1
    ival = 0
    return
  end if

  string = adjustl ( string(lchar+1:) )

  return
end
subroutine i4vec_identity ( n, a )

!*****************************************************************************80
!
!! I4VEC_IDENTITY sets an I4VEC to the identity vector A(I)=I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
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
subroutine picard ( team_num, score_matrix, rank_vector, team_games )

!*****************************************************************************80
!
!! PICARD uses Picard iteration to find an approximate nonlinear ranking.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) team_num

  real ( kind = 8 ) dx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iteration
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: picard_max = 20
  real ( kind = 8 ), parameter :: picard_tolerance = 0.001D+00
  real ( kind = 8 ) rank_vector(team_num)
  real ( kind = 8 ) rank2(team_num)
  real ( kind = 8 ) score_matrix(team_num,team_num)
  integer ( kind = 4 ) team_games(team_num)
  real ( kind = 8 ) warp2
  real ( kind = 8 ) xnrm

  iteration = 0

  do

    do i = 1, team_num
      rank2(i) = 0.0D+00
      do j = 1, team_num
        rank2(i) = rank2(i) + warp2 ( score_matrix(i,j) * rank_vector(j) )
      end do
      rank2(i) = rank2(i) / team_games(i)
    end do

    rank2(1:team_num) = rank2(1:team_num) &
      / real ( team_games(1:team_num), kind = 8 )

    dx = maxval ( abs ( rank_vector(1:team_num) - rank2(1:team_num) ) )
    xnrm = maxval ( abs ( rank_vector(1:team_num) ) )

    rank_vector(1:team_num) = rank2(1:team_num)

    iteration = iteration + 1

    if ( dx <= picard_tolerance * ( xnrm + 1.0D+00 ) ) then
      exit
    end if

    if ( picard_max < iteration ) then
      exit
    end if

  end do

  return
end
subroutine r8mat_power_method ( lda, n, a, r, v )

!*****************************************************************************80
!
!! R8MAT_POWER_METHOD applies the power method to a matrix.
!
!  Discussion:
!
!    If the power method has not converged, then calling the routine
!    again immediately with the output from the previous call will
!    continue the iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 8 ) A(LDA,N), the matrix.
!
!    Output, real ( kind = 8 ) R, V(N), the estimated eigenvalue
!    and eigenvector.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) av(n)
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: it_eps = 0.0001D+00
  integer ( kind = 4 ), parameter :: it_max = 100
  integer ( kind = 4 ), parameter :: it_min = 10
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) r_old
  real ( kind = 8 ) v(n)

  eps = sqrt ( epsilon ( 1.0D+00 ) )

  r = sqrt ( sum ( v(1:n)**2 ) )

  if ( r == 0.0D+00 ) then
    v(1:n) = 1.0D+00
    r = sqrt ( real ( n, kind = 8 ) )
  end if

  v(1:n) = v(1:n) / r

  do i = 1, it_max

    av(1:n) = matmul ( a(1:n,1:n), v(1:n) )

    r_old = r
    r = sqrt ( sum ( av(1:n)**2 ) )

    if ( it_min < i ) then
      if ( abs ( r - r_old ) <= it_eps * ( 1.0D+00 + abs ( r ) ) ) then
        exit
      end if
    end if

    v(1:n) = av(1:n)

    if ( r /= 0.0D+00 ) then
      v(1:n) = v(1:n) / r
    end if
!
!  Perturb V a bit, to avoid cases where the initial guess is exactly
!  the eigenvector of a smaller eigenvalue.
!
    if ( i < it_max / 2 ) then
      j = 1 + mod ( i-1, n )
      v(j) = v(j) + eps * ( 1.0D+00 + abs ( v(j) ) )
      r2 = sqrt ( sum ( v(1:n)**2 ) )
      v(1:n) = v(1:n) / r2
    end if

  end do

  return
end
subroutine r8vec_sort_insert_index_d ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_INSERT_INDEX_D descending index sorts an R8VEC using insertion.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Algorithm 1.1,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(N), the array to be sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sorted indices.  The array is 
!    sorted when listed from A(INDX(1)) through A(INDX(N)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x

  call i4vec_identity ( n, indx )

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( x <= a(indx(j)) ) then
        exit
      end if

      indx(j+1) = indx(j)
      j = j - 1

    end do

    indx(j+1) = i

  end do

  return
end
subroutine s_before_ss_copy ( s, ss, s2 )

!*****************************************************************************80
!
!! S_BEFORE_SS_COPY copies a string up to a given substring.
!
!  Discussion:
!
!    S and S2 can be the same object, in which case the string is
!    overwritten by a copy of itself up to the substring, followed
!    by blanks.
!
!  Example:
!
!    Input:
!
!      S = 'ABCDEFGH'
!      SS = 'EF'
!
!    Output:
!
!      S2 = 'ABCD'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be copied.
!
!    Input, character ( len = * ) SS, the substring before which the copy stops.
!
!    Output, character ( len = * ) S2, the copied portion of S.
!
  implicit none

  integer ( kind = 4 ) last
  integer ( kind = 4 ) last_s2
  character ( len = * ) s
  character ( len = * ) s2
  character ( len = * ) ss
!
!  Find the first occurrence of the substring.
!
  last = index ( s, ss )
!
!  If the substring doesn't occur at all, behave as though it begins
!  just after the string terminates.
!
!  Now redefine LAST to point to the last character to copy before
!  the substring begins.
!
  if ( last == 0 ) then
    last = len ( s )
  else
    last = last - 1
  end if
!
!  Now adjust again in case the copy holder is "short".
!
  last_s2 = len ( s2 )

  last = min ( last, last_s2 )
!
!  Copy the beginning of the string.
!  Presumably, compilers now understand that if LAST is 0, we don't
!  copy anything.
!  Clear out the rest of the copy.
!
  s2(1:last) = s(1:last)
  s2(last+1:last_s2) = ' '

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
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
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
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
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine s_token_match ( string, token_num, token, match )

!*****************************************************************************80
!
!! S_TOKEN_MATCH matches the beginning of a string and a set of tokens.
!
!  Example:
!
!    Input:
!
!      STRING = 'TOMMYGUN'
!      TOKEN = 'TOM', 'ZEBRA', 'TOMMY', 'TOMMYKNOCKER'
!
!    Output:
!
!      MATCH = 3
!
!  Discussion:
!
!    The longest possible match is taken.
!    Matching is done without regard to case or trailing blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string to be examined.
!
!    Input, integer ( kind = 4 ) TOKEN_NUM, the number of tokens to be compared.
!
!    Input, character ( len = * ) TOKEN(TOKEN_NUM), the tokens.
!
!    Output, integer ( kind = 4 ) MATCH, the index of the (longest) token that 
!    matched the string, or 0 if no match was found.
!
  implicit none

  integer ( kind = 4 ) token_num

  integer ( kind = 4 ) match
  integer ( kind = 4 ) match_length
  logical s_eqi
  character ( len = * ) string
  integer ( kind = 4 ) string_length
  integer ( kind = 4 ) token_i
  integer ( kind = 4 ) token_length
  character ( len = * ) token(token_num)

  match = 0
  match_length = 0

  string_length = len_trim ( string )

  do token_i = 1, token_num

    token_length = len_trim ( token ( token_i ) )

    if ( match_length < token_length ) then

      if ( token_length <= string_length ) then

        if ( s_eqi ( string(1:token_length), &
                     token(token_i)(1:token_length) ) ) then
          match_length = token_length
          match = token_i
        end if

      end if

    end if

  end do

  return
end
subroutine score_matrix_set ( game_num, game_score1, game_score2, game_team1, &
  game_team2, score_option, team_num, score_matrix )

!*****************************************************************************80
!
!! SCORE_MATRIX_SET sets the entries of the score matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2001
!
!  Reference:
!
!    James Keener,
!    The Perron-Frobenius Theorem and the Ranking of Football Teams,
!    SIAM Review, Volume 35, Number 1, pages 80-93, March 1993.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GAME_NUM, the number of games.
!
!    Input, integer ( kind = 4 ) GAME_SCORE1(GAME_NUM), GAME_SCORE2(GAME_NUM),
!    the scores for each game.
!
!    Input, integer ( kind = 4 ) GAME_TEAM1(GAME_NUM), GAME_TEAM2(GAME_NUM),
!    the teams for each game.
!
!    Input, integer ( kind = 4 ) SCORE_OPTION.
!    1, score = 0/0.5/1 for each loss, tie, or win.
!    2, score = Sij
!    3, score = Sij/(Sij+Sji)
!    4, score = (Sij+1)/(Sij+Sji+2)
!    5, score = h(Sij,Sji)
!
!    Input, integer ( kind = 4 ) TEAM_NUM, the number of teams.
!
!    Output, real ( kind = 8 ) SCORE_MATRIX(TEAM_NUM,TEAM_NUM),
!    the scoring matrix.
!
  implicit none

  integer ( kind = 4 ) game_num
  integer ( kind = 4 ) team_num

  real ( kind = 8 ) dij
  real ( kind = 8 ) dji
  integer ( kind = 4 ) game
  integer ( kind = 4 ) game_score1(game_num)
  integer ( kind = 4 ) game_score2(game_num)
  integer ( kind = 4 ) game_team1(game_num)
  integer ( kind = 4 ) game_team2(game_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) score_matrix(team_num,team_num)
  integer ( kind = 4 ) score_option
  integer ( kind = 4 ) sij
  integer ( kind = 4 ) sji
  integer ( kind = 4 ) team_games(team_num)
  real ( kind = 8 ) warp

  score_matrix(1:team_num,1:team_num) = 0.0D+00
  team_games(1:team_num) = 0

  do game = 1, game_num

    i = game_team1(game)
    team_games(i) = team_games(i) + 1
    j = game_team2(game)
    team_games(j) = team_games(j) + 1
    sij = game_score1(game)
    sji = game_score2(game)

    if ( score_option == 1 ) then

      if ( sji < sij ) then
        dij = 1.0D+00
        dji = 0.0D+00
      else if ( sij == sji ) then
        dij = 0.5D+00
        dji = 0.5D+00
      else if ( sij < sji ) then
        dij = 0.0D+00
        dji = 1.0D+00
      end if

    else if ( score_option == 2 ) then

      dij = real ( sij, kind = 8 )
      dji = real ( sji, kind = 8 )

    else if ( score_option == 3 ) then

      if ( sij + sji == 0 ) then
        dij = 0.0D+00
        dji = 0.0D+00
      else
        dij = real ( sij, kind = 8 ) / real ( sij + sji, kind = 8 )
        dji = real ( sji, kind = 8 ) / real ( sij + sji, kind = 8 )
      end if

    else if ( score_option == 4 ) then
      dij = real ( sij + 1, kind = 8 ) / real ( sij + sji + 2, kind = 8 )
      dji = real ( sji + 1, kind = 8 ) / real ( sij + sji + 2, kind = 8 )
    else if ( score_option == 5 ) then
      dij = real ( sij + 1, kind = 8 ) / real ( sij + sji + 2, kind = 8 )
      dij = warp ( dij )
      dji = real ( sji + 1, kind = 8 ) / real ( sij + sji + 2, kind = 8 )
      dji = warp ( dji )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCORE_MATRIX_SET - Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal score option: ', score_option
      stop
    end if

    score_matrix(i,j) = score_matrix(i,j) + dij
    score_matrix(j,i) = score_matrix(j,i) + dji

  end do
!
!  Divide each team's row by the number of games it played.
!
  do i = 1, team_num
    if ( 1 < team_games(i) ) then
      score_matrix(i,1:team_num) = score_matrix(i,1:team_num) &
        / real ( team_games(i), kind = 8 )
    end if
  end do

  return
end
subroutine score_matrix_set2 ( game_num, game_score1, game_score2, game_team1, &
  game_team2, team_num, score_matrix )

!*****************************************************************************80
!
!! SCORE_MATRIX_SET2 sets the entries of the score matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Keener,
!    The Perron-Frobenius Theorem and the Ranking of Football Teams,
!    SIAM Review, Volume 35, Number 1, pages 80-93, March 1993.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GAME_NUM, the number of games.
!
!    Input, integer ( kind = 4 ) GAME_SCORE1(GAME_NUM), GAME_SCORE2(GAME_NUM),
!    the scores for each game.
!
!    Input, integer ( kind = 4 ) GAME_TEAM1(GAME_NUM), GAME_TEAM2(GAME_NUM),
!    the teams for each game.
!
!    Input, integer ( kind = 4 ) TEAM_NUM, the number of teams.
!
!    Output, real ( kind = 8 ) SCORE_MATRIX(TEAM_NUM,TEAM_NUM),
!    the scoring matrix.
!
  implicit none

  integer ( kind = 4 ) game_num
  integer ( kind = 4 ) team_num

  real ( kind = 8 ) dij
  real ( kind = 8 ) dji
  integer ( kind = 4 ) game
  integer ( kind = 4 ) game_score1(game_num)
  integer ( kind = 4 ) game_score2(game_num)
  integer ( kind = 4 ) game_team1(game_num)
  integer ( kind = 4 ) game_team2(game_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) score_count(team_num,team_num)
  real ( kind = 8 ) score_matrix(team_num,team_num)
  integer ( kind = 4 ) sij
  integer ( kind = 4 ) sji
  integer ( kind = 4 ) team_games(team_num)
!
  score_matrix(1:team_num,1:team_num) = 0.0D+00
  team_games(1:team_num) = 0

  do game = 1, game_num

    i = game_team1(game)
    team_games(i) = team_games(i) + 1
    j = game_team2(game)
    team_games(j) = team_games(j) + 1
    sij = game_score1(game)
    sji = game_score2(game)

    dij = ( 5.0D+00 + real ( sij, kind = 8 ) &
                    + real ( sij, kind = 8 )**(2.0D+00/3.0D+00) ) &
        / ( 5.0D+00 + real ( sji, kind = 8 ) &
                    + real ( sij, kind = 8 )**(2.0D+00/3.0D+00) )

    dji = ( 5.0D+00 + real ( sji, kind = 8 ) &
                    + real ( sji )**(2.0D+00/3.0D+00) ) &
        / ( 5.0D+00 + real ( sij, kind = 8 ) &
                    + real ( sji, kind = 8 )**(2.0D+00/3.0D+00) )

    score_matrix(i,j) = score_matrix(i,j) + dij
    score_count(i,j) = score_count(i,j) + 1

    score_matrix(j,i) = score_matrix(j,i) + dji
    score_count(j,i) = score_count(j,i) + 1

  end do
!
!  A WHERE statement would be handy here.
!
  do i = 1, team_num
    do j = 1, team_num
      if ( 1 < score_count(i,j) ) then
        score_matrix(i,j) = score_matrix(i,j) &
          / real ( score_count(i,j), kind = 8 )
      end if
    end do
  end do

  return
end
subroutine team_print ( team_abbrev, team_name, team_num )

!*****************************************************************************80
!
!! TEAM_PRINT prints the team name information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 5 ) TEAM_ABBREV(TEAM_NUM), the abbreviations for
!    each team.
!
!    Input, character ( len = 30 ) TEAM_NAME(TEAM_NUM), the full name of
!    each team.
!
!    Input, integer ( kind = 4 ) TEAM_NUM, the number of teams.
!
  implicit none

  integer ( kind = 4 ) team_num

  character ( len = 5 ) team_abbrev(team_num)
  integer ( kind = 4 ) team_i
  character ( len = 30 ) team_name(team_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ABBREV  Full Name of Team'
  write ( *, '(a)' ) ' '

  do team_i = 1, team_num
    write ( *, '(i3,2x,a,2x,a)' ) team_i, team_abbrev(team_i), team_name(team_i)
  end do

  return
end
subroutine team_read ( iunit, team_abbrev, team_name, team_num )

!*****************************************************************************80
!
!! TEAM_READ reads the team name information from the game score data file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit of the file.
!
!    Output, character ( len = 5 ) TEAM_ABBREV(TEAM_NUM), the abbreviations for
!    each team.
!
!    Output, character ( len = 30 ) TEAM_NAME(TEAM_NUM), the full name of
!    each team.
!
!    Input, integer ( kind = 4 ) TEAM_NUM, the number of teams.
!
  implicit none

  integer ( kind = 4 ) team_num

  integer ( kind = 4 ) iunit
  character ( len = 100 ) line
  character ( len = 5 ) team_abbrev(team_num)
  integer ( kind = 4 ) team_i
  character ( len = 30 ) team_name(team_num)
  character ( len = 10 ) word

  do team_i = 1, team_num
    read ( iunit, '(a)' ) line
    call word_extract ( line, word )
    team_abbrev(team_i) = word
    call s_before_ss_copy ( line, '(', line )
    team_name(team_i) = line
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
  implicit none

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
subroutine token_extract ( string, token_num, token, match )

!*****************************************************************************80
!
!! TOKEN_EXTRACT "extracts" a token from the beginning of a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) STRING; on input, a string from
!    whose beginning a token is to be extracted.  On output,
!    the token, if found, has been removed.
!
!    Input, integer ( kind = 4 ) TOKEN_NUM, the number of tokens to be compared.
!
!    Input, character ( len = * ) TOKEN(TOKEN_NUM), the tokens.
!
!    Output, integer ( kind = 4 ) MATCH, the index of the (longest) token 
!    that matched the string, or 0 if no match was found.
!
  implicit none

  integer ( kind = 4 ) token_num

  integer ( kind = 4 ) left
  integer ( kind = 4 ) match
  character ( len = * ) string
  character ( len = * ) token(token_num)

  call s_token_match ( string, token_num, token, match )

  if ( match /= 0 ) then
    left = len_trim ( token(match) )
    string = adjustl ( string(left+1:) )
  end if

  return
end
function warp ( x )

!*****************************************************************************80
!
!! WARP computes a score that tries to factor out wipeouts.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Keener,
!    The Perron-Frobenius Theorem and the Ranking of Football Teams,
!    SIAM Review, Volume 35, Number 1, pages 80-93, March 1993.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the input argument, the partial score.
!
!    Output, real ( kind = 8 ) WARP, a "warped" value of X.
!
  implicit none

  real ( kind = 8 ) warp
  real ( kind = 8 ) x

  if ( x < 0.5D+00 ) then
    warp = 0.5D+00 * ( 1.0D+00 - sqrt ( abs ( 2.0D+00 * x - 1.0D+00 ) ) )
  else if ( x == 0.5D+00 ) then
    warp = 0.5D+00
  else if ( 0.5D+00 < x ) then
    warp = 0.5D+00 * ( 1.0D+00 + sqrt ( abs ( 2.0D+00 * x - 1.0D+00 ) ) )
  end if

  return
end
function warp2 ( x )

!*****************************************************************************80
!
!! WARP2 computes a score that tries to factor out wipeouts.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Keener,
!    The Perron-Frobenius Theorem and the Ranking of Football Teams,
!    SIAM Review, Volume 35, Number 1, pages 80-93, March 1993.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the input argument, the partial score.
!
!    Output, real ( kind = 8 ) WARP2, a "warped" value of X.
!
  implicit none

  real ( kind = 8 ) warp2
  real ( kind = 8 ) x

  warp2 = x * ( x + 0.05D+00 ) / ( x**2 + 0.05D+00 * x + 2.0D+00 )

  return
end
subroutine word_extract ( s, w )

!*****************************************************************************80
!
!! WORD_EXTRACT extracts the next word from a string.
!
!  Discussion:
!
!    A "word" is a string of characters terminated by a blank or
!    the end of the string.
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
!    Input/output, character ( len = * ) S, the string.  On output, the first
!    word has been removed, and the remaining string has been shifted left.
!
!    Output, character ( len = * ) W, the leading word of the string.
!
  implicit none

  integer ( kind = 4 ) iget1
  integer ( kind = 4 ) iget2
  integer ( kind = 4 ) lchar
  character ( len = * ) s
  character ( len = * ) w

  w = ' '

  lchar = len ( s )
!
!  Find the first nonblank.
!
  iget1 = 0

  do

    iget1 = iget1 + 1

    if ( lchar < iget1 ) then
      return
    end if

    if ( s(iget1:iget1) /= ' ' ) then
      exit
    end if

  end do
!
!  Look for the last contiguous nonblank.
!
  iget2 = iget1

  do

    if ( lchar <= iget2 ) then
      exit
    end if

    if ( s(iget2+1:iget2+1) == ' ' ) then
      exit
    end if

    iget2 = iget2 + 1

  end do
!
!  Copy the word.
!
  w = s(iget1:iget2)
!
!  Shift the string.
!
  s(1:iget2) = ' '
  s = adjustl ( s(iget2+1:) )

  return
end
