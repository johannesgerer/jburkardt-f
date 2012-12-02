program main

!*****************************************************************************80
!
!! DUEL_SIMULATION simulates a duel.
!
!  Discussion:
!
!    Player 1 fires at player 2, and hits with a probability of P(1).
!    If Player 2 misses, then Player 2 fires at Player 1, hitting with
!    a probability of P(2).
!
!    The duel continues with alternating shots until only one player 
!    survives.
!
!    The simulation is intended to estimate the probabilities that a
!    player will survive, and the number of turns required.
!
!    The exact probability that player 1 will survive is
!
!      P(1) / ( P(1) + P(2) - P(1) * P(2) )
!
!    Player 2's chance is
!
!     P(2) * ( 1 - P(1) ) / ( P(1) + P(2) - P(1) * P(2) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Nahin,
!    Duelling Idiots and Other Probability Puzzlers,
!    Princeton University Press, 2000,
!    ISBN13: 978-0691009797,
!    LC: QA273.N29.
!
!    Martin Shubik,
!    "Does the Fittest Necessarily Survive?",
!    in Readings in Game Theory and Political Behavior,
!    edited by Martin Shubik,
!    Doubleday, 1954,
!    LC: H61.S53.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A_ACCURACY, B_ACCURACY, the probabilities that 
!    players A and B will hit their opponent in a single shot.
!
!    Input, integer ( kind = 4 ) DUEL_NUM, the number of duels to run.
!
!    Output, real ( kind = 8 ) A_PROB, B_PROB, the estimated probablities that 
!    players A and B will survive.
!
!    Output, real ( kind = 8 ) TURN_AVERAGE, the average number of turns 
!    required to complete the duel.
!
  implicit none

  real ( kind = 8 ) a_accuracy
  real ( kind = 8 ) a_prob
  integer ( kind = 4 )  a_wins
  real ( kind = 8 ) b_accuracy
  real ( kind = 8 ) b_prob
  integer ( kind = 4 )  b_wins
  integer ( kind = 4 ) duel
  integer ( kind = 4 ) duel_num
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) winner

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'DUEL_SIMULATION:'
  write ( *, '(a)' ) '  FORTRAN90 version'

  write ( *, '(a)' ) 'Enter number of duels to run:'
  read ( *, * ) duel_num

  write ( *, '(a)' ) 'Enter player A''s accuracy between 0.0 and 1.0:'
  read ( *, * ) a_accuracy

  write ( *, '(a)' ) 'Enter player B''s accuracy between 0.0 and 1.0:'
  read ( *, * ) b_accuracy

  a_wins = 0
  b_wins = 0

  do duel = 1, duel_num

    call duel_result ( a_accuracy, b_accuracy, winner )

    if ( winner == 1 ) then
      a_wins = a_wins + 1
    else
      b_wins = b_wins + 1
    end if

  end do

  a_prob = real ( a_wins, kind = 8 ) / real ( duel_num, kind = 8 )
  b_prob = real ( b_wins, kind = 8 ) / real ( duel_num, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Player A wins with probability ', a_prob
  write ( *, '(a,g14.6)' ) '  Player B wins with probability ', b_prob

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'DUEL_SIMULATION:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine duel_result ( a_accuracy, b_accuracy, winner )

!*****************************************************************************80
!
!! DUEL_RESULT returns the outcome of a single duel.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Martin Shubik,
!    “Does the Fittest Necessarily Survive?”,
!    in Readings in Game Theory and Political Behavior,
!    edited by Martin Shubik,
!    Doubleday, 1954,
!    LC: H61.S53.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A_ACCURACY, B_ACCURACY, the probabilities that 
!    player A and B will hit their opponent in a single shot.
!
!    Output, integer ( kind = 4 ) WINNER, the survivor of the duel.
!
  implicit none

  real ( kind = 8 ) a_accuracy
  real ( kind = 8 ) b_accuracy
  real ( kind = 8 ) r
  integer ( kind = 4 ) winner

  do

    call random_number ( harvest = r )

    if ( r <= a_accuracy ) then
      winner = 1
      exit
    end if

    call random_number ( harvest = r )

    if ( r <= b_accuracy ) then
      winner = 2
      exit
    end if

  end do

  return
end
