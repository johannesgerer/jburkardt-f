program main

!*****************************************************************************80
!
!! THIRTY_ONE simulates the game of "31".
!
!  Discussion:
!
!    Two players draw alternately from a faceup tableau of cards, consisting
!    of four copies each of the numbers 1 through 6.
!
!    The total of all cards taken by both players is the score.
!
!    A player wins by picking a card which results in the score reaching 
!    exactly 31.  
!
!    If the score goes over 31 without hitting that value exactly, the game
!    is a draw.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2011
!
!  Author:
!
!    Yuen-Yick Kwan
!
  implicit none

  integer :: cards(6)
  integer :: k
  logical :: legal
  integer :: next
  integer :: total

  cards = 4
  total = 0
!
!  Each iteration simulates one move by the player and one by the computer.
!
  do

    print *, 'Cards left:'

    call print_cards ( cards )
!
!  Player selects a card.
!
    legal = .false.

    do while ( .not. legal )
      write(*,'(a)',advance='no') ' Your pick: '
      read(*,*) next
      if ( 1 <= next .and. next <= 6 ) then
        legal = ( 0 < cards(next) .and. total + next <= 31 )
      end if
    end do

    cards(next) = cards(next) - 1
    total = total + next

    if ( total == 31 ) then
      print *, 'Your win!'
      exit
    end if
!
!  Computer selects a card.
!
    next = next_play ( total, cards )

    select case ( next )
      case ( -6 : -1 )
        next = abs ( next )
      case ( -7 )
        print *, 'Draw!'
        exit
      case ( 0 )
        next = 31-total
        do while (next > 6)
          next = next-7
        end do

        if ( cards(next) == 0 ) then
          next = 1
        end if

    end select

    cards(next) = cards(next) - 1
    total = total + next
    write(*,'(a,i1,a,i2)') ' I pick ', next, ', Total = ', total
    print *, ''

    if ( total == 31 ) then
      print *, 'I win!'
      exit
    end if

    legal = .false.

     do k = 1, 6
       if ( 0 < cards(k) .and. total + k <= 31 ) then
         legal = .true.
       end if
     end do

     if ( .not. legal ) then
       print *, 'Draw!'
       exit
     end if

   end do

   stop

   contains

   subroutine print_cards ( cards )

!*****************************************************************************80
!
!! PRINT_CARDS prints the current card tableau.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2011
!
!  Author:
!
!    Yuen-Yick Kwan
!
      implicit none
      integer, intent(in) :: cards(6)
      integer :: i, k

      do k = 1, 4
         do i = 1, 6
            if (cards(i) >= k) then
               write(*,'(i3)',advance='no') i
            else
               write(*,'(a3)',advance='no') ''
            end if
         end do
         print *, ''
      end do
      print *, ''

    return

  end subroutine print_cards

  recursive function next_play ( total, cards ) result ( next )

!*****************************************************************************80
!
!! NEXT_PLAY chooses the computer's next move.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2011
!
!  Author:
!
!    Yuen-Yick Kwan
!
!  Parameters:
!
!    Input, integer TOTAL, the computer's current total.
!
!    Input, integer CARDS(6), the current card tableau.
!
  implicit none

  integer, intent(in) :: cards(6)
  integer :: k
  logical :: legal
  integer :: next
  integer :: tcards(6)
  integer, intent(in) :: total
  integer :: trial
!
!  Determine the distance to 31.
!
  next = 31 - total
!
!  If a winning card is available, pick it.
!
  if ( next <= 6 ) then
    if ( 0 < cards(next) ) then
      return
    end if
  end if
!
!  Search for a legal move that doesn't give the opponent an immediate win.
!
  tcards = cards
  legal = .false.
  next = 0

  do k = 6, 1, -1

    if ( 0 < cards(k) .and. total + k < 31 ) then

      legal = .true.
      tcards(k) = tcards(k) - 1
      trial = next_play ( total + k, tcards )
      tcards(k) = tcards(k) + 1
!
!  If the computer takes card K, and the opponent can't win on the next turn,
!  then take card K.
!
      if ( trial == 0 ) then
        next = k
        return
!
!  If the computer takes card K and the opponent is forced over 31 on the next
!  turn, consider taking card K (but keep searching!)
!
      else if ( next == 0 .and. trial < 0 ) then
        next = -k
      end if

    end if

  end do
!
!  Return NEXT = - 7 if all remaining moves take the total over 31.
!
  if ( .not. legal ) then
    next = - 7
  end if

  return

  end function next_play

end program
