subroutine simdo ( qind, qfor, iprod, kdim, jsub, ivec, ifault )

!*****************************************************************************80
!
!! SIMDO generates multi-indices, simulating nested DO-loops.
!
!  Discussion:
!
!    The loops are assumed to be nested to a depth of K.
!
!    The R-th loop is assumed to have upper limit N(R) and increment Inc(R).
!
!    The total number of executions of the innermost loop is 
!
!      N = product ( 1 <= R <= K ) N(R).
!
!    Let these executions be indexed by the single integer J, which
!    we call the index subscript.
!
!    Each value of J corresponds to a particular set of loop indices,
!    which we call the subscript vector I(J).
!
!    This routine can start with J and find I(J), or determine
!    J from I(J).
!    
!  Modified:
!
!    26 July 2008
!
!  Author:
!
!    M O'Flaherty, G MacKenzie
!
!  Reference:
!
!    M O'Flaherty, G MacKenzie,
!    Algorithm AS 172:
!    Direct Simulation of Nested Fortran DO-LOOPS,
!    Applied Statistics,
!    Volume 31, Number 1, 1982, pages 71-74.
!
!  Parameters:
!
!    Input, logical QIND.
!    TRUE to convert an index subscript J to the subscript vector I(J).
!    FALSE to convert the subscript vector I(J) to the index subscript J.
!
!    Input, logical QFOR,
!    TRUE if conversion is required in standard Fortran subscripting order,
!    FALSE otherwise.
!
!    Input, integer IPROD(KDIM), contains the partial products.
!    If QFOR is FALSE, then
!      IPROD(S) = product ( 1 <= R <= S ) N(R).
!    If QFOR is TRUE, then
!      IPROD(S) = product ( 1 <= R <= S ) N(KDIM+1-R).
!
!    Input, integer KDIM, the nesting depth of the loops.
!
!    Input/output, integer JSUB.
!    If QIND is TRUE, then JSUB is an input quantity, an index subscript J
!    to be converted into the subscript vector I(J).
!    If QIND is FALSE, then JSUB is an output quantity, the index subscript J
!    corresponding to the subscript vector I(J).
!
!    Input/output, integer IVEC(KDIM).
!    if QIND is TRUE, then IVEC is an output quantity, the subscript vector I(J)
!    corresponding to the index subscript J.
!    If QIND is FALSE, then IVEC is an input quantity, a subscript vector I(J)
!    for which the corresponding index subscript J is to be computed.
!
!    Output, integer IFAULT, error flag.
!    0, no error was detected.
!    1, if QIND is TRUE, and the input value of JSUB exceeds IPROD(KDIM).
!    2, if QIND is FALSE, and IVEC contains an illegal component.
!
  implicit none

  integer kdim

  integer i
  integer ifault
  integer ik
  integer iprod(kdim)
  integer itempv
  integer ivec(kdim)
  integer jsub
  logical qfor
  logical qind

  ifault = 0
!
!  Index subscript to subscript vector conversion.
!
  if ( qind )  then

    if ( iprod(kdim) < jsub ) then
      ifault = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIMDO - Fatal error!'
      write ( *, '(a)' ) '  JSUB is out of bounds.'
      stop
    end if

    itempv = jsub - 1

    do i = 1, kdim - 1
      ik = kdim - i
      ivec(i) = itempv / iprod(ik)
      itempv = itempv - iprod(ik) * ivec(i)
      ivec(i) = ivec(i) + 1
    end do

    ivec(kdim) = itempv + 1
    if ( qfor ) then
      ivec(1:kdim) = ivec(kdim:1:-1)
    end if
!
!  Subscript vector to index subscript conversion.
!
  else

    if ( .not. qfor ) then
      ivec(1:kdim) = ivec(kdim:1:-1)
    end if

    if ( iprod(1) < ivec(1) ) then
      ifault = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIMDO - Fatal error!'
      write ( *, '(a)' ) '  An entry of IVEC is out of bounds.'
      stop
    end if

    do i = 2, kdim
      if ( iprod(i) / iprod(i-1) < ivec(i) ) then
        ifault = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SIMDO - Fatal error!'
        write ( *, '(a)' ) '  An entry of IVEC is out of bounds.'
        stop
      end if
    end do

    jsub = ivec(1)

    do i = 2, kdim
      jsub = jsub + ( ivec(i) - 1 ) * iprod(i-1)
    end do
!
!  As a courtesy to the caller, UNREVERSE the IVEC vector
!  if you reversed it!
!
    if ( .not. qfor ) then
      ivec(1:kdim) = ivec(kdim:1:-1)
    end if

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
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
