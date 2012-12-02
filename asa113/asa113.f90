subroutine swap ( varval, class, clsize, in, ik, iv, critvl, ntrans, ifault )

!*****************************************************************************80
!
!! SWAP interchanges objects between different classes to improve a criterion.
!
!  Discussion:
!
!    This routine is given a classification of objects, including the
!    number of objects in each class, and the current value of some criterion
!    which is desired to be minimized.
!
!    The routine calculates the change in criterion for all possible swaps,
!    that is, operations in which two objects in different classes exchange
!    places. Each swap that would result in a lowering of the criterion is
!    executed, and the related quantities are updated.
!
!    When no more advantageous swaps can be found, the routine returns.
!
!    The routine relies on a user-supplied routine, CRSWAP, to report the
!    expected change in the criterion for a given swap, and to carry
!    out that transfer if requested.
!
!    The variables CLASS and CRITVL have been added to the argument list
!    of CRSWAP.
!
!    Also, the order of the two classes "L" and "M" was interchanged in
!    the call to CRSWAP.  The original order was counterintuitive.
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Banfield, Bassill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Colin Banfield, LC Bassill,
!    Algorithm AS 113:
!    A transfer for non-hierarchichal classification,
!    Applied Statistics,
!    Volume 26, Number 2, 1977, pages 206-210.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) VARVAL(IN,IV), the data values.  There are
!    IN objects, each having spatial dimension IV.
!
!    Input/output, integer ( kind = 4 ) CLASS(IN), the classification of
!    each object.
!
!    Input/output, integer ( kind = 4 ) CLSIZE(IK), the number of objects
!    in each class.
!
!    Input, integer ( kind = 4 ) IN, the number of objects.
!
!    Input, integer ( kind = 4 ) IK, the number of classes.
!
!    Input, integer ( kind = 4 ) IV, the number of spatial dimensions,
!    or variates, of the objects.
!
!    Input/output, real ( kind = 8 ) CRITVL, the current value of the criterion.
!
!    Output, integer ( kind = 4 ) NTRANS, the number of transfers executed.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no error detected.
!    1, the number of classes was less than 2.
!    2, the number of objects was less than the number of classes.
!
  implicit none

  integer ( kind = 4 ) ik
  integer ( kind = 4 ) in
  integer ( kind = 4 ) iv

  integer ( kind = 4 ) class(in)
  integer ( kind = 4 ) clsize(ik)
  real ( kind = 8 ) critvl
  real ( kind = 8 ), parameter :: eps = 1.0D-38
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ifault
  real ( kind = 8 ) inc
  integer ( kind = 4 ) iswitch
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ntrans
  real ( kind = 8 ) varval(in,iv)

  if ( ik <= 1 ) then
    ifault = 1
    return
  end if

  if ( in <= ik ) then
    ifault = 2
    return
  end if

  ifault = 0
  icount = 0
  ntrans = 0
  itop = ( in * ( in - 1 ) ) / 2

  i = 1

  do

    i = i + 1

    if ( itop <= icount ) then
      exit
    end if

    if ( in < i ) then
      i = 1
      cycle
    end if

    l = class(i)
    k = l
    it = i - 1
!
!  Test the swap of object I from class M to L,
!  and object J from class L to M.
!
    do j = 1, it

      icount = icount + 1
      m = class(j)

      if ( l /= j ) then

        if ( clsize(l) /= 1 .or. clsize(m) /= 1 ) then

          iswitch = 1
          call crswap ( varval, class, clsize, in, ik, iv, critvl, &
            i, j, l, m, iswitch, inc )

          if ( inc < - eps ) then

            critvl = critvl + inc
            icount = 0

            iswitch = 2
            call crswap ( varval, class, clsize, in, ik, iv, critvl, &
              i, j, l, m, iswitch, inc )

            ntrans = ntrans + 1
            class(i) = m
            class(j) = l
            l = m

          end if

        end if
      end if

    end do

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
!    31 May 2001   9:45:54.872 AM
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
subroutine trnsfr ( varval, class, clsize, in, ik, iv, critvl, ntrans, ifault )

!*****************************************************************************80
!
!! TRNSFR transfers objects between classes to improve a criterion.
!
!  Discussion:
!
!    This routine is given a classification of objects, including the
!    number of objects in each class, and the current value of some criterion
!    which is desired to be minimized.
!
!    The routine calculates the change in criterion for all possible transfers
!    of any object from its current class to a different class.  Each transfer
!    that would result in a lowering of the criterion is executed, and the
!    related quantities are updated.
!
!    When no more advantageous transfers can be found, the routine returns.
!
!    The routine relies on a user-supplied routine, CRTRAN, to report the
!    expected change in the criterion for a given transfer, and to carry
!    out that transfer if requested.
!
!    The variables CLASS and CRITVL have been added to the argument list
!    of CRTRAN.
!
!    Also, the order of the two classes "L" and "M" was interchanged in
!    the call to CRTRAN.  The original order was counterintuitive.
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Banfield, Bassill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Colin Banfield, LC Bassill,
!    Algorithm AS 113:
!    A transfer for non-hierarchichal classification,
!    Applied Statistics,
!    Volume 26, Number 2, 1977, pages 206-210.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) VARVAL(IN,IV), the data values.  There are IN
!    objects, each having spatial dimension IV.
!
!    Input/output, integer ( kind = 4 ) CLASS(IN), the classification of
!    each object.
!
!    Input/output, integer ( kind = 4 ) CLSIZE(IK), the number of objects in
!    each class.
!
!    Input, integer ( kind = 4 ) IN, the number of objects.
!
!    Input, integer ( kind = 4 ) IK, the number of classes.
!
!    Input, integer ( kind = 4 ) IV, the number of spatial dimensions, or
!    variates, of the objects.
!
!    Input/output, real ( kind = 8 ) CRITVL, the current value of the criterion.
!
!    Output, integer ( kind = 4 ) NTRANS, the number of transfers executed.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no error detected.
!    1, the number of classes was less than 2.
!    2, the number of objects was less than the number of classes.
!
  implicit none

  integer ( kind = 4 ) ik
  integer ( kind = 4 ) in
  integer ( kind = 4 ) iv

  integer ( kind = 4 ) class(in)
  integer ( kind = 4 ) clsize(ik)
  real ( kind = 8 ) critvl
  real ( kind = 8 ), parameter :: eps = 1.0D-38
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ifault
  real ( kind = 8 ) inc
  real ( kind = 8 ) inco
  integer ( kind = 4 ) iswitch
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ntrans
  real ( kind = 8 ) varval(in,iv)

  if ( ik <= 1 ) then
    ifault = 1
    return
  end if

  if ( in <= ik ) then
    ifault = 2
    return
  end if

  ifault = 0
  ntrans = 0
  i = 0
  icount = 0

  do

    i = i + 1

    if ( in <= icount ) then
      exit
    end if

    if ( in < i ) then
      i = 0
      icount = 0
      cycle
    end if

    m = class(i)
    if ( clsize(m) <= 1 ) then
      icount = icount + 1
      cycle
    end if

    inco = - eps
    lo = m
!
!  Test the transfer of object I from class M to class L.
!
    do l = 1, ik

      if ( l /= m ) then

        iswitch = 1
        call crtran ( varval, class, clsize, in, ik, iv, critvl, &
          i, m, l, iswitch, inc )
!
!  Remember the values of L and INC.
!
        if ( inc < inco ) then
          lo = l
          inco = inc
        end if

      end if

    end do

    icount = icount + 1
!
!  Execute the transfer of object I from class M to class LO.
!
    if ( lo /= m ) then

      l = lo
      critvl = critvl + inco
      icount = 0

      iswitch = 2
      call crtran ( varval, class, clsize, in, ik, iv, critvl, &
        i, m, l, iswitch, inc )

      ntrans = ntrans + 1
      class(i) = l
      clsize(l) = clsize(l) + 1
      clsize(m) = clsize(m) - 1

    end if

  end do

  return
end
