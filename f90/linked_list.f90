module linked_list_type

!*****************************************************************************80  
!
!! LINKED_LIST_TYPE defines a "record" of a linked list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2009
!
!  Author:
!
!    John Burkardt
!
  type record
    integer :: generation
    integer :: value
    type ( record ), pointer :: previous
    type ( record ), pointer :: next
  end type record

end module linked_list_type

module linked_list_library

  use linked_list_type

  contains

subroutine list_insert ( item, head )

!*****************************************************************************80  
!
!! LIST_INSERT inserts ITEM into the linked list pointed to by HEAD.
!
!  Discussion:
!
!    This routine requires a user defined linked-list data type "RECORD".
!
!    The items in the linked list should be ascending sorted by the VALUE component
!    of the individual records.
!
!    This routine inserts the new item to preserve this ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, type ( record ) pointer :: ITEM, a pointer to the record
!    to be inserted into the linked list.
!
!    Input/output, type ( record ) pointer :: HEAD, a pointer to the first record
!    in the linked list.
!
  implicit none

  type ( record ), pointer :: head
  integer i
  type ( record ), pointer :: item
  type ( record ), pointer :: item1
  type ( record ), pointer :: item2
!
!  In the case of an empty list.  
!
  if ( .not. associated ( head ) ) then
    head => item
    nullify ( item%previous )
    nullify ( item%next )
    return
  end if
!
!  In the special case that ITEM < HEAD, we need to update HEAD.
!
  if ( item%value < head%value ) then
    nullify ( item%previous )
    item%next => head
    head%previous => item
    head => item
    return
  end if

  nullify ( item1 )
  item2 => head

  do while ( item2%value < item%value )

    item1 => item2;
    item2 => item2%next

    if ( .not. associated ( item2 ) ) then
      exit
    end if

  end do
!
!  ITEM1 < ITEM and either ITEM <= ITEM2 or ITEM2 is null.
!
  if ( associated ( item2 ) ) then
    item2%previous => item
    item%next => item2
  else
    nullify ( item%next )
  end if
 
  item1%next => item
  item%previous => item1

  return
end subroutine
subroutine list_print ( head )

!*****************************************************************************80  
!
!! LIST_PRINT prints a linked list.
!
!  Discussion:
!
!    This routine requires a user defined linked-list data type "RECORD".
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
!  Parameters:
!
!    Input, type ( record ) pointer :: HEAD, a pointer to the first record
!    in the linked list.
!
  implicit none

  type ( record ), pointer :: head
  integer i
  type ( record ), pointer :: item

  i = 0
  item => head

  do while ( associated ( item ) )

    i = i + 1
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, item%generation, item%value

    item => item%next

  end do
 
  return
end subroutine

end module linked_list_library

program main

!*****************************************************************************80
!
!! MAIN is the main program for the linked list example.
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

  integer data_num

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINKED_LIST:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate how pointers can be used to define'
  write ( *, '(a)' ) '  and manipulate a linked list.'

  data_num = 10
  call test01 ( data_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINKED_LIST:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( data_num )

!*****************************************************************************80
!
!! TEST01 uses a linked list to store and sort random data.
!
!  Discussion:
!
!    This routine requires a user defined linked-list library.
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
  use linked_list_library

  implicit none

  integer data_num
  type ( record ), pointer :: head
  integer i
  type ( record ), pointer :: item
  real r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Create, one at a time, a sequence of integers.'
  write ( *, '(a)' ) '  As each integer is created, insert it into a sorted'
  write ( *, '(a)' ) '  list.  Print the list at the end.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We carry out this task using a linked list.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initial data generation:'
  write ( *, '(a)' ) ' '

  nullify ( head )

  do i = 1, data_num
!
!  Generate a new item.
!
    allocate ( item )

    item%generation = i

    call random_number ( harvest = r )

    item%value = int ( 1000.0 * r )
    nullify ( item%next )
    nullify ( item%previous )

    write ( *, '(2x,i8,2x,i8)' ) i, item%value
!
!  Insert the new item into the linked list.
!
    call list_insert ( item, head )

  end do
!
!  Print the linked list.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Contents of sorted linked list:'
  write ( *, '(a)' ) ' '

  call list_print ( head )

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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
