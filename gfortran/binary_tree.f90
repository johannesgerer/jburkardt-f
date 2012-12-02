module binary_tree_type

  type record
    integer :: generation
    integer :: value
    type ( record ), pointer :: parent
    type ( record ), pointer :: left
    type ( record ), pointer :: right
  end type record

end module binary_tree_type

module binary_tree_library

  use binary_tree_type

  contains

subroutine binary_tree_insert ( item, head )

!*****************************************************************************80  
!
!! BINARY_TREE_INSERT inserts ITEM into the binary tree pointed to by HEAD.
!
!  Discussion:
!
!    This routine requires a binary tree data type "RECORD".
!
!    The items already in the tree should be ascending sorted by the VALUE 
!    component of the individual records.  Lower items go to the left of
!    larger items.
!
!    This routine inserts the new item to preserve this ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, type ( record ) pointer :: ITEM, a pointer to the record
!    to be inserted into the binary tree.
!
!    Input/output, type ( record ) pointer :: HEAD, a pointer to the first 
!    record in the binary.
!
  implicit none

  type ( record ), pointer :: head
  integer i
  type ( record ), pointer :: item
  type ( record ), pointer :: item1
  type ( record ), pointer :: item2

  nullify ( item%parent )
  nullify ( item%left )
  nullify ( item%right )
!
!  In the case of an empty tree.  
!
  if ( .not. associated ( head ) ) then
    head => item
    return
  end if
!
!  ITEM is being compared to ITEM1.
!  If ITEM < ITEM1, then
!     if ITEM1%LEFT is null, ITEM1%LEFT = ITEM and return.
!     if ITEM1%LEFT is non null, set ITEM1 to ITEM1%LEFT and do again.
!  If ITEM1 < ITEM, then do the same, but to the right.
!
  item1 => head

  do

    if ( item%value <= item1%value ) then

      if ( .not. associated ( item1%left ) ) then
        item1%left => item
        item%parent => item1
        exit
      else
        item1 => item1%left
      end if

    else

      if ( .not. associated ( item1%right ) ) then
        item1%right => item
        item%parent => item1
        exit
      else
        item1 => item1%right
      end if

    end if

  end do

  return
end subroutine
recursive subroutine binary_tree_print ( head )

!*****************************************************************************80  
!
!! BINARY_TREE_PRINT prints a binary tree.
!
!  Discussion:
!
!    This routine requires a binary tree data type "RECORD".
!
!    The binary is presumed to be sorted.  For any node, all the data on
!    subnodes to the left is smaller, and all the data on subnodes to the
!    right is greater.
!
!    Therefore, we can print all the data in the tree or subtree, in order,
!    by printing the left data, the current data, and the right data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, type ( record ) pointer :: HEAD, a pointer to the first record
!    in the binary tree (or subtree).
!
  implicit none

  type ( record ), pointer :: head

  if ( associated ( head ) ) then

    call binary_tree_print ( head%left )
    write ( *, '(2x,i8,2x,i8)' ) head%generation, head%value
    call binary_tree_print ( head%right )

  end if
 
  return
end subroutine

end module binary_tree_library

program main

!*****************************************************************************80
!
!! MAIN is the main program for the binary tree example.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2007
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
  write ( *, '(a)' ) 'BINARY_TREE:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate how pointers can be used to define'
  write ( *, '(a)' ) '  and manipulate a binary tree.'

  data_num = 10
  call test01 ( data_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BINARY_TREE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( data_num )

!*****************************************************************************80
!
!! TEST01 uses a binary tree to store and sort random data.
!
!  Discussion:
!
!    This routine requires a binary tree library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2007
!
!  Author:
!
!    John Burkardt
!
  use binary_tree_library

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
  write ( *, '(a)' ) '  binary tree.  Print the binary tree at the end.'
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

    write ( *, '(2x,i8,2x,i8)' ) i, item%value
!
!  Insert the new item into the linked list.
!  The INSERT routine takes care of initializing the other fields in the
!  ITEM data.
!
    call binary_tree_insert ( item, head )

  end do
!
!  Print the binary tree.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Contents of sorted binary tree:'
  write ( *, '(a)' ) ' '

  call binary_tree_print ( head )

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
!    26 February 2005
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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

