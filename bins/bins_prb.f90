program main

!*****************************************************************************80
!
!! BINS_PRB tests routines from the BINS library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BINS_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BINS library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test05 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BINS_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BIN_TO_R8_EVEN, R8_TO_BIN_EVEN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: a = 10.0D+00
  real ( kind = 8 ), parameter :: b = 20.0D+00
  integer ( kind = 4 ) bin
  real ( kind = 8 ) c
  real ( kind = 8 ) cmax
  real ( kind = 8 ) cmin
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: nbin = 7
  real ( kind = 8 ) r8_uniform
  real ( kind = 8 ), parameter :: rmax = 23.0D+00
  real ( kind = 8 ), parameter :: rmin = 8.0D+00
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  R8_TO_BIN_EVEN puts a number into a bin.'
  write ( *, '(a)' ) '  BIN_TO_R8_EVEN returns the bin limits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The bins are equally spaced between A and B,'
  write ( *, '(a)' ) '  with two extra bins, for things less than A,'
  write ( *, '(a)' ) '  or greater than B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a,i6)' ) '  Total number of bins = ', nbin

  call get_seed ( seed )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Using random seed = ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate some random values C and put them in bins.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       C      Bin   Bin_Min  Bin_Max'
  write ( *, '(a)' ) ' '

  do i = 1, 30
    c = r8_uniform ( rmin, rmax, seed )
    call r8_to_bin_even ( nbin, a, b, c, bin )
    call bin_to_r8_even ( nbin, bin, a, b, cmin, cmax )
    write ( *, '(2x,g14.6,i4,2g14.6)' ) c, bin, cmin, cmax
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests BIN_TO_R8_EVEN2, R8_TO_BIN_EVEN2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: a = 10.0D+00
  real ( kind = 8 ), parameter :: b = 20.0D+00
  integer ( kind = 4 ) bin
  real ( kind = 8 ) c
  real ( kind = 8 ) cmax
  real ( kind = 8 ) cmin
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: nbin = 5
  real ( kind = 8 ) r8_uniform
  real ( kind = 8 ), parameter :: rmax = 23.0D+00
  real ( kind = 8 ), parameter :: rmin = 8.0D+00
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  BIN_TO_R8_EVEN2 returns the bin limits.'
  write ( *, '(a)' ) '  R8_TO_BIN_EVEN2 puts a number into a bin.'
  write ( *, '(a)' ) '  The bins are equally spaced between A and B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a,i6)' ) '  Total number of bins = ', nbin

  call get_seed ( seed )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Using random seed = ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate some random values C and put them in bins.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       C      Bin   Bin_Min  Bin_Max'
  write ( *, '(a)' ) ' '

  do i = 1, 30
    c = r8_uniform ( rmin, rmax, seed )
    call r8_to_bin_even2 ( nbin, a, b, c, bin )
    call bin_to_r8_even2 ( nbin, bin, a, b, cmin, cmax )
    write ( *, '(2x,g14.6,i4,2g14.6)' ) c, bin, cmin, cmax
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests BIN_TO_R82_EVEN, R82_TO_BIN_EVEN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 2 ) :: a = (/  5.0D+00,  0.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 2 ) :: b = (/ 15.0D+00, 20.0D+00 /)
  integer ( kind = 4 ) bin(2)
  real ( kind = 8 ) c(2)
  real ( kind = 8 ) cmax(2)
  real ( kind = 8 ) cmin(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: nbin = 7
  real ( kind = 8 ), parameter, dimension ( 2 ) :: rmin = (/  &
    3.0D+00, -2.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 2 ) :: rmax = (/ &
    23.0D+00, 21.0D+00 /)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  BIN_TO_R82_EVEN returns the bin limits.'
  write ( *, '(a)' ) '  R82_TO_BIN_EVEN puts a R82 number into a bin.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The bins are equally spaced between A and B,'
  write ( *, '(a)' ) '  with two extra bins, for things less than A,'
  write ( *, '(a)' ) '  or greater than B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A(1) = ', a(1)
  write ( *, '(a,g14.6)' ) '  B(1) = ', b(1)
  write ( *, '(a,g14.6)' ) '  A(2) = ', a(2)
  write ( *, '(a,g14.6)' ) '  B(2) = ', b(2)
  write ( *, '(a,i6)' ) '  Total number of bins = ', nbin
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate some random values C and put them in bins.'
  write ( *, '(a)' ) '  We list the X and Y components on separate lines.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       C      Bin   Bin_Min  Bin_Max'
  write ( *, '(a)' ) ' '

  call get_seed ( seed )

  do i = 1, 30
    call r82_uniform ( rmin, rmax, seed, c )
    call r82_to_bin_even ( nbin, a, b, c, bin )
    call bin_to_r82_even ( nbin, bin, a, b, cmin, cmax )
    write ( *, '(a)' ) ' '
    write ( *, '(2x,g14.6,i4,2g14.6)' ) c(1), bin(1), cmin(1), cmax(1)
    write ( *, '(2x,g14.6,i4,2g14.6)' ) c(2), bin(2), cmin(2), cmax(2)
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests R82VEC_BIN_EVEN, R82VEC_BINNED_REORDER, R82VEC_BINNED_SORT_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 30
  integer ( kind = 4 ), parameter :: nbin = 4

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ), parameter, dimension ( 2 ) :: amin = (/ &
    8.0D+00, 3.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 2 ) :: amax = (/ &
    23.0D+00, 12.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 2 ) :: bin_min = (/ &
    10.0D+00, 5.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 2 ) :: bin_max = (/ &
    20.0D+00, 10.0D+00 /)
  integer ( kind = 4 ) bin_last(nbin,nbin)
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin,nbin)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  R82VEC_BIN_EVEN constructs evenly spaced bins and'
  write ( *, '(a)' ) '    assigns each element of a R82VEC to a bin.'
  write ( *, '(a)' ) '  R82VEC_BINNED_REORDER can reorder the array'
  write ( *, '(a)' ) '    to correspond to the bin ordering.'
  write ( *, '(a)' ) '  R82VEC_BINNED_SORT_A can sort the individual bins'
  write ( *, '(a)' ) '    after the array has been reordered.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The bins are equally spaced between '
  write ( *, '(a)' ) '  BIN_MIN and BIN_MAX,'
  write ( *, '(a)' ) '  with two extra bins, for things less than BIN_MIN,'
  write ( *, '(a)' ) '  or greater than BIN_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  Component 1 range: ', bin_min(1), bin_max(1)
  write ( *, '(a,2g14.6)' ) '  Component 2 range: ', bin_min(2), bin_max(2)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of bins per row and column = ', nbin
  write ( *, '(a)' ) ' '

  call get_seed ( seed )

  call r82vec_uniform ( n, amin, amax, seed, a )

  call r82vec_print ( n, a, '  The data vector A to be binned:' )

  call r82vec_bin_even ( n, a, nbin, bin_min, bin_max, bin_start, bin_last, &
    bin_next )

  call i4mat_print ( nbin, nbin, bin_start, '  The BIN_START array:' )

  call i4mat_print ( nbin, nbin, bin_start, '  The BIN_LAST array:' )

  call i4vec_print ( n, bin_next, '  The BIN_NEXT array:' )

  do i1 = 1, nbin

    do i2 = 1, nbin

      write ( *, '(a)' ) ' '
      write ( *, '(a,2i6)' ) '  Contents of bin number ', i1, i2
      write ( *, '(a)' ) ' '

      j = bin_start(i1,i2)
      k = 0

      do while ( 0 < j )
        k = k + 1
        write ( *, '(2x,2i4,2g14.6)' ) k, j, a(1,j), a(2,j)
        j = bin_next(j)
      end do

    end do

  end do
!
!  Now reorder the data to correspond to the bins.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call R82VEC_BINNED_REORDER to reorder the array.'
  write ( *, '(a)' ) ' '

  call r82vec_binned_reorder ( n, a, nbin, bin_start, bin_last, bin_next )

  call r82vec_print ( n, a, '  The data vector, sorted by bins:' )

  call i4mat_print ( nbin, nbin, bin_start, '  The BIN_START array:' )

  call i4mat_print ( nbin, nbin, bin_last, '  The BIN_LAST array:' )

  call i4vec_print ( n, bin_next, '  The BIN_NEXT array:' )
!
!  Now sort the bins.
!
  call r82vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )

  call r82vec_print ( n, a, '  The data vector, with sorted bins:' )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests R82VEC_PART_QUICK_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ), parameter, dimension ( 2 ) :: alo = (/  0.0D+00, 2.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 2 ) :: ahi = (/ 10.0D+00, 3.0D+00 /)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) r
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  R82VEC_PART_QUICK_A reorders an R82VEC'
  write ( *, '(a)' ) '    as part of a quick sort.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call r82vec_uniform ( n, alo, ahi, seed, a )

  call r82vec_print ( n, a, '  Before rearrangment:' )

  call r82vec_part_quick_a ( n, a, l, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rearranged array'
  write ( *, '(a,i6)' ) '  Left index =  ', l
  write ( *, '(a,i6)' ) '  Key index =   ', l+1
  write ( *, '(a,i6)' ) '  Right index = ', r
  write ( *, '(a)' ) ' '

  call r82vec_print ( l,     a(1:2,1:l),   '  Left half:' )
  call r82vec_print ( 1,     a(1:2,l+1),   '  Key:' )
  call r82vec_print ( n-l-1, a(1:2,l+2:n), '  Right half:' )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests R82VEC_SORT_QUICK_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ), parameter, dimension ( 2 ) :: alo = (/  0.0D+00, 2.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 2 ) :: ahi = (/ 10.0D+00, 3.0D+00 /)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  R82VEC_SORT_QUICK_A sorts an R82VEC'
  write ( *, '(a)' ) '    using quick sort.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call r82vec_uniform ( n, alo, ahi, seed, a )
!
!  For better testing, give a few elements the same first component.
!
  a(1,3) = a(1,5)
  a(1,4) = a(1,12)
!
!  Make two entries equal.
!
  a(1:2,7) = a(1:2,11)

  call r82vec_print ( n, a, '  Before rearrangment:' )

  call r82vec_sort_quick_a ( n, a )

  call r82vec_print ( n, a, '  Sorted array:' )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests R83VEC_PART_QUICK_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ), parameter, dimension ( 3 ) :: alo = (/ &
     0.0D+00, 2.0D+00, 1.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 3 ) :: ahi = (/ &
    10.0D+00, 3.0D+00, 3.0D+00 /)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) r
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  R83VEC_PART_QUICK_A reorders an R83VEC'
  write ( *, '(a)' ) '    as part of a quick sort.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call r83vec_uniform ( n, alo, ahi, seed, a )

  call r83vec_print ( n, a, '  Before rearrangment:' )

  call r83vec_part_quick_a ( n, a, l, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rearranged array'
  write ( *, '(a,i6)' ) '  Left index =  ', l
  write ( *, '(a,i6)' ) '  Key index =   ', l+1
  write ( *, '(a,i6)' ) '  Right index = ', r
  write ( *, '(a)' ) ' '

  call r83vec_print ( l,     a(1:3,1:l),   '  Left half:' )
  call r83vec_print ( 1,     a(1:3,l+1),   '  Key:' )
  call r83vec_print ( n-l-1, a(1:3,l+2:n), '  Right half:' )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests R83VEC_SORT_QUICK_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ), parameter, dimension ( 3 ) :: alo = (/  &
     0.0D+00, 2.0D+00, 1.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 3 ) :: ahi = (/ &
    10.0D+00, 3.0D+00, 3.0D+00 /)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  R83VEC_SORT_QUICK_A sorts an R83VEC'
  write ( *, '(a)' ) '    using quick sort.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call r83vec_uniform ( n, alo, ahi, seed, a )
!
!  Give two elements the same first component.
!
  a(1,3) = a(1,5)
!
!  Give two elements the same first two components.
!
  a(1:2,4) = a(1:2,12)
!
!  Make two entries equal.
!
  a(1:3,7) = a(1:3,11)

  call r83vec_print ( n, a, '  Before rearrangment:' )

  call r83vec_sort_quick_a ( n, a )

  call r83vec_print ( n, a, '  Sorted array:' )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests R8VEC_BIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 25
  integer ( kind = 4 ), parameter :: nbin = 5

  integer ( kind = 4 ) bin(0:nbin+1)
  real ( kind = 8 ) bin_limit(0:nbin)
  real ( kind = 8 ) bin_max
  real ( kind = 8 ) bin_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  R8VEC_BIN computes bins for an R8VEC.'

  call get_seed ( seed )

  call r8vec_uniform ( n, -2.0D+00, 11.0D+00, seed, x )

  call r8vec_print ( n, x, '  The vector to be binned:' )

  bin_min =  0.0D+00
  bin_max = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of bins is ', nbin
  write ( *, '(a,g14.6)' ) '  Bin minimum is ', bin_min
  write ( *, '(a,g14.6)' ) '  Bin maximum is ', bin_max

  call r8vec_bin ( n, x, nbin, bin_min, bin_max, bin, bin_limit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lower Limit    Upper Limit    Count'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,2f8.4,i4)' ) bin_min, bin_limit(0), bin(0)
  do i = 1, nbin
    write ( *, '(2x,2f8.4,i4)' ) bin_limit(i-1), bin_limit(i), bin(i)
  end do
  write ( *, '(2x,f8.4,8x,i4)' ) bin_limit(nbin), bin(nbin+1)

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests R8VEC_BIN_EVEN, R8VEC_BINNED_REORDER, R8VEC_BINNED_SORT_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 30
  integer ( kind = 4 ), parameter :: nbin = 7

  real ( kind = 8 ) a(n)
  real ( kind = 8 ), parameter :: amax = 23.0D+00
  real ( kind = 8 ), parameter :: amin = 8.0D+00
  real ( kind = 8 ), parameter :: bin_max = 20.0D+00
  real ( kind = 8 ), parameter :: bin_min = 10.0D+00
  integer ( kind = 4 ) bin_last(nbin)
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  R8VEC_BIN_EVEN constructs evenly spaced bins and'
  write ( *, '(a)' ) '    assigns each element of a DVEC to a bin.'
  write ( *, '(a)' ) '  R8VEC_BINNED_REORDER can reorder the array'
  write ( *, '(a)' ) '    to correspond to the bin ordering.'
  write ( *, '(a)' ) '  R8VEC_BINNED_SORT_A can sort the array'
  write ( *, '(a)' ) '    once it has been reordered.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The bins are equally spaced between '
  write ( *, '(a)' ) '  BIN_MIN and BIN_MAX,'
  write ( *, '(a)' ) '  with two extra bins, for things less than BIN_MIN,'
  write ( *, '(a)' ) '  or greater than BIN_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  BIN_MIN = ', bin_min
  write ( *, '(a,g14.6)' ) '  BIN_MAX = ', bin_max
  write ( *, '(a,i6)' ) '  Total number of bins = ', nbin
  write ( *, '(a)' ) ' '

  call get_seed ( seed )

  call r8vec_uniform ( n, amin, amax, seed, a )

  call r8vec_print ( n, a, '  The data vector A to be binned:' )

  call r8vec_bin_even ( n, a, nbin, bin_min, bin_max, bin_start, &
    bin_last, bin_next )

  call i4vec_print ( nbin, bin_start, '  The BIN_START array:' )

  call i4vec_print ( nbin, bin_last, '  The BIN_LAST array:' )

  call i4vec_print ( n, bin_next, '  The BIN_NEXT array:' )

  do i = 1, nbin

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Contents of bin number ', i
    write ( *, '(a)' ) ' '

    j = bin_start(i)
    k = 0

    do while ( 0 < j )
      k = k + 1
      write ( *, '(2x,2i4,g14.6)' ) k, j, a(j)
      j = bin_next(j)
    end do

  end do
!
!  Now reorder the data to correspond to the bins.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call R8VEC_BINNED_REORDER to reorder the array.'
  write ( *, '(a)' ) ' '

  call r8vec_binned_reorder ( n, a, nbin, bin_start, bin_last, bin_next )

  call r8vec_print ( n, a, '  The data vector A:' )

  call i4vec_print ( nbin, bin_start, '  The BIN_START array:' )

  call i4vec_print ( nbin, bin_last, '  The BIN_LAST array:' )

  call i4vec_print ( n, bin_next, '  The BIN_NEXT array:' )
!
!  Now sort the data, one bin at a time
!
  call r8vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )

  call r8vec_print ( n, a, '  The sorted data vector A:' )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests POINTS_NEAREST_POINT_BINS_2D, POINTS_NEAREST_POINT_BINS_2D_2,
!    POINTS_NEAREST_POINT_BINS_2D_3, POINTS_NEAREST_POINT_NAIVE_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  integer ( kind = 4 ), parameter :: nbin = 10
  integer ( kind = 4 ), parameter, dimension ( ndim ) :: nbin2 = (/ 20, 5 /)
  integer ( kind = 4 ), parameter :: nset = 1000
  integer ( kind = 4 ), parameter :: ntest = 10

  real ( kind = 8 ), parameter, dimension ( ndim ) :: bin_min = (/  &
    0.0D+00,  0.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( ndim ) :: bin_max = (/  &
    20.0D+00, 5.0D+00 /)
  integer ( kind = 4 ) bin_last(nbin,nbin)
  integer ( kind = 4 ) bin_last2(nbin2(1),nbin2(2))
  integer ( kind = 4 ) bin_next(nset)
  integer ( kind = 4 ) bin_next2(nset)
  integer ( kind = 4 ) bin_start(nbin,nbin)
  integer ( kind = 4 ) bin_start2(nbin2(1),nbin2(2))
  integer ( kind = 4 ) compares
  real ( kind = 8 ) d_min
  logical :: debug = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min
  real ( kind = 8 ) p(ndim)
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) pset2(ndim,nset)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  Given a point in 2D, we want to find its nearest'
  write ( *, '(a)' ) '  neighbor among points in a set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POINTS_NEAREST_POINT_NAIVE_2D uses a naive algorithm.'
  write ( *, '(a)' ) '  POINTS_NEAREST_POINT_BINS_2D and'
  write ( *, '(a)' ) '  POINTS_NEAREST_POINT_BINS_2D_2 use bins, but require'
  write ( *, '(a)' ) '    the same number in each direction.'
  write ( *, '(a)' ) '  POINTS_NEAREST_POINT_BINS_2D_3 uses bins, and can use'
  write ( *, '(a)' ) '    a different number in each direction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points in the pointset is ', nset
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POINTS_NEAREST_POINT_BINS_2D and'
  write ( *, '(a,i6)' ) '  POINTS_NEAREST_POINT_BINS_2D_2 use ', nbin
  write ( *, '(a)' ) '    bins in each direction.'
  write ( *, '(a,2i6)' ) '  POINTS_NEAREST_POINT_BINS_2D_3 uses ', nbin2(1:ndim)
  write ( *, '(a)' ) '    bins in each direction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  The X coordinate range: ', bin_min(1), bin_max(1)
  write ( *, '(a,2g14.6)' ) '  The Y coordinate range: ', bin_min(2), bin_max(2)
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  Test point X range:     ', bin_min(1), bin_max(1)
  write ( *, '(a,2g14.6)' ) '  Test point Y range:     ', bin_min(2), bin_max(2)
!
!  Set the pointset.
!
  call r82vec_uniform ( nset, bin_min, bin_max, seed, pset )
!
!  We need to make a copy of the point set, because it gets sorted.
!
  pset2(1:ndim,1:nset) = pset(1:ndim,1:nset)
!
!  For the POINTS_NEAREST_POINT_BINS_2D code:
!
!    Implicitly bin the data
!    Explicitly reorder the data by bins.
!    Within each bin, sort the data.
!
  call r82vec_bin_even2 ( nset, pset, nbin, bin_min, bin_max, bin_start, &
    bin_last, bin_next )

  call r82vec_binned_reorder ( nset, pset, nbin, bin_start, bin_last, bin_next )

  call r82vec_binned_sort_a ( nset, pset, nbin, bin_start, bin_last )
!
!  For the POINTS_NEAREST_POINT_BINS_2D_3 code:
!
!    Implicitly bin the data
!    Explicitly reorder the data by bins.
!    Within each bin, sort the data.
!
  call r82vec_bin_even3 ( nset, pset2, nbin2, bin_min, bin_max, bin_start2, &
    bin_last2, bin_next2 )

  call r82vec_binned_reorder2 ( nset, pset2, nbin2, bin_start2, bin_last2, &
    bin_next2 )

  call r82vec_binned_sort_a2 ( nset, pset2, nbin2, bin_start2, bin_last2 )
!
!  Seek nearest neighbors.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Test point           Neighbor point      Distance'
  write ( *, '(a)' ) '--------------------  --------------------  ----------'
  write ( *, '(a)' ) ' '

  do i = 1, ntest

    call r82_uniform ( bin_min, bin_max, seed, p )

    write ( *, '(a)' ) ' '

    call points_nearest_point_naive_2d ( nset, pset, p, i_min, d_min )

    compares = nset
    write ( *, '(2x,2f10.4,2x,2f10.4,2x,f10.4,2x,i4)' ) p(1:ndim), &
      pset(1:ndim,i_min), d_min, compares

    call points_nearest_point_bins_2d ( nset, pset, nbin, bin_min, bin_max, &
      bin_start, bin_last, bin_next, p, i_min, d_min, compares )

    write ( *, '(2x,2f10.4,2x,2f10.4,2x,f10.4,2x,i4)' ) p(1:ndim), &
      pset(1:ndim,i_min), d_min, compares

    call points_nearest_point_bins_2d_2 ( nset, pset, nbin, bin_min, bin_max, &
      bin_start, bin_last, bin_next, p, i_min, d_min, compares )

    write ( *, '(2x,2f10.4,2x,2f10.4,2x,f10.4,2x,i4)' ) p(1:ndim), &
      pset(1:ndim,i_min), d_min, compares

    call points_nearest_point_bins_2d_3 ( nset, pset2, nbin2, bin_min, &
      bin_max, bin_start2, bin_last2, bin_next2, p, i_min, d_min, compares )

    write ( *, '(2x,2f10.4,2x,2f10.4,2x,f10.4,2x,i4)' ) p(1:ndim), &
      pset(1:ndim,i_min), d_min, compares

  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests POINTS_NEAREST_POINTS_BINS_2D, POINTS_NEAREST_POINTS_BINS_2D_2,
!    POINTS_NEAREST_POINTS_BINS_2D_3, POINTS_NEAREST_POINTS_NAIVE_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  integer ( kind = 4 ), parameter :: nbin = 10
  integer ( kind = 4 ), parameter, dimension ( ndim ) :: nbin2 = (/ 10, 10 /)
  integer ( kind = 4 ), parameter :: nset = 1000
  integer ( kind = 4 ), parameter :: ntest = 100

  real ( kind = 8 ), parameter, dimension ( ndim ) :: bin_min = (/ &
     0.0D+00,  0.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( ndim ) :: bin_max = (/ &
    10.0D+00,  10.0D+00 /)
  integer ( kind = 4 ) bin_last(nbin,nbin)
  integer ( kind = 4 ) bin_last2(nbin2(1),nbin2(2))
  integer ( kind = 4 ) bin_next(nset)
  integer ( kind = 4 ) bin_next2(nset)
  integer ( kind = 4 ) bin_start(nbin,nbin)
  integer ( kind = 4 ) bin_start2(nbin2(1),nbin2(2))
  integer ( kind = 4 ) clock_count1
  integer ( kind = 4 ) clock_count2
  integer ( kind = 4 ) clock_count3
  integer ( kind = 4 ) clock_count4
  integer ( kind = 4 ) clock_count5
  integer ( kind = 4 ) clock_count6
  integer ( kind = 4 ) clock_max
  integer ( kind = 4 ) clock_rate
  integer ( kind = 4 ) compares(ntest)
  real ( kind = 8 ) d_min0(ntest)
  real ( kind = 8 ) d_min1(ntest)
  real ( kind = 8 ) d_min2(ntest)
  real ( kind = 8 ) d_min3(ntest)
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min0(ntest)
  integer ( kind = 4 ) i_min1(ntest)
  integer ( kind = 4 ) i_min2(ntest)
  integer ( kind = 4 ) i_min3(ntest)
  integer ( kind = 4 ) n_different
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim,ntest)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  Given a point set in 2D, and a set of test points,'
  write ( *, '(a)' ) '  for each testpoint, find the nearest neighbor in'
  write ( *, '(a)' ) '  the point set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POINTS_NEAREST_POINTS_NAIVE_2D uses a naive algorithm.'
  write ( *, '(a)' ) '  POINTS_NEAREST_POINTS_BINS_2D uses equal bins.'
  write ( *, '(a)' ) '  POINTS_NEAREST_POINTS_BINS_2D_2 uses equal bins.'
  write ( *, '(a)' ) '  POINTS_NEAREST_POINTS_BINS_2D_3 uses variable bins.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points in the pointset is ', nset
  write ( *, '(a,i6)' ) '  The number of points in the test set is ', ntest
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POINTS_NEAREST_POINTS_BINS_2D and '
  write ( *, '(a,i6)' ) '  POINTS_NEAREST_POINTS_BINS_2D_2 use ', nbin
  write ( *, '(a)' ) '    bins in each direction.'
  write ( *, '(a,2i6)' ) '  POINTS_NEAREST_POINTS_BINS_2D_3 uses ', &
    nbin2(1:ndim)
  write ( *, '(a)' ) '    bins in each direction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  The X coordinate range: ', bin_min(1), bin_max(1)
  write ( *, '(a,2g14.6)' ) '  The Y coordinate range: ', bin_min(2), bin_max(2)
  write ( *, '(a)' ) ' '
!
!  Set the pointset.
!
  call r82vec_uniform ( nset, bin_min, bin_max, seed, pset )
!
!  Set the test points.
!
  call r82vec_uniform ( ntest, bin_min, bin_max, seed, ptest )
!
!  For POINTS_NEAREST_POINTS_BINS_2D and POINTS_NEAREST_POINTS_BINS_2D_2:
!
!  Implicitly bin the data.
!  Explicitly reorder the data by bins.
!  Within each bin, sort the data.
!
  call r82vec_bin_even2 ( nset, pset, nbin, bin_min, bin_max, bin_start, &
    bin_last, bin_next )

  call r82vec_binned_reorder ( nset, pset, nbin, bin_start, bin_last, bin_next )

  call r82vec_binned_sort_a ( nset, pset, nbin, bin_start, bin_last )
!
!  Seek nearest neighbors.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Print results for up to first 10 points...'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Test point		      Distance'
  write ( *, '(a)' ) '                       Naive     Bins     Bins2     Bins3'
  write ( *, '(a)' ) &
    '--------------------  ------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '

  call system_clock ( clock_count1, clock_rate, clock_max )

  call points_nearest_points_naive_2d ( nset, pset, ntest, ptest, i_min0, &
    d_min0 )

  call system_clock ( clock_count2, clock_rate, clock_max )

  call points_nearest_points_bins_2d ( nset, pset, nbin, bin_min, bin_max, &
    bin_start, bin_last, bin_next, ntest, ptest, i_min1, d_min1, compares )

  call system_clock ( clock_count3, clock_rate, clock_max )

  call points_nearest_points_bins_2d_2 ( nset, pset, nbin, bin_min, bin_max, &
    bin_start, bin_last, bin_next, ntest, ptest, i_min2, d_min2, compares )

  call system_clock ( clock_count4, clock_rate, clock_max )
!
!  We have to rework the data for BINS3, since we allow a different
!  number of bins in each direction.
!
  call r82vec_bin_even3 ( nset, pset, nbin2, bin_min, bin_max, bin_start2, &
    bin_last2, bin_next2 )

  call r82vec_binned_reorder ( nset, pset, nbin2, bin_start2, bin_last2, &
    bin_next2 )

  call r82vec_binned_sort_a ( nset, pset, nbin2, bin_start2, bin_last2 )

  call system_clock ( clock_count5, clock_rate, clock_max )

  call points_nearest_points_bins_2d_3 ( nset, pset, nbin2, bin_min, bin_max, &
    bin_start2, bin_last2, bin_next2, ntest, ptest, i_min3, d_min3, compares )

  call system_clock ( clock_count6, clock_rate, clock_max )
!
!  Print the results.
!
  do i = 1, min ( ntest, 10 )
    write ( *, '(2x,2f10.4,4x,4f10.4)' ) ptest(1:ndim,i), d_min0(i), &
      d_min1(i), d_min2(i), d_min3(i)
  end do
!
!  Check the results.
!
  eps = epsilon ( eps )

  n_different = 0
  do i = 1, ntest
    if ( eps < abs ( d_min0(i) - d_min1(i) ) ) then
      n_different = n_different + 1
    end if
  end do

  if ( n_different == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Naive and bin1 codes computed the same results.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WARNING:'
    write ( *, '(a)' ) '  Naive and bin1 codes disagreed.'
    write ( *, '(a,i6)' ) '  Number of discrepancies was ', n_different
  end if

  n_different = 0
  do i = 1, ntest
    if ( eps < abs ( d_min0(i) - d_min2(i) ) ) then
      n_different = n_different + 1
    end if
  end do

  if ( n_different == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Naive and bin2 codes computed the same results.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WARNING:'
    write ( *, '(a)' ) '  Naive and bin2 codes disagreed.'
    write ( *, '(a,i6)' ) '  Number of discrepancies was ', n_different
  end if

  n_different = 0
  do i = 1, ntest
    if ( eps < abs ( d_min0(i) - d_min3(i) ) ) then
      n_different = n_different + 1
    end if
  end do

  if ( n_different == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Naive and bin3 codes computed the same results.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WARNING:'
    write ( *, '(a)' ) '  Naive and bin3 codes disagreed.'
    write ( *, '(a,i6)' ) '  Number of discrepancies was ', n_different
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Naive code time = ', &
    real ( clock_count2 - clock_count1, kind = 8 ) &
    / real ( clock_rate, kind = 8 )

  write ( *, '(a,g14.6)' ) '  Bin code time =	', &
    real ( clock_count3 - clock_count2, kind = 8 ) &
    / real ( clock_rate, kind = 8 )

  write ( *, '(a,g14.6)' ) '  Bin2 code time =   ', &
    real ( clock_count4 - clock_count3, kind = 8 ) &
    / real ( clock_rate, kind = 8 )

  write ( *, '(a,g14.6)' ) '  Bin3 code time =   ', &
    real ( clock_count6 - clock_count5, kind = 8 ) &
    / real ( clock_rate, kind = 8 )

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests POINTS_NEAREST_POINTS_BINS_2D, POINTS_NEAREST_POINTS_BINS_2D_2,
!    POINTS_NEAREST_POINTS_BINS_2D_3, POINTS_NEAREST_POINTS_NAIVE_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  integer ( kind = 4 ), parameter :: nbin = 10
  integer ( kind = 4 ), parameter, dimension ( ndim ) :: nbin2 = (/ 4, 25 /)
  integer ( kind = 4 ), parameter :: nset = 1000
  integer ( kind = 4 ), parameter :: ntest = 100

  real ( kind = 8 ), parameter, dimension ( ndim ) :: bin_min = (/  &
    0.0D+00,  0.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( ndim ) :: bin_max = (/ &
    4.0D+00, 25.0D+00 /)
  integer ( kind = 4 ) bin_last(nbin,nbin)
  integer ( kind = 4 ) bin_last2(nbin2(1),nbin2(2))
  integer ( kind = 4 ) bin_next(nset)
  integer ( kind = 4 ) bin_next2(nset)
  integer ( kind = 4 ) bin_start(nbin,nbin)
  integer ( kind = 4 ) bin_start2(nbin2(1),nbin2(2))
  integer ( kind = 4 ) clock_count1
  integer ( kind = 4 ) clock_count2
  integer ( kind = 4 ) clock_count3
  integer ( kind = 4 ) clock_count4
  integer ( kind = 4 ) clock_count5
  integer ( kind = 4 ) clock_count6
  integer ( kind = 4 ) clock_max
  integer ( kind = 4 ) clock_rate
  integer ( kind = 4 ) compares(ntest)
  real ( kind = 8 ) d_min0(ntest)
  real ( kind = 8 ) d_min1(ntest)
  real ( kind = 8 ) d_min2(ntest)
  real ( kind = 8 ) d_min3(ntest)
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min0(ntest)
  integer ( kind = 4 ) i_min1(ntest)
  integer ( kind = 4 ) i_min2(ntest)
  integer ( kind = 4 ) i_min3(ntest)
  integer ( kind = 4 ) n_different
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim,ntest)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  Given a point set in 2D, and a set of test points,'
  write ( *, '(a)' ) '  for each testpoint, find the nearest neighbor in'
  write ( *, '(a)' ) '  the point set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, the region is RECTANGULAR.'
  write ( *, '(a)' ) &
    '  The BINS and BINS2 codes will end up using rectangular bins;'
  write ( *, '(a)' ) &
    '  We will set the BINS3 code to use the same number of bins,'
  write ( *, '(a)' ) '  but they will be square.  This should mean that BINS3'
  write ( *, '(a)' ) '  finds a match faster.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POINTS_NEAREST_POINTS_NAIVE_2D uses a naive algorithm.'
  write ( *, '(a)' ) '  POINTS_NEAREST_POINTS_BINS_2D uses bins.'
  write ( *, '(a)' ) '  POINTS_NEAREST_POINTS_BINS_2D_2 uses bins.'
  write ( *, '(a)' ) '  POINTS_NEAREST_POINTS_BINS_2D_3 uses bins.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points in the pointset is ', nset
  write ( *, '(a,i6)' ) '  The number of bins used in each direction is ', nbin
  write ( *, '(a,i6)' ) '  The number of points in the test set is ', ntest
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  The X coordinate range: ', bin_min(1), bin_max(1)
  write ( *, '(a,2g14.6)' ) '  The Y coordinate range: ', bin_min(2), bin_max(2)
  write ( *, '(a)' ) ' '
!
!  Set the pointset.
!
  call r82vec_uniform ( nset, bin_min, bin_max, seed, pset )
!
!  Set the test points.
!
  call r82vec_uniform ( ntest, bin_min, bin_max, seed, ptest )
!
!  Implicitly bin the data.
!
  call r82vec_bin_even2 ( nset, pset, nbin, bin_min, bin_max, bin_start, &
    bin_last, bin_next )
!
!  Explicitly reorder the data by bins.
!
  call r82vec_binned_reorder ( nset, pset, nbin, bin_start, bin_last, bin_next )
!
!  Within each bin, sort the data.
!
  call r82vec_binned_sort_a ( nset, pset, nbin, bin_start, bin_last )
!
!  Seek nearest neighbors.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Print results for up to first 10 points...'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Test point		      Distance'
  write ( *, '(a)' ) '                       Naive     Bins     Bins2     Bins3'
  write ( *, '(a)' ) &
    '--------------------  ------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '

  call system_clock ( clock_count1, clock_rate, clock_max )

  call points_nearest_points_naive_2d ( nset, pset, ntest, ptest, i_min0, &
    d_min0 )

  call system_clock ( clock_count2, clock_rate, clock_max )

  call points_nearest_points_bins_2d ( nset, pset, nbin, bin_min, bin_max, &
    bin_start, bin_last, bin_next, ntest, ptest, i_min1, d_min1, compares )

  call system_clock ( clock_count3, clock_rate, clock_max )

  call points_nearest_points_bins_2d_2 ( nset, pset, nbin, bin_min, bin_max, &
    bin_start, bin_last, bin_next, ntest, ptest, i_min2, d_min2, compares )

  call system_clock ( clock_count4, clock_rate, clock_max )
!
!  We have to rework the data for BINS3, since we allow a different
!  number of bins in each direction.
!
  call r82vec_bin_even3 ( nset, pset, nbin2, bin_min, bin_max, bin_start2, &
    bin_last2, bin_next2 )

  call r82vec_binned_reorder2 ( nset, pset, nbin2, bin_start2, bin_last2, &
    bin_next2 )

  call r82vec_binned_sort_a2 ( nset, pset, nbin2, bin_start2, bin_last2 )

  call system_clock ( clock_count5, clock_rate, clock_max )

  call points_nearest_points_bins_2d_3 ( nset, pset, nbin2, bin_min, bin_max, &
    bin_start2, bin_last2, bin_next2, ntest, ptest, i_min3, d_min3, compares )

  call system_clock ( clock_count6, clock_rate, clock_max )
!
!  Print the results.
!
  do i = 1, min ( ntest, 10 )
    write ( *, '(2x,2f10.4,4x,4f10.4)' ) ptest(1:ndim,i), d_min0(i), &
      d_min1(i), d_min2(i), d_min3(i)
  end do
!
!  Check the results.
!
  eps = epsilon ( eps )

  n_different = 0
  do i = 1, ntest
    if ( eps < abs ( d_min0(i) - d_min1(i) ) ) then
      n_different = n_different + 1
    end if
  end do

  if ( n_different == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Naive and bin1 codes computed the same results.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WARNING:'
    write ( *, '(a)' ) '  Naive and bin1 codes disagreed.'
    write ( *, '(a,i6)' ) '  Number of discrepancies was ', n_different
  end if

  n_different = 0
  do i = 1, ntest
    if ( eps < abs ( d_min0(i) - d_min2(i) ) ) then
      n_different = n_different + 1
    end if
  end do

  if ( n_different == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Naive and bin2 codes computed the same results.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WARNING:'
    write ( *, '(a)' ) '  Naive and bin2 codes disagreed.'
    write ( *, '(a,i6)' ) '  Number of discrepancies was ', n_different
  end if

  n_different = 0
  do i = 1, ntest
    if ( eps < abs ( d_min0(i) - d_min3(i) ) ) then
      n_different = n_different + 1
    end if
  end do

  if ( n_different == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Naive and bin3 codes computed the same results.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WARNING:'
    write ( *, '(a)' ) '  Naive and bin3 codes disagreed.'
    write ( *, '(a,i6)' ) '  Number of discrepancies was ', n_different
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Naive code time = ', &
    real ( clock_count2 - clock_count1, kind = 8 ) &
    / real ( clock_rate, kind = 8 )

  write ( *, '(a,g14.6)' ) '  Bin code time =   ', &
    real ( clock_count3 - clock_count2, kind = 8 ) &
    / real ( clock_rate, kind = 8 )

  write ( *, '(a,g14.6)' ) '  Bin2 code time =   ', &
    real ( clock_count4 - clock_count3, kind = 8 ) &
    / real ( clock_rate, kind = 8 )

  write ( *, '(a,g14.6)' ) '  Bin3 code time =   ', &
    real ( clock_count6 - clock_count5, kind = 8 ) &
    / real ( clock_rate, kind = 8 )

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests POINTS_NEAREST_POINTS_BINS_3D_2, POINTS_NEAREST_POINTS_NAIVE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nbin = 32
  integer ( kind = 4 ), parameter :: ndim = 3
  integer ( kind = 4 ), parameter :: nset = 4096
  integer ( kind = 4 ), parameter :: ntest = 1000

  real ( kind = 8 ), parameter, dimension ( ndim ) :: bin_min = (/  &
    0.0D+00,  0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( ndim ) :: bin_max = (/ &
    10.0D+00, 10.0D+00, 10.0D+00 /)
  integer ( kind = 4 ) bin_last(nbin,nbin,nbin)
  integer ( kind = 4 ) bin_next(nset)
  integer ( kind = 4 ) bin_start(nbin,nbin,nbin)
  integer ( kind = 4 ) clock_count1
  integer ( kind = 4 ) clock_count2
  integer ( kind = 4 ) clock_count3
  integer ( kind = 4 ) clock_max
  integer ( kind = 4 ) clock_rate
  integer ( kind = 4 ) compares(ntest)
  real ( kind = 8 ) d_min1(ntest)
  real ( kind = 8 ) d_min2(ntest)
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min1(ntest)
  integer ( kind = 4 ) i_min2(ntest)
  integer ( kind = 4 ) n_different
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim,ntest)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  Given a point set in 3D, and a set of test points,'
  write ( *, '(a)' ) '  for each testpoint, find the nearest neighbor in'
  write ( *, '(a)' ) '  the point set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POINTS_NEAREST_POINTS_NAIVE_3D uses a naive algorithm.'
  write ( *, '(a)' ) '  POINTS_NEAREST_POINTS_BINS_3D_2 uses bins.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points in the pointset is ', nset
  write ( *, '(a,i6)' ) '  The number of bins used in each direction is ', nbin
  write ( *, '(a,i6)' ) '  The number of points in the test set is ', ntest
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  The X coordinate range: ', bin_min(1), bin_max(1)
  write ( *, '(a,2g14.6)' ) '  The Y coordinate range: ', bin_min(2), bin_max(2)
  write ( *, '(a,2g14.6)' ) '  The Z coordinate range: ', bin_min(3), bin_max(3)
  write ( *, '(a)' ) ' '
!
!  Set the pointset.
!
  call r83vec_uniform ( nset, bin_min, bin_max, seed, pset )
!
!  Set the test points.
!
  call r83vec_uniform ( ntest, bin_min, bin_max, seed, ptest )
!
!  Implicitly bin the data.
!
  call r83vec_bin_even2 ( nset, pset, nbin, bin_min, bin_max, bin_start, &
    bin_last, bin_next )
!
!  Explicitly reorder the data by bins.
!
  call r83vec_binned_reorder ( nset, pset, nbin, bin_start, bin_last, bin_next )
!
!  Within each bin, sort the data.
!
  call r83vec_binned_sort_a ( nset, pset, nbin, bin_start, bin_last )
!
!  Seek nearest neighbors.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Print up to the first 10 points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '    Test point                       Distance        Comparisons'
  write ( *, '(a)' ) &
    '                                 Naive     Bins     Naive Bins'
  write ( *, '(a)' ) &
    '-----------------------------  --------------------  ----------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '

  call system_clock ( clock_count1, clock_rate, clock_max )

  call points_nearest_points_naive_3d ( nset, pset, ntest, ptest, i_min1, &
    d_min1 )

  call system_clock ( clock_count2, clock_rate, clock_max )

  call points_nearest_points_bins_3d_2 ( nset, pset, nbin, bin_min, bin_max, &
    bin_start, bin_last, bin_next, ntest, ptest, i_min2, d_min2, compares )

  call system_clock ( clock_count3, clock_rate, clock_max )

  do i = 1, min ( ntest, 10 )
    write ( *, '(5f10.4,2x,2i6)' ) ptest(1:ndim,i), d_min1(i), d_min2(i), &
      nset, compares(i)
  end do

  eps = epsilon ( eps )
  n_different = 0
  do i = 1, ntest
    if ( eps < abs ( d_min1(i) - d_min2(i) ) ) then
      n_different = n_different + 1
    end if
  end do

  if ( n_different == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Naive and bin codes computed the same results.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WARNING:'
    write ( *, '(a)' ) '  Naive and bin codes disagreed.'
    write ( *, '(a,i6)' ) '  Number of discrepancies was ', n_different
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Naive code time = ', &
    real ( clock_count2 - clock_count1, kind = 8 ) &
    / real ( clock_rate, kind = 8 )

  write ( *, '(a,g14.6)' ) '  Bin code time =   ', &
    real ( clock_count3 - clock_count2, kind = 8 ) &
    / real ( clock_rate, kind = 8 )

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests POINTS_NEAREST_POINT_NAIVE_2D, POINTS_NEAREST_POINT_DEL_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nabes_max = 500
  integer ( kind = 4 ), parameter :: num_pts = 13
  integer ( kind = 4 ), parameter :: ntest = 10

  real ( kind = 8 ) dnear1
  real ( kind = 8 ) dnear2
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) nabes(nabes_max)
  integer ( kind = 4 ) nabes_dim
  integer ( kind = 4 ) nabes_first(num_pts)
  integer ( kind = 4 ) nabes_num(num_pts)
  integer ( kind = 4 ) nnear1
  integer ( kind = 4 ) nnear2
  integer ( kind = 4 ) nod_tri(3,3*num_pts)
  integer ( kind = 4 ) num_tri
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) td(ntest)
  integer ( kind = 4 ) tnbr(3,3*num_pts)
  real ( kind = 8 ), dimension (2,num_pts) :: xc = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       2.0D+00, 2.0D+00, &
      -1.0D+00, 3.0D+00, &
      -2.0D+00, 2.0D+00, &
       8.0D+00, 2.0D+00, &
       9.0D+00, 5.0D+00, &
       7.0D+00, 4.0D+00, &
       5.0D+00, 6.0D+00, &
       6.0D+00, 7.0D+00, &
       8.0D+00, 8.0D+00, &
      11.0D+00, 7.0D+00, &
      10.0D+00, 4.0D+00, &
       6.0D+00, 4.0D+00 /), (/ 2, num_pts /) )
  real ( kind = 8 ) xd(2,ntest)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  Given a point set XC, and a single point XD,'
  write ( *, '(a)' ) '  find the nearest point in XC to XD.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POINTS_NEAREST_POINT_NAIVE_2D uses a naive method.'
  write ( *, '(a)' ) '  POINTS_NEAREST_POINT_DEL_2D uses the Delaunay'
  write ( *, '(a)' ) '  triangulation'
  write ( *, '(a)' ) '  TRIANGULATION_PRINT prints a triangulation.'
!
!  Set up the Delaunay triangulation.
!
  call dtris2 ( num_pts, xc, num_tri, nod_tri, tnbr )

  call triangulation_print ( num_pts, num_tri, xc, nod_tri, tnbr )
!
!  Determine the node neigbhor array.
!
  call triangulation_nabe_nodes ( num_pts, num_tri, nod_tri, nabes_first, &
    nabes_num, nabes_max, nabes_dim, nabes )

! call triangulation_nabe_nodes_print ( num_pts, nabes_first, &
!   nabes_num, nabes_dim, nabes )
!
!  Get the test points.
!
  write ( *, * ) 'DEBUG: About to call triangulation_sample.'
  call triangulation_sample ( num_pts, xc, num_tri, nod_tri, &
    ntest, seed, xd, td )
  write ( *, * ) 'DEBUG: Returned from triangulation_sample.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X        Y     Distance  Index'
  write ( *, '(a)' ) ' '

  do itest = 1, ntest

    call points_nearest_point_naive_2d ( num_pts, xc, xd(1,itest), nnear1, &
      dnear1 )

    call points_nearest_point_del_2d ( num_pts, xc, xd(1,itest), &
      nabes_first, nabes_num, nabes_dim, nabes, nnear2, dnear2 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,2f8.4   )' ) '  XD       ', xd(1:2,itest)
    write ( *, '(a,3f8.4,i6)' ) '  Naive    ', xc(1:2,nnear1), dnear1, nnear1
    write ( *, '(a,3f8.4,i6)' ) '  Delaunay ', xc(1:2,nnear2), dnear2, nnear2

  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests POINTS_NEAREST_POINT_NAIVE_ND.
!
!  !....3&11....
!  !............
!  !............
!  X..9.........
!  !.....5......
!  !...........6
!  !.4.2...10...
!  !.....8...12.
!  V............
!  !..7.........
!  !......1.....
!  !............
!  !......*.....
!  !----V----X--
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ), parameter :: ndim = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) d_min
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) i_test
  real ( kind = 8 ), dimension ( ndim, n ) :: x = &
    reshape ( &
      (/  7.0D+00,  3.0D+00, &
          4.0D+00,  7.0D+00, &
          5.0D+00, 13.0D+00, &
          2.0D+00,  7.0D+00, &
          6.0D+00,  9.0D+00, &
         12.0D+00,  8.0D+00, &
          3.0D+00,  4.0D+00, &
          6.0D+00,  6.0D+00, &
          3.0D+00, 10.0D+00, &
          8.0D+00,  7.0D+00, &
          5.0D+00, 13.0D+00, &
         10.0D+00,  6.0D+00 /), (/ ndim, n /) )
  real ( kind = 8 ), dimension ( ndim, test_num ) :: x_test = &
    reshape ( &
      (/  7.0D+00,  1.0D+00, &
          4.0D+00,  7.0D+00, &
          8.0D+00, 11.0D+00 /), (/ ndim, test_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) &
    '  POINTS_NEAREST_POINT_NAIVE_ND computes the nearest point'
  write ( *, '(a)' ) '    in a set of points, to a given point, in ND.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The spatial dimension NDIM is ', ndim
  write ( *, '(a,i6)' ) '  The number of points N is ', n

  call r8mat_print ( ndim, n, x, '  The set of points:' )

  do i_test = 1, test_num

    call points_nearest_point_naive_nd ( ndim, n, x, x_test(1,i_test), &
      i_min, d_min )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Test point is '
    write ( *, '(2x,3g14.6)' ) x_test(1:ndim,i_test)
    write ( *, '(a)' ) '  Nearest point is '
    write ( *, '(2x,3g14.6)' ) x(1:ndim,i_min)
    write ( *, '(a,g14.6)' ) '  Distance is ', d_min

  end do

  return
end
