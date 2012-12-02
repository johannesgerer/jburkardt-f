program main

!*****************************************************************************80
!
!! MAIN is the main program for F90_INTRINSICS.
!
!  Discussion:
!
!    F90_INTRINSICS tests FORTRAN intrinsic routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  interface test_present
    subroutine test_present ( a, b, ethel, fred, george )
      integer a
      integer b
      integer, optional :: ethel
      integer, optional :: fred
      integer, optional :: george
    end subroutine test_present
  end interface

  integer arg_a
  integer arg_b
  integer fred

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'F90_INTRINSICS'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FORTRAN90 intrinsic library.'
 
  call test_abs ( )
  call test_achar ( )
  call test_acos ( )
  call test_adjustl ( )
  call test_adjustr ( )
  call test_aimag ( )
  call test_aint ( )
  call test_all ( )
  call test_allocated ( )
  call test_anint ( )
  call test_any ( )
  call test_asin ( )
  call test_associated ( )
  call test_atan ( )
  call test_atan2 ( )
  call test_bit_size ( )
  call test_btest ( )
  call test_ceiling ( )
  call test_char ( )
  call test_cmplx ( )
  call test_conjg ( )
  call test_cos_c4 ( )
  call test_cos_r8 ( )
  call test_cosh ( )
  call test_count ( )
  call test_cshift ( )
  call test_date_and_time ( )
  call test_dble ( )
  call test_digits ( )
  call test_dim ( )
  call test_dot_product ( )
  call test_dprod ( )
  call test_eoshift ( )
  call test_epsilon ( )
  call test_exp ( )
  call test_exponent ( )
  call test_floor ( )
  call test_fraction ( )
  call test_huge ( )
  call test_iachar ( )
  call test_iand_i4 ( )
  call test_iand_i8 ( )
  call test_ibclr ( )
  call test_ibits ( )
  call test_ibset ( )
  call test_ichar ( )
  call test_ieor_i4 ( )
  call test_ieor_i8 ( )
  call test_index ( )
  call test_int ( )
  call test_ior_i4 ( )
  call test_ior_i8 ( )
  call test_ishft ( )
  call test_ishftc ( )
  call test_kind ( )
  call test_lbound ( )
  call test_len ( )
  call test_len_trim ( )
  call test_lge ( )
  call test_lgt ( )
  call test_lle ( )
  call test_llt ( )
  call test_log ( )
  call test_log10 ( )
  call test_logical ( )
  call test_matmul ( )
  call test_max ( )
  call test_max_vector ( )
  call test_maxexponent ( )
  call test_maxloc ( )
  call test_maxval ( )
  call test_merge ( )
  call test_min ( )
  call test_minexponent ( )
  call test_minloc ( )
  call test_minval ( )
  call test_mod_i4 ( )
  call test_mod_r4 ( )
  call test_modulo_i4 ( )
  call test_modulo_r4 ( )
  call test_mvbits ( )
  call test_nearest ( )
  call test_nint ( )
  call test_not_i4 ( )
  call test_not_i8 ( )
  call test_pack ( )
  call test_precision ( )
  call test_present ( arg_a, arg_b )
  call test_present ( arg_a, arg_b, fred = 7 )
  call test_present ( arg_a, arg_b, george = 4, fred = 3 )
  call test_present ( arg_a, arg_b, ethel = 1, fred = 2, george = 7 )
  call test_product ( )
  call test_radix ( )
  call test_random_number ( )
  call test_random_seed ( )
  call test_range ( )
  call test_real_c4 ( )
  call test_repeat ( )
  call test_reshape ( )
  call test_rrspacing ( )
  call test_scale ( )
  call test_scan ( )
  call test_selected_int_kind ( )
  call test_selected_real_kind ( )
  call test_set_exponent ( )
  call test_shape ( )
  call test_sign ( )
  call test_sin_r8 ( )
  call test_sinh ( )
  call test_size ( )
  call test_spacing ( )
  call test_spread ( )
  call test_sqrt ( )
  call test_sum ( )
  call test_sum_dim ( )
  call test_system_clock ( )
  call test_tan ( )
  call test_tanh ( )
  call test_tiny ( )
  call test_transfer ( )
  call test_transpose ( )
  call test_trim ( )
  call test_ubound ( )
  call test_unpack ( )
  call test_verify ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'F90_INTRINSICS'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test_abs ( )

!*****************************************************************************80
!
!! TEST_ABS tests ABS.
!
!  Discussion:
!
!    The FORTRAN90 function ABS returns the absolute value of a 
!    number.  For complex numbers, this is the magnitude.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) :: c4 = ( 1.0E+00, -2.0E+00 )
  complex ( kind = 8 ) :: c8 = ( 1.0D+00, -2.0D+00 )
  integer ( kind = 4 ) :: i4 = -88
  integer ( kind = 8 ) :: i8 = -88
  real ( kind = 4 ) :: r4 = 45.78E+00
  real ( kind = 8 ) :: r8 = 45.78D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ABS'
  write ( *, '(a)' ) '  ABS is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  absolute value of a numeric quantity'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Type      VALUE                ABS(VALUE)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12,12x,i12)' )     '  Integer ( kind = 4 ) ', i4, abs ( i4 )
  write ( *, '(a,i12,12x,i12)' )     '  Integer ( kind = 8 ) ', i8, abs ( i8 )
  write ( *, '(a,f12.4,12x,f12.4)' ) '  Real    ( kind = 4 ) ', r4, abs ( r4 )
  write ( *, '(a,f12.4,12x,f12.4)' ) '  Real    ( kind = 8 ) ', r8, abs ( r8 )
  write ( *, '(a,2f12.4,f12.4)' )    '  Complex ( kind = 4 ) ', c4, abs ( c4 )
  write ( *, '(a,2f12.4,f12.4)' )    '  Complex ( kind = 8 ) ', c8, abs ( c8 )

  return
end
subroutine test_achar ( )

!*****************************************************************************80
!
!! TEST_ACHAR tests ACHAR.
!
!  Discussion:
!
!    The FORTRAN90 function ACHAR returns the character corresponding
!    to the given ASCII index, between 0 and 255.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character c
  logical ch_is_printable
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ACHAR'
  write ( *, '(a)' ) '  ACHAR is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  character of given ASCII index.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I     ACHAR(I)'
  write ( *, '(a)' ) ' '

  do i = 0, 255

    c = achar ( i )

    if ( ch_is_printable ( c ) ) then
      write ( *, '(2x,i8,8x,a1)' ) i, c
    end if

  end do

  return
end
subroutine test_acos ( )

!*****************************************************************************80
!
!! TEST_ACOS tests ACOS.
!
!  Discussion:
!
!    The FORTRAN90 function ACOS returns the inverse cosine of a number X.
!    assuming -1 <= X <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 8 ) i
  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 1.0D+00
  real ( kind = 8 ) :: x_lo = -1.0D+00
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ACOS'
  write ( *, '(a)' ) '  ACOS is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  inverse cosine of a value between -1 and 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X               ACOS(X)     COS(ACOS(X))'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = acos ( x )
    z = cos ( y )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, y, z
  end do

  return
end
subroutine test_adjustl ( )

!*****************************************************************************80
!
!! TEST_ADJUSTL tests ADJUSTL.
!
!  Discussion:
!
!    The FORTRAN90 function ADJUSTL returns a left-adjusted copy of a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) s1
  character ( len = 10 ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ADJUSTL'
  write ( *, '(a)' ) '  ADJUSTL is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  left-adjusted copy of a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      S       ADJUSTL(S)'
  write ( *, '(a)' ) '  ----------  ----------'
  write ( *, '(a)' ) ' '
  s1 = '1234567890'
  s2 = adjustl ( s1 )
  write ( *, '(2x,a10,2x,a10)' ) s1, s2
  s1 = '12345     '
  s2 = adjustl ( s1 )
  write ( *, '(2x,a10,2x,a10)' ) s1, s2
  s1 = '     67890'
  s2 = adjustl ( s1 )
  write ( *, '(2x,a10,2x,a10)' ) s1, s2
  s1 = '  34 678  '
  s2 = adjustl ( s1 )
  write ( *, '(2x,a10,2x,a10)' ) s1, s2
  s1 = '    5     '
  s2 = adjustl ( s1 )
  write ( *, '(2x,a10,2x,a10)' ) s1, s2

  return
end
subroutine test_adjustr ( )

!*****************************************************************************80
!
!! TEST_ADJUSTR tests ADJUSTR.
!
!  Discussion:
!
!    The FORTRAN90 function ADJUSTR returns a right-adjusted copy of a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) s1
  character ( len = 10 ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ADJUSTR'
  write ( *, '(a)' ) '  ADJUSTR is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  right-adjusted copy of a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      S       ADJUSTR(S)'
  write ( *, '(a)' ) '  ----------  ----------'
  write ( *, '(a)' ) ' '
  s1 = '1234567890'
  s2 = adjustr ( s1 )
  write ( *, '(2x,a10,2x,a10)' ) s1, s2
  s1 = '12345     '
  s2 = adjustr ( s1 )
  write ( *, '(2x,a10,2x,a10)' ) s1, s2
  s1 = '     67890'
  s2 = adjustr ( s1 )
  write ( *, '(2x,a10,2x,a10)' ) s1, s2
  s1 = '  34 678  '
  s2 = adjustr ( s1 )
  write ( *, '(2x,a10,2x,a10)' ) s1, s2
  s1 = '    5     '
  s2 = adjustr ( s1 )
  write ( *, '(2x,a10,2x,a10)' ) s1, s2

  return
end
subroutine test_aimag ( )

!*****************************************************************************80
!
!! TEST_AIMAG tests AIMAG.
!
!  Discussion:
!
!    The FORTRAN90 function AIMAG returns the imaginary part of a
!    complex number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex c4_uniform_01
  complex c
  integer i
  integer :: seed = 123456789
  real r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_AIMAG'
  write ( *, '(a)' ) '  AIMAG is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  imaginary part of a complex number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                  X                      AIMAG(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    c = c4_uniform_01 ( seed )
    r = aimag ( c )
    write ( *, '(2x,f14.6,f14.6,6x,f14.6,f14.6)' ) c, r
  end do

  return
end
subroutine test_aint ( )

!*****************************************************************************80
!
!! TEST_AINT tests AINT.
!
!  Discussion:
!
!    The FORTRAN90 function AINT returns a real number rounded towards
!    zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 8 ) i
  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_AINT'
  write ( *, '(a)' ) '  AINT is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  value of a real number rounded towards zero.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              AINT(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = aint ( x )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, y
  end do

  return
end
subroutine test_all ( )

!*****************************************************************************80
!
!! TEST_ALL tests ALL.
!
!  Discussion:
!
!    The FORTRAN90 function ALL(MASK) returns TRUE if every MASK is true.
!    ALL(MASK,DIM) returns an array of values which are TRUE where every
!    entry in dimension DIM is TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer a_i4(4,5)
  integer i
  integer j
  real r

  do i = 1, 4
    do j = 1, 5
      call random_number ( harvest = r )
      a_i4(i,j) = int ( 100.0 * r )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ALL'
  write ( *, '(a)' ) '  ALL(MASK) is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  TRUE if every entry of MASK is true.'
  write ( *, '(a)' ) '  ALL(MASK,DIM) is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  an array of values which are TRUE where every entry'
  write ( *, '(a)' ) '  in dimension DIM is true.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer a_i4(4,5)'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(5(2x,i4))' ) a_i4(i,1:5)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  ALL ( A_I4 < 100 )      ', all ( a_i4 < 100 )
  write ( *, '(a,l1)' ) '  ALL ( A_I4 <  90 )      ', all ( a_i4 < 90 )
  write ( *, '(a,5(2x,l1))' ) &
    '  ALL ( A_I4 <  90,    1 )', all ( a_i4 < 90, 1 )
  write ( *, '(a,5(2x,l1))' ) &
    '  ALL ( A_I4 <  90,    2 )', all ( a_i4 < 90, 2 )

  write ( *, '(a,5(2x,l1))' ) &
    '  ALL ( MOD(A_I4,2)==0,1 )', all ( mod(a_i4,2)==0,1)

  write ( *, '(a,5(2x,l1))' ) &
    '  ALL ( MOD(A_I4,2)==0,2 )', all ( mod(a_i4,2)==0,2)

  return
end
subroutine test_allocated ( )

!*****************************************************************************80
!
!! TEST_ALLOCATED tests ALLOCATED.
!
!  Discussion:
!
!    The FORTRAN90 function ALLOCATED(ARRAY) returns TRUE if the allocatable
!    array ARRAY has been allocated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ), allocatable, dimension ( : ) :: a_c4
  integer ( kind = 4 ), allocatable, dimension ( : ) :: a_i4
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ALLOCATED'
  write ( *, '(a)' ) '  ALLOCATED(ARRAY) is a FORTRAN90 function which is'
  write ( *, '(a)' ) '  TRUE if the allocatable array ARRAY has actually'
  write ( *, '(a)' ) '  been allocated.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initial status:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    A_C4 <- Not allocated.'
  write ( *, '(a)' ) '    A_I4 <- Not allocated.'
  write ( *, '(a)' ) '    A_R8 <- Not allocated.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call ALLOCATED:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '    ALLOCATED(A_C4) = ', allocated ( a_c4 )
  write ( *, '(a,l1)' ) '    ALLOCATED(A_I4) = ', allocated ( a_i4 )
  write ( *, '(a,l1)' ) '    ALLOCATED(A_R8) = ', allocated ( a_r8 )

  allocate ( a_c4(1:3) )
  allocate ( a_i4(1:10) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Allocate A_C4 and A_I4:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '    ALLOCATED(A_C4) = ', allocated ( a_c4 )
  write ( *, '(a,l1)' ) '    ALLOCATED(A_I4) = ', allocated ( a_i4 )
  write ( *, '(a,l1)' ) '    ALLOCATED(A_R8) = ', allocated ( a_r8 )

  deallocate ( a_c4 )
  allocate ( a_r8(1:3,1:3) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Deallocate A_C4, allocate A_R8:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '    ALLOCATED(A_C4) = ', allocated ( a_c4 )
  write ( *, '(a,l1)' ) '    ALLOCATED(A_I4) = ', allocated ( a_i4 )
  write ( *, '(a,l1)' ) '    ALLOCATED(A_R8) = ', allocated ( a_r8 )

  deallocate ( a_i4 )
  deallocate ( a_r8 )

  return
end
subroutine test_anint ( )

!*****************************************************************************80
!
!! TEST_ANINT tests ANINT.
!
!  Discussion:
!
!    The FORTRAN90 function ANINT returns, as a real value, the nearest 
!    integer to a given real value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 8 ) i
  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ANINT'
  write ( *, '(a)' ) '  ANINT is a FORTRAN90 function which returns,'
  write ( *, '(a)' ) '  as a real value, the nearest integer to a '
  write ( *, '(a)' ) '  given real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X             ANINT(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = anint ( x )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, y
  end do

  return
end
subroutine test_any ( )

!*****************************************************************************80
!
!! TEST_ANY tests ANY.
!
!  Discussion:
!
!    The FORTRAN90 function ANY(MASK) returns TRUE if any MASK is true.
!    ANY(MASK,DIM) returns an array of values which are TRUE where any
!    entry in dimension DIM is TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer a_i4(4,5)
  integer i
  integer j
  real r
 
  do i = 1, 4
    do j = 1, 5
      call random_number ( harvest = r )
      a_i4(i,j) = int ( 100.0 * r )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ANY'
  write ( *, '(a)' ) '  ANY(MASK) is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  TRUE if any entry of MASK is true.'
  write ( *, '(a)' ) '  ANY(MASK,DIM) is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  an array of values which are TRUE where any entry'
  write ( *, '(a)' ) '  in dimension DIM is true.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer a_i4(4,5)'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(4(2x,i4))' ) a_i4(i,1:5)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  ANY ( A_I4 < 10 )      ', any ( a_i4 < 10 )
  write ( *, '(a,l1)' ) '  ANY ( A_I4 <  0 )      ', any ( a_i4 < 0 )
  write ( *, '(a,5(2x,l1))' ) &
    '  ANY ( A_I4 <  10,    1 )', any ( a_i4 < 10, 1 )
  write ( *, '(a,5(2x,l1))' ) &
    '  ANY ( A_I4 <  10,    2 )', any ( a_i4 < 10, 2 )

  write ( *, '(a,5(2x,l1))' ) &
    '  ANY ( MOD(A_I4,5)==0,1 )', any ( mod(a_i4,5)==0,1)
  write ( *, '(a,5(2x,l1))' ) &
    '  ANY ( MOD(A_I4,5)==0,2 )', any ( mod(a_i4,5)==0,2)

  return
end
subroutine test_asin ( )

!*****************************************************************************80
!
!! TEST_ASIN tests ASIN.
!
!  Discussion:
!
!    The FORTRAN90 function ASIN returns the inverse sine of a number X.
!    assuming -1 <= X <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 8 ) i
  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 1.0D+00
  real ( kind = 8 ) :: x_lo = -1.0D+00
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ASIN'
  write ( *, '(a)' ) '  ASIN is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  inverse sine of a value between -1 and 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X               ASIN(X)     SIN(ASIN(X))'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = asin ( x )
    z = sin ( y )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, y, z
  end do

  return
end
subroutine test_associated ( )

!*****************************************************************************80
!
!! TEST_ASSOCIATED tests ASSOCIATED.
!
!  Discussion:
!
!    The FORTRAN90 function ASSOCIATED(POINTER,TARGET) returns TRUE if 
!    POINTER is asscoaited with the given TARGET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ), target :: a_c4
  integer ( kind = 4 ), target :: a_i4
  real ( kind = 8 ), target :: a_r8
  complex ( kind = 4 ), pointer :: p_c4
  integer ( kind = 4 ), pointer :: p_i4
  real ( kind = 8 ), pointer :: p_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ASSOCIATED'
  write ( *, '(a)' ) '  ASSOCIATED(POINTER,TARGET) is a FORTRAN90 function'
  write ( *, '(a)' ) '  returns the association status of a pointer, or'
  write ( *, '(a)' ) '  indicates the pointer is associated with the target.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initial status:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P_C4 <- Not associated.'
  write ( *, '(a)' ) '    P_I4 <- Not associated.'
  write ( *, '(a)' ) '    P_R8 <- Not associated.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call ASSOCIATED:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) &
    '    ASSOCIATED(P_C4,A_C4) = ', associated ( p_c4, a_c4 )
  write ( *, '(a,l1)' ) &
    '    ASSOCIATED(P_I4,A_I4) = ', associated ( p_i4, a_i4 )
  write ( *, '(a,l1)' ) &
    '    ASSOCIATED(P_R8,A_R8) = ', associated ( p_r8, a_r8 )

  p_c4 => a_c4
  p_i4 => a_i4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Point to A_C4 and A_I4:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) &
    '    ASSOCIATED(P_C4,A_C4) = ', associated ( p_c4, a_c4 )
  write ( *, '(a,l1)' ) &
    '    ASSOCIATED(P_I4,A_I4) = ', associated ( p_i4, a_i4 )
  write ( *, '(a,l1)' ) &
    '    ASSOCIATED(P_R8,A_R8) = ', associated ( p_r8, a_r8 )

  nullify ( p_c4 )
  p_r8 => a_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Nullify A_C4, point to A_R8:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) &
    '    ASSOCIATED(P_C4,A_C4) = ', associated ( p_c4, a_c4 )
  write ( *, '(a,l1)' ) &
    '    ASSOCIATED(P_I4,A_I4) = ', associated ( p_i4, a_i4 )
  write ( *, '(a,l1)' ) &
    '    ASSOCIATED(P_R8,A_R8) = ', associated ( p_r8, a_r8 )

  nullify ( p_i4 )
  nullify ( p_r8 )

  return
end
subroutine test_atan ( )

!*****************************************************************************80
!
!! TEST_ATAN tests ATAN.
!
!  Discussion:
!
!    The FORTRAN90 function ATAN returns the inverse tangent of a number X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 8 ) i
  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ATAN'
  write ( *, '(a)' ) '  ATAN is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  inverse tangent of a value'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X               ATAN(X)     TAN(ATAN(X))'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = atan ( x )
    z = tan ( y )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, y, z
  end do

  return
end
subroutine test_atan2 ( )

!*****************************************************************************80
!
!! TEST_ATAN2 tests ATAN2.
!
!  Discussion:
!
!    The FORTRAN90 function ATAN2 returns the inverse tangent of a number X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ATAN2'
  write ( *, '(a)' ) '  ATAN2 is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  inverse tangent of a value'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X               Y           ATAN2(Y,X)  TAN(ATAN2(Y,X))'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = r8_uniform ( x_lo, x_hi, seed )
    z = atan2 ( y, x )
    w = tan ( z )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) x, y, z, w
  end do

  return
end
subroutine test_bit_size ( )

!*****************************************************************************80
!
!! TEST_BIT_SIZE tests BIT_SIZE.
!
!  Discussion:
!
!    The FORTRAN90 function BIT_SIZE returns the size of an integer word
!    in bits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer ( kind = 4 ) i4
  integer ( kind = 8 ) i8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_BIT_SIZE'
  write ( *, '(a)' ) '  BIT_SIZE is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  size of an integer in bits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       Type(X)             BIT_SIZE(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  INTEGER              ', bit_size ( i )
  write ( *, '(a,i8)' ) '  INTEGER ( KIND = 4 ) ', bit_size ( i4 )
  write ( *, '(a,i8)' ) '  INTEGER ( KIND = 8 ) ', bit_size ( i8 )

  return
end
subroutine test_btest ( )

!*****************************************************************************80
!
!! TEST_BTEST tests BTEST.
!
!  Discussion:
!
!    The FORTRAN90 function BTEST reports whether a given bit is 0 or 1
!    in an integer word.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i1
  integer i2
  logical l
  integer pos
  integer pos_high
  character ( len = 32 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_BTEST'
  write ( *, '(a)' ) '  BTEST(I,POS) is a FORTRAN90 function which is TRUE'
  write ( *, '(a)' ) '  if bit POS of I is 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we are only going to check the lowest 32 bits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I       POS    BTEST(I,POS)'
  write ( *, '(a)' ) ' '

  i1 = 213456
  pos_high = min ( 32, bit_size(i1) )

  do pos = 0, pos_high - 1

    l = btest ( i1, pos )

    if ( l ) then
      s(pos_high-pos:pos_high-pos) = '1'
    else
      s(pos_high-pos:pos_high-pos) = '0'
    end if

    write ( *, '(2x,i8,2x,i8,10x,l1)' ) i1, pos, l

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12,a)' ) '  The binary representation of ', i1, ' is:'
  write ( *, '(a)' ) '  "' // s //'".' 

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I       POS    BTEST(I,POS)'
  write ( *, '(a)' ) ' '

  i1 = -28

  do pos = 0, pos_high - 1

    l = btest ( i1, pos )

    if ( l ) then
      s(pos_high-pos:pos_high-pos) = '1'
    else
      s(pos_high-pos:pos_high-pos) = '0'
    end if

    write ( *, '(2x,i8,2x,i8,10x,l1)' ) i1, pos, l

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12,a)' ) '  The binary representation of ', i1, ' is:'
  write ( *, '(a)' ) '  "' // s //'".' 

  return
end
subroutine test_ceiling ( )

!*****************************************************************************80
!
!! TEST_CEILING tests CEILING.
!
!  Discussion:
!
!    The FORTRAN90 function CEILING rounds a real value up "towards infinity".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  integer ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_CEILING'
  write ( *, '(a)' ) '  CEILING is a FORTRAN90 function which rounds a'
  write ( *, '(a)' ) '  real number up "towards infinity".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              CEILING(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = ceiling ( x )
    write ( *, '(2x,g14.6,2x,i8)' ) x, y
  end do

  return
end
subroutine test_char ( )

!*****************************************************************************80
!
!! TEST_CHAR tests CHAR.
!
!  Discussion:
!
!    The FORTRAN90 function CHAR returns the character corresponding
!    to the given character index, between 0 and 255.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character c
  logical ch_is_printable
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_CHAR'
  write ( *, '(a)' ) '  CHAR is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  character of given character index.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I     CHAR(I)'
  write ( *, '(a)' ) ' '

  do i = 0, 255

    c = char ( i )

    if ( ch_is_printable ( c ) ) then
      write ( *, '(2x,i8,8x,a1)' ) i, c
    end if

  end do

  return
end
subroutine test_cmplx ( )

!*****************************************************************************80
!
!! TEST_CMPLX tests CMPLX.
!
!  Discussion:
!
!    The FORTRAN90 function CMPLX returns a complex number given its
!    real and imaginary parts.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_CMPLX'
  write ( *, '(a)' ) '  CMPLX is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  complex number formed by real and imaginary parts.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.6,f14.6)' ) '  CMPLX(1)        ', cmplx ( 1 )
  write ( *, '(a,f14.6,f14.6)' ) '  CMPLX(2,3)      ', cmplx ( 2, 3 )
  write ( *, '(a,f14.6,f14.6)' ) '  CMPLX(4.5)      ', cmplx ( 4.5 )
  write ( *, '(a,f14.6,f14.6)' ) '  CMPLX(6.7, 8.9 )', cmplx ( 6.7, 8.9 )

  return
end
subroutine test_conjg ( )

!*****************************************************************************80
!
!! TEST_CONJG tests CONJG.
!
!  Discussion:
!
!    The FORTRAN90 function CONJG returns the conjugate of a complex number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex c4_uniform_01
  complex c1
  complex c2
  integer i
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_CONJG'
  write ( *, '(a)' ) '  CONJG is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  conjugate of a complex number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                  X                             CONJG(X)'
  write ( *, '(a)' ) &
    '     --------------------------      ----------------------------'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    c1 = c4_uniform_01 ( seed )
    c2 = conjg ( c1 )
    write ( *, '(2x,f14.6,f14.6,6x,f14.6,f14.6)' ) c1, c2
  end do

  return
end
subroutine test_cos_c4 ( )

!*****************************************************************************80
!
!! TEST_COS_C4 tests COS on complex ( kind = 4 ) arguments.
!
!  Discussion:
!
!    The FORTRAN90 function COS returns the cosine of a real or complex number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex c4_uniform_01
  integer i
  integer :: seed = 123456789
  complex ( kind = 4 ) x_c4
  complex ( kind = 4 ) y_c4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_COS_C4'
  write ( *, '(a)' ) '  COS is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  cosine of a real or complex number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we use complex ( kind = 4 ) arguments.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '              X                          COS(X)'
  write ( *, '(a)' ) &
    '    --------------------------  ----------------------------'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x_c4 = c4_uniform_01 ( seed )
    y_c4 = cos ( x_c4 )
    write ( *, '(2x,2g14.6,2x,2g14.6)' ) x_c4, y_c4
  end do

  return
end
subroutine test_cos_r8 ( )

!*****************************************************************************80
!
!! TEST_COS_R8 tests COS on real ( kind = 8 ) arguments.
!
!  Discussion:
!
!    The FORTRAN90 function COS returns the cosine of a real or complex number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_COS_R8'
  write ( *, '(a)' ) '  COS is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  cosine of a real or complex number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we use real ( kind = 8 ) arguments.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              COS(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = cos ( x )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, y
  end do

  return
end
subroutine test_cosh ( )

!*****************************************************************************80
!
!! TEST_COSH tests COSH.
!
!  Discussion:
!
!    The FORTRAN90 function COSH returns the hyperbolic cosine of a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_COSH'
  write ( *, '(a)' ) '  COSH is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  hyperbolic cosine of a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              COSH(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = cosh ( x )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, y
  end do

  return
end
subroutine test_count ( )

!*****************************************************************************80
!
!! TEST_COUNT tests COUNT.
!
!  Discussion:
!
!    The FORTRAN90 function COUNT(MASK) returns the number of TRUE entries
!    in MASK.
!    COUNT(MASK,DIM) returns an array of the number of TRUE entries of
!    MASK in dimension DIM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a_i4(4,5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
 
  do i = 1, 4
    do j = 1, 5
      call random_number ( harvest = r )
      a_i4(i,j) = int ( 100.0 * r )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_COUNT'
  write ( *, '(a)' ) '  COUNT(MASK) is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the number of TRUE entries in MASK.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  COUNT(MASK,DIM) is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  an array of values which are count the number of TRUE'
  write ( *, '(a)' ) '  entries in MASK in dimension DIM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer ( kind = 4 ) a_i4(4,5)'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(5(2x,i4))' ) a_i4(i,1:5)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  COUNT ( A_I4 < 100 )      ', count ( a_i4 < 100 )
  write ( *, '(a,i4)' ) '  COUNT ( A_I4 <  90 )      ', count ( a_i4 < 90 )
  write ( *, '(a,5(2x,i4))' ) &
    '  COUNT ( A_I4 <  90,    1 )', count ( a_i4 < 90, 1 )
  write ( *, '(a,5(2x,i4))' ) &
    '  COUNT ( A_I4 <  90,    2 )', count ( a_i4 < 90, 2 )

  k = 2

  write ( *, '(a,5(2x,i4))' ) &
    '  COUNT ( MOD(A_I4,2)==0,1 )', count ( mod(a_i4,k)==0, 1 )
  write ( *, '(a,5(2x,i4))' ) &
    '  COUNT ( MOD(A_I4,2)==0,2 )', count ( mod(a_i4,k)==0, 2 )

  return
end
subroutine test_cshift ( )

!*****************************************************************************80
!
!! TEST_CSHIFT tests CSHIFT.
!
!  Discussion:
!
!    The FORTRAN90 function CSHIFT performs a circular shift on a vector
!    or array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a1_i4(4,5)
  integer ( kind = 4 ) a2_i4(4,5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer :: seed = 123456789
  integer ( kind = 4 ) v1_i4(10)
  integer ( kind = 4 ) v2_i4(10)

  do i = 1, 10
    v1_i4(i) = i
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_CSHIFT'
  write ( *, '(a)' ) '  CSHIFT is a FORTRAN90 function which applies a'
  write ( *, '(a)' ) '  circular shift to a vector or array.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer ( kind = 4 ) v1_i4(10):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,10i4)' ) v1_i4(1:10)

  v2_i4 = cshift ( v1_i4, 3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  v2_i4 = cshift ( v1_i4, 3 ):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,10i4)' ) v2_i4(1:10)

  v2_i4 = cshift ( v1_i4, -2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  v2_i4 = cshift ( v1_i4, -2 ):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,10i4)' ) v2_i4(1:10)

  do i = 1, 4
    do j = 1, 5
      a1_i4(i,j) = 10 * i + j
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer ( kind = 4 ) a1_i4(4,5):'
  write ( *, '(a)' ) ' '

  do i = 1, 4
    write ( *, '(2x,5i4)' ) a1_i4(i,1:5)
  end do

  a2_i4 = cshift ( a1_i4, shift = (/ 1, 2, 3, 4 /), dim = 2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  a2_i4 = cshift ( a1_i4, (/ 1, 2, 3, 4 /), dim = 2 ):'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(2x,5i4)' ) a2_i4(i,1:5)
  end do

  a2_i4 = cshift ( a1_i4, 2, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  a2_i4 = cshift ( a1_i4, 2, 1 ):'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(2x,5i4)' ) a2_i4(i,1:5)
  end do

  return
end
subroutine test_date_and_time ( )

!*****************************************************************************80
!
!! TEST_DATE_AND_TIME tests DATE_AND_TIME.
!
!  Discussion:
!
!    The FORTRAN90 function DATE_AND_TIME returns date and time information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 8 ) date
  character ( len = 10 ) time
  integer values(8)
  character ( len = 5 ) zone

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DATE_AND_TIME'
  write ( *, '(a)' ) '  DATE_AND_TIME is a FORTRAN90 subroutine which returns'
  write ( *, '(a)' ) '  date and time information.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  call date_and_time ( date, time, zone, values )'

  call date_and_time ( date, time, zone, values )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DATE = "' // date // '".'
  write ( *, '(a)' ) '  TIME = "' // time // '".'
  write ( *, '(a)' ) '  ZONE = "' // zone // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  VALUES(1) = year =            ', values(1)
  write ( *, '(a,i6)' ) '  VALUES(2) = month =           ', values(2)
  write ( *, '(a,i6)' ) '  VALUES(3) = day =             ', values(3)
  write ( *, '(a,i6)' ) '  VALUES(4) = UTC minute diff = ', values(4)
  write ( *, '(a,i6)' ) '  VALUES(5) = hour =            ', values(5)
  write ( *, '(a,i6)' ) '  VALUES(6) = minute =          ', values(6)
  write ( *, '(a,i6)' ) '  VALUES(7) = second =          ', values(7)
  write ( *, '(a,i6)' ) '  VALUES(8) = milliseconds =    ', values(8)

  return
end
subroutine test_dble ( )

!*****************************************************************************80
!
!! TEST_DBLE tests DBLE.
!
!  Discussion:
!
!    The FORTRAN90 function DBLE converts a numeric value to double precision.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) :: x_c4 = ( 1.1E+00, 2.2E+00 )
  complex ( kind = 8 ) :: x_c8 = ( 3.3D+00, 4.4D+00 )
  integer ( kind = 4 ) :: x_i4 = 5
  integer ( kind = 8 ) :: x_i8 = 6
  real ( kind = 4 ) :: x_r4 = 7.7E+00
  real ( kind = 8 ) :: x_r8 = 8.8D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DBLE'
  write ( *, '(a)' ) '  DBLE is a FORTRAN90 function which converts'
  write ( *, '(a)' ) '  a complex, integer or real value to double precision'
  write ( *, '(a)' ) '  real'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Type                   X             DBLE(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f6.4,2x,f6.4,2x,f6.4)' ) '  complex ( kind = 4 ) ', &
    x_c4, dble ( x_c4 )
  write ( *, '(a,f6.4,2x,f6.4,2x,f6.4)' ) '  complex ( kind = 8 ) ', &
    x_c8, dble ( x_c8 )
  write ( *, '(a,i6,10x,f6.4)'    ) '  integer ( kind = 4 ) ', &
    x_i4, dble ( x_i4 )
  write ( *, '(a,i6,10x,f6.4)'    ) '  integer ( kind = 8 ) ', &
    x_i8, dble ( x_i8 )
  write ( *, '(a,f6.4,10x,f6.4)'  ) '  real ( kind = 4 )    ', &
    x_r4, dble ( x_r4 )
  write ( *, '(a,f6.4,10x,f6.4)'  ) '  real ( kind = 8 )    ', &
    x_r8, dble ( x_r8 )

  return
end
subroutine test_digits ( )

!*****************************************************************************80
!
!! TEST_DIGITS tests DIGITS.
!
!  Discussion:
!
!    The FORTRAN90 function DIGITS returns the number of significant
!    digits in the model.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer x_i
  integer ( kind = 4 ) x_i4
  integer ( kind = 8 ) x_i8
  real x_r
  real ( kind = 4 ) x_r4
  real ( kind = 8 ) x_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DIGITS'
  write ( *, '(a)' ) '  DIGITS is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the number of significant digits for numbers'
  write ( *, '(a)' ) '  of the same kind as X.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Type         DIGITS(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  integer              ', digits ( x_i )
  write ( *, '(a,i8)' ) '  integer ( kind = 4 ) ', digits ( x_i4 )
  write ( *, '(a,i8)' ) '  integer ( kind = 8 ) ', digits ( x_i8 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  real                 ', digits ( x_r )
  write ( *, '(a,i8)' ) '  real ( kind = 4 )    ', digits ( x_r4 )
  write ( *, '(a,i8)' ) '  real ( kind = 8 )    ', digits ( x_r8 )

  return
end
subroutine test_dim ( )

!*****************************************************************************80
!
!! TEST_DIM tests DIM.
!
!  Discussion:
!
!    The FORTRAN90 function DIM(X,Y) returns the maximum of (X-Y) and 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) r
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
  integer ( kind = 4 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DIM'
  write ( *, '(a)' ) '  DIM is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the maximum of X-Y or 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Arithmetic type: integer X, Y'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X       Y  DIM(X,Y)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call random_number ( harvest = r )
    x = int ( 100.0E+00 * r )

    call random_number ( harvest = r )
    y = int ( 100.0E+00 * r )

    z = dim ( x, y )

    write ( *, '(2x,i6,2x,i6,2x,i6)' ) x, y, z

  end do

  return
end
subroutine test_dot_product ( )

!*****************************************************************************80
!
!! TEST_DOT_PRODUCT tests DOT_PRODUCT.
!
!  Discussion:
!
!    The FORTRAN90 function DOT_PRODUCT returns the dot product of
!    two complex, integer, logical or real vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 5

  complex ( kind = 4 ) c4_uniform_01
  integer i
  integer :: seed = 123456789
  complex ( kind = 4 ) x_c4(n)
  integer ( kind = 4 ) x_i4(n)
  logical ( kind = 4 ) x_l4(n)
  real ( kind = 4 ) x_r4(n)
  complex ( kind = 4 ) y_c4(n)
  integer ( kind = 4 ) y_i4(n)
  logical ( kind = 4 ) y_l4(n)
  real ( kind = 4 ) y_r4(n)

  call random_number ( harvest = x_r4(1:n) )
  call random_number ( harvest = y_r4(1:n) )

  x_c4(1:n) = cmplx ( x_r4(1:n), y_r4(1:n) )

  call random_number ( harvest = x_r4(1:n) )
  call random_number ( harvest = y_r4(1:n) )

  y_c4(1:n) = cmplx ( x_r4(1:n), y_r4(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DOT_PRODUCT'
  write ( *, '(a)' ) '  DOT_PRODUCT is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  dot product of two vectors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If the arguments are COMPLEX, then the first vector'
  write ( *, '(a)' ) '  will be conjugated to produce the result.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Two complex ( kind = 4 ) vectors:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                 X                              Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,2f14.6,2x,2f14.6)' ) x_c4(i), y_c4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,2f14.6)' ) '  DOT_PRODUCT(X,Y) = ', dot_product ( x_c4, y_c4 )

  call random_number ( harvest = x_r4(1:n) )
  call random_number ( harvest = y_r4(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Two real ( kind = 4 ) vectors:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X               Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,f14.6,2x,f14.6)' ) x_r4(i), y_r4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.6)' ) '  DOT_PRODUCT(X,Y) = ', dot_product ( x_r4, y_r4 )

  do i = 1, n
    x_i4(i) = int ( 10.0E+00 * x_r4(i) )
    y_i4(i) = int ( 10.0E+00 * y_r4(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Two integer ( kind = 4 ) vectors:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               X               Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i14,2x,i14)' ) x_i4(i), y_i4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,i14)' ) '  DOT_PRODUCT(X,Y) = ', dot_product ( x_i4, y_i4 )

  do i = 1, n
    x_l4(i) = ( 5 < x_i4(i) )
    y_l4(i) = ( 5 < y_i4(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Two logical ( kind = 4 ) vectors:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               X               Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,13x,l1,2x,13x,l1)' ) x_l4(i), y_l4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  DOT_PRODUCT(X,Y) = ', dot_product ( x_l4, y_l4 )

  return
end
subroutine test_dprod ( )

!*****************************************************************************80
!
!! TEST_DPROD tests DPROD.
!
!  Discussion:
!
!    The FORTRAN90 function DPROD(X,Y) returns the product of real
!    numbers X and Y using double precision.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  real ( kind = 4 ) x_r4
  real ( kind = 4 ) y_r4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DPROD'
  write ( *, '(a)' ) '  DPROD is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the product of real values X and Y, using'
  write ( *, '(a)' ) '  double precision.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  real ( kind = 4 ) x, y'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X       Y  DPROD(X,Y)'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    call random_number ( harvest = x_r4 )
    x_r4 = 100.0 * x_r4
    call random_number ( harvest = y_r4 )
    y_r4 = 100.0 * y_r4
    write ( *, '(2x,f12.6,2x,f12.6,2x,f12.6)' ) x_r4, y_r4, dprod ( x_r4, y_r4 )
  end do

  return
end
subroutine test_eoshift ( )

!*****************************************************************************80
!
!! TEST_EOSHIFT tests EOSHIFT.
!
!  Discussion:
!
!    The FORTRAN90 function EOSHIFT performs an "end off" shift on a vector
!    or array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer a1_i4(4,5)
  integer a2_i4(4,5)
  integer boundary
  integer i
  integer j
  integer v1_i4(10)
  integer v2_i4(10)

  do i = 1, 10
    v1_i4(i) = i
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_EOSHIFT'
  write ( *, '(a)' ) '  EOSHIFT is a FORTRAN90 function which applies an'
  write ( *, '(a)' ) '  "end off" shift to a vector or array.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer v1_i4(10):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,10i4)' ) v1_i4(1:10)

  v2_i4 = eoshift ( v1_i4, 3, 99 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  v2_i4 = eoshift ( v1_i4, 3, 99 ):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,10i4)' ) v2_i4(1:10)

  v2_i4 = eoshift ( v1_i4, -2, 88 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  v2_i4 = eoshift ( v1_i4, -2, 88 ):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,10i4)' ) v2_i4(1:10)

  do i = 1, 4
    do j = 1, 5
      a1_i4(i,j) = 10 * i + j
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer a1_i4(4,5):'
  write ( *, '(a)' ) ' '

  do i = 1, 4
    write ( *, '(2x,5i4)' ) a1_i4(i,1:5)
  end do

  a2_i4 = eoshift ( a1_i4, (/ 1, 2, 3, 4 /), 77, dim = 2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  a2_i4 = eoshift ( a1_i4, (/ 1, 2, 3, 4 /), 77, dim = 2 ):'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(2x,5i4)' ) a2_i4(i,1:5)
  end do

  a2_i4 = eoshift ( a1_i4, 2, 66, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  a2_i4 = eoshift ( a1_i4, 2, 66, 1 ):'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(2x,5i4)' ) a2_i4(i,1:5)
  end do

  return
end
subroutine test_epsilon ( )

!*****************************************************************************80
!
!! TEST_EPSILON tests EPSILON.
!
!  Discussion:
!
!    The FORTRAN90 function EPSILON returns the "machine epsilon" associated
!    with a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) x_r4
  real ( kind = 8 ) x_r8
  real ( kind = 4 ) y_r4
  real ( kind = 8 ) y_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_EPSILON'
  write ( *, '(a)' ) '  EPSILON is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  "machine epsilon" associated with a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First, some "default precision" reals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              EPSILON(X)'
  write ( *, '(a)' ) ' '
  x_r4 = 1.0
  y_r4 = epsilon ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  x_r4 = 0.0
  y_r4 = epsilon ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  x_r4 = 1000000.0
  y_r4 = epsilon ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now, some "double precision" reals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              EPSILON(X)'
  write ( *, '(a)' ) ' '
  x_r8 = 1.0D+00
  y_r8 = epsilon ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8
  x_r8 = 0.0D+00
  y_r8 = epsilon ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8
  x_r8 = 1000000.0D+00
  y_r8 = epsilon ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8

  return
end
subroutine test_exp ( )

!*****************************************************************************80
!
!! TEST_EXP tests EXP.
!
!  Discussion:
!
!    The FORTRAN90 function EXP returns the exponential of a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_EXP'
  write ( *, '(a)' ) '  EXP is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  exponential of a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              EXP(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = exp ( x )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, y
  end do

  return
end
subroutine test_exponent ( )

!*****************************************************************************80
!
!! TEST_EXPONENT tests EXPONENT.
!
!  Discussion:
!
!    The FORTRAN90 function EXPONENT returns the exponent part of
!    a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_uniform
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 4 ) x
  real ( kind = 4 ) :: x_hi = 10.0E+00
  real ( kind = 4 ) :: x_lo = -10.0E+00
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_EXPONENT'
  write ( *, '(a)' ) '  EXPONENT is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the exponent part of a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X             EXPONENT(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r4_uniform ( x_lo, x_hi, seed )
    y = exponent ( x )
    write ( *, '(2x,g14.6,2x,i8)' ) x, y
  end do

  return
end
subroutine test_floor ( )

!*****************************************************************************80
!
!! TEST_FLOOR tests FLOOR.
!
!  Discussion:
!
!    The FORTRAN90 function FLOOR rounds a real value down "towards -infinity".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  integer ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_FLOOR'
  write ( *, '(a)' ) '  FLOOR is a FORTRAN90 function which rounds a'
  write ( *, '(a)' ) '  real number down "towards -infinity".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              FLOOR(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    call random_number ( harvest = r )
    x = x_lo + r * ( x_hi - x_lo )
    y = floor ( x )
    write ( *, '(2x,g14.6,2x,i8)' ) x, y
  end do

  return
end
subroutine test_fraction ( )

!*****************************************************************************80
!
!! TEST_FRACTION tests FRACTION.
!
!  Discussion:
!
!    The FORTRAN90 function FRACTION returns the fractional part of
!    a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  real r
  real x
  real :: x_hi = 10.0
  real :: x_lo = -10.0
  real y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_FRACTION'
  write ( *, '(a)' ) '  FRACTION is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the fractional part of a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X             FRACTION(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    call random_number ( harvest = r )
    x = x_lo + r * ( x_hi - x_lo )
    y = fraction ( x )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, y
  end do

  return
end
subroutine test_huge ( )

!*****************************************************************************80
!
!! TEST_HUGE tests HUGE.
!
!  Discussion:
!
!    The FORTRAN90 function HUGE returns a "huge value" associated
!    with an integer or real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer x_i4
  real ( kind = 4 ) x_r4
  real ( kind = 8 ) x_r8
  integer y_i4
  real ( kind = 4 ) y_r4
  real ( kind = 8 ) y_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_HUGE'
  write ( *, '(a)' ) '  HUGE is a FORTRAN90 function which returns a'
  write ( *, '(a)' ) '  "huge value" associated with a real or integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First, some "default precision" reals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              HUGE(X)'
  write ( *, '(a)' ) ' '
  x_r4 = 1.0
  y_r4 = huge ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  x_r4 = 0.0
  y_r4 = huge ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  x_r4 = 1000000.0
  y_r4 = huge ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now, some "double precision" reals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              HUGE(X)'
  write ( *, '(a)' ) ' '
  x_r8 = 1.0D+00
  y_r8 = huge ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8
  x_r8 = 0.0D+00
  y_r8 = huge ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8
  x_r8 = 1000000.0D+00
  y_r8 = huge ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now, some integers:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              HUGE(X)'
  write ( *, '(a)' ) ' '
  x_i4 = 1
  y_i4 = huge ( x_i4 )
  write ( *, '(2x,i12,2x,i12)' ) x_i4, y_i4
  x_i4 = 0.0D+00
  y_i4 = huge ( x_i4 )
  write ( *, '(2x,i12,2x,i12)' ) x_i4, y_i4
  x_i4 = 1000000.0D+00
  y_i4 = huge ( x_i4 )
  write ( *, '(2x,i12,2x,i12)' ) x_i4, y_i4

  return
end
subroutine test_iachar ( )

!*****************************************************************************80
!
!! TEST_IACHAR tests IACHAR.
!
!  Discussion:
!
!    The FORTRAN90 function IACHAR returns the ASCII index (between 0 and 255)
!    of a character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character c1
  character c2
  integer i
  integer i1
  character ( len = 25 ) :: string = 'This is a string of text!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_IACHAR'
  write ( *, '(a)' ) '  IACHAR is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  ASCII index of a given character'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C        IACHAR(C)    ACHAR(IACHAR(C))'
  write ( *, '(a)' ) ' '

  do i = 1, len ( string )

    c1 = string (i:i)

    i1 = iachar ( c1 )

    c2 = achar ( i1 )

    write ( *, '(2x,a1,8x,i8,8x,a1)' ) c1, i1, c2

  end do

  return
end
subroutine test_iand_i4 ( )

!*****************************************************************************80
!
!! TEST_IAND_I4 tests IAND on integer ( kind = 4 ) arguments.
!
!  Discussion:
!
!    The FORTRAN90 function IAND returns the bitwise AND 
!    of two integers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer i4_uniform
  integer j
  integer k
  integer seed
  integer test

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_IAND_I4'
  write ( *, '(a)' ) '  IAND is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  bitwise AND of two integers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, I and J are integers of KIND = 4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J    IAND(I,J)'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    i = i4_uniform ( 0, 100, seed )
    j = i4_uniform ( 0, 100, seed )
    k = iand ( i, j )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
  end do

  return
end
subroutine test_iand_i8 ( )

!*****************************************************************************80
!
!! TEST_IAND_I8 tests IAND on integer ( kind = 8 ) arguments.
!
!  Discussion:
!
!    The FORTRAN90 function IAND returns the bitwise exclusive AND 
!    of two integers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 8 ) i
  integer i4_uniform
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  integer seed
  integer test

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_IAND_I8'
  write ( *, '(a)' ) '  IAND is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  bitwise AND of two integers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, I and J are integers of KIND = 8.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J    IAND(I,J)'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    i = i4_uniform ( 0, 100, seed )
    j = i4_uniform ( 0, 100, seed )
    k = iand ( i, j )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
  end do

  return
end
subroutine test_ibclr ( )

!*****************************************************************************80
!
!! TEST_IBCLR tests IBCLR.
!
!  Discussion:
!
!    The FORTRAN90 function IBCLR sets a given bit to zero in an integer word.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i1
  integer i2
  integer pos
!
!  Put 11 consecutive 1's into I1.
!
  i1 = 0
  do pos = 0, 10
    i1 = 2 * i1 + 1
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_IBCLR'
  write ( *, '(a)' ) '  IBCLR is a FORTRAN90 function which sets a given'
  write ( *, '(a)' ) '  bit to zero in an integer word.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I       POS    IBCLR(I,POS)'
  write ( *, '(a)' ) ' '
  do pos = 0, 10
    i2 = ibclr ( i1, pos )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i1, pos, i2
  end do

  return
end
subroutine test_ibits ( )

!*****************************************************************************80
!
!! TEST_IBITS tests IBITS.
!
!  Discussion:
!
!    The FORTRAN90 function IBITS extracts a sequence of bits from
!    an integer word.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer byte
  integer i
  integer i1
  integer i2
  integer ( kind = 4 ) i4
  integer len
  integer pos

  i1 = 1396

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_IBITS'
  write ( *, '(a)' ) '  IBITS is a FORTRAN90 function which extracts'
  write ( *, '(a)' ) '  LEN bits from word I start at position POS.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I       POS    LEN    IBITS(I,POS,LEN)'
  write ( *, '(a)' ) ' '
  len = 3
  do pos = 0, 10
    i2 = ibits ( i1, pos, len )
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) i1, pos, len, i2
  end do

  write ( *, '(a)' ) ' '

  pos = 2
  do len = 1, 10
    i2 = ibits ( i1, pos, len )
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) i1, pos, len, i2
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use IBITS to extract the 4 bytes that make up'
  write ( *, '(a)' ) '  an integer ( kind = 4 ) word.'
  write ( *, '(a)' ) ' '
!
!  3: 00111110 =  62
!
!  2: 00000100 =   4
!
!  1: 11010010 = 210
!
!  0: 00001111 =  15
!
  i4 = 2**29 + 2**28 + 2**27 + 2**26 + 2**25 &
     + 2**18 &
     + 2**15 + 2**14 + 2**12 + 2**9 &
     + 2**3 + 2**2 + 2**1 + 2**0 

  write ( *, '(a,i12)' ) '  I4 = ', i4
  write ( *, '(a)' ) ' '
  len = 8
  do i = 0, 3
    pos = i * len
    byte = ibits ( i4, pos, len )
    write ( *, '(a,i8,a,i8)' ) '  Byte ', i, ' = ', byte
  end do

  return
end
subroutine test_ibset ( )

!*****************************************************************************80
!
!! TEST_IBSET tests IBSET.
!
!  Discussion:
!
!    The FORTRAN90 function IBSET sets a given bit to one in an integer word.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i1
  integer i2
  integer pos

  i1 = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_IBSET'
  write ( *, '(a)' ) '  IBSET is a FORTRAN90 function which sets a given'
  write ( *, '(a)' ) '  bit to one in an integer word.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I       POS    IBSET(I,POS)'
  write ( *, '(a)' ) ' '
  do pos = 0, 10
    i2 = ibset ( i1, pos )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i1, pos, i2
    i1 = i2
  end do

  return
end
subroutine test_ichar ( )

!*****************************************************************************80
!
!! TEST_ICHAR tests ICHAR.
!
!  Discussion:
!
!    The FORTRAN90 function ICHAR returns the character index (between 0 
!    and 255) of a character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character c1
  character c2
  integer i
  integer i1
  character ( len = 25 ) :: string = 'This is a string of text!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ICHAR'
  write ( *, '(a)' ) '  ICHAR is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  character index of a given character'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C        ICHAR(C)    CHAR(ICHAR(C))'
  write ( *, '(a)' ) ' '

  do i = 1, len ( string )

    c1 = string (i:i)

    i1 = ichar ( c1 )

    c2 = char ( i1 )

    write ( *, '(2x,a1,8x,i8,8x,a1)' ) c1, i1, c2

  end do

  return
end
subroutine test_ieor_i4 ( )

!*****************************************************************************80
!
!! TEST_IEOR_I4 tests IEOR on integer ( kind = 4 ) arguments.
!
!  Discussion:
!
!    The FORTRAN90 function IEOR returns the bitwise exclusive OR 
!    of two integers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_IEOR_I4'
  write ( *, '(a)' ) '  IEOR is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  bitwise exclusive OR of two integers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, I and J are integers of KIND = 4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J    IEOR(I,J)'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    call random_number ( harvest = r )
    i = int ( 100.0E+00 * r )
    call random_number ( harvest = r )
    j = int ( 100.0E+00 * r )
    k = ieor ( i, j )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
  end do

  return
end
subroutine test_ieor_i8 ( )

!*****************************************************************************80
!
!! TEST_IEOR_I8 tests IEOR on integer ( kind = 8 ) arguments.
!
!  Discussion:
!
!    The FORTRAN90 function IEOR returns the bitwise exclusive OR 
!    of two integers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 8 ) i
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  real r
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_IEOR_I8'
  write ( *, '(a)' ) '  IEOR is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  bitwise exclusive OR of two integers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, I and J are integers of KIND = 8.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J    IEOR(I,J)'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    call random_number ( harvest = r )
    i = int ( 100.0 * r )
    call random_number ( harvest = r )
    j = int ( 100.0 * r )
    k = ieor ( i, j )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
  end do

  return
end
subroutine test_index ( )

!*****************************************************************************80
!
!! TEST_INDEX tests INDEX.
!
!  Discussion:
!
!    The FORTRAN90 function INDEX determines the first occurrence
!    of a substring in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INDEX'
  write ( *, '(a)' ) '  INDEX(S,SUB) is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the location of the first occurrence of '
  write ( *, '(a)' ) '  substring SUB in string S.'
  write ( *, '(a)' ) '  INDEX(S,SUB,.TRUE.) returns the location of the'
  write ( *, '(a)' ) '  LAST occurrence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  index ( ''THE CATATONIC CAT'', ''CAT'' )', &
    index ( 'THE CATATONIC CAT', 'CAT' )
  write ( *, '(a,i8)' ) '  index ( ''THE CATATONIC CAT'', ''cat'' )', &
    index ( 'THE CATATONIC CAT', 'cat' )
  write ( *, '(a,i8)' ) '  index ( ''THE CATATONIC CAT'', ''CAT'', .TRUE. )', &
    index ( 'THE CATATONIC CAT', 'CAT', .TRUE. )

  return
end
subroutine test_int ( )

!*****************************************************************************80
!
!! TEST_INT tests INT.
!
!  Discussion:
!
!    The FORTRAN90 function INT converts a numeric value to an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) :: x_c4 = ( 1.1E+00, 2.2E+00 )
  complex ( kind = 8 ) :: x_c8 = ( 3.3D+00, 4.4D+00 )
  integer ( kind = 4 ) :: x_i4 = 5
  integer ( kind = 8 ) :: x_i8 = 6
  real ( kind = 4 ) :: x_r4 = 7.7E+00
  real ( kind = 8 ) :: x_r8 = 8.8D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INT'
  write ( *, '(a)' ) '  INT is a FORTRAN90 function which converts'
  write ( *, '(a)' ) '  a complex, integer or real value to integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Type                   X             INT(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f6.4,2x,f6.4,2x,i6)' ) '  complex ( kind = 4 ) ', &
    x_c4, int ( x_c4 )
  write ( *, '(a,f6.4,2x,f6.4,2x,i6)' ) '  complex ( kind = 8 ) ', &
    x_c8, int ( x_c8 )
  write ( *, '(a,i6,10x,i6)'    ) '  integer ( kind = 4 ) ', &
    x_i4, int ( x_i4 )
  write ( *, '(a,i6,10x,i6)'    ) '  integer ( kind = 8 ) ', &
    x_i8, int ( x_i8 )
  write ( *, '(a,f6.4,10x,i6)'  ) '  real ( kind = 4 )    ', &
    x_r4, int ( x_r4 )
  write ( *, '(a,f6.4,10x,i6)'  ) '  real ( kind = 8 )    ', &
    x_r8, int ( x_r8 )

  return
end
subroutine test_ior_i4 ( )

!*****************************************************************************80
!
!! TEST_IOR_I4 tests IOR on integer ( kind = 4 ) arguments.
!
!  Discussion:
!
!    The FORTRAN90 function IOR returns the bitwise inclusive OR 
!    of two integers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer i4_uniform
  integer j
  integer k
  integer seed
  integer test

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_IOR_I4'
  write ( *, '(a)' ) '  IOR is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  bitwise inclusive OR of two integers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, I and J are integers of KIND = 4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J     IOR(I,J)'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    i = i4_uniform ( 0, 100, seed )
    j = i4_uniform ( 0, 100, seed )
    k = ior ( i, j )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
  end do

  return
end
subroutine test_ior_i8 ( )

!*****************************************************************************80
!
!! TEST_IOR_I8 tests IOR on integer ( kind = 8 ) arguments.
!
!  Discussion:
!
!    The FORTRAN90 function IEOR returns the bitwise inclusive OR 
!    of two integers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 8 ) i
  integer ( kind = 8 ) i8_0
  integer ( kind = 8 ) i8_100
  integer ( kind = 8 ) i8_uniform
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  integer ( kind = 8 ) seed
  integer ( kind = 8 ) test

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_IOR_I8'
  write ( *, '(a)' ) '  IOR is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  bitwise inclusive OR of two integers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, I and J are integers of KIND = 8.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J     IOR(I,J)'
  write ( *, '(a)' ) ' '
!
!  You may have trouble passing integer ( kind = 8 ) constants as arguments;
!  You either have to "tag" them with their type, or store them in variables.
!
  i8_0 = 0
  i8_100 = 100

  do test = 1, 10
    i = i8_uniform ( i8_0, i8_100, seed )
    j = i8_uniform ( i8_0, i8_100, seed )
    k = ior ( i, j )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
  end do

  return
end
subroutine test_ishft ( )

!*****************************************************************************80
!
!! TEST_ISHFT tests ISHFT.
!
!  Discussion:
!
!    The FORTRAN90 function ISHFT shifts the bits in an integer word.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i1
  integer i2
  integer shift

  i1 = 89

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ISHFT'
  write ( *, '(a)' ) '  ISHFT is a FORTRAN90 function which shifts'
  write ( *, '(a)' ) '  the bits in an integer word.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I       SHIFT    ISHFT(I,SHIFT)'
  write ( *, '(a)' ) ' '
  do shift = -5, 5
    i2 = ishft ( i1, shift )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i1, shift, i2
  end do

  return
end
subroutine test_ishftc ( )

!*****************************************************************************80
!
!! TEST_ISHFTC tests ISHFTC.
!
!  Discussion:
!
!    The FORTRAN90 function ISHFTC circular-shifts the bits in an integer word.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i1
  integer i2
  integer shift

  i1 = 89

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ISHFTC'
  write ( *, '(a)' ) '  ISHFTC is a FORTRAN90 function which circular-shifts'
  write ( *, '(a)' ) '  the bits in an integer word.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I       SHIFT    ISHFTC(I,SHIFT)'
  write ( *, '(a)' ) ' '
  do shift = -5, 5
    i2 = ishftc ( i1, shift )
    write ( *, '(2x,i8,2x,i8,2x,i12)' ) i1, shift, i2
  end do

  return
end
subroutine test_kind ( )

!*****************************************************************************80
!
!! TEST_KIND tests KIND.
!
!  Discussion:
!
!    The FORTRAN90 function KIND returns the "kind" associated with
!    a variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character x_ch
  complex x_c
  complex ( kind = 4 ) x_c4
  complex ( kind = 8 ) x_c8
  integer x_i
  integer ( kind = 4 ) x_i4
  integer ( kind = 8 ) x_i8
  logical x_l
  real x_r
  real ( kind = 4 ) x_r4
  real ( kind = 8 ) x_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_KIND'
  write ( *, '(a)' ) '  KIND is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the kind associated with a given variable'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Declaration              KIND(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  character              ', kind ( x_ch )
  write ( *, '(a,i8)' ) '  complex                ', kind ( x_c )
  write ( *, '(a,i8)' ) '  complex ( kind = 4 )   ', kind ( x_c4 )
  write ( *, '(a,i8)' ) '  complex ( kind = 8 )   ', kind ( x_c8 )
  write ( *, '(a,i8)' ) '  integer                ', kind ( x_i )
  write ( *, '(a,i8)' ) '  integer ( kind = 4 )   ', kind ( x_i4 )
  write ( *, '(a,i8)' ) '  integer ( kind = 8 )   ', kind ( x_i8 )
  write ( *, '(a,i8)' ) '  logical                ', kind ( x_l )
  write ( *, '(a,i8)' ) '  real                   ', kind ( x_r )
  write ( *, '(a,i8)' ) '  real ( kind = 4 )      ', kind ( x_r4 )
  write ( *, '(a,i8)' ) '  real ( kind = 8 )      ', kind ( x_r8 )

  return
end
subroutine test_lbound ( )

!*****************************************************************************80
!
!! TEST_LBOUND tests LBOUND.
!
!  Discussion:
!
!    The FORTRAN90 function LBOUND(ARRAY) returns the lower bounds
!    of all the dimensions; LBOUND(ARRAY,DIM) returns the lower bound
!    in the given dimension.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real a(5,10,17)
  integer b(4:6,-5:-1,10:20)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LBOUND'
  write ( *, '(a)' ) '  LBOUND(ARRAY) is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  lower array bounds in all dimensions;'
  write ( *, '(a)' ) '  LBOUND(ARRAY,DIM) returns the lower bound in dimension'
  write ( *, '(a)' ) '  DIM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  real a(5,10,17)'
  write ( *, '(a,3i6)' ) '  lbound(a)   = ', lbound(a)
  write ( *, '(a,3i6)' ) '  lbound(a,1) = ', lbound(a,1)
  write ( *, '(a,3i6)' ) '  lbound(a,2) = ', lbound(a,2)
  write ( *, '(a,3i6)' ) '  lbound(a,3) = ', lbound(a,3)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer b(4:6,-5:-1,10:20)'
  write ( *, '(a,3i6)' ) '  lbound(b)   = ', lbound(b)
  write ( *, '(a,3i6)' ) '  lbound(b,1) = ', lbound(b,1)
  write ( *, '(a,3i6)' ) '  lbound(b,2) = ', lbound(b,2)
  write ( *, '(a,3i6)' ) '  lbound(b,3) = ', lbound(b,3)

  return
end
subroutine test_len ( )

!*****************************************************************************80
!
!! TEST_LEN tests LEN.
!
!  Discussion:
!
!    The FORTRAN90 function LEN returns the declared length of a
!    character string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 1 ) s1
  character ( len = 2 ) s2
  character ( len = 4 ) s4
  character ( len = 8 ) s8
  character ( len = 16 ) s16

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LEN'
  write ( *, '(a)' ) '  LEN is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  declared length of a string variable, or the length of'
  write ( *, '(a)' ) '  a string constant'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      S                     LEN(S)'
  write ( *, '(a)' ) '   ----------               -----'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  character ( len =  1 ) ', len ( s1 )
  write ( *, '(a,i8)' ) '  character ( len =  2 ) ', len ( s2 )
  write ( *, '(a,i8)' ) '  character ( len =  4 ) ', len ( s4 )
  write ( *, '(a,i8)' ) '  character ( len =  8 ) ', len ( s8 )
  write ( *, '(a,i8)' ) '  character ( len = 16 ) ', len ( s16 )
  write ( *, '(a,i8)' ) ' "A STRING"              ', len ( 'A STRING' )

  return
end
subroutine test_len_trim ( )

!*****************************************************************************80
!
!! TEST_LEN_TRIM tests LEN_TRIM.
!
!  Discussion:
!
!    The FORTRAN90 function LEN_TRIM returns the "used" length of a
!    character string, up to the last non_blank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LEN_TRIM'
  write ( *, '(a)' ) '  LEN_TRIM is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  "used" length of a string variable up to the last'
  write ( *, '(a)' ) '  nonblank.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      S          LEN_TRIM(S)'
  write ( *, '(a)' ) '   ----------    ----------'
  write ( *, '(a)' ) ' '
  s = '1234567890'
  write ( *, '(a,i8)' ) '  "' // s // '"', len_trim(s)
  s = '12345     '
  write ( *, '(a,i8)' ) '  "' // s // '"', len_trim(s)
  s = '     67890'
  write ( *, '(a,i8)' ) '  "' // s // '"', len_trim(s)
  s = '    5     '
  write ( *, '(a,i8)' ) '  "' // s // '"', len_trim(s)
  s = '1 3 5 7 9 '
  write ( *, '(a,i8)' ) '  "' // s // '"', len_trim(s)

  return
end
subroutine test_lge ( )

!*****************************************************************************80
!
!! TEST_LGE tests LGE.
!
!  Discussion:
!
!    The FORTRAN90 function LGE(S1,S2) returns the value of
!    "string S1 is lexically greater than or equal to string S2".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 3 ) :: s4 = 'boy' 
  character ( len = 3 ) :: s5 = 'cat'
  character ( len = 4 ) :: s6 = 'cats'
  character ( len = 3 ) :: s7 = 'dog'
  character ( len = 3 ) :: s8 = 'CAT'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LGE'
  write ( *, '(a)' ) '  LGE is a FORTRAN90 function which returns the value'
  write ( *, '(a)' ) '  of "S1 >= S2" where S1 and S2 are strings.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    S1    S2   LGE(S1,S2)'
  write ( *, '(a)' ) '   ---   ---   ----------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,8x,l1)' ) '  "' // s4 // '"  "' // s4 // '"  ', lge ( s4, s4 )
  write ( *, '(a,8x,l1)' ) '  "' // s4 // '"  "' // s5 // '"  ', lge ( s4, s5 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s4 // '"  ', lge ( s5, s4 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s7 // '"  ', lge ( s5, s7 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s8 // '"  ', lge ( s5, s8 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s6 //  '" ', lge ( s5, s6 )

  return
end
subroutine test_lgt ( )

!*****************************************************************************80
!
!! TEST_LGT tests LGT.
!
!  Discussion:
!
!    The FORTRAN90 function LGT(S1,S2) returns the value of
!    "string S1 is lexically greater than string S2".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 3 ) :: s4 = 'boy' 
  character ( len = 3 ) :: s5 = 'cat'
  character ( len = 4 ) :: s6 = 'cats'
  character ( len = 3 ) :: s7 = 'dog'
  character ( len = 3 ) :: s8 = 'CAT'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LGT'
  write ( *, '(a)' ) '  LGT is a FORTRAN90 function which returns the value'
  write ( *, '(a)' ) '  of "S1 > S2" where S1 and S2 are strings.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    S1    S2   LGT(S1,S2)'
  write ( *, '(a)' ) '   ---   ---   ----------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,8x,l1)' ) '  "' // s4 // '"  "' // s4 // '"  ', lgt ( s4, s4 )
  write ( *, '(a,8x,l1)' ) '  "' // s4 // '"  "' // s5 // '"  ', lgt ( s4, s5 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s4 // '"  ', lgt ( s5, s4 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s7 // '"  ', lgt ( s5, s7 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s8 // '"  ', lgt ( s5, s8 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s6 //  '" ', lgt ( s5, s6 )

  return
end
subroutine test_lle ( )

!*****************************************************************************80
!
!! TEST_LLE tests LLE.
!
!  Discussion:
!
!    The FORTRAN90 function LLE(S1,S2) returns the value of
!    "string S1 is lexically less than or equal to string S2".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 3 ) :: s4 = 'boy' 
  character ( len = 3 ) :: s5 = 'cat'
  character ( len = 4 ) :: s6 = 'cats'
  character ( len = 3 ) :: s7 = 'dog'
  character ( len = 3 ) :: s8 = 'CAT'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LLE'
  write ( *, '(a)' ) '  LLE is a FORTRAN90 function which returns the value'
  write ( *, '(a)' ) '  of "S1 <= S2" where S1 and S2 are strings.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    S1    S2   LLE(S1,S2)'
  write ( *, '(a)' ) '   ---   ---   ----------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,8x,l1)' ) '  "' // s4 // '"  "' // s4 // '"  ', lle ( s4, s4 )
  write ( *, '(a,8x,l1)' ) '  "' // s4 // '"  "' // s5 // '"  ', lle ( s4, s5 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s4 // '"  ', lle ( s5, s4 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s7 // '"  ', lle ( s5, s7 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s8 // '"  ', lle ( s5, s8 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s6 //  '" ', lle ( s5, s6 )

  return
end
subroutine test_llt ( )

!*****************************************************************************80
!
!! TEST_LLT tests LLT.
!
!  Discussion:
!
!    The FORTRAN90 function LLT(S1,S2) returns the value of
!    "string S1 is lexically less than string S2".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 3 ) :: s4 = 'boy' 
  character ( len = 3 ) :: s5 = 'cat'
  character ( len = 4 ) :: s6 = 'cats'
  character ( len = 3 ) :: s7 = 'dog'
  character ( len = 3 ) :: s8 = 'CAT'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LGT'
  write ( *, '(a)' ) '  LLT is a FORTRAN90 function which returns the value'
  write ( *, '(a)' ) '  of "S1 < S2" where S1 and S2 are strings.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    S1    S2   LLT(S1,S2)'
  write ( *, '(a)' ) '   ---   ---   ----------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,8x,l1)' ) '  "' // s4 // '"  "' // s4 // '"  ', llt ( s4, s4 )
  write ( *, '(a,8x,l1)' ) '  "' // s4 // '"  "' // s5 // '"  ', llt ( s4, s5 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s4 // '"  ', llt ( s5, s4 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s7 // '"  ', llt ( s5, s7 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s8 // '"  ', llt ( s5, s8 )
  write ( *, '(a,8x,l1)' ) '  "' // s5 // '"  "' // s6 //  '" ', llt ( s5, s6 )

  return
end
subroutine test_log ( )

!*****************************************************************************80
!
!! TEST_LOG tests LOG.
!
!  Discussion:
!
!    The FORTRAN90 function LOG returns the natural logarithm of a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo =  0.0D+00
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LOG'
  write ( *, '(a)' ) '  LOG is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  natural logarithm of a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              LOG(X)     EXP(LOG(X))'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = log ( x )
    z = exp ( y )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, y, z
  end do

  return
end
subroutine test_log10 ( )

!*****************************************************************************80
!
!! TEST_LOG10 tests LOG10.
!
!  Discussion:
!
!    The FORTRAN90 function LOG10 returns the base 10 logarithm of a 
!    real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo =  0.0D+00
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LOG10'
  write ( *, '(a)' ) '  LOG10 is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  base 10 logarithm of a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              LOG10(X)     10**(LOG(X))'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = log10 ( x )
    z = 10.0D+00**y
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, y, z
  end do

  return
end
subroutine test_logical ( )

!*****************************************************************************80
!
!! TEST_LOGICAL tests LOGICAL.
!
!  Discussion:
!
!    The FORTRAN90 function LOGICAL can convert between logical kinds.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  logical ( kind = 4 ) x_l4
  logical ( kind = 1 ) x_lb

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LOGICAL'
  write ( *, '(a)' ) '  LOGICAL is a FORTRAN90 function which can convert'
  write ( *, '(a)' ) '  between logical kinds'

  x_l4 = .true.

  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  logical ( kind = 4 ) x_l4 = ', x_l4

  x_lb = logical ( x_l4, kind = 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  logical ( kind = 1 ) x_lb'
  write ( *, '(a,l1)' ) '  x_lb = logical ( x_l4, kind = 1 ) = ', x_lb

  x_lb = .not. x_lb
  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  x_lb = .not. x_lb = ', x_lb

  x_l4 = logical ( x_lb, kind = 4 )
  write ( *, '(a,l1)' ) '  x_l4 = logical ( x_lb, kind = 4 ) = ', x_l4

  return
end
subroutine test_matmul ( )

!*****************************************************************************80
!
!! TEST_MATMUL tests MATMUL.
!
!  Discussion:
!
!    The FORTRAN90 function MATMUL multiplies two matrices or a matrix
!    and a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a(4,5)
  integer ( kind = 4 ) ab(4,3)
  integer ( kind = 4 ) ac(4)
  integer ( kind = 4 ) b(5,3)
  integer ( kind = 4 ) c(5)
  integer ( kind = 4 ) cb(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) :: seed = 123456789

  do i = 1, 4
    do j = 1, 5
      a(i,j) = i4_uniform ( 0, 5, seed )
    end do
  end do

  do i = 1, 5
    do j = 1, 3
      b(i,j) = i4_uniform ( 0, 5, seed )
    end do
  end do

  do i = 1, 5
    c(i) = i4_uniform ( 0, 5, seed )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MATMUL'
  write ( *, '(a)' ) '  MATMUL is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the product of two matrices or a matrix and vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix A:'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(5i6)' ) a(i,1:5)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix B:'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(5i6)' ) a(i,1:3)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vector C:'
  write ( *, '(a)' ) ' '
  write ( *, '(5i6)' ) c(1:5)

  ab = matmul ( a, b )
  ac = matmul ( a, c )
  cb = matmul ( c, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix AB = matmul ( A, B ):'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(5i6)' ) ab(i,1:3)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vector AC = matmul ( A, C ):'
  write ( *, '(a)' ) ' '
  write ( *, '(5i6)' ) ac(1:4)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vector CB = matmul ( C, B ):'
  write ( *, '(a)' ) ' '
  write ( *, '(5i6)' ) cb(1:3)

  return
end
subroutine test_max ( )

!*****************************************************************************80
!
!! TEST_MAX tests MAX.
!
!  Discussion:
!
!    The FORTRAN90 function MAX returns the maximum value in a list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MAX'
  write ( *, '(a)' ) '  MAX is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the maximum value in a list.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  max(2,1) =     ', max ( 2, 1 )
  write ( *, '(a,i8)' ) '  max(1,3,2) =   ', max ( 1, 3, 2 )
  write ( *, '(a,i8)' ) '  max(3,2,4,1) = ', max ( 3, 2, 4, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f4.1)' ) '  max(2.1, 1.2) =           ', max ( 2.1, 1.2 )
  write ( *, '(a,f4.1)' ) '  max(1.1, 3.2, 2.3) =      ', max ( 1.1, 3.2, 2.3 )
  write ( *, '(a,f4.1)' ) &
    '  max(3.1, 2.2, 4.3, 1.4) = ', max ( 3.1, 2.2, 4.3, 1.4 )

  return
end
subroutine test_max_vector ( )

!*****************************************************************************80
!
!! TEST_MAX_VECTOR tests MAX applied to a vector.
!
!  Discussion:
!
!    The FORTRAN90 function MAX returns the maximum value in a list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) y(n)
  integer ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MAX_VECTOR'
  write ( *, '(a)' ) '  MAX is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the maximum value in a list.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We consider a vector version:'
  write ( *, '(a)' ) ' '

  seed = 123456789
  call i4vec_uniform ( n, 0, 10, seed, x )
  call i4vec_uniform ( n, 0, 10, seed, y )

  z(1:n) = max ( x(1:n),  y(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I      X(I)      Y(I)      Z(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) i, x(i), y(i), z(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now try Z(1:N) = max ( 5, Y(1:N) )'

  seed = 123456789
  call i4vec_uniform ( n, 0, 10, seed, x )
  call i4vec_uniform ( n, 0, 10, seed, y )

  z(1:n) = max ( 5,  y(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I      X(I)      Y(I)      Z(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) i, 5, y(i), z(i)
  end do
  return
end
subroutine test_maxexponent ( )

!*****************************************************************************80
!
!! TEST_MAXEXPONENT tests MAXEXPONENT.
!
!  Discussion:
!
!    The FORTRAN90 function MAXEXPONENT returns the maximum exponent
!    associated with real numbers of the given kind.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real x_r
  real ( kind = 4 ) x_r4
  real ( kind = 8 ) x_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MAXEXPONENT'
  write ( *, '(a)' ) '  MAXEXPONENT is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the maximum exponent associated with real numbers'
  write ( *, '(a)' ) '  of the same kind as X.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Type         MAXEXPONENT(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  real              ', maxexponent ( x_r )
  write ( *, '(a,i8)' ) '  real ( kind = 4 ) ', maxexponent ( x_r4 )
  write ( *, '(a,i8)' ) '  real ( kind = 8 ) ', maxexponent ( x_r8 )

  return
end
subroutine test_maxloc ( )

!*****************************************************************************80
!
!! TEST_MAXLOC tests MAXLOC.
!
!  Discussion:
!
!    The FORTRAN90 function MAXLOC returns the index of the maximum value of
!    the entries of an integer or real vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) x_i4(n)
  real ( kind = 4 ) x_r4(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MAXLOC'
  write ( *, '(a)' ) '  MAXLOC is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  index of the maximum value of the entries of a vector.'

  call random_number ( harvest = x_r4(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A real ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,f14.6)' ) x_r4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,i14)' ) '  MAXLOC(X) = ', maxloc ( x_r4 )

  do i = 1, n
    x_i4(i) = int ( 10.0E+00 * x_r4(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  An integer ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i14)' ) x_i4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,i14)' ) '  MAXLOC(X) = ', maxloc ( x_i4 )

  return
end
subroutine test_maxval ( )

!*****************************************************************************80
!
!! TEST_MAXVAL tests MAXVAL.
!
!  Discussion:
!
!    The FORTRAN90 function MAXVAL returns the maximum value of
!    the entries of an integer or real vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) x_i4(n)
  real ( kind = 4 ) x_r4(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MAXVAL'
  write ( *, '(a)' ) '  MAXVAL is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  maximum value of the entries of a vector.'

  call random_number ( harvest = x_r4(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A real ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,f14.6)' ) x_r4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.6)' ) '  MAXVAL(X) = ', maxval ( x_r4 )

  do i = 1, n
    x_i4(i) = int ( 10.0 * x_r4(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  An integer ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i14)' ) x_i4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,i14)' ) '  MAXVAL(X) = ', maxval ( x_i4 )

  return
end
subroutine test_merge ( )

!*****************************************************************************80
!
!! TEST_MERGE tests MERGE.
!
!  Discussion:
!
!    The FORTRAN90 function MERGE copies values from one array or the
!    other depending on a MASK value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer a(4,5)
  integer b(4,5)
  logical c(4,5)
  integer d(4,5)
  integer i
  integer i4_uniform
  integer j
  integer :: seed = 123456789

  do i = 1, 4
    do j = 1, 5
      a(i,j) = i4_uniform ( 0, 100, seed )
    end do
  end do

  do i = 1, 4
    do j = 1, 5
      b(i,j) = i4_uniform ( 0, 100, seed )
    end do
  end do

  do i = 1, 4
    do j = 1, 5
      c(i,j) = ( j == 3 .or. mod ( i, 2 ) == 0 )
    end do
  end do

  d = merge ( a, b, c )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MERGE'
  write ( *, '(a)' ) '  MERGE is a FORTRAN90 function which copies'
  write ( *, '(a)' ) '  entries from one array or the other depending'
  write ( *, '(a)' ) '  on a logical array.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Array A:'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(5i6)' ) a(i,1:5)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Array B:'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(5i6)' ) b(i,1:5)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Array C:'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(5(5x,l1))' ) c(i,1:5)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Array D = MERGE(A,B,C):'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(5i6)' ) d(i,1:5)
  end do

  return
end
subroutine test_min ( )

!*****************************************************************************80
!
!! TEST_MIN tests MIN.
!
!  Discussion:
!
!    The FORTRAN90 function MIN returns the minimum value in a list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MIN'
  write ( *, '(a)' ) '  MIN is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the minimum value in a list.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  min(3, 4) =       ', min ( 3, 4 )
  write ( *, '(a,i8)' ) '  min(4, 2, 3) =    ', min ( 4, 2, 3 )
  write ( *, '(a,i8)' ) '  min(2, 3, 1, 4) = ', min ( 2, 3, 1, 4 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f4.1)' ) '  min(3.1, 4.2) =           ', min ( 3.1, 4.2 )
  write ( *, '(a,f4.1)' ) '  min(4.1. 2.2, 3.3) =      ', min ( 4.1, 2.2, 3.3 )
  write ( *, '(a,f4.1)' ) &
    '  min(2.1, 3.2, 1.3, 4.4) = ', min ( 2.1, 3.2, 1.3, 4.4 )

  return
end
subroutine test_minexponent ( )

!*****************************************************************************80
!
!! TEST_MINEXPONENT tests MINEXPONENT.
!
!  Discussion:
!
!    The FORTRAN90 function MINXPONENT returns the minimum exponent
!    associated with real numbers of the given kind.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real x_r
  real ( kind = 4 ) x_r4
  real ( kind = 8 ) x_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MINEXPONENT'
  write ( *, '(a)' ) '  MINEXPONENT is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the minimum exponent associated with real numbers'
  write ( *, '(a)' ) '  of the same kind as X.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Type         MINEXPONENT(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  real              ', minexponent ( x_r )
  write ( *, '(a,i8)' ) '  real ( kind = 4 ) ', minexponent ( x_r4 )
  write ( *, '(a,i8)' ) '  real ( kind = 8 ) ', minexponent ( x_r8 )

  return
end
subroutine test_minloc ( )

!*****************************************************************************80
!
!! TEST_MINLOC tests MINLOC.
!
!  Discussion:
!
!    The FORTRAN90 function MINLOC returns the index of the minimum value of
!    the entries of an integer or real vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) x_i4(n)
  real ( kind = 4 ) x_r4(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MINLOC'
  write ( *, '(a)' ) '  MINLOC is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  index of the minimum value of the entries of a vector.'

  call random_number ( harvest = x_r4(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A real ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,f14.6)' ) x_r4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,i14)' ) '  MINLOC(X) = ', minloc ( x_r4 )

  do i = 1, n
    x_i4(i) = int ( 10.0E+00 * x_r4(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  An integer ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i14)' ) x_i4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,i14)' ) '  MINLOC(X) = ', minloc ( x_i4 )

  return
end
subroutine test_minval ( )

!*****************************************************************************80
!
!! TEST_MINVAL tests MINVAL.
!
!  Discussion:
!
!    The FORTRAN90 function MINVAL returns the minimum value of
!    the entries of an integer or real vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 5

  real ( kind = 4 ) r4_uniform_01
  integer i
  integer :: seed = 123456789
  integer ( kind = 4 ) x_i4(n)
  real ( kind = 4 ) x_r4(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MINVAL'
  write ( *, '(a)' ) '  MINVAL is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  minimum value of the entries of a vector.'

  call random_number ( harvest = x_r4(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A real ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,f14.6)' ) x_r4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.6)' ) '  MINVAL(X) = ', minval ( x_r4 )

  do i = 1, n
    x_i4(i) = int ( 10.0 * x_r4(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  An integer ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i14)' ) x_i4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,i14)' ) '  MINVAL(X) = ', minval ( x_i4 )

  return
end
subroutine test_mod_i4 ( )

!*****************************************************************************80
!
!! TEST_MOD_I4 tests MOD on integers.
!
!  Discussion:
!
!    The FORTRAN90 function MOD(A,B) returns the remainder after division.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer :: i_hi = 20
  integer :: i_lo = -10
  integer i4_uniform
  integer j
  integer k
  integer :: seed = 123456789
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MOD_I4'
  write ( *, '(a)' ) '  MOD is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  remainder after division.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, the arguments are integers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J       MOD(I,J)'
  write ( *, '(a)' ) ' '
  do test = 1, 10

    i = i4_uniform ( i_lo, i_hi, seed )
    j = i4_uniform ( i_lo, i_hi, seed )

    if ( j == 0 ) then
      write ( *, '(2x,i8,2x,i8,2x,a)' ) i, j, 'Undefined'
    else
      k = mod ( i, j )
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
    end if

  end do

  return
end
subroutine test_mod_r4 ( )

!*****************************************************************************80
!
!! TEST_MOD_R4 tests MOD on reals.
!
!  Discussion:
!
!    The FORTRAN90 function MOD(A,B) returns the remainder after division.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real r4_uniform
  integer :: seed = 123456789
  integer test
  real x
  real :: x_hi = 20.0E+00
  real :: x_lo = -10.0E+00
  real y
  real z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MOD_R4'
  write ( *, '(a)' ) '  MOD is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  remainder after division.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, the arguments are reals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X               Y              MOD(X,Y)'
  write ( *, '(a)' ) ' '
  do test = 1, 10
    x = r4_uniform ( x_lo, x_hi, seed )
    y = r4_uniform ( x_lo, x_hi, seed )
    z = mod ( x, y )
    write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6)' ) x, y, z
  end do

  return
end
subroutine test_modulo_i4 ( )

!*****************************************************************************80
!
!! TEST_MODULO_I4 tests MODULO on integers.
!
!  Discussion:
!
!    The FORTRAN90 function MODULO(A,B) returns the remainder after division.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer :: i_hi = 20
  integer :: i_lo = -10
  integer i4_uniform
  integer j
  integer k
  integer :: seed = 123456789
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MODULO_I4'
  write ( *, '(a)' ) '  MODULO is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  remainder after division.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, the arguments are integers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J    MODULO(I,J)'
  write ( *, '(a)' ) ' '

  do test = 1, 10

    i = i4_uniform ( i_lo, i_hi, seed )
    j = i4_uniform ( i_lo, i_hi, seed )

    if ( j == 0 ) then
      write ( *, '(2x,i8,2x,i8,2x,a)' ) i, j, 'Undefined'
    else
      k = modulo ( i, j )
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
    end if

  end do

  return
end
subroutine test_modulo_r4 ( )

!*****************************************************************************80
!
!! TEST_MODULO_R4 tests MODULO on reals.
!
!  Discussion:
!
!    The FORTRAN90 function MODULO(A,B) returns the remainder after division.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real r4_uniform
  integer :: seed = 123456789
  integer test
  real x
  real :: x_hi = 20.0
  real :: x_lo = -10.0
  real y
  real z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MODULO_R4'
  write ( *, '(a)' ) '  MODULO is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  remainder after division.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, the arguments are reals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X               Y           MODULO(X,Y)'
  write ( *, '(a)' ) ' '
  do test = 1, 10
    x = r4_uniform ( x_lo, x_hi, seed )
    y = r4_uniform ( x_lo, x_hi, seed )
    z = modulo ( x, y )
    write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6)' ) x, y, z
  end do

  return
end
subroutine test_mvbits ( )

!*****************************************************************************80
!
!! TEST_MVBITS tests MVBITS.
!
!  Discussion:
!
!    The FORTRAN90 function MVBITS copies a sequence of bits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2

  i1 = 1396
  i2 = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MVBITS'
  write ( *, '(a)' ) '  MVBITS is a FORTRAN90 function which extracts'
  write ( *, '(a)' ) '  bits from one place and copies them elsewhere.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CALL MVBITS(FROM,FROMPOS,LEN,TO,TOPOS)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  We will always use I1 =        ', i1
  write ( *, '(a,i8)' ) '  We will always start with I2 = ', i2

  write ( *, '(a)' ) ' '
  call mvbits ( i1, 0, 5, i2, 0 )
  write ( *, '(a,i12)' ) '  CALL MVBITS(I1,0, 5,I2,0): I2 = ', i2

  i2 = 0
  call mvbits ( i1, 0, 32, i2, 0 )
  write ( *, '(a,i12)' ) '  CALL MVBITS(I1,0,32,I2,0): I2 = ', i2

  i2 = 0
  call mvbits ( i1, 5, 5, i2, 0 )
  write ( *, '(a,i12)' ) '  CALL MVBITS(I1,5, 5,I2,0): I2 = ', i2

  i2 = 0
  call mvbits ( i1, 5, 5, i2, 5 )
  write ( *, '(a,i12)' ) '  CALL MVBITS(I1,5, 5,I2,5): I2 = ', i2

  return
end
subroutine test_nearest ( )

!*****************************************************************************80
!
!! TEST_NEAREST tests NEAREST.
!
!  Discussion:
!
!    The FORTRAN90 function NEAREST returns the nearest real number to
!    a given real number, in the given direction.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) sm_r4
  real ( kind = 4 ) sp_r4
  real ( kind = 4 ) x_r4
  real ( kind = 4 ) xm_r4
  real ( kind = 4 ) xp_r4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_NEAREST'
  write ( *, '(a)' ) '  NEAREST is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the nearest real number to a given real number, '
  write ( *, '(a)' ) '  in a given direction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X              NEAREST(X,+1.0)   NEAREST(X,-1.0)'
  write ( *, '(a)' ) ' '
  x_r4 = 1.0
  sp_r4 = +1.0
  sm_r4 = -1.0
  xp_r4 = nearest ( x_r4, sp_r4 )
  xm_r4 = nearest ( x_r4, sm_r4 )
  write ( *, '(2x,f16.8,2x,f16.8,2x,f16.8)' ) x_r4, xp_r4, xm_r4

  return
end
subroutine test_nint ( )

!*****************************************************************************80
!
!! TEST_NINT tests NINT.
!
!  Discussion:
!
!    The FORTRAN90 function NINT returns, as an integer, the nearest 
!    integer to a given real value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) j
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  integer ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_NINT'
  write ( *, '(a)' ) '  NINT is a FORTRAN90 function which returns,'
  write ( *, '(a)' ) '  as an integer, the nearest integer to a '
  write ( *, '(a)' ) '  given real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X             NINT(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    j = nint ( x )
    write ( *, '(2x,g14.6,2x,i8)' ) x, j
  end do

  return
end
subroutine test_not_i4 ( )

!*****************************************************************************80
!
!! TEST_NOT_I4 tests NOT on integer ( kind = 4 ) arguments.
!
!  Discussion:
!
!    The FORTRAN90 function NOT returns the bitwise NOT of an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_NOT_I4'
  write ( *, '(a)' ) '  NOT is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  bitwise NOT of an integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, I is an integer of KIND = 4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             I         NOT(I)'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    i = i4_uniform ( 0, 100, seed )
    j = not ( i )
    write ( *, '(2x,i12,2x,i12)' ) i, j
  end do

  return
end
subroutine test_not_i8 ( )

!*****************************************************************************80
!
!! TEST_NOT_I8 tests NOT on integer ( kind = 8 ) arguments.
!
!  Discussion:
!
!    The FORTRAN90 function NOT returns the bitwise NOT of an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 8 ) i
  integer ( kind = 8 ) i8_0
  integer ( kind = 8 ) i8_100
  integer ( kind = 8 ) i8_uniform
  integer ( kind = 8 ) j
  integer ( kind = 8 ) seed
  integer ( kind = 8 ) test

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_NOT_I8'
  write ( *, '(a)' ) '  NOT is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  bitwise NOT of an integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, I is an integer of KIND = 8.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I                  NOT(I)'
  write ( *, '(a)' ) ' '

  i8_0 = 0
  i8_100 = 100

  do test = 1, 10
    i = i8_uniform ( i8_0, i8_100, seed )
    j = not ( i )
    write ( *, '(2x,i20,2x,i20)' ) i, j
  end do

  return
end
subroutine test_pack ( )

!*****************************************************************************80
!
!! TEST_PACK tests PACK.
!
!  Discussion:
!
!    The FORTRAN90 function PACK packs an array into a vector subject
!    to a logical mask.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer a(5,4)
  integer i
  integer i4_uniform
  integer j
  logical mask(5,4)
  integer n
  integer :: seed = 123456789
  integer, allocatable, dimension ( : ) :: v

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_PACK'
  write ( *, '(a)' ) '  PACK is a FORTRAN90 function which packs the'
  write ( *, '(a)' ) '  entries of an array into a vector subject to'
  write ( *, '(a)' ) '  a logical mask.'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    do j = 1, 4
      a(i,j) = i4_uniform ( 0, 100, seed )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here is integer array a(5,4)'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    write ( *, '(4(2x,i4))' ) a(i,1:4)
  end do

  mask(1:5,1:4) = a(1:5,1:4) < 80

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here is the logical array mask(5,4), which'
  write ( *, '(a)' ) '  we will use to choose our entries of A:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    write ( *, '(4(2x,l1))' ) mask(i,1:4)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We can count the number of entries to be'
  write ( *, '(a)' ) '  copied into V by:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    N = COUNT ( MASK )'

  n = count ( mask )

  allocate ( v(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now set V = pack ( a, mask )'
  write ( *, '(a)' ) ' '

  v = pack ( a, mask )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i2,a)' ) '  Here are the entries copied into vector V:'
  write ( *, '(a)' ) ' '

  write ( *, '(10(2x,i4))' ) v(1:n)

  deallocate ( v )

  return
end
subroutine test_precision ( )

!*****************************************************************************80
!
!! TEST_PRECISION tests PRECISION.
!
!  Discussion:
!
!    The FORTRAN90 function PRECISION returns the number of decimal
!    places allocated for the storage of a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  double precision x_d
  real x_r
  real ( kind = 4 ) x_r4
  real ( kind = 8 ) x_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_PRECISION'
  write ( *, '(a)' ) '  PRECISION is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  number of decimal places available for real numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Declaration      PRECISION(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  real             ', precision ( x_r )
  write ( *, '(a,i8)' ) '  double precision ', precision ( x_d )
  write ( *, '(a,i8)' ) '  real ( kind = 4 )', precision ( x_r4 )
  write ( *, '(a,i8)' ) '  real ( kind = 8 )', precision ( x_r8 )

  return
end
subroutine test_present ( a, b, ethel, fred, george )

!*****************************************************************************80
!
!! TEST_PRESENT tests PRESENT.
!
!  Discussion:
!
!    The FORTRAN90 function PRESENT reports whether optional arguments
!    have been supplied.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer a
  integer b
  integer, optional :: ethel
  integer, optional :: fred
  integer, optional :: george

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_PRESENT'
  write ( *, '(a)' ) '  PRESENT is a FORTRAN90 function which  reports'
  write ( *, '(a)' ) '  whether optional arguments have been supplied.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This routine has optional arguments ETHEL, FRED, GEORGE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  present ( ethel ) =  ', present ( ethel )
  write ( *, '(a,l1)' ) '  present ( fred ) =   ', present ( fred )
  write ( *, '(a,l1)' ) '  present ( george ) = ', present ( george )

  return
end
subroutine test_product ( )

!*****************************************************************************80
!
!! TEST_PRODUCT tests PRODUCT.
!
!  Discussion:
!
!    The FORTRAN90 function PRODUCT returns the product of
!    the entries of a complex, integer, or real vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 5

  complex ( kind = 4 ) c4_uniform_01
  integer i
  integer :: seed = 123456789
  complex ( kind = 4 ) x_c4(n)
  integer ( kind = 4 ) x_i4(n)
  real ( kind = 4 ) x_r4(n)

  do i = 1, n
    x_c4(i) = c4_uniform_01 ( seed )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_PRODUCT'
  write ( *, '(a)' ) '  PRODUCT is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  product of the entries of a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A complex ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                 X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,2f14.6)' ) x_c4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,2f14.6)' ) '  PRODUCT(X) = ', product ( x_c4 )

  do i = 1, n
    x_r4(i) = real ( x_c4(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A real ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,f14.6)' ) x_r4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.6)' ) '  PRODUCT(X) = ', product ( x_r4 )

  do i = 1, n
    x_i4(i) = int ( 10.0 * x_r4(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  An integer ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i14)' ) x_i4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,i14)' ) '  PRODUCT(X,Y) = ', product ( x_i4 )

  return
end
subroutine test_radix ( )

!*****************************************************************************80
!
!! TEST_RADIX tests RADIX.
!
!  Discussion:
!
!    The FORTRAN90 function RADIX returns the radix or base exponent
!    associated with the representation of integer or real numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer x_i
  integer ( kind = 4 ) x_i4
  integer ( kind = 8 ) x_i8
  real x_r
  real ( kind = 4 ) x_r4
  real ( kind = 8 ) x_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_RADIX'
  write ( *, '(a)' ) '  RADIX is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  radix or base exponent associated with the'
  write ( *, '(a)' ) '  representation of numbers of a given type and kind.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  integer             ', radix ( x_i )
  write ( *, '(a,i8)' ) '  integer ( kind = 4 )', radix ( x_i4 )
  write ( *, '(a,i8)' ) '  integer ( kind = 8 )', radix ( x_i8 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  real                ', radix ( x_r )
  write ( *, '(a,i8)' ) '  real ( kind = 4 )   ', radix ( x_r4 )
  write ( *, '(a,i8)' ) '  real ( kind = 8 )   ', radix ( x_r8 )

  return
end
subroutine test_random_number ( )

!*****************************************************************************80
!
!! TEST_RANDOM_NUMBER tests RANDOM_NUMBER.
!
!  Discussion:
!
!    The FORTRAN90 subroutine RANDOM_NUMBER returns a scalar or
!    vector or array of pseudorandom values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) a_r4(3,3)
  integer i
  real ( kind = 4 ) s_r4
  real ( kind = 4 ) v_r4(5)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_RANDOM_NUMBER'
  write ( *, '(a)' ) '  RANDOM_NUMBER is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  uniformly distributed pseudorandom values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call RANDOM_NUMBER 5 times in a row:'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    call random_number ( harvest = s_r4 )
    write ( *, '(2x,g14.6)' ) s_r4
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call RANDOM_NUMBER once, with a vector argument'
  write ( *, '(a)' ) ' '
  call random_number ( harvest = v_r4 )
  do i = 1, 5
    write ( *, '(2x,g14.6)' ) v_r4(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call RANDOM_NUMBER with an array argument'
  write ( *, '(a)' ) ' '
  call random_number ( harvest = a_r4 )
  do i = 1, 3
    write ( *, '(2x,3g14.6)' ) a_r4(i,1:3)
  end do

  return
end
subroutine test_random_seed ( )

!*****************************************************************************80
!
!! TEST_RANDOM_SEED tests RANDOM_SEED.
!
!  Discussion:
!
!    The FORTRAN90 subroutine RANDOM_SEED allows the user to read or
!    set the random number seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  real r
  integer, allocatable, dimension ( : ) :: seed
  integer seed_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_RANDOM_SEED'
  write ( *, '(a)' ) '  RANDOM_SEED is a FORTRAN90 function which allows'
  write ( *, '(a)' ) '  the user to read or set the random number seed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CALL RANDOM_SEED ( SIZE = seed_size ) returns the'
  write ( *, '(a)' ) '  dimension of the random number seed.'

  call random_seed ( size = seed_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  SIZE = ', seed_size

  allocate ( seed(1:seed_size ) )

  call random_seed ( get = seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call RANDOM_SEED ( GET = SEED )'
  write ( *, '(a)' ) '  to get current contents of the SEED array:'
  write ( *, '(a)' ) ' '

  do i = 1, seed_size
    write ( *, '(i12)' ) seed(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Set the SEED array to a simple value.'
  write ( *, '(a)' ) ' '
 
  do i = 1, seed_size
    seed(i) = i
    write ( *, '(i12)' ) seed(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call RANDOM_SEED ( PUT = SEED )'
  write ( *, '(a)' ) '  to reset current contents of the SEED array:'
  write ( *, '(a)' ) ' '

  call random_seed ( put = seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call RANDOM_NUMBER 5 times.'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call random_number ( harvest = r )
    write ( *, '(2x,g14.6)' ) r
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call RANDOM_SEED ( GET = SEED )'
  write ( *, '(a)' ) '  to get current contents of the SEED array.'
  write ( *, '(a)' ) '  Notice that the seed has changed!'
  write ( *, '(a)' ) ' '

  call random_seed ( get = seed )

  do i = 1, seed_size
    write ( *, '(i12)' ) seed(i)
  end do

  do i = 1, seed_size
    seed(i) = i
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call RANDOM_SEED ( PUT = SEED )'
  write ( *, '(a)' ) '  to reset the SEED array:'
  write ( *, '(a)' ) ' '

  call random_seed ( put = seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call RANDOM_NUMBER 5 times.'
  write ( *, '(a)' ) '  Because the SEED array is the same,'
  write ( *, '(a)' ) '  the values should be repeated.'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call random_number ( harvest = r )
    write ( *, '(2x,g14.6)' ) r
  end do

  deallocate ( seed )

  return
end
subroutine test_range ( )

!*****************************************************************************80
!
!! TEST_RANGE tests RANGE.
!
!  Discussion:
!
!    The FORTRAN90 function RANGE returns an integer that is the
!    logarithm base 10 of the magnitude of the largest numbers of
!    a given kind and type.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) c4
  complex ( kind = 8 ) c8
  integer ( kind = 4 ) i4
  integer ( kind = 8 ) i8
  real ( kind = 4 ) r4
  real ( kind = 8 ) r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_RANGE'
  write ( *, '(a)' ) '  RANGE is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  base 10 logarithm of the largest magnitude objects'
  write ( *, '(a)' ) '  of a given type and kind.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2x,i8)' ) '  Integer ( kind = 4 )', range ( i4 )
  write ( *, '(a,2x,i8)' ) '  Integer ( kind = 8 )', range ( i8 )
  write ( *, '(a,2x,i8)' ) '  Real ( kind = 4 )   ', range ( r4 )
  write ( *, '(a,2x,i8)' ) '  Real ( kind = 8 )   ', range ( r8 )
  write ( *, '(a,2x,i8)' ) '  Complex ( kind = 4 )', range ( c4 )
  write ( *, '(a,2x,i8)' ) '  Complex ( kind = 8 )', range ( c8 )

  return
end
subroutine test_real_c4 ( )

!*****************************************************************************80
!
!! TEST_REAL_C4 tests REAL as applied to complex numbers.
!
!  Discussion:
!
!    The FORTRAN90 function REAL can return the real part of a complex number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex c4_uniform_01
  complex c
  integer i
  real r
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_REAL_C4'
  write ( *, '(a)' ) '  REAL is a FORTRAN90 function which can return the'
  write ( *, '(a)' ) '  real part of a complex number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                  X                         REAL(X)'
  write ( *, '(a)' ) '       ------------------------    ----------------'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    c = c4_uniform_01 ( seed )
    r = real ( c )
    write ( *, '(2x,f14.6,f14.6,6x,f14.6)' ) c, r
  end do

  return
end
subroutine test_repeat ( )

!*****************************************************************************80
!
!! TEST_REPEAT tests REPEAT.
!
!  Discussion:
!
!    The FORTRAN90 function REPEAT makes a new character string by
!    repeating a given one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer n
  character ( len = 1 ) s1
  character ( len = 2 ) s2
  character ( len = 4 ) s4
  character ( len = 5 ) s5
  character ( len = 8 ) s8
  character ( len = 12 ) s12

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_REPEAT'
  write ( *, '(a)' ) '  REPEAT(S,N) creates a new character string by repeating'
  write ( *, '(a)' ) '  a given string S N times.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S                    N  REPEAT(S,N)'
  write ( *, '(a)' ) ' '

  s1 = 'a'
  n = 5
  s5 = repeat ( s1, n )
  write ( *, '(2x,a3,9x,2x,i8,2x,a)' ) '"'//s1//'"', n, '"'//s5//'"'

  s2 = 'Ab'
  n = 4
  s8 = repeat ( s2, n )
  write ( *, '(2x,a4,8x,2x,i8,2x,a)' ) '"'//s2//'"', n, '"'//s8//'"'

  s4 = 'Abc '
  n = 3
  s12 = repeat ( s4, n )
  write ( *, '(2x,a6,6x,2x,i8,2x,a)' ) '"'//s4//'"', n, '"'//s12//'"'

  return
end
subroutine test_reshape ( )

!*****************************************************************************80
!
!! TEST_RESHAPE tests RESHAPE.
!
!  Discussion:
!
!    The FORTRAN90 function RESHAPE can reshape an array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer a(3,4)
  integer b(12)
  integer c(2,6)
  integer d(2,3,2)
  integer i
  integer j

  do i = 1, 3
    do j = 1, 4
      a(i,j) = 10 * i + j
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_RESHAPE'
  write ( *, '(a)' ) '  RESHAPE can "reshape" an array.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A is a 3 x 4 matrix, and we reshape it to a vector:'
  write ( *, '(a)' ) ' '
  do i = 1, 3
    write ( *, '(4i6)' ) a(i,1:4)
  end do

  b = reshape ( a, (/ 12 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We reshape A to a vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  B = RESHAPE ( A, (/ 12 /) )'
  write ( *, '(a)' ) ' '
  write ( *, '(12i6)' ) b(1:12)

  c = reshape ( b, (/ 2, 6 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We reshape B to a matrix C:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C = RESHAPE ( B, (/ 2, 6 /) )'
  write ( *, '(a)' ) ' '

  do i = 1, 2
    write ( *, '(6i6)' ) c(i,1:6)
  end do

  d = reshape ( b, (/ 2, 3, 2 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We reshape B to a 3D matrix D:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  D = RESHAPE ( B, (/ 2, 3, 2 /) )'
  write ( *, '(a)' ) ' '

  do i = 1, 2
    do j = 1, 3
      write ( *, '(6i6)' ) d(i,j,1:2)
    end do
    write ( *, '(a)' ) ' '
  end do

  return
end
subroutine test_rrspacing ( )

!*****************************************************************************80
!
!! TEST_RRSPACING tests RRSPACING.
!
!  Discussion:
!
!    The FORTRAN90 function RRSPACING returns the relative spacing associated
!    with a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) x_r4
  real ( kind = 8 ) x_r8
  real ( kind = 4 ) y_r4
  real ( kind = 8 ) y_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_RRSPACING'
  write ( *, '(a)' ) '  RRSPACING is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  reciprocal relative spacing associated '
  write ( *, '(a)' ) '  with a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First, some "default precision" reals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              RRSPACING(X)'
  write ( *, '(a)' ) ' '
  x_r4 = 1.0
  y_r4 = rrspacing ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  x_r4 = 0.0
  y_r4 = rrspacing ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  x_r4 = 1000000.0
  y_r4 = rrspacing ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now, some "double precision" reals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              RRSPACING(X)'
  write ( *, '(a)' ) ' '
  x_r8 = 1.0D+00
  y_r8 = rrspacing ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8
  x_r8 = 0.0D+00
  y_r8 = rrspacing ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8
  x_r8 = 1000000.0D+00
  y_r8 = rrspacing ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8

  return
end
subroutine test_scan ( )

!*****************************************************************************80
!
!! TEST_SCAN tests SCAN.
!
!  Discussion:
!
!    The FORTRAN90 function SCAN determines the first occurrence
!    of any character in a set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SCAN'
  write ( *, '(a)' ) '  SCAN(S,SET) is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the location of the first occurrence of '
  write ( *, '(a)' ) '  any character from SET in string S.'
  write ( *, '(a)' ) '  SCAN(S,SET,.TRUE.) returns the location of the'
  write ( *, '(a)' ) '  LAST occurrence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  scan ( ''THE CATATONIC CAT'', ''CAT'' )', &
    scan ( 'THE CATATONIC CAT', 'CAT' )
  write ( *, '(a,i8)' ) '  scan ( ''THE CATATONIC CAT'', ''DOG'' )', &
    scan ( 'THE CATATONIC CAT', 'DOG' )
  write ( *, '(a,i8)' ) '  scan ( ''THE CATATONIC CAT'', ''cat'' )', &
    scan ( 'THE CATATONIC CAT', 'cat' )
  write ( *, '(a,i8)' ) '  scan ( ''THE CATATONIC CAT'', ''ABC'' )', &
    scan ( 'THE CATATONIC CAT', 'ABC' )
  write ( *, '(a,i8)' ) '  scan ( ''THE CATATONIC CAT'', ''ABC'', .TRUE. )', &
    scan ( 'THE CATATONIC CAT', 'ABC', .TRUE. )

  return
end
subroutine test_set_exponent ( )

!*****************************************************************************80
!
!! TEST_SET_EXPONENT tests SET_EXPONENT.
!
!  Discussion:
!
!    The FORTRAN90 function SET_EXPONENT(X,I) returns 
!    FRACTION(X)*RADIX**I where FRACTION is the fractional
!    part of X, and RADIX is the base in the model representation of X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SET_EXPONENT'
  write ( *, '(a)' ) '  SET_EXPONENT(X,I) is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  FRACTION(X)*RADIX**I, where FRACTION(X) is the '
  write ( *, '(a)' ) '  fractional part of X, and RADIX is the base'
  write ( *, '(a)' ) '  for real number arithmetic.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X             I      SET_EXPONENT(X,I)'
  write ( *, '(a)' ) ' '
  x = 100.0
  do i = -5, 5
    y = set_exponent ( x, i )
    write ( *, '(2x,g14.6,2x,i4,2x,g14.6)' ) x, i, y
  end do

  return
end
subroutine test_scale ( )

!*****************************************************************************80
!
!! TEST_SCALE tests SCALE.
!
!  Discussion:
!
!    The FORTRAN90 function SCALE(X,I) returns X*RADIX**I where RADIX
!    is the base in the model representation of X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SCALE'
  write ( *, '(a)' ) '  SCALE(X,I) is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  value of X * RADIX**I, where RADIX is the base'
  write ( *, '(a)' ) '  for real number arithmetic.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X             I      SCALE(X,I)'
  write ( *, '(a)' ) ' '
  x = 5.0
  do i = -5, 5
    y = scale ( x, i )
    write ( *, '(2x,g14.6,2x,i4,2x,g14.6)' ) x, i, y
  end do

  return
end
subroutine test_selected_int_kind ( )

!*****************************************************************************80
!
!! TEST_SELECTED_INT_KIND tests SELECTED_INT_KIND.
!
!  Discussion:
!
!    The FORTRAN90 function SELECTED_INT_KIND(R) returns a value of the
!    KIND parameter for integers up to 10**R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer k
  integer r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SELECTED_INT_KIND'
  write ( *, '(a)' ) '  SELECTED_INT_KIND(R) is a FORTRAN90 function which'
  write ( *, '(a)' ) '  returns a value of the KIND parameter for integers'
  write ( *, '(a)' ) '  up to 10**R.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         R    SELECTED_INT_KIND(R)'
  write ( *, '(a)' ) ' '

  do r = 1, 20
    k = selected_int_kind ( r )
    write ( *, '(2x,i8,2x,i8)' ) r, k
  end do

  return
end
subroutine test_selected_real_kind ( )

!*****************************************************************************80
!
!! TEST_SELECTED_REAL_KIND tests SELECTED_REAL_KIND.
!
!  Discussion:
!
!    The FORTRAN90 function SELECTED_REAL_KIND(P,R) returns a value of the
!    KIND parameter for reals with P digits of precision and a
!    decimal exponent of up to 10**R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer k
  integer p
  integer r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SELECTED_REAL_KIND'
  write ( *, '(a)' ) '  SELECTED_INT_KIND(R) is a FORTRAN90 function which'
  write ( *, '(a)' ) '  returns a value of the KIND parameter for integers'
  write ( *, '(a)' ) '  up to 10**R.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         P         R    SELECTED_REAL_KIND(P,R)'
  write ( *, '(a)' ) ' '

  p = 5

  do r = 1, 20
    k = selected_real_kind ( p, r )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) p, r, k
  end do

  write ( *, '(a)' ) ' '

  r = 10

  do p = 1, 20
    k = selected_real_kind ( p, r )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) p, r, k
  end do


  return
end
subroutine test_shape ( )

!*****************************************************************************80
!
!! TEST_SHAPE tests SHAPE.
!
!  Discussion:
!
!    The FORTRAN90 function SHAPE(ARRAY) returns the "shape"
!    of the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real a(5,10,17)
  integer b(4:6,-5:-1,10:20)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SHAPE'
  write ( *, '(a)' ) '  SHAPE(ARRAY) is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  "shape" of the array.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  real a(5,10,17)'
  write ( *, '(a,3i6)' ) '  shape(a)   = ', shape(a)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer b(4:6,-5:-1,10:20)'
  write ( *, '(a,3i6)' ) '  shape(b)   = ', shape(b)

  return
end
subroutine test_sign ( )

!*****************************************************************************80
!
!! TEST_SIGN tests SIGN.
!
!  Discussion:
!
!    The FORTRAN90 function SIGN(X,Y) transfers the sign of Y to the
!    magnitude of X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SIGN'
  write ( *, '(a)' ) '  SIGN is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  sign of Y times the magnitude of X.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X               Y           SIGN(X,Y)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = r8_uniform ( x_lo, x_hi, seed )
    z = sign ( x, y )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, y, z
  end do

  return
end
subroutine test_sin_r8 ( )

!*****************************************************************************80
!
!! TEST_SIN_R8 tests SIN on real ( kind = 8 ) arguments.
!
!  Discussion:
!
!    The FORTRAN90 function SIN returns the sine of a real or complex number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SIN_R8'
  write ( *, '(a)' ) '  SIN is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  sine of a real or complex number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we use real ( kind = 8 ) arguments.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              SIN(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = sin ( x )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, y
  end do

  return
end
subroutine test_sinh ( )

!*****************************************************************************80
!
!! TEST_SINH tests SINH.
!
!  Discussion:
!
!    The FORTRAN90 function SINH returns the hyperbolic sine of a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SINH'
  write ( *, '(a)' ) '  SINH is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  hyperbolic sine of a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              SINH(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = sinh ( x )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, y
  end do

  return
end
subroutine test_size ( )

!*****************************************************************************80
!
!! TEST_SIZE tests SIZE.
!
!  Discussion:
!
!    The FORTRAN90 function SIZE(ARRAY) returns the "size"
!    of all the dimensions; SIZE(ARRAY,DIM) returns the size
!    in the given dimension.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real a(5,10,17)
  integer b(4:6,-5:-1,10:20)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SIZE'
  write ( *, '(a)' ) '  SIZE(ARRAY) is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  size in all dimensions;'
  write ( *, '(a)' ) '  SIZE(ARRAY,DIM) returns the size in dimension'
  write ( *, '(a)' ) '  DIM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  real a(5,10,17)'
  write ( *, '(a,i6)' ) '  size(a)   = ', size(a)
  write ( *, '(a,i6)' ) '  size(a,1) = ', size(a,1)
  write ( *, '(a,i6)' ) '  size(a,2) = ', size(a,2)
  write ( *, '(a,i6)' ) '  size(a,3) = ', size(a,3)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer b(4:6,-5:-1,10:20)'
  write ( *, '(a,i6)' ) '  size(b)   = ', size(b)
  write ( *, '(a,i6)' ) '  size(b,1) = ', size(b,1)
  write ( *, '(a,i6)' ) '  size(b,2) = ', size(b,2)
  write ( *, '(a,i6)' ) '  size(b,3) = ', size(b,3)

  return
end
subroutine test_spacing ( )

!*****************************************************************************80
!
!! TEST_SPACING tests SPACING.
!
!  Discussion:
!
!    The FORTRAN90 function SPACING returns the absolute spacing associated
!    with a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) x_r4
  real ( kind = 8 ) x_r8
  real ( kind = 4 ) y_r4
  real ( kind = 8 ) y_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SPACING'
  write ( *, '(a)' ) '  SPACING is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  absolute spacing associated with a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First, some "default precision" reals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              SPACING(X)'
  write ( *, '(a)' ) ' '
  x_r4 = 1.0
  y_r4 = spacing ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  x_r4 = 0.0
  y_r4 = spacing ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  x_r4 = 1000000.0
  y_r4 = spacing ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now, some "double precision" reals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              SPACING(X)'
  write ( *, '(a)' ) ' '
  x_r8 = 1.0D+00
  y_r8 = spacing ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8
  x_r8 = 0.0D+00
  y_r8 = spacing ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8
  x_r8 = 1000000.0D+00
  y_r8 = spacing ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8

  return
end
subroutine test_spread ( )

!*****************************************************************************80
!
!! TEST_SPREAD tests SPREAD.
!
!  Discussion:
!
!    The FORTRAN90 function SPREAD replicates an array by adding a dimension.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a1(4,3)
  integer ( kind = 4 ) a2(3,4)
  integer i
  integer ( kind = 4 ) s
  integer ( kind = 4 ) v(4)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SPREAD'
  write ( *, '(a)' ) '  SPREAD is a FORTRAN90 function which replicates'
  write ( *, '(a)' ) '  an array by adding a dimension.'
  write ( *, '(a)' ) ' '

  s = 99

  write ( *, '(a,i6)' ) '  Suppose we have a scalar S = ', s
  write ( *, '(a)' ) ' '

  v = spread ( s, 1, 4 )

  write ( *, '(a)' ) '  V = spread ( s, 1, 4 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  adds a new dimension (1) of extent 4'
  write ( *, '(a)' ) ' '
  write ( *, '(4i6)' ) v(1:4)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now first reset V to (1,2,3,4)'
  v = (/ 1, 2, 3, 4 /)

  a1 = spread ( v, 2, 3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A1 = spread ( v, 2, 3 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  adds a new dimension (2) of extent 3'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(3i6)' ) a1(i,1:3)
  end do

  a2 = spread ( v, 1, 3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A2 = spread ( v, 1, 3 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  adds a new dimension (1) of extent 3'
  write ( *, '(a)' ) ' '
  do i = 1, 3
    write ( *, '(4i6)' ) a2(i,1:4)
  end do

  return
end
subroutine test_sqrt ( )

!*****************************************************************************80
!
!! TEST_SQRT tests SQRT.
!
!  Discussion:
!
!    The FORTRAN90 function SQRT returns the square root of a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo =  0.0D+00
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SQRT'
  write ( *, '(a)' ) '  SQRT is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  square root of a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              SQRT(X)        (SQRT(X))**2'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = sqrt ( x )
    z = y * y
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, y, z
  end do

  return
end
subroutine test_sum ( )

!*****************************************************************************80
!
!! TEST_SUM tests SUM.
!
!  Discussion:
!
!    The FORTRAN90 function SUM returns the sum of
!    the entries of a complex, integer, or real vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(4,5)
  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer :: seed = 123456789
  complex ( kind = 4 ) x_c4(n)
  integer ( kind = 4 ) x_i4(n)
  real    ( kind = 4 ) x_r4(n)

  do i = 1, n
    x_c4(i) = c4_uniform_01 ( seed )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SUM'
  write ( *, '(a)' ) '  SUM is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  sum of the entries of a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A complex ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                 X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,2f14.6)' ) x_c4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,2f14.6)' ) '  SUM(X) = ', sum ( x_c4 )

  do i = 1, n
    x_r4(i) = real ( x_c4(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A real ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,f14.6)' ) x_r4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.6)' ) '  SUM(X) = ', sum ( x_r4 )

  do i = 1, n
    x_i4(i) = int ( 10.0 * x_r4(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  An integer ( kind = 4 ) vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               X'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i14)' ) x_i4(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,i14)' ) '  SUM(X) = ', sum ( x_i4 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A 4 x 5 integer array A containing 1 through 20'
  write ( *, '(a)' ) ' '
  k = 0
  do i = 1, 4
    do j = 1, 5
      k = k + 1
      a(i,j) = k
    end do
  end do

  do i = 1, 4
    write ( *, '(2x,5i4)' ) a(i,1:5)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i14)' ) '  SUM(A) = ', sum ( a )

  return
end
subroutine test_sum_dim ( )

!*****************************************************************************80
!
!! TEST_SUM_DIM tests SUM specifying a dimension.
!
!  Discussion:
!
!    The FORTRAN90 function SUM ( A, DIM=I ) returns the sum of the entries of a 
!    complex, integer, or real array along dimension I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(4,5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) sa
  integer ( kind = 4 ) sa1(5)
  integer ( kind = 4 ) sa2(4)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SUM_ARAY'
  write ( *, '(a)' ) '  SUM(A) is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  sum of the entries of a vector or array.'
  write ( *, '(a)' ) '  SUM(A,1) or SUM(A,DIM=1) returns an array of one lower '
  write ( *, '(a)' ) '  rank containing sums over index 1, and so on.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A 4 x 5 integer array A containing 1 through 20'
  write ( *, '(a)' ) ' '

  k = 0
  do i = 1, 4
    do j = 1, 5
      k = k + 1
      a(i,j) = k
    end do
  end do

  do i = 1, 4
    write ( *, '(2x,5i4)' ) a(i,1:5)
  end do

  sa = sum ( a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i14)' ) '  SUM(A) = ', sa

  sa1 = sum ( a, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SUM(A,1) = SUM(A,DIM=1) will give sums across rows:'
  write ( *, '(a)' ) ' '
  do i = 1, 5 
    write ( *, '(i4)' ) sa1(i)
  end do
!
!  DIM= can be specified.
!
  sa2 = sum ( a, dim = 2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SUM(A,2) = SUM(A,DIM=2) will give sums across columns:'
  write ( *, '(a)' ) ' '
  do i = 1, 4 
    write ( *, '(i4)' ) sa2(i)
  end do

  return
end
subroutine test_system_clock ( )

!*****************************************************************************80
!
!! TEST_SYSTEM_CLOCK tests SYSTEM_CLOCK.
!
!  Discussion:
!
!    The FORTRAN90 subroutine SYSTEM_CLOCK returns information from a
!    system clock.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer count
  integer count_max
  integer count_rate
  real t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SYSTEM_CLOCK'
  write ( *, '(a)' ) '  SYSTEM_CLOCK is a FORTRAN90 subroutine which returns'
  write ( *, '(a)' ) '  information from a system clock.'

  call system_clock ( count, count_rate, count_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  COUNT =      ', count
  write ( *, '(a,i12)' ) '  COUNT_RATE = ', count_rate
  write ( *, '(a,i12)' ) '  COUNT_MAX =  ', count_max

  t = real ( count ) / real ( count_rate )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The current value of COUNT corresponds to a '
  write ( *, '(a,g14.6,a)' ) &
    '  a time interval of COUNT / COUNT_RATE ', t, ' seconds,'

  return
end
subroutine test_tan ( )

!*****************************************************************************80
!
!! TEST_TAN tests TAN.
!
!  Discussion:
!
!    The FORTRAN90 function TAN returns the tangent of a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TAN'
  write ( *, '(a)' ) '  TAN is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  tangent of a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              TAN(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = tan ( x )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, y
  end do

  return
end
subroutine test_tanh ( )

!*****************************************************************************80
!
!! TEST_TANH tests TANH.
!
!  Discussion:
!
!    The FORTRAN90 function TANH returns the hyperbolic tangent of a 
!    real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 ) i
  integer ( kind = 8 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi = 10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TANH'
  write ( *, '(a)' ) '  TANH is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  hyperbolic tangent of a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              TANH(X)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    x = r8_uniform ( x_lo, x_hi, seed )
    y = tanh ( x )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, y
  end do

  return
end
subroutine test_tiny ( )

!*****************************************************************************80
!
!! TEST_TINY tests TINY.
!
!  Discussion:
!
!    The FORTRAN90 function TINY returns a "tiny number" associated
!    with a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) x_r4
  real ( kind = 8 ) x_r8
  real ( kind = 4 ) y_r4
  real ( kind = 8 ) y_r8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TINY'
  write ( *, '(a)' ) '  TINY is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  "tiniest number" associated with a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First, some "default precision" reals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              TINY(X)'
  write ( *, '(a)' ) ' '
  x_r4 = 1.0
  y_r4 = tiny ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  x_r4 = 0.0
  y_r4 = tiny ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  x_r4 = 1000000.0
  y_r4 = tiny ( x_r4 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r4, y_r4
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now, some "double precision" reals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              TINY(X)'
  write ( *, '(a)' ) ' '
  x_r8 = 1.0D+00
  y_r8 = tiny ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8
  x_r8 = 0.0D+00
  y_r8 = tiny ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8
  x_r8 = 1000000.0D+00
  y_r8 = tiny ( x_r8 )
  write ( *, '(2x,g14.6,2x,g14.6)' ) x_r8, y_r8

  return
end
subroutine test_transfer ( )

!*****************************************************************************80
!
!! TEST_TRANSFER tests TRANSFER.
!
!  Discussion:
!
!    The FORTRAN90 function TRANSFER allows the data stored as a variable
!    of one type to be interpreted as data defining a variable of another 
!    type.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) c
  character ( len = 1 ) ch
  integer ( kind = 4 ) i
  logical l
  real ( kind = 4 ) r
  real ( kind = 4 ) r2(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TRANSFER'
  write ( *, '(a)' ) '  TRANSFER is a FORTRAN90 function which allows'
  write ( *, '(a)' ) '  the "physical" data of a variable to be interpreted'
  write ( *, '(a)' ) '  as though it were another type or kind of variable.'

  r = 1.0
  i = transfer ( r, i )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine REAL data as though it were an INTEGER.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I = TRANSFER ( R, I )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R = ', r
  write ( *, '(a,i12)'   ) '  I = ', i

  c = ( 1.2, 3.4 )
  r2 = transfer ( c, r2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine COMPLEX data as though it were two REALS.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R2 = TRANSFER ( C, R2 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  C =  ', c
  write ( *, '(a,2g14.6)' ) '  R2 = ', r2(1:2)

  l = .true.
  i = transfer ( l, i )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine LOGICAL data as though it were an INTEGER.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I = TRANSFER ( L, I )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,l1)' ) '  L = ', l
  write ( *, '(a,i12)' ) '  I = ', i

  ch = 'a'
  i = transfer ( c, i )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine CHARACTER data as though it were an INTEGER.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I = TRANSFER ( CH, I )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '  CH = ', '"' // ch // '"'
  write ( *, '(a,i12)' ) '  I = ', i

  return
end
subroutine test_transpose ( )

!*****************************************************************************80
!
!! TEST_TRANSPOSE tests TRANSPOSE.
!
!  Discussion:
!
!    The FORTRAN90 function TRANSPOSE returns the transpose of a 
!    complex, integer or real two dimensional array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x_i4(4,4)
  real ( kind = 4 ) x_r4(2,3)
  integer ( kind = 4 ) y_i4(4,4)
  real ( kind = 4 ) y_r4(3,2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TRANSPOSE'
  write ( *, '(a)' ) '  TRANSPOSE is a FORTRAN90 function which returns'
  write ( *, '(a)' ) '  the transpose of a complex, integer or real '
  write ( *, '(a)' ) '  two dimensional array.'

  do i = 1, 4
    do j = 1, 4
      x_i4(i,j) = i4_uniform ( 0, 100, seed )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer ( kind = 4 ) x_i4(4,4)'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(4(2x,i4))' ) x_i4(i,1:4)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  y_i4 = transpose ( x_i4 )'
  write ( *, '(a)' ) ' '

  y_i4 = transpose ( x_i4 )

  do i = 1, 4
    write ( *, '(4(2x,i4))' ) y_i4(i,1:4)
  end do

  call random_number ( harvest = x_r4(1:2,1:3) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  real ( kind = 4 ) x_r4(2,3)'
  write ( *, '(a)' ) ' '
  do i = 1, 2
    write ( *, '(4(2x,f14.6))' ) x_r4(i,1:3)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  y_r4 = transpose ( x_r4 )'
  write ( *, '(a)' ) ' '

  y_r4 = transpose ( x_r4 )

  do i = 1, 3
    write ( *, '(4(2x,f14.6))' ) y_r4(i,1:2)
  end do

  return
end
subroutine test_trim ( )

!*****************************************************************************80
!
!! TEST_TRIM tests TRIM.
!
!  Discussion:
!
!    The FORTRAN90 function TRIM returns a copy of the string from which
!    the trailing blanks have been omitted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TRIM'
  write ( *, '(a)' ) '  TRIM is a FORTRAN90 function which returns a copy'
  write ( *, '(a)' ) '  of a string from which the trailing blanks have been'
  write ( *, '(a)' ) '  dropped.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      S          TRIM(S)'
  write ( *, '(a)' ) '   ----------    ----------'
  write ( *, '(a)' ) ' '
  s = '1234567890'
  write ( *, '(a)' ) '  "' // s // '"  "' // trim ( s ) // '"'
  s = '12345     '
  write ( *, '(a)' ) '  "' // s // '"  "' // trim ( s ) // '"'
  s = '     67890'
  write ( *, '(a)' ) '  "' // s // '"  "' // trim ( s ) // '"'
  s = '  34 678  '
  write ( *, '(a)' ) '  "' // s // '"  "' // trim ( s ) // '"'
  s = '    5     '
  write ( *, '(a)' ) '  "' // s // '"  "' // trim ( s ) // '"'
  s = '          '
  write ( *, '(a)' ) '  "' // s // '"  "' // trim ( s ) // '"'

  return
end
subroutine test_ubound ( )

!*****************************************************************************80
!
!! TEST_UBOUND tests UBOUND.
!
!  Discussion:
!
!    The FORTRAN90 function UBOUND(ARRAY) returns the upper bounds
!    of all the dimensions; UBOUND(ARRAY,DIM) returns the upper bound
!    in the given dimension.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real a(5,10,17)
  integer b(4:6,-5:-1,10:20)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_UBOUND'
  write ( *, '(a)' ) '  UBOUND(ARRAY) is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  upper array bounds in all dimensions;'
  write ( *, '(a)' ) '  UBOUND(ARRAY,DIM) returns the upper bound in dimension'
  write ( *, '(a)' ) '  DIM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  real a(5,10,17)'
  write ( *, '(a,3i6)' ) '  ubound(a)   = ', ubound(a)
  write ( *, '(a,3i6)' ) '  ubound(a,1) = ', ubound(a,1)
  write ( *, '(a,3i6)' ) '  ubound(a,2) = ', ubound(a,2)
  write ( *, '(a,3i6)' ) '  ubound(a,3) = ', ubound(a,3)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer b(4:6,-5:-1,10:20)'
  write ( *, '(a,3i6)' ) '  ubound(b)   = ', ubound(b)
  write ( *, '(a,3i6)' ) '  ubound(b,1) = ', ubound(b,1)
  write ( *, '(a,3i6)' ) '  ubound(b,2) = ', ubound(b,2)
  write ( *, '(a,3i6)' ) '  ubound(b,3) = ', ubound(b,3)

  return
end
subroutine test_unpack ( )

!*****************************************************************************80
!
!! TEST_UNPACK tests UNPACK.
!
!  Discussion:
!
!    The FORTRAN90 function UNPACK unpacks a vector into an array 
!    subject to a logical mask.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer a(5,4)
  integer field
  integer i
  integer i4_uniform
  integer j
  logical mask(5,4)
  integer n
  integer :: seed = 123456789
  integer v(20)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_UNPACK'
  write ( *, '(a)' ) '  UNPACK is a FORTRAN90 function which unpacks the'
  write ( *, '(a)' ) '  entries of a vector into an array subject to'
  write ( *, '(a)' ) '  a logical mask.'
  write ( *, '(a)' ) ' '

  do i = 1, 20
    v(i) = i
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i2,a)' ) '  integer v(20)'
  write ( *, '(a)' ) ' '

  write ( *, '(20(2x,i4))' ) v(1:20)

  do i = 1, 5
    do j = 1, 4
      mask(i,j) = ( mod ( i + j, 2 ) /= 0 )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  logical mask(5,4)'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    write ( *, '(4(2x,l1))' ) mask(i,1:4)
  end do

  field = -99

  a = unpack ( v, mask, field )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  integer a(5,4) = unpack ( v, mask, field )'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    write ( *, '(4(2x,i4))' ) a(i,1:4)
  end do

  return
end
subroutine test_verify ( )

!*****************************************************************************80
!
!! TEST_VERIFY tests VERIFY.
!
!  Discussion:
!
!    The FORTRAN90 function VERIFY(S1,S2,BACK) returns the location of the 
!    first character in S1 that is not in S2, or, if BACK is TRUE, the
!    last character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  logical back
  character ( len = 10 ) s1
  character ( len = 2 ) s2
  character ( len = 10 ) s3
  character ( len = 14 ) s4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_VERIFY'
  write ( *, '(a)' ) '  VERIFY(S1,S2) is a FORTRAN90 function which returns the'
  write ( *, '(a)' ) '  location of the first character in S1 not in S2,'
  write ( *, '(a)' ) '  or else 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            S1          S2  VERIFY(S1,S2)'
  write ( *, '(a)' ) ' '
  s1 = '1001010001'
  s2 = '01'
  write ( *, '(2x,a14,2x,a10,2x,i2)' ) s1, s2, verify(s1,s2)
  s1 = '1002010091'
  s2 = '01'
  write ( *, '(2x,a14,2x,a10,2x,i2)' ) s1, s2, verify(s1,s2)
  s1 = 'CAPS lower'
  s3 = 'CAPS LOWER'
  write ( *, '(2x,a14,2x,a10,2x,i2)' ) s1, s3, verify(s1,s3)
  s4 = 'Blanks count!'
  s3 = 'aBcklnstu!'
  write ( *, '(2x,a14,2x,a10,2x,i2)' ) s4, s3, verify(s4,s3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The optional parameter BACK, set to TRUE, makes'
  write ( *, '(a)' ) '  the check start at the END of S1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            S1          S2  VERIFY(S1,S2,BACK)'
  write ( *, '(a)' ) ' '
  back = .true.

  s1 = '1001010001'
  s2 = '01'
  write ( *, '(2x,a14,2x,a10,2x,i2)' ) s1, s2, verify(s1,s2,back)
  s1 = '1002010091'
  s2 = '01'
  write ( *, '(2x,a14,2x,a10,2x,i2)' ) s1, s2, verify(s1,s2,back)
  s1 = 'CAPS lower'
  s3 = 'CAPS LOWER'
  write ( *, '(2x,a14,2x,a10,2x,i2)' ) s1, s3, verify(s1,s3,back)
  s4 = 'Blanks count!'
  s3 = 'aBcklnstu!'
  write ( *, '(2x,a14,2x,a10,2x,i2)' ) s4, s3, verify(s4,s3,back)

  return
end
function c4_uniform_01 ( seed )

!*****************************************************************************80
!
!! C4_UNIFORM_01 returns a unit complex pseudorandom number.
!
!  Discussion:
!
!    The angle should be uniformly distributed between 0 and 2 * PI,
!    the square root of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, complex C4_UNIFORM_01, a pseudorandom complex value.
!
  implicit none

  complex c4_uniform_01
  integer k
  real, parameter :: pi = 3.1415926E+00
  real r
  integer seed
  real theta

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = sqrt (  real ( seed, kind = 4 ) * 4.656612875E-10 )

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  theta = 2.0E+00 * pi * ( real ( seed, kind = 4 ) * 4.6566128754E-10 )

  c4_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ) )

  return
end
function ch_is_printable ( ch )

!*****************************************************************************80
!
!! CH_IS_PRINTABLE is TRUE if C is printable.
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
!    31 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character CH, a character to check.
!
!    Output, logical CH_IS_PRINTABLE is TRUE if C is a printable character.
!
  implicit none

  character ch
  logical ch_is_printable
  integer i

  i = iachar ( ch )

  if ( 32 <= i .and. i <= 126 ) then
    ch_is_printable = .true.
  else
    ch_is_printable = .false.
  end if

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a pseudorandom integer of KIND = 4.
!
!  Discussion:
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  real ( kind = 4 ) r
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed

  r = r4_uniform_01 ( seed )

  r = ( 1.0E+00 - r ) * real ( a, kind = 4 ) &
    +             r   * real ( b, kind = 4 )

  i4_uniform = nint ( r, kind = 4 )

  return
end
subroutine i4vec_uniform ( n, a, b, seed, x )

!*****************************************************************************80
!
!! I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The pseudorandom numbers should be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) X(N), a vector of numbers between A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
    r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
      +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
    value = nint ( r, kind = 4 )

    value = max ( value, min ( a, b ) )
    value = min ( value, max ( a, b ) )

    x(i) = value

  end do

  return
end
function i8_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I8_UNIFORM returns a pseudorandom integer of KIND = 8.
!
!  Discussion:
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 8 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 8 ) I8_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 8 ) a
  integer ( kind = 8 ) b
  integer ( kind = 8 ) i8_uniform
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 8 ) seed

  r = r8_uniform_01 ( seed )

  r = ( 1.0D+00 - r ) * real ( a, kind = 8 ) &
    +             r   * real ( b, kind = 8 )

  i8_uniform = nint ( r, kind = 8 )

  return
end
function r4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R4_UNIFORM returns a scaled real ( kind = 4 ) pseudorandom number.
!
!  Discussion:
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ) k
  real ( kind = 4 ) r4_uniform
  integer ( kind = 4 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r4_uniform = a + ( b - a ) * real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit real ( kind = 4 ) pseudorandom number.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r4_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r4_uniform_01

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
function r8_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R8_UNIFORM returns a scaled real  ( kind = 8 ) pseudorandom number.
!
!  Discussion:
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 8 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 8 ) k
  real ( kind = 8 ) r8_uniform
  integer ( kind = 8 )seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r8_uniform = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit real  ( kind = 8 ) pseudorandom number.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 8 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 8 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 8 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

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
