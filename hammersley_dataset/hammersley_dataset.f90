program main

!*****************************************************************************80
!
!! MAIN is the main program for HAMMERSLEY_DATASET.
!
!  Discussion:
!
!    HAMMERSLEY_DATASET generates a Hammersley dataset and writes it to a file.
!
!    This program is meant to be used interactively.  It's also
!    possible to prepare a simple input file beforehand and use it
!    in batch mode.
!
!    The program requests input values from the user:
!
!    * DIM_NUM, the spatial dimension,
!    * N, the number of points to generate,
!    * STEP, the index of the first subsequence element to be computed.
!    * SEED(1:DIM_NUM), the sequence index corresponding to STEP = 0.
!    * LEAP(1:DIM_NUM), the successive jumps in the sequence.
!    * BASE(1:DIM_NUM), the bases (usually distinct primes or -N).
!
!    The program generates the data, writes it to the file
!
!      hammersley_DIM_NUM_N.txt
!
!    where "DIM_NUM" and "N" are the numeric values specified by the user,
!    and then asks the user for more input.   To indicate that no further
!    computations are desired, it is enough to input a nonsensical
!    value, such as -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: base
  character ( len = 80 ) :: file_out_name
  logical halham_leap_check
  logical halham_n_check
  logical halham_dim_num_check
  logical halham_seed_check
  logical halham_step_check
  logical hammersley_base_check
  integer ( kind = 4 ) ios
  integer ( kind = 4 ), allocatable, dimension ( : ) :: leap
  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r
  integer ( kind = 4 ), allocatable, dimension ( : ) :: seed
  integer ( kind = 4 ) step

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HAMMERSLEY_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate a Hammersley dataset.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program is meant to be used interactively.'
  write ( *, '(a)' ) '  It is also possible to prepare a simple input '
  write ( *, '(a)' ) '  file beforehand and use it in batch mode.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program requests input values from the user:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  * DIM_NUM, the spatial dimension,'
  write ( *, '(a)' ) '  * N, the number of points to generate,'
  write ( *, '(a)' ) '  * STEP, the index of the first subsequence element.'
  write ( *, '(a)' ) '  * SEED(1:DIM_NUM), the sequence element'
  write ( *, '(a)' ) '    corresponding to STEP = 0'
  write ( *, '(a)' ) '  * LEAP(1:DIM_NUM), the successive jumps in the sequence.'
  write ( *, '(a)' ) '  * BASE(1:DIM_NUM), the bases, usually distinct primes'
  write ( *, '(a)' ) '    or -N (to generate values like J/N).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program generates the data, writes it to the file'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    hammersley_DIM_NUM_N.txt'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where "DIM_NUM" and "N" are the numeric values specified'
  write ( *, '(a)' ) '  by the user, and then asks the user for more input.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  To indicate that no further computations are '
  write ( *, '(a)' ) '  desired, it is enough to input a nonsensical value, '
  write ( *, '(a)' ) '  such as -1.'

  do

    write ( *, '(a)' ) '  *'
    write ( *, '(a)' ) ' *'
    write ( *, '(a)' ) '*  Ready to generate a new dataset:'
    write ( *, '(a)' ) ' *'
    write ( *, '(a)' ) '  *'
    write ( *, '(a)' ) '  Enter DIM_NUM, the spatial dimension:'
    write ( *, '(a)' ) '  (Try "2" if you have no preference.)'
    write ( *, '(a)' ) '  (0 or any negative value terminates execution).'

    read ( *, *, iostat = ios ) dim_num

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAMMERSLEY_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for DIM_NUM.'
      exit
    end if

    write ( *, '(a,i12)' ) '  User input DIM_NUM = ', dim_num

    if ( .not. halham_dim_num_check ( dim_num ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAMMERSLEY_DATASET'
      write ( *, '(a,i12)' ) '  The input value of DIM_NUM = ', dim_num
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter N, the number of points to generate:'
    write ( *, '(a)' ) '  (Try "25" if you have no preference.)'
    write ( *, '(a)' ) '  (0 or any negative value terminates execution).'

    read ( *, *, iostat = ios ) n

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAMMERSLEY_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for N.'
      exit
    end if

    write ( *, '(a,i12)' ) '  User input N = ', n

    if ( .not. halham_n_check ( n ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAMMERSLEY_DATASET'
      write ( *, '(a,i12)' ) '  The input value of N = ', n
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  Enter STEP, the index of the first subsequence element:'
    write ( *, '(a)' ) '  (Try "0" or "1" if you have no preference.)'
    write ( *, '(a)' ) '  (Any negative value terminates execution).'

    read ( *, *, iostat = ios ) step

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAMMERSLEY_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for STEP.'
      exit
    end if

    write ( *, '(a,i12)' ) '  User input STEP = ', step

    if ( .not. halham_step_check ( step ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAMMERSLEY_DATASET'
      write ( *, '(a,i12)' ) '  The input value of STEP = ', step
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

    allocate ( base(1:dim_num) )
    allocate ( leap(1:dim_num) )
    allocate ( r(1:dim_num,1:n) )
    allocate ( seed(1:dim_num) )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter SEED(1:DIM_NUM), the sequence index'
    write ( *, '(a)' ) '  corresponding to STEP = 0'
    write ( *, '(a)' ) '  (Try "0 0 ... 0" if you have no preference.)'
    write ( *, '(a)' ) '  (a negative value terminates execution.)'

    read ( *, *, iostat = ios ) seed(1:dim_num)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAMMERSLEY_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for SEED.'
      exit
    end if

    call i4vec_transpose_print ( dim_num, seed, '  User input ' )

    if ( .not. halham_seed_check ( dim_num, seed ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAMMERSLEY_DATASET'
      write ( *, '(a,i12)' ) '  The input value of SEED '
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      call i4vec_transpose_print ( dim_num, seed, '  SEED:' )
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter LEAP(1:DIM_NUM), the jumps in the sequence.'
    write ( *, '(a)' ) '  (Try "1 1 ... 1" if you have no preference.)'
    write ( *, '(a)' ) '  (another choice is any prime larger than all bases.)'
    write ( *, '(a)' ) '  (any value less than 1 terminates execution.)'

    read ( *, *, iostat = ios ) leap(1:dim_num)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAMMERSLEY_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for LEAP.'
      exit
    end if

    call i4vec_transpose_print ( dim_num, leap, '  User input ' )

    if ( .not. halham_leap_check ( dim_num, leap ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAMMERSLEY_DATASET'
      write ( *, '(a,i12)' ) '  The input value of LEAP '
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      call i4vec_transpose_print ( dim_num, leap, '  LEAP:' )
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter BASE(1:DIM_NUM), the bases, usually distinct'
    write ( *, '(a)' ) '  primes, but any NEGATIVE base generates '
    write ( *, '(a)' ) '  values J/|BASE|.' 
    write ( *, '(a)' ) '  (Try "-N 2 3 5 7 11 ..." if you have no preference.)'
    write ( *, '(a)' ) '  (any value of 0 or 1 terminates execution.)'

    read ( *, *, iostat = ios ) base(1:dim_num)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAMMERSLEY_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for BASE.'
      exit
    end if

    call i4vec_transpose_print ( dim_num, base, '  User input ' )

    if ( .not. hammersley_base_check ( dim_num, base ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAMMERSLEY_DATASET'
      write ( *, '(a,i12)' ) '  The input value of BASE '
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      call i4vec_transpose_print ( dim_num, base, '  BASE:' )
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

    call i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base, r )

    write ( file_out_name, '(a,i2.2,a,i5.5,a)' ) &
      'hammersley_', dim_num, '_', n, '.txt'

    call halham_write ( dim_num, n, step, seed, leap, base, r, file_out_name )

    deallocate ( base )
    deallocate ( leap )
    deallocate ( seed )
    deallocate ( r )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The data was written to the file "' &
       // trim ( file_out_name ) // '".'

  end do

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
