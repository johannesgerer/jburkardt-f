program main

!*****************************************************************************80
!
!! MAIN generates colored noise data for a sequence of values of ALPHA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) q_d
  integer ( kind = 4 ) seed_init

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COLORED_NOISE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the COLORED_NOISE library.'

  n = 128
  q_d = 1.0D+00
  alpha = 0.00D+00
  seed_init = 123456789

  do i = 0, 8
    alpha = 0.25D+00 * real ( i, kind = 8 )
    call test01 ( n, q_d, alpha, seed_init )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COLORED_NOISE_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, * ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( n, q_d, alpha, seed_init )

!*****************************************************************************80
!
!! TEST01 calls F_ALPHA with particular parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence 
!    to generate.
!
!    Input, real ( kind = 8 ) Q_D, the variance of the sequence.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the power law.
!
!    Input, integer ( kind = 4 ) SEED_INIT, the initial seed for the 
!    random number generator.
!
  implicit none

  integer   ( kind = 4 )  n

  real      ( kind = 8 )  alpha
  integer   ( kind = 4 )  i
  integer   ( kind = 4 )  ios
  character ( len = 255 ) output_filename
  integer   ( kind = 4 )  output_unit
  real      ( kind = 8 )  q_d
  integer   ( kind = 4 )  seed
  integer   ( kind = 4 )  seed_init
  real      ( kind = 8 )  x(n)

  write ( output_filename, '(a,f4.2,a)' ) "alpha_", alpha, '.txt'
!
!  Report parameters.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a,i8,a)' ) '  Generating ', n, ' sample points.'
  write ( *, '(a,g14.6)' ) '  1/F^ALPHA noise has ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  Variance is ', q_d
  write ( *, '(a,i12)' ) '  Initial random number seed = ', seed_init

  seed = seed_init

  call f_alpha ( n, q_d, alpha, seed, x )
!
!  Print no more than 10 entries of the data.
!
  call r8vec_print_part ( n, x, 10, '  Noise sample:' )
!
!  Write the data to a file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  do i = 1, n
    write ( output_unit, '(g14.6)' ) x(i)
  end do

  close ( unit = output_unit )

  write ( *, '(a)' ) '  Data written to file "' &
    // trim ( output_filename ) // '."'

  return
end
