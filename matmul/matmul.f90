program main

!*****************************************************************************80
!
!! MAIN is the main program for MATMUL.
!
!  Discussion:
!
!    MATMUL is an interactive program that sets up, carries out, and times 
!    a matrix multiplication.  
!
!    At the user's command, MATMUL can use many different algorithms,
!    matrix sizes, arithmetic types.
!
!  Complaints:
!
!    NREP doesn't do what I think I want.  Ah, I think I see.
!    NREP1 does a number of multiplications inside the timing loop.
!    NREP2 simply does the multiplication several times, timing
!      each one.
!
!    So I want the ability to 
!    *A) Fatten the loop being timed, with multiple invocations, and
!        one timed value at the end; (average would be nice)
!    OR
!    *B) Run multiple copies of the loop, printing the times.
!
!  Making a version for a given machine:
!
!    1) Set the value of "MACHINE" to be the name of your machine.
!
!    This occurs in routine INIT.
!
!    2) Make some special routines available for special machines:
!
!     For the Cray,
!       uncomment the calls to
!         MXMA
!         SAXPY
!         SDOT
!         SGEMM
!         TIMEF.
!
!     For the SGI/IRIS,
!       uncomment the calls to
!         SAXPY
!         SDOT
!         SECNDS
!         SGEMM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 100

  character ( len = 20 ) command
  logical c_show
  logical header_show
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) lda
  logical lda_show
  integer ( kind = 4 ) lenc
  character ( len = 7 ) lingo
  logical lingo_show
  logical lval
  character ( len = 10 ) machine
  logical machine_show
  logical mflops_show
  integer ( kind = 4 ) n
  logical n_show
  integer ( kind = 4 ) nhi
  integer ( kind = 4 ) ninc
  integer ( kind = 4 ) nlo
  integer ( kind = 4 ) nmult
  logical noshow
  integer ( kind = 4 ) nrep
  logical nrep_show
  character ( len = 10 ) order
  character ( len = 10 ) order_list(order_max)
  integer ( kind = 4 ) order_num
  character ( len = 82 ) output
  logical s_eqi
  logical time_show

  call timestamp ( )
!
!  Initialize the data.
!
  call init ( command, c_show, lda, lingo, lingo_show, lda_show, &
    machine, machine_show, mflops_show,  n, nhi, ninc, nlo, nmult, &
    noshow, nrep, nrep_show, n_show, order, order_list, order_max, &
    order_num, time_show )
!
!  Say hello to the user.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATMUL'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  An interactive demonstration of the speed'
  write ( *, '(a)' ) '  of matrix multiplication.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This is version 1.20'
  write ( *, '(a)' ) '  Last modified on 16 December 2002.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This is the version for ' // trim ( machine )

  do
!
!  Print the prompt.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Command?  (Type H for help)'
!
!  Read the command, but if it's blank or a comment, keep reading.
!
    do

      read ( *, '(a)', iostat = ios ) command

      if ( ios /= 0 ) then
        command = 'QUIT'
      end if

      if ( len_trim ( command ) /= 0 .and. &
        command(1:1) /= '#' ) then
        exit
      end if

    end do
!
!  Capitalize the command.
!
    call s_cap ( command )
!
!  Remove all blanks to make interpretation easier.
!
    call s_blank_delete ( command )

    ierror = 0
!
!  "H" means print out the help message.
!
    if ( s_eqi ( command(1:1), 'H' ) ) then
 
      call help
!
!  "LDA = " means the user wants to set lda.
!
    else if ( s_eqi ( command(1:4), 'LDA=' ) ) then
 
      call s_to_i4 ( command(5:), itemp, ierror, lchar )
 
      if ( ierror /= 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  MATMUL did not understand your definition of LDA.'

      else if ( itemp <= 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The assignment of LDA was not acceptable!'
        write ( *, '(a)' ) '  LDA must be positive.'

      else

        lda = itemp
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) 'LDA has been set to ', lda

      end if
 
      if ( lda < n ) then
        n = lda
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Note:'
        write ( *, '(a)' ) '  Since N must be no greater than LDA,'
        write ( *, '(a)' ) '  MATMUL has decreased the value of N.'
        write ( *, '(a,i6)' ) '  N has been set to ', n
      end if
!
!  "M" means the user wants the multiplication to be carried out.
!
    else if ( s_eqi ( command(1:1), 'M' ) ) then
!
!  Carry out multiplication for one, or many values of N.
!
      n = nlo
 
      header_show = .true.

      do
 
        call mult ( c_show, mflops_show, header_show, lda, lingo, lingo_show, &
          lda_show, machine, machine_show, n, noshow, nrep, nrep_show, &
          n_show, order, order_list, order_num, output, time_show )

        call n_step ( ido, n, nhi, ninc, nlo, nmult )

        if ( ido /= 1 ) then
          exit
        end if

        header_show = .false.

      end do
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix multiplication has been carried out.'
!
!  "N=" means the user wants to set the matrix size N;
!  It is also possible to set a range of sizes, and an increment or multiplier.
!
    else if ( s_eqi ( command(1:2), 'N=' ) ) then
 
      call n_get ( command(3:), ierror, lda, n, nhi, ninc, nlo, nmult )

      if ( lda < n ) then
        lda = n
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Because N was changed, LDA is automatically'
        write ( *, '(a,i6)' ) '  increased to LDA = ', lda
      end if
!
!  "NOSHOW" means the user wants to turn off the display of all quantities.
!
    else if ( s_eqi ( command, 'NOSHOW' ) ) then
 
      command = 'NOSHOW=ALL'
      lval = .false.

      call getsho ( command(8:), c_show, lingo_show, lda_show, &
        lval, machine_show, mflops_show, noshow, nrep_show, n_show, &
        time_show )
!
!  "NOSHOW=" means the user wants to turn off the display of a
!  particular quantity.
!
    else if ( s_eqi ( command(1:7), 'NOSHOW=' ) ) then
 
      lval = .false.
      call getsho ( command(8:), c_show, lingo_show, lda_show, &
        lval, machine_show, mflops_show, noshow, nrep_show, n_show, &
        time_show )
!
!  "NREP=" sets the number of repetitions.
!
    else if ( s_eqi ( command(1:5), 'NREP=' ) ) then
 
      call s_to_i4 ( command(6:), nrep, ierror, lchar )
      write ( *, '(a,i6)' ) '  The repetition factor is now NREP = ', nrep
 
      if ( nrep == 1 ) then
        nrep_show = .false.
      else
        nrep_show = .true.
      end if
!
!  "ORDER=" means the user wants to set the method.
!
    else if ( s_eqi ( command(1:5), 'ORDER' ) ) then
 
      i = index ( command, '=' ) 

      if ( i == 0 ) then
        call order_list_print ( order_list, order_num )
        write ( *, '(a)' ) 'Enter your choice for the order'
        read ( *, '(a)' ) command
        call order_get ( command, ierror, order, order_list, order_num )
      else
        command = adjustl ( command(i+1:) )
        call order_get ( command, ierror, order, order_list, order_num )
      end if
!
!  "Q" means the user wants to quit.
!
    else if ( s_eqi ( command(1:1), 'Q' ) ) then
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Type "Y" to confirm that you want to quit.'

      read ( *, '(a)', iostat = ios ) command

      if ( ios /= 0 ) then
        command = 'Y'
      end if

      call s_cap ( command )

      if ( s_eqi ( command(1:1), 'Y' ) ) then
        exit
      end if
!
!  "P" means the user wants to print out the current settings.
!
    else if ( s_eqi ( command(1:1), 'P' ) ) then
 
      call printr ( lda, lingo, n, nhi, ninc, nlo, nmult, nrep, order )
!
!  "SHOW" means the user wants all items to be displayed.
!
    else if ( s_eqi ( command, 'SHOW' ) ) then
 
      command = 'SHOW=ALL'
      lval = .true.

      call getsho ( command(6:), c_show, lingo_show, lda_show, &
        lval, machine_show, mflops_show, noshow, nrep_show, n_show, &
        time_show )
!
!  "SHOW=" means the user wants a particular item displayed.
!
    else if ( s_eqi(command(1:5), 'SHOW=' ) ) then
 
      lval = .true.

      call getsho ( command(6:), c_show, lingo_show, lda_show, &
        lval, machine_show, mflops_show, noshow, nrep_show, n_show, &
        time_show )
!
!  The user's input did not match any acceptable command.
!
    else

      if ( 0 < len_trim ( command ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Your command was not recognized.'
        write ( *, '(a)' ) '  You typed ' // trim ( command)
        write ( *, '(a)' ) '  Type HELP for a list of commands.'
      end if

    end if
 
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATMUL:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine c4_ijk ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! C4_IJK uses C4 arithmetic and IJK order.
!
!  Discussion:
!
!    C4 arithmetic uses "single precision" complex numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) b(lda,n)
  complex ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call c4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do j = 1, n
        do k = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine c4_matmul ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! C4_MATMUL uses C4 arithmetic and FORTRAN90 MATMUL.
!
!  Discussion:
!
!    C4 arithmetic uses "single precision" complex numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) b(lda,n)
  complex ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call c4_set ( a, b, c, n, n )
 
    call matmul_cpu_timer ( time1 )
 
    a(1:n,1:n) = matmul ( b(1:n,1:n), c(1:n,1:n) )
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine c4_set ( a, b, c, lda, n )

!*****************************************************************************80
!
!! C4_SET initializes the A, B and C matrices using C4 arithmetic.
!
!  Discussion:
!
!    C4 arithmetic uses "single precision" complex numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, complex ( kind = 4 ) A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
!    used in the multiplication.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) b(lda,n)
  complex ( kind = 4 ) c(lda,n)

  a(1:n,1:n) = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )
  b(1:n,1:n) = cmplx ( 2.0E+00, 1.0E+00, kind = 4 )
  c(1:n,1:n) = cmplx ( 1.0E+00, 1.0E+00, kind = 4 )
 
  return
end
subroutine domethod ( lda, n, nrep, order, ttime )

!*****************************************************************************80
!
!! DOMETHOD calls a specific multiplication routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Input, character ( len = 10 ) ORDER, specifies the method to be used.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrep
  character ( len = 10 ) order
  logical s_eqi
  real ( kind = 8 ) ttime

  if ( s_eqi ( order, 'C4_IJK' ) ) then

    call c4_ijk ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'C4_MATMUL' ) ) then

    call c4_matmul ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R8_IJK' ) ) then

    call r8_ijk ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R8_MATMUL' ) ) then

    call r8_matmul ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'I4_IJK' ) ) then

    call i4_ijk ( lda, n, nrep, ttime )
 
  else if ( s_eqi ( order, 'I4_MATMUL' ) ) then

    call i4_matmul ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'L4_IJK' ) ) then

    call l4_ijk ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_DOT_PRODUCT' ) ) then

    call r4_dot_product ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_IJ' ) ) then
 
    call r4_ij ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_IJK' ) ) then
 
    call r4_ijk ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_IJK_IMPLICIT' ) ) then
 
    call r4_ijk_implicit ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_IJK_S' ) ) then

    call r4_ijk_s ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_IJK_I2' ) ) then

    call r4_ijk_i2 ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_IJK_I4' ) ) then

    call r4_ijk_i4 ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_IJK_I8' ) ) then

    call r4_ijk_i8 ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_IJK_J4' ) ) then
 
    call r4_ijk_j4 ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_IJK_K4' ) ) then

    call r4_ijk_k4 ( lda, n, nrep, ttime )
 
  else if ( s_eqi ( order, 'R4_IJK_M' ) ) then

    call r4_ijk_m ( lda, n, nrep, ttime )
 
  else if ( s_eqi ( order, 'R4_IKJ' ) ) then
 
    call r4_ikj ( lda, n, nrep, ttime )
 
  else if ( s_eqi ( order, 'R4_IKJ_DOT' ) ) then
 
    call r4_ikj_dot ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_JIK' ) ) then
 
    call r4_jik ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_JIK_IMPLICIT' ) ) then
 
    call r4_jik_implicit ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_JKI' ) ) then
 
    call r4_jki ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_JKI_IMPLICIT' ) ) then
 
    call r4_jki_implicit ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_KIJ' ) ) then

    call r4_kij ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_KIJ_DOT' ) ) then

    call r4_kij_dot ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_KJI' ) ) then

    call r4_kji ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_KJI_IMPLICIT' ) ) then

    call r4_kji_implicit ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_KJI_M' ) ) then

    call r4_kji_m ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_MATMUL' ) ) then
 
    call r4_matmul ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_MXMA' ) ) then

    call r4_mxma ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_SAXPYC' ) ) then

    call r4_saxpyc ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_SAXPYR' ) ) then

    call r4_saxpyr ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_SDOT' ) ) then

    call r4_sdot ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_SGEMM' ) ) then

    call r4_sgemm ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_TAXPYC' ) ) then

    call r4_taxpyc ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_TAXPYR' ) ) then

    call r4_taxpyr ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_TDOT' ) ) then

    call r4_tdot ( lda, n, nrep, ttime )

  else if ( s_eqi ( order, 'R4_TGEMM' ) ) then

    call r4_tgemm ( lda, n, nrep, ttime )

  end if

  return
end
subroutine getsho ( string, c_show, lingo_show, lda_show, lval, &
  machine_show, mflops_show, noshow, nrep_show, n_show, time_show )

!*****************************************************************************80
!
!! GETSHO determines what items the user wishes to print out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character STRING*(*), a string which specifies which
!    variable is to be assigned the input value of LVAL.
!
!    Output, logical C_SHOW, is TRUE if the multiplication method
!    is to be shown.
!
!    Output, logical LINGO_SHOW, is TRUE if the software's programming
!    language, stored in LINGO, is to be shown.
!
!    Output, logical LDA_SHOW, is TRUE if the variable LDA is to be shown.
!
!    Input, logical LVAL, a logical value, which is to be assigned
!    to one of the variables.
!
!    Output, logical MACHINE_SHOW, is TRUE if the machine name is to be shown.
!
!    Output, logical MFLOPS_SHOW, is TRUE if the MegaFLOP rate is to be shown.
!
!    Output, logical NOSHOW, is TRUE if the number of operations is
!    to be shown.
!
!    Output, logical NREP_SHOW, is TRUE if the number of repetitions is
!    to be shown.
!
!    Output, logical N_SHOW, is TRUE if the matrix size N is to be shown.
!
!    Output, logical TIME_SHOW, is TRUE if the execution time is to be shown.
!
  implicit none

  logical c_show
  logical lingo_show
  logical lda_show
  logical lval
  logical machine_show
  logical mflops_show
  logical noshow
  logical nrep_show
  logical n_show
  logical s_eqi
  character ( len = * ) string
  logical time_show

  if ( s_eqi ( string(1:3), 'ALL' ) ) then
    c_show = lval
    lingo_show = lval
    lda_show = lval
    machine_show = lval
    mflops_show = lval
    noshow = lval
    nrep_show = lval
    n_show = lval
    time_show = lval
  else if ( s_eqi ( string(1:3), 'CPU' ) ) then
    time_show = lval
  else if ( s_eqi ( string(1:8), 'LANGUAGE' ) ) then
    lingo_show = lval
  else if ( s_eqi ( string(1:3), 'LDA' ) ) then
    lda_show = lval
  else if ( s_eqi ( string(1:7), 'MACHINE' ) ) then
    machine_show = lval
  else if ( s_eqi ( string(1:6), 'MFLOPS' ) ) then
    mflops_show = lval
  else if ( s_eqi ( string(1:1), 'N' ) ) then
    n_show = lval
  else if ( s_eqi ( string(1:4), 'NREP' ) ) then
    nrep_show = lval
  else if ( s_eqi ( string(1:3), 'OPS' ) ) then
    noshow = lval
  else if ( s_eqi ( string(1:5), 'ORDER' ) ) then
    c_show = lval
  else if ( s_eqi ( string(1:4), 'TIME' ) ) then
    time_show = lval
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  That is not a legal name!'
    write ( *, '(a)' ) '  Legal names are CPU, LANGUAGE, LDA, MACHINE, MFLOPS,'
    write ( *, '(a)' ) '  N, NREP, OPS, ORDER, and TIME.'
  end if
 
  return
end
subroutine header ( c_show, lingo_show, lda_show, machine_show, &
  mflops_show, noshow, nrep_show, n_show, output, time_show )

!*****************************************************************************80
!
!! HEADER prints out a header for the results.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical C_SHOW, is TRUE if the multiplication method
!    is to be shown.
!
!    Input, logical LINGO_SHOW, is TRUE if the software's programming
!    language, stored in LINGO, is to be shown.
!
!    Input, logical LDA_SHOW, is TRUE if the variable LDA is to be shown.
!
!    Input, logical MACHINE_SHOW, is TRUE if the machine name is to be shown.
!
!    Input, logical MFLOPS_SHOW, is TRUE if the MegaFLOP rate is to be shown.
!
!    Input, logical NOSHOW, is TRUE if the number of operations is
!    to be shown.
!
!    Input, logical NREP_SHOW, is TRUE if the number of repetitions is
!    to be shown.
!
!    Input, logical N_SHOW, is TRUE if the matrix size N is to be shown.
!
!    Output, character ( len = 82 ) OUTPUT, a string containing a header for
!    all the variables which are to be printed.
!
!    Input, logical TIME_SHOW, is TRUE if the execution time is to be shown.
!
  implicit none

  logical c_show
  logical lingo_show
  logical lda_show
  logical machine_show
  logical mflops_show
  integer ( kind = 4 ) next
  logical noshow
  logical nrep_show
  logical n_show
  character ( len = 82 ) output
  logical time_show
!
!  Prepare the header string.  
!
  next = 1
  
  if ( c_show ) then
    output(next:) = '     Order'
    next = len_trim(output) + 1
  end if

  if ( lda_show ) then
    output(next:) = ' LDA'
    next = len_trim(output) + 1
  end if

  if ( n_show ) then
    output(next:) = '   N'
    next = len_trim(output) + 1
  end if

  if ( c_show ) then
    output(next:) = '      CPU'
    next = len_trim(output) + 1
  end if

  if ( time_show ) then
    output(next:) = ' Secs'
    next = len_trim(output) + 1
  end if

  if ( noshow ) then
    output(next:) = '       Ops'
    next = len_trim(output) + 1
  end if

  if ( nrep_show ) then
    output(next:) = ' NREP'
    next = len_trim(output) + 1
  end if

  if ( mflops_show ) then
    output(next:) = '    MFLOPS'
    next = len_trim(output) + 1
  end if

  if ( machine_show ) then
    output(next:) = '  Machine'
    next = len_trim(output) + 1
  end if
 
  if ( lingo_show ) then
    output(next:) = '  Language'
    next = len_trim(output) + 1
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) output(1:next-1)
  write ( *, '(a)' ) ' '
 
  return
end
subroutine help

!*****************************************************************************80
!
!! HELP prints a list of the available commands.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 May 2000
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This is the list of legal commands.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  H                 Help. List the legal commands.'
  write ( *, '(a)' ) '  LDA=value         Set leading dimension of arrays.'
  write ( *, '(a)' ) '  M                 Multiply two matrices.'
  write ( *, '(a)' ) '  N=value           Assigns the size of the arrays.'
  write ( *, '(a)' ) '  N=nlo,nhi,ninc    Sets N=nlo, N=nlo+ninc, ....'
  write ( *, '(a)' ) '  N=nlo,nhi,*nmult  Sets N=nlo, N=nlo*nmult, ....'
  write ( *, '(a)' ) '  NREP=nrep         Sets the repetition factor.'
  write ( *, '(a)' ) '  ORDER=name        Chooses the algorithm.'
  write ( *, '(a)' ) '  P                 Prints out current results.'
  write ( *, '(a)' ) '  Q                 Quit.'
  write ( *, '(a)' ) '  SHOW=name         Include "name" in output.'
  write ( *, '(a)' ) '  		  "name" = CPU, LANGUAGE, LDA,'
  write ( *, '(a)' ) '                    MACHINE, MFLOPS, N, NREP, '
  write ( *, '(a)' ) '                    OPS, ORDER, or TIME.'
  write ( *, '(a)' ) '                    "SHOW=ALL" means all items are shown.'
  write ( *, '(a)' ) '  NOSHOW=name       removes an item from output list.'
 
  return
end
subroutine init ( command, c_show, lda, lingo, lingo_show, lda_show, &
  machine, machine_show, mflops_show, n, nhi, ninc, nlo, nmult, noshow, &
  nrep, nrep_show, n_show, order, order_list, order_max, order_num, time_show )

!*****************************************************************************80
!
!! INIT initializes data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 6 ) COMMAND, the most recent user command.
!
!    Output, logical C_SHOW, is TRUE if the multiplication method
!    is to be shown.
!
!    Output, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Output, character ( len = 7 ) LINGO, the language in which MATMUL 
!    is written.
!
!    Output, logical LINGO_SHOW, is TRUE if the software's programming
!    language, stored in LINGO, is to be shown.
!
!    Output, logical LDA_SHOW, is TRUE if the variable LDA is to be shown.
!
!    Output, character ( len = 10 ) MACHINE, the computer where MATMUL 
!    is running.
!
!    Output, logical MACHINE_SHOW, is TRUE if the machine name is to be shown.
!
!    Output, logical MFLOPS_SHOW, is TRUE if the MegaFLOP rate is to be shown.
!
!    Output, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Output, integer ( kind = 4 ) NHI, the maximum value of N to use.
!
!    Output, integer ( kind = 4 ) NINC, the additive increment to use, if additive
!    steps are being taken.
!
!    Output, integer ( kind = 4 ) NLO, the smallest value of N to use.
!
!    Output, integer ( kind = 4 ) NMULT, the multiplier to use, if multiplicative
!    steps are being taken.
!
!    Output, logical NOSHOW, is TRUE if the number of operations is
!    to be shown.
!
!    Output, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, logical NREP_SHOW, is TRUE if the number of repetitions is
!    to be shown.
!
!    Output, logical N_SHOW, is TRUE if the matrix size N is to be shown.
!
!    Output, character ( len = 10 ) ORDER, specifies the method to be used.
!
!    Output, character ( len = 10 ) ORDER_LIST(ORDER_NUM), the list of
!    available methods.
!
!    Output, integer ( kind = 4 ) ORDER_NUM, the number of available methods.
!
!    Output, logical TIME_SHOW, is TRUE if the execution time is to be shown.
!
  implicit none

  integer ( kind = 4 ) order_max

  character ( len = 6 ) command
  logical c_show
  integer ( kind = 4 ) lda
  character ( len = 7 ) lingo
  logical lingo_show
  logical lda_show
  character ( len = 10 ) machine
  logical machine_show
  logical mflops_show
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nhi
  integer ( kind = 4 ) ninc
  integer ( kind = 4 ) nlo
  integer ( kind = 4 ) nmult
  logical noshow
  integer ( kind = 4 ) nrep
  logical nrep_show
  logical n_show
  character ( len = 10 ) order
  character ( len = 10 ) order_list(order_max)
  integer ( kind = 4 ) order_num
  logical s_eqi
  logical time_show

  command = ' '
  c_show = .true.
  lda = 16
  lda_show = .true.
  lingo = 'F90'
  lingo_show = .true.

! machine = 'Alpha'
! machine = 'CM-2'
! machine = 'Cray C90'
! machine = 'Cray YMP'
! machine = 'DECstation'
! machine = 'IBM PC'
! machine = 'Macintosh'
! machine = 'Mac/6881'
  machine = 'Mac/G5'
! machine = 'SGI/IRIS'
! machine = 'SUN'
! machine = 'VAX/VMS'
! machine = 'VECTOR/VAX'

  machine_show = .true.
  mflops_show = .true.
  n = 16
  n_show = .true.
  nhi = n
  ninc = 1
  nlo = n
  nmult = 1
  noshow = .true.
  nrep = 1
  nrep_show = .false.

  order = 'R_IJK'

  order_num = 0

  order_num = order_num + 1
  order_list(order_num) = 'C4_IJK'
  order_num = order_num + 1
  order_list(order_num) = 'C4_MATMUL'
  order_num = order_num + 1
  order_list(order_num) = 'I4_IJK'
  order_num = order_num + 1
  order_list(order_num) = 'I4_MATMUL'
  order_num = order_num + 1
  order_list(order_num) = 'L4_IJK'
  order_num = order_num + 1
  order_list(order_num) = 'R4_DOT_PRODUCT'
  order_num = order_num + 1
  order_list(order_num) = 'R4_IJ'
  order_num = order_num + 1
  order_list(order_num) = 'R4_IJK'
  order_num = order_num + 1
  order_list(order_num) = 'R4_IJK_IMPLICIT'
  order_num = order_num + 1
  order_list(order_num) = 'R4_IJK_I2'
  order_num = order_num + 1
  order_list(order_num) = 'R4_IJK_I4'
  order_num = order_num + 1
  order_list(order_num) = 'R4_IJK_I8'
  order_num = order_num + 1
  order_list(order_num) = 'R4_IJK_J4'
  order_num = order_num + 1
  order_list(order_num) = 'R4_IJK_K4'

  if ( s_eqi ( machine(1:4), 'CRAY') ) then
    order_num = order_num + 1
    order_list(order_num) = 'R4_IJK_M'
  end if

  if ( s_eqi ( machine(1:4), 'CRAY') ) then
    order_num = order_num + 1
    order_list(order_num) = 'R4_IJK_S'
  end if

  order_num = order_num + 1
  order_list(order_num) = 'R4_IKJ'
  order_num = order_num + 1
  order_list(order_num) = 'R4_IKJ_DOT'
  order_num = order_num + 1
  order_list(order_num) = 'R4_JIK'
  order_num = order_num + 1
  order_list(order_num) = 'R4_JIK_IMPLICIT'
  order_num = order_num + 1
  order_list(order_num) = 'R4_JKI'
  order_num = order_num + 1
  order_list(order_num) = 'R4_JKI_IMPLICIT'
  order_num = order_num + 1
  order_list(order_num) = 'R4_KIJ'
  order_num = order_num + 1
  order_list(order_num) = 'R4_KIJ_DOT'
  order_num = order_num + 1
  order_list(order_num) = 'R4_KJI'
  order_num = order_num + 1
  order_list(order_num) = 'R4_KJI_IMPLICIT'

  if ( s_eqi ( machine, 'SGI/IRIS' ) ) then
    order_num = order_num + 1
    order_list(order_num) = 'R4_KJI_M'
  end if

  order_num = order_num + 1
  order_list(order_num) = 'R4_MATMUL'

  if ( s_eqi ( machine(1:4), 'CRAY') ) then
    order_num = order_num + 1
    order_list(order_num) = 'R4_MXMA'
  end if

  if ( s_eqi ( machine(1:4), 'CRAY') .or. s_eqi ( machine, 'SGI/IRIS' ) ) then
    order_num = order_num + 1
    order_list(order_num) = 'R4_SAXPYC' 
  end if

  if ( s_eqi ( machine(1:4), 'CRAY') .or. s_eqi ( machine, 'SGI/IRIS' ) ) then
    order_num = order_num + 1
    order_list(order_num) = 'R4_SAXPYR' 
  end if

  if ( s_eqi ( machine(1:4), 'CRAY') .or. s_eqi ( machine, 'SGI/IRIS' ) ) then
    order_num = order_num + 1
    order_list(order_num) = 'R4_SDOT'
  end if

  if ( s_eqi ( machine(1:4), 'CRAY') .or. s_eqi ( machine, 'SGI/IRIS' ) ) then
    order_num = order_num + 1
    order_list(order_num) = 'R4_SGEMM'
  end if

  order_num = order_num + 1
  order_list(order_num) = 'R4_TAXPYC'
  order_num = order_num + 1
  order_list(order_num) = 'R4_TAXPYR'
  order_num = order_num + 1
  order_list(order_num) = 'R4_TDOT'
  order_num = order_num + 1
  order_list(order_num) = 'R4_TGEMM'

  order_num = order_num + 1
  order_list(order_num) = 'R8_IJK'
  order_num = order_num + 1
  order_list(order_num) = 'R8_MATMUL'

  time_show = .true.
 
  return
end
subroutine i4_ijk ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! I4_IJK uses I4 arithmetic and IJK order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) b(lda,n)
  integer ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call i4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do j = 1, n
        do k = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine i4_matmul ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! I4_MATMUL uses I4 arithmetic and FORTRAN90 MATMUL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) b(lda,n)
  integer ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call i4_set ( a, b, c, n, n )
 
    call matmul_cpu_timer ( time1 )
 
    a(1:n,1:n) = matmul ( b(1:n,1:n), c(1:n,1:n) )
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine i4_set ( a, b, c, lda, n )

!*****************************************************************************80
!
!! I4_SET initializes the A, B and C matrices using I4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Workspace, integer ( kind = 4 ) A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
!    used in the multiplication.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) b(lda,n)
  integer ( kind = 4 ) c(lda,n)

  a(1:n,1:n) = 0
  b(1:n,1:n) = 1
  c(1:n,1:n) = 1
 
  return
end
subroutine l4_ijk ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! L4_IJK uses L4 arithmetic and IJK order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  logical a(lda,n)
  logical b(lda,n)
  logical c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call l4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do j = 1, n
        do k = 1, n
          a(i,k) = a(i,k) .or. ( b(i,j) .and. c(j,k) )
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine l4_set ( a, b, c, lda, n )

!*****************************************************************************80
!
!! L4_SET initializes the A, B and C matrices using L4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
!    used in the multiplication.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  logical a(lda,n)
  logical b(lda,n)
  logical c(lda,n)

  a(1:n,1:n) = .false.
  b(1:n,1:n) = .true.
  c(1:n,1:n) = .true.
 
  return
end
subroutine matmul_cpu_timer ( cpu )

!*****************************************************************************80
!
!! MATMUL_CPU_TIMER computes total CPU seconds.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) CPU, the total CPU time, in seconds, since the
!    program began running.
!
  implicit none

  real ( kind = 8 ) cpu

  call cpu_time ( cpu )

  return
end
subroutine matmul_real_timer ( seconds )

!*****************************************************************************80
!
!! MATMUL_REAL_TIMER returns a reading of the real time clock.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) SECONDS, a current real time in seconds.
!
  implicit none

  integer ( kind = 4 ) clock_count
  integer ( kind = 4 ) clock_max
  integer ( kind = 4 ) clock_rate
  real ( kind = 4 ) seconds

  call system_clock ( clock_count, clock_rate, clock_max )

  seconds = real ( clock_count ) / real ( clock_rate )

  return
end
subroutine mult ( c_show, mflops_show, header_show, lda, lingo, lingo_show, &
  lda_show, machine, machine_show, n, noshow, nrep, nrep_show, n_show, &
  order, order_list, order_num, output, time_show )

!*****************************************************************************80
!
!! MULT carries out the matrix multiplication, using the requested method.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical C_SHOW, is TRUE if the multiplication method
!    is to be shown.
!
!    Input, logical MFLOPS_SHOW, is TRUE if the MegaFLOP rate is to be shown.
!
!    Input, logical HEADER_SHOW, is TRUE if the header is to be printed.
!
!    Input/output, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!    Automatically reset to N if it is less than N.
!
!    Input, character ( len = 7 ) LINGO, the language in which MATMUL is
!    written.
!
!    Input, logical LINGO_SHOW, is TRUE if the software's programming
!    language, stored in LINGO, is to be shown.
!
!    Input, logical LDA_SHOW, is TRUE if the variable LDA is to be shown.
!
!    Input, character ( len = 10 ) MACHINE, the computer where MATMUL is
!    running.
!
!    Input, logical MACHINE_SHOW, is TRUE if the machine name is to be shown.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, logical NOSHOW, is TRUE if the number of operations is
!    to be shown.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Input, logical NREP_SHOW, is TRUE if the number of repetitions is
!    to be shown.
!
!    Input, logical N_SHOW, is TRUE if the matrix size N is to be shown.
!
!    Input, character ( len = 10 ) ORDER, specifies the method to be used.
!
!    Input, character ( len = 10 ) ORDER_LIST(ORDER_NUM), the list of
!    available methods.
!
!    Input, integer ( kind = 4 ) ORDER_NUM, the number of available methods.
!
!    Output, character ( len = 82 ) OUTPUT, a string containing the data for
!    this multiplication.
!
!    Input, logical TIME_SHOW, is TRUE if the execution time is to be shown.
!
  implicit none

  integer ( kind = 4 ) order_num

  logical c_show
  logical mflops_show
  logical header_show
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lda
  logical lda_show
  character ( len = 7 ) lingo
  logical lingo_show
  character ( len = 10 ) machine
  logical machine_show
  integer ( kind = 4 ) n
  logical noshow
  integer ( kind = 4 ) nrep
  logical nrep_show
  logical n_show
  character ( len = 10 ) order
  character ( len = 10 ) order_list(order_num)
  character ( len = 82 ) output
  logical s_eqi
  logical time_show
  real ( kind = 8 ) ttime

  if ( lda < n ) then
    lda = n
  end if

  if ( header_show ) then

    call header ( c_show, lingo_show, lda_show, machine_show, &
      mflops_show, noshow, nrep_show, n_show, output, time_show )

  end if

  if ( s_eqi ( order, 'ALL' ) ) then

    do i = 1, order_num

      order = order_list(i)

      call domethod ( lda, n, nrep, order, ttime )

      call report ( c_show, mflops_show, lda, lda_show, lingo, lingo_show, machine, &
        machine_show, n, noshow, nrep, nrep_show, n_show, order, output, &
        time_show, ttime )
 
    end do

    order = 'ALL'

  else

    call domethod ( lda, n, nrep, order, ttime )

    call report ( c_show, mflops_show, lda, lda_show, lingo, lingo_show, &
      machine, machine_show, n, noshow, nrep, nrep_show, n_show, order, &
      output, time_show, ttime )

  end if

  return
end
subroutine n_get ( string, ierror, lda, n, nhi, ninc, nlo, nmult )

!*****************************************************************************80
!
!! N_GET determines the problem sizes desired by the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, a string containing the user's
!    data for an "N" command.  This command might have one of the
!    following forms, where "n", "nlo", "nhi", "ninc" and "nmult"
!    are actually numeric values:
!      n                   Solve a single problem of size N.
!      nlo, nhi            Solve every problem from N = NLO to N = NHI.
!      nlo, nhi, ninc      Solve from N = NLO to N = NHI, incrementing by NINC.
!      nlo, nhi, *nmult    Solve from N = NLO to N = NHI, multiplying by NMULT.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error occurred.
!    nonzero, an error occurred, and the operation could not be done.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Output, integer ( kind = 4 ) N, the new value for the number of rows and columns 
!    to use in the matrices.
!
!    Output, integer ( kind = 4 ) NHI, the maximum value of N to use.
!
!    Output, integer ( kind = 4 ) NINC, the additive increment to use, if additive
!    steps are being taken, or 0.
!
!    Output, integer ( kind = 4 ) NLO, the smallest value of N to use.
!
!    Output, integer ( kind = 4 ) NMULT, the multiplier to use, if multiplicative
!    steps are being taken, or 0.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nhi
  integer ( kind = 4 ) ninc
  integer ( kind = 4 ) nlo
  integer ( kind = 4 ) nmult
  character ( len = * ) string

  nhi = 0
  ninc = 0
  nmult = 0
  nlo = 0
!
!  Read the first number, N or NLO.
!
  call s_to_i4 ( string, itemp, ierror, lchar )
  next = lchar + 1

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I could not understand your definition of N.'
    return
  end if

  if ( itemp <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This value of N is not acceptable!'
    write ( *, '(a)' ) '  N must be positive.'
    return
  end if

  n = itemp
  nlo = itemp
  nhi = itemp
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  N has been set to ', n
  write ( *, '(a,i6)' ) '  NLO has been set to ', nlo

  if ( len_trim ( string ) < next ) then
    write ( *, '(a,i6)' ) '  NHI has been set to ', nhi
    write ( *, '(a,i6)' ) '  NINC has been set to ', ninc
    write ( *, '(a,i6)' ) '  NMULT has been set to ', nmult
    return
  end if
!
!  Assume the next character is a comma, so wipe it out.
!
  if ( string(next:next) == ',' ) then
    next = next + 1
  end if
!
!  Read the second number, NHI.
!
  call s_to_i4 ( string(next:), itemp, ierror, lchar )
  next = next + lchar

  if ( ierror /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  I could not understand your definition of NHI.'
    write ( *, '(a)' ) '  "' // trim ( string(next:) ) // '".'
    return

  else if ( itemp <= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This value of NHI is not acceptable!'
    write ( *, '(a)' ) '  NHI must be positive.'
    return

  else

    nhi = itemp
    write ( *, '(a,i6)' ) '  NHI has been set to ', nhi

  end if

  if ( len_trim ( string ) < next ) then
    ninc = 1
    write ( *, '(a,i6)' ) '  NINC has been set to ', ninc
    write ( *, '(a,i6)' ) '  NMULT has been set to ', nmult
    return
  end if
!
!  Assume the next character is a comma, so wipe it out.
!
  if ( string(next:next) == ',' ) then
    next = next + 1
  end if

  if ( string(next:next) == '*' ) then

    next = next + 1
    call s_to_i4 ( string(next:), itemp, ierror, lchar )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I could not understand your definition of NMULT.'
      write ( *, '(a)' ) '  "' // trim ( string(next:) ) // '".'
    else
      ninc = 0
      nmult = itemp
      write ( *, '(a,i6)' ) '  NINC has been set to ', ninc
      write ( *, '(a,i6)' ) '  NMULT has been set to ', nmult
    end if

  else

    call s_to_i4 ( string(next:), itemp, ierror, lchar )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I could not understand your definition of NINC.'
      write ( *, '(a)' ) '  "' // trim ( string(next:) ) // '".'
    else
      ninc = itemp
      nmult = 0
      write ( *, '(a,i6)' ) '  NINC has been set to ', ninc
      write ( *, '(a,i6)' ) '  NMULT has been set to ', nmult
    end if

  end if
 
  return
end
subroutine n_step ( ido, n, nhi, ninc, nlo, nmult )

!*****************************************************************************80
!
!! N_STEP is used when a set of values of N is being generated.
!
!  Discussion:
!
!    The routine checks whether addition or multiplication is being used,
!    and increases the value of N.  It also checks whether the set of
!    values is done, or whether the input values are inconsistent.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IDO.
!    * 0, N has reached or surpassed NHI, and no further
!    increase should be carried out.  N has been reset to NLO.
!    * 1, N had not reached or surpassed NHI, and so N has been
!    incremented by NINC, or multiplied by NMULT.
!    * 2, N had not reached or surpassed NHI, but we're never
!    going to get there!  Either:
!      NINC is negative but NHI is greater than NLO, or
!      NINC is positive but NHI is less than NLO, or
!      NMULT is positive but NLO is greater than NHI.
!    * 3, NINC is 0 and NMULT is less than or equal to 1.
!
!    Input/output, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!    This routine modifies the value of N as appropriate.
!
!    Input, integer ( kind = 4 ) NHI, the maximum value of N to use.
!
!    Input, integer ( kind = 4 ) NINC, the additive increment to use, if additive
!    steps are being taken.
!
!    Input, integer ( kind = 4 ) NLO, the smallest value of N to use.
!
!    Input, integer ( kind = 4 ) NMULT, the multiplier to use, if multiplicative
!    steps are being taken.
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nhi
  integer ( kind = 4 ) ninc
  integer ( kind = 4 ) nlo
  integer ( kind = 4 ) nmult
!
!  If NINC is not 0, then
!    if it's pointing in the right direction, then
!      add NINC to N,
!      set a continuation flag
!      if N+NINC exceeds NHI, then
!        reset N, and
!        set a completion flag
!    else
!      set an error flag.
!
  if ( ninc /= 0 ) then

    if ( ( nlo < nhi .and. 0 < ninc ) .or. &
         ( nhi < nlo .and. ninc < 0 ) ) then
      n = n + ninc
      ido = 1
      if ( ( nhi < n .and. nlo <= nhi ) .or. &
           ( n < nhi .and. nhi <= nlo ) ) then
        ido = 0
        n = nlo
      end if
    else
      ido = 2
    end if
!
!  If NMULT is greater than 1, then
!    if it's pointing in the right direction, then
!      multiply N by NMULT,
!      set a continuation flag
!      if N*NMULT exceeds NHI, then
!        reset N, and
!        set a completion flag
!    else
!      set an error flag.
!
  else if ( 1 < nmult ) then

    if ( nlo < nhi ) then
      n = n * nmult
      ido = 1
      if ( ( nhi < n .and. nlo <= nhi ) .or. &
         ( n < nhi .and. nhi <= nlo ) ) then
        ido = 0
        n = nlo
      end if
    else
      ido = 2
    end if
!
!  NINC was 0, and NMULT wasn't greater than 1.
!
  else

    ido = 3

  end if

  return
end
subroutine order_get ( string, ierror, order, order_list, order_num )

!*****************************************************************************80
!
!! ORDER_GET reads a new value of order from the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, a string of characters, containing
!    the user's choice for the matrix multiplication method to employ.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error occurred.
!    nonzero, an error occurred and the operation could not be completed.
!
!    Output, character ( len = 10 ) ORDER, specifies the method to be used.
!
!    Input, character ( len = 10 ) ORDER_LIST(ORDER_NUM), the list of
!    available methods.
!
!    Input, integer ( kind = 4 ) ORDER_NUM, the number of available methods.
!
  implicit none

  integer ( kind = 4 ) order_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  character ( len = 10 ) order
  character ( len = 10 ) order_list(order_num)
  logical s_eqi
  character ( len = * ) string

  ierror = 0

  if ( s_eqi ( string, 'ALL' ) ) then
    order = string
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The order has been set to '// trim ( order )
    return
  end if

  do i = 1, order_num
    if ( s_eqi ( string, order_list(i) ) ) then
      order = string
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The order has been set to ' // trim ( order )
      return
    end if
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The order you requested was not valid:'
  write ( *, '(a)' ) '  "' // trim ( string ) // '"'

  ierror = 1

  return
end
subroutine order_list_print ( order_list, order_num )

!*****************************************************************************80
!
!! ORDER_LIST_PRINT prints the list of choices for the algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 10 ) ORDER_LIST(ORDER_NUM), the list of
!    available methods.
!
!    Input, integer ( kind = 4 ) ORDER_NUM, the number of available methods.
!
  implicit none

  integer ( kind = 4 ) order_num

  integer ( kind = 4 ) i
  character ( len = * ) order_list(order_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Valid choices for the order are:'
  write ( *, '(a)' ) ' '
  write ( *, '(5a10)' ) ( order_list(i), i = 1, order_num )  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'or "ALL" to try them all.'

  return
end
subroutine printr ( lda, lingo, n, nhi, ninc, nlo, nmult, nrep, order )

!*****************************************************************************80
!
!! PRINTR prints out those parameters the user wants to see.
!
!  Discussion:
!
!    These parameters include:
!
!      the language MATMUL is written in,
!      the algorithm,
!      the leading dimension,
!      the maximum allowable dimension,
!      the actual size of arrays,
!      the number of multiplications carried out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, character ( len = 7 ) LINGO, the language in which MATMUL 
!    was written.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NHI, the maximum value of N to use.
!
!    Input, integer ( kind = 4 ) NINC, the additive increment to use, if additive
!    steps are being taken.
!
!    Input, integer ( kind = 4 ) NLO, the smallest value of N to use.
!
!    Input, integer ( kind = 4 ) NMULT, the multiplier to use, if multiplicative
!    steps are being taken.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Input, character ( len = 10 ) ORDER, specifies the method to be used.
!
  implicit none

  integer ( kind = 4 ) lda
  character ( len = 7 ) lingo
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nhi
  integer ( kind = 4 ) ninc
  integer ( kind = 4 ) nlo
  integer ( kind = 4 ) nmult
  integer ( kind = 4 ) nrep
  character ( len = 10 ) order

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This version of MATMUL is written in ' // trim ( lingo )
  write ( *, '(a)' ) '  The algorithm chosen is ' // trim ( order )
  write ( *, '(a,i6)' ) '  The leading dimension of arrays, LDA, is ', lda
  write ( *, '(a,i6)' ) '  The actual size of the arrays, N, is ', n

  if ( nhi /= nlo ) then
    write ( *, '(a)' ) '  Several problem sizes will be solved in order.'
    write ( *, '(a,i6)' ) '  The final size of arrays, NHI, will be ', nhi
  end if
 
  if ( ninc /= 0 ) then
    write ( *, '(a,i6)' ) '  Array size will be incremented by NINC = ', ninc
  end if
 
  if ( nmult /= 1 ) then
    write ( *, '(a,i6)' ) '  Array size will be multiplied by NMULT = ', nmult
  end if
 
  if ( nrep /= 1 ) then
    write ( *, '(a,i6)' ) '  Multiplications repeated NREP = ', nrep, ' times.'
  end if
 
  return
end
subroutine r4_dot_product ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_DOT_PRODUCT multiplies A = B*C using DOT_PRODUCT and R4 arithmetic.
!
!  Discussion:
!
!    This is equivalent to the following "IKJ" code:
!
!       do i = 1, n
!         do k = 1, n
!           do j = 1, n
!             a(i,k) = a(i,k) + b(i,j) * c(j,k)
!          end do
!        end do
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00

  do irep = 1, nrep

    call r4_set ( a, b, c, lda, n )

    call matmul_cpu_timer ( time1 )

    do i = 1, n
      do k = 1, n
        a(i,k) = dot_product ( b(i,1:n), c(1:n,k) )
      end do
    end do

    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do

  return
end
subroutine r4_matmul ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_MATMUL computes A = B*C using FORTRAN90 MATMUL and R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, n, n )
 
    call matmul_cpu_timer ( time1 )
 
    a(1:n,1:n) = matmul ( b(1:n,1:n), c(1:n,1:n) )
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r4_ij ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_IJ sets A = B*C using index order IJ with implicit K and R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00

  do irep = 1, nrep

    call r4_set ( a, b, c, lda, n )

    call matmul_cpu_timer ( time1 )

    do i = 1, n
      do j = 1, n
        a(i,1:n) = a(i,1:n) + b(i,j) * c(j,1:n)
      end do
    end do

    call matmul_cpu_timer ( time2 )

    ttime = ttime + time2 - time1

  end do

  return
end
subroutine r4_ijk ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_IJK multiplies A = B*C using index order IJK and R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do j = 1, n
        do k = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
 
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_ijk_implicit ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_IJK_IMPLICIT sets A = B*C, index order IJK, implicit loops, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do j = 1, n
        a(i,1:n) = a(i,1:n) + b(i,j) * c(j,1:n)
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
 
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_ijk_m ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_IJK_M multiplies A = B*C using index order IJK and R4 arithmetic.
!
!  Discussion:
!
!    The routine uses a Cray directive to run the triple loop using
!    multitasking.  The benefit of such a directive depends on the
!    algorithm and the load on the machine.
!
!    Except on the Cray, this routine should not be used, and in
!    particular, the call to TIMEF should be commented out.
!
!    Note that the Cray routine TIMEF must be called, rather than
!    SECOND.  TIMEF reports elapsed "real" time or "wallclock" time,
!    which should go down with multitasking, whereas CPU time should
!    remain roughly constant.
!
!    In order for parallel processing to occur, this routine must be
!    compiled on the Cray with the directive "-Zu"; moreover, the user must
!    set the environment variable NCPUS to the number of processors the
!    user would like.  For instance, a C shell user would type:
!
!      setenv NCPUS 8
!
!    while a Bourne shell user would type
!
!      NCPUS = 8
!      export NCPUS
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  time1 = 0.0E+00
  time2 = 0.0E+00
  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )

    do i = 1, n
      do j = 1, n
        do k = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    ttime = ttime + ( time2 - time1 ) / 1000.0E+00
 
  end do
 
  return
end
subroutine r4_ijk_s ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_IJK_S sets A = B*C, index order IJK, no Cray vectorization, R4 arithmetic.
!
!  Discussion:
!
!    The routine uses a Cray directive to run the inner do loop WITHOUT
!    vectorization.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do j = 1, n
        do k = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_ijk_i2 ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_IJK_I2 sets A = B*C, index order IJK, unrolling I 2 times, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  integer ( kind = 4 ), parameter :: nroll = 2
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    ihi = ( n / nroll ) * nroll
    do i = 1, ihi, nroll
      do j = 1, n
        do k = 1, n
          a(i,k)   = a(i,k)   + b(i,j)   * c(j,k)
          a(i+1,k) = a(i+1,k) + b(i+1,j) * c(j,k)
        end do
      end do
    end do
!
!  Take care of the few cases we missed if N is not a multiple of NROLL.
!
    do i = ihi+1, n
      do j = 1, n
        do k = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r4_ijk_i4 ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_IJK_I4 sets A = B*C, index order IJK, unrolling I 4 times, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  integer ( kind = 4 ), parameter :: nroll = 4
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    ihi = ( n / nroll ) * nroll
    do i = 1, ihi, nroll
      do j = 1, n
        do k = 1, n
          a(i,k)   = a(i,k)   + b(i,j)   * c(j,k)
          a(i+1,k) = a(i+1,k) + b(i+1,j) * c(j,k)
          a(i+2,k) = a(i+2,k) + b(i+2,j) * c(j,k)
          a(i+3,k) = a(i+3,k) + b(i+3,j) * c(j,k)
        end do
      end do
    end do
!
!  Take care of the few cases we missed if N is not a multiple of NROLL.
!
    do i = ihi+1, n
      do j = 1, n
        do k = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r4_ijk_i8 ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_IJK_I8 sets A = B*C, index order IJK, unrolling I 8 times, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  integer ( kind = 4 ), parameter :: nroll = 8
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    ihi = ( n / nroll ) * nroll
    do i = 1, ihi, nroll
      do j = 1, n
        do k = 1, n
          a(i,k)   = a(i,k)   + b(i,j)   * c(j,k)
          a(i+1,k) = a(i+1,k) + b(i+1,j) * c(j,k)
          a(i+2,k) = a(i+2,k) + b(i+2,j) * c(j,k)
          a(i+3,k) = a(i+3,k) + b(i+3,j) * c(j,k)
          a(i+4,k) = a(i+4,k) + b(i+4,j) * c(j,k)
          a(i+5,k) = a(i+5,k) + b(i+5,j) * c(j,k)
          a(i+6,k) = a(i+6,k) + b(i+6,j) * c(j,k)
          a(i+7,k) = a(i+7,k) + b(i+7,j) * c(j,k)
        end do
      end do
    end do
!
!  Take care of the few cases we missed if N is not a multiple of NROLL.
!
    do i = ihi+1, n
      do j = 1, n
        do k = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r4_ijk_j4 ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_IJK_J4 sets A = B*C, index order IJK, unrolling on J, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  integer ( kind = 4 ), parameter :: nroll = 4
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
    jhi = ( n / nroll ) * nroll
 
    do i = 1, n
      do j = 1, jhi, nroll
        do k = 1, n
          a(i,k) = a(i,k) + b(i,j  ) * c(j,k)   + b(i,j+1) * c(j+1,k) &
                          + b(i,j+2) * c(j+2,k) + b(i,j+3) * c(j+3,k)
        end do
      end do
    end do
!
!  Take care of the few cases we missed if N is not a multiple of NROLL.
!
    do i = 1, n
      do j = jhi+1, n
        do k = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_ijk_k4 ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_IJK_K4 sets A = B*C, index order IJK, unrolling on K, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) nrep
  integer ( kind = 4 ), parameter :: nroll = 4
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    khi = ( n / nroll ) * nroll

    do i = 1, n
      do j = 1, n
        do k = 1, khi, nroll
          a(i,k)   = a(i,k)   + b(i,j) * c(j,k)
          a(i,k+1) = a(i,k+1) + b(i,j) * c(j,k+1)
          a(i,k+2) = a(i,k+2) + b(i,j) * c(j,k+2)
          a(i,k+3) = a(i,k+3) + b(i,j) * c(j,k+3)
        end do
      end do
    end do
!
!  Take care of the few cases we missed if N is not a multiple of NROLL.
!
    do i = 1, n
      do j = 1, n
        do k = khi+1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_ikj ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_IKJ multiplies A = B*C using index order IKJ, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do k = 1, n
        do j = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_ikj_dot ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_IKJ multiplies A = B*C, index order IKJ, DOT_PRODUCT, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do k = 1, n
        a(i,k) = dot_product ( b(i,1:n), c(1:n,k) )
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_jik ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_JIK multiplies A = B*C using index order JIK, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do j = 1, n
      do i = 1, n
        do k = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_jik_implicit ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_JIK_implicit sets A = B*C, index order JIK, implicit loops, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do j = 1, n
      do i = 1, n
        a(i,1:n) = a(i,1:n) + b(i,j) * c(j,1:n)
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_jki ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_JKI multiplies A = B*C using index order JKI, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Workspace, real ( kind = 4 ) A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
!    used in the multiplication.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do j = 1, n
      do k = 1, n
        do i = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_jki_implicit ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_JKI_IMPLICIT sets A = B*C, index order JKI, implicit loops, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Workspace, real ( kind = 4 ) A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
!    used in the multiplication.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do j = 1, n
      do k = 1, n
        a(1:n,k) = a(1:n,k) + b(1:n,j) * c(j,k)
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_kij ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_KIJ multiplies A = B*C using index order KIJ, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do k = 1, n
      do i = 1, n
        do j = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_kij_dot ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_KIJ_DOT sets A = B*C using index order KIJ and DOT_PRODUCT, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do k = 1, n
      do i = 1, n
        a(i,k) = dot_product ( b(i,1:n), c(1:n,k) )
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_kji ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_KJI multiplies A = B*C using index order KJI, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do k = 1, n
      do j = 1, n
        do i = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_kji_implicit ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_KJI_IMPLICIT sets A = B*C, index order KJI, implicit loops, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do k = 1, n
      do j = 1, n
        a(1:n,k) = a(1:n,k) + b(1:n,j) * c(j,k)
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_kji_m ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_KJI_M sets A = B*C using index order KJI and multitasking, R4 arithmetic.
!
!  Discussion:
!
!    The routine uses an SGI/IRIS parallel processing directive to run the
!    triple loop using multitasking.
!
!    The benefit of such a directive depends on the algorithm and the
!    load on the machine.
!
!    Except on the SGI/IRIS, this routine should not be used, and in
!    particular, the call to SECNDS should be commented out.
!
!    Note that the SGI/IRIS routine SECNDS must be called, rather than
!    SECOND.  SECNDS reports elapsed "real" time or "wallclock" time,
!    which should go down with multitasking, whereas CPU time should
!    remain roughly constant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  time1 = 0.0E+00
  time2 = 0.0E+00
  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_real_timer ( time1 )

    do k = 1, n
      do j = 1, n
        do i = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_real_timer ( time2 )

    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_mxma ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_MXMA multiplies A = B*C using optimized MXMA, R4 arithmetic.
!
!  Discussion:
!
!    Since the routine MXMA is only available on the Cray, in the
!    SCILIB library, the statement
!
!      call mxma(...)
!
!    should be commented out in versions of MATMUL that are to run
!    on other machines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ), parameter :: inca = 1
  integer ( kind = 4 ), parameter :: incb = 1
  integer ( kind = 4 ), parameter :: incc = 1
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )

!       call mxma ( b, incb, lda, c, incc, lda, a, inca, lda, n, n, n )

    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r4_saxpyc ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_SAXPYC multiplies A = B*C columnwise, using optimized SAXPY, R4 arithmetic.
!
!  Discussion:
!
!    SAXPY is used to carry out the multiplication "columnwise".
!
!    This is equivalent to the following "JKI" code:
!
!       do j = 1, n
!         do k = 1, n
!           do i = 1, n
!              a(i,k) = a(i,k)+b(i,j)*c(j,k)
!            end do
!          end do
!        end do
!
!    Except on the Cray and SGI/IRIS, the statement "call saxpy" below
!    should be commented out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ), parameter :: inca = 1
  integer ( kind = 4 ), parameter :: incb = 1
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do j = 1, n
      do k = 1, n
        call saxpy ( n, c(j,k), b(1,j), incb, a(1,k), inca )
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r4_saxpyr ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_SAXPYR multiplies A = B*C "rowwise", using optimized SAXPY, R4 arithmetic.
!
!  Discussion:
!
!    This is equivalent to the following "IJK" code:
!
!     do i = 1, n
!       do j = 1, n
!          do k = 1, n
!            a(i,k) = a(i,k)+b(i,j)*c(j,k)
!          end do
!        end do
!      end do
!
!    Except on the Cray and SGI/IRIS, the statement "call saxpy" below
!    should be commented out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do j = 1, n
        call saxpy ( n, b(i,j), c(j,1), lda, a(i,1), lda )
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r4_sdot ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_SDOT multiplies A = B*C using optimized SDOT, R4 arithmetic.
!
!  Discussion:
!
!    This is equivalent to the following "IKJ" code:
!
!       do i = 1, n
!         do k = 1, n
!           do j = 1, n
!             a(i,k) = a(i,k)+b(i,j)*c(j,k)
!          end do
!        end do
!      end do
!
!    Except on the Cray or SGI/IRIS, the call to SDOT should be commented out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 4 ) sdot
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do k = 1, n
        a(i,k) = sdot ( n, b(i,1), lda, c(1,k), 1 )
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r4_set ( a, b, c, lda, n )

!*****************************************************************************80
!
!! R4_SET initializes the A, B and C matrices using R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
!    used in the multiplication.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)

  a(1:n,1:n) = 0.0E+00
  b(1:n,1:n) = 1.0E+00
  c(1:n,1:n) = 1.0E+00
 
  return
end
subroutine r4_sgemm ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_SGEMM multiplies A = B*C using optimized SGEMM, R4 arithmetic.
!
!  Discussion:
!
!    Except on the Cray or SGI/IRIS, the call to SGEMM should be commented out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ), parameter :: alpha = 1.0E+00
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ), parameter :: beta = 0.0E+00
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )

    call sgemm ( 'n', 'n', n, n, n, alpha, b, lda, c, lda, beta, a, lda )

    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r4_taxpyc ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_TAXPYC sets A = B*C columnwise, using source code SAXPY, R4 arithmetic.
!
!  Discussion:
!
!    This is equivalent to the following "JKI" code:
!
!       do j = 1, n
!         do k = 1, n
!           do i = 1, n
!             a(i,k) = a(i,k) + b(i,j) * c(j,k)
!           end do
!         end do
!       end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do j = 1, n
      do k = 1, n
        call taxpy ( n, c(j,k), b(1,j), 1, a(1,k), 1 )
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r4_taxpyr ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_TAXPYR multiplies A = B*C rowwise using source code SAXPY, R4 arithmetic
!
!  Discussion:
!
!    This is equivalent to the following "IJK" code:
!
!        do i = 1, n
!          do j = 1, n
!            do k = 1, n
!              a(i,k) = a(i,k) + b(i,j) * c(j,k)
!            end do
!          end do
!        end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do j = 1, n
        call taxpy ( n, b(i,j), c(j,1), lda, a(i,1), lda )
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r4_tdot ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_TDOT multiplies A = B * C using source code SDOT, R4 arithmetic.
!
!  Discussion:
!
!    This is equivalent to the following "IKJ" code:
!
!      do i = 1, n
!        do k = 1, n
!          do j = 1, n
!            a(i,k) = a(i,k) + b(i,j) * c(j,k)
!          end do
!        end do
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 4 ) tdot
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do k = 1, n
        a(i,k) = tdot ( n, b(i,1), lda, c(1,k), 1 )
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r4_tgemm ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R4_TGEMM multiplies A = B*C using source code SGEMM, R4 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ), parameter :: alpha = 1.0E+00
  real ( kind = 4 ) b(lda,n)
  real ( kind = 4 ), parameter :: beta = 0.0E+00
  real ( kind = 4 ) c(lda,n)
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r4_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )

    call tgemm ( 'n', 'n', n, n, n, alpha, b, lda, c, lda, beta, a, lda )

    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1
 
  end do
 
  return
end
subroutine r8_matmul ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R8_MATMUL computes A = B*C using FORTRAN90 MATMUL and R8 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(lda,n)
  real ( kind = 8 ) c(lda,n)
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r8_set ( a, b, c, n, n )
 
    call matmul_cpu_timer ( time1 )
 
    a(1:n,1:n) = matmul ( b(1:n,1:n), c(1:n,1:n) )
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r8_ijk ( lda, n, nrep, ttime )

!*****************************************************************************80
!
!! R8_IJK multiplies A = B*C using index order IJK and R8 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  double precision a(lda,n)
  double precision b(lda,n)
  double precision c(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrep
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) ttime

  ttime = 0.0D+00
 
  do irep = 1, nrep
 
    call r8_set ( a, b, c, lda, n )
 
    call matmul_cpu_timer ( time1 )
 
    do i = 1, n
      do j = 1, n
        do k = 1, n
          a(i,k) = a(i,k) + b(i,j) * c(j,k)
        end do
      end do
    end do
 
    call matmul_cpu_timer ( time2 )
    ttime = ttime + time2 - time1

  end do
 
  return
end
subroutine r8_set ( a, b, c, lda, n )

!*****************************************************************************80
!
!! R8_SET initializes the matrices A, B and C using R8 arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
!    used in the multiplication.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  double precision a(lda,n)
  double precision b(lda,n)
  double precision c(lda,n)

  a(1:n,1:n) = 0.0D+00
  b(1:n,1:n) = 1.0D+00
  c(1:n,1:n) = 1.0D+00
 
  return
end
subroutine report ( c_show, mflops_show, lda, lda_show, lingo, lingo_show, &
  machine, machine_show, n, noshow, nrep, nrep_show, n_show, order, &
  output, time_show, ttime )

!*****************************************************************************80
!
!! REPORT reports the results for each multiplication experiment.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical C_SHOW, is TRUE if the multiplication method
!    is to be shown.
!
!    Input, logical MFLOPS_SHOW, is TRUE if the MegaFLOP rate is to be shown.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension used for arrays.
!
!    Input, logical LDA_SHOW, is TRUE if the variable LDA is to be shown.
!
!    Input, character ( len = 7 ) LINGO, the language in which MATMUL 
!    is written.
!
!    Input, logical LINGO_SHOW, is TRUE if the software's programming
!    language, stored in LINGO, is to be shown.
!
!    Input, character ( len = 10 ) MACHINE, the computer where MATMUL 
!    is running.
!
!    Input, logical MACHINE_SHOW, is TRUE if the machine name is to be shown.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the matrices.
!
!    Input, logical NOSHOW, is TRUE if the number of operations is
!    to be shown.
!
!    Input, integer ( kind = 4 ) NREP, the number of repetitions to carry out.
!
!    Input, logical NREP_SHOW, is TRUE if the number of repetitions is
!    to be shown.
!
!    Input, logical N_SHOW, is TRUE if the matrix size N is to be shown.
!
!    Input, character ( len = 10 ) ORDER, specifies the method to be used.
!
!    Output, character ( len = 82 ) OUTPUT, a string containing the data for
!    this multiplication.
!
!    Input, logical TIME_SHOW, is TRUE if the execution time is to be shown.
!
!    Output, real ( kind = 8 ) TTIME, an estimate of the CPU time in seconds required
!    for the matrix multiplications.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  logical c_show
  logical mflops_show
  real ( kind = 4 ) ftemp
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  logical lda_show
  character ( len = 7 ) lingo
  logical lingo_show
  character ( len = 10 ) machine
  logical machine_show
  logical noshow
  integer ( kind = 4 ) nrep
  logical nrep_show
  logical n_show
  integer ( kind = 4 ) ntemp
  character ( len = 10 ) order
  character ( len = 82 ) output
  logical time_show
  real ( kind = 8 ) ttime

  output = ' '
  ihi = 0

  if ( c_show ) then
    ilo = ihi + 1
    ihi = ilo + 9
    output(ilo:ihi) = order
  end if

  if ( lda_show ) then
    ilo = ihi + 1
    ihi = ilo + 3
    write ( output(ilo:ihi), '(i4)' ) lda
  end if

  if ( n_show ) then
    ilo = ihi + 1
    ihi = ilo + 3
    write ( output(ilo:ihi), '(i4)' ) n
  end if

  if ( time_show ) then
    ilo = ihi + 1
    ihi = ilo + 13
    write ( output(ilo:ihi), '(f14.6)' ) ttime
  end if

  if ( noshow ) then
    ntemp = 2 * n * n * n
    ilo = ihi + 1
    ihi = ilo + 9
    write ( output(ilo:ihi), '(i10)' ) ntemp
  end if

  if ( nrep_show ) then
    ilo = ihi + 1
    ihi = ilo + 4
    write ( output(ilo:ihi), '(i5)' ) nrep
  end if

  if ( mflops_show ) then

    if ( ttime == 0.0E+00 ) then
      ftemp = 0.0E+00
    else
      ftemp = real ( ntemp * nrep ) / ( 1.0E+06 * ttime )
    end if

    ilo = ihi + 1
    ihi = ilo + 9
    write ( output(ilo:ihi), '(f10.4)' ) ftemp

  end if

  if ( machine_show ) then
    ilo = ihi + 1
    ilo = ilo + 1
    ihi = ilo + 9
    output(ilo:ihi) = machine
  end if

  if ( lingo_show ) then
    ilo = ihi + 1
    ilo = ilo + 1
    ihi = ilo + 6
    output(ilo:ihi) = lingo
  end if

  if ( 0 < ihi ) then
    write ( *, '(a)' ) output(1:ihi)
  end if
 
  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  iput = 0

  do iget = 1, len ( s )

    c = s(iget:iget)

    if ( c /= ' ' .and. c /= TAB ) then
      iput = iput + 1
      s(iput:iput) = c
    end if

  end do

  s(iput+1:) = ' '

  return
end
subroutine s_cap ( string )

!*****************************************************************************80
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) STRING, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) string

  nchar = len ( string )

  do i = 1, nchar

    c = string(i:i)
    call ch_cap ( c )
    string(i:i) = c

  end do

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine taxpy ( n, sa, sx, incx, sy, incy )

!*****************************************************************************80
!
!! TAXPY is unoptimized standard BLAS routine SAXPY.
!
!  Discussion:
!
!    Roughly, TAXPY adds SA * SX(I) to SY(I) for I = 1 to N.
!    However, the increments INCX and INCY allow this to be done
!    even when SX or SY is a row or column of a matrix.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the "logical" number of items in the vectors.
!
!    Input, real ( kind = 4 ) SA, the multiplier.
!
!    Input, real ( kind = 4 ) SX(*), a vector, a multiple of which is to be added to SY.
!
!    Input, integer ( kind = 4 ) INCX, the increment in SX between the successive
!    elements that we will use.
!
!    Input/output, real ( kind = 4 ) SY(*), a vector, to which is to be added SA 
!    times an entry of SX.
!
!    Input, integer ( kind = 4 ) INCY, the increment in SY between the successive
!    elements that we will use.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ) sa
  real ( kind = 4 ) sx(*)
  real ( kind = 4 ) sy(*)

  if ( n <= 0 ) then
    return
  end if
 
  if ( sa == 0.0E+00 ) then
    return
  end if
 
  if ( incx /= 1 .or. incy /= 1 ) then
 
    if ( incx < 0 ) then
      ix = ( - n + 1 ) * incx + 1
    else
      ix = 1
    end if
 
    if ( incy < 0 ) then
      iy = ( - n + 1 ) * incy + 1
    else
      iy = 1
    end if
 
    do i = 1, n
      sy(iy) = sy(iy) + sa * sx(ix)
      ix = ix + incx
      iy = iy + incy
    end do
 
  else
 
    m = mod ( n, 4 )
 
    do i = 1, m
      sy(i) = sy(i) + sa * sx(i)
    end do
 
    do i = m+1, n, 4
      sy(i) =   sy(i)   + sa * sx(i)
      sy(i+1) = sy(i+1) + sa * sx(i+1)
      sy(i+2) = sy(i+2) + sa * sx(i+2)
      sy(i+3) = sy(i+3) + sa * sx(i+3)
    end do
 
  end if
 
  return
end
function tdot ( n, sx, incx, sy, incy )

!*****************************************************************************80
!
!! TDOT computes the inner product of two vectors.
!
!  Discussion:
!
!    TDOT is a source code version of the BLAS level 1 routine SDOT,
!    which can be used to compare performance with optimized versions
!    of SDOT supplied by the compiler or operating system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the "logical" number of items in the vectors.
!
!    Input, real ( kind = 4 ) SX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment in SX between the successive
!    elements that we will use.
!
!    Input/output, real ( kind = 4 ) SY(*), the second vector.
!
!    Input, integer ( kind = 4 ) INCY, the increment in SY between the successive
!    elements that we will use.
!
!    Output, real ( kind = 4 ) TDOT, the sum of the products of the appropriate
!    entries of SX and SY.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ) stemp
  real ( kind = 4 ) sx(*)
  real ( kind = 4 ) sy(*)
  real ( kind = 4 ) tdot

  tdot = 0.0E+00
 
  if ( n <= 0 ) then
    return
  end if
 
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( incx < 0 ) then
      ix = ( - n + 1 ) * incx + 1
    else
      ix = 1
    end if
 
    if ( incy < 0 ) then
      iy = ( - n + 1 ) * incy + 1
    else
      iy = 1
    end if
 
    stemp = 0.0E+00
    do i = 1, n
      stemp = stemp + sx(ix) * sy(iy)
      ix = ix + incx
      iy = iy + incy
    end do
 
    tdot = stemp
 
  else
 
    m = mod ( n, 5)
 
    stemp = 0.0E+00
 
    do i = 1, m
      stemp = stemp + sx(i) * sy(i)
    end do
 
    do i = m+1, n, 5
      stemp = stemp + sx(i)   * sy(i)   + sx(i+1) * sy(i+1) &
                    + sx(i+2) * sy(i+2) + sx(i+3) * sy(i+3) &
                    + sx(i+4) * sy(i+4)
    end do
 
  end if
 
  tdot = stemp
 
  return
end
subroutine tgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, &
  ldc )

!*****************************************************************************80
!
!! TGEMM is a source code copy of SGEMM, a BLAS matrix * matrix routine.
!
!  Discussion:
!
!    TGEMM performs one of the matrix-matrix operations
!
!       C : =  ALPHA * op( A ) * op( B ) + BETA * C,
!
!    where  op( X ) is one of
!
!       op( X ) = X   or   op( X ) = X',
!
!    ALPHA and BETA are scalars, and A, B and C are matrices, with op( A )
!    an M by K matrix, op( B ) a K by N matrix and C an M by N matrix.
!
!  Parameters:
!
!  transa - character.
!           On entry, transa specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              transa = 'N' or 'n',  op( A ) = A.
!
!              transa = 'T' or 't',  op( A ) = A'.
!
!              transa = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  transb - character.
!           On entry, transb specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              transa = 'N' or 'n',  op( B ) = B.
!
!              transa = 'T' or 't',  op( B ) = B'.
!
!              transa = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  m      - integer.
!           On entry,  m  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  m  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - integer.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - integer.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  alpha  - real .
!           On entry, alpha specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real array of DIMENSION ( lda, ka ), where ka is
!           k  when  transa = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  transa = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix a,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  lda    - integer.
!           On entry, lda specifies the first dimension of A as declared
!           in the calling (sub) program. When  transa = 'N' or 'n' then
!           lda must be at least  max( 1, m ), otherwise  lda must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - real array of DIMENSION ( ldb, kb ), where kb is
!           n  when  transb = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  transb = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  ldb    - integer.
!           On entry, ldb specifies the first dimension of B as declared
!           in the calling (sub) program. When  transb = 'N' or 'n' then
!           ldb must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  beta   - real .
!           On entry,  beta  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  c      - real array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - integer.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldb
  integer ( kind = 4 ) ldc

  real ( kind = 4 ) a(lda,*)
  real ( kind = 4 ) alpha
  real ( kind = 4 ) b(ldb,*)
  real ( kind = 4 ) beta
  real ( kind = 4 ) c(ldc,*)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncola
  logical nota
  logical notb
  integer ( kind = 4 ) nrowa
  integer ( kind = 4 ) nrowb
  logical tlsame
  character transa
  character transb
!
!  Set nota and notb as true if A and B respectively are not
!  transposed, and set nrowa, ncola and nrowb as the number of rows
!  and columns of A and the number of rows of B respectively.
!
  nota = tlsame ( transa, 'N' )
  notb = tlsame ( transb, 'N' )

  if ( nota  ) then
     nrowa = m
     ncola = k
  else
     nrowa = k
     ncola = m
  end if

  if ( notb  ) then
     nrowb = k
  else
     nrowb = n
  end if
!
!  Test the input parameters.
!
  info = 0
  if ( ( .not. nota                 ) .and. &
       ( .not. tlsame( transa, 'T' ) ).and. &
       ( .not. tlsame( transa, 'C' ) )       ) then
     info = 1
  else if ( ( .not. notb                 ) .and. &
    ( .not. tlsame( transb, 'T' ) ).and. &
    ( .not. tlsame( transb, 'C' ) )       ) then
     info = 2
  else if ( m < 0  ) then
     info = 3
  else if ( n < 0  ) then
     info = 4
  else if ( k < 0  ) then
     info = 5
  else if ( lda < max ( 1,nrowa )  ) then
     info = 8
  else if ( ldb < max ( 1,nrowb )  ) then
     info = 10
  else if ( ldc < max ( 1, m )  ) then
     info = 13
  end if
 
  if ( info /= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TGEMM - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value for input argument #', info
    stop
  end if
!
!  Quick return if possible.
!
  if ( m == 0 .or.  n==0 .or. &
    ( ( alpha == 0.0E+00 .or. k == 0 ) .and. beta == 1.0E+00 ) ) then
    return
  end if
!
!  Start the operations.
!
  if ( k == 0  ) then
!
!  Form  C : =  beta*C.
!
     if ( beta == 0.0E+00  ) then
       c(1:m,1:n) = 0.0E+00
     else
       c(1:m,1:n) = beta * c(1:m,1:n)
     end if

  else if ( notb  ) then
!
!  Form  C : =  alpha*op( A )*B + beta*C.
!
    do j = 1, n
      call tgemvf ( transa, nrowa, ncola, alpha, a, lda, b(1,j), 1, beta, &
        c(1,j), 1 )
    end do
 
  else
!
!  Form  C : =  alpha*op( A )*B' + beta*C.
!
    do j = 1, n
      call tgemvf ( transa, nrowa, ncola, alpha, a, lda, b(j,1), ldb, &
        beta, c(1,j), 1 )
    end do
 
  end if
 
  return
end
subroutine tgemvf ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

!*****************************************************************************80
!
!! TGEMVF is a source code copy of BLAS SGEMVF, a matrix * vector routine.
!
!  Discussion:
!
!    TGEMVF performs one of the matrix-vector operations
!
!      Y := alpha * A  * X + beta * Y,   
!
!    or
!   
!      Y := alpha * A' * X + beta * Y,
!
!    where ALPHA and BETA are scalars, X and Y are vectors and A is an
!    M by N matrix.
!
!  Parameters:
!
!  trans  - character.
!           On entry, trans specifies the operation to be performed as
!           follows:
!
!              trans = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              trans = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              trans = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  m      - integer.
!           On entry, m specifies the number of rows of the matrix A.
!           m must be at least zero.
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  alpha  - real .
!           On entry, alpha specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real array of DIMENSION ( lda, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  lda    - integer.
!           On entry, lda specifies the first dimension of A as declared
!           in the calling (sub) program. lda must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - real array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
!  beta   - real .
!           On entry, beta specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - real array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
!           Before entry with beta non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  incy   - integer.
!           On entry, incy specifies the increment for the elements of
!           Y. incy must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer ( kind = 4 ) lda

  real ( kind = 4 ) a(lda,*)
  real ( kind = 4 ) alpha
  real ( kind = 4 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jx
  integer ( kind = 4 ) jy
  integer ( kind = 4 ) kx
  integer ( kind = 4 ) ky
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) leny
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ) temp
  logical tlsame
  character trans
  real ( kind = 4 ) x(*)
  real ( kind = 4 ) y(*)
!
!  Test the input parameters.
!
  info = 0

  if ( .not. tlsame ( trans, 'N' ) .and. &
    .not. tlsame ( trans, 'T' ) .and. &
    .not. tlsame( trans, 'C' )       ) then
    info = 1
  else if ( m < 0  ) then
    info = 2
  else if ( n < 0  ) then
    info = 3
  else if ( lda < max ( 1, m )  ) then
    info = 6
  else if ( incx == 0  ) then
    info = 8
  else if ( incy == 0  ) then
    info = 11
  end if

  if ( info /= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TGEMVF - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value for input argument #', info
    stop
  end if
!
!  Quick return if possible.
!
  if ( ( m == 0 ) .or. ( n == 0 ).or. &
     ( ( alpha == 0.0E+00 ) .and. ( beta==1.0E+00 ) ) ) then
    return
  end if
!
!  Set  lenx  and  leny, the lengths of the vectors x and y, and set
!  up the start points in  X  and  Y.
!
  if ( tlsame ( trans, 'N' )  ) then
    lenx = n
    leny = m
  else
    lenx = m
    leny = n
  end if

  if ( 0 < incx  ) then
    kx = 1
  else
    kx = 1 - ( lenx - 1 )*incx
  end if

  if ( 0 < incy ) then
    ky = 1
  else
    ky = 1 - ( leny - 1 ) * incy
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
!  First form  y : =  beta*y.
!
  if ( beta /= 1.0E+00  ) then
     if ( incy == 1  ) then
        if ( beta == 0.0E+00  ) then
          y(1:leny) = 0.0E+00
        else
          y(1:leny) = beta * y(1:leny)
        end if
     else
        iy = ky
        if ( beta==0.0E+00 ) then
          do i = 1, leny
            y(iy) = 0.0E+00
            iy      = iy   + incy
          end do
        else
          do i = 1, leny
            y(iy) = beta * y(iy)
            iy      = iy           + incy
           end do
        end if
     end if
  end if

  if ( alpha==0.0E+00 ) then
    return
  end if

  if ( tlsame( trans, 'N' )  ) then
!
!  Form  y : =  alpha*A*x + y.
!
     jx = kx
     if ( incy == 1 ) then
        do j = 1, n
           if ( x( jx ) /= 0.0E+00  ) then
              temp = alpha * x(jx)
              do i = 1, m
                 y(i) = y(i) + temp * a(i,j)
              end do
           end if
           jx = jx + incx
        end do
     else
        do j = 1, n
           if ( x( jx ) /= 0.0E+00 ) then
              temp = alpha*x( jx )
              iy   = ky
              do i = 1, m
                 y(iy) = y(iy) + temp * a(i,j)
                 iy      = iy      + incy
             end do
           end if
           jx = jx + incx
        end do
     end if
  else
!
!  Form  y : =  alpha*A'*x + y.
!
     jy = ky
     if ( incx == 1  ) then
       do j = 1, n
           temp = 0.0E+00
           do i = 1, m
              temp = temp + a(i,j) * x(i)
           end do
           y( jy ) = y( jy ) + alpha * temp
           jy      = jy      + incy
       end do
     else
        do j = 1, n
           temp = 0.0E+00
           ix   = kx
           do i = 1,m
              temp = temp + a(i,j) * x(ix)
              ix   = ix   + incx
           end do
           y(jy) = y(jy) + alpha * temp
           jy      = jy      + incy
      end do
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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2001
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
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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
function tlsame ( ca, cb )

!*****************************************************************************80
!
!! TLSAME is a source code copy of BLAS LSAME, testing character equality.
!
!  Discussion:
!
!    CB is assumed to be an upper case letter.  The routine returns TRUE if
!    CA is either the same as CB or the equivalent lower case letter.
!
!    This version of the routine is only correct for ASCII code.
!    Installers must modify the routine for other character-codes.
!
!    For EBCDIC systems the constant ioff must be changed to -64.
!
!  Parameters:
!
!    Input, character CA, CB, the characters to be compared.
!
!    Output, logical TLSAME, is TRUE if CA is the same as CB, or
!    if CA is the lower-case version of CB.
!
  implicit none

  character ca
  character cb
  integer ( kind = 4 ), parameter :: ioff = 32
  logical tlsame
!
!  Test if the characters are equal.
!
  tlsame = ( ca == cb )
!
!  Test for equivalence.
!
  if ( .not. tlsame ) then
    tlsame = ichar ( ca ) - ioff == ichar ( cb )
  end if
 
  return
end
