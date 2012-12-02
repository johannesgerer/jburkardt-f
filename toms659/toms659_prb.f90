program main

!*******************************************************************************
!
!! MAIN is the main program for TOMS659_PRB.
!
!       THIS PROGRAM TESTS ACCURACY OF
!       NUMERICAL INTEGRATION USING "GOSOBL"
!       AND INTEGRAND (2) OF DAVIS AND
!       RABINOWITZ, PAGE 406
!
      integer, parameter :: maxdim = 1111

      integer atmost
      integer dimen
      double precision f
      logical flag(2)
      integer i
      integer j
      double precision quasi(maxdim)
      double precision sum
      integer taus

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS659_PRB'
      write ( *, '(a)' ) '  FORTRAN90 version'
      write ( *, '(a)' ) '  Test the TOMS659 library.'

   10 continue

      read (*,*) dimen,atmost
      if (dimen.eq.0) stop ' Run ends normally'
      write ( *,fmt='(1h1)')
      write ( *,fmt=*) 'Test SOBOL'
      write ( *,fmt=*) 'DIMENSION = ',DIMEN
      write ( *,fmt=*) 'ATMOST = ',ATMOST

      call cpu_time ( t1 )

      call insobl ( dimen, atmost )

      write ( *,*) 'Start time = ',T1
      write ( *,fmt=*) 'I = Iteration number'
      write ( *,fmt=*) 'EI = Estimate of integral'
      write ( *,fmt='(1h )')

      sum = 0.0
      do i = 1, atmost

          call i4_sobol ( dimen, quasi )

          f = 1.0
          do j = 1, dimen
              f = f * abs ( 4.0 * quasi(j) - 2.0 )
          end do
          sum = sum + f

          if ( mod ( i, 5000 ) .eq. 0 ) then
            write ( *,*) 'i = ', i
            write ( *,*) 'ei = ',sum / i
            call cpu_time ( t2 )
            write ( *,*) 'TIMEI = ', t2 - t1
            write ( *,'(1h )')
          end if

      end do

      write ( *,fmt=*) ' EI = ',sum / atmost

      go to 10

      end
