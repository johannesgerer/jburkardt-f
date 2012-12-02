program main

!*****************************************************************************80
!
!! MAIN is the main program for STARPAC_PRB.
!
!  Discussion:
!
!    STARPAC_PRB calls the STARPAC test problems.
!
!  Modified:
!
!    03 December 2006
!
  implicit none

  integer ldstak
  parameter ( ldstak = 3000 )

  double precision dstak(ldstak)
  integer ierr
  integer iflag
  integer mbo
  integer mbol
  integer mspect
  integer nfact
  integer nparar
  integer npardf
  integer nparma
  integer nrests
  integer parar
  integer pardf
  integer parma
  real q
  integer t
  integer temp

  save / cstak /
  save / errchk /
  save / mdltsc /
  save / notopt /

  common / cstak / dstak
  common / errchk / ierr
  common / mdltsc / mspect,nfact,pardf,npardf,parar,nparar,parma, &
     nparma,mbo,mbol,t,temp,nrests,iflag
  common / notopt / q

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STARPAC_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tests for the STARPAC statistical package.'

  call xacf ( ldstak )
  call xaimd ( ldstak )
  call xaimt ( ldstak )
  call xaov1 ( ldstak )
  call xbfs ( ldstak )
  call xccf ( ldstak )
  call xcorr ( ldstak )
  call xdckld ( ldstak )
  call xdckle ( ldstak )
  call xdcklt ( ldstak )
  call xdemod ( ldstak )
  call xdflt ( ldstak )
  call xhist ( ldstak )
  call xlls ( ldstak )
  call xnlsd ( ldstak )
  call xnlse ( ldstak )
  call xnlst ( ldstak )
  call xnrand ( ldstak )
  call xpgm ( ldstak )
  call xpp ( )
  call xstat ( ldstak )
  call xstpld ( ldstak )
  call xstple ( ldstak )
  call xstplt ( ldstak )
  call xuas ( ldstak )
  call xufs ( ldstak )
  call xvp ( )
  call xxch1 ( ldstak )
  call xxch2 ( )
  call xxch3 ( )
  call xxch4 ( ldstak )
  call xxch5 ( ldstak )
  call xxch6 ( ldstak )
  call xxch7 ( ldstak )
  call xxch8 ( ldstak )
  call xxch9 ( ldstak )
  call xxch10 ( )
  call xxch11 ( ldstak )
  call xxch12 ( ldstak )
  call xxch13 ( ldstak )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STARPAC_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine xacf ( lds )

!*****************************************************************************80
!
!! XACF tests the time series correlation subroutines.
!
!  Discussion:
!
!    Series y is listed as series x1 on page 362 in jenkins and watts.
!
!    Series yd is listed as series g on page 531 of box and jenkins.
!
!  Modified:
!
!    01 December 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    George Box, Gwilym Jenkins,
!    Time Series Analysis: Forecasting and Control,
!    Holden-Day, 1970,
!    QA280.B67.
!
!    Gwilym Jenkins, Donald Watts,
!    Spectral Analysis and its Applications,
!    Holden-Day 1968.
!
!  Parameters:
!
!     real acov(21)
!        the autocovariance vector.
!     real amiss
!        the missing value code for the returned acvf estimates
!        (vector acov).
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an indexing variable.
!     integer iar
!        the order of the autoregressive process chosen.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list
!        if ierr == 0, no errors were detected
!     integer iod(2)
!        the order of each of the difference factors.
!     integer itest
!        the number of the test being run
!     integer lacov
!        the length of the acvf related vectors.
!     integer lagmax
!        the maximum lag value requested.
!     integer ldstak
!        the length of the array dstak.
!     integer lyfft
!        the length of the arrays used when the computations are
!        performed by the fft.
!     integer n
!        the integer number of observations in each series
!     integer nd(2)
!        the array containing the number of times the difference
!        factors are to be applied.
!     integer nfac
!        the number of difference factors.
!     integer nlppa(21)
!        the array containing the number of lagged product pairs
!        used to compute each acvf estimate.
!     integer nprt
!        the indicator variable used to specify whether or not
!        printed output is to be given, where if the value of
!        nprt is zero, no output is made.
!     integer nyd
!        the number of observations in the series to be differenced.
!     real phi(21)
!        the array of autoregressive coefficients for the selected
!        order.
!     real y(100), yd(150)
!        the vector containing the observed time series
!     real yfft(150)
!        the vectors used for storing the series for the routines
!        using the fft.
!     real ymiss
!        the missing value codes for series y and ym.
!
  implicit none

  real acov(21)
  real amiss
  double precision dstak(12)
  integer i
  integer iar
  integer ierr
  integer iod(2)
  integer itest
  integer lacov
  integer lagmax
  integer lds
  integer ldstak
  integer lyfft
  integer n
  integer nd(2)
  integer nfac
  integer nlppa(21)
  integer nprt
  integer nyd
  real phi(21)
  real y(100)
  real yd(150)
  real yfft(150)
  real ymiss

  external acf,acfd,acff,acffs,acfm,acfms,acfs,iprint

  common /cstak/ dstak
  common /errchk/ ierr

  data    y(  1),   y(  2),   y(  3),   y(  4),   y(  5),   y(  6) &
      / -2.07e0, -1.15e0,  0.69e0, -0.46e0, -1.49e0, -0.70e0/
  data    y(  7),   y(  8),   y(  9),   y( 10),   y( 11),   y( 12) &
      / -1.07e0, -0.69e0, -0.68e0,  1.27e0, -1.05e0, -0.05e0/
  data    y( 13),   y( 14),   y( 15),   y( 16),   y( 17),   y( 18) &
      / -0.84e0, -0.62e0, -0.49e0, -1.29e0, -0.49e0, -1.06e0/
  data    y( 19),   y( 20),   y( 21),   y( 22),   y( 23),   y( 24) &
      / -0.38e0, -0.52e0, -0.13e0,  1.30e0, -1.51e0, -0.43e0/
  data    y( 25),   y( 26),   y( 27),   y( 28),   y( 29),   y( 30) &
      / -1.33e0, -0.78e0,  0.31e0, -0.95e0, -0.90e0, -0.30e0/
  data    y( 31),   y( 32),   y( 33),   y( 34),   y( 35),   y( 36) &
      / -1.02e0, -0.53e0,  0.15e0,  1.40e0,  1.22e0,  0.59e0/
  data    y( 37),   y( 38),   y( 39),   y( 40),   y( 41),   y( 42) &
      /  0.70e0,  1.70e0,  2.78e0,  1.98e0,  1.39e0,  1.85e0/
  data    y( 43),   y( 44),   y( 45),   y( 46),   y( 47),   y( 48) &
      /  2.60e0,  0.51e0,  2.77e0,  1.16e0,  1.07e0, -0.48e0/
  data    y( 49),   y( 50),   y( 51),   y( 52),   y( 53),   y( 54) &
      / -0.52e0,  0.37e0,  0.00e0, -1.99e0, -1.75e0,  0.70e0/
  data    y( 55),   y( 56),   y( 57),   y( 58),   y( 59),   y( 60) &
      /  0.73e0,  1.16e0,  0.06e0, -0.02e0,  1.10e0, -0.35e0/
  data    y( 61),   y( 62),   y( 63),   y( 64),   y( 65),   y( 66) &
      / -1.67e0, -1.57e0,  1.16e0,  1.84e0,  3.35e0,  0.40e0/
  data    y( 67),   y( 68),   y( 69),   y( 70),   y( 71),   y( 72) &
      /  0.45e0,  1.30e0,  0.93e0,  1.17e0, -1.74e0, -1.28e0/
  data    y( 73),   y( 74),   y( 75),   y( 76),   y( 77),   y( 78) &
      / -0.07e0,  1.50e0,  0.53e0,  0.20e0, -0.42e0,  1.18e0/
  data    y( 79),   y( 80),   y( 81),   y( 82),   y( 83),   y( 84) &
      /  0.82e0,  1.50e0,  2.92e0,  1.18e0,  1.23e0,  3.16e0/
  data    y( 85),   y( 86),   y( 87),   y( 88),   y( 89),   y( 90) &
      /  0.79e0,  0.68e0,  1.14e0,  1.02e0,  1.02e0, -0.71e0/
  data    y( 91),   y( 92),   y( 93),   y( 94),   y( 95),   y( 96) &
      / -0.17e0, -1.50e0, -0.26e0, -0.38e0,  0.93e0, -0.33e0/
  data    y( 97),   y( 98),   y( 99),   y(100) &
      / -1.12e0, -2.95e0, -2.09e0, -1.11e0                    /
!
  data   yd(  1),  yd(  2),  yd(  3),  yd(  4),  yd(  5),  yd(  6) &
      /  112.0e0, 118.0e0, 132.0e0, 129.0e0, 121.0e0, 135.0e0/
  data   yd(  7),  yd(  8),  yd(  9),  yd( 10),  yd( 11),  yd( 12) &
      /  148.0e0, 148.0e0, 136.0e0, 119.0e0, 104.0e0, 118.0e0/
  data   yd( 13),  yd( 14),  yd( 15),  yd( 16),  yd( 17),  yd( 18) &
      /  115.0e0, 126.0e0, 141.0e0, 135.0e0, 125.0e0, 149.0e0/
  data   yd( 19),  yd( 20),  yd( 21),  yd( 22),  yd( 23),  yd( 24) &
      /  170.0e0, 170.0e0, 158.0e0, 133.0e0, 114.0e0, 140.0e0/
  data   yd( 25),  yd( 26),  yd( 27),  yd( 28),  yd( 29),  yd( 30) &
      /  145.0e0, 150.0e0, 178.0e0, 163.0e0, 172.0e0, 178.0e0/
  data   yd( 31),  yd( 32),  yd( 33),  yd( 34),  yd( 35),  yd( 36) &
      /  199.0e0, 199.0e0, 184.0e0, 162.0e0, 146.0e0, 166.0e0/
  data   yd( 37),  yd( 38),  yd( 39),  yd( 40),  yd( 41),  yd( 42) &
      /  171.0e0, 180.0e0, 193.0e0, 181.0e0, 183.0e0, 218.0e0/
  data   yd( 43),  yd( 44),  yd( 45),  yd( 46),  yd( 47),  yd( 48) &
      /  230.0e0, 242.0e0, 209.0e0, 191.0e0, 172.0e0, 194.0e0/
  data   yd( 49),  yd( 50),  yd( 51),  yd( 52),  yd( 53),  yd( 54) &
      /  196.0e0, 196.0e0, 236.0e0, 235.0e0, 229.0e0, 243.0e0/
  data   yd( 55),  yd( 56),  yd( 57),  yd( 58),  yd( 59),  yd( 60) &
      /  264.0e0, 272.0e0, 237.0e0, 211.0e0, 180.0e0, 201.0e0/
  data   yd( 61),  yd( 62),  yd( 63),  yd( 64),  yd( 65),  yd( 66) &
      /  204.0e0, 188.0e0, 235.0e0, 227.0e0, 234.0e0, 264.0e0/
  data   yd( 67),  yd( 68),  yd( 69),  yd( 70),  yd( 71),  yd( 72) &
      /  302.0e0, 293.0e0, 259.0e0, 229.0e0, 203.0e0, 229.0e0/
  data   yd( 73),  yd( 74),  yd( 75),  yd( 76),  yd( 77),  yd( 78) &
      /  242.0e0, 233.0e0, 267.0e0, 269.0e0, 270.0e0, 315.0e0/
  data   yd( 79),  yd( 80),  yd( 81),  yd( 82),  yd( 83),  yd( 84) &
      /  364.0e0, 347.0e0, 312.0e0, 274.0e0, 237.0e0, 278.0e0/
  data   yd( 85),  yd( 86),  yd( 87),  yd( 88),  yd( 89),  yd( 90) &
      /  284.0e0, 277.0e0, 317.0e0, 313.0e0, 318.0e0, 374.0e0/
  data   yd( 91),  yd( 92),  yd( 93),  yd( 94),  yd( 95),  yd( 96) &
      /  413.0e0, 405.0e0, 355.0e0, 306.0e0, 271.0e0, 306.0e0/
  data   yd( 97),  yd( 98),  yd( 99),  yd(100),  yd(101),  yd(102) &
      /  315.0e0, 301.0e0, 356.0e0, 348.0e0, 355.0e0, 422.0e0/
  data   yd(103),  yd(104),  yd(105),  yd(106),  yd(107),  yd(108) &
      /  465.0e0, 467.0e0, 404.0e0, 347.0e0, 305.0e0, 336.0e0/
  data   yd(109),  yd(110),  yd(111),  yd(112),  yd(113),  yd(114) &
      /  340.0e0, 318.0e0, 362.0e0, 348.0e0, 363.0e0, 435.0e0/
  data   yd(115),  yd(116),  yd(117),  yd(118),  yd(119),  yd(120) &
      /  491.0e0, 505.0e0, 404.0e0, 359.0e0, 310.0e0, 337.0e0/
  data   yd(121),  yd(122),  yd(123),  yd(124),  yd(125),  yd(126) &
      /  360.0e0, 342.0e0, 406.0e0, 396.0e0, 420.0e0, 472.0e0/
  data   yd(127),  yd(128),  yd(129),  yd(130),  yd(131),  yd(132) &
      /  548.0e0, 559.0e0, 463.0e0, 407.0e0, 362.0e0, 405.0e0/
  data   yd(133),  yd(134),  yd(135),  yd(136),  yd(137),  yd(138) &
      /  417.0e0, 391.0e0, 419.0e0, 461.0e0, 472.0e0, 535.0e0/
  data   yd(139),  yd(140),  yd(141),  yd(142),  yd(143),  yd(144) &
      /  622.0e0, 606.0e0, 508.0e0, 461.0e0, 390.0e0, 432.0e0/

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XACF'
  write ( *, '(a)' ) '  Test the time series correlation routines.'
  write ( *, '(a)' ) ' '

  itest = 1
  ldstak = lds

  n = 100
  lagmax = 20
  nprt = 1
  lyfft = 150
  lacov = 21
  nyd = 144
  nfac = 2
  nd(1) = 1
  nd(2) = 1
  iod(1) = 12
  iod(2) = 1
  ymiss = 1.16e0
!
!  Test of acf
!
5 continue

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACF:'
  write ( *, '(a)' ) ' '

  call acf ( y, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  Test of acfs
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACFS:'
  write ( *, '(a)' ) ' '

  call acfs ( y, n, lagmax, lacov, acov, iar, phi, nprt, ldstak )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  print storage from acfs
!
  if (ierr==0) then
    write ( *,'(9f10.5)') (acov(i),i=1,lagmax+1)
    write ( *,'(9f10.5)') (phi(i),i=1,iar)
  end if
!
!  Test of acfd
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACFD:'
  write ( *, '(a)' ) ' '

  call acfd(yd, nyd, lagmax, nfac, nd, iod, ldstak)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  Test of acfm
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACFM:'
  write ( *, '(a)' ) ' '

  call acfm ( y, ymiss, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  Test of acfms
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACFMS:'
  write ( *, '(a)' ) ' '

  call acfms(y, ymiss, n, lagmax, lacov, acov, amiss, nlppa, nprt, &
     ldstak)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  print storage from acfms
!
  if (ierr==0) then
    write ( *,'(9f10.5)') (acov(i),i=1,lagmax+1)
    write ( *,'(9i10)') (nlppa(i),i=1,lagmax+1)
  end if
!
!  copy data into yfft for acff
!
  yfft(1:n) = y(1:n)
!
!  Test of acff
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACFF:'
  write ( *, '(a)' ) ' '

  call acff ( yfft, n, lyfft, ldstak )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  copy data into yfft for acffs
!
  yfft(1:n) = y(1:n)
!
!  Test of acffs
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACFFS:'
  write ( *, '(a)' ) ' '

  call acffs ( yfft, n, lyfft, ldstak, lagmax, lacov, acov, iar, phi, &
     nprt )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  print storage from acffs
!
  if ( ierr == 0 ) then
    write ( *,'(9f10.5)') (acov(i),i=1,lagmax+1)
    write ( *,'(9f10.5)') (phi(i),i=1,iar)
  end if
!
!  Test minimum problem size
!
  if ( itest == 1 ) then

    itest = itest + 1
    n = 13
    lagmax = 1
    nfac = 1
    nd(1) = 1
    iod(1) = 1
!
!  Check error handling
!
  else if ( itest == 2 ) then

    itest = itest + 1
    n = 0
    lagmax = 20
    lyfft = 0
    lacov = 0
    nyd = 0
    nfac = 1
    nd(1) = 0
    iod(1) = 0
    go to 5
!
!  Check error handling
!
  else if ( itest == 3 ) then

    itest = itest + 1
    n = 100
    lagmax = 0
    lyfft = 0
    lacov = 0
    nyd = 144
    nfac = 0
    ldstak = 0
    go to 5

  end if

  return
end
subroutine xaimd ( ldstak )

!*****************************************************************************80
!
!! XAIMD demonstrates the user callable routines in the ARIMA family.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     real delta
!        the maximum change allowed in the model parameters at the
!        first iteration.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fcst(200, 10)
!        the forecasts.
!     real fcstsd(200, 10)
!        the standard deviations of the forecasts.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr .ge. 1, errors were detected.
!     integer ifcst
!        the first dimension of the array fcst.
!     integer ifcsto(10)
!        the indices of the origins for the forecasts.
!     integer ifixed(10)
!        the indicator values used to designate whether the
!        parameters are to be optimized or are to be held fixed.  if
!        ifixed(i).ne.0, then par(i) will be optimized.  if
!        ifixed(i)==0, then par(i) will be held fixed.
!        ifixed(i).lt.0, then all par(i),i=1,npar, will be optimized..
!     integer ivaprx
!        an indicator value used to designate which option is to be used
!        to compute the variance covariance matrix (vcv), where
!        ivaprx le 0 indicates the the default option will be used
!        ivaprx eq 1 indicates the vcv is to be computed by
!                       inverse(hessian)*transpose(jacobian)*jacobian
!                          *inverse(hessian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 2 indicates the vcv is to be computed by
!                       inverse(hessian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 3 indicates the vcv is to be computed by
!                       inverse(transpose(jacobian)*jacobian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 4 indicates the vcv is to be computed by
!                       inverse(hessian)*transpose(jacobian)*jacobian
!                          *inverse(hessian)
!                    using only the model subroutine
!        ivaprx eq 5 indicates the vcv is to be computed by
!                       inverse(hessian)
!                    using only the model subroutine
!        ivaprx eq 6 indicates the vcv is to be computed by
!                       inverse(transpose(jacobian)*jacobian)
!                    using only the model subroutine
!        ivaprx ge 7 indicates the default option will be used
!     integer ivcv
!        the first dimension of the variance covariance matrix vcv.
!     integer ldstak
!        the length of the array dstak.
!     integer mit
!        the maximum number of iterations allowed.
!     integer mspec(4,10)
!        the array containing the values of p, d, q, and s for each
!        factor.
!     integer n
!        the number of observations.
!     integer nfac
!        the number of factors in the model
!     integer nfcst
!        the number of forecasts.
!     integer nfcsto
!        the number of the origins.
!     integer npar
!        the number of unknown parameters in the model.
!     integer npare
!        the number of parameters estimated by the routine.
!     integer nprt
!        the parameter used to indicate how much printed output is
!        to be provided.
!     integer ntest
!        the number of the current test.
!     real par(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real pv(200)
!        the predicted value based on the current parameter estimates
!     real res(200)
!        the residuals from the fit.
!     real rsd
!        the value of the residual standard deviation at the solution.
!     real scale(10)
!        a value to indicate use of the default values of
!        the typical size of the unknown parameters.
!     real sdpv(200)
!        the standard deviation of the predicted value.
!     real sdres(200)
!        the standard deviations of the residuals.
!     real stopp
!        the stopping criterion for the test based on the maximum scaled
!        relative change in the elements of the model parameter vector
!     real stopss
!        the stopping criterion for the test based on the ratio of the
!        predicted decrease in the residual sum of squares (computed
!        by starpac) to the current residual sum of squares estimate.
!     real stp(10)
!        the rcstep size array.
!     real vcv(10,10)
!        the covariance matrix.
!     real y(200)
!        the array of the dependent variable.
!
  implicit none

  real delta
  double precision dstak(12)
  integer ierr
  integer ldstak
  integer mspec(4,10)
  integer n
  integer nfac
  integer npar
  real par(10)
  real res(200)
  real rsd
  real stopp
  real stopss
  real y(200)

  integer &
     ifcst,ivaprx,ivcv,mit,mxfac,mxfc,mxfco,mxn,mxpar, &
     nfcst,nfcsto,npare,nprt,ntest
!
!  local arrays
  real &
     fcst(200,10),fcstsd(200,10),pv(200), &
     scale(10),sdpv(200),sdres(200),stp(10),vcv(10,10)
  integer &
     ifcsto(10),ifixed(10)
!
!  external subroutines
  external aime,aimec,aimes,aimf,aimfs,aimx1,fitxsp,iprint
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!  define constants
!
  data    y(  1),   y(  2),   y(  3),   y(  4),   y(  5),   y(  6) &
      / 112.0e0, 118.0e0, 132.0e0, 129.0e0, 121.0e0, 135.0e0/
  data    y(  7),   y(  8),   y(  9),   y( 10),   y( 11),   y( 12) &
      / 148.0e0, 148.0e0, 136.0e0, 119.0e0, 104.0e0, 118.0e0/
  data    y( 13),   y( 14),   y( 15),   y( 16),   y( 17),   y( 18) &
      / 115.0e0, 126.0e0, 141.0e0, 135.0e0, 125.0e0, 149.0e0/
  data    y( 19),   y( 20),   y( 21),   y( 22),   y( 23),   y( 24) &
      / 170.0e0, 170.0e0, 158.0e0, 133.0e0, 114.0e0, 140.0e0/
  data    y( 25),   y( 26),   y( 27),   y( 28),   y( 29),   y( 30) &
      / 145.0e0, 150.0e0, 178.0e0, 163.0e0, 172.0e0, 178.0e0/
  data    y( 31),   y( 32),   y( 33),   y( 34),   y( 35),   y( 36) &
      / 199.0e0, 199.0e0, 184.0e0, 162.0e0, 146.0e0, 166.0e0/
  data    y( 37),   y( 38),   y( 39),   y( 40),   y( 41),   y( 42) &
      / 171.0e0, 180.0e0, 193.0e0, 181.0e0, 183.0e0, 218.0e0/
  data    y( 43),   y( 44),   y( 45),   y( 46),   y( 47),   y( 48) &
      / 230.0e0, 242.0e0, 209.0e0, 191.0e0, 172.0e0, 194.0e0/
  data    y( 49),   y( 50),   y( 51),   y( 52),   y( 53),   y( 54) &
      / 196.0e0, 196.0e0, 236.0e0, 235.0e0, 229.0e0, 243.0e0/
  data    y( 55),   y( 56),   y( 57),   y( 58),   y( 59),   y( 60) &
      / 264.0e0, 272.0e0, 237.0e0, 211.0e0, 180.0e0, 201.0e0/
  data    y( 61),   y( 62),   y( 63),   y( 64),   y( 65),   y( 66) &
      / 204.0e0, 188.0e0, 235.0e0, 227.0e0, 234.0e0, 264.0e0/
  data    y( 67),   y( 68),   y( 69),   y( 70),   y( 71),   y( 72) &
      / 302.0e0, 293.0e0, 259.0e0, 229.0e0, 203.0e0, 229.0e0/
  data    y( 73),   y( 74),   y( 75),   y( 76),   y( 77),   y( 78) &
      / 242.0e0, 233.0e0, 267.0e0, 269.0e0, 270.0e0, 315.0e0/
  data    y( 79),   y( 80),   y( 81),   y( 82),   y( 83),   y( 84) &
      / 364.0e0, 347.0e0, 312.0e0, 274.0e0, 237.0e0, 278.0e0/
  data    y( 85),   y( 86),   y( 87),   y( 88),   y( 89),   y( 90) &
      / 284.0e0, 277.0e0, 317.0e0, 313.0e0, 318.0e0, 374.0e0/
  data    y( 91),   y( 92),   y( 93),   y( 94),   y( 95),   y( 96) &
      / 413.0e0, 405.0e0, 355.0e0, 306.0e0, 271.0e0, 306.0e0/
  data    y( 97),   y( 98),   y( 99),   y(100),   y(101),   y(102) &
      / 315.0e0, 301.0e0, 356.0e0, 348.0e0, 355.0e0, 422.0e0/
  data    y(103),   y(104),   y(105),   y(106),   y(107),   y(108) &
      / 465.0e0, 467.0e0, 404.0e0, 347.0e0, 305.0e0, 336.0e0/
  data    y(109),   y(110),   y(111),   y(112),   y(113),   y(114) &
      / 340.0e0, 318.0e0, 362.0e0, 348.0e0, 363.0e0, 435.0e0/
  data    y(115),   y(116),   y(117),   y(118),   y(119),   y(120) &
      / 491.0e0, 505.0e0, 404.0e0, 359.0e0, 310.0e0, 337.0e0/
  data    y(121),   y(122),   y(123),   y(124),   y(125),   y(126) &
      / 360.0e0, 342.0e0, 406.0e0, 396.0e0, 420.0e0, 472.0e0/
  data    y(127),   y(128),   y(129),   y(130),   y(131),   y(132) &
      / 548.0e0, 559.0e0, 463.0e0, 407.0e0, 362.0e0, 405.0e0/
  data    y(133),   y(134),   y(135),   y(136),   y(137),   y(138) &
      / 417.0e0, 391.0e0, 419.0e0, 461.0e0, 472.0e0, 535.0e0/
  data    y(139),   y(140),   y(141),   y(142),   y(143),   y(144) &
      / 622.0e0, 606.0e0, 508.0e0, 461.0e0, 390.0e0, 432.0e0/

  y(1:144) = log(y(1:144))
!
!  set dimensions
!
  mxn = 200
  mxpar = 10
  mxfc = 200
  mxfco = 10
  mxfac = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XAIMD'
  write ( *, '(a)' ) '  Demonstrate codes in the ARIMA family.'

  ntest = 0
!
!  Test on normal statement.
!
  ntest = ntest + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ARIMA test number ', ntest
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) '  Normal problem'
  write ( *, '(a)' ) '  Test of AIM'
!
!  Set the starting parameters for AIMX.
!
  call aimx1 ( mxn, mxpar, mxfc, mxfco, mxfac, &
     1, n, mspec, nfac, par, npar, res, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv, &
     nfcst, nfcsto, ifcsto, fcst, ifcst, fcstsd )

  call aime ( y, n, mspec, nfac, par, npar, res, ldstak )

  write ( *,1120) ierr

  call fitxsp ( par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, &
     n, npare, rsd )

  ntest = ntest + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ARIMA test number ', ntest
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal problem'
  write ( *, '(a)' ) '  Test of AIMC'
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt

  call aimx1 ( mxn, mxpar, mxfc, mxfco, mxfac, &
     1, n, mspec, nfac, par, npar, res, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv, &
     nfcst, nfcsto, ifcsto, fcst, ifcst, fcstsd )

  call aimec ( y, n, mspec, nfac, par, npar, res, ldstak, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt )

  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  call fitxsp ( par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, &
     n, npare, rsd )

  ntest = ntest + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ARIMA test number ', ntest
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) '  Normal problem'
  write ( *, '(a)' ) '  Test of AIMS'
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt

  call aimx1 ( mxn, mxpar, mxfc, mxfco, mxfac, &
     1, n, mspec, nfac, par, npar, res, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv, &
     nfcst, nfcsto, ifcsto, fcst, ifcst, fcstsd)

  call aimes ( y, n, mspec, nfac, par, npar, res, ldstak, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv )
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp ( par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, &
     n, npare, rsd )

  ntest = ntest + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ARIMA test number ', ntest
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) '  Normal problem'
  write ( *, '(a)' ) '  Test of AIMF'

  call aimx1 ( mxn, mxpar, mxfc, mxfco, mxfac, &
     1, n, mspec, nfac, par, npar, res, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv, &
     nfcst, nfcsto, ifcsto, fcst, ifcst, fcstsd )

  call aimf ( y, n, mspec, nfac, par, npar, ldstak )
  write ( *,1120) ierr

  ntest = ntest + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ARIMA test number ', ntest
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) '  Normal problem'
  write ( *, '(a)' ) '  Test of AIMFS'
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt

  call aimx1 ( mxn, mxpar, mxfc, mxfco, mxfac, &
     1, n, mspec, nfac, par, npar, res, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv, &
     nfcst, nfcsto, ifcsto, fcst, ifcst, fcstsd)

  call aimfs ( y, n, mspec, nfac, par, npar, ldstak, &
     nfcst, nfcsto, ifcsto, nprt, fcst, ifcst, fcstsd )
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  call fitxsp ( par, fcst(1,1), fcst(1,2), fcst(1,3), fcstsd, vcv, &
    n, npar, ivcv, n, npare, rsd )

  return

!1010 format ('test of aimc'  )
!1020 format ('test of aims'  )
!1030 format ('test of aimf' )
!1040 format ('test of aimfs' )

 1120 format (/' ***** returned results *****', 5x, '(-1 indicates ', &
     'value not changed by called subroutine)'//' ierr is ', i3)
 1340 format (//24h input   -  ifixed(1) = , i6, 9x, ', stp(1) = ', &
     g15.8, ',    mit = ',i5, ', stopss = ', g15.8, 10h, stopp = , &
     g15.8/13x, 'scale(1) = ', g15.8, ',  delta = ', g15.8, &
     ', ivaprx = ', i5, ',   nprt = ', i5//)
 1350 format (//24h output  -  ifixed(1) = , i6, 9x, ', stp(1) = ', &
     g15.8, ',    mit = ',i5, ', stopss = ', g15.8, 10h, stopp = , &
     g15.8/13x, 'scale(1) = ', g15.8, ',  delta = ', g15.8, &
     ', ivaprx = ', i5, ',   nprt = ', i5//)
end
subroutine xaimt ( ldstak )

!*****************************************************************************80
!
!! XAIMT tests the time series model estimation routines.
!
!  Discussion:
!
!    Series y is the airline data listed on page 531 of box and jenkins
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    George Box, Gwilym Jenkins,
!    Time Series Analysis: Forecasting and Control,
!    Holden-Day, 1970,
!    QA280.B67.
!
!  Parameters:
!
!     real delta
!        the maximum change allowed in the model parameters at the
!        first iteration.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fcst(50,5)
!        the forecasts.
!     real fcstsd(50)
!        the standard deviations of the forecasts.
!     integer i
!        *
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr .ge. 1, errors were detected.
!     integer ifcst
!        the first dimension of the array fcst.
!     integer ifixed(50)
!        the indicator values used to designate whether the
!        parameters are to be optimized or are to be held fixed.  if
!        ifixed(i).ne.0, then par(i) will be optimized.  if
!        ifixed(i)==0, then par(i) will be held fixed.
!        ifixed(i).lt.0, then all par(i),i=1,npar, will be optimized..
!     integer ivaprx
!        an indicator value used to designate which option is to be used
!        to compute the variance covariance matrix (vcv), where
!        ivaprx le 0 indicates the the default option will be used
!        ivaprx eq 1 indicates the vcv is to be computed by
!                       inverse(hessian)*transpose(jacobian)*jacobian
!                          *inverse(hessian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 2 indicates the vcv is to be computed by
!                       inverse(hessian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 3 indicates the vcv is to be computed by
!                       inverse(transpose(jacobian)*jacobian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 4 indicates the vcv is to be computed by
!                       inverse(hessian)*transpose(jacobian)*jacobian
!                          *inverse(hessian)
!                    using only the model subroutine
!        ivaprx eq 5 indicates the vcv is to be computed by
!                       inverse(hessian)
!                    using only the model subroutine
!        ivaprx eq 6 indicates the vcv is to be computed by
!                       inverse(transpose(jacobian)*jacobian)
!                    using only the model subroutine
!        ivaprx ge 7 indicates the default option will be used
!     integer ivcv
!        the first dimension of the variance covariance matrix vcv.
!     integer ldstak
!        the length of the array dstak.
!     integer mit
!        the maximum number of iterations allowed.
!     integer mspec(4,50)
!        the array containing the values of p, d, q, and s for each
!        factor.
!     integer nfac
!        the number of factors in the model
!     integer npar
!        the number of unknown parameters in the model.
!     integer npare
!        the number of parameters estimated by the routine.
!     integer nprt
!        the parameter used to indicate how much printed output is
!        to be provided.
!     integer ny
!        the number of observations.
!     real par(50)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real pv(200)
!        the predicted value based on the current parameter estimates
!     real res(200)
!        the residuals from the fit.
!     real rsd
!        the value of the residual standard deviation at the solution.
!     real scale(50)
!        a value to indicate use of the default values of
!        the typical size of the unknown parameters.
!     real sdpv(200)
!        the standard deviation of the predicted value.
!     real sdres(200)
!        the standard deviations of the residuals.
!     real stopp
!        the stopping criterion for the test based on the maximum scaled
!        relative change in the elements of the model parameter vector
!     real stopss
!        the stopping criterion for the test based on the ratio of the
!        predicted decrease in the residual sum of squares (computed
!        by starpac) to the current residual sum of squares estimate.
!     real stp(50)
!        the rcstep size array.
!     real vcv(10,10)
!        the covariance matrix.
!     real y(200),ylog(200),yt(200)
!        the array of the dependent variable.
!
  implicit none

  double precision dstak(12)
  integer ierr
  integer ldstak
!
!  local scalars
  real &
     delta,rsd,stopp,stopss
  integer ifcst,ivaprx,ivcv,mit,nfac,npar,npare,nprt,ny
!
!  local arrays
  real &
     fcst(50,5),fcstsd(50),par(50),pv(200),res(200),scale(50), &
     sdpv(200),sdres(200),stp(50),vcv(10,10),y(200),ylog(200), &
     yt(200)
  integer &
     ifixed(50),mspec(4,50)
!
!  external subroutines
  external aime,aimec,aimes,aimf,aimfs,iprint
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data    y(  1),   y(  2),   y(  3),   y(  4),   y(  5),   y(  6) &
      / 112.0e0, 118.0e0, 132.0e0, 129.0e0, 121.0e0, 135.0e0/
  data    y(  7),   y(  8),   y(  9),   y( 10),   y( 11),   y( 12) &
      / 148.0e0, 148.0e0, 136.0e0, 119.0e0, 104.0e0, 118.0e0/
  data    y( 13),   y( 14),   y( 15),   y( 16),   y( 17),   y( 18) &
      / 115.0e0, 126.0e0, 141.0e0, 135.0e0, 125.0e0, 149.0e0/
  data    y( 19),   y( 20),   y( 21),   y( 22),   y( 23),   y( 24) &
      / 170.0e0, 170.0e0, 158.0e0, 133.0e0, 114.0e0, 140.0e0/
  data    y( 25),   y( 26),   y( 27),   y( 28),   y( 29),   y( 30) &
      / 145.0e0, 150.0e0, 178.0e0, 163.0e0, 172.0e0, 178.0e0/
  data    y( 31),   y( 32),   y( 33),   y( 34),   y( 35),   y( 36) &
      / 199.0e0, 199.0e0, 184.0e0, 162.0e0, 146.0e0, 166.0e0/
  data    y( 37),   y( 38),   y( 39),   y( 40),   y( 41),   y( 42) &
      / 171.0e0, 180.0e0, 193.0e0, 181.0e0, 183.0e0, 218.0e0/
  data    y( 43),   y( 44),   y( 45),   y( 46),   y( 47),   y( 48) &
      / 230.0e0, 242.0e0, 209.0e0, 191.0e0, 172.0e0, 194.0e0/
  data    y( 49),   y( 50),   y( 51),   y( 52),   y( 53),   y( 54) &
      / 196.0e0, 196.0e0, 236.0e0, 235.0e0, 229.0e0, 243.0e0/
  data    y( 55),   y( 56),   y( 57),   y( 58),   y( 59),   y( 60) &
      / 264.0e0, 272.0e0, 237.0e0, 211.0e0, 180.0e0, 201.0e0/
  data    y( 61),   y( 62),   y( 63),   y( 64),   y( 65),   y( 66) &
      / 204.0e0, 188.0e0, 235.0e0, 227.0e0, 234.0e0, 264.0e0/
  data    y( 67),   y( 68),   y( 69),   y( 70),   y( 71),   y( 72) &
      / 302.0e0, 293.0e0, 259.0e0, 229.0e0, 203.0e0, 229.0e0/
  data    y( 73),   y( 74),   y( 75),   y( 76),   y( 77),   y( 78) &
      / 242.0e0, 233.0e0, 267.0e0, 269.0e0, 270.0e0, 315.0e0/
  data    y( 79),   y( 80),   y( 81),   y( 82),   y( 83),   y( 84) &
      / 364.0e0, 347.0e0, 312.0e0, 274.0e0, 237.0e0, 278.0e0/
  data    y( 85),   y( 86),   y( 87),   y( 88),   y( 89),   y( 90) &
      / 284.0e0, 277.0e0, 317.0e0, 313.0e0, 318.0e0, 374.0e0/
  data    y( 91),   y( 92),   y( 93),   y( 94),   y( 95),   y( 96) &
      / 413.0e0, 405.0e0, 355.0e0, 306.0e0, 271.0e0, 306.0e0/
  data    y( 97),   y( 98),   y( 99),   y(100),   y(101),   y(102) &
      / 315.0e0, 301.0e0, 356.0e0, 348.0e0, 355.0e0, 422.0e0/
  data    y(103),   y(104),   y(105),   y(106),   y(107),   y(108) &
      / 465.0e0, 467.0e0, 404.0e0, 347.0e0, 305.0e0, 336.0e0/
  data    y(109),   y(110),   y(111),   y(112),   y(113),   y(114) &
      / 340.0e0, 318.0e0, 362.0e0, 348.0e0, 363.0e0, 435.0e0/
  data    y(115),   y(116),   y(117),   y(118),   y(119),   y(120) &
      / 491.0e0, 505.0e0, 404.0e0, 359.0e0, 310.0e0, 337.0e0/
  data    y(121),   y(122),   y(123),   y(124),   y(125),   y(126) &
      / 360.0e0, 342.0e0, 406.0e0, 396.0e0, 420.0e0, 472.0e0/
  data    y(127),   y(128),   y(129),   y(130),   y(131),   y(132) &
      / 548.0e0, 559.0e0, 463.0e0, 407.0e0, 362.0e0, 405.0e0/
  data    y(133),   y(134),   y(135),   y(136),   y(137),   y(138) &
      / 417.0e0, 391.0e0, 419.0e0, 461.0e0, 472.0e0, 535.0e0/
  data    y(139),   y(140),   y(141),   y(142),   y(143),   y(144) &
      / 622.0e0, 606.0e0, 508.0e0, 461.0e0, 390.0e0, 432.0e0/
!
!  test against published results
!
  ny = 144
  ylog(1:144) = log(y(1:144))

  nfac = 2
  mspec(1,1) = 0
  mspec(2,1) = 1
  mspec(3,1) = 1
  mspec(4,1) = 1

  mspec(1,2) = 0
  mspec(2,2) = 1
  mspec(3,2) = 1
  mspec(4,2) = 12

  npar = 3
  par(1) = 0.0e0
  par(2) = 0.40e0
  par(3) = 0.60e0

  ifixed(1) = 1
  ifixed(2) = 0
  ifixed(3) = 0

  stopss = -1.0e0
  stopp = -1.0e0
  scale(1) = -1.0e0
  scale(2) = 1.0e-7
  scale(3) = 1.0e-7
  stp(1) = -1.0e0
  stp(2) = 1.0e-7
  stp(3) = 1.0e-7
  mit = 0
  nprt = -1
  delta = -1.0e0
  ivaprx = -1

  write ( *, 1000)
  call aimec ( ylog, ny, mspec, nfac, &
     par, npar, res,ldstak, ifixed,stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )

  write ( *, 1005)
  par(1) = 0.0e0
  par(2) = 0.395e0
  par(3) = 0.615e0
  call aimfs ( ylog, ny, mspec, nfac, &
     par, npar, ldstak, ny/10+1, 1, ny, nprt, fcst, 50, fcstsd )

  scale(1) = 1.0e-7
  scale(2) = 1.0e-7
  scale(3) = 1.0e-7

  nfac = 2
  mspec(1,1) = 0
  mspec(2,1) = 1
  mspec(3,1) = 1
  mspec(4,1) = 1

  mspec(1,2) = 0
  mspec(2,2) = 0
  mspec(3,2) = 1
  mspec(4,2) = 12

  write ( *, 1000)
  call aimec ( ylog, ny, mspec, nfac, &
     par, npar, res,ldstak, ifixed,stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
  ny = 20
  write ( *, 1000)
  call aimec ( ylog, ny, mspec, nfac, &
     par, npar, res,ldstak, ifixed,stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )

  nfac = 2
  mspec(1,1) = 0
  mspec(2,1) = 0
  mspec(3,1) = 1
  mspec(4,1) = 1

  mspec(1,2) = 0
  mspec(2,2) = 0
  mspec(3,2) = 1
  mspec(4,2) = 12

  ny = 144
  write ( *, 1000)
  call aimec ( ylog, ny, mspec, nfac, &
     par, npar, res,ldstak, ifixed,stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
!
!  example from page 212 of box and jenkins (1970)
!  add print statements to mdlts2 to check computations
!  at first call against those listed on page 214.
!
  write ( *, 1000)
  ny = 10
  yt(1) = 460.0e0
  yt(2) = 457.0e0
  yt(3) = 452.0e0
  yt(4) = 459.0e0
  yt(5) = 462.0e0
  yt(6) = 459.0e0
  yt(7) = 463.0e0
  yt(8) = 479.0e0
  yt(9) = 493.0e0
  yt(10) = 490.0e0

  nfac = 1
  mspec(1,1) = 0
  mspec(2,1) = 1
  mspec(3,1) = 1
  mspec(4,1) = 1

  npar = 2
  par(1) = 0.0e0
  par(2) = 0.5e0

  ifixed(1) = 1
  ifixed(2) = 0

  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res,ldstak, ifixed,stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
!
!  example from page 216 of box and jenkins (1970)
!  add print statements to mdlts2 to check computations
!  at first call against those listed on page 218.
!
  write ( *, 1000)
  ny = 12
  yt(1) = 2.0e0
  yt(2) = 0.8e0
  yt(3) = -0.3e0
  yt(4) = -0.3e0
  yt(5) = -1.9e0
  yt(6) = 0.3e0
  yt(7) = 3.2e0
  yt(8) = 1.6e0
  yt(9) = -0.7e0
  yt(10) = 3.0e0
  yt(11) = 4.3e0
  yt(12) = 1.1e0

  nfac = 1
  mspec(1,1) = 1
  mspec(2,1) = 0
  mspec(3,1) = 1
  mspec(4,1) = 1

  npar = 3
  par(1) = 0.3e0
  par(2) = 0.0e0
  par(3) = 0.7e0

  ifixed(1) = 0
  ifixed(2) = 1
  ifixed(3) = 0

  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
!
!  test error messages
!
  write ( *, 1010)
  ny = 0
  nfac = 0
  call aime ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak )
  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
  call aimes ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt, npare, rsd, pv, sdpv, sdres, vcv, &
     ivcv )
  call aimf ( y, ny, mspec, nfac, par, npar, ldstak )
  call aimfs ( y, ny, mspec, nfac, &
     par, npar, ldstak, ny/10+1, 1, ny, nprt, fcst, 50, fcstsd )

  ny = 144
  nfac = 2
  mspec(1,1) = -1
  call aime ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak )
  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
  call aimes ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt, npare, rsd, pv, sdpv, sdres, vcv, &
     ivcv )
  call aimf ( y, ny, mspec, nfac, par, npar, ldstak )
  call aimfs ( y, ny, mspec, nfac, &
     par, npar, ldstak, ny/10+1, 1, ny, nprt, fcst, 50, fcstsd )
  ny = 144
  nfac = 2
  mspec(1,1) = 0
  npar = 1
  call aime ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak )
  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
  call aimes ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt, npare, rsd, pv, sdpv, sdres, vcv, &
     ivcv )
  call aimf ( y, ny, mspec, nfac, par, npar, ldstak )
  call aimfs ( y, ny, mspec, nfac, &
     par, npar, ldstak, ny/10+1, 1, ny, nprt, fcst, 50, fcstsd )
  ny = 144
  nfac = 2
  mspec(1,1) = 0
  npar = 3
  ifixed(1:npar) = 1
  ivcv = 0
  ifcst = 0
  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
  call aimes ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt, npare, rsd, pv, sdpv, sdres, vcv, &
     ivcv )
  call aimfs ( y, ny, mspec, nfac, &
     par, npar, ldstak, ny/10+1, 1, ny, nprt, fcst, ifcst, fcstsd )
  ifixed(1:npar) = 1
  ivcv = 0
  stp(2) = -1.0e0
  scale(2) = -1.0e0
  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
  return

 1000 format ('1test of arima estimation routines')
 1005 format ('1test of arima forecasting routines')
 1010 format ('1test of error checking facilities')
end
subroutine xaov1 ( ldstak )

!*****************************************************************************80
!
!! XAOV1 exercises the oneway family routines.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Linda Mitchell,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fi
!        float of index i
!     real gstat(10,4)
!        the group statistics.  columns correspond to the tag
!        value, sample size, group mean, and group standard deviation.
!     integer i
!        index variable
!     integer ierr
!        common flag indicating whether or not there were any errors
!     integer igstat
!        the first dimension of gstat.
!     integer ldsmin
!        the smallest acceptable size of the common cstak
!     integer ldstak
!        the size of the common cstak
!     integer n
!        the number of observations
!     integer ng
!        the number of different groups
!     real tag(20)
!        the tag values for each observation
!     real y(20)
!        the vector of observations
!     real z(10)
!        test vector
!     real ztag(10)
!        test tag vector
!
  implicit none

  integer ierr
  integer ldstak
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     fi
  integer &
     i,igstat,ldsmin,n,ng
!
!  local arrays
  real &
     gstat(10,4),tag(20),y(20),z(10),ztag(10)
!
!  external subroutines
  external aov1,aov1s,aov1xp,iprint,ldscmp,msgx,setrv
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data   y( 1),   y( 2),   y( 3),   y( 4),   y( 5) &
     /    61.0e0,    61.0e0,    67.0e0,    67.0e0,    64.0e0/
  data   y( 6),   y( 7),   y( 8),   y( 9),   y(10) &
     /    78.0e0,    71.0e0,    75.0e0,    72.0e0,    74.0e0/
  data   y(11),   y(12),   y(13),   y(14),   y(15) &
     /    83.0e0,    81.0e0,    76.0e0,    78.0e0,    79.0e0/
  data   y(16),   y(17) &
     /    72.0e0,   72.0e0/

  data tag( 1), tag( 2), tag( 3), tag( 4), tag( 5) &
     /    11.5e0,    11.5e0,    11.5e0,    11.5e0,    11.5e0/
  data tag( 6), tag( 7), tag( 8), tag( 9), tag(10) &
     /    12.0e0,    12.0e0,    12.0e0,    12.0e0,    12.0e0/
  data tag(11), tag(12), tag(13), tag(14), tag(15) &
     /    11.0e0,    11.0e0,    11.0e0,    11.0e0,    11.0e0/
  data tag(16), tag(17) &
     /   -11.0e0,   11.0e0/
!
!  set various dimensions and program variables
!
  n = 17
  igstat = 10
!
! test with correct call statements.
!
  write ( *,1060)
!
!  test aov1
!
  write ( *,1000)
  call aov1( y, tag, n, ldstak )
  call msgx ( 0 )
  write ( *,1010)
  write ( *,1020) (tag(i),i=1,n)
!
!  test of aov1s
!
!  printout not supressed
!
  write ( *,1040)
  write ( *,1030)
  call aov1s(y, tag, n, ldstak, 1, gstat, igstat, ng)
  call msgx(0 )
!
!  print storage and zero vectors before using again
!
  call aov1xp(gstat, igstat, ng)
!
!  printout supressed
!
  write ( *,1050)
  call aov1s(y, tag, n, ldstak, 0, gstat, igstat, ng)
  call msgx(0 )
!
!  print storage and zero vectors before using again
!
  call aov1xp(gstat, igstat, ng)
!
!  number of observations less than 2
!
  write ( *,1090)
  call aov1(y, tag, 1, ldstak)
  call msgx(1 )

  call aov1s(y, tag, -14, ldstak, 0, gstat, igstat, ng)
  call msgx(1 )
!
!  all observations the same value
!
  write ( *,1100)
  write ( *,1000)
  call setrv(z, 10, 0.0e0)
  call aov1(z, tag, 10, ldstak)
  call msgx(0 )

  call setrv(z, 10, 2.0e0)
  write ( *,1100)
  write ( *,1030)
  call aov1s(z, tag, 10, ldstak, 1, gstat, igstat, ng)
  call msgx(0 )
!
!  test work area size handling
!
  call ldscmp(11, 0, 33, 0, 0, 0, 's', 40, ldsmin)
  write ( *,1070)
  write ( *,1000)
  call aov1(y, tag, n, 1)
  call msgx(1 )
  write ( *,1070)
  write ( *,1000)
  call aov1(y, tag, n, ldsmin-1)
  call msgx(1 )
  write ( *,1010)
  write ( *,1020) (tag(i),i=1,n)
  write ( *,1080)
  write ( *,1000)
  call aov1(y, tag, n, ldsmin)
  call msgx(0 )

  call ldscmp(11, 0, 33, 0, 0, 0, 's', 28, ldsmin)
  write ( *,1070)
  write ( *,1030)
  call aov1s(y, tag, n, 1, 0, gstat, igstat, ng)
  call msgx(1 )
  write ( *,1070)
  write ( *,1030)
  call aov1s(y, tag, n, ldsmin-1, 0, gstat, igstat, ng)
  call msgx(1 )
  write ( *,1010)
  write ( *,1020) (tag(i),i=1,n)
  write ( *,1080)
  write ( *,1030)
  call aov1s(y, tag, n, ldsmin, 1, gstat, igstat, ng)
  call msgx(0 )
!
!  same number of groups as non-zero tags
!
  write ( *,1120)
  do i=1,10
     fi = i
     ztag(i) = fi
     z(i) = 13.0e0 - fi
  end do
  call aov1(z, ztag, 10, ldstak)
  call msgx(1 )
  call aov1s(z, ztag, 10, ldstak, 1, gstat, igstat, ng)
  call msgx(1 )
!
!  less than 2 different tag groups
!
  write ( *,1130)
  call setrv(ztag, 10, 1.0e0)
  call aov1(z, ztag, 10, ldstak)
  call msgx(1 )
  call aov1s(z, ztag, 10, ldstak, 1, gstat, igstat, ng)
  call msgx(1 )
!
!  less than 2 tags
!
  call setrv(ztag, 9, 0.0e0)
  write ( *,1140)
  call aov1(z, ztag, 10, ldstak)
  call msgx(1 )
  call aov1s(z, ztag, 10, ldstak, 1, gstat, igstat, ng)
  call msgx(1 )
!
!  incorrect dimension of gstat
!
  write ( *,1150)
  call aov1s(y, tag, n, ldstak, 1, gstat, 2, ng)
  call msgx(1 )
!
!  all observations within a group same value
!
  z(1) = 53.0e0
  ztag(1) = 1.0e0
  z(2) = 62.0e0
  ztag(2) = 3.0e0
  z(3) = 53.0e0
  ztag(3) = 1.0e0
  z(4) = 71.0e0
  ztag(4) = 4.0e0
  z(5) = 89.0e0
  ztag(5) = 2.0e0
  z(6) = 71.0e0
  ztag(6) = 4.0e0
  z(7) = 89.0e0
  ztag(7) = 2.0e0
  z(8) = 62.0e0
  ztag(8) = 3.0e0
  z(9) = 71.0e0
  ztag(9) = 4.0e0
  z(10) = 62.0e0
  ztag(10) = 3.0e0
  write ( *,1160)
  call aov1(z, ztag, 10, ldstak)
  call msgx(0 )
!
!  2 tags
!
  write ( *,1170)
  call aov1(z, ztag, 3, ldstak)
  call msgx(0 )
!
!  all groups(except for 1) with 1 observation
!
  write ( *,1180)
  call aov1(z, ztag, 5, ldstak)
  call msgx(0 )

  return

 1000 format(' test of aov1 ')
 1010 format(' check to see if tags have been changed')
 1020 format(4f12.6)
 1030 format(' test of aov1s ')
 1040 format('1printout not supressed.')
 1050 format(' printout supressed.')
 1060 format('1****test routines with correct call****')
 1070 format('1****test with insufficient work area****')
 1080 format('1****test with exactly the right amount of work area****')
 1090 format('1****number of observations less than 2****')
 1100 format('1****all observations with same value****')
 1120 format('1****same number of groups as non-zero tags****')
 1130 format(' ****less than 2 different tag groups****')
 1140 format(' ****less than 2 tags****')
 1150 format('1****incorrect dimension of gstat****')
 1160 format('1****all observations within a group same value****')
 1170 format('1****test with 2 tags****')
 1180 format('1****all groups except for 1 with 1 observation ****')
end
subroutine xbfs ( lds )

!*****************************************************************************80
!
!! XBFS tests the Fourier spectrum analysis routines.
!
!  Discussion:
!
!    series y1 and y2 are listed as series x1 and x2 on page of 361 of
!    jenkins and watts.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    Gwilym Jenkins, Donald Watts,
!    Spectral Analysis and its Applications,
!    Holden-Day 1968.
!
!  Parameters:
!
!     real cmiss
!         the missing value code for the returned ccvf estimates.
!     real ccov(101,2,2)
!        the covariances.
!     real cspc2(300,2)
!        the squared coherency component of the spectrum.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fmax, fmin
!        the maximum and minimum frequencies at which the
!        spectrum is to be computed.
!     real freq(300)
!        the vector of frequencies at which the spectrum is computed.
!     integer i
!        an index variable
!     integer iccov
!        the first dimension of the array ccov.
!     integer icspc2
!        the first dimension of the array cspc2.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list
!        if ierr == 0, no errors were detected
!        if ierr == 1, errors have been detected
!     integer iphas
!        the first dimension of the array phas.
!     integer ispcf
!         the actual dimension for the spectrum arrays.
!     integer j
!        index variable.
!     integer lacov
!        the length of the vector acov.
!     integer lagmax
!        the indexing variable indicating the lag value of the
!        autocovariance being computed and the maximum lag to be used,
!        respectively.
!     integer lags(4)
!        the array used to store the lag window truccation
!        points used for each set of spectrum values.
!     integer lds, ldstak
!        the length of the vector dstak in common cstak.
!     integer lyfft
!        the length of the vectors yfft and zfft, respectively..
!     integer nf
!        the number of frequencies at which the spectrum is
!        to be computed.
!     integer nlppc(101,2,2)
!        the numbers of lagged product pairs used for each acvf.
!     integer nprt
!        a code used to specify the type of plot.
!        if nprt < 0 the plot is decibles/linear
!        if nprt = 0 the plot is suppressed.
!        if nprt > 0 the plot is log/linear
!     integer nw
!        the number of different lag window truncation points specified,
!        and therefore, the number of plots.
!     integer n
!        the number of observations in the series y.
!     real phas(300,2)
!        the phase component of the spectrum.
!     real yfft1(400), yfft2(400)
!        the vectors of the observed time series to be analyzed using
!        the fft.
!     real ymiss, ymiss1, ymiss2, ymmiss(4)
!        the user supplied code which is used to determine whether or
!        not an observation in the series is missing.  if y(i) = ymiss,
!        the value is assumed missing, otherwise it is not.
!     real ym(150,2)
!        the multivariate representation of the data
!     real y1(150), y2(150)
!         the vectors containing the time series from jenkins and watts.
!
  implicit none

  integer &
     lds
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     cmiss,fmax,fmin,ymiss,ymiss1,ymiss2
  integer &
     i,iccov,icspc2,index1,index2,inlppc,iphas,ispcf,j, &
     jccov,jnlppc,lacov,lagmax,ldstak,lyfft,n,nf,nprt,nw
!
!  local arrays
  real &
     ccov(101,2,2),cspc2(300,2),freq(300),phas(300,2),y1(150), &
     y2(150),yfft1(400),yfft2(400),ym(150,2),ymmiss(4)
  integer &
     lags(4),nlppc(101,2,2)
!
!  external subroutines
  external bfs,bfsf,bfsfs,bfsm,bfsms,bfsmv,bfsmvs,bfss,bfsv,bfsvs, &
     ccfms,ccfs,iprint,nrand,scopy,setiv,setrv
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!  equivalences
  equivalence (ym(1,1),y1(1))
  equivalence (ym(1,2),y2(1))

  data   y1(  1),  y1(  2),  y1(  3),  y1(  4),  y1(  5),  y1(  6) &
      /-0.88e0, -0.16e0, -1.87e0, -1.12e0,  1.38e0,  2.13e0/
  data   y1(  7),  y1(  8),  y1(  9),  y1( 10),  y1( 11),  y1( 12) &
      / 2.76e0,  0.56e0, -0.69e0, -1.79e0, -3.82e0, -2.38e0/
  data   y1( 13),  y1( 14),  y1( 15),  y1( 16),  y1( 17),  y1( 18) &
      / 1.00e0,  0.70e0, -0.15e0,  0.98e0,  0.11e0, -0.35e0/
  data   y1( 19),  y1( 20),  y1( 21),  y1( 22),  y1( 23),  y1( 24) &
      /-0.73e0,  0.89e0, -1.63e0, -0.44e0, -1.37e0, -1.71e0/
  data   y1( 25),  y1( 26),  y1( 27),  y1( 28),  y1( 29),  y1( 30) &
      /-1.22e0, -2.00e0, -0.22e0,  0.38e0,  1.31e0,  0.71e0/
  data   y1( 31),  y1( 32),  y1( 33),  y1( 34),  y1( 35),  y1( 36) &
      / 0.32e0,  0.48e0, -1.88e0, -0.94e0, -1.54e0, -0.13e0/
  data   y1( 37),  y1( 38),  y1( 39),  y1( 40),  y1( 41),  y1( 42) &
      / 1.02e0,  0.02e0, -0.77e0,  0.11e0, -0.60e0, -0.52e0/
  data   y1( 43),  y1( 44),  y1( 45),  y1( 46),  y1( 47),  y1( 48) &
      /-0.09e0,  1.23e0,  1.46e0,  0.61e0,  0.42e0,  2.16e0/
  data   y1( 49),  y1( 50),  y1( 51),  y1( 52),  y1( 53),  y1( 54) &
      / 3.18e0,  2.10e0,  0.37e0, -0.24e0,  0.57e0, -0.53e0/
  data   y1( 55),  y1( 56),  y1( 57),  y1( 58),  y1( 59),  y1( 60) &
      / 2.44e0,  1.02e0, -0.53e0, -2.49e0, -2.12e0, -1.04e0/
  data   y1( 61),  y1( 62),  y1( 63),  y1( 64),  y1( 65),  y1( 66) &
      /-0.12e0, -1.88e0, -1.50e0,  1.54e0,  3.33e0,  3.08e0/
  data   y1( 67),  y1( 68),  y1( 69),  y1( 70),  y1( 71),  y1( 72) &
      / 1.71e0,  0.79e0,  1.55e0,  0.89e0, -0.89e0, -1.18e0/
  data   y1( 73),  y1( 74),  y1( 75),  y1( 76),  y1( 77),  y1( 78) &
      / 0.89e0,  1.71e0,  3.05e0,  0.15e0, -1.04e0,  0.12e0/
  data   y1( 79),  y1( 80),  y1( 81),  y1( 82),  y1( 83),  y1( 84) &
      / 0.08e0,  0.11e0, -2.62e0, -1.28e0,  1.07e0,  3.20e0/
  data   y1( 85),  y1( 86),  y1( 87),  y1( 88),  y1( 89),  y1( 90) &
      / 1.92e0,  0.53e0, -1.08e0,  0.49e0, -0.58e0,  0.17e0/
  data   y1( 91),  y1( 92),  y1( 93),  y1( 94),  y1( 95),  y1( 96) &
      / 1.15e0, -0.97e0, -1.63e0,  1.14e0, -0.67e0, -0.88e0/
  data   y1( 97),  y1( 98),  y1( 99),  y1(100) &
      /-0.07e0,  0.24e0,  0.55e0, -2.16e0/
  data   y2(  1),  y2(  2),  y2(  3),  y2(  4),  y2(  5),  y2(  6) &
      / 0.79e0,  1.12e0, -1.10e0, -2.39e0, -1.75e0, -0.82e0/
  data   y2(  7),  y2(  8),  y2(  9),  y2( 10),  y2( 11),  y2( 12) &
      /-0.36e0,  1.27e0,  1.75e0,  2.44e0,  0.36e0, -2.10e0/
  data   y2( 13),  y2( 14),  y2( 15),  y2( 16),  y2( 17),  y2( 18) &
      /-1.93e0, -1.30e0, -1.75e0, -0.34e0,  0.74e0,  0.49e0/
  data   y2( 19),  y2( 20),  y2( 21),  y2( 22),  y2( 23),  y2( 24) &
      / 0.70e0,  0.71e0,  0.09e0,  0.59e0,  1.54e0,  0.14e0/
  data   y2( 25),  y2( 26),  y2( 27),  y2( 28),  y2( 29),  y2( 30) &
      / 0.55e0, -1.40e0, -2.55e0, -1.66e0, -0.43e0,  0.58e0/
  data   y2( 31),  y2( 32),  y2( 33),  y2( 34),  y2( 35),  y2( 36) &
      / 2.18e0, -0.24e0,  0.58e0, -0.18e0, -1.55e0, -0.64e0/
  data   y2( 37),  y2( 38),  y2( 39),  y2( 40),  y2( 41),  y2( 42) &
      /-1.09e0,  0.90e0, -0.66e0, -0.35e0,  0.48e0,  0.50e0/
  data   y2( 43),  y2( 44),  y2( 45),  y2( 46),  y2( 47),  y2( 48) &
      / 0.05e0, -0.68e0,  0.24e0,  0.58e0, -1.26e0, -0.25e0/
  data   y2( 49),  y2( 50),  y2( 51),  y2( 52),  y2( 53),  y2( 54) &
      / 0.25e0,  2.18e0,  2.96e0,  1.56e0, -0.36e0, -0.59e0/
  data   y2( 55),  y2( 56),  y2( 57),  y2( 58),  y2( 59),  y2( 60) &
      /-0.12e0,  3.03e0,  2.11e0,  0.78e0,  0.89e0, -1.45e0/
  data   y2( 61),  y2( 62),  y2( 63),  y2( 64),  y2( 65),  y2( 66) &
      /-0.36e0, -0.37e0, -1.39e0, -4.19e0, -0.73e0, -0.98e0/
  data   y2( 67),  y2( 68),  y2( 69),  y2( 70),  y2( 71),  y2( 72) &
      / 0.36e0,  0.06e0, -1.94e0, -0.08e0,  0.17e0,  1.00e0/
  data   y2( 73),  y2( 74),  y2( 75),  y2( 76),  y2( 77),  y2( 78) &
      /-0.05e0,  0.43e0,  0.15e0,  2.69e0,  0.57e0,  0.29e0/
  data   y2( 79),  y2( 80),  y2( 81),  y2( 82),  y2( 83),  y2( 84) &
      / 1.10e0,  0.48e0, -1.06e0, -2.28e0, -2.03e0, -0.75e0/
  data   y2( 85),  y2( 86),  y2( 87),  y2( 88),  y2( 89),  y2( 90) &
      / 1.00e0,  1.71e0,  0.58e0,  1.97e0,  0.99e0,  1.94e0/
  data   y2( 91),  y2( 92),  y2( 93),  y2( 94),  y2( 95),  y2( 96) &
      / 2.18e0,  3.14e0,  0.60e0,  0.51e0,  1.35e0,  0.56e0/
  data   y2( 97),  y2( 98),  y2( 99),  y2(100) &
      / 0.11e0,  0.00e0,  2.34e0,  1.88e0/

  call setrv(ymmiss, 4, 0.89e0)
!
!  check error handling
!
!  test 1  -  miscellaneous error checking
!
  write ( *, 2000)
  lagmax = -1
  n = -10
  index1 = 0
  index2 = 0
  iccov = 0
  jccov = 0
  inlppc = 0
  jnlppc = 0
  icspc2 = -10
  iphas = -10
  lacov = -11
  lyfft = -11
  nw = -1
  nf = -5
  fmin = 0.5e0
  fmax = 0.0e0
  nprt = -1
  ispcf = -20
  ldstak = 0
  ymiss1 = 0.89e0
  ymiss2 = 0.89e0

  write ( *, 1001)
  call bfs (y1, y2, n)
  write ( *, 1002) ierr

  write ( *, 1003)
  call bfss(y1, y2, n, nw, lags, nf, fmin, fmax, nprt, &
     cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1019)
  call bfsf (yfft1, yfft2, n, lyfft, ldstak)
  write ( *, 1002) ierr

  write ( *, 1020)
  call bfsfs(yfft1, yfft2, n, lyfft, ldstak, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq)
  write ( *, 1002) ierr

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of BFSM'
  write ( *, '(a)' ) ' '

  call bfsm (y1, ymiss1, y2, ymiss2, n)
  write ( *, 1002) ierr

  write ( *, 1006)
  call bfsms(y1, ymiss1, y2, ymiss2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1007)
  call bfsv(ccov, index1, index2, n, lagmax, iccov, jccov)
  write ( *, 1002) ierr

  write ( *, 1008)
  call bfsvs (ccov, index1, index2, n, iccov, jccov, &
     nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, phas, iphas, &
     freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1021)
  call bfsmv(ccov, nlppc, index1, index2, n, lagmax, iccov, jccov, &
     inlppc, jnlppc)
  write ( *, 1002) ierr

  write ( *, 1022)
  call bfsmvs (ccov, nlppc, index1, index2, n, iccov, jccov, &
     inlppc, jnlppc, nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, &
     phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  test 2  -  miscellaneous error checking (continued)
!
  write ( *, 2010)
  n = 100
  lagmax = 40
  index1 = 0
  index2 = 0
  iccov = 0
  jccov = 0
  inlppc = 0
  jnlppc = 0
  icspc2 = 300
  iphas = 300
  lacov = 101
  lyfft = -11
  nw = 2
  lags(1) = 0
  lags(2) = 100
  nf = 202
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 3
  ispcf = 101
  ldstak = 0

  write ( *, 1003)
  call bfss(y1, y2, n, nw, lags, nf, fmin, fmax, nprt, &
     cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1019)
  call bfsf (yfft1, yfft2, n, lyfft, ldstak)
  write ( *, 1002) ierr

  write ( *, 1020)
  call bfsfs(yfft1, yfft2, n, lyfft, ldstak, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq)
  write ( *, 1002) ierr

  write ( *, 1006)
  call bfsms(y1, ymiss1, y2, ymiss2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1008)
  call bfsvs (ccov, index1, index2, n, iccov, jccov, &
     nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, phas, iphas, &
     freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1022)
  call bfsmvs (ccov, nlppc, index1, index2, n, iccov, jccov, &
     inlppc, jnlppc, nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, &
     phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  test 3  -  ldstak too small
!
  write ( *, 2030)
  n = 100
  index1 = 2
  index2 = 1
  iccov = 101
  jccov = 2
  inlppc = 101
  jnlppc = 2
  icspc2 = 300
  iphas = 300
  lagmax = 99
  lacov = 101
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 26
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = 0

  write ( *, 1003)
  call bfss(y1, y2, n, nw, lags, nf, fmin, fmax, nprt, &
     cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1019)
  call bfsf (yfft1, yfft2, n, lyfft, ldstak)
  write ( *, 1002) ierr

  write ( *, 1020)
  call bfsfs(yfft1, yfft2, n, lyfft, ldstak, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq)
  write ( *, 1002) ierr

  write ( *, 1006)
  call bfsms(y1, ymiss1, y2, ymiss2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1008)
  call bfsvs (ccov, index1, index2, n, iccov, jccov, &
     nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, phas, iphas, &
     freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1022)
  call bfsmvs (ccov, nlppc, index1, index2, n, iccov, jccov, &
     inlppc, jnlppc, nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, &
     phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  test 4  -  all data and covariances missing
!
  write ( *, 2040)
  n = 100
  lagmax = 99
  icspc2 = 300
  iphas = 300
  lacov = 101
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 26
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = lds
  call setrv(yfft1, n, ymiss1)
  call setrv(yfft2, n, ymiss2)
  call setrv(ccov, 404, 0.0e0)
  call setiv(nlppc, 404, 0)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of BFSM'
  write ( *, '(a)' ) ' '

  call bfsm (yfft1, ymiss1, yfft2, ymiss2, n)
  write ( *, 1002) ierr

  write ( *, 1006)
  call bfsms(yfft1, ymiss1, yfft2, ymiss2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1021)
  call bfsmv(ccov, nlppc, index1, index2, n, lagmax, iccov, jccov, &
     inlppc, jnlppc)
  write ( *, 1002) ierr

  write ( *, 1022)
  call bfsmvs (ccov, nlppc, index1, index2, n, iccov, jccov, &
     inlppc, jnlppc, nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, &
     phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  test 5  -  every other value missing
!
  write ( *, '(a)' ) '  Every other data value missing.'
  n = 100
  lagmax = 99
  icspc2 = 300
  iphas = 300
  lacov = 101
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 26
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = lds
  call setrv(yfft1, n, ymiss1)
  call setrv(yfft2, n, ymiss2)
  yfft1(1:n:2) = y1(1:n:2)
  yfft2(1:n:2) = y2(1:n:2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of BFSM'
  write ( *, '(a)' ) ' '

  call bfsm ( yfft1, ymiss1, yfft2, ymiss2, n )
  write ( *, 1002) ierr

  write ( *, 1006)
  call bfsms(yfft1, ymiss1, yfft2, ymiss2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  check results from valid call
!
  write ( *, 2020)
  ymiss = 1.16e0
  n = 100
  lagmax = 99
  icspc2 = 300
  iphas = 300
  lacov = 101
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 26
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = lds
!
!  test of bfs
!
  write ( *, 1001)
  call bfs (y1, y2, n)
  write ( *, 1002) ierr
!
!  test of bfss
!
  write ( *, 2020)
  write ( *, 1003)
  call bfss(y1, y2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  print returned values from bfss
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  test of bfsf
!
  write ( *, 2020)
  write ( *, 1019)
  call scopy(n, y1, 1, yfft1, 1)
  call scopy(n, y2, 1, yfft2, 1)
  call bfsf (yfft1, yfft2, n, lyfft, ldstak)
  write ( *, 1002) ierr
!
!  test of bfsfs
!
  write ( *, 2020)
  write ( *, 1020)
  call scopy(n, y1, 1, yfft1, 1)
  call scopy(n, y2, 1, yfft2, 1)
  call bfsfs(yfft1, yfft2, n, lyfft, ldstak, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq)
  write ( *, 1002) ierr
!
!  print returned values from bfsfs
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  test of bfsm
!
  write ( *, 2020)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of BFSM'
  write ( *, '(a)' ) ' '

  call bfsm (y1, ymiss1, y2, ymiss2, n)
  write ( *, 1002) ierr
!
!  test of bfsms
!
  write ( *, 2020)
  write ( *, 1006)
  call bfsms(y1, ymiss1, y2, ymiss2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  print returned values from bfsms
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  test of bfsv
!
  write ( *, 2020)
  call ccfs (ym, n, 2, 150, lagmax, ccov, iccov, jccov, 0, &
     ldstak)
  write ( *, 1007)
  call bfsv(ccov, index1, index2, n, lagmax, iccov, jccov)
  write ( *, 1002) ierr
!
!  test of bfsvs
!
  write ( *, 2020)
  write ( *, 1008)
  call bfsvs(ccov, index1, index2, n, iccov, jccov, nw, lags, nf, &
     fmin, fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  print returned values from bfsvs
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  test of bfsmv
!
  write ( *, 2020)
  call ccfms (ym, ymmiss, n, 2, 150, lagmax, ccov, cmiss, iccov, &
    jccov, nlppc, inlppc, jnlppc, 0, ldstak)
  write ( *, 1021)
  call bfsmv(ccov, nlppc, index1, index2, n, lagmax, iccov, &
    jccov, inlppc, jnlppc)
  write ( *, 1002) ierr
!
!  test of bfsmvs
!
  write ( *, 2020)
  write ( *, 1022)
  call bfsmvs(ccov, nlppc, index1, index2, n, iccov, &
    jccov, inlppc, jnlppc, nw, lags, nf, fmin, fmax, nprt, &
    cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  print returned values from bfsmvs
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  minimum problem size
!
  ymiss = 1.16e0
  n = 17
  lagmax = 1
  icspc2 = 1
  iphas = 1
  lacov = 101
  lyfft = 400
  nw = 1
  lags(1) = 1
  lags(2) = 16
  nf = 1
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = lds
!
!  test of bfs
!
  write ( *, 2060)
  write ( *, 1001)
  call bfs(y1, y2, n)
  write ( *, 1002) ierr
!
!  test of bfss
!
  write ( *, 2060)
  write ( *, 1003)
  call bfss(y1, y2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  print returned values from bfss
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  check handling of fmin and fmax
!
  ymiss = 1.16e0
  n = 100
  lagmax = 99
  icspc2 = 300
  iphas = 300
  lacov = 101
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 26
  fmin = 0.45e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = lds
!
!  test of bfss
!
  write ( *, 2070)
  write ( *, 1003)
  call bfss(y1, y2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  print returned values from bfss
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  check results for white noise spectrum
!
  ymiss = 1.16e0
  call nrand(yfft1, n, 12343)
  call nrand(yfft2, n, 34523)
  n = 100
  lagmax = 99
  icspc2 = 300
  iphas = 300
  lacov = 101
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 26
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = lds
!
!  test of bfsfs
!
  write ( *, 2080)
  write ( *, 1003)
  call bfsfs(yfft1, yfft2, n, lyfft, ldstak, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq)
  write ( *, 1002) ierr
!
!  print returned values from bfss
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)

  return

 1001 format (' test of bfs')
 1002 format (' ierr is', i5/)
 1003 format (' test of bfss')
 1004 format (5(1x, e15.7))
 1006 format (' test of bfsms')
 1007 format (' test of bfsv')
 1008 format (' test of bfsvs')
 1019 format (' test of bfsf')
 1020 format (' test of bfsfs')
 1021 format (' test of bfsmv')
 1022 format (' test of bfsmvs')
 2000 format ('1check error handling  -  test 1')
 2010 format ('1check error handling  -  test 2')
 2020 format ('1valid problem')
 2030 format ('1lds too small')
 2040 format ('1all data and covariances missing')
 2060 format ('1minimum problem size')
 2070 format ('1check handling of fmin and fmax')
 2080 format ('1white noise spectrum')
end
subroutine xccf ( lds )

!*****************************************************************************80
!
!! XCCF tests the time series correlation subroutines.
!
!  Discussion:
!
!    series y1 and y2 are listed as series x1 and x2 on page of 361 of
!    jenkins and watts.  ccf for series y1 and y2 are plotted on page 3
!    and listed on page 420.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    Gwilym Jenkins, Donald Watts,
!    Spectral Analysis and its Applications,
!    Holden-Day 1968.
!
!  Parameters:
!
!     real ccov(30,5,5)
!        the cross covariance array.
!     real cmiss
!        the missing value code for the returned ccvf estimates
!        (vector ccov).
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer iccov
!        the first dimension of the array ccov.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list
!        if ierr == 0, no errors were detected
!     integer inlppc
!        the first dimension of the array nlppc.
!     integer itest
!        the number of the test being performed
!     integer iym, iymfft
!        the first dimension of the arrays ym and ymfft, respectively.
!     integer jccov, jnlppc
!        the second dimensions of the arrays ccov and nlppc,
!        respectively.
!     integer lagmax
!        the maximum lag value requested.
!     integer lds, ldstak
!        the length of the array dstak.
!     integer lyfft
!        the length of the arrays used when the computations are
!        performed by the fft.
!     integer m
!        the number of series in the multivariate time series ym.
!     integer n
!        the integer number of observations in each series
!     integer nlag
!        the number of lags at which the acvf was computed.
!     integer nlppc(30,5,5)
!        the array containing the number of lagged product pairs
!        used to compute each acvf estimate.
!     integer nprt
!        the indicator variable used to specify whether or not
!        printed output is to be given, where if the value of
!        nprt is zero, no output is made.
!     integer nyd
!        the number of observations in the series to be differenced.
!     real yfft1(150), yfft2(150)
!        the vectors used for storing the series for the routines
!        using the fft.
!     real ym(150,5), ymfft(150,5)
!        the arrays used for multivariate time series.
!     real ymiss0, ymmiss(5)
!        the missing value codes for series y and ym.
!     real y1(100), y1miss
!        the first series, and its missing value code.
!     real y2(100), y2miss
!        the second series, and its missing value code.
!
  implicit none

  integer &
     lds
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     cmiss,y1miss,y2miss,ymiss0
  integer &
     iccov,inlppc,itest,iym,iymfft,jccov,jnlppc,lagmax, &
     ldstak,lyfft,m,n,nlag,nprt,nyd
!
!  local arrays
  real &
     ccov(30,5,5),y1(100),y2(100),yfft1(150),yfft2(150),ym(150,5), &
     ymfft(150,5),ymmiss(5)
  integer &
     nlppc(30,5,5)
!
!  external subroutines
  external ccf,ccff,ccffs,ccfm,ccfms,ccfs,ccfxp,iprint,scopy,setra, &
     setrv
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data   y1(  1),  y1(  2),  y1(  3),  y1(  4),  y1(  5),  y1(  6) &
      /-0.88e0, -0.16e0, -1.87e0, -1.12e0,  1.38e0,  2.13e0/
  data   y1(  7),  y1(  8),  y1(  9),  y1( 10),  y1( 11),  y1( 12) &
      / 2.76e0,  0.56e0, -0.69e0, -1.79e0, -3.82e0, -2.38e0/
  data   y1( 13),  y1( 14),  y1( 15),  y1( 16),  y1( 17),  y1( 18) &
      / 1.00e0,  0.70e0, -0.15e0,  0.98e0,  0.11e0, -0.35e0/
  data   y1( 19),  y1( 20),  y1( 21),  y1( 22),  y1( 23),  y1( 24) &
      /-0.73e0,  0.89e0, -1.63e0, -0.44e0, -1.37e0, -1.71e0/
  data   y1( 25),  y1( 26),  y1( 27),  y1( 28),  y1( 29),  y1( 30) &
      /-1.22e0, -2.00e0, -0.22e0,  0.38e0,  1.31e0,  0.71e0/
  data   y1( 31),  y1( 32),  y1( 33),  y1( 34),  y1( 35),  y1( 36) &
      / 0.32e0,  0.48e0, -1.88e0, -0.94e0, -1.54e0, -0.13e0/
  data   y1( 37),  y1( 38),  y1( 39),  y1( 40),  y1( 41),  y1( 42) &
      / 1.02e0,  0.02e0, -0.77e0,  0.11e0, -0.60e0, -0.52e0/
  data   y1( 43),  y1( 44),  y1( 45),  y1( 46),  y1( 47),  y1( 48) &
      /-0.09e0,  1.23e0,  1.46e0,  0.61e0,  0.42e0,  2.16e0/
  data   y1( 49),  y1( 50),  y1( 51),  y1( 52),  y1( 53),  y1( 54) &
      / 3.18e0,  2.10e0,  0.37e0, -0.24e0,  0.57e0, -0.53e0/
  data   y1( 55),  y1( 56),  y1( 57),  y1( 58),  y1( 59),  y1( 60) &
      / 2.44e0,  1.02e0, -0.53e0, -2.49e0, -2.12e0, -1.04e0/
  data   y1( 61),  y1( 62),  y1( 63),  y1( 64),  y1( 65),  y1( 66) &
      /-0.12e0, -1.88e0, -1.50e0,  1.54e0,  3.33e0,  3.08e0/
  data   y1( 67),  y1( 68),  y1( 69),  y1( 70),  y1( 71),  y1( 72) &
      / 1.71e0,  0.79e0,  1.55e0,  0.89e0, -0.89e0, -1.18e0/
  data   y1( 73),  y1( 74),  y1( 75),  y1( 76),  y1( 77),  y1( 78) &
      / 0.89e0,  1.71e0,  3.05e0,  0.15e0, -1.04e0,  0.12e0/
  data   y1( 79),  y1( 80),  y1( 81),  y1( 82),  y1( 83),  y1( 84) &
      / 0.08e0,  0.11e0, -2.62e0, -1.28e0,  1.07e0,  3.20e0/
  data   y1( 85),  y1( 86),  y1( 87),  y1( 88),  y1( 89),  y1( 90) &
      / 1.92e0,  0.53e0, -1.08e0,  0.49e0, -0.58e0,  0.17e0/
  data   y1( 91),  y1( 92),  y1( 93),  y1( 94),  y1( 95),  y1( 96) &
      / 1.15e0, -0.97e0, -1.63e0,  1.14e0, -0.67e0, -0.88e0/
  data   y1( 97),  y1( 98),  y1( 99),  y1(100) &
      /-0.07e0,  0.24e0,  0.55e0, -2.16e0/
  data   y2(  1),  y2(  2),  y2(  3),  y2(  4),  y2(  5),  y2(  6) &
      / 0.79e0,  1.12e0, -1.10e0, -2.39e0, -1.75e0, -0.82e0/
  data   y2(  7),  y2(  8),  y2(  9),  y2( 10),  y2( 11),  y2( 12) &
      /-0.36e0,  1.27e0,  1.75e0,  2.44e0,  0.36e0, -2.10e0/
  data   y2( 13),  y2( 14),  y2( 15),  y2( 16),  y2( 17),  y2( 18) &
      /-1.93e0, -1.30e0, -1.75e0, -0.34e0,  0.74e0,  0.49e0/
  data   y2( 19),  y2( 20),  y2( 21),  y2( 22),  y2( 23),  y2( 24) &
      / 0.70e0,  0.71e0,  0.09e0,  0.59e0,  1.54e0,  0.14e0/
  data   y2( 25),  y2( 26),  y2( 27),  y2( 28),  y2( 29),  y2( 30) &
      / 0.55e0, -1.40e0, -2.55e0, -1.66e0, -0.43e0,  0.58e0/
  data   y2( 31),  y2( 32),  y2( 33),  y2( 34),  y2( 35),  y2( 36) &
      / 2.18e0, -0.24e0,  0.58e0, -0.18e0, -1.55e0, -0.64e0/
  data   y2( 37),  y2( 38),  y2( 39),  y2( 40),  y2( 41),  y2( 42) &
      /-1.09e0,  0.90e0, -0.66e0, -0.35e0,  0.48e0,  0.50e0/
  data   y2( 43),  y2( 44),  y2( 45),  y2( 46),  y2( 47),  y2( 48) &
      / 0.05e0, -0.68e0,  0.24e0,  0.58e0, -1.26e0, -0.25e0/
  data   y2( 49),  y2( 50),  y2( 51),  y2( 52),  y2( 53),  y2( 54) &
      / 0.25e0,  2.18e0,  2.96e0,  1.56e0, -0.36e0, -0.59e0/
  data   y2( 55),  y2( 56),  y2( 57),  y2( 58),  y2( 59),  y2( 60) &
      /-0.12e0,  3.03e0,  2.11e0,  0.78e0,  0.89e0, -1.45e0/
  data   y2( 61),  y2( 62),  y2( 63),  y2( 64),  y2( 65),  y2( 66) &
      /-0.36e0, -0.37e0, -1.39e0, -4.19e0, -0.73e0, -0.98e0/
  data   y2( 67),  y2( 68),  y2( 69),  y2( 70),  y2( 71),  y2( 72) &
      / 0.36e0,  0.06e0, -1.94e0, -0.08e0,  0.17e0,  1.00e0/
  data   y2( 73),  y2( 74),  y2( 75),  y2( 76),  y2( 77),  y2( 78) &
      /-0.05e0,  0.43e0,  0.15e0,  2.69e0,  0.57e0,  0.29e0/
  data   y2( 79),  y2( 80),  y2( 81),  y2( 82),  y2( 83),  y2( 84) &
      / 1.10e0,  0.48e0, -1.06e0, -2.28e0, -2.03e0, -0.75e0/
  data   y2( 85),  y2( 86),  y2( 87),  y2( 88),  y2( 89),  y2( 90) &
      / 1.00e0,  1.71e0,  0.58e0,  1.97e0,  0.99e0,  1.94e0/
  data   y2( 91),  y2( 92),  y2( 93),  y2( 94),  y2( 95),  y2( 96) &
      / 2.18e0,  3.14e0,  0.60e0,  0.51e0,  1.35e0,  0.56e0/
  data   y2( 97),  y2( 98),  y2( 99),  y2(100) &
      / 0.11e0,  0.00e0,  2.34e0,  1.88e0/

  itest = 1
  ldstak = lds

  n = 100
  lagmax = 20
  nlag = 30
  nprt = 1
  lyfft = 150
  iccov = 30
  jccov = 5
  iym = 150
  m = 4
  iymfft = 150
  inlppc = 30
  jnlppc = 5
  nyd = 144
  ymiss0 = 1.16e0
  y1miss = 0.89e0
  y2miss = 0.89e0
!
!  copy data into ym for ccfs and ccfms
!
  call scopy(n, y1, 1, ym(1,1), 1)
  call scopy(n, y2, 1, ym(1,2), 1)
  call scopy(n, y1, 1, ym(1,3), 1)
  call scopy(n, y2, 1, ym(1,4), 1)
  call setrv(ymmiss, 4, ymiss0)
!
!  test of ccf
!
  write ( *,1060)
  call ccf(y1, y2, n)
!
!  print returned results
!
  call ccfxp (.false., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test of ccfs
!
  write ( *,1080)
  call ccfs(ym, n, m, iym, lagmax, ccov, iccov, jccov, nprt, &
     ldstak)
!
!  print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test of ccfm without missing values
!
  write ( *,1070)
  write ( *, 1050)
  call ccfm(y1, ymiss0, y2, ymiss0, n)
!
!  print returned results
!
  call ccfxp ( .false., lagmax, m, ccov, iccov, jccov, .true., &
     nlppc,  inlppc, jnlppc, cmiss )
!
!  test of ccfms without missing values
!
  write ( *,1140)
  write ( *, 1050)
  call ccfms ( ym, ymmiss, n, m, iym, lagmax, ccov, cmiss, &
     iccov, jccov, nlppc, inlppc, jnlppc, nprt, ldstak )
!
!  print returned results
!
  call ccfxp ( .true., lagmax, m, ccov, iccov, jccov, .true., &
     nlppc,  inlppc, jnlppc, cmiss )
!
!  copy data into yfft1, yfft2 and ymfft for ccff and ccffs
!
  call scopy(n, y1, 1, yfft1, 1)
  call scopy(n, y2, 1, yfft2, 1)
  call scopy(n, y1, 1, ymfft(1,1), 1)
  call scopy(n, y2, 1, ymfft(1,2), 1)
  call scopy(n, y1, 1, ymfft(1,3), 1)
  call scopy(n, y2, 1, ymfft(1,4), 1)
!
!  test of ccff
!
  write ( *,1100)
  call ccff(yfft1, yfft2, n, lyfft, ldstak)
!
!  print returned results
!
  call ccfxp (.false., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test of ccffs
!
  write ( *,1150)
  call ccffs(ymfft, n, m, iymfft, lagmax, ccov, &
     iccov, jccov, nprt, ldstak)
!
!  print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  reset ymmiss
!
  ymmiss(1) = y1miss
  ymmiss(2) = y2miss
  ymmiss(3) = y1miss
  ymmiss(4) = y2miss
!
!  test of ccfm with missing values
!
  write ( *,1070)
  write ( *, 1040)
  call ccfm(y1, y1miss, y2, y2miss, n)
!
!  print returned results
!
  call ccfxp (.false., lagmax, m, ccov, iccov, jccov, .true., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test of ccfms with missing values
!
  write ( *,1140)
  write ( *, 1040)
  call ccfms(ym, ymmiss, n, m, iym, lagmax, ccov, cmiss, &
     iccov, jccov, nlppc, inlppc, jnlppc, nprt, ldstak)
!
!     print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .true., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test print control
!
  nprt = 0
!
!  test of ccfs
!
  write ( *,1080)
  write ( *, 1020)
  call ccfs(ym, n, m, lagmax, iym, ccov, iccov, jccov, nprt, &
     ldstak)
!
!  print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test of ccfms with missing values
!
  write ( *,1140)
  write ( *, 1040)
  write ( *, 1020)
  call ccfms(ym, ymmiss, n, m, iym, lagmax, ccov, cmiss, &
     iccov, jccov, nlppc, inlppc, jnlppc, nprt, ldstak)
!
!  print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .true., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  copy data into ymfft for ccffs
!
  call scopy(n, y1, 1, ymfft(1,1), 1)
  call scopy(n, y2, 1, ymfft(1,2), 1)
  call scopy(n, y1, 1, ymfft(1,3), 1)
  call scopy(n, y2, 1, ymfft(1,4), 1)
!
!  test of ccffs
!
  write ( *,1150)
  write ( *, 1020)
  call ccffs(ymfft, n, m, iymfft, lagmax, ccov, &
     iccov, jccov, nprt, ldstak)
!
!  print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test lead/lag message
!
  nprt = 1

  call setra(ymfft, iymfft, m, n, 0.0e0)
  ymfft(5,1) = 1.0e0
  ymfft(15,2) = 1.0e0
  ymfft(5,3) = ymfft(5,1)
  ymfft(15,4) = ymfft(15,2)
!
!  test of ccffs
!
  write ( *,1150)
  write ( *, 1020)
  call ccffs(ymfft, n, m, iymfft, lagmax, ccov, &
     iccov, jccov, nprt, ldstak)
!
!  print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
  go to (100, 200, 300, 400), itest
!
!     test minimum problem size
!
  100 itest = itest + 1
  n = 3
  lagmax = 1
  lyfft = 150
  iccov = 30
  jccov = 5
  iym = 150
  m = 1
  iymfft = 150
  inlppc = 30
  jnlppc = 5
  nyd = 144
  ymiss0 = 1.16e0
  y1miss = 0.89e0
  y2miss = 0.89e0
!
!  test error handling
!
  200 itest = itest + 1
  n = 0
  lagmax = 1
  lyfft = 0
  iccov = 0
  jccov = 0
  iym = 0
  m = 0
  iymfft = 0
  inlppc = 0
  jnlppc = 0
  nyd = 0
!
!  test error handling
!
  300 itest = itest + 1
  n = 100
  lagmax = 100
  lyfft = 0
  iccov = 0
  jccov = 0
  iym = 0
  m = 0
  iymfft = 0
  inlppc = 0
  jnlppc = 0
  nyd = 144
  ldstak = 0

  400 return

 1020 format (' output suppressed')
 1040 format (' with missing values')
 1050 format (' without missing values')
 1060 format ('test of ccf')
 1070 format ('test of ccfm')
 1080 format ('test of ccfs')
 1100 format ('test of ccff')
 1140 format ('test of ccfms')
 1150 format ('test of ccffs')
end
subroutine xcorr ( ldstak )

!*****************************************************************************80
!
!! XCORR exercises all aspects of the correlation family routines.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Linda Mitchell,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    Draper, Smith,
!    Applied Regression Analysis
!    page 216
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!         an index variable.
!     integer ierr
!        common flag indicating if any errors were detected
!        if ierr = 0, then no errors were found
!     integer ivcv
!        the row dimension of vcv
!     integer iym
!        the row dimension of ym
!     integer j
!        an index variable.
!     integer ldsmin
!        the smallest acceptable size of common area cstak
!     integer ldstak
!        the size of the common area cstak
!     integer m
!        the number of variables
!     integer n
!        the number of observations
!     real vcv(4,4)
!        the variance covariance matrix
!     real ym(10,4)
!        general data set, from draper and smith
!     real z(10,4)
!        test observation matrix
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  integer &
     i,ivcv,iym,j,ldsmin,m,n
!
!  local arrays
  real &
     vcv(4,4),ym(10,4),z(10,4)
!
!  external subroutines
  external corr,corrs,corrxp,genr,iprint,ldscmp,msgx,setra
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data     ym(1,1),   ym(1,2),   ym(1,3),   ym(1,4) &
      /      42.2e0,  11.2e0,  31.9e0, 167.1e0/
  data     ym(2,1),   ym(2,2),   ym(2,3),   ym(2,4) &
      /      48.6e0,  10.6e0,  13.2e0, 174.4e0/
  data     ym(3,1),   ym(3,2),   ym(3,3),   ym(3,4) &
      /      42.6e0,  10.6e0,  28.7e0, 160.8e0/
  data     ym(4,1),   ym(4,2),   ym(4,3),   ym(4,4) &
      /      39.0e0,  10.4e0,  26.1e0, 162.0e0/
  data     ym(5,1),   ym(5,2),   ym(5,3),   ym(5,4) &
      /      34.7e0,   9.3e0,  30.1e0, 140.8e0/
  data     ym(6,1),   ym(6,2),   ym(6,3),   ym(6,4) &
      /      44.5e0,  10.8e0,   8.5e0, 174.6e0/
  data     ym(7,1),   ym(7,2),   ym(7,3),   ym(7,4) &
      /      39.1e0,  10.7e0,  24.3e0, 163.7e0/
  data     ym(8,1),   ym(8,2),   ym(8,3),   ym(8,4) &
      /      40.1e0,  10.0e0,  18.6e0, 174.5e0/
  data     ym(9,1),   ym(9,2),   ym(9,3),   ym(9,4) &
      /      45.9e0,  12.0e0,  20.4e0, 185.7e0/

  ivcv = 4
  iym = 10
  m = 4
  n = 9
  ierr = 0
!
!  test routines with correct call statement.
!
  write ( *,1000)
  write ( *,1010)
!
!  test corr
!
  write ( *,1020)
  write ( *,1060)
  call corr(ym, n, m, iym, ldstak)
  call msgx(0 )
!
!  test corrs
!
!  printout suppressed
!
  write ( *,1030)
  write ( *,1040)
  write ( *,1060)
  call corrs(ym, n, m, iym, ldstak, 0, vcv, ivcv)
  call msgx ( 0 )
!
!  print stored output and zero arrays
!
  call corrxp ( m, vcv, ivcv )
!
!  with printout
!
  write ( *,1050)
  write ( *,1060)
  call corrs(ym, n, m, iym, ldstak, 1, vcv, ivcv)
  call msgx ( 0 )
!
!  print stored output
!
  call corrxp ( m, vcv, ivcv )
!
!  special 2 column matrix.
!
  write ( *,1070)
  write ( *,1060)
  call corr(ym, n, 2, iym, ldstak)
  call msgx ( 0 )
!
!  test work area requirements.
!
!  test corr
!
  call ldscmp(12, 0, max(n,m), 0, 0, 0, 's', &
     m*m + (max(n,m)+m+n*(m+3)+6*m*m), ldsmin)
  write ( *,1090)
  call corr(ym, n, m, iym, ldsmin-1)
  call msgx ( 1 )
  write ( *,1100)
  call corr(ym, n, m, iym, ldsmin)
  call msgx ( 0 )
!
!  test corrs with printout
!
  call ldscmp(12, 0, max(n,m), 0, 0, 0, 's', &
     max(n,m)+m+n*(m+3)+6*m*m, ldsmin)
  write ( *,1090)
  call corrs(ym, n, m, iym, ldsmin-1, 1, vcv, ivcv)
  call msgx ( 1 )
  write ( *,1100)
  call corrs(ym, n, m, iym, ldsmin, 1, vcv, ivcv)
  call corrxp(m, vcv, ivcv )
  call msgx ( 0 )
!
!  test corrs without printout
!
  call ldscmp(12, 0, 0, 0, 0, 0, 's', 0, ldsmin)
  write ( *,1090)
  call corrs(ym, n, m, iym, ldsmin-1, 0, vcv, ivcv)
  call msgx ( 1 )
  write ( *,1100)
  call corrs(ym, n, m, iym, ldsmin, 0, vcv, ivcv)
  call corrxp(m, vcv, ivcv )
  call msgx ( 0 )
!
!  number of variables less than 2.
!
  write ( *,1110)
!
!  test corr
!
  call corr(ym, n, 1, iym, ldstak)
  call msgx ( 1 )
!
!  test corrs
!
  call corrs(ym, n, 1, iym, ldstak, 1, vcv, ivcv)
  call msgx ( 1 )
!
!  number of observations less than 3.
!
  write ( *,1120)
!
!  test corr
!
  call corr(ym, 2, 4, iym, ldstak)
  call msgx ( 1 )
!
!  test corrs
!
  call corrs(ym, 2, 4, iym, ldstak, 1, vcv, ivcv)
  call msgx ( 1 )
!
!  observation matrix dimensioned less than n.
!
  write ( *,1150)
!
!  test corr
!
  call corr(ym, n, m, 8, ldstak)
  call msgx ( 1 )
!
!  test corrs
!
  call corrs(ym, n, m, 8, ldstak, 1, vcv, ivcv)
  call msgx ( 1 )
!
!  vcv matrix dimensioned less than m.
!
  write ( *,1130)
  call corrs(ym, n, m, iym, ldstak, 1, vcv, 2)
  call msgx ( 1 )
!
!  all observations on a single variable equal to zero.
!
  write ( *,1140)
  call setra(z, 10, 4, 10, 0.0e0)
  call corr(z, 9, 4, 10, ldstak)
  call msgx ( 1 )
  call corrs(z, 9, 4, 10, ldstak, 1, vcv, ivcv)
  call corrxp(m, vcv, ivcv )
  call msgx ( 1 )

  z(1:10,1) = i
  z(1:10,2) = 0.0e0
  call corr(z, 10, 4, 10, ldstak)
  call msgx ( 1 )
!
!  array filled with a single value.
!
  write ( *,1160)
  call setra(z, 10, 4, 10, 4.0e0)
  call corr(z, 4, 10, 4, ldstak)
  call msgx ( 1 )
!
!  2 columns the same.
!
  do i=1,3
     call genr(z(1,i), 5, 5.0e0*i, 5.0e0*i)
  end do
  z(1:5,4) = z(1:5,3)

  write ( *,1170)
  call corr(z, 5, 4, 10, ldstak)
  call msgx ( 1 )
!
!  2 columns inversely related.
!
  j = 5
  do i=1,5
     j = j - 1
     z(j,4) = z(i,3)
  end do
  write ( *,1170)
  call corr(z, 5, 4, 10, ldstak)
  call msgx ( 1 )

  return

 1000 format('1')
 1010 format(' ****test routines with correct call****')
 1020 format(' test of corr')
 1030 format('1test of corrs')
 1040 format(' printout supressed.')
 1050 format('1printout not supressed.')
 1060 format(' draper and smith data set (page 216).')
 1070 format('1****special case 2 column matrix****')
 1090 format('1****test with insufficient work area****')
 1100 format('1****test with exactly the right amount of work area****')
 1110 format('1****number of variables less than 2****')
 1120 format(' ****number of observations less than 3****')
 1130 format(' ****inadequate space in storage arrays****')
 1140 format('1****all observations on a variable equal to zero****')
 1150 format(' ****observation matrix dimensioned less than number', &
         ' of observations designated****')
 1160 format('1****array containing a single value****')
 1170 format('1****2 columns related****')
end
subroutine xdckld ( ldstak )

!*****************************************************************************80
!
!! XDCKLD tests the derivative checking routines.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     real delta
!        *
!     external drv4a
!        the name of the user supplied subroutine which computes the
!        analytic derivatives (jacobian matrix) of the model.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index variable.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ixm
!        the first dimension of the independent variable array xm.
!     integer ldsmin
!        the minimum length of the array dstak allowed.
!     integer ldstak
!        the length of the array dstak.
!     integer m
!        the number of independent variables.
!     external mdl4
!        the name of the user supplied subroutine which computes the
!        predicted values based on the current parameter estimates.
!     integer n
!        the number of observations of data.
!     integer neta
!        the number of reliable digits in the model.
!     integer npar
!        the number of unknown parameters in the model.
!     integer nprt
!        the indicator variable used to specify whether or not
!        printed output is to be provided, where if the value of
!        nprt is zero, no printed output is given.
!     integer nrow
!        the number of the row of the independent variable array at
!        which the derivative is to be checked.
!     integer ntau
!        the number of digits of agreement required between the
!        numerical derivatives and the user supplied derivatives.
!     integer ntest
!        the number of the current test.
!     real par(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real scale(10)
!        a dummy array, indicating use of default values for
!        the typical size of the parameters.
!     real xm(200,2)
!        the array in which one row of the independent variable array
!        is stored.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     delta
  integer &
     i,ixm,ldsmin,m,n,neta,npar,nprt,nrow,ntau,ntest
!
!  local arrays
  real &
     par(10),scale(10),xm(200,2)
!
!  external subroutines
  external dckls,dckls1,dcklsc,drv4a,iprint,ldscmp,mdl4
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!  set parameter values
!
  call dckls1(n, m, ixm, par, npar, neta, ntau, nrow, scale, nprt)
  call ldscmp(5, 0, 2*npar+1, 0, 0, 0, 's', &
              n*npar+npar+n, ldsmin)

  if (ldsmin.gt.ldstak) then
     write ( *, 1020) ldsmin
     return
  end if
!
!  create independent variable
!
  delta = 0.0625e0
  xm(1,1) = 0.0e0
  do i=2,n
     xm(i,1) = xm(i-1,1) + delta
  end do

  ntest = 0
!
!  check results from valid calls
!
!  simple example
!
!  check result for correctly computed derivative
!
  call dckls1(n, m, ixm, par, npar, neta, ntau, nrow, scale, nprt)

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1100)
  write ( *,1040)
  write ( *,1000)
  ierr = -1
  call dckls(xm, n, m, ixm, mdl4, drv4a, par, npar, ldsmin)
  write ( *,1050) ierr

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1100)
  write ( *,1040)
  write ( *,1060) neta, ntau, scale(1), nrow, nprt
  write ( *,1010)
  ierr = -1
  call dcklsc(xm, n, m, ixm, mdl4, drv4a, par, npar, ldsmin, neta, &
     ntau, scale, nrow, nprt)
  write ( *,1050) ierr
  write ( *,1140) neta, ntau, scale(1), nrow, nprt

  return

 1000 format (' test of dckls' )
 1010 format ('test of dcklsc')
 1020 format ('ldstak must be greater than or equal to ', i6)
 1040 format ('simple example')
 1050 format ('***** returned results *****', 5x, ' (-1 indicates ', &
     'value not changed by called subroutine)'//' ierr is ', i3)
 1060 format (' input   -  neta = ', i5, ', ntau = ', i5, &
     ', scale(1) = ', g15.8, ', nrow = ', i5, ', nprt = ', i5)
 1100 format (' correctly coded derivative')
 1130 format ('derivative checking subroutine test number', i5/)
 1140 format ('output  -  neta = ', i5, ', ntau = ', i5, &
     ', scale(1) = ', g15.8, ', nrow = ', i5, ', nprt = ', i5//)
end
subroutine xdckle ( ldstak )

!*****************************************************************************80
!
!! XDCKLE tests the derivative checking routines.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     real delta
!        *
!     external drv4a, drv4b
!        the name of the user supplied subroutine which computes the
!        analytic derivatives (jacobian matrix) of the model.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index variable.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ixm
!        the first dimension of the independent variable array xm.
!     integer ldsmin
!        the minimum length of the array dstak allowed.
!     integer ldstak
!        the length of the array dstak.
!     integer m
!        the number of independent variables.
!     external mdl4
!        the name of the user supplied subroutine which computes the
!        predicted values based on the current parameter estimates.
!     integer n
!        the number of observations of data.
!     integer neta
!        the number of reliable digits in the model.
!     integer npar
!        the number of unknown parameters in the model.
!     integer nprt
!        the indicator variable used to specify whether or not
!        printed output is to be provided, where if the value of
!        nprt is zero, no printed output is given.
!     integer nrow
!        the number of the row of the independent variable array at
!        which the derivative is to be checked.
!     integer ntau
!        the number of digits of agreement required between the
!        numerical derivatives and the user supplied derivatives.
!     integer ntest
!        the number of the current test.
!     real par(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real scale(10)
!        a dummy array, indicating use of default values for
!        the typical size of the parameters.
!     real xm(200,2)
!        the array in which one row of the independent variable array
!        is stored.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     delta
  integer &
     i,ixm,ldsmin,m,n,neta,npar,nprt,nrow,ntau,ntest
!
!  local arrays
  real &
     par(10),scale(10),xm(200,2)
!
!  external subroutines
  external dckls,dckls1,dcklsc,drv4a,drv4b,iprint,ldscmp,mdl4
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!     set parameter values
!
  call dckls1(n, m, ixm, par, npar, neta, ntau, nrow, scale, nprt)
  call ldscmp(5, 0, 2*npar+1, 0, 0, 0, 's', &
              n*npar+npar+n, ldsmin)

  if (ldsmin.le.ldstak) go to 5

  write ( *, 1040) ldsmin
  return

    5 continue
!
!  create independent variable
!
  delta = 0.0625e0
  xm(1,1) = 0.0e0
  do i=2,n
     xm(i,1) = xm(i-1,1) + delta
  end do

  ntest = 0
!
!  check error handling
!
!  test 1  -  miscellaneous error checking
!
  n = -5
  m = -5
  ixm = -10
  npar = -10

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1020)
  write ( *,1000)
  ierr = -1
  call dckls(xm, n, m, ixm, mdl4, drv4a, par, npar, ldstak)
  write ( *,1050) ierr
  write ( *,1010)
  ierr = -1
  call dcklsc(xm, n, m, ixm, mdl4, drv4a, par, npar, ldstak, neta, &
     ntau, scale, nrow, nprt)
  write ( *,1050) ierr
!
!  test 2  -  miscellaneous error checking (continued)
!
  call dckls1(n, m, ixm, par, npar, neta, ntau, nrow, scale, nprt)
  scale(2) = 0.0e0

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1030)
  write ( *,1000)
  ierr = -1
  call dckls(xm, n, m, ixm, mdl4, drv4a, par, npar, ldsmin-1)
  write ( *,1050) ierr
  write ( *,1010)
  ierr = -1
  call dcklsc(xm, n, m, ixm, mdl4, drv4a, par, npar, ldsmin-1, &
     neta, ntau, scale, nrow, nprt)
  write ( *,1050) ierr

  return

 1000 format (' test of dckls' )
 1010 format (' test of dcklsc')
 1020 format (' check error handling  -  test 1')
 1030 format (' check error handling  -  test 2')
 1040 format ('1 *** ldstak must be greater than or equal to ', i6)
 1050 format (' ***** returned results *****', 5x, &
    '(-1 indicates value not changed by called subroutine)'//' ierr is ', i3)
 1130 format ('1derivative checking subroutine test number', i5)
end
subroutine xdcklt ( ldstak )

!*****************************************************************************80
!
!! XDCKLT tests the derivative checking routines.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     real delta
!        *
!     external drv4a, drv4b
!        the name of the user supplied subroutine which computes the
!        analytic derivatives (jacobian matrix) of the model.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index variable.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ixm
!        the first dimension of the independent variable array xm.
!     integer j, jstop
!        *
!     integer ldsmin
!        the minimum length of the array dstak allowed.
!     integer ldstak
!        the length of the array dstak.
!     integer m
!        the number of independent variables.
!     external mdl4
!        the name of the user supplied subroutine which computes the
!        predicted values based on the current parameter estimates.
!     integer n
!        the number of observations of data.
!     integer neta
!        the number of reliable digits in the model.
!     integer nettst(6)
!        various test values for neta.
!     integer npar
!        the number of unknown parameters in the model.
!     integer nprt
!        the indicator variable used to specify whether or not
!        printed output is to be provided, where if the value of
!        nprt is zero, no printed output is given.
!     integer nrotst(5)
!        various test values for nrow.
!     integer nrow
!        the number of the row of the independent variable array at
!        which the derivative is to be checked.
!     integer ntatst(6)
!         various test values for ntau.
!     integer ntau
!        the number of digits of agreement required between the
!        numerical derivatives and the user supplied derivatives.
!     integer ntest
!        the number of the current test.
!     real par(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real scale(10)
!        a dummy array, indicating use of default values for
!        the typical size of the parameters.
!     real xm(200,2)
!        the array in which one row of the independent variable array
!        is stored.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     delta
  integer &
     i,ixm,j,jstop,ldsmin,m,n,neta,npar,nprt,nrow,ntau, &
     ntest
!
!  local arrays
  real &
     par(10),scale(10),xm(200,2)
  integer &
     nettst(6),nrotst(5),ntatst(6)
!
!  external functions
  real &
     r1mach
  external r1mach
!
!  external subroutines
  external dckls,dckls1,dcklsc,drv4a,drv4b,iprint,ldscmp,mdl4,setrv
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!     set parameter values
!
  call dckls1(n, m, ixm, par, npar, neta, ntau, nrow, scale, nprt)
  call ldscmp(5, 0, 2*npar+1, 0, 0, 0, 's', &
              n*npar+npar+n, ldsmin)

  if (ldsmin.le.ldstak) go to 5

  write ( *, 1020) ldsmin
  return

    5 continue
!
!  create independent variable
!
  delta = 0.0625e0
  xm(1,1) = 0.0e0
  do i=2,n
     xm(i,1) = xm(i-1,1) + delta
  end do

  ntest = 0
!
!  test various values of neta and ntau
!
  scale(1) = 0.0e0

  nettst(1) = -1
  nettst(2) = 0
  nettst(3) = 1
  nettst(4) = 2

  nettst(5) = -log10(r1mach(4))
  nettst(6) = nettst(5) + 1

  ntatst(1) = -1
  ntatst(2) = 0
  ntatst(3) = 1

  jstop = 3

  do i=1,6

     ntatst(4) = nettst(i)/4
     if (i.le.5) then
        ntatst(5) = (nettst(i)-1)/2
        ntatst(6) = ntatst(5) + 1
     end if

     if (i==5) jstop = 6

     do j=1,jstop

        ntest = ntest + 1
        write ( *,1130) ntest
        write ( *,1100)
        write ( *,1040)
        write ( *,1060) nettst(i), ntatst(j), scale(1), nrow, nprt
        write ( *,1000)
        ierr = -1
        call dcklsc(xm, n, m, ixm, mdl4, drv4a, par, npar, ldstak, &
           nettst(i), ntatst(j), scale, nrow, nprt)
        write ( *,1050) ierr
        write ( *,1140) nettst(i), ntatst(j), scale(1), nrow, nprt

      end do

   end do
!
!  suppress output
!
  nprt = 0

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1100)
  write ( *,1040)
  write ( *,1060) neta, ntau, scale(1), nrow, nprt
  write ( *,1010)
  ierr = -1
  call dcklsc(xm, n, m, ixm, mdl4, drv4a, par, npar, ldstak, neta, &
     ntau, scale, nrow, nprt)
  write ( *,1050) ierr
  write ( *,1140) neta, ntau, scale(1), nrow, nprt
!
!  large calculation error problem
!
  call dckls1(n, m, ixm, par, npar, neta, ntau, nrow, scale, nprt)
  par(3) = 10.0e0**ntatst(5)
  scale(1) = 0.0e0
  nrow = 51

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1100)
  write ( *,1070)
  write ( *,1080)
  write ( *,1060) neta, ntau, scale(1), nrow, nprt
  write ( *,1010)
  ierr = -1
  call dcklsc(xm, n, m, ixm, mdl4, drv4a, par, npar, ldstak, neta, &
     ntau, scale, nrow, nprt)
  write ( *,1050) ierr
  write ( *,1140) neta, ntau, scale(1), nrow, nprt
!
!  nearly zero derivative
!
  nrow = 50

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1100)
  write ( *,1070)
  write ( *,1090)
  write ( *,1060) neta, ntau, scale(1), nrow, nprt
  write ( *,1010)
  ierr = -1
  call dcklsc(xm, n, m, ixm, mdl4, drv4a, par, npar, ldstak, neta, &
     ntau, scale, nrow, nprt)
  write ( *,1050) ierr
  write ( *,1140) neta, ntau, scale(1), nrow, nprt
!
!  incorrectly coded derivative
!
!  simple example
!
!  set parameter values
!
  call dckls1(n, m, ixm, par, npar, neta, ntau, nrow, scale, nprt)

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1110)
  write ( *,1040)
  write ( *,1000)
  ierr = -1
  call dckls(xm, n, m, ixm, mdl4, drv4b, par, npar, ldstak)
  write ( *,1050) ierr

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1110)
  write ( *,1040)
  write ( *,1060) neta, ntau, scale(1), nrow, nprt
  write ( *,1010)
  ierr = -1
  call dcklsc(xm, n, m, ixm, mdl4, drv4b, par, npar, ldstak, neta, &
     ntau, scale, nrow, nprt)
  write ( *,1050) ierr
  write ( *,1140) neta, ntau, scale(1), nrow, nprt
!
!  suppress output
!
  nprt = 0

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1110)
  write ( *,1040)
  write ( *,1060) neta, ntau, scale(1), nrow, nprt
  write ( *,1010)
  ierr = -1
  call dcklsc(xm, n, m, ixm, mdl4, drv4b, par, npar, ldstak, neta, &
     ntau, scale, nrow, nprt)
  write ( *,1050) ierr
  write ( *,1140) neta, ntau, scale(1), nrow, nprt
!
!  large calculation error problem
!
  call dckls1(n, m, ixm, par, npar, neta, ntau, nrow, scale, nprt)

  par(3) = 10.0e0**ntatst(5)
  nrow = 26

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1110)
  write ( *,1070)
  write ( *,1060) neta, ntau, scale(1), nrow, nprt
  write ( *,1010)
  ierr = -1
  call dcklsc(xm, n, m, ixm, mdl4, drv4b, par, npar, ldstak, neta, &
     ntau, scale, nrow, nprt)
  write ( *,1050) ierr
  write ( *,1140) neta, ntau, scale(1), nrow, nprt

  par(4) = 0.75e0
  nrow = 1

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1110)
  write ( *,1070)
  write ( *,1060) neta, ntau, scale(1), nrow, nprt
  write ( *,1010)
  ierr = -1
  call dcklsc(xm, n, m, ixm, mdl4, drv4b, par, npar, ldstak, neta, &
     ntau, scale, nrow, nprt)
  write ( *,1050) ierr
  write ( *,1140) neta, ntau, scale(1), nrow, nprt
!
!  check various values of nrow
!
  call dckls1(n, m, ixm, par, npar, neta, ntau, nrow, scale, nprt)

  call setrv(xm(1,1), n, 0.0e0)
  nrotst(1) = -1
  nrotst(2) = 0
  nrotst(3) = 1
  nrotst(4) = n
  nrotst(5) = n + 1

  do i=1,5

     ntest = ntest + 1
     write ( *,1130) ntest
     write ( *,1110)
     write ( *,1120)
     write ( *,1060) neta, ntau, scale(1), nrotst(i), nprt
     write ( *,1010)
     ierr = -1
     call dcklsc(xm, n, m, ixm, mdl4, drv4b, par, npar, ldstak, &
        neta, ntau, scale, nrotst(i), nprt)
     write ( *,1050) ierr
     write ( *,1140) neta, ntau, scale(1), nrotst(i), nprt

  end do

  return

 1000 format ('test of dckls' )
 1010 format ('test of dcklsc')
 1020 format ('1 *** ldstak must be greater than or equal to ', i6)
 1040 format ('simple example')
 1050 format ('***** returned results *****', 5x, '(-1 indicates ', &
     'value not changed by called subroutine)'//' ierr is ', i3)
 1060 format ('input   -  neta = ', i5, ', ntau = ', i5, &
     ', scale(1) = ', g15.8, ', nrow = ', i5, ', nprt = ', i5)
 1070 format ('large calculation error problem')
 1080 format ('zero derivative')
 1090 format ('nearly zero derivative')
 1100 format ('correctly coded derivative')
 1110 format (' incorrectly coded derivative for parameters 1, 2 and 4')
 1120 format (' all independent variables equal to zero')
 1130 format ('derivative checking subroutine test number', i5)
 1140 format ('output  -  neta = ', i5, ', ntau = ', i5, &
     ', scale(1) = ', g15.8, ', nrow = ', i5, ', nprt = ', i5//)
end
subroutine xdemod ( lds )

!*****************************************************************************80
!
!! XDEMOD tests the time series complex demodulation routines.
!
!  Discussion:
!
!    Series Y is the Wolf sunspot data from 1700 to 1960 as
!    tabulated by Waldmeier.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    Max Waldmeier,
!    The Sunspot-Activity in the Years 1610-1960,
!    Shulthess, Zurich, 1961.
!
!  Parameters:
!
!     real ampl(300)
!        the array in which the amplitudes are stored.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fc
!        the cutoff frequency used for the low pass filter.
!     real fd
!        the demodulation frequency.
!     integer i
!        an indexing variable.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors have been detected.
!     integer itest
!        the number of the test being run
!     integer k
!        the number of terms in the symetric linear filter.
!     integer lds, ldstak
!        the length of the array dstak.
!     integer n
!        the number of observations in the input series.
!     integer ndem
!        the number of values in the demodulated series, i. e., it
!        is the number of values in the amplitude and phase arrays.
!     integer nprt
!        a code used to specify the type of plot, where if
!        nprt == 0 the plot is suppressed
!        nprt .ne. 1 the plot is provided
!     real phas(300)
!        the array in which the primary phase estimates are returned.
!     real y(300)
!        the vector containing the observed time series.
!
  implicit none

  integer &
     lds
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     fc,fd
  integer &
     i,itest,k,ldstak,n,ndem,nprt
!
!  local arrays
  real &
     ampl(300),phas(300),y(300)
!
!  external subroutines
  external demod,demods,iprint
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data   y(  1),  y(  2),  y(  3),  y(  4),  y(  5),  y(  6) &
      /     5.0e0, 11.0e0, 16.0e0, 23.0e0, 36.0e0, 58.0e0/
  data   y(  7),  y(  8),  y(  9),  y( 10),  y( 11),  y( 12) &
      /    29.0e0, 20.0e0, 10.0e0,  8.0e0,  3.0e0,  0.0e0/
  data   y( 13),  y( 14),  y( 15),  y( 16),  y( 17),  y( 18) &
      /     0.0e0, 2.0e0, 11.0e0, 27.0e0, 47.0e0, 63.0e0/
  data   y( 19),  y( 20),  y( 21),  y( 22),  y( 23),  y( 24) &
      /    60.0e0, 39.0e0, 28.0e0, 26.0e0, 22.0e0, 11.0e0/
  data   y( 25),  y( 26),  y( 27),  y( 28),  y( 29),  y( 30) &
      /    21.0e0, 40.0e0, 78.0e0,122.0e0,103.0e0, 73.0e0/
  data   y( 31),  y( 32),  y( 33),  y( 34),  y( 35),  y( 36) &
      /    47.0e0, 35.0e0, 11.0e0,  5.0e0, 16.0e0, 34.0e0/
  data   y( 37),  y( 38),  y( 39),  y( 40),  y( 41),  y( 42) &
      /    70.0e0, 81.0e0,111.0e0,101.0e0, 73.0e0, 40.0e0/
  data   y( 43),  y( 44),  y( 45),  y( 46),  y( 47),  y( 48) &
      /    20.0e0, 16.0e0,  5.0e0, 11.0e0, 22.0e0, 40.0e0/
  data   y( 49),  y( 50),  y( 51),  y( 52),  y( 53),  y( 54) &
      /    60.0e0, 80.9e0, 83.4e0, 47.7e0, 47.8e0, 30.7e0/
  data   y( 55),  y( 56),  y( 57),  y( 58),  y( 59),  y( 60) &
      /    12.2e0,  9.6e0, 10.2e0, 32.4e0, 47.6e0, 54.0e0/
  data   y( 61),  y( 62),  y( 63),  y( 64),  y( 65),  y( 66) &
      /    62.9e0, 85.9e0, 61.2e0, 45.1e0, 36.4e0, 20.9e0/
  data   y( 67),  y( 68),  y( 69),  y( 70),  y( 71),  y( 72) &
      /    11.4e0, 37.8e0, 69.8e0,106.1e0,100.8e0, 81.6e0/
  data   y( 73),  y( 74),  y( 75),  y( 76),  y( 77),  y( 78) &
      /    66.5e0, 34.8e0, 30.6e0,  7.0e0, 19.8e0, 92.5e0/
  data   y( 79),  y( 80),  y( 81),  y( 82),  y( 83),  y( 84) &
      /   154.4e0,125.9e0, 84.8e0, 68.1e0, 38.5e0, 22.8e0/
  data   y( 85),  y( 86),  y( 87),  y( 88),  y( 89),  y( 90) &
      /    10.2e0, 24.1e0, 82.9e0,132.0e0,130.9e0,118.1e0/
  data   y( 91),  y( 92),  y( 93),  y( 94),  y( 95),  y( 96) &
      /    89.9e0, 66.6e0, 60.0e0, 46.9e0, 41.0e0, 21.3e0/
  data   y( 97),  y( 98),  y( 99),  y(100),  y(101),  y(102) &
      /    16.0e0,  6.4e0,  4.1e0,  6.8e0, 14.5e0, 34.0e0/
  data   y(103),  y(104),  y(105),  y(106),  y(107),  y(108) &
      /    45.0e0, 43.1e0, 47.5e0, 42.2e0, 28.1e0, 10.1e0/
  data   y(109),  y(110),  y(111),  y(112),  y(113),  y(114) &
      /     8.1e0,  2.5e0,  0.0e0,  1.4e0,  5.0e0, 12.2e0/
  data   y(115),  y(116),  y(117),  y(118),  y(119),  y(120) &
      /    13.9e0, 35.4e0, 45.8e0, 41.1e0, 30.1e0, 23.9e0/
  data   y(121),  y(122),  y(123),  y(124),  y(125),  y(126) &
      /    15.6e0,  6.6e0,  4.0e0,  1.8e0,  8.5e0, 16.6e0/
  data   y(127),  y(128),  y(129),  y(130),  y(131),  y(132) &
      /    36.3e0, 49.6e0, 64.2e0, 67.0e0, 70.9e0, 47.8e0/
  data   y(133),  y(134),  y(135),  y(136),  y(137),  y(138) &
      /    27.5e0,  8.5e0, 13.2e0, 56.9e0,121.5e0,138.3e0/
  data   y(139),  y(140),  y(141),  y(142),  y(143),  y(144) &
      /   103.2e0, 85.7e0, 64.6e0, 36.7e0, 24.2e0, 10.7e0/
  data   y(145),  y(146),  y(147),  y(148),  y(149),  y(150) &
      /    15.0e0, 40.1e0, 61.5e0, 98.5e0,124.7e0, 96.3e0/
  data   y(151),  y(152),  y(153),  y(154),  y(155),  y(156) &
      /    66.6e0, 64.5e0, 54.1e0, 39.0e0, 20.6e0,  6.7e0/
  data   y(157),  y(158),  y(159),  y(160),  y(161),  y(162) &
      /     4.3e0, 22.7e0, 54.8e0, 93.8e0, 95.8e0, 77.2e0/
  data   y(163),  y(164),  y(165),  y(166),  y(167),  y(168) &
      /    59.1e0, 44.0e0, 47.0e0, 30.5e0, 16.3e0,  7.3e0/
  data   y(169),  y(170),  y(171),  y(172),  y(173),  y(174) &
      /    37.6e0, 74.0e0,139.0e0,111.2e0,101.6e0, 66.2e0/
  data   y(175),  y(176),  y(177),  y(178),  y(179),  y(180) &
      /    44.7e0, 17.0e0, 11.3e0, 12.4e0,  3.4e0,  6.0e0/
  data   y(181),  y(182),  y(183),  y(184),  y(185),  y(186) &
      /    32.3e0, 54.3e0, 59.7e0, 63.7e0, 63.5e0, 52.2e0/
  data   y(187),  y(188),  y(189),  y(190),  y(191),  y(192) &
      /    25.4e0, 13.1e0,  6.8e0,  6.3e0,  7.1e0, 35.6e0/
  data   y(193),  y(194),  y(195),  y(196),  y(197),  y(198) &
      /    73.0e0, 85.1e0, 78.0e0, 64.0e0, 41.8e0, 26.2e0/
  data   y(199),  y(200),  y(201),  y(202),  y(203),  y(204) &
      /    26.7e0, 12.1e0,  9.5e0,  2.7e0,  5.0e0, 24.4e0/
  data   y(205),  y(206),  y(207),  y(208),  y(209),  y(210) &
      /    42.0e0, 63.5e0, 53.8e0, 62.0e0, 48.5e0, 43.9e0/
  data   y(211),  y(212),  y(213),  y(214),  y(215),  y(216) &
      /    18.6e0,  5.7e0,  3.6e0,  1.4e0,  9.6e0, 47.4e0/
  data   y(217),  y(218),  y(219),  y(220),  y(221),  y(222) &
      /    57.1e0,103.9e0, 80.6e0, 63.6e0, 37.6e0, 26.1e0/
  data   y(223),  y(224),  y(225),  y(226),  y(227),  y(228) &
      /    14.2e0,  5.8e0, 16.7e0, 44.3e0, 63.9e0, 69.0e0/
  data   y(229),  y(230),  y(231),  y(232),  y(233),  y(234) &
      /    77.8e0, 64.9e0, 35.7e0, 21.2e0, 11.1e0,  5.7e0/
  data   y(235),  y(236),  y(237),  y(238),  y(239),  y(240) &
      /     8.7e0, 36.1e0, 79.7e0,114.4e0,109.6e0, 88.8e0/
  data   y(241),  y(242),  y(243),  y(244),  y(245),  y(246) &
      /    67.8e0, 47.5e0, 30.6e0, 16.3e0,  9.6e0, 33.2e0/
  data   y(247),  y(248),  y(249),  y(250),  y(251),  y(252) &
      /    92.6e0,151.6e0,136.3e0,134.7e0, 83.9e0, 69.4e0/
  data   y(253),  y(254),  y(255),  y(256),  y(257),  y(258) &
      /    31.5e0, 13.9e0,  4.4e0, 38.0e0,141.7e0,190.2e0/
  data   y(259),  y(260),  y(261) &
      /   184.8e0,159.0e0,112.3e0/

  itest = 1
  ldstak = lds

  n = 261
  nprt = 1
  fd = 1.0e0/11.0e0
  fc = 1.0e0/22.0e0
  k = 41
!
!  test of demod
!
5 continue

  write ( *, 1016)
  call demod (y, n, fd, fc, k, ldstak)
  write ( *, 1002) ierr
!
!  test of demods
!
  write ( *, 1017)
  call demods (y, n, fd, fc, k, ampl, phas, ndem, nprt, ldstak)
  write ( *, 1002) ierr
!
!  print storage from demods
!
  if (ierr==0) then
    write ( *, 1004) (ampl(i), i = 1, ndem)
    write ( *, 1004) (phas(i), i = 1, ndem)
  end if

  go to (100, 200, 300), itest
!
!  test minimum problem specifications
!
  100 continue

  itest = itest + 1
  n = 17
  k = 15
  nprt = -1
  go to 5
!
!  test error conditions
!
  200 continue

  itest = itest + 1
  n = 0
  fd = 0.5e0
  fc = 0.3e0
  k = 1
  go to 5

  300 return

 1002 format (' ierr is', i5)
 1004 format (10f10.5)
 1016 format ('test of demod')
 1017 format ('test of demods')
end
subroutine xdflt ( lds )

!*****************************************************************************80
!
!! XDFLT tests time series digital filtering and complex demodulation routines.
!
!  Discussion:
!
!    series y is the wolf sunspot data from 1700 to 1960 as
!    tabulated by waldmeier
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    Max Waldmeier,
!    The Sunspot-Activity in the Years 1610-1960,
!    Shulthess, Zurich, 1961.
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fc
!        the cutoff frequency used for the low pass filter.
!     real fmax, fmin
!        the minimum and maximum frequency for which the gain
!        function is to be estimated.
!     real freq(101)
!        the vector of frequencies at which the gain function
!        has been estimated.
!     real gain(101)
!        the vector in which the gain function estimates are
!        stored.
!     real hhp(50)
!        the array in which the -ideal- high pass filter coefficients
!        will be returned.
!     real hlp(50)
!        the array in which the input low pass filter coefficients
!        are stored.
!     integer i
!        an indexing variable.
!     integer iar
!        the number of filter coefficients.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors have been detected.
!     integer iod(10)
!        the order of each of the difference factors.
!     integer itest
!        the number of the test being performed
!     integer k
!        the number of terms in the symetric linear filter.
!     integer lds, ldstak
!        the length of the array dstak.
!     integer lphi
!        the length of the vector phi.
!     integer n
!        the number of observations in the input series.
!     integer nd(10)
!        the array containing the number of times the difference
!        factors are to be applied.
!     integer nfac
!        the number of difference factors.
!     integer nf
!        the number of frequencies at which the gain function
!        is to be estimated.
!     integer nprt
!        a code used to specify the type of plot, where if
!        nprt = 0 the plot is suppressed
!        nprt = 1 the plot is decibels/linear
!        nprt = 2 the plot is log/linear
!     integer nyf
!        the number of values in the filtered series.
!     integer nys
!        the number of values in the sampled series.
!     real phas(300)
!        the array in which the primary phase estimates are returned.
!     real phi(50)
!        the vector containing the filter coefficients.
!     real y(300)
!        the vector containing the observed time series.
!     real yf(300)
!        the vector in which the filtered series is returned.
!     real yfmiss
!        the missing value code used in the filtered series.
!     real ymiss
!        the missing value code used in the input time series.
!     real ys(300)
!        the array containing the sampled series.
!
  implicit none

  integer &
     lds
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     fc,fmax,fmin,yfmiss,ymiss
  integer &
     i,iar,itest,k,ldstak,lphi,n,nf,nfac,nprt,nyf,nys
!
!  local arrays
  real &
     freq(101),gain(101),hhp(50),hlp(50),phas(300),phi(50),y(300), &
     yf(300),ys(300)
  integer &
     iod(10),nd(10)
!
!  external subroutines
  external arflt,dif,difc,difm,difmc,gfarf,gfarfs,gfslf,gfslfs, &
     hipass,hpcoef,iprint,lopass,lpcoef,maflt,sample,slflt
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data   y(  1),  y(  2),  y(  3),  y(  4),  y(  5),  y(  6) &
      /     5.0e0, 11.0e0, 16.0e0, 23.0e0, 36.0e0, 58.0e0/
  data   y(  7),  y(  8),  y(  9),  y( 10),  y( 11),  y( 12) &
      /    29.0e0, 20.0e0, 10.0e0,  8.0e0,  3.0e0,  0.0e0/
  data   y( 13),  y( 14),  y( 15),  y( 16),  y( 17),  y( 18) &
      /     0.0e0, 2.0e0, 11.0e0, 27.0e0, 47.0e0, 63.0e0/
  data   y( 19),  y( 20),  y( 21),  y( 22),  y( 23),  y( 24) &
      /    60.0e0, 39.0e0, 28.0e0, 26.0e0, 22.0e0, 11.0e0/
  data   y( 25),  y( 26),  y( 27),  y( 28),  y( 29),  y( 30) &
      /    21.0e0, 40.0e0, 78.0e0,122.0e0,103.0e0, 73.0e0/
  data   y( 31),  y( 32),  y( 33),  y( 34),  y( 35),  y( 36) &
      /    47.0e0, 35.0e0, 11.0e0,  5.0e0, 16.0e0, 34.0e0/
  data   y( 37),  y( 38),  y( 39),  y( 40),  y( 41),  y( 42) &
      /    70.0e0, 81.0e0,111.0e0,101.0e0, 73.0e0, 40.0e0/
  data   y( 43),  y( 44),  y( 45),  y( 46),  y( 47),  y( 48) &
      /    20.0e0, 16.0e0,  5.0e0, 11.0e0, 22.0e0, 40.0e0/
  data   y( 49),  y( 50),  y( 51),  y( 52),  y( 53),  y( 54) &
      /    60.0e0, 80.9e0, 83.4e0, 47.7e0, 47.8e0, 30.7e0/
  data   y( 55),  y( 56),  y( 57),  y( 58),  y( 59),  y( 60) &
      /    12.2e0,  9.6e0, 10.2e0, 32.4e0, 47.6e0, 54.0e0/
  data   y( 61),  y( 62),  y( 63),  y( 64),  y( 65),  y( 66) &
      /    62.9e0, 85.9e0, 61.2e0, 45.1e0, 36.4e0, 20.9e0/
  data   y( 67),  y( 68),  y( 69),  y( 70),  y( 71),  y( 72) &
      /    11.4e0, 37.8e0, 69.8e0,106.1e0,100.8e0, 81.6e0/
  data   y( 73),  y( 74),  y( 75),  y( 76),  y( 77),  y( 78) &
      /    66.5e0, 34.8e0, 30.6e0,  7.0e0, 19.8e0, 92.5e0/
  data   y( 79),  y( 80),  y( 81),  y( 82),  y( 83),  y( 84) &
      /   154.4e0,125.9e0, 84.8e0, 68.1e0, 38.5e0, 22.8e0/
  data   y( 85),  y( 86),  y( 87),  y( 88),  y( 89),  y( 90) &
      /    10.2e0, 24.1e0, 82.9e0,132.0e0,130.9e0,118.1e0/
  data   y( 91),  y( 92),  y( 93),  y( 94),  y( 95),  y( 96) &
      /    89.9e0, 66.6e0, 60.0e0, 46.9e0, 41.0e0, 21.3e0/
  data   y( 97),  y( 98),  y( 99),  y(100),  y(101),  y(102) &
      /    16.0e0,  6.4e0,  4.1e0,  6.8e0, 14.5e0, 34.0e0/
  data   y(103),  y(104),  y(105),  y(106),  y(107),  y(108) &
      /    45.0e0, 43.1e0, 47.5e0, 42.2e0, 28.1e0, 10.1e0/
  data   y(109),  y(110),  y(111),  y(112),  y(113),  y(114) &
      /     8.1e0,  2.5e0,  0.0e0,  1.4e0,  5.0e0, 12.2e0/
  data   y(115),  y(116),  y(117),  y(118),  y(119),  y(120) &
      /    13.9e0, 35.4e0, 45.8e0, 41.1e0, 30.1e0, 23.9e0/
  data   y(121),  y(122),  y(123),  y(124),  y(125),  y(126) &
      /    15.6e0,  6.6e0,  4.0e0,  1.8e0,  8.5e0, 16.6e0/
  data   y(127),  y(128),  y(129),  y(130),  y(131),  y(132) &
      /    36.3e0, 49.6e0, 64.2e0, 67.0e0, 70.9e0, 47.8e0/
  data   y(133),  y(134),  y(135),  y(136),  y(137),  y(138) &
      /    27.5e0,  8.5e0, 13.2e0, 56.9e0,121.5e0,138.3e0/
  data   y(139),  y(140),  y(141),  y(142),  y(143),  y(144) &
      /   103.2e0, 85.7e0, 64.6e0, 36.7e0, 24.2e0, 10.7e0/
  data   y(145),  y(146),  y(147),  y(148),  y(149),  y(150) &
      /    15.0e0, 40.1e0, 61.5e0, 98.5e0,124.7e0, 96.3e0/
  data   y(151),  y(152),  y(153),  y(154),  y(155),  y(156) &
      /    66.6e0, 64.5e0, 54.1e0, 39.0e0, 20.6e0,  6.7e0/
  data   y(157),  y(158),  y(159),  y(160),  y(161),  y(162) &
      /     4.3e0, 22.7e0, 54.8e0, 93.8e0, 95.8e0, 77.2e0/
  data   y(163),  y(164),  y(165),  y(166),  y(167),  y(168) &
      /    59.1e0, 44.0e0, 47.0e0, 30.5e0, 16.3e0,  7.3e0/
  data   y(169),  y(170),  y(171),  y(172),  y(173),  y(174) &
      /    37.6e0, 74.0e0,139.0e0,111.2e0,101.6e0, 66.2e0/
  data   y(175),  y(176),  y(177),  y(178),  y(179),  y(180) &
      /    44.7e0, 17.0e0, 11.3e0, 12.4e0,  3.4e0,  6.0e0/
  data   y(181),  y(182),  y(183),  y(184),  y(185),  y(186) &
      /    32.3e0, 54.3e0, 59.7e0, 63.7e0, 63.5e0, 52.2e0/
  data   y(187),  y(188),  y(189),  y(190),  y(191),  y(192) &
      /    25.4e0, 13.1e0,  6.8e0,  6.3e0,  7.1e0, 35.6e0/
  data   y(193),  y(194),  y(195),  y(196),  y(197),  y(198) &
      /    73.0e0, 85.1e0, 78.0e0, 64.0e0, 41.8e0, 26.2e0/
  data   y(199),  y(200),  y(201),  y(202),  y(203),  y(204) &
      /    26.7e0, 12.1e0,  9.5e0,  2.7e0,  5.0e0, 24.4e0/
  data   y(205),  y(206),  y(207),  y(208),  y(209),  y(210) &
      /    42.0e0, 63.5e0, 53.8e0, 62.0e0, 48.5e0, 43.9e0/
  data   y(211),  y(212),  y(213),  y(214),  y(215),  y(216) &
      /    18.6e0,  5.7e0,  3.6e0,  1.4e0,  9.6e0, 47.4e0/
  data   y(217),  y(218),  y(219),  y(220),  y(221),  y(222) &
      /    57.1e0,103.9e0, 80.6e0, 63.6e0, 37.6e0, 26.1e0/
  data   y(223),  y(224),  y(225),  y(226),  y(227),  y(228) &
      /    14.2e0,  5.8e0, 16.7e0, 44.3e0, 63.9e0, 69.0e0/
  data   y(229),  y(230),  y(231),  y(232),  y(233),  y(234) &
      /    77.8e0, 64.9e0, 35.7e0, 21.2e0, 11.1e0,  5.7e0/
  data   y(235),  y(236),  y(237),  y(238),  y(239),  y(240) &
      /     8.7e0, 36.1e0, 79.7e0,114.4e0,109.6e0, 88.8e0/
  data   y(241),  y(242),  y(243),  y(244),  y(245),  y(246) &
      /    67.8e0, 47.5e0, 30.6e0, 16.3e0,  9.6e0, 33.2e0/
  data   y(247),  y(248),  y(249),  y(250),  y(251),  y(252) &
      /    92.6e0,151.6e0,136.3e0,134.7e0, 83.9e0, 69.4e0/
  data   y(253),  y(254),  y(255),  y(256),  y(257),  y(258) &
      /    31.5e0, 13.9e0,  4.4e0, 38.0e0,141.7e0,190.2e0/
  data   y(259),  y(260),  y(261) &
      /   184.8e0,159.0e0,112.3e0/

  itest = 1
  ldstak = lds

  n = 261
  nprt = 2
  fc = 1.0e0/22.0e0
  nf = 101
  fmin = 0.0e0
  fmax = 0.2e0
  lphi = 50
  nfac = 1
  nd(1) = 1
  iod(1) = 1
  iar = 1
  phi(1) = 0.6e0
  k = 41
  ymiss = 11.0e0
!
!  test of lpcoef
!
   10 write ( *, 1001)
  call lpcoef (fc, k, hlp)
  write ( *, 1002) ierr
!
!  print storage from lpcoef
!
  if (ierr==0) write ( *, 1004) hlp(1:k)
!
!  test of lopass
!
  write ( *, 1007)
  call lopass (y, n, fc, k, hlp, yf, nyf)
  write ( *, 1002) ierr
!
!  print storage from lopass
!
  if (ierr==0) then
    write ( *, 1004) (hlp(i), i = 1, k)
    write ( *, 1004) yf(1:nyf)
  end if
!
!  test of hipass
!
  write ( *, 1008)
  call hipass (y, n, fc, k, hhp, yf, nyf)
  write ( *, 1002) ierr
!
!  print storage from hipass
!
  if (ierr==0) then
    write ( *, 1004) (hhp(i), i = 1, k)
    write ( *, 1004) yf(1:nyf)
  end if
!
!  test of hpcoef
!
   20 write ( *, 1003)
  call hpcoef (hlp, k, hhp)
  write ( *, 1002) ierr
!
!  print storage from hpcoef
!
  if (ierr==0) write ( *, 1004) (hhp(i), i = 1, k)
!
!  test of maflt
!
  write ( *, 1020)
  call maflt (y, n, k, yf, nyf)
  write ( *, 1002) ierr
!
!  print storage from maflt
!
  if (ierr==0) write ( *, 1004) yf(1:nyf)
!
!  test of slflt
!
  write ( *, 1005)
  call slflt (y, n, k, hlp, yf, nyf)
  write ( *, 1002) ierr
!
!  print storage from slflt
!
  if (ierr==0) write ( *, 1004) yf(1:nyf)
!
!  test of sample
!
  write ( *, 1006)
  call sample (yf, n, k, ys, nys)
  write ( *, 1002) ierr
!
!  print storage from sample
!
  if (ierr==0) write ( *, 1004) yf(1:nys)
!
!  test of arflt
!
  write ( *, 1009)
  call arflt (y, n,  iar, phi, yf, nyf)
  write ( *, 1002) ierr
!
!  print storage from arflt
!
  if (ierr==0) write ( *, 1004) yf(1:nyf)
!
!  test of dif
!
  write ( *, 1015)
  call dif (y, n, yf, nyf)
  write ( *, 1002) ierr
!
!  print storage from dif
!
  if (ierr==0) write ( *, 1004) yf(1:nyf)
!
!  test of difm
!
  write ( *, 1018)
  call difm (y, ymiss, n, yf, yfmiss, nyf)
  write ( *, 1002) ierr
!
!  print storage from difm
!
  if (ierr==0) then
     write ( *, 1004) yf(1:nyf)
     write ( *, 1004) yfmiss
  end if
!
!  test of gfslf
!
  write ( *, 1011)
  call gfslf (hlp, k)
  write ( *, 1002) ierr
!
!  test of gfarf
!
  write ( *, 1013)
  call gfarf (phi, iar)
  write ( *, 1002) ierr
!
!  test of difc
!
   30 write ( *, 1010)
  call difc (y, n, nfac, nd, iod, iar, phi, lphi, yf, nyf, ldstak)
  write ( *, 1002) ierr
!
!  print storage from difc
!
  if (ierr==0) then
    write ( *, 1004) (phi(i), i = 1, k)
    write ( *, 1004) yf(1:nyf)
  end if
!
!  test of difmc
!
  write ( *, 1019)
  call difmc (y, ymiss, n, nfac, nd, iod, iar, phi, lphi, yf, &
     yfmiss, nyf, ldstak)
  write ( *, 1002) ierr
!
!  print storage from difmc
!
  if (ierr==0) then
    write ( *, 1004) (phi(i), i = 1, k)
    write ( *, 1004) yf(1:nyf)
    write ( *, 1004) yfmiss
  end if
!
!  test of gfslfs
!
  write ( *, 1012)
  call gfslfs (hlp, k, nf, fmin, fmax, gain, freq, nprt, ldstak)
  write ( *, 1002) ierr
!
!  print storage from gfslfs
!
  if (ierr==0) then
    write ( *, 1004) (gain(i), i = 1, nf)
    write ( *, 1004) (freq(i), i = 1, nf)
  end if
!
!  test of gfarfs
!
  write ( *, 1014)
  call gfarfs (phi, iar, nf, fmin, fmax, gain, phas, freq, nprt, &
     ldstak)
  write ( *, 1002) ierr
!
!  print storage from gfarfs
!
  if (ierr==0) then
    write ( *, 1004) (gain(i), i = 1, nf)
    write ( *, 1004) (phas(i), i = 1, nf)
    write ( *, 1004) (freq(i), i = 1, nf)
  end if

  go to (100, 200, 300, 400), itest
!
!  test special cases
!
  100 itest = itest + 1
!
!  test of gfslfs
!
  fmin = 0.4e0
  fmax = 0.1e0
  nprt = 1
  write ( *, 1012)
  call gfslfs (hlp, k, nf, fmin, fmax, gain, freq, nprt, ldstak)
  write ( *, 1002) ierr
!
!  print storage from gfslfs
!
  if (ierr==0) then
    write ( *, 1004) (gain(i), i = 1, nf)
    write ( *, 1004) (freq(i), i = 1, nf)
  end if
!
!  test of gfarfs
!
  nprt = -1
  write ( *, 1014)
  call gfarfs (phi, iar, nf, fmin, fmax, gain, phas, freq, nprt, &
     ldstak)
  write ( *, 1002) ierr
!
!  print storage from gfarfs
!
  if (ierr==0) then
    write ( *, 1004) (gain(i), i = 1, nf)
    write ( *, 1004) (phas(i), i = 1, nf)
    write ( *, 1004) (freq(i), i = 1, nf)
  end if
!
!  test minimum problem size
!
  n = 3
  k = 1
  nprt = -1
  iar = 1
  nf = 1
  go to 20
!
!  test error conditions
!
  200 itest = itest + 1
  n = -5
  fc = 1.0e0
  nf = 0
  lphi = 0
  nfac = 1
  nd(1) = -1
  iod(1) = -1
  iar = 0
  k = -1
  go to 10
!
!  test ldstak
!
  300 itest = itest + 1
  n = 261
  nprt = 2
  fc = 1.0e0/22.0e0
  nf = 101
  fmin = 0.0e0
  fmax = 0.2e0
  lphi = 50
  nfac = 1
  nd(1) = 1
  iod(1) = 1
  iar = 1
  phi(1) = 0.6e0
  k = 41
  ymiss = 11.0e0
  ldstak = 0
  go to 30
!
  400 return

 1001 format ('test of lpcoef')
 1002 format (/' ierr is ', i5)
 1003 format ('test of hpcoef')
 1004 format (10e10.3)
 1005 format ('test of slflt')
 1006 format ('test of sample')
 1007 format ('test of lopass')
 1008 format ('test of hipass')
 1009 format ('test of arflt')
 1010 format ('test of difc')
 1011 format ('test of gfslf')
 1012 format ('test of gfslfs')
 1013 format ('test of gfarf')
 1014 format ('test of gfarfs')
 1015 format ('test of dif')
 1018 format ('test of difm')
 1019 format ('test of difmc')
 1020 format ('test of maflt')
end
subroutine xhist ( ldstak )

!*****************************************************************************80
!
!! XHIST tests the HIST family.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson and John Koontz,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        a loop index.
!     integer ierr
!        flag to indicate presence of error detected by preceding
!        starpac call.  (0 is ok, 1 is error)
!     integer ldsmin
!        the minimum amount of work area needed for a given problem.
!     integer ldstak
!        amount of work area.  size of dstak.
!     integer n
!        the length of the vector y.
!     integer ncell
!        the user supplied value for the number of cells in the
!        histogram.  if ncell is less than or equal to zero, the
!        number of cells to be used (ncells) will be calculated from
!        the number of observations.
!     integer nconst
!        length of the vector yconst.
!     integer nprtof
!        flag for no output (except error messages).
!     integer nprton
!        flag for full printout.
!     real y(84)
!        data vector for tests.
!     real yconst(10)
!        vector of constant data.
!     real ylb
!        the lower bound for selecting data from y for the histogram.
!     real ylong(200)
!        long vector of data
!     real ypath(10)
!        a vector of y values designed to force different paths
!        through the summation routines.
!     real yub
!        the upper bound for selecting data from y for the histogram.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     ylb,yub
  integer ldsmin,n,ncell,nconst
!
!  local arrays
  real &
     y(84),yconst(10),ylong(200),ypath(10)
!
!  external subroutines
  external hist,histc,iprint,ldscmp,nrand
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

!
!     data initializations.
!
  data n /84/
  data nconst /10/
  data ncell/10/
  data ylb/0.60e0/, yub/0.63e0/
!
!     davis-harrison r.h. data, pikes peak.
!
!     this is an arbitrarily chosen data set.
!
  data y( 1), y( 2), y( 3), y( 4) &
      / 0.6067e0, 0.6087e0, 0.6086e0, 0.6134e0/
  data y( 5), y( 6), y( 7) &
      / 0.6108e0, 0.6138e0, 0.6125e0/
  data y( 8), y( 9), y(10), y(11) &
      / 0.6122e0, 0.6110e0, 0.6104e0, 0.7213e0/
  data y(12), y(13), y(14) &
      / 0.7078e0, 0.7021e0, 0.7004e0/
  data y(15), y(16), y(17), y(18) &
      / 0.6981e0, 0.7242e0, 0.7268e0, 0.7418e0/
  data y(19), y(20), y(21) &
      / 0.7407e0, 0.7199e0, 0.6225e0/
  data y(22), y(23), y(24), y(25) &
      / 0.6254e0, 0.6252e0, 0.6267e0, 0.6218e0/
  data y(26), y(27), y(28) &
      / 0.6178e0, 0.6216e0, 0.6192e0/
  data y(29), y(30), y(31), y(32) &
      / 0.6191e0, 0.6250e0, 0.6188e0, 0.6233e0/
  data y(33), y(34), y(35) &
      / 0.6225e0, 0.6204e0, 0.6207e0/
  data y(36), y(37), y(38), y(39) &
      / 0.6168e0, 0.6141e0, 0.6291e0, 0.6231e0/
  data y(40), y(41), y(42) &
      / 0.6222e0, 0.6252e0, 0.6308e0/
  data y(43), y(44), y(45), y(46) &
      / 0.6376e0, 0.6330e0, 0.6303e0, 0.6301e0/
  data y(47), y(48), y(49) &
      / 0.6390e0, 0.6423e0, 0.6300e0/
  data y(50), y(51), y(52), y(53) &
      / 0.6260e0, 0.6292e0, 0.6298e0, 0.6290e0/
  data y(54), y(55), y(56) &
      / 0.6262e0, 0.5952e0, 0.5951e0/
  data y(57), y(58), y(59), y(60) &
      / 0.6314e0, 0.6440e0, 0.6439e0, 0.6326e0/
  data y(61), y(62), y(63) &
      / 0.6392e0, 0.6417e0, 0.6412e0/
  data y(64), y(65), y(66), y(67) &
      / 0.6530e0, 0.6411e0, 0.6355e0, 0.6344e0/
  data y(68), y(69), y(70) &
      / 0.6623e0, 0.6276e0, 0.6307e0/
  data y(71), y(72), y(73), y(74) &
      / 0.6354e0, 0.6197e0, 0.6153e0, 0.6340e0/
  data y(75), y(76), y(77) &
      / 0.6338e0, 0.6284e0, 0.6162e0/
  data y(78), y(79), y(80), y(81) &
      / 0.6252e0, 0.6349e0, 0.6344e0, 0.6361e0/
  data y(82), y(83), y(84) &
      / 0.6373e0, 0.6337e0, 0.6383e0/
!
!  check for sufficient work area length.
!
  if (ldstak.lt.300) then
    write ( *, 1000)
     return
  end if

  yconst(1:nconst) = 1.0e0
!
!  heading.
!
  write ( *,1150)
!
!  test 1.  check all error messages.
!
  write ( *,1160)
!
!  error 1, zero or fewer elements.
!
  write ( *,1180)
  call hist(y, 0, ldstak)
  write ( *, 1350)
  write ( *, 1360) y(1:n)
  write ( *,1170) ierr
  call histc(y, 0, ncell, ylb, yub, ldstak)
  write ( *, 1350)
  write ( *, 1360) y(1:n)
  write ( *,1170) ierr
!
!  error 2, not enough space in cstak.
!
  write ( *,1190)
  call ldscmp(2, 0, n, 0, 0, 0, 's', &
              min(nint(5.5+1.5*anint(log10(real(n)))),25),ldsmin)
  call hist(y, n, ldsmin-1)
  write ( *,1170) ierr
  write ( *,1195)
  call hist(y, n, ldsmin)
  write ( *,1170) ierr
  write ( *,1190)
  call ldscmp(2, 0, n, 0, 0, 0, 's', ncell, ldsmin)
  call histc(y, n, ncell, ylb, yub, ldsmin-1)
  write ( *,1170) ierr
  write ( *,1195)
  call histc(y, n, ncell, ylb, yub, ldsmin)
  write ( *,1170) ierr
!
!  constant y. (not an error)
!
  write ( *,1200)
  call hist(yconst, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1200)
  call histc(yconst, nconst, ncell, ylb, yub, ldstak)
  write ( *,1170) ierr
!
!  error 4, no data within user supplied limits
!
  write ( *, 1110)
  call histc(y, n, 0, 4.0e0, 10.0e0, ldstak)
  write ( *, 1170) ierr
!
!  test 2.  make a working run of each routine to check the output.
!
  write ( *,1300)
  write ( *,1310)
  call hist(y, n, ldstak)
  write ( *, 1350)
  write ( *, 1360) y(1:n)
  write ( *,1170) ierr

  write ( *,1340)
  call histc(y, n, ncell, ylb, yub, ldstak)
  write ( *, 1350)
  write ( *, 1360) y(1:n)
  write ( *,1170) ierr
!
!  run data set 6.7.
!
  ypath(1:10) = 0.0e0
  ypath(1) = -1.0e0
  ypath(10) = 1.0e0
  write ( *,1130)
  call hist(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *, 1130)
  call histc(ypath, nconst, 0, 0.0e0, 0.0e0, ldstak)
  write ( *, 1130)
  call histc(ypath, nconst, 1, 0.0e0, 0.0e0, ldstak)
  write ( *, 1130)
  call histc(ypath, nconst, 0, -0.5e0, 0.5e0, ldstak)
  write ( *, 1130)
  call histc(ypath, nconst, 0, 1.0e0, 4.0e0, ldstak)
!
!  run data set 6.8
!
  write ( *, 1120)
  call nrand (ylong, 200, 3254767)
  call hist (ylong, 200, ldstak)
  return

 1000 format ('1the dimension of dstak and the value of ldstak needed'/ &
    ' for histx must equal or exceed 300.  change driver'/ &
    ' and recall histx.')
 1110 format ('try no data within user supplied limits.')
 1120 format ('run hist on 200 pseudo-randon numbers')
 1130 format('run hist on -1, 8*0, 1.')
 1150 format ('test runs for the histogram family of routines.')
 1160 format('test 1.  generate one of each of the possible ', &
     'error messages.')
 1170 format('the value of ierr is ', i4)
 1180 format('try zero or fewer elements.')
 1190 format('1test with insufficient work area')
 1195 format(' test with exactly the right amount of work area.')
 1200 format('1try constant y. (not an error)')
 1300 format('test 4.  make working runs of all routines to check', &
     ' the output.')
 1310 format('run hist on the davis-harrison pikes peak data.')
 1340 format('run histc on the davis-harrison pikes peak data.')
 1350 format(/' print the data to insure the original order has', &
     ' been restored.')
 1360 format (7f10.5)
end
subroutine xlls ( lds )

!*****************************************************************************80
!
!! XLLS tests the linear least squares subroutines.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index.
!     integer ierr
!        the integer value designating whether any errors were
!        detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ivcv
!        the first dimension of the matrix vcv.
!     integer ixm
!        the first dimension of the matrix x.
!     integer j
!        an index.
!     integer lds
!       ..
!     integer ldsmin
!        the minimum length allowed for the array dstak.
!     integer ldstak
!        the length of the vector dstak in common cstak.
!     integer lpar
!        the actual length of the parameter array.
!     integer n
!        the number of observations.
!     integer ndeg
!        the degree of the polynomial model to be fit.
!     integer npar
!        the number of parameters.
!     integer nprt
!        the indicator variable used to designate the amount of
!        printed output.
!     real par(10)
!        the parameters  to be estimated.
!     real pv(50)
!        the predicted values.
!     real rand(1)
!        *
!     real res(50)
!        the residuals.
!     real rsd
!        the residual standard deviation.
!     real sdpv(50)
!        the standard deviations of the predicted values.
!     real sdres(50)
!        the standardized residuals.
!     real vcv(10,10)
!        the variance covariance matrix.
!     real wt(50)
!        the weights (a dummy vector in the unweighted case).
!     real x(50,9)
!        the independent variable.
!     real xm(50,10)
!        the independent variable.
!     real xm1(50,10)
!        the independent variable.
!     real y(50)
!        the dependent variable.
!     real y1(50)
!        the dependent variable.
!
  implicit none

  integer &
     lds
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     rsd,sum,term
  integer &
     i,ivcv,ixm,j,ldsmin,ldstak,lpar,n,ndeg,npar,nprt
!
!  local arrays
  real &
     par(10),pv(50),rand(1),res(50),sdpv(50),sdres(50),vcv(10,10), &
     wt(50),x(50,9),xm(50,10),xm1(50,10),y(50),y1(50)
!
!  external subroutines
  external fitsxp,genr,iprint,ldscmp,lls,llsp,llsps,llspw,llspws, &
     llss,llsw,llsws,nrand,setrv
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!  equivalences
  equivalence (xm(1,2),x(1,1))
!
  data      xm(1,1),  xm(1,2),  xm(1,3),  xm(1,4) &
      /      1.0e0, 42.2e0, 11.2e0, 31.9e0/
  data      xm(2,1),  xm(2,2),  xm(2,3),  xm(2,4) &
      /      1.0e0, 48.6e0, 10.6e0, 13.2e0/
  data      xm(3,1),  xm(3,2),  xm(3,3),  xm(3,4) &
      /      1.0e0, 42.6e0, 10.6e0, 28.7e0/
  data      xm(4,1),  xm(4,2),  xm(4,3),  xm(4,4) &
      /      1.0e0, 39.0e0, 10.4e0, 26.1e0/
  data      xm(5,1),  xm(5,2),  xm(5,3),  xm(5,4) &
      /      1.0e0, 34.7e0,  9.3e0, 30.1e0/
  data      xm(6,1),  xm(6,2),  xm(6,3),  xm(6,4) &
      /      1.0e0, 44.5e0, 10.8e0,  8.5e0/
  data      xm(7,1),  xm(7,2),  xm(7,3),  xm(7,4) &
      /      1.0e0, 39.1e0, 10.7e0, 24.3e0/
  data      xm(8,1),  xm(8,2),  xm(8,3),  xm(8,4) &
      /      1.0e0, 40.1e0, 10.0e0, 18.6e0/
  data      xm(9,1),  xm(9,2),  xm(9,3),  xm(9,4) &
      /      1.0e0, 45.9e0, 12.0e0, 20.4e0/
  data         y(1),     y(2),     y(3) &
      /    167.1e0,174.4e0,160.8e0/
  data         y(4),     y(5),     y(6) &
      /    162.0e0,140.8e0,174.6e0/
  data         y(7),     y(8),     y(9) &
      /    163.7e0,174.5e0,185.7e0/
!
!  set parameters necessary for the computations
!
  n = 9
  npar = 4
  ndeg = 3
  nprt = 2
  lpar = 10
  ivcv = 10
  ixm = 50
  ldstak = lds

  call setrv(wt, n, 1.0e0)
!
!  check error handling
!
!  error 1  -  non positive number of observations and parameter
!  number of parameters greater than n
!  ixm less than number of observations
!  ivcv less than number of parameters
!  lpar too small
!
  n = -5
  npar = 0
  ndeg = -1
  ixm = -10
  lpar = -1
  ivcv = -10
  nprt = -1
  write ( *,1200)
  write ( *,1000)
  call lls(y, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1020)
  call llsw(y, wt, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1040)
  call llsp(y, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1050)
  call llsps(y, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1060)
  call llspw(y, wt, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  n = 9
  npar = 4
  ndeg = 3
  ixm = 50
  lpar = -10
  ivcv = 10
!
!  error 2  -  lds too small, lpar too small
!
  ldstak = 0
  write ( *,1220)
  write ( *,1000)
  call lls(y, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1020)
  call llsw(y, wt, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1040)
  call llsp(y, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1050)
  call llsps(y, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1060)
  call llspw(y, wt, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  ldstak = lds
  nprt = 2
  lpar = 10
!
!  error 3  -  negative weights
!
  wt(1) = -1.0e0
  write ( *,1240)
  write ( *,1020)
  call llsw(y, wt, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1060)
  call llspw(y, wt, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  wt(1) = 1.0e0
!
!  error 4  -  too few positive weights
!
  call setrv(wt(2), n-1, 0.0e0)
  write ( *,1250)
  write ( *,1020)
  call llsw(y, wt, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1060)
  call llspw(y, wt, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call setrv(wt(2), n-1, 1.0e0)
!
!  check results from valid call
!
  write ( *,1260)
  write ( *,1000)
  call lls(y, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1260)
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1260)
  write ( *,1020)
  call llsw(y, wt, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1260)
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1260)
  write ( *,1040)
  call llsp(y, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1260)
  write ( *,1050)
  call llsps(y, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1260)
  write ( *,1060)
  call llspw(y, wt, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1260)
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
!
!  check results from exact fit
!
  n = npar
  ndeg = npar-1

  write ( *,1270)
  write ( *,1000)
  call lls(y, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1270)
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1270)
  write ( *,1040)
  call llsp(y, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1270)
  write ( *,1050)
  call llsps(y, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
!
  n = 9
!
  call setrv(wt(npar+1), n-npar, 0.0e0)
!
  write ( *,1270)
  write ( *,1020)
  call llsw(y, wt, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1270)
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1270)
  write ( *,1060)
  call llspw(y, wt, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1270)
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
!
  call setrv(wt(npar+1), n-npar, 1.0e0)
!
!     check results from rank deficient fit
!
  xm(1:n,5) = xm(1:n,4)
  write ( *,1280)
  write ( *,1000)
  call lls(y, xm, n, ixm, npar+1, res, ldstak)
  write ( *,1500) ierr
!
!     check results from a poorly scaled problem.
!
  do i = 1, n
     y1(i) = y(i) * 1.0e-8
     xm1(i,1:4) = xm(i,1:4)
     xm1(i,3) = xm1(i,3) * 1.0e+8
  end do

  write ( *,1290)
  write ( *,1000)
  call lls(y1, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1290)
  write ( *,1000)
  call lls(y, xm1, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1290)
  write ( *,1000)
  call lls(y1, xm1, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
!
!  minimum amount of work area.
!
  call ldscmp(15, 0, 0, 0, 0, 0, 's', &
              6*n + npar*(n+2*npar+5) + 1, ldsmin)

  write ( *,1300)
  write ( *,1000)
  call lls(y, xm, n, ixm, npar, res, ldsmin)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
!
!  check results for weighted analysis
!
  nprt = 1111
  call setrv(wt, n, 100.0e0)
  write ( *,1310)
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  wt(1) = 0.0e0
  wt(5) = 0.0e0
  wt(9) = 0.0e0

  write ( *,1310)
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  call setrv(wt, n, 100.0e0)

  call genr(wt, n, 1.0e0, 1.0e0)
  write ( *,1310)
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  call setrv(wt, n, 100.0e0)
!
!  check print control
!
  nprt = 1000
  write ( *,1320) nprt
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  nprt = 2000
  write ( *,1320) nprt
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  nprt = 200
  write ( *,1320) nprt
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  nprt = 20
  write ( *,1320) nprt
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  nprt = 2
  write ( *,1320) nprt
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  nprt = 0
  write ( *,1320) nprt
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
!
!  check results for n = 2, npar = id+1 = 1
!
  nprt = 2222
  n = 2
  npar = 1
  ndeg = 0
  write ( *,1330)
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1330)
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
!
!  check results for n = 1, npar = id+1 = 1
!
  nprt = 2222
  n = 1
  npar = 1
  ndeg = 0
  write ( *,1330)
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1330)
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  n = 9
  npar = 4
  ndeg = 3
!
!  ill-conditioned
!
  do i = 1, 50
     term = 1.0e0
     sum = 0.0e0
     do j = 1, 6
        xm1(i,j) = term
        sum = sum + term
        term = (i-1)*term
     end do
     y1(i) = sum
  end do

  n = 21
  npar = 6
  ndeg = 5
  write ( *,1340)
  write ( *,1010)
  call llss(y1, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1340)
  write ( *,1050)
  call llsps(y1, xm1(1,2), n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  n = 50
  npar = 6
  ndeg = 5
  call nrand(rand, 1, 223)
  do i = 1, n
     call nrand(rand, 1, 0)
     y1(i) = y1(i) + rand(1)
  end do
  write ( *,1340)
  write ( *,1010)
  call llss(y1, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1340)
  write ( *,1050)
  call llsps(y1, xm1(1,2), n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
!
  nprt = 1000
  write ( *,1340)
  write ( *,1010)
  call llss(y1, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1340)
  write ( *,1050)
  call llsps(y1, xm1(1,2), n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  n = 45
  call setrv(wt, n, 1.0e0)
  write ( *,1340)
  write ( *,1010)
  call llsws(y1, wt, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  n = 44
  write ( *,1340)
  write ( *,1010)
  call llsws(y1, wt, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  n = 41
  write ( *,1340)
  write ( *,1010)
  call llsws(y1, wt, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  n = 40
  write ( *,1340)
  write ( *,1010)
  call llsws(y1, wt, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  return

 1000 format (' call to lls   ')
 1010 format (' call to llss  ')
 1020 format (' call to llsw  ')
 1030 format (' call to llsws ')
 1040 format (' call to llsp  ')
 1050 format (' call to llsps ')
 1060 format (' call to llspw ')
 1070 format (' call to llspws')
 1200 format ('1miscellaneous errors  -  test 1')
 1220 format ('1miscellaneous errors  -  test 2')
 1240 format ('1negative weights')
 1250 format ('1too few positive weights')
 1260 format ('1valid problem')
 1270 format ('1zero residual problem')
 1280 format ('1rank deficient problem')
 1290 format ('1poorly scaled problem')
 1300 format ('1minimum work area size')
 1310 format ('1weighted analysis')
 1320 format ('1check print control  -  nprt = ', i5)
 1330 format ('1check minimum problem size')
 1340 format ('1ill-conditioned problem')
 1430 format (//4h res/ (1x, e22.14))
 1500 format (/' ierr = ', i5)
end
subroutine xnlsd ( ldstak )

!*****************************************************************************80
!
!! XNLSD demonstrates the nonlinear least squares family.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     real delta
!        the maximum change allowed in the model parameters at the
!        first iteration.
!     external drv1a
!        the name of the user supplied subroutine which computes the
!        derivative (jacobian) matrix of the model.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer idrvck
!        the variable used to indicate whether the derivatives are
!        to be checked (idrvck = 1) or not (idrvck = 0).
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr .ge. 1, errors were detected.
!     integer ifixed(10)
!        the indicator values used to designate whether the
!        parameters are to be optimized or are to be held fixed.  if
!        ifixed(i).ne.0, then par(i) will be optimized.  if
!        ifixed(i)==0, then par(i) will be held fixed.
!        ifixed(i).lt.0, then all par(i),i=1,npar, will be optimized..
!     integer ivaprx
!        an indicator value used to designate which option is to be used
!        to compute the variance covariance matrix (vcv), where
!        ivaprx le 0 indicates the the default option will be used
!        ivaprx eq 1 indicates the vcv is to be computed by
!                       inverse(hessian)*transpose(jacobian)*jacobian
!                          *inverse(hessian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 2 indicates the vcv is to be computed by
!                       inverse(hessian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 3 indicates the vcv is to be computed by
!                       inverse(transpose(jacobian)*jacobian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 4 indicates the vcv is to be computed by
!                       inverse(hessian)*transpose(jacobian)*jacobian
!                          *inverse(hessian)
!                    using only the model subroutine
!        ivaprx eq 5 indicates the vcv is to be computed by
!                       inverse(hessian)
!                    using only the model subroutine
!        ivaprx eq 6 indicates the vcv is to be computed by
!                       inverse(transpose(jacobian)*jacobian)
!                    using only the model subroutine
!        ivaprx ge 7 indicates the default option will be used
!     integer ivcv
!        the first dimension of the variance covariance matrix vcv.
!     integer ixm1
!        the first dimension of the independent variable array.
!     integer ldsa1, ldsn1a, ldsn1b
!        the minimum length allowed for the array dstak
!        for the routines with analytic derivatives and
!        numerical derivatives, respectively.
!     integer ldsmin
!        the minimum length of the array dstak allowed.
!     integer ldstak
!        the length of the array dstak.
!     integer m1
!        the number of independent variables.
!     integer mit
!        the maximum number of iterations allowed.
!     external mdl1
!        the name of the user supplied subroutine which computes the
!        predicted values based on the current parameter estimate.
!     integer n1
!        the number of observations.
!     integer npare
!        the number of parameters estimated by the routine.
!     integer npar1
!        the number of unknown parameters in the model.
!     integer nnzw
!        the number of non zero weights.
!     integer nprt
!        the parameter used to indicate how much printed output is
!        to be provided.
!     integer ntest
!        the number of the current test.
!     real par1(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real pv(100)
!        the starting location in rstak/dstak of
!        the predicted value based on the current parameter estimates
!     real res(100)
!        the residuals from the fit.
!     real rsd
!        the value of the residual standard deviation at the solution.
!     real scale(10)
!        a value to indicate use of the default values of
!        the typical size of the unknown parameters.
!     real sdpv(100)
!        the starting location in rstak/dstak of
!        the standard deviation of the predicted value.
!     real sdres(100)
!        the starting location in rstak/dstak of the
!        the standard deviations of the residuals.
!     real stopp
!        the stopping criterion for the test based on the maximum scaled
!        relative change in the elements of the model parameter vector
!     real stopss
!        the stopping criterion for the test based on the ratio of the
!        predicted decrease in the residual sum of squares (computed
!        by starpac) to the current residual sum of squares estimate.
!     real stp(10)
!        the rcstep size array.
!     real vcv(6,6)
!        the covariance matrix.
!     real wt(100)
!        the user supplied weights.
!     real xm1(10,2)
!        the array in which one row of the independent variable array
!        is stored.
!     real y1(10)
!        the array of the dependent variable.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     delta,rsd,stopp,stopss
  integer &
     idrvck,ivaprx,ivcv,ixm1,ldsa1,ldsmin,ldsn1a,ldsn1b, &
     m1,mit,n1,nnzw,npar1,npare,nprt,ntest
!
!  local arrays
  real &
     par1(10),pv(100),res(100),scale(10),sdpv(100),sdres(100), &
     stp(10),vcv(6,6),wt(100),xm1(10,2),y1(10)
  integer &
     ifixed(10)
!
!  external subroutines
  external drv1a,fitxsp,iprint,ldscmp,mdl1,nl2x,nls,nlsc,nlsd,nlsdc, &
     nlsds,nlss,nlsw,nlswc,nlswd,nlswdc,nlswds,nlsws,nlsx1, &
     nlsx2,setrv
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!     define constants
!
  data y1(1), y1(2), y1(3), y1(4), y1(5), y1(6) &
     /2.138e0, 3.421e0, 3.597e0, 4.340e0, 4.882e0, 5.660e0/

  data xm1(1,1), xm1(2,1), xm1(3,1), xm1(4,1), xm1(5,1), xm1(6,1) &
     /1.309e0, 1.471e0, 1.490e0, 1.565e0, 1.611e0, 1.680e0/

  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)

  call setrv(wt, n1, 1.0e0)

  call ldscmp(6, 0, 60+2*npar1, 0, 0, 0, &
     's', 94+n1*(3+npar1)+npar1*(3*npar1+35)/2, ldsa1)
  call ldscmp(14, 0, max(2*(n1+npar1),60+2*npar1), 0, 0, 0, &
     's', max(10*n1,94+n1*(3+npar1)+npar1*(3*npar1+37)/2), ldsn1a)
  call ldscmp(14, 0, 60+2*npar1, 0, 0, 0, &
     's', 94+n1*(3+npar1)+npar1*(3*npar1+37)/2, ldsn1b)

  ldsmin = max(ldsa1, ldsn1a, ldsn1b)

  if (ldsmin.le.ldstak) go to 5

  write ( *, 1140) ldsmin
  return

    5 continue

  ntest = 0
!
!  test on normal statement
!
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1000)
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nls(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldsn1a)
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1010)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsc(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldsn1b, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1020)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlss(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldsn1b, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1030)
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsw(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldsn1a)
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1040)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswc(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldsn1b, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1050)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsws(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldsn1b, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1060)
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsd(y1, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, res, &
     ldsa1)
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1070)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsdc(y1, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, res, &
     ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1080)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsds(y1, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, res, &
     ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1090)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswd(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldsa1)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1100)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)
!
!     test nl2sol and nl2sno directly
!
  write ( *,1320)
  call nl2x

  return

 1000 format ('test of nls'  )
 1010 format ('test of nlsc'  )
 1020 format ('test of nlss'  )
 1030 format ('test of nlsw' )
 1040 format ('test of nlswc' )
 1050 format ('test of nlsws' )
 1060 format ('test of nlsd' )
 1070 format ('test of nlsdc' )
 1080 format ('test of nlsds' )
 1090 format ('test of nlswd')
 1100 format ('test of nlswdc')
 1110 format ('test of nlswds')
 1120 format (/'returned results (-1 indicates ', &
     'value not changed by called subroutine)'//' ierr is ', i3)
 1130 format ('normal problem')
 1140 format ('ldstak must be greater than or equal to ', i6)
 1320 format ('test of nl2sol and nl2sno called directly')
 1330 format ('nonlinear least squares estimation subroutine test number', &
     i5/)
 1340 format (' input   -  ifixed(1) = ', i6, 9x, ', stp(1) = ', &
     g15.8, ',    mit = ',i5, ', stopss = ', g15.8, ', stopp = ', &
     g15.8/13x, 'scale(1) = ', g15.8, ',  delta = ', g15.8, &
     ', ivaprx = ', i5, ',   nprt = ', i5)
 1350 format (//24h output  -  ifixed(1) = , i6, 9x, 11h, stp(1) = , &
     g15.8, 11h,    mit = ,i5, 11h, stopss = , g15.8, 10h, stopp = , &
     g15.8/13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i5//)
 1360 format (24h input   -  ifixed(1) = , i6, 9x, 11h, idrvck = , &
     i5, 10x, 11h,    mit = ,i5, 11h, stopss = , g15.8, &
     10h, stopp = , g15.8/ &
     13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i6)
 1370 format (//24h output  -  ifixed(1) = , i6, 9x, 11h, idrvck = , &
     i5, 10x, 11h,    mit = ,i5, 11h, stopss = , g15.8, &
     10h, stopp = , g15.8/ &
     13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i5//)
end
subroutine xnlse ( ldstak )

!*****************************************************************************80
!
!! XNLSE demonstrate the nonlinear least squares family.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     real delta
!        the maximum change allowed in the model parameters at the
!        first iteration.
!     external drv1a, drv1b
!        the name of the user supplied subroutine which computes the
!        derivative (jacobian) matrix of the model.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index value.
!     integer idrvck
!        the variable used to indicate whether the derivatives are
!        to be checked (idrvck = 1) or not (idrvck = 0).
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr .ge. 1, errors were detected.
!     integer ifixed(10)
!        the indicator values used to designate whether the
!        parameters are to be optimized or are to be held fixed.  if
!        ifixed(i).ne.0, then par(i) will be optimized.  if
!        ifixed(i)==0, then par(i) will be held fixed.
!        ifixed(i).lt.0, then all par(i),i=1,npar, will be optimized..
!     integer ivaprx
!        an indicator value used to designate which option is to be used
!        to compute the variance covariance matrix (vcv), where
!        ivaprx le 0 indicates the the default option will be used
!        ivaprx eq 1 indicates the vcv is to be computed by
!                       inverse(hessian)*transpose(jacobian)*jacobian
!                          *inverse(hessian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 2 indicates the vcv is to be computed by
!                       inverse(hessian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 3 indicates the vcv is to be computed by
!                       inverse(transpose(jacobian)*jacobian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 4 indicates the vcv is to be computed by
!                       inverse(hessian)*transpose(jacobian)*jacobian
!                          *inverse(hessian)
!                    using only the model subroutine
!        ivaprx eq 5 indicates the vcv is to be computed by
!                       inverse(hessian)
!                    using only the model subroutine
!        ivaprx eq 6 indicates the vcv is to be computed by
!                       inverse(transpose(jacobian)*jacobian)
!                    using only the model subroutine
!        ivaprx ge 7 indicates the default option will be used
!     integer ivcv
!        the first dimension of the variance covariance matrix vcv.
!     integer ixm1, ixm3
!        the first dimension of the independent variable array.
!     integer ldsa1, ldsn3a
!        the minimum length allowed for the array dstak
!        for the routines with analytic derivatives and
!        numerical derivatives, respectively.
!     integer ldsmin
!        the minimum length of the array dstak allowed.
!     integer ldstak
!        the length of the array dstak.
!     integer m1, m3
!        the number of independent variables.
!     integer mit
!        the maximum number of iterations allowed.
!     external mdl1, mdl3
!        the name of the user supplied subroutine which computes the
!        predicted values based on the current parameter estimate.
!     integer n1, n3
!        the number of observations.
!     integer npare
!        the number of parameters estimated by the routine.
!     integer npar1, npar3
!        the number of unknown parameters in the model.
!     integer nnzw
!        the number of non zero weights.
!     integer nprt
!        the parameter used to indicate how much printed output is
!        to be provided.
!     integer ntest
!        the number of the current test.
!     real par1(10), par3(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real pv(100)
!        the starting location in rstak/dstak of
!        the predicted value based on the current parameter estimates
!     real res(100)
!        the residuals from the fit.
!     real rsd
!        the value of the residual standard deviation at the solution.
!     real scale(10)
!        a value to indicate use of the default values of
!        the typical size of the unknown parameters.
!     real sdpv(100)
!        the starting location in rstak/dstak of
!        the standard deviation of the predicted value.
!     real sdres(100)
!        the starting location in rstak/dstak of the
!        the standard deviations of the residuals.
!     real stopp
!        the stopping criterion for the test based on the maximum scaled
!        relative change in the elements of the model parameter vector
!     real stopss
!        the stopping criterion for the test based on the ratio of the
!        predicted decrease in the residual sum of squares (computed
!        by starpac) to the current residual sum of squares estimate.
!     real stp(10)
!        the rcstep size array.
!     real vcv(6,6)
!        the covariance matrix.
!     real wt(100)
!        the user supplied weights.
!     real xm1(10,2), xm3(101,5)
!        the array in which one row of the independent variable array
!        is stored.
!     real y1(10), y3(100)
!        the array of the dependent variable.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     delta,rsd,stopp,stopss
  integer &
     i,idrvck,ivaprx,ivcv,ixm1,ixm3,ldsa1,ldsmin,ldsn3a, &
     m1,m3,mit,n1,n3,nnzw,npar1,npar3,npare,nprt,ntest
!
!  local arrays
  real &
     par1(10),par3(10),pv(100),res(100),scale(10),sdpv(100), &
     sdres(100),stp(10),vcv(6,6),wt(100),xm1(10,2),xm3(101,5), &
     y1(10),y3(100)
  integer &
     ifixed(10)
!
!  external subroutines
  external drv1a,drv1b,iprint,ldscmp,mdl1,mdl3,nls,nlsc,nlsd,nlsdc, &
     nlsds,nlss,nlsw,nlswc,nlswd,nlswdc,nlswds,nlsws,nlsx1, &
     nlsx2,setrv
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!     define constants
!
  data y1(1), y1(2), y1(3), y1(4), y1(5), y1(6) &
     /2.138e0, 3.421e0, 3.597e0, 4.340e0, 4.882e0, 5.660e0/

  data xm1(1,1), xm1(2,1), xm1(3,1), xm1(4,1), xm1(5,1), xm1(6,1) &
     /1.309e0, 1.471e0, 1.490e0, 1.565e0, 1.611e0, 1.680e0/

  data n3 /50/, m3 /5/, ixm3 /101/, npar3 /5/

  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)

  call setrv(wt, n3, 1.0e0)

  call ldscmp(6, 0, 60+2*npar1, 0, 0, 0, &
     's', 94+n1*(3+npar1)+npar1*(3*npar1+35)/2, ldsa1)
  call ldscmp(14, 0, max(2*(n3+npar3),60+2*npar3), 0, 0, 0, &
     's', max(10*n3,94+n3*(3+npar3)+npar3*(3*npar3+37)/2), ldsn3a)

  ldsmin = max(ldsa1, ldsn3a)

  if (ldsmin.le.ldstak) go to 5

  write ( *, 1190) ldsmin
  return

    5 continue

  do i=1,n3
     xm3(i,1) = 1.0e0
     xm3(i,2) = i
     xm3(i,3) = xm3(i,2)*xm3(i,2)
     xm3(i,4) = xm3(i,3)*xm3(i,2)
     xm3(i,5) = xm3(i,4)*xm3(i,2)
     y3(i) = xm3(i,1) + xm3(i,2) + xm3(i,3) + xm3(i,4) + xm3(i,5)
  end do

  ntest = 0
!
!  check error handling
!
!  test 1  -  problem specification
!
  n1 = -5
  m1 = -1
  ixm1 = -10
  npar1 = 0
  ivcv = -10

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1140)

  write ( *,1000)
  call nls(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldstak)
  write ( *,1120) ierr

  write ( *,1010)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsc(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldstak, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  npar1 = 8
  n1 = 2
  write ( *,1020)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlss(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldstak, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  write ( *,1030)
  call nlsw(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldstak)
  write ( *,1120) ierr

  write ( *,1040)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlswc(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldstak, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  write ( *,1050)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsws(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldstak, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  write ( *,1060)
  call nlsd(y1, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, res, &
     ldstak)
  write ( *,1120) ierr

  write ( *,1070)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsdc(y1, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, res, &
     ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  n1 = 15
  write ( *,1080)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsds(y1, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, res, &
     ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  write ( *,1090)
  call nlswd(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak)
  write ( *,1120) ierr

  write ( *,1100)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
!
!  test 2  -  weights and control values
!
  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)

  wt(n1) = -1.0e0
  stp(1) = 1.0e0
  stp(2) = 0.0e0
  scale(1) = 1.0e0
  scale(2) = 0.0e0
  ifixed(1:npar1) = 1

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1150)

  write ( *,1050)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsws(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldstak, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  write ( *,1090)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlswd(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak)
  call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
!
!  test 3  -  too few positive weights
!
  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)
!
  call setrv(wt(2), n1-1, 0.0e0)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1160)

  write ( *,1030)
  call nlsw(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldstak)
  write ( *,1120) ierr

  write ( *,1110)
  call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1120) ierr
!
!
!  test 4  -  definite error in derivative
!
  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)

  call setrv(wt, n1, 1.0e0)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1170)

  write ( *,1060)
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsd(y1, xm1, n1, m1, ixm1, mdl1, drv1b, par1, npar1, res, &
     ldstak)
  write ( *,1120) ierr
!
!  test 5  -  possible error in derivative
!
  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)
  idrvck = 1
  nprt = 10000

  call setrv(wt, n1, 1.0e0)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1180)
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  par1(1) = 0.0e0
  write ( *,1070)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsdc(y1, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, res, &
     ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
!
!  test 6 -  insufficient work area length
!
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1230)
  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)

  write ( *,1000)
  call nlsx1(3, par3, npar3, pv, sdpv, res, sdres, vcv, n3, ivcv, &
     nnzw, npare, rsd)
  call nls(y3, xm3, n3, m3, ixm3, mdl3, par3, npar3, res, ldsn3a-1)
  write ( *,1120) ierr

  write ( *,1090)
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n3, ivcv, &
     nnzw, npare, rsd)
  call nlswd(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldsa1-1)
  write ( *,1120) ierr

  return

 1000 format ('test of nls'  )
 1010 format ('test of nlsc'  )
 1020 format ('test of nlss'  )
 1030 format ('test of nlsw' )
 1040 format ('test of nlswc' )
 1050 format ('test of nlsws' )
 1060 format ('test of nlsd' )
 1070 format ('test of nlsdc' )
 1080 format ('test of nlsds' )
 1090 format ('test of nlswd')
 1100 format ('test of nlswdc')
 1110 format ('test of nlswds')
 1120 format (/'returned results (-1 indicates ', &
     'value not changed by called subroutine)'//'ierr is ', i3)
 1140 format (46h error handling test 1 - problem specification)
 1150 format (51h error handling test 2 - weights and control values)
 1160 format (49h error handling test 3 - too few positive weights)
 1170 format (53h error handling test 4 - definite error in derivative)
 1180 format (53h error handling test 5 - possible error in derivative)
 1190 format (45h1 *** ldstak must be greater than or equal to , i6)
 1230 format (' error handling test 6 - insufficient work area length')
 1330 format (54h1nonlinear least squares estimation subroutine test nu, &
     4hmber, i5/)
 1340 format (24h input   -  ifixed(1) = , i6, 9x, 11h, stp(1) = , &
     g15.8, 11h,    mit = ,i5, 11h, stopss = , g15.8, 10h, stopp = , &
     g15.8/13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i5)
 1350 format (//24h output  -  ifixed(1) = , i6, 9x, 11h, stp(1) = , &
     g15.8, 11h,    mit = ,i5, 11h, stopss = , g15.8, 10h, stopp = , &
     g15.8/13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i5//)
 1360 format (24h input   -  ifixed(1) = , i6, 9x, 11h, idrvck = , &
     i5, 10x, 11h,    mit = ,i5, 11h, stopss = , g15.8, &
     10h, stopp = , g15.8/ &
     13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i6)
 1370 format (//24h output  -  ifixed(1) = , i6, 9x, 11h, idrvck = , &
     i5, 10x, 11h,    mit = ,i5, 11h, stopss = , g15.8, &
     10h, stopp = , g15.8/ &
     13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i5//)
end
subroutine xnlst ( ldstak )

!*****************************************************************************80
!
!! XNLST demonstrates the nonlinear least squares family.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     real delta
!        the maximum change allowed in the model parameters at the
!        first iteration.
!     external drv1a, drv2, drv3
!        the name of the user supplied subroutine which computes the
!        derivative (jacobian) matrix of the model.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index value.
!     integer idrvck
!        the variable used to indicate whether the derivatives are
!        to be checked (idrvck = 1) or not (idrvck = 0).
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr .ge. 1, errors were detected.
!     integer ifixed(10)
!        the indicator values used to designate whether the
!        parameters are to be optimized or are to be held fixed.  if
!        ifixed(i).ne.0, then par(i) will be optimized.  if
!        ifixed(i)==0, then par(i) will be held fixed.
!        ifixed(i).lt.0, then all par(i),i=1,npar, will be optimized..
!     integer ivctst(9)
!        variance-covariance code test values.
!     integer ivaprx
!        an indicator value used to designate which option is to be used
!        to compute the variance covariance matrix (vcv), where
!        ivaprx le 0 indicates the the default option will be used
!        ivaprx eq 1 indicates the vcv is to be computed by
!                       inverse(hessian)*transpose(jacobian)*jacobian
!                          *inverse(hessian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 2 indicates the vcv is to be computed by
!                       inverse(hessian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 3 indicates the vcv is to be computed by
!                       inverse(transpose(jacobian)*jacobian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 4 indicates the vcv is to be computed by
!                       inverse(hessian)*transpose(jacobian)*jacobian
!                          *inverse(hessian)
!                    using only the model subroutine
!        ivaprx eq 5 indicates the vcv is to be computed by
!                       inverse(hessian)
!                    using only the model subroutine
!        ivaprx eq 6 indicates the vcv is to be computed by
!                       inverse(transpose(jacobian)*jacobian)
!                    using only the model subroutine
!        ivaprx ge 7 indicates the default option will be used
!     integer ivcv
!        the first dimension of the variance covariance matrix vcv.
!     integer ixm1, ixm2, ixm3
!        the first dimension of the independent variable array.
!     integer ldsa1, ldsn1b
!        the minimum length allowed for the array dstak
!        for the routines with analytic derivatives and
!        numerical derivatives, respectively.
!     integer ldstak
!        the length of the array dstak.
!     integer m1, m2, m3
!        the number of independent variables.
!     integer mit
!        the maximum number of iterations allowed.
!     external mdl1, mdl2, mdl3
!        the name of the user supplied subroutine which computes the
!        predicted values based on the current parameter estimate.
!     integer n1, n2, n3
!        the number of observations.
!     integer npare
!        the number of parameters estimated by the routine.
!     integer npar1, npar2, npar3
!        the number of unknown parameters in the model.
!     integer nnzw
!        the number of non zero weights.
!     integer nprt
!        the parameter used to indicate how much printed output is
!        to be provided.
!     integer ntest
!        the number of the current test.
!     real par1(10), par2(10), par3(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real pv(100)
!        the starting location in rstak/dstak of
!        the predicted value based on the current parameter estimates
!     real res(100)
!        the residuals from the fit.
!     real rsd
!        the value of the residual standard deviation at the solution.
!     real scale(10)
!        a value to indicate use of the default values of
!        the typical size of the unknown parameters.
!     real sdpv(100)
!        the starting location in rstak/dstak of
!        the standard deviation of the predicted value.
!     real sdres(100)
!        the starting location in rstak/dstak of the
!        the standard deviations of the residuals.
!     real stop(8)
!        stopping criteria test variable.
!     real stopp
!        the stopping criterion for the test based on the maximum scaled
!        relative change in the elements of the model parameter vector
!     real stopss
!        the stopping criterion for the test based on the ratio of the
!        predicted decrease in the residual sum of squares (computed
!        by starpac) to the current residual sum of squares estimate.
!     real stp(10)
!        the rcstep size array.
!     real vcv(6,6)
!        the covariance matrix.
!     real wt(100)
!        the user supplied weights.
!     real xm1(10,2), xm2(10,3), xm3(101,5)
!        the array in which one row of the independent variable array
!        is stored.
!     real y1(10), y2(10), y3(100)
!        the array of the dependent variable.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     delta,rsd,stopp,stopss
  integer &
     i,idrvck,ivaprx,ivcv,ixm1,ixm2,ixm3,ldsa1,ldsmin, &
     ldsn1b,m1,m2,m3,mit,n1,n2,n3,nnzw,npar1,npar2,npar3,npare, &
     nprt,ntest
!
!  local arrays
  real &
     par1(10),par2(10),par3(10),pv(100),res(100),scale(10), &
     sdpv(100),sdres(100),stop(8),stp(10),vcv(6,6),wt(100), &
     xm1(10,2),xm2(10,3),xm3(101,5),y1(10),y2(10),y3(100)
  integer &
     ifixed(10),ivctst(9)
!
!  external functions
  real &
     rmdcon
  external rmdcon
!
!  external subroutines
  external drv1a,drv2,drv3,fitxsp,iprint,ldscmp,mdl1,mdl2,mdl3,nlsc, &
     nlsdc,nlsds,nlss,nlswc,nlswdc,nlswds,nlsws,nlsx1,nlsx2, &
     setrv

!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!     define constants
!
  data y1(1), y1(2), y1(3), y1(4), y1(5), y1(6) &
     /2.138e0, 3.421e0, 3.597e0, 4.340e0, 4.882e0, 5.660e0/

  data xm1(1,1), xm1(2,1), xm1(3,1), xm1(4,1), xm1(5,1), xm1(6,1) &
     /1.309e0, 1.471e0, 1.490e0, 1.565e0, 1.611e0, 1.680e0/

  data n2 /10/, m2 /3/, ixm2 /10/, npar2 /3/

  data n3 /50/, m3 /5/, ixm3 /101/, npar3 /5/

  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)

  call setrv(wt, n3, 1.0e0)

  call ldscmp(6, 0, 60+2*npar1, 0, 0, 0, &
     's', 94+n1*(3+npar1)+npar1*(3*npar1+35)/2, ldsa1)
  call ldscmp(14, 0, 60+2*npar1, 0, 0, 0, &
     's', 94+n1*(3+npar1)+npar1*(3*npar1+37)/2, ldsn1b)
  call ldscmp(14, 0, max(2*(n3+npar3),60+2*npar3), 0, 0, 0, &
     's', max(10*n3,94+n3*(3+npar3)+npar3*(3*npar3+37)/2), ldsmin)

  if (ldsmin.le.ldstak) go to 5

  write ( *, 1000) ldsmin
  return

    5 continue

  do i=1,n2
     y2(i) = 0.0e0
     xm2(i,1) = i
     xm2(i,2) = i + 0.125e0
     xm2(i,3) = i + 0.25e0
  end do

  do i=1,n3
     xm3(i,1) = 1.0e0
     xm3(i,2) = i
     xm3(i,3) = xm3(i,2)*xm3(i,2)
     xm3(i,4) = xm3(i,3)*xm3(i,2)
     xm3(i,5) = xm3(i,4)*xm3(i,2)
     y3(i) = xm3(i,1) + xm3(i,2) + xm3(i,3) + xm3(i,4) + xm3(i,5)
  end do

  ntest = 0
!
!  test checking of control criteria
!
  write ( *,1240)

  stop(1) = rmdcon(3)
  stop(2) = 0.1e0
  stop(3) = 0.9e0*rmdcon(3)
  stop(4) = 0.11e0
  stop(5) = 0.0e0
  stop(6) = 1.0e0
  stop(7) = -1.0e0
  stop(8) = 1.1e0

  nprt = 11000
  mit = 0
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1250) mit
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  write ( *,1100)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, scale(1), &
     delta, ivaprx, nprt
  call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, scale(1), &
     delta, ivaprx, nprt
  write ( *,1120) ierr

  mit = 1
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1250) mit
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  write ( *,1100)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, scale(1), &
     delta, ivaprx, nprt
  call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, scale(1), &
     delta, ivaprx, nprt
  write ( *,1120) ierr

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1250) mit
  call nlsx1(2, par2, npar2, pv, sdpv, res, sdres, vcv, n2, ivcv, &
     nnzw, npare, rsd)
  write ( *,1050)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsws(y2, wt, xm2, n2, m2, ixm2, mdl2, par2, npar2, res, &
     ldstak, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par2, pv, sdpv, res, sdres, vcv, n2, npar2, ivcv, &
     nnzw, npare, rsd)

  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)

  do i=1,4
     ntest = ntest + 1
     write ( *,1330) ntest
     write ( *,1260) stop(i)
     call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, &
        ivcv, nnzw, npare, rsd)
     write ( *,1100)
     write ( *,1360) ifixed(1), idrvck, mit, stop(i), stopp, &
        scale(1), delta, ivaprx, nprt
     call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, &
        npar1, res, ldstak, ifixed, idrvck, mit, stop(i), stopp, &
        scale, delta, ivaprx, nprt)
     write ( *,1370) ifixed(1), idrvck, mit, stop(i), stopp, &
        scale(1), delta, ivaprx, nprt
     write ( *,1120) ierr
  end do

  do i=5,8
     ntest = ntest + 1
     write ( *,1330) ntest
     write ( *,1270) stop(i)
     call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, &
        ivcv, nnzw, npare, rsd)
     write ( *,1100)
     write ( *,1360) ifixed(1), idrvck, mit, stopss, stop(i), &
        scale(1), delta, ivaprx, nprt
     call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, &
        npar1, res, ldstak, ifixed, idrvck, mit, stopss, stop(i), &
        scale, delta, ivaprx, nprt)
     write ( *,1370) ifixed(1), idrvck, mit, stopss, stop(i), &
        scale(1), delta, ivaprx, nprt
     write ( *,1120) ierr
  end do

  nprt = 100000

  do i=1,6
     ntest = ntest + 1
     write ( *,1330) ntest
     write ( *,1280) nprt
     call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, &
        ivcv, nnzw, npare, rsd)
     write ( *,1110)
     write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
        scale(1), delta, ivaprx, nprt
     call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, &
        npar1, res, ldsa1, ifixed, idrvck, mit, stopss, stopp, &
        scale, delta, ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, &
        sdres, vcv, ivcv)
     write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
        scale(1), delta, ivaprx, nprt
     write ( *,1120) ierr
     call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
        nnzw, npare, rsd)
     nprt = nprt/10
  end do

  nprt = 11000

  do i=1,2
     ntest = ntest + 1
     write ( *,1330) ntest
     write ( *,1280) nprt
     call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, &
        ivcv, nnzw, npare, rsd)
     write ( *,1110)
     write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
        scale(1), delta, ivaprx, nprt
     call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, &
        npar1, res, ldsa1, ifixed, idrvck, mit, stopss, stopp, &
        scale, delta, ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, &
        sdres, vcv, ivcv)
     write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
        scale(1), delta, ivaprx, nprt
     write ( *,1120) ierr
     call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
        nnzw, npare, rsd)
     nprt = 11001
  end do

  nprt = 0
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt

  write ( *,1010)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsc(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldsn1b, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1020)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlss(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldsn1b, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1040)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswc(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldsn1b, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1050)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsws(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldsn1b, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1070)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsdc(y1, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, res, &
     ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1080)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsds(y1, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, res, &
     ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1100)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  nprt = -1

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)
!
!  test parameter handling**
!
  write ( *,1190)
!
!  all zero
!
  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)
  stp(1) = -1.0e0
  scale(1) = -1.0e0
  delta = -1.0e0
  nprt = 11000

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1200)
  call nlsx1(4, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  write ( *,1010)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsc(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldstak, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  call nlsx1(4, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1200)
  write ( *,1100)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)
!
!  test with constant y
!
!  constant y=0
!
  nprt = 21222
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1210)
  call nlsx1(2, par2, npar2, pv, sdpv, res, sdres, vcv, n2, ivcv, &
     nnzw, npare, rsd)
  write ( *,1050)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsws(y2, wt, xm2, n2, m2, ixm2, mdl2, par2, npar2, res, &
     ldstak, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par2, pv, sdpv, res, sdres, vcv, n2, npar2, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1210)
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(2, par2, npar2, pv, sdpv, res, sdres, vcv, n2, ivcv, &
     nnzw, npare, rsd)
  call nlswds(y2, wt, xm2, n2, m2, ixm2, mdl2, drv2, par2, npar2, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par2, pv, sdpv, res, sdres, vcv, n2, npar2, ivcv, &
     nnzw, npare, rsd)
!
!  test with linear model
!
  nprt = 11212
  ivaprx = 1

  ifixed(1:npar3) = 0
  ifixed(1) = 1

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1220)
  write ( *,1010)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(3, par3, npar3, pv, sdpv, res, sdres, vcv, n3, ivcv, &
     nnzw, npare, rsd)
  call nlsc(y3, xm3, n3, m3, ixm3, mdl3, par3, npar3, res, ldstak, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1220)
  nprt = 11111
  write ( *,1070)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(3, par3, npar3, pv, sdpv, res, sdres, vcv, n3, ivcv, &
     nnzw, npare, rsd)
  call nlsdc(y3, xm3, n3, m3, ixm3, mdl3, drv3, par3, npar3, res, &
     ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
!
!  test xm
!
!  first column zero
!
  call setrv(y2, n2, 2.0e0)
  call setrv(xm2(1,1), n2, 0.0e0)

  nprt = 11000

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1310)
  call nlsx1(2, par2, npar2, pv, sdpv, res, sdres, vcv, n2, ivcv, &
     nnzw, npare, rsd)
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlswds(y2, wt, xm2, n2, m2, ixm2, mdl2, drv2, par2, npar2, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, &
     delta, ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, &
     ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par2, pv, sdpv, res, sdres, vcv, n2, npar2, ivcv, &
     nnzw, npare, rsd)
!
!  2 columns zero
!
  call setrv(xm2(1,2), n2, 0.0e0)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1300)
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(2, par2, npar2, pv, sdpv, res, sdres, vcv, n2, ivcv, &
     nnzw, npare, rsd)
  call nlswds(y2, wt, xm2, n2, m2, ixm2, mdl2, drv2, par2, npar2, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par2, pv, sdpv, res, sdres, vcv, n2, npar2, ivcv, &
     nnzw, npare, rsd)
!
!  test variance-covariance matrix computations
!
  ivctst(1) = -1
  ivctst(2) = 0
  ivctst(3) = 1
  ivctst(4) = 2
  ivctst(5) = 3
  ivctst(6) = 4
  ivctst(7) = 5
  ivctst(8) = 6
  ivctst(9) = 7
  nprt = 2
  do i=1,9
     ntest = ntest + 1
     write ( *,1330) ntest
     write ( *,1380)
     call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, &
        ivcv, nnzw, npare, rsd)
     write ( *,1110)
     write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
        scale(1), delta, ivctst(i), nprt
     call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, &
        npar1, res, ldstak, ifixed, idrvck, mit, stopss, stopp, &
        scale, delta, ivctst(i), nprt, nnzw, npare, rsd, pv, sdpv, &
        sdres, vcv, ivcv)
     write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
        scale(1), delta, ivctst(i), nprt
     write ( *,1120) ierr
     call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
        nnzw, npare, rsd)
  end do
!
!  test with 2 zero weights
!
  nprt = 22222
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1290)
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  wt(3) = 0.0e0
  wt(5) = 0.0e0
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  return

 1000 format ('ldstak must be greater than or equal to ', i6)
 1010 format ('test of nlsc'  )
 1020 format ('test of nlss'  )
 1040 format ('test of nlswc' )
 1050 format ('test of nlsws' )
 1070 format ('test of nlsdc' )
 1080 format ('test of nlsds' )
 1100 format ('test of nlswdc')
 1110 format ('test of nlswds')
 1120 format (/'returned results (-1 indicates ', &
     'value not changed by called subroutine)'//'ierr is ', i3)
 1190 format ('test parameter handling')
 1200 format ('all parameters zero')
 1210 format ('test with constant zero y')
 1220 format ('test linear model')
 1240 format ('test control criteria')
 1250 format ('maximum number of iterations = ', i5)
 1260 format (12h --stopss = , g14.8)
 1270 format (11h --stopp = , g14.8)
 1280 format (10h --nprt = , i6)
 1290 format (29h **test with 2 zero weights**)
 1300 format (19h **2 columns zero**)
 1310 format (18h **1 column zero**)
 1330 format (54h1nonlinear least squares estimation subroutine test nu, &
     4hmber, i5/)
 1340 format (24h input   -  ifixed(1) = , i6, 9x, 11h, stp(1) = , &
     g15.8, 11h,    mit = ,i5, 11h, stopss = , g15.8, 10h, stopp = , &
     g15.8/13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i6)
 1350 format (//24h output  -  ifixed(1) = , i6, 9x, 11h, stp(1) = , &
     g15.8, 11h,    mit = ,i5, 11h, stopss = , g15.8, 10h, stopp = , &
     g15.8/13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i6//)
 1360 format (24h input   -  ifixed(1) = , i6, 9x, 11h, idrvck = , &
     i5, 10x, 11h,    mit = ,i5, 11h, stopss = , g15.8, &
     10h, stopp = , g15.8/ &
     13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i6)
 1370 format (//24h output  -  ifixed(1) = , i6, 9x, 11h, idrvck = , &
     i5, 10x, 11h,    mit = ,i5, 11h, stopss = , g15.8, &
     10h, stopp = , g15.8/ &
     13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i6//)
 1380 format (54h test handling of variance-covariance computation code, &
     's')
end
subroutine xnrand ( ldstak )

!*****************************************************************************80
!
!! XNRAND tests the NRAND family.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson, John Koontz,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index variable
!     integer ierr
!        flag to indicate presence of error detected by preceding
!        starpac call.  (0 is ok, 1 is error)
!     integer iseed
!        the seed for the random number generator.
!     integer ldstak
!        amount of work area.  size of dstak.
!     integer n
!        the length of the vector y.
!     real sigma
!        the s.d. of the sample.
!     real y(1000)
!        data vector for tests.
!     real ymean
!        the mean of the sample.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     sigma,ymean
  integer &
     i,iseed,n
!
!  local arrays
  real &
     y(1000)
!
!  external subroutines
  external hist,iprint,nrand,nrandc
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!  check for sufficient work area length.
!
  if (ldstak.lt.1000) then
    write ( *, 1000)
     return
  end if
!
!  heading.
!
  write ( *,1150)
!
!  test 1.  check all error messages.
!
  write ( *,1160)
!
!  error 1, zero or fewer elements or negative standard deviation.
!
  iseed = 0
  sigma = -1
  n = 0
  ymean = 0.0e0
  call nrand(y, n, iseed)
  write ( *,1170) ierr
  call nrandc(y, n, iseed, ymean, sigma)
  write ( *,1170) ierr
!
!  compare results
!
  iseed = 334
  n = 10
  ymean = 0.0e0
  sigma = 1.0e0

  write ( *, 1120) n, iseed
  call nrand (y, n, iseed)
  write ( *, 1100) (y(i),i=1,n)

  iseed = 333
  write ( *, 1130) n, ymean, sigma, iseed
  call nrandc (y, n, iseed, ymean, sigma)
  write ( *, 1100) (y(i),i=1,n)

  iseed = 13531
  n = 1000
  ymean = 0.0e0
  sigma = 1.0e0

  write ( *, 1120) n, iseed
  call nrand (y, n, iseed)
  call hist (y, n, ldstak)

  write ( *, 1130) n, ymean, sigma, iseed
  call nrandc (y, n, iseed, ymean, sigma)
  call hist (y, n, ldstak)

  iseed = 99999
  n = 1000
  ymean = 4.0e0
  sigma = 4.0e0

  write ( *, 1120) n, iseed
  call nrand (y, n, iseed)
  call hist (y, n, ldstak)

  write ( *, 1130) n, ymean, sigma, iseed
  call nrandc (y, n, iseed, ymean, sigma)
  call hist (y, n, ldstak)

  return

 1000 format ('1the dimension of dstak and the value of ldstak needed'/ &
    ' for nrandx must equal or exceed 1000.  change driver'/ &
    ' and recall nrandx.')
 1100 format (5e15.8)
 1120 format ('1generate ', i4, &
    ' standard normal numbers using iseed = ', i5)
 1130 format ('1generate ', i4, &
    ' normally distributed numbers with ymean = ', f5.2, &
    ' and sigma = ', f5.2, ' using iseed = ', i5)
 1150 format ('1test runs for the nrand family of routines.')
 1160 format(' test 1.  generate each of the possible ', &
     15herror messages.)
 1170 format(/22h the value of ierr is , i4//)
end
subroutine xpgm ( ldstak )

!*****************************************************************************80
!
!! XPGM tests the time series periodogram subroutines.
!
!  Discussion:
!
!    series y is the first 50 values of the series listed on page
!    318 of jenkins and watts.  the spectrum of this series is shown
!    for various bandwidth on page 270 of jenkins and watts.
!
!    series z is the wolf sunspot numbers from 1700 to 1960 as
!    tabulated by waldmeier.  the raw and smoothed periodograms of
!    tapered series are shown on pages 95 and 176, respectively, of
!    bloomfield.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    Gwilym Jenkins, Donald Watts,
!    Spectral Analysis and its Applications,
!    Holden-Day 1968.
!
!    Max Waldmeier,
!    The Sunspot-Activity in the Years 1610-1960,
!    Shulthess, Zurich, 1961.
!
!  Parameters:
!
!     real ab(600)
!        the vector of the nf real and imaginary components of the
!        fourier coefficients.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real freq(300)
!        the vector of frequencies at which the spectrum is computed.
!     integer i
!        an index variable
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list
!        if ierr == 0, no errors were detected
!        if ierr == 1, errors have been detected
!     integer iextnd
!        the indicator variable used to designate whether zero
!        (iextnd == 0) or the series mean (iextnd .ne. 0) is to be
!        used to extend the series.
!     integer itest
!        the number of the test set being performed
!     integer k(10)
!        the vector of the modified daniel filter lengths.
!     integer lab
!        the length of the vector ab.
!     integer ldstak
!        the length of the vector dstak in common cstak.
!     integer lfreq
!        the length of the vector freq.
!     integer lper
!        the length of the vector per.
!     integer lperi
!        the length of the vector peri.
!     integer lzfft
!        the length of the vectors yfft and zfft, respectively..
!     integer nf
!        the number of frequencies at which the spectrum is
!        to be computed.
!     integer nfft
!        the extended series length.
!     integer nk
!        the number of modified daniel filters to be applied.
!     integer nprt
!        the variable controling printed output, where
!        for the periodogram routines
!        if nprt <= -2, the output consists of a page plot of the
!                         periodogram on a log-linear scale,
!        if nprt == -1, the output consists of a page plot of the
!                         periodogram in decibels on a linear scale,
!        if nprt ==  0, the output is suppressed,
!        if nprt ==  1, the output consists of a vertical plot of the
!                         periodogram in decibels on a linear scale.
!        if nprt .ge.  2, the output consists of a vertical plot of the
!                         periodogram on a log-linear scale,
!        and for the integrated periodogram routines
!        if nprt ==  0, the output is suppressed,
!        if nprt .ne.  0, the output consists of a page plot of the
!                         integrated periodogram
!     integer ntemp
!        a temporary storage location
!     integer ny
!        the number of observations in the series y.
!     integer nz
!        the number of observations in the series z.
!     real per(300)
!        the series periodogram.
!     real perf(300)
!        the filtered (smoothed) periodogram.
!     real peri(300)
!        the series integrated periodogram.
!     real taperp
!        the percent of the series to be tapered.
!     real y(150)
!         the array containing the time series from jenkins and watts.
!     real yfft(400)
!        the vector of the observed time series to be analyzed using
!        the fft.
!     real z(275)
!        the array of the wolf sunspot numbers.
!     real zfft(600)
!        the vector of the tapered wolf sunspot numbers to be
!        analyzed using the fft.
!     real zt(275)
!        the array of the tapered sunspot numbers.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     taperp
  integer &
     i,iextnd,itest,lab,lfreq,lper,lperi,lzfft,nf,nfft,nk, &
     nprt,ntemp,ny,nz
!
!  local arrays
  real &
     ab(600),freq(300),per(300),perf(300),peri(300),y(150), &
     yfft(400),z(275),zfft(600),zt(275)
  integer &
     k(10)
!
!  external subroutines
  external center,fftlen,fftr,ipgm,ipgmp,ipgmps,ipgms,iprint,mdflt, &
     pgm,pgms,ppl,scopy,taper
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data   y(  1),  y(  2),  y(  3),  y(  4),  y(  5),  y(  6) &
      / -.88e0,  -.12e0,  -.89e0, -1.38e0,  -.07e0,  1.03e0/
  data   y(  7),  y(  8),  y(  9),  y( 10),  y( 11),  y( 12) &
      / 2.14e0,   .35e0, -1.10e0, -1.78e0, -2.76e0, -1.77e0/
  data   y( 13),  y( 14),  y( 15),  y( 16),  y( 17),  y( 18) &
      /  .98e0,  1.00e0,  -.70e0, -1.01e0, -1.30e0,  -.85e0/
  data   y( 19),  y( 20),  y( 21),  y( 22),  y( 23),  y( 24) &
      / -.46e0,  1.63e0,   .06e0,  -.17e0, -1.01e0, -1.04e0/
  data   y( 25),  y( 26),  y( 27),  y( 28),  y( 29),  y( 30) &
      / -.66e0, -1.12e0,  -.51e0,  -.71e0,  -.20e0,  -.13e0/
  data   y( 31),  y( 32),  y( 33),  y( 34),  y( 35),  y( 36) &
      /  .14e0,  1.59e0,  -.76e0, -1.08e0, -1.77e0, -1.20e0/
  data   y( 37),  y( 38),  y( 39),  y( 40),  y( 41),  y( 42) &
      /  .45e0,  -.07e0,  -.63e0,  -.35e0,  -.87e0,  -.62e0/
  data   y( 43),  y( 44),  y( 45),  y( 46),  y( 47),  y( 48) &
      /  .28e0,  1.90e0,  2.14e0,  1.05e0,   .31e0,  1.07e0/
  data   y( 49),  y( 50) &
      / 2.67e0,  2.44e0/
!
  data   z(  1),  z(  2),  z(  3),  z(  4),  z(  5),  z(  6) &
      /  5.0e0,  11.0e0,  16.0e0,  23.0e0,  36.0e0,  58.0e0/
  data   z(  7),  z(  8),  z(  9),  z( 10),  z( 11),  z( 12) &
      / 29.0e0,  20.0e0,  10.0e0,   8.0e0,   3.0e0,   0.0e0/
  data   z( 13),  z( 14),  z( 15),  z( 16),  z( 17),  z( 18) &
      /  0.0e0,   2.0e0,  11.0e0,  27.0e0,  47.0e0,  63.0e0/
  data   z( 19),  z( 20),  z( 21),  z( 22),  z( 23),  z( 24) &
      / 60.0e0,  39.0e0,  28.0e0,  26.0e0,  22.0e0,  11.0e0/
  data   z( 25),  z( 26),  z( 27),  z( 28),  z( 29),  z( 30) &
      / 21.0e0,  40.0e0,  78.0e0, 122.0e0, 103.0e0,  73.0e0/
  data   z( 31),  z( 32),  z( 33),  z( 34),  z( 35),  z( 36) &
      / 47.0e0,  35.0e0,  11.0e0,   5.0e0,  16.0e0,  34.0e0/
  data   z( 37),  z( 38),  z( 39),  z( 40),  z( 41),  z( 42) &
      / 70.0e0,  81.0e0, 111.0e0, 101.0e0,  73.0e0,  40.0e0/
  data   z( 43),  z( 44),  z( 45),  z( 46),  z( 47),  z( 48) &
      / 20.0e0,  16.0e0,   5.0e0,  11.0e0,  22.0e0,  40.0e0/
  data   z( 49),  z( 50),  z( 51),  z( 52),  z( 53),  z( 54) &
      / 60.0e0,  80.9e0,  83.4e0,  47.7e0,  47.8e0,  30.7e0/
  data   z( 55),  z( 56),  z( 57),  z( 58),  z( 59),  z( 60) &
      / 12.2e0,   9.6e0,  10.2e0,  32.4e0,  47.6e0,  54.0e0/
  data   z( 61),  z( 62),  z( 63),  z( 64),  z( 65),  z( 66) &
      / 62.9e0,  85.9e0,  61.2e0,  45.1e0,  36.4e0,  20.9e0/
  data   z( 67),  z( 68),  z( 69),  z( 70),  z( 71),  z( 72) &
      / 11.4e0,  37.8e0,  69.8e0, 106.1e0, 100.8e0,  81.6e0/
  data   z( 73),  z( 74),  z( 75),  z( 76),  z( 77),  z( 78) &
      / 66.5e0,  34.8e0,  30.6e0,   7.0e0,  19.8e0,  92.5e0/
  data   z( 79),  z( 80),  z( 81),  z( 82),  z( 83),  z( 84) &
      /154.4e0, 125.9e0,  84.8e0,  68.1e0,  38.5e0,  22.8e0/
  data   z( 85),  z( 86),  z( 87),  z( 88),  z( 89),  z( 90) &
      / 10.2e0,  24.1e0,  82.9e0, 132.0e0, 130.9e0, 118.1e0/
  data   z( 91),  z( 92),  z( 93),  z( 94),  z( 95),  z( 96) &
      / 89.9e0,  66.6e0,  60.0e0,  46.9e0,  41.0e0,  21.3e0/
  data   z( 97),  z( 98),  z( 99),  z(100),  z(101),  z(102) &
      / 16.0e0,   6.4e0,   4.1e0,   6.8e0,  14.5e0,  34.0e0/
  data   z(103),  z(104),  z(105),  z(106),  z(107),  z(108) &
      / 45.0e0,  43.1e0,  47.5e0,  42.2e0,  28.1e0,  10.1e0/
  data   z(109),  z(110),  z(111),  z(112),  z(113),  z(114) &
      /  8.1e0,   2.5e0,   0.0e0,   1.4e0,   5.0e0,  12.2e0/
  data   z(115),  z(116),  z(117),  z(118),  z(119),  z(120) &
      / 13.9e0,  35.4e0,  45.8e0,  41.1e0,  30.1e0,  23.9e0/
  data   z(121),  z(122),  z(123),  z(124),  z(125),  z(126) &
      / 15.6e0,   6.6e0,   4.0e0,   1.8e0,   8.5e0,  16.6e0/
  data   z(127),  z(128),  z(129),  z(130),  z(131),  z(132) &
      / 36.3e0,  49.6e0,  64.2e0,  67.0e0,  70.9e0,  47.8e0/
  data   z(133),  z(134),  z(135),  z(136),  z(137),  z(138) &
      / 27.5e0,   8.5e0,  13.2e0,  56.9e0, 121.5e0, 138.3e0/
  data   z(139),  z(140),  z(141),  z(142),  z(143),  z(144) &
      /103.2e0,  85.7e0,  64.6e0,  36.7e0,  24.2e0,  10.7e0/
  data   z(145),  z(146),  z(147),  z(148),  z(149),  z(150) &
      / 15.0e0,  40.1e0,  61.5e0,  98.5e0, 124.7e0,  96.3e0/
  data   z(151),  z(152),  z(153),  z(154),  z(155),  z(156) &
      / 66.6e0,  64.5e0,  54.1e0,  39.0e0,  20.6e0,   6.7e0/
  data   z(157),  z(158),  z(159),  z(160),  z(161),  z(162) &
      /  4.3e0,  22.7e0,  54.8e0,  93.8e0,  95.8e0,  77.2e0/
  data   z(163),  z(164),  z(165),  z(166),  z(167),  z(168) &
      / 59.1e0,  44.0e0,  47.0e0,  30.5e0,  16.3e0,   7.3e0/
  data   z(169),  z(170),  z(171),  z(172),  z(173),  z(174) &
      / 37.6e0,  74.0e0, 139.0e0, 111.2e0, 101.6e0,  66.2e0/
  data   z(175),  z(176),  z(177),  z(178),  z(179),  z(180) &
      / 44.7e0,  17.0e0,  11.3e0,  12.4e0,   3.4e0,   6.0e0/
  data   z(181),  z(182),  z(183),  z(184),  z(185),  z(186) &
      / 32.3e0,  54.3e0,  59.7e0,  63.7e0,  63.5e0,  52.2e0/
  data   z(187),  z(188),  z(189),  z(190),  z(191),  z(192) &
      / 25.4e0,  13.1e0,   6.8e0,   6.3e0,   7.1e0,  35.6e0/
  data   z(193),  z(194),  z(195),  z(196),  z(197),  z(198) &
      / 73.0e0,  85.1e0,  78.0e0,  64.0e0,  41.8e0,  26.2e0/
  data   z(199),  z(200),  z(201),  z(202),  z(203),  z(204) &
      / 26.7e0,  12.1e0,   9.5e0,   2.7e0,   5.0e0,  24.4e0/
  data   z(205),  z(206),  z(207),  z(208),  z(209),  z(210) &
      / 42.0e0,  63.5e0,  53.8e0,  62.0e0,  48.5e0,  43.9e0/
  data   z(211),  z(212),  z(213),  z(214),  z(215),  z(216) &
      / 18.6e0,   5.7e0,   3.6e0,   1.4e0,   9.6e0,  47.4e0/
  data   z(217),  z(218),  z(219),  z(220),  z(221),  z(222) &
      / 57.1e0, 103.9e0,  80.6e0,  63.6e0,  37.6e0,  26.1e0/
  data   z(223),  z(224),  z(225),  z(226),  z(227),  z(228) &
      / 14.2e0,   5.8e0,  16.7e0,  44.3e0,  63.9e0,  69.0e0/
  data   z(229),  z(230),  z(231),  z(232),  z(233),  z(234) &
      / 77.8e0,  64.9e0,  35.7e0,  21.2e0,  11.1e0,   5.7e0/
  data   z(235),  z(236),  z(237),  z(238),  z(239),  z(240) &
      /  8.7e0,  36.1e0,  79.7e0, 114.4e0, 109.6e0,  88.8e0/
  data   z(241),  z(242),  z(243),  z(244),  z(245),  z(246) &
      / 67.8e0,  47.5e0,  30.6e0,  16.3e0,   9.6e0,  33.2e0/
  data   z(247),  z(248),  z(249),  z(250),  z(251),  z(252) &
      / 92.6e0, 151.6e0, 136.3e0, 134.7e0,  83.9e0,  69.4e0/
  data   z(253),  z(254),  z(255),  z(256),  z(257),  z(258) &
      / 31.5e0,  13.9e0,   4.4e0,  38.0e0, 141.7e0, 190.2e0/
  data   z(259),  z(260),  z(261) &
      /184.8e0, 159.0e0, 112.3e0/

  itest = 1
!
!  make calls with valid data
!
  ny = 50
  nz = 261
  nfft = 514
  nprt = 2
  lzfft = 600
  lper = 514
  lperi = 514
  taperp = 0.10e0
  lfreq = 300
  iextnd = 0
  lab = 600
  nk = 3
  k(1) = 8
  k(2) = 8
  k(3) = 8
!
!  test of center
!
    5 write ( *, 1018)
  call center (z, nz, zt)
  write ( *, 1002) ierr
!
!  print returned variables from center
!
  if (ierr==0) write ( *, 1004) (zt(i), i = 1, nz)
!
!  test of taper
!
  write ( *, 1015)
  call taper (z, nz, taperp, zt)
  write ( *, 1002) ierr
!
!  print returned variables from taper
!
  if (ierr==0) write ( *, 1004) (zt(i), i = 1, nz)
!
!  test of pgm
!
  write ( *, 1013)
  call scopy (nz, zt, 1, zfft, 1)
  call pgm (zfft, nz, lzfft, ldstak)
  write ( *, 1002) ierr
!
!  test of fftlen
!
  write ( *, 1026)
  call fftlen(nfft-2, 2, ntemp)
  write ( *, 1002) ierr
!
!  print returned variables from fftlen
!
  if (ierr==0) write ( *, 1027) ntemp
!
!  test of pgms
!
  ntemp = nfft-1
  write ( *, 1025)
  call scopy (nz, zt, 1, zfft, 1)
  call pgms (zfft, nz, ntemp, lzfft, iextnd, nf, per, lper, freq, &
     lfreq, -2)
  write ( *, 1027) ntemp
  write ( *, 1002) ierr
!
!  print returned variables from pgms
!
  if (ierr==0) then
    write ( *, 1004) (freq(i), i = 1, nf)
    write ( *, 1004) (per(i), i = 1, nf)
  end if
!
!  test of mdflt
!
  write ( *, 1016)
  call mdflt (per, nf, nk, k, perf, ldstak)
  write ( *, 1002) ierr
  if (ierr==0) then
!
!  print returned variables from mdflt
!
    write ( *, 1004) (perf(i), i = 1, nf)
!
! display smoothed periodogram on a log plot
!
    write ( *, 1028)
    call ppl (perf, freq, nf, 1)
  end if
!
!  test of ipgmp
!
  write ( *, 1029)
  call ipgmp (per, freq, nf, nz, ldstak)
  write ( *, 1002) ierr
!
!  test of ipgmps
!
  write ( *, 1030)
  call ipgmps (per, freq, nf, nz, ldstak, peri, nprt)
  write ( *, 1002) ierr
  if (ierr==0) then
!
!  print returned variables from ipgmps
!
     write ( *, 1004) (freq(i), i = 1, nf)
     write ( *, 1004) (peri(i), i = 1, nf)
  end if
!
!  test of ipgm
!
  write ( *, 1017)
  call scopy (nz, zt, 1, zfft, 1)
  call ipgm (zfft, nz, lzfft, ldstak)
  write ( *, 1002) ierr
!
!  test of ipgms
!
  write ( *, 1014)
  call scopy (nz, zt, 1, zfft, 1)
  call ipgms (zfft, nz, lzfft, ldstak, nf, peri, lperi, freq, lfreq, &
     nprt)
  write ( *, 1002) ierr
  if (ierr==0) then
!
!  print returned variables from ipgms
!
     write ( *, 1004) (freq(i), i = 1, nf)
     write ( *, 1004) (peri(i), i = 1, nf)
  end if
!
!  test of fftr (centered data - o percent taper)
!
  taperp = -1.0e0
  write ( *, 1031)
  if (ny.ge.1) call taper (y, ny, taperp, yfft)
  call fftr (yfft, ny, nfft, iextnd, nf, ab, lab)
  write ( *, 1002) ierr
!
!  print returned variables from fftr
!
  if (ierr==0) write ( *, 1004) (ab(i), i = 1, nf)
!
  go to (10, 20, 30, 40) itest
!
!  check minimum problem size
!
   10 itest = itest + 1
  nz = 17
  ny = 17
  go to 5
!
!  check various options
!
   20 itest = itest + 1
!
!  test of mdflt (elements of k not even)
!
  k(1) = 7
  write ( *, 1016)
  call mdflt (per, nf, nk, k, perf, ldstak)
  write ( *, 1002) ierr
!
!  print returned variables from mdflt
!
  write ( *, 1004) (perf(i), i = 1, nf)
!
!  test of pgms (uncentered data)
!
  iextnd = 1
  nprt = 1
  write ( *, 1025)
  call scopy (nz, z, 1, zfft, 1)
  call pgms (zfft, nz, nfft, lzfft, iextnd, nf, per, lper, freq, &
     lfreq, nprt)
  write ( *, 1002) ierr
!
!  print returned variables from pgms
!
  write ( *, 1004) (freq(i), i = 1, nf)
  write ( *, 1004) (per(i), i = 1, nf)

  nprt = 2
  write ( *, 1025)
  call scopy (nz, z, 1, zfft, 1)
  call pgms (zfft, nz, nfft, lzfft, iextnd, nf, per, lper, freq, &
     lfreq, nprt)
  write ( *, 1002) ierr
!
!  print returned variables from pgms
!
  write ( *, 1004) (freq(i), i = 1, nf)
  write ( *, 1004) (per(i), i = 1, nf)
!
!  test of fftr (centered data - 100 percent taper)
!
  taperp = 1.1e0
  write ( *, 1031)
  call taper (y, ny, taperp, yfft)
  call fftr (yfft, ny, nfft, iextnd, nf, ab, lab)
  write ( *, 1002) ierr
!
!  print returned variables from fftr
!
  write ( *, 1004) (ab(i), i = 1, nf)
!
!  test of fftr (centered data - o percent taper)
!
  taperp = -1.0e0
  write ( *, 1031)
  call taper (y, ny, taperp, yfft)
  call fftr (yfft, ny, nfft-1, iextnd, nf, ab, nfft-2)
  write ( *, 1002) ierr
!
!  print returned variables from fftr
!
  write ( *, 1004) (ab(i), i = 1, nf)
!
!  perform error checking
!
  ny = -1
  nz = -1
  nfft = 0
  nprt = 2
  lzfft = 0
  lper = 0
  lperi = 0
  taperp = 0.10e0
  lfreq = 0
  iextnd = 0
  lab = 0
  nk = 0
  k(1) = 0
  k(2) = 0
  k(3) = 0
  go to 5
!
!     perform more error checking
!
   30 itest = itest + 1
  ny = 50
  nz = 261
  nprt = 2
  lzfft = 0
  lper = 0
  lperi = 0
  taperp = 0.10e0
  lfreq = 0
  iextnd = 0
  lab = 0
  nk = 3
  k(1) = 0
  k(2) = 0
  k(3) = 0
  go to 5

   40 return

 1002 format ('ierr is ', i5)
 1004 format (3(1x, e16.8))
 1013 format ('test of pgm')
 1014 format ('test of ipgms')
 1015 format ('test of taper')
 1016 format ('test of mdflt')
 1017 format ('test of ipgm')
 1018 format ('test of center')
 1025 format ('test of pgms')
 1026 format ('test of fftlen')
 1027 format (/'nfft is ', i6)
 1028 format ('display of periodogram smoothed with modified', &
     ' daniel filter')
 1029 format ('test of ipgmp')
 1030 format ('test of ipgmps')
 1031 format ('test of fftr')
end
subroutine xpp ( )

!*****************************************************************************80
!
!! XPP tests the plotting subroutines.
!
!  Discussion:
!
!    series y is the airline data listed on page 531 of box and jenkins.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    George Box, Gwilym Jenkins,
!    Time Series Analysis: Forecasting and Control,
!    Holden-Day, 1970,
!    QA280.B67.
!
!  Parameters:
!
!     real air(144)
!        the airline data.
!     integer ierr
!        a common variable used as a flag to indicate whether
!        or not there are any errors, if =0 then no errors.
!     integer ilog
!        the two digit integer, pq, used to select axis scale, where
!        p designates the x-axis and q designates the y-axis.
!        if p==0 (q==0), then the x-axis (y-axis) is linear.
!        if p.ne.0 (q.ne.0), then the x-axis (y-axis) is log.
!     integer isize
!        the two digit integer, pq, used to select axis size, where
!        p designates the x-axis and q designates the y-axis.
!        if p==0 (q==0), then the x-axis (y-axis) is the maximum.
!        if p.ne.0 (q.ne.0), then the x-axis (y-axis) is half the maximu
!     integer isym(144)
!        vector containing symbol designations for plotting
!     integer itest
!        the number of the test.
!     integer iym
!        actual dimension of ym in users main program
!     integer m
!        the number of vectors in ym
!     integer nout
!        used to indicate how many of the points outside the bounds
!        of the plot are to be listed.
!     integer ny, nym
!        the number of observations in arrays y and ym, respectively.
!     real time(144)
!        the time values for the airline data.
!     real x(144)
!        vector of observations for x(horizontal) coordinates
!     real xlb
!        the lower bound for the x-axis.  (xlb=xub indicates limits are
!        to be determined from the range of the data.)
!     real xmiss
!        the missing value code for the x-axis.
!     real xub
!        the upper bound for the x-axis.  (xlb=xub indicates limits are
!        to be determined from the range of the data.)
!     real y(144)
!        vector of observations for the y (vertical) coordinates
!     real ylb
!        the lower bound for the y-axis.  (ylb=yub indicates limits are
!        to be determined from the range of the data.)
!     real ym(12,12)
!        multivariate observations for the y (vertical) coordinates.
!     real ymiss
!        the missing value code for the y-axis.
!     real ymmiss(144)
!        the missing value codes for each column of ym.
!     real yub
!        the upper bound for the y-axis.  (ylb=yub indicates limits are
!        to be determined from the range of the data.)
!
  implicit none
!
!  scalars in common
  integer &
     ierr
!
!  local scalars
  real &
     xlb,xmiss,xub,ylb,ymiss,yub
  integer &
     ilog,isize,itest,iym,m,nout,ny,nym
!
!  local arrays
  real &
     air(144),time(144),x(144),y(144),ym(12,12),ymmiss(144)
  integer &
     isym(144)
!
!  external subroutines
  external iprint,mpp,mppc,mppl,mppm,mppmc,mppml,pp,ppc,ppl,ppm, &
     ppmc,ppml,scopy,setrv,spp,sppc,sppl,sppm,sppmc,sppml
!
!  common blocks
  common /errchk/ierr
!
!  equivalences
  equivalence (y(1),ym(1,1))

  data     xmiss,    ymiss &
      /      7.0e0,    180.0e0/

  data isym(  1),isym(  2),isym(  3),isym(  4),isym(  5),isym(  6) &
      /    -5000,     6000,        7,        8,        9,       10/
  data isym(  7),isym(  8),isym(  9),isym( 10),isym( 11),isym( 12) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 13),isym( 14),isym( 15),isym( 16),isym( 17),isym( 18) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 19),isym( 20),isym( 21),isym( 22),isym( 23),isym( 24) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 25),isym( 26),isym( 27),isym( 28),isym( 29),isym( 30) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 31),isym( 32),isym( 33),isym( 34),isym( 35),isym( 36) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 37),isym( 38),isym( 39),isym( 40),isym( 41),isym( 42) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 43),isym( 44),isym( 45),isym( 46),isym( 47),isym( 48) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 49),isym( 50),isym( 51),isym( 52),isym( 53),isym( 54) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 55),isym( 56),isym( 57),isym( 58),isym( 59),isym( 60) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 61),isym( 62),isym( 63),isym( 64),isym( 65),isym( 66) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 67),isym( 68),isym( 69),isym( 70),isym( 71),isym( 72) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 73),isym( 74),isym( 75),isym( 76),isym( 77),isym( 78) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 79),isym( 80),isym( 81),isym( 82),isym( 83),isym( 84) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 85),isym( 86),isym( 87),isym( 88),isym( 89),isym( 90) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 91),isym( 92),isym( 93),isym( 94),isym( 95),isym( 96) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 97),isym( 98),isym( 99),isym(100),isym(101),isym(102) &
      /        5,        6,        7,        8,        9,       10/
  data isym(103),isym(104),isym(105),isym(106),isym(107),isym(108) &
      /       11,       12,       13,       14,       15,       16/
  data isym(109),isym(110),isym(111),isym(112),isym(113),isym(114) &
      /        5,        6,        7,        8,        9,       10/
  data isym(115),isym(116),isym(117),isym(118),isym(119),isym(120) &
      /       11,       12,       13,       14,       15,       16/
  data isym(121),isym(122),isym(123),isym(124),isym(125),isym(126) &
      /        5,        6,        7,        8,        9,       10/
  data isym(127),isym(128),isym(129),isym(130),isym(131),isym(132) &
      /       11,       12,       13,       14,       15,       16/
  data isym(133),isym(134),isym(135),isym(136),isym(137),isym(138) &
      /        5,        6,        7,        8,        9,       10/
  data isym(139),isym(140),isym(141),isym(142),isym(143),isym(144) &
      /       11,       12,       13,       14,       15,       16/
!
  data time(  1),time(  2),time(  3),time(  4),time(  5),time(  6) &
      /   1.0e0,    2.0e0,    3.0e0,    4.0e0,    5.0e0,    6.0e0/
  data time(  7),time(  8),time(  9),time( 10),time( 11),time( 12) &
      /   7.0e0,    8.0e0,    9.0e0,   10.0e0,   11.0e0,   12.0e0/
  data time( 13),time( 14),time( 15),time( 16),time( 17),time( 18) &
      /  13.0e0,   14.0e0,   15.0e0,   16.0e0,   17.0e0,   18.0e0/
  data time( 19),time( 20),time( 21),time( 22),time( 23),time( 24) &
      /  19.0e0,   20.0e0,   21.0e0,   22.0e0,   23.0e0,   24.0e0/
  data time( 25),time( 26),time( 27),time( 28),time( 29),time( 30) &
      /  25.0e0,   26.0e0,   27.0e0,   28.0e0,   29.0e0,   30.0e0/
  data time( 31),time( 32),time( 33),time( 34),time( 35),time( 36) &
      /  31.0e0,   32.0e0,   33.0e0,   34.0e0,   35.0e0,   36.0e0/
  data time( 37),time( 38),time( 39),time( 40),time( 41),time( 42) &
      /  37.0e0,   38.0e0,   39.0e0,   40.0e0,   41.0e0,   42.0e0/
  data time( 43),time( 44),time( 45),time( 46),time( 47),time( 48) &
      /  43.0e0,   44.0e0,   45.0e0,   46.0e0,   47.0e0,   48.0e0/
  data time( 49),time( 50),time( 51),time( 52),time( 53),time( 54) &
      /  49.0e0,   50.0e0,   51.0e0,   52.0e0,   53.0e0,   54.0e0/
  data time( 55),time( 56),time( 57),time( 58),time( 59),time( 60) &
      /  55.0e0,   56.0e0,   57.0e0,   58.0e0,   59.0e0,   60.0e0/
  data time( 61),time( 62),time( 63),time( 64),time( 65),time( 66) &
      /  61.0e0,   62.0e0,   63.0e0,   64.0e0,   65.0e0,   66.0e0/
  data time( 67),time( 68),time( 69),time( 70),time( 71),time( 72) &
      /  67.0e0,   68.0e0,   69.0e0,   70.0e0,   71.0e0,   72.0e0/
  data time( 73),time( 74),time( 75),time( 76),time( 77),time( 78) &
      /  73.0e0,   74.0e0,   75.0e0,   76.0e0,   77.0e0,   78.0e0/
  data time( 79),time( 80),time( 81),time( 82),time( 83),time( 84) &
      /  79.0e0,   80.0e0,   81.0e0,   82.0e0,   83.0e0,   84.0e0/
  data time( 85),time( 86),time( 87),time( 88),time( 89),time( 90) &
      /  85.0e0,   86.0e0,   87.0e0,   88.0e0,   89.0e0,   90.0e0/
  data time( 91),time( 92),time( 93),time( 94),time( 95),time( 96) &
      /  91.0e0,   92.0e0,   93.0e0,   94.0e0,   95.0e0,   96.0e0/
  data time( 97),time( 98),time( 99),time(100),time(101),time(102) &
      /  97.0e0,   98.0e0,   99.0e0,  100.0e0,  101.0e0,  102.0e0/
  data time(103),time(104),time(105),time(106),time(107),time(108) &
      / 103.0e0,  104.0e0,  105.0e0,  106.0e0,  107.0e0,  108.0e0/
  data time(109),time(110),time(111),time(112),time(113),time(114) &
      / 109.0e0,  110.0e0,  111.0e0,  112.0e0,  113.0e0,  114.0e0/
  data time(115),time(116),time(117),time(118),time(119),time(120) &
      / 115.0e0,  116.0e0,  117.0e0,  118.0e0,  119.0e0,  120.0e0/
  data time(121),time(122),time(123),time(124),time(125),time(126) &
      / 121.0e0,  122.0e0,  123.0e0,  124.0e0,  125.0e0,  126.0e0/
  data time(127),time(128),time(129),time(130),time(131),time(132) &
      / 127.0e0,  128.0e0,  129.0e0,  130.0e0,  131.0e0,  132.0e0/
  data time(133),time(134),time(135),time(136),time(137),time(138) &
      / 133.0e0,  134.0e0,  135.0e0,  136.0e0,  137.0e0,  138.0e0/
  data time(139),time(140),time(141),time(142),time(143),time(144) &
      / 139.0e0,  140.0e0,  141.0e0,  142.0e0,  143.0e0,  144.0e0/
!
  data  air(  1), air(  2), air(  3), air(  4), air(  5), air(  6) &
      / 112.0e0,  118.0e0,  132.0e0,  129.0e0,  121.0e0,  135.0e0/
  data  air(  7), air(  8), air(  9), air( 10), air( 11), air( 12) &
      / 148.0e0,  148.0e0,  136.0e0,  119.0e0,  104.0e0,  118.0e0/
  data  air( 13), air( 14), air( 15), air( 16), air( 17), air( 18) &
      / 115.0e0,  126.0e0,  141.0e0,  135.0e0,  125.0e0,  149.0e0/
  data  air( 19), air( 20), air( 21), air( 22), air( 23), air( 24) &
      / 170.0e0,  170.0e0,  158.0e0,  133.0e0,  114.0e0,  140.0e0/
  data  air( 25), air( 26), air( 27), air( 28), air( 29), air( 30) &
      / 145.0e0,  150.0e0,  178.0e0,  163.0e0,  172.0e0,  178.0e0/
  data  air( 31), air( 32), air( 33), air( 34), air( 35), air( 36) &
      / 199.0e0,  199.0e0,  184.0e0,  162.0e0,  146.0e0,  166.0e0/
  data  air( 37), air( 38), air( 39), air( 40), air( 41), air( 42) &
      / 171.0e0,  180.0e0,  193.0e0,  181.0e0,  183.0e0,  218.0e0/
  data  air( 43), air( 44), air( 45), air( 46), air( 47), air( 48) &
      / 230.0e0,  242.0e0,  209.0e0,  191.0e0,  172.0e0,  194.0e0/
  data  air( 49), air( 50), air( 51), air( 52), air( 53), air( 54) &
      / 196.0e0,  196.0e0,  236.0e0,  235.0e0,  229.0e0,  243.0e0/
  data  air( 55), air( 56), air( 57), air( 58), air( 59), air( 60) &
      / 264.0e0,  272.0e0,  237.0e0,  211.0e0,  180.0e0,  201.0e0/
  data  air( 61), air( 62), air( 63), air( 64), air( 65), air( 66) &
      / 204.0e0,  188.0e0,  235.0e0,  227.0e0,  234.0e0,  264.0e0/
  data  air( 67), air( 68), air( 69), air( 70), air( 71), air( 72) &
      / 302.0e0,  293.0e0,  259.0e0,  229.0e0,  203.0e0,  229.0e0/
  data  air( 73), air( 74), air( 75), air( 76), air( 77), air( 78) &
      / 242.0e0,  233.0e0,  267.0e0,  269.0e0,  270.0e0,  315.0e0/
  data  air( 79), air( 80), air( 81), air( 82), air( 83), air( 84) &
      / 364.0e0,  347.0e0,  312.0e0,  274.0e0,  237.0e0,  278.0e0/
  data  air( 85), air( 86), air( 87), air( 88), air( 89), air( 90) &
      / 284.0e0,  277.0e0,  317.0e0,  313.0e0,  318.0e0,  374.0e0/
  data  air( 91), air( 92), air( 93), air( 94), air( 95), air( 96) &
      / 413.0e0,  405.0e0,  355.0e0,  306.0e0,  271.0e0,  306.0e0/
  data  air( 97), air( 98), air( 99), air(100), air(101), air(102) &
      / 315.0e0,  301.0e0,  356.0e0,  348.0e0,  355.0e0,  422.0e0/
  data  air(103), air(104), air(105), air(106), air(107), air(108) &
      / 465.0e0,  467.0e0,  404.0e0,  347.0e0,  305.0e0,  336.0e0/
  data  air(109), air(110), air(111), air(112), air(113), air(114) &
      / 340.0e0,  318.0e0,  362.0e0,  348.0e0,  363.0e0,  435.0e0/
  data  air(115), air(116), air(117), air(118), air(119), air(120) &
      / 491.0e0,  505.0e0,  404.0e0,  359.0e0,  310.0e0,  337.0e0/
  data  air(121), air(122), air(123), air(124), air(125), air(126) &
      / 360.0e0,  342.0e0,  406.0e0,  396.0e0,  420.0e0,  472.0e0/
  data  air(127), air(128), air(129), air(130), air(131), air(132) &
      / 548.0e0,  559.0e0,  463.0e0,  407.0e0,  362.0e0,  405.0e0/
  data  air(133), air(134), air(135), air(136), air(137), air(138) &
      / 417.0e0,  391.0e0,  419.0e0,  461.0e0,  472.0e0,  535.0e0/
  data  air(139), air(140), air(141), air(142), air(143), air(144) &
      / 622.0e0,  606.0e0,  508.0e0,  461.0e0,  390.0e0,  432.0e0/

  call setrv(ymmiss, 144, ymiss)
  call scopy(144, air, 1, y, 1)
  call scopy(144, time, 1, x, 1)
!
!     define constants
!

  itest = 0
!
!  short calls
!
  ny = 144
  nym = 12
  iym = 12
  m = 12
  ilog = -1
  isize = -1
  nout = -1
  ylb = 0.0e0
  yub = 0.0e0
  xlb = 0.0e0
  xub = 0.0e0

   10 continue
!
!  test of pp
!
  write ( *, 1000)
  write ( *, 3100) itest
  write ( *, 3010) ny
  call pp(y, x, ny)
  write ( *, 3000) ierr
!
!  test of ppm
!
  write ( *, 1030)
  write ( *, 3100) itest
  write ( *, 3010) ny
  call ppm(y, ymiss, x, xmiss, ny)
  write ( *, 3000) ierr
!
!  test of spp
!
  write ( *, 1120)
  write ( *, 3100) itest
  write ( *, 3010) ny
  call spp(y, x, ny, isym)
  write ( *, 3000) ierr
!
!  test of sppm
!
  write ( *, 1150)
  write ( *, 3100) itest
  write ( *, 3010) ny
  call sppm(y, ymiss, x, xmiss, ny, isym)
  write ( *, 3000) ierr
!
!  test of mpp
!
  write ( *, 1060)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  call mpp(ym, x, nym, m, iym)
  write ( *, 3000) ierr
!
!  test of mppm
!
  write ( *, 1090)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  call mppm(ym, ymmiss, x, xmiss, nym, m, iym)
  write ( *, 3000) ierr
!
!  log option calls
!
   20 continue
!
!  test of ppl
!
  write ( *, 1010)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  call ppl(y, x, ny, ilog)

  write ( *, 3000) ierr
!
!  test of ppml
!
  write ( *, 1040)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  call ppml(y, ymiss, x, xmiss, ny, ilog)
  write ( *, 3000) ierr
!
!  test of sppl
!
  write ( *, 1130)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  call sppl(y, x, ny, isym, ilog)
  write ( *, 3000) ierr
!
!  test of sppml
!
  write ( *, 1160)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  call sppml(y, ymiss, x, xmiss, ny, isym, ilog)
  write ( *, 3000) ierr
!
!  test of mppl
!
  write ( *, 1070)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3040) ilog
  call mppl(ym, x, nym, m, iym, ilog)
  write ( *, 3000) ierr
!
!  test of mppml
!
  write ( *, 1100)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3040) ilog
  call mppml(ym, ymmiss, x, xmiss, nym, m, iym, ilog)
  write ( *, 3000) ierr
!
!  test of long calls
!
   30 continue
!
!  test of ppc
!
  write ( *, 1020)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  write ( *, 3050) isize, nout
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3080) xub
  call ppc(y, x, ny, ilog, isize, nout, ylb, yub, xlb, xub)
  write ( *, 3000) ierr
!
!  test of ppmc
!
  write ( *, 1050)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  write ( *, 3050) isize, nout
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3080) xub
  call ppmc(y, ymiss, x, xmiss, ny, ilog, isize, nout, ylb, yub, &
     xlb, xub)
  write ( *, 3000) ierr
!
!  test of sppc
!
  write ( *, 1140)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  write ( *, 3050) isize, nout
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3080) xub
  call sppc(y, x, ny, isym, ilog, isize, nout, ylb, yub, xlb, &
     xub)
  write ( *, 3000) ierr
!
!  test of sppmc
!
  write ( *, 1170)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  write ( *, 3050) isize, nout
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3080) xub
  call sppmc(y, ymiss, x, xmiss, ny, isym, ilog, isize, nout, &
     ylb, yub, xlb, xub)
  write ( *, 3000) ierr
!
!  test of mppc
!
   40 write ( *, 1080)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3040) ilog
  write ( *, 3050) isize, nout
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3080) xub
  call mppc(ym, x, nym, m, iym, ilog, isize, nout, ylb, yub, &
     xlb, xub)
  write ( *, 3000) ierr
!
!  test of mppmc
!
   50 write ( *, 1110)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3040) ilog
  write ( *, 3050) isize, nout
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3080) xub
  call mppmc(ym, ymmiss, x, xmiss, nym, m, iym, ilog, isize, nout, &
     ylb, yub, xlb, xub)
  write ( *, 3000) ierr

  itest = itest + 1

  go to (110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300), &
    itest
!
!  test valid options
!
  110 ilog = 0
  isize = 0
  nout = 0
  ylb = 100.0e0
  yub = 700.0e0
  xlb = 4.0e0
  xub = 16.0e0
  go to 20

  120 ilog = 2
  isize = 2
  nout = 5
  go to 20

  130 ilog = 20
  isize = 20
  nout = 55
  yub = 300.0e0
  go to 30

  140 ilog = 22
  isize = 22
  go to 40

  150 ny = 1
  nym = 1
  m = 144
  iym = 1
  x(1) = 10.0e0
  go to 40
!
  160 call setrv(y, 144, 1.0e0)
  call setrv(x, 144, 1.0e0)
  nym = 6
  iym = 12
  m = 6
  ny = 36
  ylb = 0.0e0
  yub = 0.0e0
  xlb = 0.0e0
  xub = 0.0e0
  go to 30
!
!     test error response
!
  170 ny = 0
  nym = 0
  m = 0
  iym = -1
  go to 10

  180 ny = 144
  nym = 12
  m = 12
  iym = -1
  xlb = -1.0e0
  ylb = -1.0e0
  go to 40

  190 iym = 12
  x(1) = 0.0e0
  y(1) = 0.0e0
  go to 50

  200 call setrv(x, 144, xmiss)
  call setrv(y, 144, ymiss)
  xlb = xub
  ylb = yub
  go to 50

  300 continue

  return

 1000 format ('1', 10htest of pp)
 1010 format ('1', 11htest of ppl)
 1020 format ('1', 11htest of ppc)
 1030 format ('1', 11htest of ppm)
 1040 format ('1', 12htest of ppml)
 1050 format ('1', 12htest of ppmc)
 1060 format ('1', 11htest of mpp)
 1070 format ('1', 12htest of mppl)
 1080 format ('1', 12htest of mppc)
 1090 format ('1', 12htest of mppm)
 1100 format ('1', 13htest of mppml)
 1110 format ('1', 13htest of mppmc)
 1120 format ('1', 11htest of spp)
 1130 format ('1', 12htest of sppl)
 1140 format ('1', 12htest of sppc)
 1150 format ('1', 12htest of sppm)
 1160 format ('1', 13htest of sppml)
 1170 format ('1', 13htest of sppmc)
 3000 format (/8h ierr = , i4)
 3010 format (' ', 5x, 10h   n     =, i5)
 3020 format ('+', 20x, 10h / m     =, i5, 10h / iym   =, i5)
 3040 format ('+', 65x, 10h / ilog  =, i5)
 3050 format (' ',  5x, 10h   isize =, i5, 10h / nout  =, i5)
 3070 format ('+', 50x, 10h / ylb   =, f10.4, 10h / yub   =, f10.4, &
     10h / xlb   =, f10.4)
 3080 format ('+', 110x, 10h / xub   =, f10.4)
 3100 format (' ', 13h test number , i5)
end
subroutine xstat ( ldstak )

!*****************************************************************************80
!
!! XSTAT tests the STAT family.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson, John Koontz,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fplm
!        the floating point largest magnitude.
!     integer i
!        a loop index.
!     integer ierr
!        flag to indicate presence of error detected by preceding
!        starpac call.  (0 is ok, 1 is error)
!     integer ldstak
!        amount of work area.  size of dstak.
!     integer n
!        the length of the vector y.
!     integer nconst
!        length of the vector yconst.
!     integer nprtof
!        flag for no output (except error messages).
!     integer nprton
!        flag for full printout.
!     real sts(53)
!        vector of statistics.
!     real wt(84)
!        weights vector.
!     real wtall0(10)
!        n vector of 0 weights.
!     real wtall1(84)
!        n vector of 1 weights.
!     real wtemp
!        temporary storage for one of the weights.
!     real y(84)
!        data vector for tests.
!     real yconst(10)
!        vector of constant data.
!     real ypath(10)
!        a vector of y values designed to force different paths
!        through the summation routines.
!     real ytempn, ytemp1
!        temporary storage for the first and last y value.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     fplm,wtemp,ytemp1,ytempn
  integer &
     i,n,nconst,nprtof,nprton
!
!  local arrays
  real &
     sts(53),wt(84),wtall0(10),wtall1(84),y(84),yconst(10), &
     ypath(10)
!
!  external functions
  real &
     r1mach
  external r1mach
!
!  external subroutines
  external iprint,stat,stats,statw,statws
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

!
!     data initializations.
!
  data n /84/
  data nconst /10/
  data nprton /1/
  data nprtof /0/
!
!     davis-harrison r.h. data, pikes peak.
!
!     this is an arbitrarily chosen data set.
!
  data y( 1), y( 2), y( 3), y( 4) &
      / 0.6067e0, 0.6087e0, 0.6086e0, 0.6134e0/
  data y( 5), y( 6), y( 7) &
      / 0.6108e0, 0.6138e0, 0.6125e0/
  data y( 8), y( 9), y(10), y(11) &
      / 0.6122e0, 0.6110e0, 0.6104e0, 0.7213e0/
  data y(12), y(13), y(14) &
      / 0.7078e0, 0.7021e0, 0.7004e0/
  data y(15), y(16), y(17), y(18) &
      / 0.6981e0, 0.7242e0, 0.7268e0, 0.7418e0/
  data y(19), y(20), y(21) &
      / 0.7407e0, 0.7199e0, 0.6225e0/
  data y(22), y(23), y(24), y(25) &
      / 0.6254e0, 0.6252e0, 0.6267e0, 0.6218e0/
  data y(26), y(27), y(28) &
      / 0.6178e0, 0.6216e0, 0.6192e0/
  data y(29), y(30), y(31), y(32) &
      / 0.6191e0, 0.6250e0, 0.6188e0, 0.6233e0/
  data y(33), y(34), y(35) &
      / 0.6225e0, 0.6204e0, 0.6207e0/
  data y(36), y(37), y(38), y(39) &
      / 0.6168e0, 0.6141e0, 0.6291e0, 0.6231e0/
  data y(40), y(41), y(42) &
      / 0.6222e0, 0.6252e0, 0.6308e0/
  data y(43), y(44), y(45), y(46) &
      / 0.6376e0, 0.6330e0, 0.6303e0, 0.6301e0/
  data y(47), y(48), y(49) &
      / 0.6390e0, 0.6423e0, 0.6300e0/
  data y(50), y(51), y(52), y(53) &
      / 0.6260e0, 0.6292e0, 0.6298e0, 0.6290e0/
  data y(54), y(55), y(56) &
      / 0.6262e0, 0.5952e0, 0.5951e0/
  data y(57), y(58), y(59), y(60) &
      / 0.6314e0, 0.6440e0, 0.6439e0, 0.6326e0/
  data y(61), y(62), y(63) &
      / 0.6392e0, 0.6417e0, 0.6412e0/
  data y(64), y(65), y(66), y(67) &
      / 0.6530e0, 0.6411e0, 0.6355e0, 0.6344e0/
  data y(68), y(69), y(70) &
      / 0.6623e0, 0.6276e0, 0.6307e0/
  data y(71), y(72), y(73), y(74) &
      / 0.6354e0, 0.6197e0, 0.6153e0, 0.6340e0/
  data y(75), y(76), y(77) &
      / 0.6338e0, 0.6284e0, 0.6162e0/
  data y(78), y(79), y(80), y(81) &
      / 0.6252e0, 0.6349e0, 0.6344e0, 0.6361e0/
  data y(82), y(83), y(84) &
      / 0.6373e0, 0.6337e0, 0.6383e0/
  data wt( 1), wt( 2), wt( 3), wt( 4), wt( 5), wt( 6), wt( 7) &
     / 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0/
  data wt( 8), wt( 9), wt(10), wt(11), wt(12), wt(13), wt(14) &
     / 0.5e0, 0.5e0, 0.5e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0/
  data wt(15), wt(16), wt(17), wt(18), wt(19), wt(20), wt(21) &
     / 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.5e0/
  data wt(22), wt(23), wt(24), wt(25), wt(26), wt(27), wt(28) &
     / 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0/
  data wt(29), wt(30), wt(31), wt(32), wt(33), wt(34), wt(35) &
     / 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0/
  data wt(36), wt(37), wt(38), wt(39), wt(40), wt(41), wt(42) &
     / 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0/
  data wt(43), wt(44), wt(45), wt(46), wt(47), wt(48), wt(49) &
     / 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0/
  data wt(50), wt(51), wt(52), wt(53), wt(54), wt(55), wt(56) &
     / 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 0.0e0, 0.0e0/
  data wt(57), wt(58), wt(59), wt(60), wt(61), wt(62), wt(63) &
     / 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0/
  data wt(64), wt(65), wt(66), wt(67), wt(68), wt(69), wt(70) &
     / 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0/
  data wt(71), wt(72), wt(73), wt(74), wt(75), wt(76), wt(77) &
     / 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0/
  data wt(78), wt(79), wt(80), wt(81), wt(82), wt(83), wt(84) &
     / 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0/

  fplm = r1mach(2)
!
!  set up the weights vectors.
!
  wtall1(1:n) = 1.0e0

  yconst(1:nconst) = 1.0e0
  wtall0(1:nconst) = 0.0e0
!
!  heading.
!
  write ( *,1150)
!
!  test 1.  check all error messages.
!
!  error 1, two or fewer elements.
!
  write ( *,1180)
  write ( *,1230)
  write ( *,1240)
  call stat(y, 2, ldstak)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1250)
  call stats(y, 2, ldstak, sts, nprton)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1400)
  call statw(y, wt, 2, ldstak)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1410)
  call statws(y, wt, 2, ldstak, sts, nprton)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
!
!  error 2, not enough space in cstak.
!
  write ( *,1190)
  write ( *,1230)
  write ( *,1240)
  call stat(y, n, n/4)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1250)
  call stats(y, n, n/4, sts, nprton)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1400)
  call statw(y, wt, n, n/4)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1410)
  call statws(y, wt, n, n/4, sts, nprton)
  write ( *,1170) ierr
!
!  error 4, negative weights.
!
  write ( *,1210)
  wtemp = wt(2)
  wt(2) = -1.0e0
  write ( *,1230)
  write ( *,1400)
  call statw(y, wt, n, ldstak)
  write ( *,1390) y(1:10)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1410)
  call statws(y, wt, n, ldstak, sts, nprton)
  write ( *,1390) y(1:10)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
  wt(2) = wtemp
!
!     error 5, all weights zero (plus constant y).
!
  write ( *,1220)
  write ( *,1230)
  write ( *,1400)
  call statw(yconst, wtall0, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1410)
  call statws(yconst, wtall0, nconst, ldstak, sts, nprton)
  write ( *,1170) ierr
!
!  test 2.  check for reading outside of data array.
!
  write ( *,1160)
  ytemp1 = yconst(1)
  yconst(1) = fplm
  ytempn = yconst(nconst)
  yconst(nconst) = fplm
  write ( *,1440)
  write ( *,1240)
  call stat(yconst(2), nconst-2, ldstak)
  write ( *,1170) ierr
  write ( *,1440)
  write ( *,1250)
  call stats(yconst(2), nconst-2, ldstak, sts, nprton)
  write ( *,1170) ierr
  write ( *,1440)
  write ( *,1400)
  call statw(yconst(2), wt, nconst-2, ldstak)
  write ( *,1170) ierr
  write ( *,1440)
  write ( *,1410)
  call statws(yconst(2), wt, nconst-2, ldstak, sts, nprton)
  write ( *,1170) ierr
  yconst(1) = ytemp1
  yconst(nconst) = ytempn
!
!  test 3.  constant y.
!
  write ( *,1200)
  write ( *,1440)
  write ( *,1240)
  call stat(yconst, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1440)
  write ( *,1250)
  call stats(yconst, nconst, ldstak, sts, nprton)
  write ( *,1170) ierr
  write ( *,1440)
  write ( *,1400)
  call statw(yconst, wt, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1440)
  write ( *,1410)
  call statws(yconst, wt, nconst, ldstak, sts, nprton)
  write ( *,1170) ierr
!
!  test 4.  see if turning off the printout works.
!
  write ( *,1260)
  write ( *,1270)
  write ( *,1230)
  write ( *,1250)
  call stats(y, n, ldstak, sts, nprtof)
  write ( *,1390) y(1:10)
  write ( *,1280)
  write ( *,1230)
  write ( *,1410)
  call statws(y, wt, n, ldstak, sts, nprtof)
  write ( *,1390) y(1:10)
!
!     test 5.  make a working run of each routine  first with
!              n=2 (the minimun valid value) and then for the whole
!              data set to check the output.
!
  write ( *,1300)
  write ( *,1310)
  write ( *,1240)
  call stat(y, 3, ldstak)
  write ( *,1310)
  write ( *,1240)
  call stat(y, n, ldstak)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr

  write ( *,1320)
  write ( *,1400)
  call statw(y, wt, 3, ldstak)
  write ( *,1320)
  write ( *,1400)
  call statw(y, wt, n, ldstak)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr

  write ( *,1340)
  write ( *,1250)
  call stats(y, 3, ldstak, sts, nprton)
  write ( *,1340)
  write ( *,1250)
  call stats(y, n, ldstak, sts, nprton)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr

  write ( *,1350)
  write ( *,1410)
  call statws(y, wt, 3, ldstak, sts, nprton)
  write ( *,1350)
  write ( *,1410)
  call statws(y, wt, n, ldstak, sts, nprton)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
!
!  test 5.  check results of weighting all observations
!  with 1.0e0.  compare with stat execution.
!
  write ( *,1370)
  write ( *,1400)
  call statw(y, wtall1, n, ldstak)
  write ( *,1170) ierr
!
!     test 6.  check results of forcing difference paths through
!              the summation routines, using small, simple data sets.
!
  write ( *,1000)
!
!  run data set 6.1
!
  do i=1,10
     ypath(i) = i
  end do
  write ( *,1010)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1020)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!  run data set 6.2
!
  do i=1,10
     ypath(i) = -i
  end do
  write ( *,1030)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1040)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!     run data set 6.3
!
  do i=1,10
     ypath(i) = i-1
  end do
  write ( *,1050)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1060)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!     run data set 6.4
!
  do i=1,10
     ypath(i) = 1-i
  end do
  write ( *,1070)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1080)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!     run data set 6.5
!
  do i=1,10
     ypath(i) = i-6
  end do
  write ( *,1090)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1100)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!     run data set 6.6
!
  do i=1,10
     ypath(i) = i-5
  end do
  write ( *,1110)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1120)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!     run data set 6.7
!
  ypath(1:10) = 0.0e0
  ypath(1) = -5.0e0
  ypath(10) = 5.0e0
  write ( *,1130)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1140)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!     run data set 6.8
!
  ypath(1:10) = 0.0e0
  ypath(1) = -5.0e0
  wtall1(1) = 0.0e0
  ypath(10) = 5.0e0
  wtall1(10) = 0.0e0
  write ( *,1380)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
  return

 1000 format('test 6.  try different paths through the summation code.')
 1010 format('1run stat on 1, ..., 10.')
 1020 format('1run statw on 1, ..., 10.  weights are all 1.')
 1030 format('1run stat on -1, ..., -10.')
 1040 format('1run statw on -1, ..., -10.  weights are all 1.')
 1050 format('1run stat on 0, ..., 9.')
 1060 format('1run statw on 0, ..., 9.  weights are all 1.')
 1070 format('1run stat on 0, ..., -9.')
 1080 format('1run statw on 0, ..., -9.  weights are all 1.')
 1090 format('1stat on -5, ..., 4.')
 1100 format('1run statw on -5, ..., 4.  weights are all 1.')
 1110 format('1run stat on -4, ..., 5.')
 1120 format('1run statw on -4, ..., 5.  weights are all 1.')
 1130 format('1run stat on -1, 8*0, 1.')
 1140 format('1run statw on -1, 8*0, 1.  weights are all 1.')
 1150 format('1test runs for the statistical analysis family routines.')
 1160 format('1test runs to be sure code is not reading outside', &
         ' data array.')
 1170 format(/' the value of ierr is ', i4)
 1180 format('1try two or fewer elements.')
 1190 format('1try insufficient work area.')
 1200 format('1try constant y.')
 1210 format('1try negative weights.')
 1220 format('1try all weights zero (and constant y).')
 1230 format (///)
 1240 format (' call to stat')
 1250 format (' call to stats')
 1260 format(45h1test3.  try turning off the print for those , &
     24hroutines which allow it.)
 1270 format(37h try turning the print off for stats.)
 1280 format(38h try turning the print off for statws.)
 1300 format(52h1test 4.  make working runs of all routines to check, &
     16h the statistics.)
 1310 format('1run stat on the davis-harrison pikes peak data.')
 1320 format('1run statw on the davis-harrison pikes peak data.')
 1340 format('1run stats on the davis-harrison pikes peak data.')
 1350 format('1run statws on the davis-harrison pikes peak data.')
 1370 format('1run statw on the davis-harrison pikes peak data.', &
    '  weights all equal to one.  compare to stat above, not to', &
    ' statw.')
 1380 format('series with nonzero values weighted zero.')
 1390 format(/' data = ', 10f7.4)
 1400 format ('call to statw')
 1410 format (' call to statws')
 1440 format ('1')
end
subroutine xstpld ( ldstak )

!*****************************************************************************80
!
!! XSTPLD tests the derivative checking routines.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     real delta
!        *
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real exmpt
!        the proportion of observations for which the computed
!        numerical derivatives wrt a given parameter are exempted
!        from meeting the derivative acceptance criteria.
!     integer i
!        an index variable.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ixm
!        the first dimension of the independent variable array xm.
!     integer ldsmin
!        the minimum length of the array dstak allowed.
!     integer ldstak
!        the length of the array dstak.
!     integer m
!        the number of independent variables.
!     external mdl4
!        the name of the user supplied subroutine which computes the
!        predicted values based on the current parameter estimates.
!     integer n
!        the number of observations of data.
!     integer neta
!        the number of reliable digits in the model.
!     integer npar
!        the number of unknown parameters in the model.
!     integer nprt
!        the indicator variable used to specify whether or not
!        printed output is to be provided, where if the value of
!        nprt is zero, no printed output is given.
!     integer ntest
!        the number of the current test.
!     real par(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real scale(10)
!        a dummy array, indicating use of default values for
!        the typical size of the parameters.
!     real stp(10)
!        the selected step sizes for each parameter.
!     real xm(200,2)
!        the array in which one row of the independent variable array
!        is stored.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     delta,exmpt
  integer &
     i,ixm,ldsmin,m,n,neta,npar,nprt,ntest
!
!  local arrays
  real &
     par(10),scale(10),stp(10),xm(200,2)
!
!  external subroutines
  external iprint,ldscmp,lstvec,mdl4,stpls,stpls1,stpls2,stplsc
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!     set parameter values
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)
  call stpls2(npar, stp)
  call ldscmp(14, 0, 2*(n+npar), 0, 0, 0, 's', 10*n, ldsmin)

  if ( ldstak < ldsmin ) then
    write ( *, 1020) ldsmin
    return
  end if
!
!  create independent variable
!
  delta = 0.0625e0
  xm(1,1) = 0.0e0
  do i=2,n
     xm(i,1) = xm(i-1,1) + delta
  end do

  ntest = 0
!
!  check results from valid calls
!
!  simple example
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)

  ntest = ntest + 1
  write ( *,1090) ntest
  write ( *,1040)
  write ( *,1000)
  call stpls2(npar, stp)
  call stpls(xm, n, m, ixm, mdl4, par, npar, ldsmin, stp)
  write ( *,1050) ierr
  write ( *,1080)
  call lstvec(4, stp)

  ntest = ntest + 1
  write ( *,1090) ntest
  write ( *,1040)
  write ( *,1060) neta, exmpt, scale(1), nprt
  write ( *,1010)
  call stpls2(npar, stp)
  call stplsc(xm, n, m, ixm, mdl4, par, npar, ldsmin, stp, neta, &
     exmpt, scale, nprt)
  write ( *,1100) neta, exmpt, scale(1), nprt
  write ( *,1050) ierr
  write ( *,1080)
  call lstvec(4, stp)

  return

 1000 format ('test of stpls' )
 1010 format ('test of stplsc')
 1020 format ('ldstak must be greater than or equal to ', i6)
 1040 format ('simple example')
 1050 format (/'***** returned results *****', 5x, '(-1 indicates ', &
     'value not changed by called subroutine)'//'ierr is ', i3)
 1060 format (' input   -  neta = ', i5, ', exmpt = ', g15.8, &
     ', scale(1) = ', g15.8, ', nprt = ', i5)
 1080 format (//'returned values of stp')
 1090 format ('derivative step size selection subroutine test number', &
     i5)
 1100 format (//'output  -  neta = ', i5, ', exmpt = ', g15.8, &
     ', scale(1) = ', g15.8, ', nprt = ', i5//)
end
subroutine xstple ( ldstak )

!*****************************************************************************80
!
!! XSTPLE tests the derivative checking routines.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     real delta
!        *
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real exmpt
!        the proportion of observations for which the computed
!        numerical derivatives wrt a given parameter are exempted
!        from meeting the derivative acceptance criteria.
!     integer i
!        an index variable.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ixm
!        the first dimension of the independent variable array xm.
!     integer ldsmin
!        the minimum length of the array dstak allowed.
!     integer ldstak
!        the length of the array dstak.
!     integer m
!        the number of independent variables.
!     external mdl4
!        the name of the user supplied subroutine which computes the
!        predicted values based on the current parameter estimates.
!     integer n
!        the number of observations of data.
!     integer neta
!        the number of reliable digits in the model.
!     integer npar
!        the number of unknown parameters in the model.
!     integer nprt
!        the indicator variable used to specify whether or not
!        printed output is to be provided, where if the value of
!        nprt is zero, no printed output is given.
!     integer ntest
!        the number of the current test.
!     real par(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real scale(10)
!        a dummy array, indicating use of default values for
!        the typical size of the parameters.
!     real stp(10)
!        the selected step sizes for each parameter.
!     real xm(200,2)
!        the array in which one row of the independent variable array
!        is stored.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     delta,exmpt
  integer &
     i,ixm,ldsmin,m,n,neta,npar,nprt,ntest
!
!  local arrays
  real &
     par(10),scale(10),stp(10),xm(200,2)
!
!  external subroutines
  external iprint,ldscmp,mdl4,stpls,stpls1,stpls2,stplsc
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!  set parameter values
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)
  call stpls2(npar, stp)
  call ldscmp(14, 0, 2*(n+npar), 0, 0, 0, 's', 10*n, ldsmin)

  if (ldsmin.le.ldstak) go to 5

  write ( *, 1040) ldsmin
  return

    5 continue
!
!  Create independent variable.
!
  delta = 0.0625e0
  xm(1,1) = 0.0e0
  do i=2,n
     xm(i,1) = xm(i-1,1) + delta
  end do

  ntest = 0
!
!  check error handling
!
!  test 1  -  miscellaneous error checking
!
  n = -5
  m = -5
  ixm = -10
  npar = -10

  ntest = ntest + 1
  write ( *,1090) ntest
  write ( *,1020)
  write ( *,1000)
  ierr = -1
  call stpls(xm, n, m, ixm, mdl4, par, npar, ldstak, stp)
  write ( *,1050) ierr
  write ( *,1010)
  ierr = -1
  call stplsc(xm, n, m, ixm, mdl4, par, npar, ldstak, stp, neta, &
     exmpt, scale, nprt)
  write ( *,1050) ierr
!
!  test 2  -  miscellaneous error checking (continued)
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)
  scale(2) = 0.0e0

  ntest = ntest + 1
  write ( *,1090) ntest
  write ( *,1030)
  write ( *,1000)
  ierr = -1
  call stpls(xm, n, m, ixm, mdl4, par, npar, ldsmin-1, stp)
  write ( *,1050) ierr
  write ( *,1010)
  ierr = -1
  call stplsc(xm, n, m, ixm, mdl4, par, npar, ldsmin-1, stp, neta, &
     exmpt, scale, nprt)
  write ( *,1050) ierr

  return

 1000 format (15h test of stpls )
 1010 format (15h test of stplsc)
 1020 format (32h check error handling  -  test 1)
 1030 format (32h check error handling  -  test 2)
 1040 format (45h1 *** ldstak must be greater than or equal to , i6)
 1050 format (/29h ***** returned results *****, 5x, 15h (-1 indicates , &
     39hvalue not changed by called subroutine)//9h ierr is , i3)
 1090 format (54h1derivative step size selection subroutine test number, &
     i5)
end
subroutine xstplt ( ldstak )

!*****************************************************************************80
!
!! XSTPLT tests the derivative checking routines.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     real delta
!        *
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real exmpt
!        the proportion of observations for which the computed
!        numerical derivatives wrt a given parameter are exempted
!        from meeting the derivative acceptance criteria.
!     real exmtst(5)
!        various test values for exmpt.
!     integer i
!        an index variable.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ixm
!        the first dimension of the independent variable array xm.
!     integer ldsmin
!        the minimum length of the array dstak allowed.
!     integer ldstak
!        the length of the array dstak.
!     integer m
!        the number of independent variables.
!     external mdl4
!        the name of the user supplied subroutine which computes the
!        predicted values based on the current parameter estimates.
!     integer n
!        the number of observations of data.
!     integer neta
!        the number of reliable digits in the model.
!     integer nettst(6)
!        various test values for neta.
!     integer npar
!        the number of unknown parameters in the model.
!     integer nprt
!        the indicator variable used to specify whether or not
!        printed output is to be provided, where if the value of
!        nprt is zero, no printed output is given.
!     integer ntest
!        the number of the current test.
!     real par(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real scale(10)
!        a dummy array, indicating use of default values for
!        the typical size of the parameters.
!     real stp(10)
!        the selected step sizes for each parameter.
!     real xm(200,2)
!        the array in which one row of the independent variable array
!        is stored.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     delta,exmpt
  integer &
     i,ixm,ldsmin,m,n,neta,npar,nprt,ntest
!
!  local arrays
  real &
     exmtst(5),par(10),scale(10),stp(10),xm(200,2)
  integer &
     nettst(6)
!
!  external functions
  real &
     r1mach
  external r1mach
!
!  external subroutines
  external iprint,ldscmp,lstvec,mdl4,stpls1,stpls2,stplsc
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

!
!  set parameter values
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)
  call stpls2(npar, stp)
  call ldscmp(14, 0, 2*(n+npar), 0, 0, 0, 's', 10*n, ldsmin)

  if ( ldstak < ldsmin ) then
    write ( *, 1000) ldsmin
    return
  end if
!
!     create independent variable
!
  delta = 0.0625e0
  xm(1,1) = 0.0e0
  do i=2,n
     xm(i,1) = xm(i-1,1) + delta
  end do

  ntest = 0
!
!  test various values of exmpt
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)
  exmtst(1) = -1.0e0
  exmtst(2) = 0.0001e0
  exmtst(3) = 0.5e0
  exmtst(4) = 1.0e0
  exmtst(5) = 1.1e0

  do i=1,5

     ntest = ntest + 1
     write ( *,1090) ntest
     write ( *,1040)
     write ( *,1060) neta, exmtst(i), scale(1), nprt
     write ( *,1010)
     call stpls2(npar, stp)
     call stplsc(xm, n, m, ixm, mdl4, par, npar, ldstak, stp, neta, &
        exmtst(i), scale, nprt)
     write ( *,1100) neta, exmtst(i), scale(1), nprt
     write ( *,1050) ierr
     write ( *,1080)
     call lstvec(4, stp)

  end do
!
!  test various values of neta
!
  nettst(1) = -1
  nettst(2) = 0
  nettst(3) = 1
  nettst(4) = 2

  nettst(5) = -log10(r1mach(4))
  nettst(6) = nettst(5) + 1

  scale(1) = 0.0e0

  do i=1,6

     ntest = ntest + 1
     write ( *,1090) ntest
     write ( *,1040)
     write ( *,1060) nettst(i), exmpt, scale(1), nprt
     write ( *,1010)
     call stpls2(npar, stp)
     call stplsc(xm, n, m, ixm, mdl4, par, npar, ldstak, stp, &
        nettst(i), exmpt, scale, nprt)
     write ( *,1100) nettst(i), exmpt, scale(1), nprt
     write ( *,1050) ierr
     write ( *,1080)
     call lstvec(4, stp)

   end do
!
!  suppress output
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)
  nprt = 0

  ntest = ntest + 1
  write ( *,1090) ntest
  write ( *,1040)
  write ( *,1060) neta, exmpt, scale(1), nprt
  write ( *,1010)
  call stpls2(npar, stp)
  call stplsc(xm, n, m, ixm, mdl4, par, npar, ldstak, stp, neta, &
     exmpt, scale, nprt)
  write ( *,1100) neta, exmpt, scale(1), nprt
  write ( *,1050) ierr
  write ( *,1080)
  call lstvec(4, stp)
!
!  large calculation error problem
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)
  par(3) = 10.0e0**((nettst(5)-1)/2)
  scale(1) = -1.0e0
!
  ntest = ntest + 1
  write ( *,1090) ntest
  write ( *,1070)
  write ( *,1060) neta, exmpt, scale(1), nprt
  write ( *,1010)
  call stpls2(npar, stp)
  call stplsc(xm, n, m, ixm, mdl4, par, npar, ldstak, stp, neta, &
     exmpt, scale, nprt)
  write ( *,1100) neta, exmpt, scale(1), nprt
  write ( *,1050) ierr
  write ( *,1080)
  call lstvec(4, stp)
!
  exmpt = 0.11e0
  nprt = 0
!
  ntest = ntest + 1
  write ( *,1090) ntest
  write ( *,1070)
  write ( *,1060) neta, exmpt, scale(1), nprt
  write ( *,1010)
  call stpls2(npar, stp)
  call stplsc(xm, n, m, ixm, mdl4, par, npar, ldstak, stp, neta, &
     exmpt, scale, nprt)
  write ( *,1100) neta, exmpt, scale(1), nprt
  write ( *,1050) ierr
  write ( *,1080)
  call lstvec(4, stp)
!
  return

 1000 format ('ldstak must be greater than or equal to ', i6)
 1010 format ('test of stplsc')
 1040 format ('simple example')
 1050 format (/'returned results (-1 indicates ', &
     'value not changed by called subroutine)'//' ierr is ', i3)
 1060 format (19h input   -  neta = , i5, 10h, exmpt = , g15.8, &
     13h, scale(1) = , g15.8, 9h, nprt = , i5)
 1070 format (32h large calculation error problem)
 1080 format (//23h returned values of stp)
 1090 format (54h1derivative step size selection subroutine test number, &
     i5)
 1100 format (//19h output  -  neta = , i5, 10h, exmpt = , g15.8, &
     13h, scale(1) = , g15.8, 9h, nprt = , i5//)
end
subroutine xuas ( ldstak )

!*****************************************************************************80
!
!! XUAS tests the autoregressive spectrum analysis routines.
!
!  Discussion:
!
!    Series y is the first 50 values of the series listed on page
!    318 of jenkins and watts.  the spectrum of this series is shown
!    for various bandwidth on page 270 of jenkins and watts.
!
!  Modified:
!
!    02 December 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    Gwilym Jenkins, Donald Watts,
!    Spectral Analysis and its Applications,
!    Holden-Day 1968.
!
!  Parameters:
!
!     real acov(101)
!        the autocovariance vector.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fmax, fmin
!        the maximum and minimum frequencies at which the
!        spectrum is to be computed.
!     real freq(300)
!        the vector of frequencies at which the spectrum is computed.
!     integer i
!        an index variable
!     integer iar
!        the order of the autoregressive model to be used.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list
!        if ierr == 0, no errors were detected
!        if ierr == 1, errors have been detected
!     integer lacov
!        the length of the vector acov.
!     integer lag, lagmax
!        the indexing variable indicating the lag value of the
!        autocovariance being computed and the maximum lag to be used,
!        respectively.
!     integer lds
!        the length of the vector dstak in common cstak.
!     integer lyfft
!        the length of the vector yfft.
!     integer nf
!        the number of frequencies at which the spectrum is
!        to be computed.
!     integer nprt
!        a code used to specify the type of plot, where if
!        nprt < 0 the plot is decibles/linear
!        nprt = 0 the plot is suppressed
!        nprt > 0 the plot is log/linear
!     integer ny
!        the number of observations in the series y.
!     real phi(101)
!        the vector of the order iar autoregressive model coefficients.
!     real spca(101)
!        the arrays in which the autoregressive spectrum is stored
!        for each lag window.
!     real spcf(101)
!        the arrays in which the fourier spectrum is stored
!        for each lag window.
!     real y(150)
!         the array containing the time series from jenkins and watts.
!     real yfft(400)
!        the vector of the observed time series to be analyzed using
!        the fft.
!     real ymiss
!        the user supplied code which is used to determine whether or
!        not an observation in the series is missing.  if y(i) = ymiss,
!        the value is assumed missing, otherwise it is not.
!
  implicit none

  real acov(101)
  double precision dstak(12)
  real fmax
  real fmin
  real freq(300)
  integer i
  integer iar
  integer ierr
  integer lacov
  integer lag
  integer lagmax
  integer lds
  integer ldstak
  integer lyfft
  integer nf
  integer nprt
  integer ny
  real phi(101)
  real spca(101)
  real spcf(101)
  real y(150)
  real yfft(400)
  real ymiss

  external acfs,iprint,scopy,setrv,uas,uasf,uasfs,uass,uasv,uasvs

  common /cstak/dstak
  common /errchk/ierr

  data   y(  1), y(  2), y(  3), y(  4), y(  5), y(  6) &
      /-0.88e0, -0.12e0, -0.89e0, -1.38e0, -0.07e0,  1.03e0/
  data   y(  7), y(  8), y(  9), y( 10), y( 11), y( 12) &
      / 2.14e0,  0.35e0, -1.10e0, -1.78e0, -2.76e0, -1.77e0/
  data   y( 13), y( 14), y( 15), y( 16), y( 17), y( 18) &
      / 0.98e0,  1.00e0, -0.70e0, -1.01e0, -1.30e0, -0.85e0/
  data   y( 19), y( 20), y( 21), y( 22), y( 23), y( 24) &
      /-0.46e0,  1.63e0,  0.06e0, -0.17e0, -1.01e0, -1.04e0/
  data   y( 25), y( 26), y( 27), y( 28), y( 29), y( 30) &
      /-0.66e0, -1.12e0, -0.51e0, -0.71e0, -0.20e0, -0.13e0/
  data   y( 31), y( 32), y( 33), y( 34), y( 35), y( 36) &
      / 0.14e0,  1.59e0, -0.76e0, -1.08e0, -1.77e0, -1.20e0/
  data   y( 37), y( 38), y( 39), y( 40), y( 41), y( 42) &
      / 0.45e0, -0.07e0, -0.63e0, -0.35e0, -0.87e0, -0.62e0/
  data   y( 43), y( 44), y( 45), y( 46), y( 47), y( 48) &
      / 0.28e0,  1.90e0,  2.14e0,  1.05e0,  0.31e0,  1.07e0/
  data   y( 49), y( 50) &
      / 2.67e0,  2.44e0/

!
!  Check error handling
!
!  Test 1  -  miscellaneous error checking
!
  write ( *, 2000)
  ymiss = 1.16e0
  lagmax = -1
  ny = -10
  lacov = 101
  lag = -2
  iar = -2
  lyfft = -11
  nf = -5
  fmin = 0.5e0
  fmax = 0.0e0
  nprt = -1
  lds = 0
  write ( *, 1001)
  call uas(y, ny)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 1 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt,spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1019)
  call uasf (yfft, ny, lyfft, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1020)
  call uasfs(yfft, ny, lyfft, lds, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1007)
  call uasv(acov, lagmax, ny)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1008)
  call uasvs(acov, lagmax, y, ny, iar, phi, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Test 2  -  miscellaneous error checking (continued)
!
  write ( *, 2010)
  ymiss = 1.16e0
  ny = 50
  lagmax = 50
  lag = 101
  iar = 101
  call setrv(phi, iar, 2.0e0)
  call setrv(acov, lagmax+1, 2.0e0)
  acov(1) = 1.0e0
  lyfft = -11
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 5
  lds = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 2 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1019)
  call uasf (yfft, ny, lyfft, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1020)
  call uasfs(yfft, ny, lyfft, lds, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1007)
  call uasv(acov, lagmax, ny)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1008)
  call uasvs(acov, lagmax, y, ny, iar, phi, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Test 3  -  lds too small
!
  write ( *, 2030)
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  lag = 16
  iar = 2
  phi(1) = 1.0e0
  phi(2) = -0.5e0
  call acfs (y, ny, lagmax, lacov, acov, iar, phi, 0, 700)
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 3
  lds = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 3 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1019)
  call uasf (yfft, ny, lyfft, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1020)
  call uasfs(yfft, ny, lyfft, lds, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1008)
  call uasvs(acov, lagmax, y, ny, iar, phi, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Check results from valid call
!
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  lag = 16
  iar = 2
  phi(1) = 1.0e0
  phi(2) = -0.5e0
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 2
  lds = 700
!
!  Test of UAS.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VALID PROBLEM:'
  write ( *, '(a)' ) ' '

  write ( *, 1001)
  call uas ( y, ny )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Test of UASS.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VALID PROBLEM:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 4 of UASS:'
  write ( *, '(a)' ) ' '

  call uass ( y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from UASS.
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))
!
!  Test of uasf
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VALID PROBLEM:'
  write ( *, '(a)' ) ' '

  write ( *, 1019)
  call scopy(ny, y, 1, yfft, 1)
  call uasf (yfft, ny, lyfft, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Test of uasfs
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VALID PROBLEM:'
  write ( *, '(a)' ) ' '

  write ( *, 1020)
  call scopy(ny, y, 1, yfft, 1)
  call uasfs(yfft, ny, lyfft, lds, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from uasfs
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))
!
!  Test of uasv
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VALID PROBLEM:'
  write ( *, '(a)' ) ' '

  write ( *, 1007)
  call uasv(acov, lagmax, ny)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Test of uasvs
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VALID PROBLEM:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 5 of UASS:'
  write ( *, '(a)' ) ' '

  call uasvs(acov, lagmax, y, ny, iar, phi, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from uasvs
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))
!
!  Minimum problem size
!
  ymiss = 1.16e0
  ny = 17
  lagmax = 1
  lyfft = 400
  lag = 1
  iar = -1
  phi(1) = 1.0e0
  phi(2) = -0.5e0
  nf = 1
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 2
  lds = ldstak
!
!  Test of uas
!
  write ( *, 2060)
  write ( *, 1001)
  call uas(y, ny)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Test of UASS.
!
  write ( *, 2060)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 6 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from UASS.
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))
!
!  Check handling of fmin and fmax, and lag==0 and iar==0
!
  ny = 50
  lagmax = 49
  lyfft = 400
  lag = 0
  iar = 0
  phi(1) = 1.0e0
  phi(2) = -0.5e0
  nf = 51
  fmin = 0.45e0
  fmax = 0.5e0
  nprt = 2
!
!  Test of UASS.
!
  write ( *, 2070)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 7 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from UASS.
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))
!
!  White noise spectrum
!
  ymiss = 1.16e0
  call setrv(yfft, ny, 0.0e0)
  ny = 50
  lagmax = 49
  lyfft = 400
  lag = 16
  iar = 2
  phi(1) = 1.0e0
  phi(2) = -0.5e0
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 2
!
!  Test of uass
!
  write ( *, 2080)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 8 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(yfft, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from UASS.
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))
!
!  Suppress output and check handling of lag .lt.0 and iar .lt. 0
!
  ny = 50
  lagmax = 49
  lyfft = 400
  lag = 0
  iar = 0
  phi(1) = 1.0e0
  phi(2) = -0.5e0
  nf = 51
  fmin = 0.45e0
  fmax = 0.5e0
  nprt = 0
!
!  Test of UASS.
!
  write ( *, 2090)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 9 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from UASS.
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))

  return

 1001 format ('test of uas')
 1004 format (3(1x, e16.8))
 1005 format (/' iar = ', i5/' lag = ', i5)
 1006 format (/' phi = ', (1x, 5e21.8))
 1007 format ('test of uasv')
 1008 format ('test of uasvs')
 1019 format ('test of uasf')
 1020 format ('test of uasfs')
 2000 format ('check error handling  -  test 1')
 2010 format ('check error handling  -  test 2')
 2030 format ('lds too small')
 2060 format ('minimum problem size')
 2070 format ( &
     'check handling of fmin and fmax, lag and iar equal to zero')
 2080 format ('white noise spectrum')
 2090 format ('suppress output, lag and iar less than zero')
end
subroutine xufs ( ldstak )

!*****************************************************************************80
!
!! XUFS tests the Fourier spectrum analysis routines.
!
!  Discussion:
!
!    series y is the first 50 values of the series listed on page
!    318 of jenkins and watts.  the spectrum of this series is shown
!    for various bandwidth on page 270 of jenkins and watts.
!
!    series z is the wolf sunspot numbers from 1700 to 1960 as
!    tabulated by waldmeier.  the raw and smoothed periodograms of
!    tapered series are shown on pages 95 and 176, respectively, of
!    bloomfield.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    Gwilym Jenkins, Donald Watts,
!    Spectral Analysis and its Applications,
!    Holden-Day 1968.
!
!    Max Waldmeier,
!    The Sunspot-Activity in the Years 1610-1960,
!    Shulthess, Zurich, 1961.
!
!  Parameters:
!
!     real acov(101)
!        the autocovariance vector.
!     real amiss
!         the missing value code for the returned acvf estimates.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fmax, fmin
!        the maximum and minimum frequencies at which the
!        spectrum is to be computed.
!     real freq(300)
!        the vector of frequencies at which the spectrum is computed.
!     integer i
!        an index variable
!     integer iar
!        the order of the autoregressive model to be used.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list
!        if ierr == 0, no errors were detected
!        if ierr == 1, errors have been detected
!     integer ispcf
!         the actual dimension for the spectrum arrays.
!     integer j
!        index variable.
!     integer lacov
!        the length of the vector acov.
!     integer lagmax
!        the indexing variable indicating the lag value of the
!        autocovariance being computed and the maximum lag to be used,
!        respectively.
!     integer lags(4)
!        the array used to store the lag window truccation
!        points used for each set of spectrum values.
!     integer lds
!        the length of the vector dstak in common cstak.
!     integer lyfft
!        the length of the vector yfft.
!     integer nf
!        the number of frequencies at which the spectrum is
!        to be computed.
!     integer nlppa(101)
!        the numbers of lagged product pairs used for each acvf.
!     integer nprt
!        a code used to specify the type of plot, where if
!        nprt = 0 the plot is suppressed, if
!        nprt = 2 the plot is decibels/linear, if
!        nprt = 2 the plot is log/linear, if
!        nprt = 3 the plot is decibels/log, and if
!        nprt = 4 the plot is log/log.
!     integer nw
!        the number of different lag window truncation points specified,
!        and therefore, the number of plots.
!     integer ny
!        the number of observations in the series y.
!     real phi(100)
!        the vector of the order iar autoregressive model coefficients.
!     real spcf(101, 4)
!        the arrays in which the fourier spectrum is stored
!        for each lag window.
!     real y(150)
!         the array containing the time series from jenkins and watts.
!     real yfft(400)
!        the vector of the observed time series to be analyzed using
!        the fft.
!     real ymiss
!        the user supplied code which is used to determine whether or
!        not an observation in the series is missing.  if y(i) = ymiss,
!        the value is assumed missing, otherwise it is not.
!
  implicit none

  integer &
     ldstak
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     amiss,fmax,fmin,ymiss
  integer &
     i,iar,ispcf,j,lacov,lagmax,lds,lyfft,nf,nprt,nw,ny
!
!  local arrays
  real &
     acov(101),freq(300),phi(100),spcf(101,4),y(150),yfft(400)
  integer &
     lags(4),nlppa(101)
!
!  external subroutines
  external acfms,acfs,iprint,nrand,scopy,setrv,ufs,ufsf,ufsfs,ufsm, &
     ufsms,ufsmv,ufsmvs,ufss,ufsv,ufsvs
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data   y(  1), y(  2), y(  3), y(  4), y(  5), y(  6) &
      /-0.88e0, -0.12e0, -0.89e0, -1.38e0, -0.07e0,  1.03e0/
  data   y(  7), y(  8), y(  9), y( 10), y( 11), y( 12) &
      / 2.14e0,  0.35e0, -1.10e0, -1.78e0, -2.76e0, -1.77e0/
  data   y( 13), y( 14), y( 15), y( 16), y( 17), y( 18) &
      / 0.98e0,  1.00e0, -0.70e0, -1.01e0, -1.30e0, -0.85e0/
  data   y( 19), y( 20), y( 21), y( 22), y( 23), y( 24) &
      /-0.46e0,  1.63e0,  0.06e0, -0.17e0, -1.01e0, -1.04e0/
  data   y( 25), y( 26), y( 27), y( 28), y( 29), y( 30) &
      /-0.66e0, -1.12e0, -0.51e0, -0.71e0, -0.20e0, -0.13e0/
  data   y( 31), y( 32), y( 33), y( 34), y( 35), y( 36) &
      / 0.14e0,  1.59e0, -0.76e0, -1.08e0, -1.77e0, -1.20e0/
  data   y( 37), y( 38), y( 39), y( 40), y( 41), y( 42) &
      / 0.45e0, -0.07e0, -0.63e0, -0.35e0, -0.87e0, -0.62e0/
  data   y( 43), y( 44), y( 45), y( 46), y( 47), y( 48) &
      / 0.28e0,  1.90e0,  2.14e0,  1.05e0,  0.31e0,  1.07e0/
  data   y( 49), y( 50) &
      / 2.67e0,  2.44e0/

!     check error handling
!
!        test 1  -  miscellaneous error checking
!
  write ( *, 2000)
  ymiss = 1.16e0
  lagmax = -1
  ny = -10
  lacov = 101
  lyfft = -11
  nw = -1
  nf = -5
  fmin = 0.5e0
  fmax = 0.0e0
  nprt = -1
  ispcf = -20
  lds = 0
  write ( *, 1001)
  call ufs (y, ny)
  write ( *, 1002) ierr
  write ( *, 1003)
  call ufss(y, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1019)
  call ufsf (yfft, ny, lyfft, lds)
  write ( *, 1002) ierr
  write ( *, 1020)
  call ufsfs(yfft, ny, lyfft, lds, nw, lags, nf, fmin, fmax, nprt, &
     spcf, ispcf, freq)
  write ( *, 1002) ierr
  write ( *, 1005)
  call ufsm (y, ymiss, ny)
  write ( *, 1002) ierr
  write ( *, 1006)
  call ufsms(y, ymiss, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1007)
  call ufsv(acov, lagmax, ny)
  write ( *, 1002) ierr
  write ( *, 1008)
  call ufsvs (acov, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1021)
  call ufsmv(acov, nlppa, lagmax, ny)
  write ( *, 1002) ierr
  write ( *, 1022)
  call ufsmvs (acov, nlppa, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  test 2  -  miscellaneous error checking (continued)
!
  write ( *, 2010)
  ymiss = 1.16e0
  ny = 50
  lagmax = 55
  lyfft = -11
  nw = 2
  lags(1) = 0
  lags(2) = 50
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 3
  ispcf = 20
  lds = 0
  write ( *, 1003)
  call ufss(y, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1019)
  call ufsf (yfft, ny, lyfft, lds)
  write ( *, 1002) ierr
  write ( *, 1020)
  call ufsfs(yfft, ny, lyfft, lds, nw, lags, nf, fmin, fmax, nprt, &
     spcf, ispcf, freq)
  write ( *, 1002) ierr
  write ( *, 1006)
  call ufsms(y, ymiss, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1008)
  call ufsvs (acov, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1022)
  call ufsmvs (acov, nlppa, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  test 3  -  lds too small
!
  write ( *, 2030)
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  nw = 2
  lags(1) = 0
  lags(2) = 50
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 3
  ispcf = 101
  lds = 0
  write ( *, 1003)
  call ufss(y, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1019)
  call ufsf (yfft, ny, lyfft, lds)
  write ( *, 1002) ierr
  write ( *, 1020)
  call ufsfs(yfft, ny, lyfft, lds, nw, lags, nf, fmin, fmax, nprt, &
     spcf, ispcf, freq)
  write ( *, 1002) ierr
  write ( *, 1006)
  call ufsms(y, ymiss, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1008)
  call ufsvs (acov, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1022)
  call ufsmvs (acov, nlppa, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  test 4  -  all data and covariances missing
!
  write ( *, 2040)
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 3
  ispcf = 101
  lds = 700
  call setrv(yfft, ny, ymiss)
  call setrv(acov, lagmax, 0.0e0)
  nlppa(1:lagmax) = 0
  write ( *, 1005)
  call ufsm(yfft, ymiss, ny)
  write ( *, 1002) ierr
  write ( *, 1006)
  call ufsms(yfft, ymiss, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1021)
  call ufsmv (acov, nlppa, lagmax, ny)
  write ( *, 1002) ierr
  write ( *, 1022)
  call ufsmvs (acov, nlppa, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  test 5  -  every other value missing
!
  write ( *, 2050)
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 3
  ispcf = 101
  lds = 700
  call setrv(yfft, ny, ymiss)
  yfft(1:ny:2) = y(1:ny:2)
  write ( *, 1005)
  call ufsm(yfft, ymiss, ny)
  write ( *, 1002) ierr
  write ( *, 1006)
  call ufsms(yfft, ymiss, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  check results from valid call
!
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 2
  ispcf = 101
  lds = ldstak
!
!  test of ufs
!
  write ( *, 2020)
  write ( *, 1001)
  call ufs (y, ny)
  write ( *, 1002) ierr
!
!  test of ufss
!
  write ( *, 2020)
  write ( *, 1003)
  call ufss(y, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  print returned values from ufss
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!  test of ufsf
!
  write ( *, 2020)
  write ( *, 1019)
  call scopy(ny, y, 1, yfft, 1)
  call ufsf (yfft, ny, lyfft, lds)
  write ( *, 1002) ierr
!
!  test of ufsfs
!
  write ( *, 2020)
  write ( *, 1020)
  call scopy(ny, y, 1, yfft, 1)
  call ufsfs(yfft, ny, lyfft, lds, nw, lags, nf, fmin, fmax, nprt, &
     spcf, ispcf, freq)
  write ( *, 1002) ierr
!
!  print returned values from ufsfs
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!  test of ufsm
!
  write ( *, 2020)
  write ( *, 1005)
  call ufsm (y, ymiss, ny)
  write ( *, 1002) ierr
!
!     test of ufsms
!
  write ( *, 2020)
  write ( *, 1006)
  call ufsms(y, ymiss, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
!
!     print returned values from ufsms
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!     test of ufsv
!
  write ( *, 2020)
  call acfs (y, ny, lagmax, lacov, acov, iar, phi, 0, lds)
  write ( *, 1007)
  call ufsv(acov, lagmax, ny)
  write ( *, 1002) ierr
!
!     test of ufsvs
!
  write ( *, 2020)
  write ( *, 1008)
  call ufsvs (acov, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
!
!     print returned values from ufsvs
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!     test of ufsmv
!
  write ( *, 2020)
  call acfms (y, ymiss, ny, lagmax, lacov, acov, amiss, nlppa, &
     0, lds)
  write ( *, 1021)
  call ufsmv(acov, nlppa, lagmax, ny)
  write ( *, 1002) ierr
!
!     test of ufsmvs
!
  write ( *, 2020)
  write ( *, 1022)
  call ufsmvs (acov, nlppa, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  print returned values from ufsmvs
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!  minimum problem size
!
  ymiss = 1.16e0
  ny = 17
  lagmax = 1
  lyfft = 400
  nw = 2
  lags(1) = 1
  lags(2) = 16
  nf = 1
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 2
  ispcf = 101
  lds = 700
!
!  test of ufs
!
  write ( *, 2060)
  write ( *, 1001)
  call ufs (y, ny)
  write ( *, 1002) ierr
!
!  test of ufss
!
  write ( *, 2060)
  write ( *, 1003)
  call ufss(y, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  print returned values from ufss
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!  check handling of fmin and fmax
!
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 51
  fmin = 0.45e0
  fmax = 0.5e0
  nprt = 2
  ispcf = 101
  lds = 700
!
!  test of ufss
!
  write ( *, 2070)
  write ( *, 1003)
  call ufss(y, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
!
!     print returned values from ufss
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!  white noise spectrum
!
  ymiss = 1.16e0
  call nrand(yfft, ny, 12345)
  ny = 50
  lagmax = 49
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 2
  ispcf = 101
  lds = 700
!
!  test of ufss
!
  write ( *, 2080)
  write ( *, 1003)
  call ufss(yfft, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  print returned values from ufss
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)

  return

 1001 format (' test of ufs')
 1002 format (/' ierr is', i5/)
 1003 format (' test of ufss')
 1004 format (3(1x, e16.8))
 1005 format (' test of ufsm')
 1006 format (' test of ufsms')
 1007 format (' test of ufsv')
 1008 format (' test of ufsvs')
 1019 format (' test of ufsf')
 1020 format (' test of ufsfs')
 1021 format (' test of ufsmv')
 1022 format (' test of ufsmvs')
 2000 format ('1check error handling  -  test 1')
 2010 format ('1check error handling  -  test 2')
 2020 format ('1valid problem')
 2030 format ('1lds too small')
 2040 format ('1all data and covariances missing')
 2050 format ('1every other data value missing')
 2060 format ('1minimum problem size')
 2070 format ('1check handling of fmin and fmax')
 2080 format ('1white noise spectrum')
end
subroutine xvp ( )

!*****************************************************************************80
!
!! XVP tests the plotting subroutines.
!
!  Discussion:
!
!    series y is the airline data listed on page 531 of box and jenkins.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    George Box, Gwilym Jenkins,
!    Time Series Analysis: Forecasting and Control,
!    Holden-Day, 1970,
!    QA280.B67.
!
!  Parameters:
!
!     real air(144)
!        the airline data.
!     integer ibar
!        the variable used to determine if single points (ibar .ne. 0)
!        or bars (ibar == 0) are to be plotted.
!     integer ierr
!        a common variable used as a flag to indicate whether
!        or not there are any errors, if =0 then no errors.
!     integer ilog
!        the two digit integer, pq, used to select axis scale, where
!        p designates the x-axis and q designates the y-axis.
!        if p==0 (q==0), then the x-axis (y-axis) is linear.
!        if p.ne.0 (q.ne.0), then the x-axis (y-axis) is log.
!     integer irlin
!        the indicator variable used to designate whether zero or the
!        series mean is to be plotted as a reference line, or whether
!        no reference line is to be plotted.
!        if irlin <= -1, no reference line is plotted.
!        if irlin ==  0, zero is plotted as the reference line.
!        if irlin .ge.  1, the series mean is plotted.
!     integer isize
!        the two digit integer, pq, used to select axis size, where
!        p designates the x-axis and q designates the y-axis.
!        if p==0 (q==0), then the x-axis (y-axis) is the maximum.
!        if p.ne.0 (q.ne.0), then the x-axis (y-axis) is half the maximu
!     integer isym(144)
!        vector containing symbol designations for plotting
!     integer itest
!        the number of the test.
!     integer iym
!        actual dimension of ym in users main program
!     integer ldstak
!        *
!     integer m
!        the number of vectors in ym
!     integer ns
!        the sampling frequency,
!        where if ns <= 1, every point is plotted,
!                       = 2, every other point is plotted,
!                       = 3, every third point is plotted, etc.
!     integer ny, nym
!        the number of observations in arrays y and ym, respectively.
!     real xinc
!        the increment for the x axis.
!     real xlb
!        the lower bound for the x-axis.
!     real y(144)
!        vector of observations for the y (vertical) coordinates
!     real ylb
!        the lower bound for the y-axis.  (ylb=yub indicates limits are
!        to be determined from the range of the data.)
!     real ym(12,12)
!        multivariate observations for the y (vertical) coordinates.
!     real ymiss
!        the missing value code for the y-axis.
!     real ymmiss(144)
!        the missing value codes for each column of ym.
!     real yub
!        the upper bound for the y-axis.  (ylb=yub indicates limits are
!        to be determined from the range of the data.)
!
!
  implicit none
!
!  scalars in common
  integer &
     ierr
!
!  local scalars
  real &
     xinc,xlb,ylb,ymiss,yub
  integer &
     ibar,ilog,irlin,isize,itest,iym,m,nout,ns,ny,nym
!
!  local arrays
  real &
     air(144),y(144),ym(12,12),ymmiss(144)
  integer &
     isym(144)
!
!  external subroutines
  external iprint,mvp,mvpc,mvpl,mvpm,mvpmc,mvpml,scopy,setrv,svp, &
     svpc,svpl,svpm,svpmc,svpml,vp,vpc,vpl,vpm,vpmc,vpml
!
!  common blocks
  common /errchk/ierr
!
!  equivalences
  equivalence (y(1),ym(1,1))

  data ymiss/180.0e0/
!
  data isym(  1),isym(  2),isym(  3),isym(  4),isym(  5),isym(  6) &
      /    -5000,     6000,        7,        8,        9,       10/
  data isym(  7),isym(  8),isym(  9),isym( 10),isym( 11),isym( 12) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 13),isym( 14),isym( 15),isym( 16),isym( 17),isym( 18) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 19),isym( 20),isym( 21),isym( 22),isym( 23),isym( 24) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 25),isym( 26),isym( 27),isym( 28),isym( 29),isym( 30) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 31),isym( 32),isym( 33),isym( 34),isym( 35),isym( 36) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 37),isym( 38),isym( 39),isym( 40),isym( 41),isym( 42) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 43),isym( 44),isym( 45),isym( 46),isym( 47),isym( 48) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 49),isym( 50),isym( 51),isym( 52),isym( 53),isym( 54) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 55),isym( 56),isym( 57),isym( 58),isym( 59),isym( 60) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 61),isym( 62),isym( 63),isym( 64),isym( 65),isym( 66) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 67),isym( 68),isym( 69),isym( 70),isym( 71),isym( 72) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 73),isym( 74),isym( 75),isym( 76),isym( 77),isym( 78) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 79),isym( 80),isym( 81),isym( 82),isym( 83),isym( 84) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 85),isym( 86),isym( 87),isym( 88),isym( 89),isym( 90) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 91),isym( 92),isym( 93),isym( 94),isym( 95),isym( 96) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 97),isym( 98),isym( 99),isym(100),isym(101),isym(102) &
      /        5,        6,        7,        8,        9,       10/
  data isym(103),isym(104),isym(105),isym(106),isym(107),isym(108) &
      /       11,       12,       13,       14,       15,       16/
  data isym(109),isym(110),isym(111),isym(112),isym(113),isym(114) &
      /        5,        6,        7,        8,        9,       10/
  data isym(115),isym(116),isym(117),isym(118),isym(119),isym(120) &
      /       11,       12,       13,       14,       15,       16/
  data isym(121),isym(122),isym(123),isym(124),isym(125),isym(126) &
      /        5,        6,        7,        8,        9,       10/
  data isym(127),isym(128),isym(129),isym(130),isym(131),isym(132) &
      /       11,       12,       13,       14,       15,       16/
  data isym(133),isym(134),isym(135),isym(136),isym(137),isym(138) &
      /        5,        6,        7,        8,        9,       10/
  data isym(139),isym(140),isym(141),isym(142),isym(143),isym(144) &
      /       11,       12,       13,       14,       15,       16/
!
  data  air(  1), air(  2), air(  3), air(  4), air(  5), air(  6) &
      / 112.0e0, 118.0e0, 132.0e0, 129.0e0, 121.0e0, 135.0e0/
  data  air(  7), air(  8), air(  9), air( 10), air( 11), air( 12) &
      / 148.0e0, 148.0e0, 136.0e0, 119.0e0, 104.0e0, 118.0e0/
  data  air( 13), air( 14), air( 15), air( 16), air( 17), air( 18) &
      / 115.0e0, 126.0e0, 141.0e0, 135.0e0, 125.0e0, 149.0e0/
  data  air( 19), air( 20), air( 21), air( 22), air( 23), air( 24) &
      / 170.0e0, 170.0e0, 158.0e0, 133.0e0, 114.0e0, 140.0e0/
  data  air( 25), air( 26), air( 27), air( 28), air( 29), air( 30) &
      / 145.0e0, 150.0e0, 178.0e0, 163.0e0, 172.0e0, 178.0e0/
  data  air( 31), air( 32), air( 33), air( 34), air( 35), air( 36) &
      / 199.0e0, 199.0e0, 184.0e0, 162.0e0, 146.0e0, 166.0e0/
  data  air( 37), air( 38), air( 39), air( 40), air( 41), air( 42) &
      / 171.0e0, 180.0e0, 193.0e0, 181.0e0, 183.0e0, 218.0e0/
  data  air( 43), air( 44), air( 45), air( 46), air( 47), air( 48) &
      / 230.0e0, 242.0e0, 209.0e0, 191.0e0, 172.0e0, 194.0e0/
  data  air( 49), air( 50), air( 51), air( 52), air( 53), air( 54) &
      / 196.0e0, 196.0e0, 236.0e0, 235.0e0, 229.0e0, 243.0e0/
  data  air( 55), air( 56), air( 57), air( 58), air( 59), air( 60) &
      / 264.0e0, 272.0e0, 237.0e0, 211.0e0, 180.0e0, 201.0e0/
  data  air( 61), air( 62), air( 63), air( 64), air( 65), air( 66) &
      / 204.0e0, 188.0e0, 235.0e0, 227.0e0, 234.0e0, 264.0e0/
  data  air( 67), air( 68), air( 69), air( 70), air( 71), air( 72) &
      / 302.0e0, 293.0e0, 259.0e0, 229.0e0, 203.0e0, 229.0e0/
  data  air( 73), air( 74), air( 75), air( 76), air( 77), air( 78) &
      / 242.0e0, 233.0e0, 267.0e0, 269.0e0, 270.0e0, 315.0e0/
  data  air( 79), air( 80), air( 81), air( 82), air( 83), air( 84) &
      / 364.0e0, 347.0e0, 312.0e0, 274.0e0, 237.0e0, 278.0e0/
  data  air( 85), air( 86), air( 87), air( 88), air( 89), air( 90) &
      / 284.0e0, 277.0e0, 317.0e0, 313.0e0, 318.0e0, 374.0e0/
  data  air( 91), air( 92), air( 93), air( 94), air( 95), air( 96) &
      / 413.0e0, 405.0e0, 355.0e0, 306.0e0, 271.0e0, 306.0e0/
  data  air( 97), air( 98), air( 99), air(100), air(101), air(102) &
      / 315.0e0, 301.0e0, 356.0e0, 348.0e0, 355.0e0, 422.0e0/
  data  air(103), air(104), air(105), air(106), air(107), air(108) &
      / 465.0e0, 467.0e0, 404.0e0, 347.0e0, 305.0e0, 336.0e0/
  data  air(109), air(110), air(111), air(112), air(113), air(114) &
      / 340.0e0, 318.0e0, 362.0e0, 348.0e0, 363.0e0, 435.0e0/
  data  air(115), air(116), air(117), air(118), air(119), air(120) &
      / 491.0e0, 505.0e0, 404.0e0, 359.0e0, 310.0e0, 337.0e0/
  data  air(121), air(122), air(123), air(124), air(125), air(126) &
      / 360.0e0, 342.0e0, 406.0e0, 396.0e0, 420.0e0, 472.0e0/
  data  air(127), air(128), air(129), air(130), air(131), air(132) &
      / 548.0e0, 559.0e0, 463.0e0, 407.0e0, 362.0e0, 405.0e0/
  data  air(133), air(134), air(135), air(136), air(137), air(138) &
      / 417.0e0, 391.0e0, 419.0e0, 461.0e0, 472.0e0, 535.0e0/
  data  air(139), air(140), air(141), air(142), air(143), air(144) &
      / 622.0e0, 606.0e0, 508.0e0, 461.0e0, 390.0e0, 432.0e0/
!
  call setrv(ymmiss, 144, ymiss)
  call scopy(144, air, 1, y, 1)

  itest = 0
!
!  short calls
!
  ny = 144
  nym = 12
  iym = 12
  m = 12
  ns = 1
  ilog = -1
  isize = -1
  isize = -1
  irlin = -1
  ibar = -1
  ylb = 0.0e0
  yub = 0.0e0
  xlb = 0.0e0
  xinc = 0.0e0

   10 continue
!
!  test of vp
!
  write ( *, 2000)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  call vp ( y, ny, ns )
  write ( *, 3000) ierr
!
!  test of vpm
!
  write ( *, 2030)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  call vpm ( y, ymiss, ny, ns )
  write ( *, 3000) ierr
!
!  test of svp
!
  write ( *, 2120)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  call svp (y, ny, ns, isym)
  write ( *, 3000) ierr
!
!  test of svpm
!
  write ( *, 2150)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  call svpm (y, ymiss, ny, ns, isym)
  write ( *, 3000) ierr
!
!  test of mvp
!
  write ( *, 2060)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3030) ns
  call mvp (ym, nym, m, iym, ns)
  write ( *, 3000) ierr
!
!  test of mvpm
!
  write ( *, 2090)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3030) ns
  call mvpm (ym, ymmiss, nym, m, iym, ns)
  write ( *, 3000) ierr
!
!     log option calls
!
   20 continue
!
!     test of vpl
!
  write ( *, 2010)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  call vpl (y, ny, ns, ilog)
  write ( *, 3000) ierr
!
!     test of vpml
!
  write ( *, 2040)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  call vpml (y, ymiss, ny, ns, ilog)
  write ( *, 3000) ierr
!
!     test of svpl
!
  write ( *, 2130)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  call svpl (y, ny, ns, isym, ilog)
  write ( *, 3000) ierr
!
!     test of svpml
!
  write ( *, 2160)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  call svpml (y, ymiss, ny, ns, isym, ilog)
  write ( *, 3000) ierr
!
!     test of mvpl
!
  write ( *, 2070)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3030) ns
  write ( *, 3040) ilog
  call mvpl (ym, nym, m, iym, ns, ilog)
  write ( *, 3000) ierr
!
!     test of mvpml
!
  write ( *, 2100)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3030) ns
  write ( *, 3040) ilog
  call mvpml(ym, ymmiss, nym, m, iym, ns, ilog)
  write ( *, 3000) ierr
!
!     test of long calls
!
   30 continue
!
!     test of vpc
!
  write ( *, 2020)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  write ( *, 3060) isize, irlin, ibar
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3090) xinc
  call vpc (y, ny, ns, ilog, isize, irlin, ibar, ylb, &
     yub, xlb, xinc)
  write ( *, 3000) ierr
!
!     test of vpmc
!
  write ( *, 2050)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  write ( *, 3060) isize, irlin, ibar
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3090) xinc
  call vpmc (y, ymiss, ny, ns, ilog, isize, irlin, ibar, ylb, &
     yub, xlb, xinc)
  write ( *, 3000) ierr
!
!     test of svpc
!
  write ( *, 2140)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  write ( *, 3060) isize, irlin, ibar
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3090) xinc
  call svpc (y, ny, ns, isym, ilog, isize, irlin, ibar, ylb, &
     yub, xlb, xinc)
  write ( *, 3000) ierr
!
!     test of svpmc
!
  write ( *, 2170)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  write ( *, 3060) isize, irlin, ibar
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3090) xinc
  call svpmc(y, ymiss, ny, ns, isym, ilog, isize, irlin, ibar, &
     ylb, yub, xlb, xinc)
  write ( *, 3000) ierr
!
!     test of mvpc
!
   40 write ( *, 2080)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3030) ns
  write ( *, 3040) ilog
  write ( *, 3060) isize, irlin, ibar
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3090) xinc
  call mvpc(ym, nym, m, iym, ns, ilog, isize, ylb, &
     yub, xlb, xinc)
  write ( *, 3000) ierr
!
!     test of mvpmc
!
   50 write ( *, 2110)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3030) ns
  write ( *, 3040) ilog
  write ( *, 3060) isize, irlin, ibar
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3090) xinc
  call mvpmc(ym, ymmiss, nym, m, iym, ns, ilog, isize, ylb, &
     yub, xlb, xinc)
  write ( *, 3000) ierr
!
  itest = itest + 1
!
  go to (110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300), &
    itest
!
!     test valid options
!
  110 ilog = 0
  isize = 0
  ylb = 100.0e0
  yub = 700.0e0
  xlb = 4.0e0
  xinc = 16.0e0
  go to 20
!
  120 ilog = 2
  isize = 2
  nout = 5
  xinc = -1.0e0
  go to 20
!
  130 ilog = 20
  isize = 20
  nout = 55
  yub = 300.0e0
  go to 30
!
  140 ilog = 22
  isize = 22
  go to 40
!
  150 ny = 1
  nym = 1
  m = 144
  iym = 1
  go to 40
!
  160 call setrv(y, 144, 1.0e0)
  nym = 6
  iym = 12
  m = 6
  ny = 36
  ylb = 0.0e0
  yub = 0.0e0
  xlb = 0.0e0
  xinc = 0.0e0
  go to 30
!
!     test error response
!
  170 ny = 0
  nym = 0
  m = 0
  iym = -1
  go to 10

  180 ny = 144
  nym = 12
  m = 12
  iym = -1
  xlb = -1.0e0
  ylb = -1.0e0
  go to 40

  190 iym = 12
  y(1) = 0.0e0
  go to 50

  200 call setrv(y, 144, ymiss)
  xlb = xinc
  ylb = yub
  go to 50

  300 continue

  return

 2000 format ('test of vp')
 2010 format ('test of vpl')
 2020 format ('test of vpc')
 2030 format ('test of vpm')
 2040 format ('test of vpml')
 2050 format ('test of vpmc')
 2060 format ('test of mvp')
 2070 format ('test of mvpl')
 2080 format ('test of mvpc')
 2090 format ('test of mvpm')
 2100 format ('test of mvpml')
 2110 format ('test of mvpmc')
 2120 format ('test of svp')
 2130 format ('test of svpl')
 2140 format ('test of svpc')
 2150 format ('test of svpm')
 2160 format ('test of svpml')
 2170 format ('test of svpmc')
 3000 format (/' ierr = ', i4)
 3010 format (' ', 5x, '   n     =', i5)
 3020 format ('+', 20x, 10h / m     =, i5, 10h / iym   =, i5)
 3030 format ('+', 50x, 10h / ns    =, i5)
 3040 format ('+', 65x, 10h / ilog  =, i5)
 3060 format (' ',  5x, '   isize=', i5, ' / irlin=', i5, &
     10h / ibar  =, i5)
 3070 format ('+', 50x, 10h / ylb   =, f10.4, 10h / yub   =, f10.4, &
     10h / xlb   =, f10.4)
 3090 format ('+', 110x, 10h / xinc  =, f10.4)
 3100 format (' test number ', i5)
end
subroutine xxch1 ( ldstak )

!*****************************************************************************80
!
!! XXCH1 tests the page plot and statistical analysis families of routines.
!
!  Discussion:
!
!    The data set is 84 relative humidity measurememts from pikes peak.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson and John Koontz,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index variable.
!     integer ierr
!        error flag
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer n
!        the length of the vector y.
!     real x(100)
!        the order indices of the data.
!     real y(100)
!        data vector for tests.
!
!
  implicit none

      integer &
         ldstak
!
!  scalars in common
      integer &
         ierr
!
!  arrays in common
      double precision dstak(12)
!
!  local scalars
      integer &
         i,n
!
!  local arrays
      real &
         x(100),y(100)
!
!  external subroutines
      external iprint,pp,stat
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr

      data n /84/
!
      data y( 1), y( 2), y( 3), y( 4) &
          / 0.6067e0, 0.6087e0, 0.6086e0, 0.6134e0/
      data y( 5), y( 6), y( 7) &
          / 0.6108e0, 0.6138e0, 0.6125e0/
      data y( 8), y( 9), y(10), y(11) &
          / 0.6122e0, 0.6110e0, 0.6104e0, 0.7213e0/
      data y(12), y(13), y(14) &
          / 0.7078e0, 0.7021e0, 0.7004e0/
      data y(15), y(16), y(17), y(18) &
          / 0.6981e0, 0.7242e0, 0.7268e0, 0.7418e0/
      data y(19), y(20), y(21) &
          / 0.7407e0, 0.7199e0, 0.6225e0/
      data y(22), y(23), y(24), y(25) &
          / 0.6254e0, 0.6252e0, 0.6267e0, 0.6218e0/
      data y(26), y(27), y(28) &
          / 0.6178e0, 0.6216e0, 0.6192e0/
      data y(29), y(30), y(31), y(32) &
          / 0.6191e0, 0.6250e0, 0.6188e0, 0.6233e0/
      data y(33), y(34), y(35) &
          / 0.6225e0, 0.6204e0, 0.6207e0/
      data y(36), y(37), y(38), y(39) &
          / 0.6168e0, 0.6141e0, 0.6291e0, 0.6231e0/
      data y(40), y(41), y(42) &
          / 0.6222e0, 0.6252e0, 0.6308e0/
      data y(43), y(44), y(45), y(46) &
          / 0.6376e0, 0.6330e0, 0.6303e0, 0.6301e0/
      data y(47), y(48), y(49) &
          / 0.6390e0, 0.6423e0, 0.6300e0/
      data y(50), y(51), y(52), y(53) &
          / 0.6260e0, 0.6292e0, 0.6298e0, 0.6290e0/
      data y(54), y(55), y(56) &
          / 0.6262e0, 0.5952e0, 0.5951e0/
      data y(57), y(58), y(59), y(60) &
          / 0.6314e0, 0.6440e0, 0.6439e0, 0.6326e0/
      data y(61), y(62), y(63) &
          / 0.6392e0, 0.6417e0, 0.6412e0/
      data y(64), y(65), y(66), y(67) &
          / 0.6530e0, 0.6411e0, 0.6355e0, 0.6344e0/
      data y(68), y(69), y(70) &
          / 0.6623e0, 0.6276e0, 0.6307e0/
      data y(71), y(72), y(73), y(74) &
          / 0.6354e0, 0.6197e0, 0.6153e0, 0.6340e0/
      data y(75), y(76), y(77) &
          / 0.6338e0, 0.6284e0, 0.6162e0/
      data y(78), y(79), y(80), y(81) &
          / 0.6252e0, 0.6349e0, 0.6344e0, 0.6361e0/
      data y(82), y(83), y(84) &
          / 0.6373e0, 0.6337e0, 0.6383e0/

      do i=1,n
         x(i) = i
      end do
!
!     print heading
!
      write ( *,1000)
!
!     perform simple test of pp
!
      write ( *,1100)
      call pp(y, x, n)
      write ( *,2000) ierr
!
!     perform simple test of stat
!
      write ( *,1200)
      call stat(y, n, ldstak)
      write ( *,2000) ierr

      return
!
!     formats
!
 1000 format ('1*ch1')
 1100 format (' simple test of pp')
 1200 format ('1simple test of stat')
 2000 format (/' the value of ierr is ', i4)
end
subroutine xxch2 ( )

!*****************************************************************************80
!
!! XXCH2 tests the page plot family of routines.
!
!  Discussion:
!
!    Data is the airline data listed on page 531 of box and jenkins.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    George Box, Gwilym Jenkins,
!    Time Series Analysis: Forecasting and Control,
!    Holden-Day, 1970,
!    QA280.B67.
!
!  Parameters:
!
!     real air(144)
!        the airline data.
!     integer ierr
!        a common variable used as a flag to indicate whether
!        or not there are any errors, if =0 then no errors.
!     integer iym
!        the exact value of the first dimension of the matrix ym.
!     integer ldstak
!        a dummy variable for this test subprogram.
!     integer m
!        the number of columns of data in ym.
!     integer n
!        the number of observations in each column of ym.
!     real x(12)
!        vector of observations for x(horizontal) coordinates
!     real ym(12,12)
!        multivariate observations for the y (vertical) coordinates.
!
  implicit none
!
!  scalars in common
      integer &
         ierr
!
!  local scalars
      integer iym,m,n
!
!  local arrays
      real &
         air(144),x(12),ym(12,12)
!
!  external subroutines
      external iprint,mpp
!
!  common blocks
      common /errchk/ierr
!
!  equivalences
      equivalence (air(1),ym(1,1))

!
      data    x(  1),   x(  2),   x(  3),   x(  4),   x(  5),   x(  6) &
          /   1.0e0,    2.0e0,    3.0e0,    4.0e0,    5.0e0,    6.0e0/
      data    x(  7),   x(  8),   x(  9),   x( 10),   x( 11),   x( 12) &
          /   7.0e0,    8.0e0,    9.0e0,   10.0e0,   11.0e0,   12.0e0/
!
      data  air(  1), air(  2), air(  3), air(  4), air(  5), air(  6) &
          / 112.0e0,  118.0e0,  132.0e0,  129.0e0,  121.0e0,  135.0e0/
      data  air(  7), air(  8), air(  9), air( 10), air( 11), air( 12) &
          / 148.0e0,  148.0e0,  136.0e0,  119.0e0,  104.0e0,  118.0e0/
      data  air( 13), air( 14), air( 15), air( 16), air( 17), air( 18) &
          / 115.0e0,  126.0e0,  141.0e0,  135.0e0,  125.0e0,  149.0e0/
      data  air( 19), air( 20), air( 21), air( 22), air( 23), air( 24) &
          / 170.0e0,  170.0e0,  158.0e0,  133.0e0,  114.0e0,  140.0e0/
      data  air( 25), air( 26), air( 27), air( 28), air( 29), air( 30) &
          / 145.0e0,  150.0e0,  178.0e0,  163.0e0,  172.0e0,  178.0e0/
      data  air( 31), air( 32), air( 33), air( 34), air( 35), air( 36) &
          / 199.0e0,  199.0e0,  184.0e0,  162.0e0,  146.0e0,  166.0e0/
      data  air( 37), air( 38), air( 39), air( 40), air( 41), air( 42) &
          / 171.0e0,  180.0e0,  193.0e0,  181.0e0,  183.0e0,  218.0e0/
      data  air( 43), air( 44), air( 45), air( 46), air( 47), air( 48) &
          / 230.0e0,  242.0e0,  209.0e0,  191.0e0,  172.0e0,  194.0e0/
      data  air( 49), air( 50), air( 51), air( 52), air( 53), air( 54) &
          / 196.0e0,  196.0e0,  236.0e0,  235.0e0,  229.0e0,  243.0e0/
      data  air( 55), air( 56), air( 57), air( 58), air( 59), air( 60) &
          / 264.0e0,  272.0e0,  237.0e0,  211.0e0,  180.0e0,  201.0e0/
      data  air( 61), air( 62), air( 63), air( 64), air( 65), air( 66) &
          / 204.0e0,  188.0e0,  235.0e0,  227.0e0,  234.0e0,  264.0e0/
      data  air( 67), air( 68), air( 69), air( 70), air( 71), air( 72) &
          / 302.0e0,  293.0e0,  259.0e0,  229.0e0,  203.0e0,  229.0e0/
      data  air( 73), air( 74), air( 75), air( 76), air( 77), air( 78) &
          / 242.0e0,  233.0e0,  267.0e0,  269.0e0,  270.0e0,  315.0e0/
      data  air( 79), air( 80), air( 81), air( 82), air( 83), air( 84) &
          / 364.0e0,  347.0e0,  312.0e0,  274.0e0,  237.0e0,  278.0e0/
      data  air( 85), air( 86), air( 87), air( 88), air( 89), air( 90) &
          / 284.0e0,  277.0e0,  317.0e0,  313.0e0,  318.0e0,  374.0e0/
      data  air( 91), air( 92), air( 93), air( 94), air( 95), air( 96) &
          / 413.0e0,  405.0e0,  355.0e0,  306.0e0,  271.0e0,  306.0e0/
      data  air( 97), air( 98), air( 99), air(100), air(101), air(102) &
          / 315.0e0,  301.0e0,  356.0e0,  348.0e0,  355.0e0,  422.0e0/
      data  air(103), air(104), air(105), air(106), air(107), air(108) &
          / 465.0e0,  467.0e0,  404.0e0,  347.0e0,  305.0e0,  336.0e0/
      data  air(109), air(110), air(111), air(112), air(113), air(114) &
          / 340.0e0,  318.0e0,  362.0e0,  348.0e0,  363.0e0,  435.0e0/
      data  air(115), air(116), air(117), air(118), air(119), air(120) &
          / 491.0e0,  505.0e0,  404.0e0,  359.0e0,  310.0e0,  337.0e0/
      data  air(121), air(122), air(123), air(124), air(125), air(126) &
          / 360.0e0,  342.0e0,  406.0e0,  396.0e0,  420.0e0,  472.0e0/
      data  air(127), air(128), air(129), air(130), air(131), air(132) &
          / 548.0e0,  559.0e0,  463.0e0,  407.0e0,  362.0e0,  405.0e0/
      data  air(133), air(134), air(135), air(136), air(137), air(138) &
          / 417.0e0,  391.0e0,  419.0e0,  461.0e0,  472.0e0,  535.0e0/
      data  air(139), air(140), air(141), air(142), air(143), air(144) &
          / 622.0e0,  606.0e0,  508.0e0,  461.0e0,  390.0e0,  432.0e0/
!
!  define constants
!
      iym = 12
      n = 12
      m = 12
!
!  write header
!
      write ( *, 1000)
!
!  run simple test of mpp
!
      write ( *, 1100)
      call mpp(ym, x, n, m, iym)
      write ( *,2000) ierr

      return

 1000 format ('1*ch2')
 1100 format (' simple test of mpp')
 2000 format (/' the value of ierr is ', i4)
end
subroutine xxch3 ( )

!*****************************************************************************80
!
!! XXCH3 tests the normal random number generator family of routines.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson and John Koontz,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     integer i
!        an index variable.
!     integer ierr
!        flag to indicate presence of error detected by preceding
!        starpac call.  (0 is ok, 1 is error)
!     integer iseed
!        the seed for the random number generator.
!     integer iym
!        the exact value of the first dimension of array ym.
!     integer ldstak
!        a dummy variable for this test subprogram.
!     integer m
!        the number of sets of numbers to be generated
!     integer n
!        the number of observations to be generated.
!     real sigma
!        the s.d. of the sample.
!     real ym(50,2)
!        data vector for tests.
!     real ymean
!        the mean of the sample.
!
  implicit none

  integer ierr
!
!  local scalars
      real &
         sigma,ymean
      integer &
         i,iseed,iym,m,n
!
!  local arrays
      real &
         ym(50,2)
!
!  external subroutines
      external iprint,mvp,nrand,nrandc
!
!  common blocks
      common /errchk/ierr

!
!     data initialization
!
      iym = 50
      iseed = 531
      n = 50
      m = 2
      ymean = 4.0e0
      sigma = 0.5e0
!
!
!  write heading
!
      write ( *,1000)
!
!     generate standard normal pseudo-random numbers into column 1 of ym
!
      write ( *,1100)
      call nrand(ym(1,1), n, iseed)
      write ( *,2000) ierr
      write ( *, 1400) (ym(i,1),i=1,n)
!
!     generate normal pseudo-random numbers
!     with mean 4.0 and standard deviation 0.5 into column 2 of ym
!
      write ( *,1200)
      call nrandc(ym(1,2), n, iseed, ymean, sigma)
      write ( *,2000) ierr
      write ( *, 1400) (ym(i,2),i=1,n)
!
!     plot results, sampling every observation
!
      write ( *,1500)
      call mvp (ym, n, m, iym, 1)
!
      return

 1000 format ('1*ch3')
 1100 format (' simple test of nrand')
 1200 format ('1simple test of nrandc')
 2000 format (/' the value of ierr is ', i4)
 1400 format (/' generated results = '//(5e15.8))
 1500 format ('1mvp display of generated results')
end
subroutine xxch4 ( ldstak )

!*****************************************************************************80
!
!! XXCH4 tests the histogram family of routines.
!
!  Discussion:
!
!    Data set is from page 39 of mandel [1964]
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson and John Koontz,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer ierr
!        error flag
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer n
!        the length of the vector y.
!     real y(40)
!        data vector for tests.
!
  implicit none

      integer &
         ldstak
!
!  scalars in common
      integer &
         ierr
!
!  arrays in common
      double precision dstak(12)
!
!  local scalars
      integer n
!
!  local arrays
      real &
         y(40)
!
!  external subroutines
      external hist,iprint
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr

!
!     data initializations.
!
      data n /39/
!
      data y( 1), y( 2), y( 3), y( 4) &
          / 0.4, 0.6, 1.0, 1.0/
      data y( 5), y( 6), y( 7), y( 8) &
          / 1.0, 0.5, 0.6, 0.7/
      data y( 9), y(10), y(11), y(12) &
          / 1.0, 0.6, 0.2, 1.9/
      data y(13), y(14), y(15), y(16) &
          / 0.2, 0.4, 0.0, -0.4/
      data y(17), y(18), y(19), y(20) &
          / -0.3, 0.0, -0.4, -0.3/
      data y(21), y(22), y(23), y(24) &
          / 0.1, -0.1, 0.2, -0.5/
      data y(25), y(26), y(27), y(28) &
          / 0.3, -0.1, 0.2, -0.2/
      data y(29), y(30), y(31), y(32) &
          / 0.8, 0.5, 0.6, 0.8/
      data y(33), y(34), y(35), y(36) &
          / 0.7, 0.7, 0.2, 0.5/
      data y(37), y(38), y(39) &
          / 0.7, 0.8, 1.1/
!
!     print heading
!
      write ( *,1000)
!
!     perform simple test of hist
!
      write ( *,1100)
      call hist(y, n, ldstak)
      write ( *,2000) ierr
!
      return
!
!     formats
!
 1000 format ('1*ch4')
 1100 format (' simple test of hist')
 2000 format (/' the value of ierr is ', i4)
end
subroutine xxch5 ( ldstak )

!*****************************************************************************80
!
!! XXCH5 tests the statistical analysis family of routines.
!
!  Discussion:
!
!    Data set is from page 39 of mandel [1964]
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson, John Koontz,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer ierr
!        error flag
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer n
!        the length of the vector y.
!     real y(40)
!        data vector for tests.
!
  implicit none

      integer &
         ldstak
!
!  scalars in common
      integer &
         ierr
!
!  arrays in common
      double precision dstak(12)
!
!  local scalars
      integer n
!
!  local arrays
      real &
         y(40)
!
!  external subroutines
      external iprint,stat
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr
!
!     data initializations.
!
      data n /39/
!
      data y( 1), y( 2), y( 3), y( 4) &
          / 0.4, 0.6, 1.0, 1.0/
      data y( 5), y( 6), y( 7), y( 8) &
          / 1.0, 0.5, 0.6, 0.7/
      data y( 9), y(10), y(11), y(12) &
          / 1.0, 0.6, 0.2, 1.9/
      data y(13), y(14), y(15), y(16) &
          / 0.2, 0.4, 0.0, -0.4/
      data y(17), y(18), y(19), y(20) &
          / -0.3, 0.0, -0.4, -0.3/
      data y(21), y(22), y(23), y(24) &
          / 0.1, -0.1, 0.2, -0.5/
      data y(25), y(26), y(27), y(28) &
          / 0.3, -0.1, 0.2, -0.2/
      data y(29), y(30), y(31), y(32) &
          / 0.8, 0.5, 0.6, 0.8/
      data y(33), y(34), y(35), y(36) &
          / 0.7, 0.7, 0.2, 0.5/
      data y(37), y(38), y(39) &
          / 0.7, 0.8, 1.1/
!
!     print heading
!
      write ( *,1000)
!
!     perform simple test of stat
!
      write ( *,1100)
      call stat(y, n, ldstak)
      write ( *,2000) ierr
!
      return
!
!     formats
!
 1000 format ('1*ch5')
 1100 format (' simple test of stat')
 2000 format (/' the value of ierr is ', i4)
end
subroutine xxch6 ( ldstak )

!*****************************************************************************80
!
!! XXCH6 tests the oneway analysis of variance family of routines.
!
!  Discussion:
!
!    Data set is from pages 314-316 of Brownlee [1965].
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson and John Koontz,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer ierr
!        error flag
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer n
!        the length of the vector y.
!     real tag(20)
!        the tag values for each observation
!     real y(20)
!        data vector for tests.
!
  implicit none

      integer &
         ldstak
!
!  scalars in common
      integer &
         ierr
!
!  arrays in common
      double precision dstak(12)
!
!  local scalars
      integer n
!
!  local arrays
      real &
         tag(20),y(20)
!
!  external subroutines
      external aov1,iprint
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr
!
!     data initializations.
!
      data n /16/

      data y( 1), y( 2), y( 3), y( 4) &
          / 83.0, 81.0, 76.0, 78.0/
      data y( 5), y( 6), y( 7), y( 8) &
          / 79.0, 72.0, 61.0, 61.0/
      data y( 9), y(10), y(11), y(12) &
          / 67.0, 67.0, 64.0, 78.0/
      data y(13), y(14), y(15), y(16) &
          / 71.0, 75.0, 72.0, 74.0/

      data tag( 1), tag( 2), tag( 3), tag( 4) &
          / 1.0, 1.0, 1.0, 1.0/
      data tag( 5), tag( 6), tag( 7), tag( 8) &
          / 1.0, 1.0, 2.0, 2.0/
      data tag( 9), tag(10), tag(11), tag(12) &
          / 2.0, 2.0, 2.0, 3.0/
      data tag(13), tag(14), tag(15), tag(16) &
          / 3.0, 3.0, 3.0, 3.0/
!
!
!  print heading
!
      write ( *,1000)
!
!     perform simple test of aov1
!
      write ( *,1100)
      call aov1(y, tag, n, ldstak)
      write ( *,2000) ierr
!
      return
!
!     formats
!
 1000 format ('1*ch6')
 1100 format (' simple test of aov1')
 2000 format (/' the value of ierr is ', i4)
end
subroutine xxch7 ( ldstak )

!*****************************************************************************80
!
!! XXCH7 tests the correlation analysis family of routines.
!
!  Discussion:
!
!    Data is from draper and smith [1968], page 216.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    Draper, Smith,
!    Applied Regression Analysis
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer ierr
!        the integer value designating whether any errors were
!        detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer iym
!        the first dimension of the array ym.
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer m
!        the number of variables measured for each observation.
!     integer n
!        the number of observations.
!     real ym(10,5)
!        the observed multivariate data.
!
!
  implicit none

      integer &
         ldstak
!
!  scalars in common
      integer &
         ierr
!
!  arrays in common
      double precision dstak(12)
!
!  local scalars
      integer iym,m,n
!
!  local arrays
      real &
         ym(10,5)
!
!  external subroutines
      external corr,iprint
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr

      data     ym(1,1),   ym(1,2),   ym(1,3),   ym(1,4) &
          /      42.2e0,  11.2e0,  31.9e0, 167.1e0/
      data     ym(2,1),   ym(2,2),   ym(2,3),   ym(2,4) &
          /      48.6e0,  10.6e0,  13.2e0, 174.4e0/
      data     ym(3,1),   ym(3,2),   ym(3,3),   ym(3,4) &
          /      42.6e0,  10.6e0,  28.7e0, 160.8e0/
      data     ym(4,1),   ym(4,2),   ym(4,3),   ym(4,4) &
          /      39.0e0,  10.4e0,  26.1e0, 162.0e0/
      data     ym(5,1),   ym(5,2),   ym(5,3),   ym(5,4) &
          /      34.7e0,   9.3e0,  30.1e0, 140.8e0/
      data     ym(6,1),   ym(6,2),   ym(6,3),   ym(6,4) &
          /      44.5e0,  10.8e0,   8.5e0, 174.6e0/
      data     ym(7,1),   ym(7,2),   ym(7,3),   ym(7,4) &
          /      39.1e0,  10.7e0,  24.3e0, 163.7e0/
      data     ym(8,1),   ym(8,2),   ym(8,3),   ym(8,4) &
          /      40.1e0,  10.0e0,  18.6e0, 174.5e0/
      data     ym(9,1),   ym(9,2),   ym(9,3),   ym(9,4) &
          /      45.9e0,  12.0e0,  20.4e0, 185.7e0/
!
!
!     set parameters necessary for the computations
!
      iym = 10
      n = 9
      m = 4
!
!     print header
!
      write ( *,1000)
!
!     run simple example of corr
!
      write ( *,1100)
      call corr(ym, n, m, iym, ldstak)
      write ( *,2000) ierr

      return
!
!     format statements
!
 1000 format ('1*ch7')
 1100 format (' simple test of corr')
 2000 format (/' the value of ierr is ', i4)

end
subroutine xxch8 ( ldstak )

!*****************************************************************************80
!
!! XXCH8 tests the linear least squares family of routines.
!
!  Discussion:
!
!    LLS problem is from daniel and wood [1971], pages 61-65.
!
!    LLSP problem is from miller and freund [1977], page 311.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer ierr
!        the integer value designating whether any errors were
!        detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ixm
!        the first dimension of the matrix x.
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer ndeg
!        the degree of the polynomial model to be fit.
!     integer npar
!        the number of parameters to be estimated.
!     integer n1, n2
!        the number of observations in each problem.
!     real res(25)
!        the residuals.
!     real x(25)
!        the independent variable.
!     real xm(25,5)
!        the independent variable.
!     real y1(25), y2(25)
!        the dependent variable.
!
  implicit none

      integer &
         ldstak
!
!  scalars in common
      integer &
         ierr
!
!  arrays in common
      double precision dstak(12)
!
!  local scalars
      integer ixm,n1,n2,ndeg,npar
!
!  local arrays
      real &
         res(25),x(25),xm(25,5),y1(25),y2(25)
!
!  external subroutines
      external iprint,lls,llsp
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr

      data      xm(1,1),  xm(1,2),  xm(1,3),  xm(1,4) &
          /      1.0e0, 80.0e0, 27.0e0, 89.0e0/
      data      xm(2,1),  xm(2,2),  xm(2,3),  xm(2,4) &
          /      1.0e0, 80.0e0, 27.0e0, 88.0e0/
      data      xm(3,1),  xm(3,2),  xm(3,3),  xm(3,4) &
          /      1.0e0, 75.0e0, 25.0e0, 90.0e0/
      data      xm(4,1),  xm(4,2),  xm(4,3),  xm(4,4) &
          /      1.0e0, 62.0e0, 24.0e0, 87.0e0/
      data      xm(5,1),  xm(5,2),  xm(5,3),  xm(5,4) &
          /      1.0e0, 62.0e0, 22.0e0, 87.0e0/
      data      xm(6,1),  xm(6,2),  xm(6,3),  xm(6,4) &
          /      1.0e0, 62.0e0, 23.0e0, 87.0e0/
      data      xm(7,1),  xm(7,2),  xm(7,3),  xm(7,4) &
          /      1.0e0, 62.0e0, 24.0e0, 93.0e0/
      data      xm(8,1),  xm(8,2),  xm(8,3),  xm(8,4) &
          /      1.0e0, 62.0e0, 24.0e0, 93.0e0/
      data      xm(9,1),  xm(9,2),  xm(9,3),  xm(9,4) &
          /      1.0e0, 58.0e0, 23.0e0, 87.0e0/
      data     xm(10,1), xm(10,2), xm(10,3), xm(10,4) &
          /      1.0e0, 58.0e0, 18.0e0, 80.0e0/
      data     xm(11,1), xm(11,2), xm(11,3), xm(11,4) &
          /      1.0e0, 58.0e0, 18.0e0, 89.0e0/
      data     xm(12,1), xm(12,2), xm(12,3), xm(12,4) &
          /      1.0e0, 58.0e0, 17.0e0, 88.0e0/
      data     xm(13,1), xm(13,2), xm(13,3), xm(13,4) &
          /      1.0e0, 58.0e0, 18.0e0, 82.0e0/
      data     xm(14,1), xm(14,2), xm(14,3), xm(14,4) &
          /      1.0e0, 58.0e0, 19.0e0, 93.0e0/
      data     xm(15,1), xm(15,2), xm(15,3), xm(15,4) &
          /      1.0e0, 50.0e0, 18.0e0, 89.0e0/
      data     xm(16,1), xm(16,2), xm(16,3), xm(16,4) &
          /      1.0e0, 50.0e0, 18.0e0, 86.0e0/
      data     xm(17,1), xm(17,2), xm(17,3), xm(17,4) &
          /      1.0e0, 50.0e0, 19.0e0, 72.0e0/
      data     xm(18,1), xm(18,2), xm(18,3), xm(18,4) &
          /      1.0e0, 50.0e0, 19.0e0, 79.0e0/
      data     xm(19,1), xm(19,2), xm(19,3), xm(19,4) &
          /      1.0e0, 50.0e0, 20.0e0, 80.0e0/
      data     xm(20,1), xm(20,2), xm(20,3), xm(20,4) &
          /      1.0e0, 56.0e0, 20.0e0, 82.0e0/
      data     xm(21,1), xm(21,2), xm(21,3), xm(21,4) &
          /      1.0e0, 70.0e0, 20.0e0, 91.0e0/
!
      data        y1(1),    y1(2),    y1(3) &
          /     42.0e0, 37.0e0, 37.0e0/
      data        y1(4),    y1(5),    y1(6) &
          /     28.0e0, 18.0e0, 18.0e0/
      data        y1(7),    y1(8),    y1(9) &
          /     19.0e0, 20.0e0, 15.0e0/
      data       y1(10),   y1(11),   y1(12) &
          /     14.0e0, 14.0e0, 13.0e0/
      data       y1(13),   y1(14),   y1(15) &
          /     11.0e0, 12.0e0,  8.0e0/
      data       y1(16),   y1(17),   y1(18) &
          /      7.0e0,  8.0e0,  8.0e0/
      data       y1(19),   y1(20),   y1(21) &
          /      9.0e0, 15.0e0, 15.0e0/
!
      data         x(1),     x(2),     x(3) &
          /      0.0e0,  1.0e0,  2.0e0/
      data         x(4),     x(5),     x(6) &
          /      3.0e0,  4.0e0,  5.0e0/
      data         x(7),     x(8),     x(9) &
          /      6.0e0,  7.0e0,  8.0e0/

      data        y2(1),    y2(2),    y2(3) &
          /     12.0e0, 10.5e0, 10.0e0/
      data        y2(4),    y2(5),    y2(6) &
          /      8.0e0,  7.0e0,  8.0e0/
      data        y2(7),    y2(8),    y2(9) &
          /      7.5e0,  8.5e0,  9.0e0/
!
!     set parameters necessary for the computations
!
      ixm = 25
      n1 = 21
      n2 = 9
      npar = 4
      ndeg = 2
!
!     print header
!
      write ( *,1000)
!
!     run simple example of lls
!
      write ( *,1100)
      call lls(y1, xm, n1, ixm, npar, res, ldstak)
      write ( *,2000) ierr
!
!     run simple example of llsp
!
      write ( *,1200)
      call llsp(y2, x, n2, ndeg, res, ldstak)
      write ( *,2000) ierr

      return

 1000 format ('1*ch8')
 1100 format (' simple test of lls')
 1200 format ('1simple test of llsp')
 2000 format (/' the value of ierr is ', i4)
end
subroutine xxch9 ( ldstak )

!*****************************************************************************80
!
!! XXCH9 tests the nonlinear least squares family of routines.
!
!  Discussion:
!
!    Data is from daniel and wood [1980], pages 428-441.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     external drv1a, drv1b
!        the name of the ''user supplied'' derivative routines.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer ierr
!        the integer value designating whether any errors were
!        detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ixm
!        the first dimension of the matrix x.
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer m
!        the number of independent variables.
!     external mdl1
!        the name of the ''user supplied'' model routines.
!     integer n
!        the number of observations in each problem.
!     integer npar
!        the number of parameters to be estimated.
!     real par(5)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real res(10)
!        the residuals.
!     real stp(5)
!        the step sizes selected for generating finite difference
!        derivatives.
!     real xm(10,2)
!        the independent variable.
!     real y(10)
!        the dependent variable.
!
  implicit none
!
      integer &
         ldstak
!
!  scalars in common
      integer &
         ierr
!
!  arrays in common
      double precision dstak(12)
!
!  local scalars
      integer ixm,m,n,npar
!
!  local arrays
      real &
         par(5),res(10),stp(5),xm(10,2),y(10)
!
!  external subroutines
      external dckls,drv1a,drv1b,iprint,mdl1,nls,nlsd,stpls
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr

      data y(1), y(2), y(3), y(4), y(5), y(6) &
         /2.138e0, 3.421e0, 3.597e0, 4.340e0, 4.882e0, 5.660e0/

      data xm(1,1), xm(2,1), xm(3,1), xm(4,1), xm(5,1), xm(6,1) &
         /1.309e0, 1.471e0, 1.490e0, 1.565e0, 1.611e0, 1.680e0/
!
!  set parameters necessary for the computations
!
      ixm = 10
      n = 6
      m = 1
      npar = 2
!
!  print header
!
      write ( *,1000)
!
!  run simple example of nls
!
      write ( *,1100)
      par(1) = 0.725
      par(2) = 4.000
      call nls(y, xm, n, m, ixm, mdl1, par, npar, res, ldstak)
      write ( *,2000) ierr
!
!  run simple example of nlsd
!
      write ( *,1200)
      par(1) = 0.725
      par(2) = 4.000
      call nlsd(y, xm, n, m, ixm, mdl1, drv1a, par, npar, res, ldstak)
      write ( *,2000) ierr
!
!  run simple example of stpls
!
      write ( *,1300)
      par(1) = 0.725
      par(2) = 4.000
      call stpls(xm, n, m, ixm, mdl1, par, npar, ldstak, stp)
      write ( *,2000) ierr
!
!  run simple example of dckls
!
      write ( *,1400)
      par(1) = 0.000
      par(2) = 4.000
      call dckls(xm, n, m, ixm, mdl1, drv1b, par, npar, ldstak)
      write ( *,2000) ierr

      return

 1000 format ('1*ch9')
 1100 format (' simple test of nls')
 1200 format ('1simple test of nlsd')
 1300 format ('1simple test of stpls')
 1400 format ('1simple test of dckls')
 2000 format (/' the value of ierr is ', i4)
end
subroutine xxch10 ( )

!*****************************************************************************80
!
!! XXCH10 tests the histogram family of routines.
!
!  Discussion:
!
!    Data is the airline data listed on page 531 of box and jenkins.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    George Box, Gwilym Jenkins,
!    Time Series Analysis: Forecasting and Control,
!    Holden-Day, 1970,
!    QA280.B67.
!
!  Parameters:
!
!     real air(144)
!        the airline data.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index variable.
!     integer iar
!        the number of coefficients in the difference filter.
!     integer ierr
!        a common variable used as a flag to indicate whether
!        or not there are any errors, if =0 then no errors.
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer n
!        the number of observations.
!     integer nyf
!        the number of observations in the filtered series.
!     real phi(5)
!        the filter coefficients.
!     real y(150)
!        the log of the airline data.
!     real yf(150)
!        the filtered data.
!
!
  implicit none
!
!  scalars in common
      integer &
         ierr
!
!  arrays in common
      double precision dstak(12)
!
!  local scalars
      integer iar,n,nyf
!
!  local arrays
      real &
         air(144),phi(5),y(150),yf(150)
!
!  external subroutines
      external dif,gfarf,iprint,vp
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr

      data  air(  1), air(  2), air(  3), air(  4), air(  5), air(  6) &
          / 112.0e0,  118.0e0,  132.0e0,  129.0e0,  121.0e0,  135.0e0/
      data  air(  7), air(  8), air(  9), air( 10), air( 11), air( 12) &
          / 148.0e0,  148.0e0,  136.0e0,  119.0e0,  104.0e0,  118.0e0/
      data  air( 13), air( 14), air( 15), air( 16), air( 17), air( 18) &
          / 115.0e0,  126.0e0,  141.0e0,  135.0e0,  125.0e0,  149.0e0/
      data  air( 19), air( 20), air( 21), air( 22), air( 23), air( 24) &
          / 170.0e0,  170.0e0,  158.0e0,  133.0e0,  114.0e0,  140.0e0/
      data  air( 25), air( 26), air( 27), air( 28), air( 29), air( 30) &
          / 145.0e0,  150.0e0,  178.0e0,  163.0e0,  172.0e0,  178.0e0/
      data  air( 31), air( 32), air( 33), air( 34), air( 35), air( 36) &
          / 199.0e0,  199.0e0,  184.0e0,  162.0e0,  146.0e0,  166.0e0/
      data  air( 37), air( 38), air( 39), air( 40), air( 41), air( 42) &
          / 171.0e0,  180.0e0,  193.0e0,  181.0e0,  183.0e0,  218.0e0/
      data  air( 43), air( 44), air( 45), air( 46), air( 47), air( 48) &
          / 230.0e0,  242.0e0,  209.0e0,  191.0e0,  172.0e0,  194.0e0/
      data  air( 49), air( 50), air( 51), air( 52), air( 53), air( 54) &
          / 196.0e0,  196.0e0,  236.0e0,  235.0e0,  229.0e0,  243.0e0/
      data  air( 55), air( 56), air( 57), air( 58), air( 59), air( 60) &
          / 264.0e0,  272.0e0,  237.0e0,  211.0e0,  180.0e0,  201.0e0/
      data  air( 61), air( 62), air( 63), air( 64), air( 65), air( 66) &
          / 204.0e0,  188.0e0,  235.0e0,  227.0e0,  234.0e0,  264.0e0/
      data  air( 67), air( 68), air( 69), air( 70), air( 71), air( 72) &
          / 302.0e0,  293.0e0,  259.0e0,  229.0e0,  203.0e0,  229.0e0/
      data  air( 73), air( 74), air( 75), air( 76), air( 77), air( 78) &
          / 242.0e0,  233.0e0,  267.0e0,  269.0e0,  270.0e0,  315.0e0/
      data  air( 79), air( 80), air( 81), air( 82), air( 83), air( 84) &
          / 364.0e0,  347.0e0,  312.0e0,  274.0e0,  237.0e0,  278.0e0/
      data  air( 85), air( 86), air( 87), air( 88), air( 89), air( 90) &
          / 284.0e0,  277.0e0,  317.0e0,  313.0e0,  318.0e0,  374.0e0/
      data  air( 91), air( 92), air( 93), air( 94), air( 95), air( 96) &
          / 413.0e0,  405.0e0,  355.0e0,  306.0e0,  271.0e0,  306.0e0/
      data  air( 97), air( 98), air( 99), air(100), air(101), air(102) &
          / 315.0e0,  301.0e0,  356.0e0,  348.0e0,  355.0e0,  422.0e0/
      data  air(103), air(104), air(105), air(106), air(107), air(108) &
          / 465.0e0,  467.0e0,  404.0e0,  347.0e0,  305.0e0,  336.0e0/
      data  air(109), air(110), air(111), air(112), air(113), air(114) &
          / 340.0e0,  318.0e0,  362.0e0,  348.0e0,  363.0e0,  435.0e0/
      data  air(115), air(116), air(117), air(118), air(119), air(120) &
          / 491.0e0,  505.0e0,  404.0e0,  359.0e0,  310.0e0,  337.0e0/
      data  air(121), air(122), air(123), air(124), air(125), air(126) &
          / 360.0e0,  342.0e0,  406.0e0,  396.0e0,  420.0e0,  472.0e0/
      data  air(127), air(128), air(129), air(130), air(131), air(132) &
          / 548.0e0,  559.0e0,  463.0e0,  407.0e0,  362.0e0,  405.0e0/
      data  air(133), air(134), air(135), air(136), air(137), air(138) &
          / 417.0e0,  391.0e0,  419.0e0,  461.0e0,  472.0e0,  535.0e0/
      data  air(139), air(140), air(141), air(142), air(143), air(144) &
          / 622.0e0,  606.0e0,  508.0e0,  461.0e0,  390.0e0,  432.0e0/
!
!  define constants
!
      n = 144
!
!  take log of data
!
      y(1:n) = log ( air(1:n) )
!
!  write header
!
      write ( *, 1000)
!
!  run simple test of dif
!
      write ( *, 1100)
      call dif (y, n, yf, nyf)
      write ( *,2000) ierr
!
!  plot original series
!
      write ( *, 1200)
      call vp (y, n, 1)
      write ( *,2000) ierr
!
!  plot differenced series
!
      write ( *, 1300)
      call vp (yf, nyf, 1)
      write ( *,2000) ierr
!
!  run simple test of gfarf
!
      write ( *, 1400)
      phi(1) = 1.0
      iar = 1
      call gfarf (phi, iar)
      write ( *,2000) ierr
!
      return
!
!     format statements
!
 1000 format ('1*ch10')
 1100 format (' simple test of dif (no output unless error found)')
 1200 format ('1plot of original series')
 1300 format ('1plot of differenced series')
 1400 format ('1simple test of gfarf')
 2000 format (/' the value of ierr is ', i4)
end
subroutine xxch11 ( ldstak )

!*****************************************************************************80
!
!! XXCH10 tests the complex demodulation family of routines.
!
!  Discussion:
!
!    Data is the wolf sunspot numbers for the years 1700 to 1960 as
!    tabulated by waldmeier [1961].
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    Max Waldmeier,
!    The Sunspot-Activity in the Years 1610-1960,
!    Shulthess, Zurich, 1961.
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fc
!        the cutoff frequency used for the low pass filter.
!     real fd
!        the demodulation frequency.
!     integer ierr
!        a common variable used as a flag to indicate whether
!        or not there are any errors, if =0 then no errors.
!     integer k
!        the number of terms in the symetric linear filter.
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer n
!        the number of observations.
!     real y(300)
!        the log of the airline data.
!
!
  implicit none

      integer &
         ldstak
!
!  scalars in common
      integer &
         ierr
!
!  arrays in common
      double precision dstak(12)
!
!  local scalars
      real &
         fc,fd
      integer k,n
!
!  local arrays
      real &
         y(300)
!
!  external subroutines
      external demod,iprint
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr

      data   y(  1),  y(  2),  y(  3),  y(  4),  y(  5),  y(  6) &
          /     5.0e0, 11.0e0, 16.0e0, 23.0e0, 36.0e0, 58.0e0/
      data   y(  7),  y(  8),  y(  9),  y( 10),  y( 11),  y( 12) &
          /    29.0e0, 20.0e0, 10.0e0,  8.0e0,  3.0e0,  0.0e0/
      data   y( 13),  y( 14),  y( 15),  y( 16),  y( 17),  y( 18) &
          /     0.0e0, 2.0e0, 11.0e0, 27.0e0, 47.0e0, 63.0e0/
      data   y( 19),  y( 20),  y( 21),  y( 22),  y( 23),  y( 24) &
          /    60.0e0, 39.0e0, 28.0e0, 26.0e0, 22.0e0, 11.0e0/
      data   y( 25),  y( 26),  y( 27),  y( 28),  y( 29),  y( 30) &
          /    21.0e0, 40.0e0, 78.0e0,122.0e0,103.0e0, 73.0e0/
      data   y( 31),  y( 32),  y( 33),  y( 34),  y( 35),  y( 36) &
          /    47.0e0, 35.0e0, 11.0e0,  5.0e0, 16.0e0, 34.0e0/
      data   y( 37),  y( 38),  y( 39),  y( 40),  y( 41),  y( 42) &
          /    70.0e0, 81.0e0,111.0e0,101.0e0, 73.0e0, 40.0e0/
      data   y( 43),  y( 44),  y( 45),  y( 46),  y( 47),  y( 48) &
          /    20.0e0, 16.0e0,  5.0e0, 11.0e0, 22.0e0, 40.0e0/
      data   y( 49),  y( 50),  y( 51),  y( 52),  y( 53),  y( 54) &
          /    60.0e0, 80.9e0, 83.4e0, 47.7e0, 47.8e0, 30.7e0/
      data   y( 55),  y( 56),  y( 57),  y( 58),  y( 59),  y( 60) &
          /    12.2e0,  9.6e0, 10.2e0, 32.4e0, 47.6e0, 54.0e0/
      data   y( 61),  y( 62),  y( 63),  y( 64),  y( 65),  y( 66) &
          /    62.9e0, 85.9e0, 61.2e0, 45.1e0, 36.4e0, 20.9e0/
      data   y( 67),  y( 68),  y( 69),  y( 70),  y( 71),  y( 72) &
          /    11.4e0, 37.8e0, 69.8e0,106.1e0,100.8e0, 81.6e0/
      data   y( 73),  y( 74),  y( 75),  y( 76),  y( 77),  y( 78) &
          /    66.5e0, 34.8e0, 30.6e0,  7.0e0, 19.8e0, 92.5e0/
      data   y( 79),  y( 80),  y( 81),  y( 82),  y( 83),  y( 84) &
          /   154.4e0,125.9e0, 84.8e0, 68.1e0, 38.5e0, 22.8e0/
      data   y( 85),  y( 86),  y( 87),  y( 88),  y( 89),  y( 90) &
          /    10.2e0, 24.1e0, 82.9e0,132.0e0,130.9e0,118.1e0/
      data   y( 91),  y( 92),  y( 93),  y( 94),  y( 95),  y( 96) &
          /    89.9e0, 66.6e0, 60.0e0, 46.9e0, 41.0e0, 21.3e0/
      data   y( 97),  y( 98),  y( 99),  y(100),  y(101),  y(102) &
          /    16.0e0,  6.4e0,  4.1e0,  6.8e0, 14.5e0, 34.0e0/
      data   y(103),  y(104),  y(105),  y(106),  y(107),  y(108) &
          /    45.0e0, 43.1e0, 47.5e0, 42.2e0, 28.1e0, 10.1e0/
      data   y(109),  y(110),  y(111),  y(112),  y(113),  y(114) &
          /     8.1e0,  2.5e0,  0.0e0,  1.4e0,  5.0e0, 12.2e0/
      data   y(115),  y(116),  y(117),  y(118),  y(119),  y(120) &
          /    13.9e0, 35.4e0, 45.8e0, 41.1e0, 30.1e0, 23.9e0/
      data   y(121),  y(122),  y(123),  y(124),  y(125),  y(126) &
          /    15.6e0,  6.6e0,  4.0e0,  1.8e0,  8.5e0, 16.6e0/
      data   y(127),  y(128),  y(129),  y(130),  y(131),  y(132) &
          /    36.3e0, 49.6e0, 64.2e0, 67.0e0, 70.9e0, 47.8e0/
      data   y(133),  y(134),  y(135),  y(136),  y(137),  y(138) &
          /    27.5e0,  8.5e0, 13.2e0, 56.9e0,121.5e0,138.3e0/
      data   y(139),  y(140),  y(141),  y(142),  y(143),  y(144) &
          /   103.2e0, 85.7e0, 64.6e0, 36.7e0, 24.2e0, 10.7e0/
      data   y(145),  y(146),  y(147),  y(148),  y(149),  y(150) &
          /    15.0e0, 40.1e0, 61.5e0, 98.5e0,124.7e0, 96.3e0/
      data   y(151),  y(152),  y(153),  y(154),  y(155),  y(156) &
          /    66.6e0, 64.5e0, 54.1e0, 39.0e0, 20.6e0,  6.7e0/
      data   y(157),  y(158),  y(159),  y(160),  y(161),  y(162) &
          /     4.3e0, 22.7e0, 54.8e0, 93.8e0, 95.8e0, 77.2e0/
      data   y(163),  y(164),  y(165),  y(166),  y(167),  y(168) &
          /    59.1e0, 44.0e0, 47.0e0, 30.5e0, 16.3e0,  7.3e0/
      data   y(169),  y(170),  y(171),  y(172),  y(173),  y(174) &
          /    37.6e0, 74.0e0,139.0e0,111.2e0,101.6e0, 66.2e0/
      data   y(175),  y(176),  y(177),  y(178),  y(179),  y(180) &
          /    44.7e0, 17.0e0, 11.3e0, 12.4e0,  3.4e0,  6.0e0/
      data   y(181),  y(182),  y(183),  y(184),  y(185),  y(186) &
          /    32.3e0, 54.3e0, 59.7e0, 63.7e0, 63.5e0, 52.2e0/
      data   y(187),  y(188),  y(189),  y(190),  y(191),  y(192) &
          /    25.4e0, 13.1e0,  6.8e0,  6.3e0,  7.1e0, 35.6e0/
      data   y(193),  y(194),  y(195),  y(196),  y(197),  y(198) &
          /    73.0e0, 85.1e0, 78.0e0, 64.0e0, 41.8e0, 26.2e0/
      data   y(199),  y(200),  y(201),  y(202),  y(203),  y(204) &
          /    26.7e0, 12.1e0,  9.5e0,  2.7e0,  5.0e0, 24.4e0/
      data   y(205),  y(206),  y(207),  y(208),  y(209),  y(210) &
          /    42.0e0, 63.5e0, 53.8e0, 62.0e0, 48.5e0, 43.9e0/
      data   y(211),  y(212),  y(213),  y(214),  y(215),  y(216) &
          /    18.6e0,  5.7e0,  3.6e0,  1.4e0,  9.6e0, 47.4e0/
      data   y(217),  y(218),  y(219),  y(220),  y(221),  y(222) &
          /    57.1e0,103.9e0, 80.6e0, 63.6e0, 37.6e0, 26.1e0/
      data   y(223),  y(224),  y(225),  y(226),  y(227),  y(228) &
          /    14.2e0,  5.8e0, 16.7e0, 44.3e0, 63.9e0, 69.0e0/
      data   y(229),  y(230),  y(231),  y(232),  y(233),  y(234) &
          /    77.8e0, 64.9e0, 35.7e0, 21.2e0, 11.1e0,  5.7e0/
      data   y(235),  y(236),  y(237),  y(238),  y(239),  y(240) &
          /     8.7e0, 36.1e0, 79.7e0,114.4e0,109.6e0, 88.8e0/
      data   y(241),  y(242),  y(243),  y(244),  y(245),  y(246) &
          /    67.8e0, 47.5e0, 30.6e0, 16.3e0,  9.6e0, 33.2e0/
      data   y(247),  y(248),  y(249),  y(250),  y(251),  y(252) &
          /    92.6e0,151.6e0,136.3e0,134.7e0, 83.9e0, 69.4e0/
      data   y(253),  y(254),  y(255),  y(256),  y(257),  y(258) &
          /    31.5e0, 13.9e0,  4.4e0, 38.0e0,141.7e0,190.2e0/
      data   y(259),  y(260),  y(261) &
          /   184.8e0,159.0e0,112.3e0/
!
!  define constants
!
      n = 261
      fd = 1.0/11.0
      fc = 1.0/22.0
      k = 41
!
!  write header
!
      write ( *, 1000)
!
!  run simple test of dif
!
      write ( *, 1100)
      call demod (y, n, fd, fc, k, ldstak)
      write ( *,2000) ierr

      return

 1000 format ('1*ch11')
 1100 format (' simple test of demod')
 2000 format (/' the value of ierr is ', i4)
end
subroutine xxch12 ( ldstak )

!*****************************************************************************80
!
!! XXCH12 tests time series correlation and spectrum analysis routines.
!
!  Discussion:
!
!    data for acf is taken from p. 362 of jenkins and watts [1968]
!
!    data for ccf is taken from p. 361 of jenkins and watts [1968]
!
!    data for ufs is taken from p. 318 of jenkins and watts [1968]
!
!    data for uas is taken from p. 318 of jenkins and watts [1968]
!
!    data for taper, pgms, mdflt and ppl is
!    the wolf sunspot numbers for the years 1700 to 1960 as
!    tabulated by waldmeier [1961].
!
!    data for bfs is taken from pp. 387-388 of jenkins and watts [1968]
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    Gwilym Jenkins, Donald Watts,
!    Spectral Analysis and its Applications,
!    Holden-Day 1968.
!
!    Max Waldmeier,
!    The Sunspot-Activity in the Years 1610-1960,
!    Shulthess, Zurich, 1961.
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real freq(300)
!        the frequencies at which the periodogram is computed.
!     integer ierr
!        a common variable used as a flag to indicate whether
!        or not there are any errors, if =0 then no errors.
!     integer iextnd
!        the indicator variable used to designate whether zero or the
!        series mean is to be used to extend the series.
!     integer ilog
!        the indicator variable used to designate whether the plot is
!        to have logarithmic axis or not.
!     integer kmd(10)
!        the vector of modified daniel filter lengths.
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer lfreq
!        the length of vector freq.
!     integer lper
!        the length of vector per.
!     integer lyfft
!        the length of vector yfft.
!     integer nf
!        the number of frequencies.
!     integer nfft
!        the extended series length for the fft.
!     integer nk
!        the number of daniel filters to apply.
!     integer nprt
!        the print control variable.
!     integer ny1, ny2, ny3, ny4, ny5, ny6
!        the number of observations.
!     real per(300)
!        the periodogram.
!     real perf(300)
!        the filtered periodogram.
!     real taperp
!        the percentage of the series to be tapered.
!     real yfft(600)
!        an array for the fft computations.
!     real y1(100)
!        the data from page 362 of jenkins and watts.
!     real y2a(100), y2b(100)
!        the data from page 361 of jenkins and watts.
!     real y3(50), y4(50)
!        the data from page 318 of jenkins and watts.
!     real y5(300)
!        the wolf sunspot data.
!     real y6a(100), y6b(100)
!        the data from page 387 and 388 of jenkins and watts.
!
!
  implicit none

      integer &
         ldstak
!
!  scalars in common
      integer &
         ierr
!
!  arrays in common
      double precision dstak(12)
!
!  local scalars
      real &
         taperp
      integer &
         iextnd,ilog,lfreq,lper,lyfft,nf,nfft,nk,nprt,ny1,ny2, &
         ny3,ny4,ny5,ny6
!
!  local arrays
      real &
         freq(300),per(300),perf(300),y1(100),y2a(100),y2b(100), &
         y3(50),y4(50),y5(300),y6a(100),y6b(100),yfft(600)
      integer &
         kmd(10)
!
!  external subroutines
      external acf,bfs,ccf,iprint,mdflt,pgms,ppl,taper,uas,ufs
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr
!
!  equivalences
      equivalence (y3(1),y4(1))
      equivalence (y2a(1),y6a(1))
      equivalence (y2b(1),y6b(1))

      data   y1(  1),  y1(  2),  y1(  3),  y1(  4),  y1(  5),  y1(  6) &
          / -2.07e0, -1.15e0,  0.69e0, -0.46e0, -1.49e0, -0.70e0/
      data   y1(  7),  y1(  8),  y1(  9),  y1( 10),  y1( 11),  y1( 12) &
          / -1.07e0, -0.69e0, -0.68e0,  1.27e0, -1.05e0, -0.05e0/
      data   y1( 13),  y1( 14),  y1( 15),  y1( 16),  y1( 17),  y1( 18) &
          / -0.84e0, -0.62e0, -0.49e0, -1.29e0, -0.49e0, -1.06e0/
      data   y1( 19),  y1( 20),  y1( 21),  y1( 22),  y1( 23),  y1( 24) &
          / -0.38e0, -0.52e0, -0.13e0,  1.30e0, -1.51e0, -0.43e0/
      data   y1( 25),  y1( 26),  y1( 27),  y1( 28),  y1( 29),  y1( 30) &
          / -1.33e0, -0.78e0,  0.31e0, -0.95e0, -0.90e0, -0.30e0/
      data   y1( 31),  y1( 32),  y1( 33),  y1( 34),  y1( 35),  y1( 36) &
          / -1.02e0, -0.53e0,  0.15e0,  1.40e0,  1.22e0,  0.59e0/
      data   y1( 37),  y1( 38),  y1( 39),  y1( 40),  y1( 41),  y1( 42) &
          /  0.70e0,  1.70e0,  2.78e0,  1.98e0,  1.39e0,  1.85e0/
      data   y1( 43),  y1( 44),  y1( 45),  y1( 46),  y1( 47),  y1( 48) &
          /  2.60e0,  0.51e0,  2.77e0,  1.16e0,  1.07e0, -0.48e0/
      data   y1( 49),  y1( 50),  y1( 51),  y1( 52),  y1( 53),  y1( 54) &
          / -0.52e0,  0.37e0,  0.00e0, -1.99e0, -1.75e0,  0.70e0/
      data   y1( 55),  y1( 56),  y1( 57),  y1( 58),  y1( 59),  y1( 60) &
          /  0.73e0,  1.16e0,  0.06e0, -0.02e0,  1.10e0, -0.35e0/
      data   y1( 61),  y1( 62),  y1( 63),  y1( 64),  y1( 65),  y1( 66) &
          / -1.67e0, -1.57e0,  1.16e0,  1.84e0,  3.35e0,  0.40e0/
      data   y1( 67),  y1( 68),  y1( 69),  y1( 70),  y1( 71),  y1( 72) &
          /  0.45e0,  1.30e0,  0.93e0,  1.17e0, -1.74e0, -1.28e0/
      data   y1( 73),  y1( 74),  y1( 75),  y1( 76),  y1( 77),  y1( 78) &
          / -0.07e0,  1.50e0,  0.53e0,  0.20e0, -0.42e0,  1.18e0/
      data   y1( 79),  y1( 80),  y1( 81),  y1( 82),  y1( 83),  y1( 84) &
          /  0.82e0,  1.50e0,  2.92e0,  1.18e0,  1.23e0,  3.16e0/
      data   y1( 85),  y1( 86),  y1( 87),  y1( 88),  y1( 89),  y1( 90) &
          /  0.79e0,  0.68e0,  1.14e0,  1.02e0,  1.02e0, -0.71e0/
      data   y1( 91),  y1( 92),  y1( 93),  y1( 94),  y1( 95),  y1( 96) &
          / -0.17e0, -1.50e0, -0.26e0, -0.38e0,  0.93e0, -0.33e0/
      data   y1( 97),  y1( 98),  y1( 99),  y1(100) &
          / -1.12e0, -2.95e0, -2.09e0, -1.11e0                    /
!
      data  y2a(  1), y2a(  2), y2a(  3), y2a(  4), y2a(  5), y2a(  6) &
          /-0.88e0, -0.16e0, -1.87e0, -1.12e0,  1.38e0,  2.13e0/
      data  y2a(  7), y2a(  8), y2a(  9), y2a( 10), y2a( 11), y2a( 12) &
          / 2.76e0,  0.56e0, -0.69e0, -1.79e0, -3.82e0, -2.38e0/
      data  y2a( 13), y2a( 14), y2a( 15), y2a( 16), y2a( 17), y2a( 18) &
          / 1.00e0,  0.70e0, -0.15e0,  0.98e0,  0.11e0, -0.35e0/
      data  y2a( 19), y2a( 20), y2a( 21), y2a( 22), y2a( 23), y2a( 24) &
          /-0.73e0,  0.89e0, -1.63e0, -0.44e0, -1.37e0, -1.71e0/
      data  y2a( 25), y2a( 26), y2a( 27), y2a( 28), y2a( 29), y2a( 30) &
          /-1.22e0, -2.00e0, -0.22e0,  0.38e0,  1.31e0,  0.71e0/
      data  y2a( 31), y2a( 32), y2a( 33), y2a( 34), y2a( 35), y2a( 36) &
          / 0.32e0,  0.48e0, -1.88e0, -0.94e0, -1.54e0, -0.13e0/
      data  y2a( 37), y2a( 38), y2a( 39), y2a( 40), y2a( 41), y2a( 42) &
          / 1.02e0,  0.02e0, -0.77e0,  0.11e0, -0.60e0, -0.52e0/
      data  y2a( 43), y2a( 44), y2a( 45), y2a( 46), y2a( 47), y2a( 48) &
          /-0.09e0,  1.23e0,  1.46e0,  0.61e0,  0.42e0,  2.16e0/
      data  y2a( 49), y2a( 50), y2a( 51), y2a( 52), y2a( 53), y2a( 54) &
          / 3.18e0,  2.10e0,  0.37e0, -0.24e0,  0.57e0, -0.53e0/
      data  y2a( 55), y2a( 56), y2a( 57), y2a( 58), y2a( 59), y2a( 60) &
          / 2.44e0,  1.02e0, -0.53e0, -2.49e0, -2.12e0, -1.04e0/
      data  y2a( 61), y2a( 62), y2a( 63), y2a( 64), y2a( 65), y2a( 66) &
          /-0.12e0, -1.88e0, -1.50e0,  1.54e0,  3.33e0,  3.08e0/
      data  y2a( 67), y2a( 68), y2a( 69), y2a( 70), y2a( 71), y2a( 72) &
          / 1.71e0,  0.79e0,  1.55e0,  0.89e0, -0.89e0, -1.18e0/
      data  y2a( 73), y2a( 74), y2a( 75), y2a( 76), y2a( 77), y2a( 78) &
          / 0.89e0,  1.71e0,  3.05e0,  0.15e0, -1.04e0,  0.12e0/
      data  y2a( 79), y2a( 80), y2a( 81), y2a( 82), y2a( 83), y2a( 84) &
          / 0.08e0,  0.11e0, -2.62e0, -1.28e0,  1.07e0,  3.20e0/
      data  y2a( 85), y2a( 86), y2a( 87), y2a( 88), y2a( 89), y2a( 90) &
          / 1.92e0,  0.53e0, -1.08e0,  0.49e0, -0.58e0,  0.17e0/
      data  y2a( 91), y2a( 92), y2a( 93), y2a( 94), y2a( 95), y2a( 96) &
          / 1.15e0, -0.97e0, -1.63e0,  1.14e0, -0.67e0, -0.88e0/
      data  y2a( 97), y2a( 98), y2a( 99), y2a(100) &
          /-0.07e0,  0.24e0,  0.55e0, -2.16e0/
!
      data  y2b(  1), y2b(  2), y2b(  3), y2b(  4), y2b(  5), y2b(  6) &
          / 0.79e0,  1.12e0, -1.10e0, -2.39e0, -1.75e0, -0.82e0/
      data  y2b(  7), y2b(  8), y2b(  9), y2b( 10), y2b( 11), y2b( 12) &
          /-0.36e0,  1.27e0,  1.75e0,  2.44e0,  0.36e0, -2.10e0/
      data  y2b( 13), y2b( 14), y2b( 15), y2b( 16), y2b( 17), y2b( 18) &
          /-1.93e0, -1.30e0, -1.75e0, -0.34e0,  0.74e0,  0.49e0/
      data  y2b( 19), y2b( 20), y2b( 21), y2b( 22), y2b( 23), y2b( 24) &
          / 0.70e0,  0.71e0,  0.09e0,  0.59e0,  1.54e0,  0.14e0/
      data  y2b( 25), y2b( 26), y2b( 27), y2b( 28), y2b( 29), y2b( 30) &
          / 0.55e0, -1.40e0, -2.55e0, -1.66e0, -0.43e0,  0.58e0/
      data  y2b( 31), y2b( 32), y2b( 33), y2b( 34), y2b( 35), y2b( 36) &
          / 2.18e0, -0.24e0,  0.58e0, -0.18e0, -1.55e0, -0.64e0/
      data  y2b( 37), y2b( 38), y2b( 39), y2b( 40), y2b( 41), y2b( 42) &
          /-1.09e0,  0.90e0, -0.66e0, -0.35e0,  0.48e0,  0.50e0/
      data  y2b( 43), y2b( 44), y2b( 45), y2b( 46), y2b( 47), y2b( 48) &
          / 0.05e0, -0.68e0,  0.24e0,  0.58e0, -1.26e0, -0.25e0/
      data  y2b( 49), y2b( 50), y2b( 51), y2b( 52), y2b( 53), y2b( 54) &
          / 0.25e0,  2.18e0,  2.96e0,  1.56e0, -0.36e0, -0.59e0/
      data  y2b( 55), y2b( 56), y2b( 57), y2b( 58), y2b( 59), y2b( 60) &
          /-0.12e0,  3.03e0,  2.11e0,  0.78e0,  0.89e0, -1.45e0/
      data  y2b( 61), y2b( 62), y2b( 63), y2b( 64), y2b( 65), y2b( 66) &
          /-0.36e0, -0.37e0, -1.39e0, -4.19e0, -0.73e0, -0.98e0/
      data  y2b( 67), y2b( 68), y2b( 69), y2b( 70), y2b( 71), y2b( 72) &
          / 0.36e0,  0.06e0, -1.94e0, -0.08e0,  0.17e0,  1.00e0/
      data  y2b( 73), y2b( 74), y2b( 75), y2b( 76), y2b( 77), y2b( 78) &
          /-0.05e0,  0.43e0,  0.15e0,  2.69e0,  0.57e0,  0.29e0/
      data  y2b( 79), y2b( 80), y2b( 81), y2b( 82), y2b( 83), y2b( 84) &
          / 1.10e0,  0.48e0, -1.06e0, -2.28e0, -2.03e0, -0.75e0/
      data  y2b( 85), y2b( 86), y2b( 87), y2b( 88), y2b( 89), y2b( 90) &
          / 1.00e0,  1.71e0,  0.58e0,  1.97e0,  0.99e0,  1.94e0/
      data  y2b( 91), y2b( 92), y2b( 93), y2b( 94), y2b( 95), y2b( 96) &
          / 2.18e0,  3.14e0,  0.60e0,  0.51e0,  1.35e0,  0.56e0/
      data  y2b( 97), y2b( 98), y2b( 99), y2b(100) &
          / 0.11e0,  0.00e0,  2.34e0,  1.88e0/
!
      data  y3(  1),y3(  2),y3(  3),y3(  4),y3(  5),y3(  6) &
          /-0.88e0, -0.12e0, -0.89e0, -1.38e0, -0.07e0,  1.03e0/
      data  y3(  7),y3(  8),y3(  9),y3( 10),y3( 11),y3( 12) &
          / 2.14e0,  0.35e0, -1.10e0, -1.78e0, -2.76e0, -1.77e0/
      data  y3( 13),y3( 14),y3( 15),y3( 16),y3( 17),y3( 18) &
          / 0.98e0,  1.00e0, -0.70e0, -1.01e0, -1.30e0, -0.85e0/
      data  y3( 19),y3( 20),y3( 21),y3( 22),y3( 23),y3( 24) &
          /-0.46e0,  1.63e0,  0.06e0, -0.17e0, -1.01e0, -1.04e0/
      data  y3( 25),y3( 26),y3( 27),y3( 28),y3( 29),y3( 30) &
          /-0.66e0, -1.12e0, -0.51e0, -0.71e0, -0.20e0, -0.13e0/
      data  y3( 31),y3( 32),y3( 33),y3( 34),y3( 35),y3( 36) &
          / 0.14e0,  1.59e0, -0.76e0, -1.08e0, -1.77e0, -1.20e0/
      data  y3( 37),y3( 38),y3( 39),y3( 40),y3( 41),y3( 42) &
          / 0.45e0, -0.07e0, -0.63e0, -0.35e0, -0.87e0, -0.62e0/
      data  y3( 43),y3( 44),y3( 45),y3( 46),y3( 47),y3( 48) &
          / 0.28e0,  1.90e0,  2.14e0,  1.05e0,  0.31e0,  1.07e0/
      data  y3( 49),y3( 50) &
          / 2.67e0,  2.44e0/
!
      data  y5(  1), y5(  2), y5(  3), y5(  4), y5(  5), y5(  6) &
          /     5.0e0, 11.0e0, 16.0e0, 23.0e0, 36.0e0, 58.0e0/
      data  y5(  7), y5(  8), y5(  9), y5( 10), y5( 11), y5( 12) &
          /    29.0e0, 20.0e0, 10.0e0,  8.0e0,  3.0e0,  0.0e0/
      data  y5( 13), y5( 14), y5( 15), y5( 16), y5( 17), y5( 18) &
          /     0.0e0, 2.0e0, 11.0e0, 27.0e0, 47.0e0, 63.0e0/
      data  y5( 19), y5( 20), y5( 21), y5( 22), y5( 23), y5( 24) &
          /    60.0e0, 39.0e0, 28.0e0, 26.0e0, 22.0e0, 11.0e0/
      data  y5( 25), y5( 26), y5( 27), y5( 28), y5( 29), y5( 30) &
          /    21.0e0, 40.0e0, 78.0e0,122.0e0,103.0e0, 73.0e0/
      data  y5( 31), y5( 32), y5( 33), y5( 34), y5( 35), y5( 36) &
          /    47.0e0, 35.0e0, 11.0e0,  5.0e0, 16.0e0, 34.0e0/
      data  y5( 37), y5( 38), y5( 39), y5( 40), y5( 41), y5( 42) &
          /    70.0e0, 81.0e0,111.0e0,101.0e0, 73.0e0, 40.0e0/
      data  y5( 43), y5( 44), y5( 45), y5( 46), y5( 47), y5( 48) &
          /    20.0e0, 16.0e0,  5.0e0, 11.0e0, 22.0e0, 40.0e0/
      data  y5( 49), y5( 50), y5( 51), y5( 52), y5( 53), y5( 54) &
          /    60.0e0, 80.9e0, 83.4e0, 47.7e0, 47.8e0, 30.7e0/
      data  y5( 55), y5( 56), y5( 57), y5( 58), y5( 59), y5( 60) &
          /    12.2e0,  9.6e0, 10.2e0, 32.4e0, 47.6e0, 54.0e0/
      data  y5( 61), y5( 62), y5( 63), y5( 64), y5( 65), y5( 66) &
          /    62.9e0, 85.9e0, 61.2e0, 45.1e0, 36.4e0, 20.9e0/
      data  y5( 67), y5( 68), y5( 69), y5( 70), y5( 71), y5( 72) &
          /    11.4e0, 37.8e0, 69.8e0,106.1e0,100.8e0, 81.6e0/
      data  y5( 73), y5( 74), y5( 75), y5( 76), y5( 77), y5( 78) &
          /    66.5e0, 34.8e0, 30.6e0,  7.0e0, 19.8e0, 92.5e0/
      data  y5( 79), y5( 80), y5( 81), y5( 82), y5( 83), y5( 84) &
          /   154.4e0,125.9e0, 84.8e0, 68.1e0, 38.5e0, 22.8e0/
      data  y5( 85), y5( 86), y5( 87), y5( 88), y5( 89), y5( 90) &
          /    10.2e0, 24.1e0, 82.9e0,132.0e0,130.9e0,118.1e0/
      data  y5( 91), y5( 92), y5( 93), y5( 94), y5( 95), y5( 96) &
          /    89.9e0, 66.6e0, 60.0e0, 46.9e0, 41.0e0, 21.3e0/
      data  y5( 97), y5( 98), y5( 99), y5(100), y5(101), y5(102) &
          /    16.0e0,  6.4e0,  4.1e0,  6.8e0, 14.5e0, 34.0e0/
      data  y5(103), y5(104), y5(105), y5(106), y5(107), y5(108) &
          /    45.0e0, 43.1e0, 47.5e0, 42.2e0, 28.1e0, 10.1e0/
      data  y5(109), y5(110), y5(111), y5(112), y5(113), y5(114) &
          /     8.1e0,  2.5e0,  0.0e0,  1.4e0,  5.0e0, 12.2e0/
      data  y5(115), y5(116), y5(117), y5(118), y5(119), y5(120) &
          /    13.9e0, 35.4e0, 45.8e0, 41.1e0, 30.1e0, 23.9e0/
      data  y5(121), y5(122), y5(123), y5(124), y5(125), y5(126) &
          /    15.6e0,  6.6e0,  4.0e0,  1.8e0,  8.5e0, 16.6e0/
      data  y5(127), y5(128), y5(129), y5(130), y5(131), y5(132) &
          /    36.3e0, 49.6e0, 64.2e0, 67.0e0, 70.9e0, 47.8e0/
      data  y5(133), y5(134), y5(135), y5(136), y5(137), y5(138) &
          /    27.5e0,  8.5e0, 13.2e0, 56.9e0,121.5e0,138.3e0/
      data  y5(139), y5(140), y5(141), y5(142), y5(143), y5(144) &
          /   103.2e0, 85.7e0, 64.6e0, 36.7e0, 24.2e0, 10.7e0/
      data  y5(145), y5(146), y5(147), y5(148), y5(149), y5(150) &
          /    15.0e0, 40.1e0, 61.5e0, 98.5e0,124.7e0, 96.3e0/
      data  y5(151), y5(152), y5(153), y5(154), y5(155), y5(156) &
          /    66.6e0, 64.5e0, 54.1e0, 39.0e0, 20.6e0,  6.7e0/
      data  y5(157), y5(158), y5(159), y5(160), y5(161), y5(162) &
          /     4.3e0, 22.7e0, 54.8e0, 93.8e0, 95.8e0, 77.2e0/
      data  y5(163), y5(164), y5(165), y5(166), y5(167), y5(168) &
          /    59.1e0, 44.0e0, 47.0e0, 30.5e0, 16.3e0,  7.3e0/
      data  y5(169), y5(170), y5(171), y5(172), y5(173), y5(174) &
          /    37.6e0, 74.0e0,139.0e0,111.2e0,101.6e0, 66.2e0/
      data  y5(175), y5(176), y5(177), y5(178), y5(179), y5(180) &
          /    44.7e0, 17.0e0, 11.3e0, 12.4e0,  3.4e0,  6.0e0/
      data  y5(181), y5(182), y5(183), y5(184), y5(185), y5(186) &
          /    32.3e0, 54.3e0, 59.7e0, 63.7e0, 63.5e0, 52.2e0/
      data  y5(187), y5(188), y5(189), y5(190), y5(191), y5(192) &
          /    25.4e0, 13.1e0,  6.8e0,  6.3e0,  7.1e0, 35.6e0/
      data  y5(193), y5(194), y5(195), y5(196), y5(197), y5(198) &
          /    73.0e0, 85.1e0, 78.0e0, 64.0e0, 41.8e0, 26.2e0/
      data  y5(199), y5(200), y5(201), y5(202), y5(203), y5(204) &
          /    26.7e0, 12.1e0,  9.5e0,  2.7e0,  5.0e0, 24.4e0/
      data  y5(205), y5(206), y5(207), y5(208), y5(209), y5(210) &
          /    42.0e0, 63.5e0, 53.8e0, 62.0e0, 48.5e0, 43.9e0/
      data  y5(211), y5(212), y5(213), y5(214), y5(215), y5(216) &
          /    18.6e0,  5.7e0,  3.6e0,  1.4e0,  9.6e0, 47.4e0/
      data  y5(217), y5(218), y5(219), y5(220), y5(221), y5(222) &
          /    57.1e0,103.9e0, 80.6e0, 63.6e0, 37.6e0, 26.1e0/
      data  y5(223), y5(224), y5(225), y5(226), y5(227), y5(228) &
          /    14.2e0,  5.8e0, 16.7e0, 44.3e0, 63.9e0, 69.0e0/
      data  y5(229), y5(230), y5(231), y5(232), y5(233), y5(234) &
          /    77.8e0, 64.9e0, 35.7e0, 21.2e0, 11.1e0,  5.7e0/
      data  y5(235), y5(236), y5(237), y5(238), y5(239), y5(240) &
          /     8.7e0, 36.1e0, 79.7e0,114.4e0,109.6e0, 88.8e0/
      data  y5(241), y5(242), y5(243), y5(244), y5(245), y5(246) &
          /    67.8e0, 47.5e0, 30.6e0, 16.3e0,  9.6e0, 33.2e0/
      data  y5(247), y5(248), y5(249), y5(250), y5(251), y5(252) &
          /    92.6e0,151.6e0,136.3e0,134.7e0, 83.9e0, 69.4e0/
      data  y5(253), y5(254), y5(255), y5(256), y5(257), y5(258) &
          /    31.5e0, 13.9e0,  4.4e0, 38.0e0,141.7e0,190.2e0/
      data  y5(259), y5(260), y5(261) &
          /   184.8e0,159.0e0,112.3e0/
!
!  define constants
!
      lper = 300
      lfreq = 300
      lyfft = 600

      ny1 = 100
      ny2 = 50
      ny3 = 50
      ny4 = 50
      ny5 = 261
      ny6 = 100

      nk = 3
      kmd(1) = 8
      kmd(2) = 8
      kmd(3) = 8

      taperp = 0.10

      nfft = 514
      iextnd = 0
      nprt = -1

      ilog = 1
!
!  write header
!
      write ( *, 1000)
!
!  Run simple test of acf
!
      write ( *, 1100)
      call acf ( y1, ny1 )
      write ( *,2000) ierr
!
!  Run simple test of ccf
!
      write ( *, 1200)
      call ccf ( y2a, y2b, ny2 )
      write ( *,2000) ierr
!
!  Run simple test of ufs
!
      write ( *, 1300)
      call ufs ( y3, ny3 )
      write ( *,2000) ierr
!
!  Run simple test of uas
!
      write ( *, 1400)
      call uas ( y4, ny4 )
      write ( *,2000) ierr
!
!  Run simple test of taper
!
      write ( *, 1510)
      call taper ( y5, ny5, taperp, yfft )
      write ( *,2000) ierr
!
!  Run simple test of pgms
!
      write ( *, 1520)
      call pgms ( yfft, ny5, nfft, lyfft, &
                iextnd, nf, per, lper, freq, lfreq, nprt )
      write ( *,2000) ierr
!
!  Run simple test of mdflt
!
      write ( *, 1530)
      call mdflt ( per, nf, nk, kmd, perf, ldstak )
      write ( *,2000) ierr
!
!  display results of mdflt
!
      write ( *, 1540)
      call ppl ( perf, freq, nf, ilog )
      write ( *,2000) ierr
!
!  Run simple test of bfs
!
      write ( *, 1600)
      call bfs ( y6a, y6b, ny6 )
      write ( *,2000) ierr

      return

 1000 format ('1*ch12')
 1100 format (' simple test of acf')
 1200 format ('1simple test of ccf')
 1300 format ('1simple test of ufs')
 1400 format ('1simple test of uas')
 1510 format ('1simple test of taper (no output unless error found)')
 1520 format ('1simple test of pgms')
 1530 format ('1simple test of mdflt (no output unless error found)')
 1540 format ('1display results of mdflt')
 1600 format ('1simple test of bfs')
 2000 format (/' the value of ierr is ', i4)
end
subroutine xxch13 ( ldstak )

!*****************************************************************************80
!
!! XXCH13 tests the arima modeling and forecasting family of routines.
!
!  Discussion:
!
!    Data is the airline data listed on page 531 of box and jenkins.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    George Box, Gwilym Jenkins,
!    Time Series Analysis: Forecasting and Control,
!    Holden-Day, 1970,
!    QA280.B67.
!
!  Parameters:
!
!     real air(200)
!        the airline data.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index variable.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr .ge. 1, errors were detected.
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer mspec(4,10)
!        the array containing the values of p, d, q, and s for each
!        factor.
!     integer n
!        the number of observations.
!     integer nfac
!        the number of factors in the model
!     integer npar
!        the number of unknown parameters in the model.
!     real par(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real res(200)
!        the residuals from the fit.
!     real y(200)
!        the array of the dependent variable.
!
  implicit none

      integer &
         ldstak
!
!  scalars in common
      integer &
         ierr
!
!  arrays in common
      double precision dstak(12)
!
!  local scalars
      integer n,nfac,npar
!
!  local arrays
      real &
         air(200),par(10),res(200),y(200)
      integer &
         mspec(4,10)
!
!  external subroutines
      external aime,aimf,iprint
!
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr
!
!  define constants
!
      data mspec(1,1), mspec(2,1), mspec(3,1), mspec(4,1) &
         /          0,          1,          1,          1/
      data mspec(1,2), mspec(2,2), mspec(3,2), mspec(4,2) &
         /          0,          1,          1,         12/

      data  air(  1), air(  2), air(  3), air(  4), air(  5), air(  6) &
          / 112.0e0, 118.0e0, 132.0e0, 129.0e0, 121.0e0, 135.0e0/
      data  air(  7), air(  8), air(  9), air( 10), air( 11), air( 12) &
          / 148.0e0, 148.0e0, 136.0e0, 119.0e0, 104.0e0, 118.0e0/
      data  air( 13), air( 14), air( 15), air( 16), air( 17), air( 18) &
          / 115.0e0, 126.0e0, 141.0e0, 135.0e0, 125.0e0, 149.0e0/
      data  air( 19), air( 20), air( 21), air( 22), air( 23), air( 24) &
          / 170.0e0, 170.0e0, 158.0e0, 133.0e0, 114.0e0, 140.0e0/
      data  air( 25), air( 26), air( 27), air( 28), air( 29), air( 30) &
          / 145.0e0, 150.0e0, 178.0e0, 163.0e0, 172.0e0, 178.0e0/
      data  air( 31), air( 32), air( 33), air( 34), air( 35), air( 36) &
          / 199.0e0, 199.0e0, 184.0e0, 162.0e0, 146.0e0, 166.0e0/
      data  air( 37), air( 38), air( 39), air( 40), air( 41), air( 42) &
          / 171.0e0, 180.0e0, 193.0e0, 181.0e0, 183.0e0, 218.0e0/
      data  air( 43), air( 44), air( 45), air( 46), air( 47), air( 48) &
          / 230.0e0, 242.0e0, 209.0e0, 191.0e0, 172.0e0, 194.0e0/
      data  air( 49), air( 50), air( 51), air( 52), air( 53), air( 54) &
          / 196.0e0, 196.0e0, 236.0e0, 235.0e0, 229.0e0, 243.0e0/
      data  air( 55), air( 56), air( 57), air( 58), air( 59), air( 60) &
          / 264.0e0, 272.0e0, 237.0e0, 211.0e0, 180.0e0, 201.0e0/
      data  air( 61), air( 62), air( 63), air( 64), air( 65), air( 66) &
          / 204.0e0, 188.0e0, 235.0e0, 227.0e0, 234.0e0, 264.0e0/
      data  air( 67), air( 68), air( 69), air( 70), air( 71), air( 72) &
          / 302.0e0, 293.0e0, 259.0e0, 229.0e0, 203.0e0, 229.0e0/
      data  air( 73), air( 74), air( 75), air( 76), air( 77), air( 78) &
          / 242.0e0, 233.0e0, 267.0e0, 269.0e0, 270.0e0, 315.0e0/
      data  air( 79), air( 80), air( 81), air( 82), air( 83), air( 84) &
          / 364.0e0, 347.0e0, 312.0e0, 274.0e0, 237.0e0, 278.0e0/
      data  air( 85), air( 86), air( 87), air( 88), air( 89), air( 90) &
          / 284.0e0, 277.0e0, 317.0e0, 313.0e0, 318.0e0, 374.0e0/
      data  air( 91), air( 92), air( 93), air( 94), air( 95), air( 96) &
          / 413.0e0, 405.0e0, 355.0e0, 306.0e0, 271.0e0, 306.0e0/
      data  air( 97), air( 98), air( 99), air(100), air(101), air(102) &
          / 315.0e0, 301.0e0, 356.0e0, 348.0e0, 355.0e0, 422.0e0/
      data  air(103), air(104), air(105), air(106), air(107), air(108) &
          / 465.0e0, 467.0e0, 404.0e0, 347.0e0, 305.0e0, 336.0e0/
      data  air(109), air(110), air(111), air(112), air(113), air(114) &
          / 340.0e0, 318.0e0, 362.0e0, 348.0e0, 363.0e0, 435.0e0/
      data  air(115), air(116), air(117), air(118), air(119), air(120) &
          / 491.0e0, 505.0e0, 404.0e0, 359.0e0, 310.0e0, 337.0e0/
      data  air(121), air(122), air(123), air(124), air(125), air(126) &
          / 360.0e0, 342.0e0, 406.0e0, 396.0e0, 420.0e0, 472.0e0/
      data  air(127), air(128), air(129), air(130), air(131), air(132) &
          / 548.0e0, 559.0e0, 463.0e0, 407.0e0, 362.0e0, 405.0e0/
      data  air(133), air(134), air(135), air(136), air(137), air(138) &
          / 417.0e0, 391.0e0, 419.0e0, 461.0e0, 472.0e0, 535.0e0/
      data  air(139), air(140), air(141), air(142), air(143), air(144) &
          / 622.0e0, 606.0e0, 508.0e0, 461.0e0, 390.0e0, 432.0e0/

      nfac = 2
      n = 144

      npar = 3
      par(1) = 0.000
      par(2) = 0.395
      par(3) = 0.615

      y(1:n) = log(air(1:n))
!
!  run simple test of aime
!
      write ( *,1000)
      write ( *,1100)
      call aime (y, n, mspec, nfac, par, npar, res, ldstak)
      write ( *,2000) ierr
!
!  run simple test of aimf
!
      write ( *,1200)
      call aimf (y, n, mspec, nfac, par, npar, ldstak)
      write ( *,2000) ierr

      return
 1000 format ('1*ch13')
 1100 format (' simple test of aime')
 1200 format ('1simple test of aimf')
 2000 format (/' the value of ierr is ', i4)
      end
