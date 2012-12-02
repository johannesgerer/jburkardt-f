! CLASS = W
!
!
!  This file is generated automatically by the setparams utility.
!  It sets the number of processors and the class of the NPB
!  in this directory. Do not modify it by hand.
!
        integer            na, nonzer, niter
        double precision   shift, rcond
        parameter(  na=7000, &
                    nonzer=8,&
                    niter=15,&
                    shift=12.,&
                    rcond=1.0d-1 )
        logical  convertdouble
        parameter (convertdouble = .false.)
        character compiletime*11
        parameter (compiletime='07 Apr 2009')
        character npbversion*3
        parameter (npbversion='3.3')
        character cs1*23
        parameter (cs1='/usr/local/bin/gfortran')
        character cs2*6
        parameter (cs2='$(F77)')
        character cs3*6
        parameter (cs3='(none)')
        character cs4*6
        parameter (cs4='(none)')
        character cs5*2
        parameter (cs5='-O')
        character cs6*2
        parameter (cs6='-O')
        character cs7*6
        parameter (cs7='randi8')
