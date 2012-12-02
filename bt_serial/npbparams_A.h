! CLASS = A
!
!
!  This file is generated automatically by the setparams utility.
!  It sets the number of processors and the class of the NPB
!  in this directory. Do not modify it by hand.
!
        integer problem_size, niter_default
        parameter (problem_size=64, niter_default=200)
        double precision dt_default
        parameter (dt_default = 0.0008d0)
        logical  convertdouble
        parameter (convertdouble = .false.)
        character compiletime*11
        parameter (compiletime='28 Mar 2009')
        character npbversion*3
        parameter (npbversion='3.3')
        character cs1*23
        parameter (cs1='/usr/local/bin/gfortran')
        character cs2*6
        parameter (cs2='$(F90)')
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
