subroutine assess ( d, iv, p, step, stlstg, v, x, x0 )

!*****************************************************************************80
!
!! ASSESS assesses a candidate step.
!
!  Discussion:
!
!    This subroutine is called by an unconstrained minimization
!    routine to assess the next candidate step.  It may recommend one
!    of several courses of action, such as accepting the step,
!    recomputing it using the same or a new quadratic model, or
!    halting due to convergence or false convergence.  See the return
!    code listing below.
!
!    This routine is called as part of the NL2SOL (nonlinear
!    least-squares) package.  It may be used in any unconstrained
!    minimization solver that uses dogleg, Goldfeld-Quandt-Trotter,
!    or Levenberg-Marquardt steps.
!
!    See Dennis, Gay and Welsch for further discussion of the assessing
!    and model switching strategies.  While NL2SOL considers only two
!    models, ASSESS is designed to handle any number of models.
!
!    On the first call of an iteration, only the I/O variables
!    step, X, IV(IRC), IV(MODEL), V(F), V(DSTNRM), V(GTSTEP), and
!    V(PREDUC) need have been initialized.  Between calls, no I/O
!    values execpt STEP, X, IV(MODEL), V(G) and the stopping tolerances
!    should be changed.
!
!    After a return for convergence or false convergence, one can
!    change the stopping tolerances and call ASSESS again, in which
!    case the stopping tests will be repeated.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Reference:
!
!    John Dennis, David Gay, Roy Welsch,
!    An Adaptive Nonlinear Least Squares Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, 1981.
!
!    M J D Powell,
!    A FORTRAN Subroutine for Solving Systems of Nonlinear Algebraic Equations,
!    in Numerical Methods for Nonlinear Algebraic Equations,
!    edited by Philip Rabinowitz,
!    Gordon and Breach, London, 1970.
!
!  Parameters:
!
!     iv (i/o) integer parameter and scratch vector -- see description
!             below of iv values referenced.
!
!    Input, real D(P), a scale vector used in computing V(RELDX).
!
!    Input, integer P, the number of parameters being optimized.
!
!   step (i/o) on input, step is the step to be assessed.  it is un-
!             changed on output unless a previous step achieved a
!             better objective function reduction, in which case stlstg
!             will have been copied to step.
!
! stlstg (i/o) when assess recommends recomputing step even though the
!             current (or a previous) step yields an objective func-
!             tion decrease, it saves in stlstg the step that gave the
!             best function reduction seen so far (in the current itera-
!             tion).  if the recomputed step yields a larger function
!             value, then step is restored from stlstg and
!             x = x0 + step is recomputed.
!
!      v (i/o) real parameter and scratch vector -- see description
!             below of v values referenced.
!
!      x (i/o) on input, x = x0 + step is the point at which the objec-
!             tive function has just been evaluated.  if an earlier
!             step yielded a bigger function decrease, then x is
!             restored to the corresponding earlier value.  otherwise,
!             if the current step does not give any function decrease,
!             then x is restored to x0.
!
!     x0 (in)  initial objective function parameter vector (at the
!             start of the current iteration).
!
!  iv values referenced
!
!    iv(irc) (i/o) on input for the first step tried in a new iteration,
!             iv(irc) should be set to 3 or 4 (the value to which it is
!             set when step is definitely to be accepted).  on input
!             after step has been recomputed, iv(irc) should be
!             unchanged since the previous return of assess.
!                on output, iv(irc) is a return code having one of the
!             following values...
!                  1 = switch models or try smaller step.
!                  2 = switch models or accept step.
!                  3 = accept step and determine v(radfac) by gradient
!                       tests.
!                  4 = accept step, v(radfac) has been determined.
!                  5 = recompute step (using the same model).
!                  6 = recompute step with radius = v(lmax0) but do not
!                       evaulate the objective function.
!                  7 = x-convergence (see v(xctol)).
!                  8 = relative function convergence (see v(rfctol)).
!                  9 = both x- and relative function convergence.
!                 10 = absolute function convergence (see v(afctol)).
!                 11 = singular convergence (see v(lmax0)).
!                 12 = false convergence (see v(xftol)).
!                 13 = iv(irc) was out of range on input.
!             return code i has precdence over i+1 for i = 9, 10, 11.
! iv(mlstgd) (i/o) saved value of iv(model).
!  iv(model) (i/o) on input, iv(model) should be an integer identifying
!             the current quadratic model of the objective function.
!             if a previous step yielded a better function reduction,
!             then iv(model) will be set to iv(mlstgd) on output.
! iv(nfcall) (in)  invocation count for the objective function.
! iv(nfgcal) (i/o) value of iv(nfcall) at step that gave the biggest
!             function reduction this iteration.  iv(nfgcal) remains
!             unchanged until a function reduction is obtained.
! iv(radinc) (i/o) the number of radius increases (or minus the number
!             of decreases) so far this iteration.
! iv(restor) (out) set to 0 unless x and v(f) have been restored, in
!             which case assess sets iv(restor) = 1.
!  iv(stage) (i/o) count of the number of models tried so far in the
!             current iteration.
! iv(stglim) (in)  maximum number of models to consider.
! iv(switch) (out) set to 0 unless a new model is being tried and it
!             gives a smaller function value than the previous model,
!             in which case assess sets iv(switch) = 1.
! iv(toobig) (in)  is nonzero if step was too big (e.g. if it caused
!             overflow).
!   iv(xirc) (i/o) value that iv(irc) would have in the absence of
!             convergence, false convergence, and oversized steps.
!
!  v values referenced
!
! v(afctol) (in)  absolute function convergence tolerance.  if the
!             absolute value of the current function value v(f) is less
!             than v(afctol), then assess returns with iv(irc) = 10.
! v(decfac) (in)  factor by which to decrease radius when iv(toobig) is
!             nonzero.
! v(dstnrm) (in)  the 2-norm of d * step.
! v(dstsav) (i/o) value of v(dstnrm) on saved step.
!   v(dst0) (in)  the 2-norm of d times the Newton step (when defined,
!             i.e., for v(nreduc) >= 0).
!      v(f) (i/o) on both input and output, v(f) is the objective func-
!             tion value at x.  if x is restored to a previous value,
!             then v(f) is restored to the corresponding value.
!   v(fdif) (out) the function reduction v(f0) - v(f) (for the output
!             value of v(f) if an earlier step gave a bigger function
!             decrease, and for the input value of v(f) otherwise).
! v(flstgd) (i/o) saved value of v(f).
!     v(f0) (in)  objective function value at start of iteration.
! v(gtslst) (i/o) value of v(gtstep) on saved step.
! v(gtstep) (in)  inner product between step and gradient.
! v(incfac) (in)  minimum factor by which to increase radius.
!  v(lmax0) (in)  maximum reasonable step size (and initial step bound).
!             if the actual function decrease is no more than twice
!             what was predicted, if a return with iv(irc) = 7, 8, 9,
!             or 10 does not occur, if v(dstnrm) > v(lmax0), and if
!             v(preduc) <= v(rfctol) * abs(v(f0)), then assess re-
!             turns with iv(irc) = 11.  if so doing appears worthwhile,
!             then assess repeats this test with v(preduc) computed for
!             a step of length v(lmax0) (by a return with iv(irc) = 6).
! v(nreduc) (i/o)  function reduction predicted by quadratic model for
!             Newton step.  if assess is called with iv(irc) = 6, i.e.,
!             if v(preduc) has been computed with radius = v(lmax0) for
!             use in the singular convervence test, then v(nreduc) is
!             set to -v(preduc) before the latter is restored.
! v(plstgd) (i/o) value of v(preduc) on saved step.
! v(preduc) (i/o) function reduction predicted by quadratic model for
!             current step.
! v(radfac) (out) factor to be used in determining the new radius,
!             which should be v(radfac)*dst, where  dst  is either the
!             output value of v(dstnrm) or the 2-norm of
!             diag(newd) * step  for the output value of step and the
!             updated version, newd, of the scale vector d.  for
!             iv(irc) = 3, v(radfac) = 1.0 is returned.
! v(rdfcmn) (in)  minimum value for v(radfac) in terms of the input
!             value of v(dstnrm) -- suggested value = 0.1.
! v(rdfcmx) (in)  maximum value for v(radfac) -- suggested value = 4.0.
!  v(reldx) (out) scaled relative change in x caused by step, computed
!             by function  reldst  as
!                 max (d(i)*abs(x(i)-x0(i)), 1 <= i <= p) /
!                    max (d(i)*(abs(x(i))+abs(x0(i))), 1 <= i <= p).
!             if an acceptable step is returned, then v(reldx) is com-
!             puted using the output (possibly restored) values of x
!             and step.  otherwise it is computed using the input
!             values.
! v(rfctol) (in)  relative function convergence tolerance.  if the
!             actual function reduction is at most twice what was pre-
!             dicted and  v(nreduc) <= v(rfctol)*abs(v(f0)),  then
!             assess returns with iv(irc) = 8 or 9.  see also v(lmax0).
! v(STPPAR) (in)  Marquardt parameter -- 0 means full Newton step.
! v(tuner1) (in)  tuning constant used to decide if the function
!             reduction was much less than expected.  suggested
!             value = 0.1.
! v(tuner2) (in)  tuning constant used to decide if the function
!             reduction was large enough to accept step.  suggested
!             value = 10**-4.
! v(tuner3) (in)  tuning constant used to decide if the radius
!             should be increased.  suggested value = 0.75.
!  v(xctol) (in)  x-convergence criterion.  if step is a Newton step
!             (v(STPPAR) = 0) having v(reldx) <= v(xctol) and giving
!             at most twice the predicted function decrease, then
!             assess returns iv(irc) = 7 or 9.
!  v(xftol) (in)  false convergence tolerance.  if step gave no or only
!             a small function decrease and v(reldx) <= v(xftol),
!             then assess returns with iv(irc) = 12.
!
  implicit none

  integer p

  integer, parameter :: afctol = 31
  real d(p)
  integer, parameter :: decfac = 22
  real emax
  logical goodx
  real gts
  integer i
  integer, parameter :: irc = 3
  integer iv(13)
  integer, parameter :: lmax0 = 35
  integer nfc
  integer, parameter :: nreduc = 6
  integer, parameter :: plstgd = 15
  integer, parameter :: preduc = 7
  integer, parameter :: radfac = 16
  integer, parameter :: rdfcmn = 24
  integer, parameter :: rdfcmx = 25
  real reldst
  integer, parameter :: reldx = 17
  real reldx1
  real rfac1
  integer, parameter :: rfctol = 32
  real step(p)
  real stlstg(p)
  integer, parameter :: stppar = 5
  integer, parameter :: tuner1 = 26
  integer, parameter :: tuner2 = 27
  integer, parameter :: tuner3 = 28
  real v(35)
  real x(p)
  real x0(p)
  integer, parameter :: xctol = 33
  integer, parameter :: xftol = 34
  integer, parameter :: xirc = 13
  real xmax
!
!  subscripts for iv and v
!
      integer dstnrm, dstsav, dst0, f, fdif, flstgd, f0
      integer gtslst, gtstep, incfac, mlstgd, model, nfcall, &
              nfgcal, radinc, &
              restor, stage, stglim, &
              switch, toobig
     parameter ( mlstgd=4, model=5, nfcall=6 )
     parameter ( nfgcal=7, radinc=8, restor=9, stage=10 )
     parameter ( stglim=11, switch=12, toobig=2 )
     parameter ( dstnrm=2, dst0=3 )
     parameter ( dstsav=18, f=10, fdif=11, flstgd=12, f0=13 )
     parameter ( gtslst=14, gtstep=4, incfac=23 )

  nfc = iv(nfcall)
  iv(switch) = 0
  iv(restor) = 0
  rfac1 = 1.0E+00
  goodx = .true.
  i = iv(irc)

  if ( i < 1 .or. 12 < i ) then
    iv(irc) = 13
    return
  end if

  go to (20,30,10,10,40,360,290,290,290,290,290,140), i
!
!  Initialize for new iteration.
!
 10   continue

      iv(stage) = 1
      iv(radinc) = 0
      v(flstgd) = v(f0)

      if ( iv(toobig) /= 0 ) then
         iv(stage) = -1
         iv(xirc) = i
         v(radfac) = v(decfac)
         iv(radinc) = iv(radinc) - 1
         iv(irc) = 5
         return
       end if

       go to 90
!
!  Step was recomputed with new model or smaller radius.
!  First decide which.
!
 20   continue
!
!  Old model retained, smaller radius tried.
!  Do not consider any more new models this iteration.
!
      if (iv(model) == iv(mlstgd)) then
         iv(stage) = iv(stglim)
         iv(radinc) = -1
         go to 90
      end if
!
!  A new model is being tried.  Decide whether to keep it.
!
 30   iv(stage) = iv(stage) + 1
!
!  Now we add the possibiltiy that step was recomputed with
!  the same model, perhaps because of an oversized step.
!
 40   if ( 0 < iv(stage) ) go to 50
!
!  Step was recomputed because it was too big.
!
         if (iv(toobig) /= 0) then
           v(radfac) = v(decfac)
           iv(radinc) = iv(radinc) - 1
           iv(irc) = 5
           return
         end if
!
!  Restore IV(STAGE) and pick up where we left off.
!
         iv(stage) = -iv(stage)
         i = iv(xirc)
         go to (20, 30, 90, 90, 70), i

 50   if (iv(toobig) == 0) go to 70
!
!  Handle oversize step.
!
      if ( iv(radinc) <= 0 ) then
         iv(stage) = -iv(stage)
         iv(xirc) = iv(irc)
         v(radfac) = v(decfac)
         iv(radinc) = iv(radinc) - 1
         iv(irc) = 5
         return
      end if

      go to 80

 70   if (v(f) < v(flstgd)) go to 90
!
!  The new step is a loser.  Restore old model.
!
      if ( iv(model) /= iv(mlstgd) ) then
        iv(model) = iv(mlstgd)
        iv(switch) = 1
      end if
!
!  Restore step, etc. only if a previous step decreased V(F).
!
 80   continue

      if ( v(flstgd) < v(f0) ) then
        iv(restor) = 1
        v(f) = v(flstgd)
        v(preduc) = v(plstgd)
        v(gtstep) = v(gtslst)
        if (iv(switch) == 0) rfac1 = v(dstnrm) / v(dstsav)
        v(dstnrm) = v(dstsav)
        nfc = iv(nfgcal)
        goodx = .false.
      end if
!
!  Compute relative change in X by current step.
!
 90   continue

      reldx1 = reldst ( p, d, x, x0 )
!
!  Restore X and STEP if necessary.
!
      if ( .not. goodx ) then

        step(1:p) = stlstg(1:p)
        x(1:p) = x0(1:p) + stlstg(1:p)

      end if

      v(fdif) = v(f0) - v(f)
!
!  No (or only a trivial) function decrease,
!  so try new model or smaller radius.
!
      if ( v(fdif) <= v(tuner2) * v(preduc) ) then

         v(reldx) = reldx1

         if ( v(f0) <= v(f) ) then
           iv(mlstgd) = iv(model)
           v(flstgd) = v(f)
           v(f) = v(f0)
           x(1:p) = x0(1:p)
           iv(restor) = 1
         else
           iv(nfgcal) = nfc
         end if

         iv(irc) = 1

         if ( iv(stglim) <= iv(stage) ) then
           iv(irc) = 5
           iv(radinc) = iv(radinc) - 1
         end if

      else
!
!  Nontrivial function decrease achieved.
!
      iv(nfgcal) = nfc
      rfac1 = 1.0E+00
      if ( goodx ) then
        v(reldx) = reldx1
      end if
      v(dstsav) = v(dstnrm)

      if ( v(preduc) * v(tuner1) < v(fdif) ) then
        go to 200
      end if
!
!  Decrease was much less than predicted: either change models
!  or accept step with decreased radius.
!
      if ( iv(stage) < iv(stglim) ) then
        iv(irc) = 2
      else
        iv(irc) = 4
      end if

    end if
!
!  Set V(RADFAC) to Fletcher's decrease factor.
!
      iv(xirc) = iv(irc)
      emax = v(gtstep) + v(fdif)
      v(radfac) = 0.5E+00 * rfac1

      if (emax < v(gtstep) ) then
        v(radfac) = rfac1 * max ( v(rdfcmn), 0.5E+00 * v(gtstep) / emax )
      end if
!
!  Do a false convergence test.
!
 140  continue

      if ( v(reldx) <= v(xftol)) go to 160
         iv(irc) = iv(xirc)
         if (v(f) < v(f0)) go to 230
              go to 300

 160  iv(irc) = 12
      go to 310
!
!  Handle good function decrease,
!
 200  if (v(fdif) < (-v(tuner3) * v(gtstep))) go to 260
!
!  Increasing radius looks worthwhile.  See if we just
!  recomputed step with a decreased radius or restored step
!  after recomputing it with a larger radius.
!
      if (iv(radinc) < 0) go to 260
      if (iv(restor) == 1) go to 260
!
!  We did not.  Try a longer step unless this was a Newton step.
!
         v(radfac) = v(rdfcmx)
         gts = v(gtstep)
         if ( v(fdif) < ( 0.5E+00 / v(radfac) - 1.0E+00) * gts ) then
           v(radfac) = max ( v(incfac), 0.5E+00 * gts / ( gts + v(fdif) ) )
         end if
         iv(irc) = 4

         if ( v(stppar) == 0.0E+00 ) then
           go to 300
         end if
!
!  Step was not a Newton step.  Recompute it with a larger radius.
!
              iv(irc) = 5
              iv(radinc) = iv(radinc) + 1
!
!  Save values corresponding to good step.
!
 230  v(flstgd) = v(f)
      iv(mlstgd) = iv(model)
      stlstg(1:p) = step(1:p)
      v(dstsav) = v(dstnrm)
      iv(nfgcal) = nfc
      v(plstgd) = v(preduc)
      v(gtslst) = v(gtstep)
      go to 300
!
!  Accept step with radius unchanged.
!
 260  v(radfac) = 1.0E+00
      iv(irc) = 3
      go to 300
!
!  Come here for a restart after convergence.
!
 290  iv(irc) = iv(xirc)

      if ( v(dstsav) < 0.0E+00 ) then
        iv(irc) = 12
      end if
      go to 310
!
!  Perform convergence tests.
!
 300  iv(xirc) = iv(irc)

 310  continue

      if ( abs ( v(f) ) < v(afctol) ) then
        iv(irc) = 10
      end if

      if ( 0.5E+00 * v(fdif) > v(preduc)) then
        return
      end if

      emax = v(rfctol) * abs(v(f0))

      if ( v(dstnrm) > v(lmax0) .and. v(preduc) <= emax) then
        iv(irc) = 11
      end if

      if ( 0.0E+00 <= v(dst0) ) then

        if ((v(nreduc) > 0.0E+00 .and. v(nreduc) <= emax) .or. &
            (v(nreduc) == 0.0E+00 .and. v(preduc) == 0.0))  then
          i = 2
        else
          i = 0
        end if

        if ( v(stppar) == 0.0E+00 .and. &
             v(reldx) <= v(xctol) .and. &
             goodx ) then
          i = i + 1
        end if

        if ( 0 < i ) then
          iv(irc) = i + 6
        end if

      end if
!
!  Consider recomputing step of length V(LMAX0) for singular
!  convergence test.
!
      if ( abs (iv(irc)-3) > 2 .and. iv(irc) /= 12 ) then
        return
      end if

      if ( v(lmax0) < v(dstnrm) ) then

        if ( 0.5E+00 * v(dstnrm) <= v(lmax0) ) then
          return
        end if

        xmax = v(lmax0) / v(dstnrm)

        if ( emax <= xmax * ( 2.0E+00 - xmax) * v(preduc) ) then
          return
        end if

      else

         if ( emax <= v(preduc) ) then
           return
         end if

         if ( 0.0E+00 < v(dst0) ) then

           if ( 0.5E+00 * v(dst0) <= v(lmax0) ) then
             return
           end if

         end if

      end if

  if ( v(nreduc) < 0.0E+00 ) then
    if ( -v(nreduc) <= v(rfctol) * abs ( v(f0) ) ) then
      iv(irc) = 11
    end if
    return
  end if
!
!  Recompute V(PREDUC) for use in singular convergence test.
!
  v(gtslst) = v(gtstep)
  v(dstsav) = v(dstnrm)
  if ( iv(irc) == 12 ) then
    v(dstsav) = -v(dstsav)
  end if
  v(plstgd) = v(preduc)
  iv(irc) = 6
  stlstg(1:p) = step(1:p)
  return
!
!  Perform singular convergence test with recomputed V(PREDUC).
!
 360  continue

  v(gtstep) = v(gtslst)
  v(dstnrm) = abs(v(dstsav))
  step(1:p) = stlstg(1:p)

  if ( v(dstsav) <= 0.0E+00 ) then
    iv(irc) = 12
  else
    iv(irc) = iv(xirc)
  end if

  v(nreduc) = -v(preduc)
  v(preduc) = v(plstgd)

  if ( -v(nreduc) <= v(rfctol) * abs(v(f0)) ) then
    iv(irc) = 11
  end if

  return
end
subroutine covclc ( covirc, d, iv, j, n, nn, p, r, v, x )

!*****************************************************************************80
!
!! COVCLC computes the covariance matrix for NL2ITR.
!
!  Discussion:
!
!    Let K = abs ( IV(COVREQ) ).
!
!    For K <= 2, a finite-difference hessian H is computed,
!    * using function and gradient values if IV(COVREQ) is nonnegative,
!    * using only function values if IV(COVREQ) is negative).
!
!    Let
!      SCALE = 2 * F(X) / max ( 1, N - P ),
!    where 2 * F(X) is the residual sum of squares.
!
!    COVCLC computes:
!      K = 0 or 1:  SCALE * inverse ( H ) * ( J' * J ) * inverse ( H ).
!      K = 2:       SCALE * inverse ( H );
!      K >= 3:      SCALE * inverse ( J' * J ).
!
!  Modified:
!
!    13 April 2006
!
!  Parameters:
!
!    ?, integer COVIRC, ?
!
!    Input, real D(P), the scaling vector.
!
!    Input/output, integer IV(*), the NL2SOL integer parameter vector.
!
!    Input, real J(NN,P), the N by P Jacobian matrix.
!
!    Input, integer N, the number of functions.
!
!    Input, integer NN, the leading dimension of J.
!
!    Input, integer P, the number of variables.
!
!    ?, real R(N), ?
!
!    Input, real V(*), the NL2SOL real parameter array.
!
!    ?, real X(P), ?
!
  implicit none

  integer n
  integer nn
  integer p

  integer cov
  integer covirc
  integer, parameter :: covmat = 26
  integer, parameter :: covreq = 15
  real d(p)
  real del
  integer, parameter :: delta = 50
  integer, parameter :: delta0 = 44
  integer, parameter :: dltfdc = 40
  integer, parameter :: f = 10
  integer, parameter :: fx = 46
  integer, parameter :: g = 28
  integer g1
  integer gp
  integer gsave1
  integer, parameter :: h = 44
  logical havej
  integer hc
  integer hmi
  integer hpi
  integer hpm
  integer i
  integer, parameter :: ierr = 32
  integer ip1
  integer, parameter :: ipiv0 = 60
  integer ipivi
  integer ipivk
  integer, parameter :: ipivot = 61
  integer irc
  integer iv(*)
  real j(nn,p)
  integer k
  integer, parameter :: kagqt = 35
  integer, parameter :: kalm = 36
  integer kind
  integer kl
  integer l
  integer, parameter :: lmat = 58
  integer m
  integer mm1
  integer mm1o2
  integer, parameter :: mode = 38
  integer, parameter :: nfgcal = 7
  integer pp1o2
  integer, parameter :: qtr = 49
  integer qtr1
  real r(n)
  integer, parameter :: rd = 51
  integer rd1
  integer, parameter :: rsave = 52
  integer, parameter :: savei = 54
  integer stp0
  integer stpi
  integer stpm
  integer, parameter :: switch = 12
  real t
  integer, parameter :: toobig = 2
  real v(*)
  integer, parameter :: w = 59
  integer w0
  integer w1
  real wk
  integer wl
  real x(p)
  integer, parameter :: xmsave = 49

  cov = iv(lmat)
  covirc = 4
  kind = iv(covreq)
  m = iv(mode)

  if ( m <= 0 ) then

    iv(kagqt) = -1

    if ( 0 < iv(kalm) ) then
      iv(kalm) = 0
    end if

    if ( 3 <= abs ( kind ) ) then

      rd1 = iv(rd)

      if ( iv(kalm) == -1 ) then
        qtr1 = iv(qtr)
        v(qtr1:qtr1+n-1) = r(1:n)
        w1 = iv(w) + p
        call qrfact ( nn, n, p, j, v(rd1), iv(ipivot), iv(ierr), 0, v(w1) )
        iv(kalm) = -2
      end if

      iv(covmat) = -1

      if (iv(ierr) /= 0) then
        return
      end if

      cov = iv(lmat)
      hc = abs ( iv(h) )
      iv(h) = -hc
!
!  Set HC = R matrix from QRFACT.
!
      l = hc
      do i = 1, p
        if ( 1 < i ) then
          call vcopy ( i-1, v(l), j(1,i) )
        end if
        l = l + i - 1
        v(l) = v(rd1)
        l = l + 1
        rd1 = rd1 + 1
      end do

      go to 350

    end if

    v(fx) = v(f)
    k = iv(rsave)
    v(k:k+n-1) = r(1:n)

  end if

  if ( m <= p ) then

    if ( kind < 0 ) go to 100
!
!  Compute finite-difference hessian using both function and
!  gradient values.
!
    gsave1 = iv(w) + p
    g1 = iv(g)
!
!  First call on COVCLC.  Set GSAVE = G, take first step.
!
    if ( m <= 0 ) then

      v(gsave1:gsave1+p-1) = v(g1:g1+p-1)
      iv(switch) = iv(nfgcal)

    else

      del = v(delta)
      x(m) = v(xmsave)
!
!  Handle oversize V(DELTA).
!
      if ( iv(toobig) /= 0 ) then

        if ( 0.0E+00 < del * x(m) ) then
          del = -0.5E+00 * del
          x(m) = x(m) + del
          v(delta) = del
          covirc = 2
          return
        end if

        iv(covmat) = -2
        go to 190

      end if

      cov = iv(lmat)
      gp = g1 + p - 1
!
!  Set G = ( G - GSAVE ) / DEL.
!
      do i = g1, gp
        v(i) = (v(i) - v(gsave1)) / del
        gsave1 = gsave1 + 1
      end do
!
!  Add G as new column to finite-difference hessian matrix.
!
      k = cov + ( m * ( m - 1 ) ) / 2
      l = k + m - 2
!
!  Set H(1:M-1,M) = 0.5 * (H(1:M-1,m) + G(1:M-1)).
!
      if ( m /= 1 ) then

        do i = k, l
          v(i) = 0.5E+00 * ( v(i) + v(g1) )
          g1 = g1 + 1
        end do

      end if
!
!  Add H(M:P,M) = G(M:P).
!
      l = l + 1
      do i = m, p
        v(l) = v(g1)
        l = l + i
        g1 = g1 + 1
      end do

    end if

    m = m + 1
    iv(mode) = m

    if ( p < m ) then
      go to 190
    end if
!
!  Choose next finite-difference step, return to get G there.
!
    del = v(delta0) * max ( 1.0E+00 / d(m), abs ( x(m) ) )
    if ( x(m) < 0.0E+00 ) then
      del = -del
    end if
    v(xmsave) = x(m)
    x(m) = x(m) + del
    v(delta) = del
    covirc = 2
    return
!
!  Compute finite-difference hessian using function values only.
!
 100  continue

    stp0 = iv(w) + p - 1
    mm1 = m - 1
    mm1o2 = m * mm1 / 2
!
!  First call on COVCLC.
!
    if ( m <= 0 ) then

      iv(savei) = 0

    else

      i = iv(savei)

      if ( i <= 0 ) then
!
!  Handle oversize step.
!
        if ( iv(toobig) /= 0 ) then

          stpm = stp0 + m
          del = v(stpm)
!
!  We already tried shrinking the step, so quit.
!
          if ( del * v(xmsave) <= 0.0E+00 ) then
            iv(covmat) = -2
            return
          end if
!
!  Try shrinking the step.
!
          del = -0.5E+00 * del
          x(m) = v(xmsave) + del
          v(stpm) = del
          covirc = 1
          return

        end if
!
!  Save F(X + STP(M)*E(M)) in H(P,M).
!
        pp1o2 = ( p * ( p - 1 ) ) / 2
        cov = iv(lmat)
        hpm = cov + pp1o2 + mm1
        v(hpm) = v(f)
!
!  Start computing row M of the finite-difference hessian H.
!
        hmi = cov + mm1o2
        hpi = cov + pp1o2

        do i = 1, mm1
          v(hmi) = v(fx) - (v(f) + v(hpi))
          hmi = hmi + 1
          hpi = hpi + 1
        end do

        v(hmi) = v(f) - 2.0E+00 * v(fx)
!
!  Compute function values needed to complete row M of H.
!
        i = 1

        iv(savei) = i
        stpi = stp0 + i
        v(delta) = x(i)
        x(i) = x(i) + v(stpi)
        if ( i == m ) then
          x(i) = v(xmsave) - v(stpi)
        end if
        covirc = 1
        return

      end if

      x(i) = v(delta)
!
!  Punt in the event of an oversize step.
!
      if ( iv(toobig) /= 0 ) then
        iv(covmat) = -2
        return
      end if
!
!  Finish computing H(M,I).
!
      stpi = stp0 + i
      hmi = cov + mm1o2 + i - 1
      stpm = stp0 + m
      v(hmi) = ( v(hmi) + v(f) ) / ( v(stpi) * v(stpm) )
      i = i + 1

      if ( i <= m ) then
        iv(savei) = i
        stpi = stp0 + i
        v(delta) = x(i)
        x(i) = x(i) + v(stpi)
        if ( i == m ) then
          x(i) = v(xmsave) - v(stpi)
        end if
        covirc = 1
        return
      end if

      iv(savei) = 0
      x(m) = v(xmsave)

    end if

    m = m + 1
    iv(mode) = m
!
!  Prepare to compute row M of the finite-difference hessian H.
!  Compute the M-th step size STP(M), then return to obtain
!  F(X + STP(M)*E(M)), where E(M) = M-th standard unit vector.
!
    if ( m <= p ) then
      del = v(dltfdc) * max ( 1.0E+00 / d(m), abs(x(m)) )
      if (x(m) < 0.0E+00 ) then
        del = -del
      end if
      v(xmsave) = x(m)
      x(m) = x(m) + del
      stpm = stp0 + m
      v(stpm) = del
      covirc = 1
      return
    end if
!
!  Restore R, V(F), etc.
!
 190  continue

    k = iv(rsave)
    r(1:n) = v(k:k+n-1)
    v(f) = v(fx)

    if ( 0 <= kind ) then

      iv(nfgcal) = iv(switch)
      qtr1 = iv(qtr)
      v(qtr1:qtr1+n-1) = r(1:n)

      if ( 0 <= iv(covmat) ) then
        covirc = 3
      end if

      return
    end if

  end if

  cov = iv(lmat)
!
!  The complete finite-difference hessian is now stored at V(COV).
!  Use it to compute the requested covariance matrix.
!
!  Compute Cholesky factor C of H = C * C' and store it at V(HC).
!
  hc = cov

  if ( abs ( kind ) /= 2 ) then
    hc = abs ( iv(h) )
    iv(h) = -hc
  end if

  call lsqrt ( 1, p, v(hc), v(cov), irc )
  iv(covmat) = -1

  if ( irc /= 0 ) then
    return
  end if

  w1 = iv(w) + p

  if ( 1 < abs ( kind ) ) then
    go to 350
  end if
!
!  Covariance = SCALE * inverse ( H ) * (J' * J) * inverse ( H ).
!
  v(cov:cov+(p*(p+1))/2) = 0.0E+00
  havej = iv(kalm) == (-1)
!
!  HAVEJ = .true. means J is in its original form, while
!  HAVEJ = .false. means QRFACT has been applied to J.
!
  if ( havej ) then
    m = n
  else
    m = p
  end if
  w0 = w1 - 1
  rd1 = iv(rd)

  do i = 1, m
!
!  Set W = IPIVOT * (row I of R matrix from QRFACT).
!
    if ( .not. havej ) then

      v(w1:w1+p-1) = 0.0E+00
      ipivi = ipiv0 + i
      l = w0 + iv(ipivi)
      v(l) = v(rd1)
      rd1 = rd1 + 1

      do k = i+1, p
        ipivk = ipiv0 + k
        l = w0 + iv(ipivk)
        v(l) = j(i,k)
      end do
!
!  Set W = (row I of J).
!
    else

      l = w0
      do k = 1, p
        l = l + 1
        v(l) = j(i,k)
      end do

    end if
!
!  Set W = inverse ( H ) * W.
!
    call livmul ( p, v(w1), v(hc), v(w1) )
    call litvmu ( p, v(w1), v(hc), v(w1) )
!
!  Add W * W' to covariance matrix.
!
    kl = cov
    do k = 1, p
      l = w0 + k
      wk = v(l)
      do l = 1, k
        wl = w0 + l
        v(kl) = v(kl)  +  wk * v(wl)
        kl = kl + 1
      end do
    end do

  end do

  go to 380
!
!  The Cholesky factor C of the unscaled inverse covariance matrix
!  (or permutation thereof) is stored at V(HC).
!
!  Set C = inverse ( C ).
!
 350  continue

  call linvrt ( p, v(hc), v(hc) )
!
!  Set C = C' * C.
!
  call ltsqar ( p, v(hc), v(hc) )
!
!  C = permuted, unscaled covariance.
!  Set COV = IPIVOT * C * IPIVOT'.
!
  if ( hc /= cov ) then

    do i = 1, p
      m = ipiv0 + i
      ipivi = iv(m)
      kl = cov-1 + ( ipivi * (ipivi-1) ) / 2
      do k = 1, i
        m = ipiv0 + k
        ipivk = iv(m)
        l = kl + ipivk
        if ( ipivi < ipivk ) then
          l = l + ( (ipivk-ipivi) * (ipivk+ipivi-3) ) / 2
        end if
        v(l) = v(hc)
        hc = hc + 1
      end do
    end do

  end if

380  continue

  iv(covmat) = cov
!
!  Apply scale factor = (residual sum of squares) / max(1,n-p).
!
  t = v(f) / ( 0.5E+00 * real ( max ( 1, n-p ) ) )
  k = cov - 1 + ( p * ( p + 1 ) ) / 2

  v(cov:k) = t * v(cov:k)

  return
end
subroutine dfault ( iv, v )

!*****************************************************************************80
!
!! DFAULT supplies default values to IV and V.
!
!  Discussion:
!
!    Only entries in the first 25 positions of IV and the first 45
!    positions of V are reset.
!
!  Modified:
!
!    05 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Output, integer IV(25), contains default values for specific entries.
!
!    Output, real V(45), contains default values for specific values.
!
  implicit none

  integer :: afctol = 31
  integer :: cosmin = 43
  integer :: covprt = 14
  integer :: covreq = 15
  integer :: d0init = 37
  integer :: decfac = 22
  integer :: delta0 = 44
  integer :: dfac = 41
  integer :: dinit = 38
  integer :: dltfdc = 40
  integer :: dltfdj = 36
  integer :: dtype = 16
  integer :: inits = 25
  integer :: epslon = 19
  integer :: fuzz = 45
  integer :: incfac = 23
  integer iv(25)
  integer :: jtinit = 39
  integer :: lmax0 = 35
  real machep
  real mepcrt
  integer :: mxfcal = 17
  integer :: mxiter = 18
  integer :: outlev = 19
  integer :: parprt = 20
  integer :: phmnfc = 20
  integer :: phmxfc = 21
  integer :: prunit = 21
  integer :: rdfcmn = 24
  integer :: rdfcmx = 25
  integer :: rfctol = 32
  integer :: rlimit = 42
  integer :: solprt = 22
  real sqteps
  integer :: statpr = 23
  integer :: tuner1 = 26
  integer :: tuner2 = 27
  integer :: tuner3 = 28
  integer :: tuner4 = 29
  integer :: tuner5 = 30
  real v(45)
  integer :: x0prt = 24
  integer :: xctol = 33
  integer :: xftol = 34

  iv(1) = 12
  iv(covprt) = 1
  iv(covreq) = 1
  iv(dtype) = 1
  iv(inits) = 0
  iv(mxfcal) = 200
  iv(mxiter) = 150
  iv(outlev) = -1
  iv(parprt) = 1
  iv(prunit) = 6
  iv(solprt) = 1
  iv(statpr) = 1
  iv(x0prt) = 1

  machep = epsilon ( machep )
  v(afctol) = 1.0e-20
  if ( 1.0e-10 < machep ) then
    v(afctol) = machep**2
  end if
  v(cosmin) = max ( 1.0e-06, 1.0e+02 * machep )
  v(decfac) = 0.5e+00
  sqteps = sqrt ( epsilon ( sqteps ) )
  v(delta0) = sqteps
  v(dfac) = 0.6e+00
  v(dinit) = 0.0e+00
  mepcrt = machep ** ( 1.0E+00 / 3.0E+00 )
  v(dltfdc) = mepcrt
  v(dltfdj) = sqteps
  v(d0init) = 1.0e+00
  v(epslon) = 0.1e+00
  v(fuzz) = 1.5e+00
  v(incfac) = 2.0e+00
  v(jtinit) = 1.0e-06
  v(lmax0) = 100.e+00
  v(phmnfc) = -0.1e+00
  v(phmxfc) = 0.1e+00
  v(rdfcmn) = 0.1e+00
  v(rdfcmx) = 4.e+00
  v(rfctol) = max ( 1.0E-10, mepcrt**2 )
  v(rlimit) = sqrt ( 0.999E+00 * huge ( v(rlimit) ) )
  v(tuner1) = 0.1e+00
  v(tuner2) = 1.0e-04
  v(tuner3) = 0.75e+00
  v(tuner4) = 0.5e+00
  v(tuner5) = 0.75e+00
  v(xctol) = sqteps
  v(xftol) = 1.0e+02 * machep

  return
end
function dotprd ( p, x, y )

!*****************************************************************************80
!
!! DOTPRD returns the inner product of two vectors.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input, integer P, the number of entries in the vectors.
!
!    Input, real X(P), Y(P), the vectors.
!
!    Output, real DOTPRD, the dot product of X and Y.
!
  implicit none

  integer p

  real dotprd
  integer i
  real, save :: sqteta = 0.0E+00
  real t
  real x(p)
  real y(p)

  dotprd = 0.0E+00

  if ( p <= 0 ) then
    return
  end if

  if ( sqteta == 0.0E+00 ) then
    sqteta = sqrt ( 1.001E+00 * tiny ( sqteta ) )
  end if

  do i = 1, p

    t = max ( abs ( x(i) ), abs ( y(i) ) )

    if ( t < sqteta ) then

    else if ( 1.0E+00 < t ) then

      dotprd = dotprd + x(i) * y(i)

    else

      t = ( x(i) / sqteta ) * y(i)

      if ( sqteta <= abs ( t ) ) then
        dotprd = dotprd + x(i) * y(i)
      end if

    end if

  end do

  return
end
subroutine dupdat ( d, iv, j, n, nn, p, v )

!*****************************************************************************80
!
!! DUPDAT updates the scale vector for NL2ITR.
!
!  Modified:
!
!    05 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input/output, real D(P), the scale vector.
!
!    Input, integer IV(*), the NL2SOL integer array.
!
!    Input, real J(NN,P), the N by P Jacobian matrix.
!
!    Input, integer N, the number of functions.
!
!    Input, integer NN, the leading dimension of J.
!
!    Input, integer P, the number of variables.
!
!    Input, real V(*), the NL2SOL real array.
!
  implicit none

  integer nn
  integer p

  real d(p)
  integer d0
  integer :: dfac = 41
  integer :: dtype = 16
  integer i
  integer iv(*)
  real j(nn,p)
  integer :: jtol0 = 86
  integer jtoli
  integer n
  integer :: niter = 31
  integer :: s = 53
  integer s1
  real sii
  real t
  real v(*)
  real v2norm
  real vdfac

  i = iv(dtype)

  if ( i /= 1 ) then

    if ( 0 < iv(niter) ) then
      return
    end if

  end if

  vdfac = v(dfac)
  d0 = jtol0 + p
  s1 = iv(s) - 1

  do i = 1, p

    s1 = s1 + i
    sii = v(s1)
    t = v2norm ( n, j(1,i) )

    if ( 0.0E+00 < sii ) then
      t = sqrt ( t * t + sii )
    end if

    jtoli = jtol0 + i
    d0 = d0 + 1

    if ( t < v(jtoli) ) then
      t = max ( v(d0), v(jtoli) )
    end if

    d(i) = max ( vdfac * d(i), t )

  end do

  return
end
subroutine gqtstp ( d, dig, dihdi, ka, l, p, step, v, w )

!*****************************************************************************80
!
!! GQTSTP computes the Goldfeld-Quandt-Trotter step by More-Hebden technique.
!
!  Discussion:
!
!    Given the compactly stored lower triangle of a scaled
!    hessian approximation and a nonzero scaled gradient vector,
!    this subroutine computes a Goldfeld-Quandt-Trotter step of
!    approximate length V(RADIUS) by the More-Hebden technique.
!
!    In other words, STEP is computed to approximately minimize
!      PSI(STEP) = G' * STEP + 0.5 * STEP' * H * STEP
!    such that the 2-norm of D * STEP is at most approximately V(RADIUS),
!    where G is the gradient, H is the hessian, and D is a diagonal
!    scale matrix whose diagonal is stored in the parameter D.
!
!    GQTSTP assumes:
!
!      DIG = inverse ( D ) * G,
!      DIHDI = inverse ( D ) * H * inverse ( D ).
!
!    If G = 0, however, STEP = 0 is returned, even at a saddle point.
!
!    If it is desired to recompute STEP using a different value of
!    V(RADIUS), then this routine may be restarted by calling it
!    with all parameters unchanged except V(RADIUS).  This explains
!    why STEP and W are listed as I/O.  On an initial call, with
!    KA < 0, STEP and W need not be initialized and only components
!    V(EPSLON), V(STPPAR), V(PHMNFC), V(PHMXFC), V(RADIUS), and
!    V(RAD0) of V must be initialized.  To compute STEP from a saddle
!    point, where the true gradient vanishes and H has a negative
!    eigenvalue, a nonzero G with small components should be passed.
!
!    This routine is called as part of the NL2SOL package, but it could
!    be used in solving any unconstrained minimization problem.
!
!    The desired G-Q-T step (references 2, 3, 4) satisfies
!    (H + ALPHA*D**2) * STEP = -G  for some nonnegative ALPHA such that
!    H + ALPHA*D**2 is positive semidefinite.  ALPHA and STEP are
!    computed by a scheme analogous to the one described in reference 5.
!    Estimates of the smallest and largest eigenvalues of the hessian
!    are obtained from the Gerschgorin circle theorem enhanced by a
!    simple form of the scaling described in reference 6.
!
!    Cases in which H + ALPHA*D**2 is nearly or exactly singular are
!    handled by the technique discussed in reference 2.  In these
!    cases, a step of exact length V(RADIUS) is returned for which
!    PSI(STEP) exceeds its optimal value by less than
!    -V(EPSLON)*PSI(STEP).
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Reference:
!
!    John Dennis, David Gay, Roy Welsch,
!    An Adaptive Nonlinear Least Squares Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, 1981.
!
!    David Gay,
!    Computing Optimal Locally Constrained Steps,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 2, Number 2, pages 186-197, 1981.
!
!    S M Goldfeld, R E Quandt, H F Trotter,
!    Maximization by Quadratic Hill-climbing,
!    Econometrica,
!    Volume 34, pages 541-551, 1966.
!
!    M D Hebden,
!    An Algorithm for Minimization using Exact Second Derivatives,
!    Report TP515,
!    Theoretical Physics Division,
!    AERE, Harwell, Oxon., England, 1973.
!
!    Jorge More,
!    The Levenberg-Marquardt Algorithm, Implementation and Theory,
!    in Springer Lecture Notes in Mathematics, Number 630,
!    edited by G A Watson,
!    Springer Verlag, Berlin and New York, pages 105-116, 1978.
!
!    Richard Varga,
!    Minimal Gerschgorin Sets,
!    Pacific Journal of Mathematics,
!    Volume 15, pages 719-729, 1965.
!
!  Parameters:
!
!    Input, real D(P), the scale vector, that is, the diagonal of the scale
!    matrix D mentioned above.
!
!    Input, real DIG(P), the scaled gradient vector, inverse ( D ) * G.
!    If G = 0, then STEP = 0 and V(STPPAR) = 0 are returned.
!
!    Input, real DIHDI((P*(P+1))/2), the lower triangle of the scaled
!    hessian approximation, that is,
!      inverse ( D ) * H * inverse ( D ),
!    stored compactly by rows, in the order (1,1), (2,1), (2,2), (3,1),
!    (3,2), and so on.
!
!    Input/output, integer KA, the number of Hebden iterations taken so
!    far to determine STEP.  KA < 0 on input means this is the first
!    attempt to determine STEP for the present DIG and DIHDI.
!    KA is initialized to 0 in this case.  Output with KA = 0  or
!    V(STPPAR) = 0 means STEP = -inverse(H)*G.
!
!     l (i/o) = workspace of length p*(p+1)/2 for cholesky factors.
!
!     p (in)  = number of parameters -- the hessian is a  p x p  matrix.
!
!  step (i/o) = the step computed.
!
!     v (i/o) contains various constants and variables described below.
!
!     w (i/o) = workspace of length 4*p + 6.
!
!  entries in v
!
! v(dgnorm) (i/o) = 2-norm of (d**-1)*g.
! v(dstnrm) (output) = 2-norm of d * step.
! v(dst0)   (i/o) = 2-norm of d*(h**-1)*g (for pos. def. h only), or
!             overestimate of smallest eigenvalue of (d**-1)*h*(d**-1).
! v(epslon) (in)  = max. relative error allowed for psi(step).  for the
!             step returned, psi(step) will exceed its optimal value
!             by less than -v(epslon)*psi(step).  suggested value = 0.1.
! v(gtstep) (out) = inner product between g and step.
! v(nreduc) (out) = psi(-(h**-1)*g) = psi(Newton step)  (for pos. def.
!             h only -- v(nreduc) is set to zero otherwise).
! v(phmnfc) (in)  = tol. (together with v(phmxfc)) for accepting step
!             (More's sigma).  the error v(dstnrm) - v(radius) must lie
!             between v(phmnfc)*v(radius) and v(phmxfc)*v(radius).
! v(phmxfc) (in)  (see v(phmnfc).)
!             suggested values -- v(phmnfc) = -0.25, v(phmxfc) = 0.5.
! v(preduc) (out) = psi(step) = predicted obj. func. reduction for step.
! v(radius) (in)  = radius of current (scaled) trust region.
! v(rad0)   (i/o) = value of v(radius) from previous call.
! v(STPPAR) (i/o) is normally the Marquardt parameter, i.e. the alpha
!             described below under algorithm notes.  if h + alpha*d**2
!             (see algorithm notes) is (nearly) singular, however,
!             then v(STPPAR) = -alpha.
!
  implicit none

  integer p

  real aki
  real akk
  real alphak
  real d(p)
  real delta
  integer dggdmx
  integer, parameter :: dgnorm = 1
  real, save :: dgxfac = 0.0E+00
  integer diag
  integer diag0
  real dig(p)
  real dihdi((p*(p+1))/2)
  real dotprd
  real dst
  integer, parameter :: dst0 = 3
  integer, parameter :: dstnrm = 2
  integer dstsav
  integer emax
  integer emin
  real, parameter :: epsfac = 50.0E+00
  integer, parameter :: epslon = 19
  real epso6
  integer, parameter :: gtstep = 4
  integer i
  integer inc
  integer irc
  integer j
  integer k
  integer k1
  integer ka
  integer kalim
  real, parameter :: kappa = 2.0E+00
  real l((p*(p+1))/2)
  real lk
  integer lk0
  real lsvmin
  integer, parameter :: nreduc = 6
  real oldphi
  real phi
  real phimax
  real phimin
  integer phipin
  integer, parameter :: phmnfc = 20
  integer, parameter :: phmxfc = 21
  integer, parameter :: preduc = 7
  real psifac
  integer q
  integer q0
  real rad
  integer, parameter :: rad0 = 9
  integer, parameter :: radius = 8
  logical restrt
  real root
  real si
  real sk
  real step(p)
  integer, parameter :: stppar = 5
  real sw
  real t
  real t1
  real twopsi
  real uk
  integer uk0
  real v(21)
  real v2norm
  real w(4*p+7)
  real wi
  integer x
  integer x0
!
!  Store largest absolute entry in inverse(D)*H*inverse(D) at W(DGGDMX).
!
  dggdmx = p + 1
!
!  Store Gerschgorin over- and underestimates of the largest
!  and smallest eigenvalues of inverse(D)*H*inverse(D) at W(EMAX)
!  and W(EMIN) respectively.
!
  emax = dggdmx + 1
  emin = emax + 1
!
!  For use in recomputing step, the final values of LK, UK, DST,
!  and the inverse derivative of More's PHI at 0, for positive definite
!  H, are stored in W(LK0), W(UK0), W(DSTSAV), and W(PHIPIN)
!  respectively.
!
  lk0 = emin + 1
  phipin = lk0 + 1
  uk0 = phipin + 1
  dstsav = uk0 + 1
!
!  Store diagonal of inverse(D)*H*inverse(D) in W(DIAG:DIAG+P-1).
!
  diag0 = dstsav
  diag = diag0 + 1
!
!  Store -D * STEP in W(Q:Q+P-1).
!
  q0 = diag0 + p
  q = q0 + 1
  rad = v(radius)
!
!  PHITOL = maximum error allowed in DST = V(DSTNRM) = 2-norm of
!  D * STEP.
!
  phimax = v(phmxfc) * rad
  phimin = v(phmnfc) * rad
!
!  EPSO6 and PSIFAC are used in checking for the special case
!  of nearly singular H + ALPHA*D**2.  See reference 2.
!
  psifac = 2.0E+00 * v(epslon) &
    / ( 3.0E+00 * ( 4.0E+00 * ( v(phmnfc) + 1.0E+00 ) * &
    ( kappa + 1.0E+00 )  +  kappa  +  2.0E+00 ) * rad**2 )
!
!  OLDPHI is used to detect limits of numerical accuracy.  If
!  we recompute step and it does not change, then we accept it.
!
  oldphi = 0.0E+00
  epso6 = v(epslon) / 6.0E+00
  irc = 0
  restrt = .false.
  kalim = ka + 50
!
!  Start or restart, depending on KA.
!
  if ( 0 <= ka ) then
!
!  Restart with new radius.
!
!  Prepare to return Newton step.
!
    if ( 0.0E+00 < v(dst0) .and. v(dst0) - rad <= phimax ) then

      restrt = .true.
      ka = ka + 1
      k = 0
      do i = 1, p
        k = k + i
        j = diag0 + i
        dihdi(k) = w(j)
      end do
      uk = -1.0E+00
      go to 40

    end if

    if ( ka == 0 ) then
      go to 60
    end if

    dst = w(dstsav)
    alphak = abs ( v(stppar) )
    phi = dst - rad
    t = v(dgnorm) / rad
!
!  Smaller radius.
!
    if ( rad <= v(rad0) ) then

      uk = t - w(emin)
      lk = 0.0E+00
      if ( 0.0E+00 < alphak ) then
        lk = w(lk0)
      end if
      lk = max ( lk, t - w(emax) )
      if ( 0.0E+00 < v(dst0) ) then
        lk = max ( lk, ( v(dst0) - rad ) * w(phipin) )
      end if
!
!  Bigger radius.
!
    else

      uk = t - w(emin)
      if ( 0.0E+00 < alphak ) then
        uk = min ( uk, w(uk0) )
      end if
      lk = max ( 0.0E+00, -v(dst0), t - w(emax) )
      if ( 0.0E+00 < v(dst0) ) then
        lk = max ( lk, (v(dst0)-rad)*w(phipin) )
      end if

    end if

    go to 260

  end if
!
!  Fresh start.
!
  k = 0
  uk = -1.0E+00
  ka = 0
  kalim = 50
!
!  Store diagonal of DIHDI in W(DIAG0+1:DIAG0+P).
!
  j = 0
  do i = 1, p
    j = j + i
    k1 = diag0 + i
    w(k1) = dihdi(j)
  end do
!
!  Determine W(DGGDMX), the largest element of DIHDI.
!
  t1 = 0.0E+00
  j = p * (p + 1) / 2
  do i = 1, j
    t = abs(dihdi(i))
    t1 = max ( t1, t )
  end do
  w(dggdmx) = t1
!
!  Try ALPHA = 0.
!
40 continue

  call lsqrt ( 1, p, l, dihdi, irc )
!
!  Indefinite H.  Underestimate smallest eigenvalue, use this
!  estimate to initialize lower bound LK on ALPHA.
!
  if ( irc == 0 ) go to 60

  j = ( irc * ( irc + 1 ) ) / 2
  t = l(j)
  l(j) = 1.0E+00
  w(1:irc) = 0.0E+00
  w(irc) = 1.0E+00
  call litvmu(irc, w, l, w)
  t1 = v2norm(irc, w)
  lk = -t / t1 / t1
  v(dst0) = -lk

  if (restrt) go to 210

  v(nreduc) = 0.0E+00
  go to 70
!
!  Positive definite H.  Compute unmodified Newton step.
!
60 continue

  lk = 0.0E+00
  call livmul(p, w(q), l, dig)
  v(nreduc) = 0.5E+00 * dotprd(p, w(q), w(q))
  call litvmu(p, w(q), l, w(q))
  dst = v2norm(p, w(q))
  v(dst0) = dst
  phi = dst - rad

  if ( phi <= phimax ) then
    alphak = 0.0E+00
    go to 290
  end if

  if (restrt) go to 210
!
!  Prepare to compute Gerschgorin estimates of largest and
!  smallest eigenvalues.
!
70 continue

  v(dgnorm) = v2norm ( p, dig )

  if ( v(dgnorm) == 0.0E+00 ) then
    v(stppar) = 0.0E+00
    v(preduc) = 0.0E+00
    v(dstnrm) = 0.0E+00
    v(gtstep) = 0.0E+00
    step(1:p) = 0.0E+00
    return
  end if

  k = 0
  do i = 1, p
    wi = 0.0E+00
    do j = 1, i - 1
      k = k + 1
      t = abs ( dihdi(k) )
      wi = wi + t
      w(j) = w(j) + t
    end do
    w(i) = wi
    k = k + 1
  end do
!
!  Underestimate smallest eigenvalue of inverse(D)*H*inverse(D).
!
  k = 1
  t1 = w(diag) - w(1)

  do i = 2, p
    j = diag0 + i
    t = w(j) - w(i)
    if ( t < t1 ) then
      t1 = t
      k = i
    end if
  end do

  sk = w(k)
  j = diag0 + k
  akk = w(j)
  k1 = ( k * ( k - 1 ) ) / 2 + 1
  inc = 1
  t = 0.0E+00

  do i = 1, p

    if ( i == k ) then
      inc = i
      k1 = k1 + inc
    else
      aki = abs(dihdi(k1))
      si = w(i)
      j = diag0 + i
      t1 = 0.5E+00 * (akk - w(j) + si - aki)
      t1 = t1 + sqrt(t1*t1 + sk*aki)
      if (t < t1) t = t1
      if ( k <= i ) then
        inc = i
      end if
      k1 = k1 + inc
    end if

  end do

  w(emin) = akk - t
  uk = v(dgnorm) / rad - w(emin)
!
!  Compute Gerschgorin overestimate of largest eigenvalue.
!
  k = 1
  t1 = w(diag) + w(1)

  do i = 2, p
    j = diag0 + i
    t = w(j) + w(i)
    if ( t1 < t ) then
      t1 = t
      k = i
    end if
  end do

  sk = w(k)
  j = diag0 + k
  akk = w(j)
  k1 = ( k * ( k - 1 ) ) / 2 + 1
  inc = 1
  t = 0.0E+00

  do i = 1, p
    if (i == k) then
      inc = i
      k1 = k1 + inc
    else
      aki = abs ( dihdi(k1) )
      si = w(i)
      j = diag0 + i
      t1 = 0.5E+00 * ( w(j) + si - aki - akk )
      t1 = t1 + sqrt ( t1 * t1 + sk * aki )
      if (t < t1) t = t1
      if ( k <= i ) then
        inc = i
      end if
      k1 = k1 + inc
    end if
  end do

  w(emax) = akk + t
  lk = max ( lk, v(dgnorm) / rad - w(emax) )
!
!  ALPHAK = current value of ALPHA.  We
!  use More's scheme for initializing it.
!
  alphak = abs ( v(stppar) ) * v(rad0) / rad
!
!  Compute L0 for positive definite H.
!
  if ( irc == 0 ) then

    call livmul(p, w, l, w(q))
    t = v2norm(p, w)
    w(phipin) = dst / t / t
    lk = max ( lk, phi * w(phipin) )

  end if
!
!  Safeguard ALPHAK and add ALPHAK*IDENTITY to inverse(D)*H*inverse(D).
!
210 continue

  ka = ka + 1

  if ( -v(dst0) >= alphak .or. alphak < lk .or. alphak >= uk ) then
    alphak = uk * max ( 0.001E+00, sqrt ( lk / uk ) )
  end if

  k = 0
  do i = 1, p
    k = k + i
    j = diag0 + i
    dihdi(k) = w(j) + alphak
  end do
!
!  Try computing Cholesky decomposition.
!
  call lsqrt(1, p, l, dihdi, irc)
!
!  inverse(D)*H*inverse(D) + ALPHAK*IDENTITY  is indefinite.  Overestimate
!  smallest eigenvalue for use in updating LK.
!
  if ( irc /= 0 ) then

    j = ( irc * ( irc + 1 ) ) / 2
    t = l(j)
    l(j) = 1.0E+00
    w(1:irc) = 0.0E+00
    w(irc) = 1.0E+00
    call litvmu ( irc, w, l, w )
    t1 = v2norm ( irc, w )
    lk = alphak - t / t1 / t1
    v(dst0) = -lk
    go to 210

  end if
!
!  ALPHAK makes inverse(D)*H*inverse(D) positive definite.
!  Compute Q = -D * STEP, check for convergence.
!
  call livmul(p, w(q), l, dig)
  call litvmu(p, w(q), l, w(q))
  dst = v2norm(p, w(q))
  phi = dst - rad

  if (phi <= phimax .and. phi >= phimin) go to 290

  if (phi == oldphi) go to 290

  oldphi = phi

  if ( phi > 0.0E+00 ) go to 260
!
!  Check for the special case of H + ALPHA*D**2  (nearly)
!  singular.  delta is >= the smallest eigenvalue of
!  inverse(D)*H*inverse(D) + ALPHAK*IDENTITY.
!
  if ( v(dst0) > 0.0E+00 ) go to 260

  delta = alphak + v(dst0)
  twopsi = alphak * dst * dst + dotprd ( p, dig, w(q) )

  if ( delta < psifac*twopsi ) then
    go to 270
  end if
!
!  Unacceptable ALPHAK.  Update LK, UK, ALPHAK.
!
260 continue

  if (ka >= kalim) go to 290

  call livmul(p, w, l, w(q))
  t1 = v2norm(p, w)
!
!  The following min is necessary because of restarts.
!
  if ( phi < 0.0E+00 ) then
    uk = min ( uk, alphak )
  end if

  alphak = alphak + ( phi / t1 ) * ( dst / t1 ) * ( dst / rad )
  lk = max ( lk, alphak )
  go to 210
!
!  Decide how to handle nearly singular H + ALPHA*D**2.
!
!  If not yet available, obtain machine dependent value dgxfac.
!
270 continue

  if ( dgxfac == 0.0E+00 ) then
    dgxfac = epsfac * epsilon ( dgxfac )
  end if
!
!  Is DELTA so small we cannot handle the special case in
!  the available arithmetic?  If so, accept STEP as it is.
!
  if ( dgxfac * w(dggdmx) < delta ) then
!
!  Handle nearly singular H + ALPHA*D**2.
!  Negate ALPHAK to indicate special case.
!
    alphak = -alphak
!
!  Allocate storage for scratch vector X.
!
    x0 = q0 + p
    x = x0 + 1
!
!  Use inverse power method with start from LSVMIN to obtain
!  approximate eigenvector corresponding to smallest eigenvalue
!  of inverse ( D ) * H * inverse ( D ).
!
    delta = kappa * delta
    t = lsvmin(p, l, w(x), w)
    k = 0
!
!  Normalize W.
!
    do

      w(1:p) = t * w(1:p)
!
!  Complete current inverse power iteration.
!  Replace W by inverse ( L' ) * W.
!
      call litvmu ( p, w, l, w )
      t1 = 1.0E+00 / v2norm(p, w)
      t = t1 * t

      if ( t <= delta ) then
        exit
      end if

      if ( 30 < k ) then
        go to 290
      end if

      k = k + 1
!
!  Start next inverse power iteration by storing normalized W in X.
!
      do i = 1, p
        j = x0 + i
        w(j) = t1 * w(i)
      end do
!
!  Compute W = inverse ( L ) * X.
!
      call livmul(p, w, l, w(x))
      t = 1.0E+00 / v2norm(p, w)

    end do

    w(1:p) = t1 * w(1:p)
!
!  Now W is the desired approximate unit eigenvector and
!  T * X = ( inverse(D) * H * inverse(D) + ALPHAK * I ) * W.
!
    sw = dotprd ( p, w(q), w )
    t1 = ( rad + dst ) * ( rad - dst )
    root = sqrt ( sw * sw + t1 )
    if ( sw < 0.0E+00 ) then
      root = -root
    end if
    si = t1 / (sw + root)
!
!  Accept current step if adding SI * W would lead to a
!  further relative reduction in PSI of less than V(EPSLON) / 3.
!
    v(preduc) = 0.5E+00 * twopsi
    t1 = 0.0E+00
    t = si * ( alphak * sw &
      - 0.5E+00 * si * ( alphak + t * dotprd ( p, w(x), w ) ) )

    if ( epso6 * twopsi <= t ) then
      v(preduc) = v(preduc) + t
      dst = rad
      t1 = -si
    end if

    do i = 1, p
      j = q0 + i
      w(j) = t1 * w(i) - w(j)
      step(i) = w(j) / d(i)
    end do

    v(gtstep) = dotprd ( p, dig, w(q) )
!
!  Save values for use in a possible restart.
!
    v(dstnrm) = dst
    v(stppar) = alphak
    w(lk0) = lk
    w(uk0) = uk
    v(rad0) = rad
    w(dstsav) = dst
!
!  Restore diagonal of DIHDI.
!
    j = 0
    do i = 1, p
      j = j + i
      k = diag0 + i
      dihdi(j) = w(k)
    end do

    return

  end if
!
!  Successful step.  Compute STEP = - inverse ( D ) * Q.
!
290 continue

  do i = 1, p
    j = q0 + i
    step(i) = -w(j) / d(i)
  end do
  v(gtstep) = -dotprd(p, dig, w(q))
  v(preduc) = 0.5E+00 * ( abs ( alphak ) *dst*dst - v(gtstep))
!
!  Save values for use in a possible restart.
!
  v(dstnrm) = dst
  v(stppar) = alphak
  w(lk0) = lk
  w(uk0) = uk
  v(rad0) = rad
  w(dstsav) = dst
!
!  Restore diagonal of DIHDI.
!
  j = 0
  do i = 1, p
    j = j + i
    k = diag0 + i
    dihdi(j) = w(k)
  end do

  return
end
subroutine itsmry ( d, iv, p, v, x )

!*****************************************************************************80
!
!! ITSMRY prints an iteration summary.
!
!  Modified:
!
!    06 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input, real D(P), the scale vector.
!
!    Input/output, integer IV(*), the NL2SOL integer parameter array.
!
!    Input, integer P, the number of variables.
!
!    Input, real V(*), the NL2SOL real array.
!
!    Input, real X(P), the current estimate of the minimizer.
!
  implicit none

  integer p

  integer cov1
  integer :: covmat = 26
  integer :: covprt = 14
  integer :: covreq = 15
  real d(p)
  integer :: dstnrm = 2
  integer :: f = 10
  integer :: f0 = 13
  integer :: fdif = 11
  integer :: g = 28
  integer g1
  integer i
  integer i1
  integer ii
  integer iv(*)
  integer iv1
  integer j
  integer m
  character*7 model(6)
  integer :: needhd = 39
  integer nf
  integer :: nfcall = 6
  integer :: nfcov = 40
  integer ng
  integer :: ngcall = 30
  integer :: ngcov = 41
  integer :: niter = 31
  integer :: nreduc = 6
  real nreldf
  integer ol
  real oldf
  integer :: outlev = 19
  integer :: preduc = 7
  real preldf
  integer :: prntit = 48
  integer :: prunit = 21
  integer pu
  real reldf
  integer :: reldx = 17
  integer :: size = 47
  integer :: solprt = 22
  integer :: statpr = 23
  integer :: stppar = 5
  integer :: sused = 57
  real v(*)
  real x(p)
  integer :: x0prt = 24

  data model / &
    '      G', &
    '      S', &
    '    G-S', &
    '    S-G', &
    '  G-S-G', &
    '  S-G-S' /

  pu = iv(prunit)

  if ( pu == 0 ) then
    return
  end if

  iv1 = iv(1)
  ol = iv(outlev)

  if ( iv1 < 2 .or. 15 < iv1 ) then
    write ( pu, '(a,i5)' ) 'IV(1) = ', iv1
    return
  end if

  if ( ol == 0 ) go to 20

  if ( iv1 >= 12 ) go to 20

  if (iv1 >= 10 .and. iv(prntit) == 0) go to 20

  if ( iv1 <= 2 ) then
    iv(prntit) = iv(prntit) + 1
    if (iv(prntit) < abs ( ol ) ) then
      return
    end if
  end if

 10   continue

      nf = iv(nfcall) - abs ( iv(nfcov) )
      iv(prntit) = 0
      reldf = 0.0E+00
      preldf = 0.0E+00
      oldf = v(f0)

      if ( 0.0E+00 < oldf ) then
         reldf = v(fdif) / oldf
         preldf = v(preduc) / oldf
      end if
!
!  Print short summary line.
!
      if ( ol <= 0 ) then

         if ( iv(needhd) == 1 ) then
           write ( pu, * ) ' '
           write ( pu, '(a)' ) &
           '    it    nf      f        reldf      preldf     reldx'
         end if

         iv(needhd) = 0
         write(pu,1017) iv(niter), nf, v(f), reldf, preldf, v(reldx)
!
!  Print long summary line.
!
      else

        if (iv(needhd) == 1) then
          write ( pu, * ) ' '
          write ( pu, * ) &
            '    it    nf      f        reldf      preldf     reldx' // &
            '    model    STPPAR      size      d*step     npreldf'
        end if

      iv(needhd) = 0
      m = iv(sused)
      if ( 0.0E+00 < oldf ) then
        nreldf = v(nreduc) / oldf
      else
        nreldf = 0.0E+00
      end if

      write(pu,1017) iv(niter), nf, v(f), reldf, preldf, v(reldx), &
                     model(m), v(stppar), v(size), &
                     v(dstnrm), nreldf
 1017 format(1x,i5,i6,4e11.3,a7,4e11.3)

  end if

 20   continue

  if ( iv1 == 1 ) then

    return

  else if ( iv1 == 2 ) then

    return

  else if ( iv1 == 3 ) then

    write ( pu, * ) ' '
    write ( pu, '(a)' ) 'X-convergence.'

  else if ( iv1 == 4 ) then

    write ( pu, * ) ' '
    write ( pu, '(a)' ) 'Relative function convergence.'

  else if ( iv1 == 5 ) then

    write ( pu, * ) ' '
    write ( pu, '(a)' ) 'X- and relative function convergence.'

  else if ( iv1 == 6 ) then

    write ( pu, * ) ' '
    write ( pu, '(a)' ) 'Absolute function convergence.'

  else if ( iv1 == 7 ) then

    write ( pu, * ) ' '
    write ( pu, '(a)' ) 'Singular convergence.'

  else if ( iv1 == 8 ) then

    write ( pu, * ) ' '
    write ( pu, '(a)' ) 'False convergence.'

  else if ( iv1 == 9 ) then

    write ( pu, * ) ' '
    write ( pu, '(a)' ) 'Function evaluation limit.'

  else if ( iv1 == 10 ) then

    write ( pu, * ) ' '
    write ( pu, '(a)' ) 'Iteration limit.'

  else if ( iv1 == 11 ) then

    write ( pu, * ) ' '
    write ( pu, '(a)' ) 'Stopx.'

  else if ( iv1 == 14 ) then

    write ( pu, * ) ' '
    write ( pu, '(a)' ) 'Bad parameters to ASSESS.'
    return
!
!  Initial call on ITSMRY.
!
  else if ( iv1 == 12 .or. iv1 == 13 .or. iv1 == 15 ) then

    if ( iv1 == 15 ) then
      write ( pu, * ) ' '
      write ( pu, '(a)' ) 'J could not be computed.'
      if ( 0 < iv(niter) ) then
        go to 190
      end if
    end if

    if ( iv1 == 13 ) then
      write ( pu, * ) ' '
      write ( pu, '(a)' ) 'Initial sum of squares overflows.'
    end if

    if ( iv(x0prt) /= 0 ) then
      write ( pu, * ) ' '
      write ( pu, * ) '     I     Initial X(i)      D(i)'
      write ( pu, * ) ' '
      write(pu,1150) (i, x(i), d(i), i = 1, p)
    end if

 1150 format((1x,i5,e17.6,e14.3))

    if ( iv1 == 13 ) then
      return
    end if

    iv(needhd) = 0
    iv(prntit) = 0

    if ( ol == 0 ) then
      return
    else if ( ol < 0 ) then
      write ( pu, '(a)' ) ' '
      write ( pu, '(a)' ) &
        '    it    nf      f        reldf      preldf     reldx'
    else if ( 0 < ol ) then
      write ( pu, '(a)' ) ' '
      write ( pu, '(a)' ) &
        '    it    nf      f        reldf      preldf     reldx' // &
        '    model    STPPAR      size      d*step     npreldf'
    end if

    write ( pu, * ) ' '
    write(pu,1160) v(f)
 1160 format('     0     1',e11.3,11x,e11.3)
    return

  else

    return

  end if
!
!  Print various information requested on solution.
!
180 continue

      iv(needhd) = 1

      if ( iv(statpr) /= 0 ) then

         oldf = v(f0)

         if ( 0.0E+00 < oldf ) then
           preldf = v(preduc) / oldf
           nreldf = v(nreduc) / oldf
         else
           preldf = 0.0E+00
           nreldf = 0.0E+00
         end if

         nf = iv(nfcall) - iv(nfcov)
         ng = iv(ngcall) - iv(ngcov)
         write ( pu, * ) ' '
         write(pu,1180) v(f), v(reldx), nf, ng, preldf, nreldf
 1180 format(' function',e17.6,'   reldx',e20.6/' func. evals', &
         i8,9x,'grad. evals',i8/' preldf',e19.6,3x,'npreldf',e18.6)

         if ( 0 < iv(nfcov) ) then
           write ( pu, * ) ' '
           write ( pu, '(i5,a)' ) iv(nfcov), &
             ' extra function evaluations for covariance.'
         end if

         if ( 0 < iv(ngcov) ) then
           write ( pu, '(i5,a)' ) iv(ngcov), &
             ' extra gradient evaluations for covariance.'
         end if
      end if

 190  continue

      if ( iv(solprt) /= 0 ) then

         iv(needhd) = 1
         g1 = iv(g)

         write ( pu, '(a)' ) ' '
         write ( pu, '(a)' ) &
           '     I      Final X(I)        D(I)          G(I)'
         write ( pu, '(a)' ) ' '

         do i = 1, p
           write ( pu, '(i5,e17.6,2e14.3)' ) i, x(i), d(i), v(g1)
           g1 = g1 + 1
         end do

      end if

      if ( iv(covprt) == 0 ) then
        return
      end if

      cov1 = iv(covmat)
      iv(needhd) = 1

      if ( cov1 < 0 ) then

        if ( -1 == cov1 ) then
          write ( pu, '(a)' ) 'Indefinite covariance matrix'
        else if (-2 == cov1) then
          write ( pu, '(a)' ) 'Oversize steps in computing covariance'
        end if

      else if ( cov1 == 0 ) then

        write ( pu, '(a)' ) 'Covariance matrix not computed'

      else if ( 0 < cov1 ) then

        write ( pu, * ) ' '
        i = abs ( iv(covreq) )
        if ( i <= 1 ) then
          write ( pu, '(a)' ) 'Covariance = scale * H**-1 * (J'' * J) * H**-1'
        else if ( i == 2 ) then
          write ( pu, '(a)' ) 'Covariance = scale * inverse ( H )'
        else if ( 3 <= i ) then
          write ( pu, '(a)' ) 'Covariance = scale * inverse ( J'' * J )'
        end if
        write ( pu, * ) ' '

        ii = cov1 - 1
        if ( ol <= 0 ) then
          do i = 1, p
            i1 = ii + 1
            ii = ii + i
            write(pu,1270) i, v(i1:ii)
          end do
 1270 format(' row',i3,2x,5e12.4/(9x,5e12.4))
        else

          do i = 1, p
            i1 = ii + 1
            ii = ii + i
            write(pu,1250) i, v(i1:ii)
          end do

 1250 format(' row',i3,2x,9e12.4/(9x,9e12.4))

    end if

  end if

  return
end
subroutine linvrt ( n, lin, l )

!*****************************************************************************80
!
!! LINVRT computes the inverse of a lower triangular matrix.
!
!  Discussion:
!
!    LIN = inverse ( L ), both N by N lower triangular matrices stored
!    compactly by rows.  LIN and L may share the same storage.
!
!  Modified:
!
!    05 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input, integer N, the order of L and LIN.
!
!    Output, real LIN((N*(N+1))/2), the inverse of L, a lower triangular
!    matrix stored by rows.
!
!    Input, real L((N*(N+1))/2), a lower triangular matrix stored by rows.
!
  implicit none

  integer n

  integer i
  integer ii
  integer j0
  integer j1
  integer jj
  integer k
  integer k0
  real l((n*(n+1))/2)
  real lin((n*(n+1))/2)
  real t

  j0 = ( n * ( n + 1 ) ) / 2

  do ii = 1, n

    i = n + 1 - ii
    lin(j0) = 1.0E+00 / l(j0)

    if ( i <= 1 ) then
      return
    end if

    j1 = j0

    do jj = 1, i - 1

      t = 0.0E+00
      j0 = j1
      k0 = j1 - jj

      do k = 1, jj
        t = t - l(k0) * lin(j0)
        j0 = j0 - 1
        k0 = k0 + k - i
      end do

      lin(j0) = t / l(k0)

    end do

    j0 = j0 - 1

  end do

  return
end
subroutine litvmu ( n, x, l, y )

!*****************************************************************************80
!
!! LITVMU solves L' * X = Y, where L is a lower triangular matrix.
!
!  Discussion:
!
!    This routine solves L' * X = Y, where L is an N by N lower
!    triangular matrix stored compactly by rows.  X and Y may occupy
!    the same storage.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input, integer N, the order of L.
!
!    Output, real X(N), the solution.
!
!    Input, real L((N*(N+1))/2), the lower triangular matrix, stored
!    by rows.
!
!    Input, real Y(N), the right hand side.
!
  implicit none

  integer n

  integer i
  integer i0
  integer ii
  integer ij
  integer j
  real l((n*(n+1))/2)
  real x(n)
  real xi
  real y(n)

  x(1:n) = y(1:n)
  i0 = ( n * ( n + 1 ) ) / 2

  do ii = 1, n

    i = n + 1 - ii
    xi = x(i) / l(i0)
    x(i) = xi

    if ( i <= 1 ) then
      return
    end if

    i0 = i0 - i

    if ( xi /= 0.0E+00 ) then

      do j = 1, i - 1
        ij = i0 + j
        x(j) = x(j) - xi * l(ij)
      end do

    end if

  end do

  return
end
subroutine livmul ( n, x, l, y )

!*****************************************************************************80
!
!! LIVMUL solves L * X = Y, where L is a lower triangular matrix.
!
!  Discussion:
!
!    This routine solves L * X = Y, where L is an N by N lower
!    triangular matrix stored compactly by rows.  X and Y may occupy
!    the same storage.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input, integer N, the order of L.
!
!    Output, real X(N), the solution.
!
!    Input, real L((N*(N+1))/2), the lower triangular matrix, stored
!    by rows.
!
!    Input, real Y(N), the right hand side.
!
  implicit none

  integer n

  real dotprd
  integer i
  integer j
  real l((n*(n+1))/2)
  real t
  real x(n)
  real y(n)

  x(1) = y(1) / l(1)

  j = 1

  do i = 2, n
    t = dotprd ( i-1, l(j+1), x )
    j = j + i
    x(i) = ( y(i) - t ) / l(j)
  end do

  return
end
subroutine lmstep ( d, g, ierr, ipivot, ka, p, qtr, r, step, v, w )

!*****************************************************************************80
!
!! LMSTEP computes a Levenberg-Marquardt step by More-Hebden techniques.
!
!  Discussion:
!
!    Given the R matrix from the QR decomposition of a jacobian
!    matrix, J, as well as Q' times the corresponding
!    residual vector, RESID, this subroutine computes a Levenberg-
!    Marquardt step of approximate length V(RADIUS) by the More
!    technique.
!
!    If it is desired to recompute step using a different value of
!    V(RADIUS), then this routine may be restarted by calling it
!    with all parameters unchanged except V(RADIUS).  This explains
!    why many parameters are listed as I/O.  On an initial call
!    with KA = -1, the caller need only have initialized D, G, KA, P,
!    QTR, R, V(EPSLON), V(PHMNFC), V(PHMXFC), V(RADIUS), and V(RAD0).
!
!    This code implements the step computation scheme described in
!    refs. 2 and 4.  Fast Givens transformations (see reference 3,
!    pages 60-62) are used to compute step with a nonzero Marquardt
!    parameter.
!
!    A special case occurs if J is nearly singular and V(RADIUS)
!    is sufficiently large.  In this case the step returned is such
!    that  twonorm(R)**2 - twonorm(R - J * STEP)**2  differs from its
!    optimal value by less than V(EPSLON) times this optimal value,
!    where J and R denote the original jacobian and residual.  See
!    reference 2 for more details.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Reference:
!
!    John Dennis, David Gay, Roy Welsch,
!    An Adaptive Nonlinear Least Squares Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, 1981.
!
!    David Gay,
!    Computing Optimal Locally Constrained Steps,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 2, Number 2, pages 186-197, 1981.
!
!    Charles Lawson and Richard Hanson,
!    Solving Least Squares Problems,
!    Prentice Hall, 1974.
!
!    Jorge More,
!    The Levenberg-Marquardt Algorithm, Implementation and Theory,
!    in Springer Lecture Notes in Mathematics, Number 630,
!    edited by G A Watson,
!    Springer Verlag, Berlin and New York, pages 105-116, 1978.
!
!  Parameters:
!
!    Input, real D(P), the scale vector.
!
!    Input, real G(P), the gradient vector J'*R.
!
!   ierr (i/o) = return code from QRFACT or QRFGS -- 0 means r has
!             full rank.
!
! ipivot (i/o) = permutation array from QRFACT or QRFGS, which compute
!             qr decompositions with column pivoting.
!
!     ka (i/o).  ka < 0 on input means this is the first call on
!             lmstep for the current r and qtr.  on output ka con-
!             tains the number of Hebden iterations needed to determine
!             step.  ka = 0 means a Gauss-Newton step.
!
!      p (in)  = number of parameters.
!
!    qtr (in)  = Q' * residual.
!
!      r (in)  = the R matrix, stored compactly by columns.
!
!   step (out) = the Levenberg-Marquardt step computed.
!
!      v (i/o) contains various constants and variables described below.
!
!      w (i/o) = workspace of length p*(p+5)/2 + 4.
!
!  entries in v
!
! v(dgnorm) (i/o) = 2-norm of (d**-1)*g.
! v(dstnrm) (i/o) = 2-norm of d * step.
! v(dst0)   (i/o) = 2-norm of Gauss-Newton step (for nonsing. j).
! v(epslon) (in) = max. relative error allowed in twonorm(r)**2 minus
!             twonorm(r - j * step)**2.  (see algorithm notes below.)
! v(gtstep) (out) = inner product between G and STEP.
! v(nreduc) (out) = half the reduction in the sum of squares predicted
!             for a Gauss-Newton step.
! v(phmnfc) (in)  = tol. (together with v(phmxfc)) for accepting step
!             (More's sigma).  the error v(dstnrm) - v(radius) must lie
!             between v(phmnfc)*v(radius) and v(phmxfc)*v(radius).
! v(phmxfc) (in)  (see v(phmnfc).)
! v(preduc) (out) = half the reduction in the sum of squares predicted
!             by the step returned.
! v(radius) (in)  = radius of current (scaled) trust region.
! v(rad0)   (i/o) = value of v(radius) from previous call.
! v(STPPAR) (i/o) = Marquardt parameter (or its negative if the special
!             case mentioned below in the algorithm notes occurs).
!
  implicit none

  integer p

  real a
  real adi
  real alphak
  real b
  real d(p)
  real d1
  real d2
  real, parameter :: dfac = 256.0E+00
  real dfacsq
  integer, parameter :: dgnorm = 1
  real dotprd
  real dst
  integer, parameter :: dst0 = 3
  integer, parameter ::dstnrm = 2
  integer dstsav
  real dtol
  integer, parameter :: epslon = 19
  real g(p)
  integer i
  integer i1
  integer ierr
  integer ip1
  integer ipivot(p)
  integer j1
  integer k
  integer ka
  integer kalim
  integer l
  real lk
  integer lk0
  real oldphi
  real phi
  real phimax
  real phimin
  integer phipin
  integer pp1o2
  real psifac
  real qtr(p)
  real r(1)
  real rad
  integer, parameter :: rad0 = 9
  integer res
  integer res0
  integer rmat
  integer rmat0
  real si
  real sj
  real sqrtak
  real step(p)
  integer, parameter :: stppar = 5
  real t
  real twopsi
  real uk
  integer uk0
  real v(21)
  real v2norm
  real w(p*(p+5)/2 + 4)
  real wl
!
!  subscripts for v
!
      integer gtstep, nreduc, phmnfc, &
              phmxfc, preduc, radius

     parameter ( gtstep=4, nreduc=6, phmnfc=20 )
    parameter ( phmxfc=21, preduc=7, radius=8 )
!
!  For use in recomputing STEP, the final values of LK and UK,
!  the inverse derivative of More's PHI at 0 (for nonsingular J)
!  and the value returned as V(DSTNRM) are stored at W(LK0),
!  W(UK0), W(PHIPIN), and W(DSTSAV) respectively.
!
  lk0 = p + 1
  phipin = lk0 + 1
  uk0 = phipin + 1
  dstsav = uk0 + 1
  rmat0 = dstsav
!
!  A copy of the R matrix from the QR decomposition of J is
!  stored in W starting at W(RMAT), and a copy of the residual
!  vector is stored in W starting at W(RES).  The loops below
!  that update the QR decomposition for a nonzero Marquardt parameter
!  work on these copies.
!
  rmat = rmat0 + 1
  pp1o2 = ( p * ( p + 1 ) ) / 2
  res0 = pp1o2 + rmat0
  res = res0 + 1
  rad = v(radius)
  if ( 0.0E+00 < rad ) then
    psifac = v(epslon) &
      / ( ( 8.0E+00 * ( v(phmnfc) + 1.0E+00 ) + 3.0E+00 ) * rad**2)
  end if
  phimax = v(phmxfc) * rad
  phimin = v(phmnfc) * rad
!
!  DTOL, DFAC, and DFACSQ are used in rescaling the fast Givens
!  representation of the updated QR decomposition.
!
  dtol = 1.0E+00 / dfac
  dfacsq = dfac * dfac
!
!  OLDPHI is used to detect limits of numerical accuracy.  If
!  we recompute STEP and it does not change, then we accept it.
!
  oldphi = 0.0E+00
  lk = 0.0E+00
  uk = 0.0E+00
  kalim = ka + 12
!
!  Start or restart, depending on KA.
!
  if ( 0 < ka ) then
    go to 370
  end if
!
!  Fresh start.  Compute V(NREDUC).
!
  if ( ka < 0 ) then
    ka = 0
    kalim = 12
    k = p
    if ( ierr /= 0 ) then
      k = abs ( ierr ) - 1
    end if
    v(nreduc) = 0.5E+00 * dotprd ( k, qtr, qtr )
  end if
!
!  Set up to try initial Gauss-Newton step.
!
 20   continue

  v(dst0) = -1.0E+00
!
!  Compute Gauss-Newton step.
!
!  Note that the R matrix is stored compactly by columns in
!  R(1), R(2), R(3), ...  It is the transpose of a
!  lower triangular matrix stored compactly by rows, and we
!  treat it as such when using LITVMU and LIVMUL.
!
  if ( ierr == 0 ) then

    call litvmu ( p, w, r, qtr )
!
!  Temporarily store permuted -D * STEP in STEP.
!
    do i = 1, p
      j1 = ipivot(i)
      step(i) = d(j1) * w(i)
    end do

    dst = v2norm(p, step)
    v(dst0) = dst
    phi = dst - rad

    if ( phi <= phimax ) then
      go to 410
    end if
!
!  If this is a restart, go to 110.
!
    if ( 0 < ka ) then
      go to 110
    end if
!
!  Gauss-Newton step was unacceptable.  Compute L0.
!
    do i = 1, p
      j1 = ipivot(i)
      step(i) = d(j1) * ( step(i) / dst )
    end do

    call livmul ( p, step, r, step )
    t = 1.0E+00 / v2norm(p, step)
    w(phipin) = ( t / dst ) * t
    lk = phi * w(phipin)

  end if
!
!  Compute U0.
!
  w(1:p) = g(1:p) / d(1:p)
  v(dgnorm) = v2norm(p, w)
  uk = v(dgnorm) / rad
!
!  Special case.  RAD <= 0 or (G = 0 and J is singular).
!
  if ( uk <= 0.0E+00 ) then
    v(stppar) = 0.0E+00
    dst = 0.0E+00
    lk = 0.0E+00
    uk = 0.0E+00
    v(gtstep) = 0.0E+00
    v(preduc) = 0.0E+00
    step(1:p) = 0.0E+00
    v(dstnrm) = dst
    w(dstsav) = dst
    w(lk0) = lk
    w(uk0) = uk
    v(rad0) = rad
    return
  end if
!
!  ALPHAK will be used as the current Marquardt parameter.  We
!  use More's scheme for initializing it.
!
  alphak = abs ( v(stppar) ) * v(rad0) / rad
!
!  Top of loop.  Increment KA, copy R to RMAT, QTR to RES.
!
110  continue

  ka = ka + 1
  w(rmat:rmat+pp1o2-1) = r(1:pp1o2)
  w(res:res+p-1) = qtr(1:p)
!
!  Safeguard ALPHAK and initialize fast Givens scale vector.
!
      if (alphak <= 0.0E+00 .or. alphak < lk .or. alphak >= uk ) then
        alphak = uk * max ( 0.001E+00, sqrt ( lk / uk ) )
      end if

      sqrtak = sqrt(alphak)
      w(1:p) = 1.0E+00
!
!  Add ALPHAK * D and update QR decomposition using fast Givens transform.
!
      do i = 1, p
!
!  Generate, apply first Givens transformation for row I of ALPHAK * D.
!  Use STEP to store temporary row.
!
         l = ( i * ( i + 1 ) ) / 2 + rmat0
         wl = w(l)
         d2 = 1.0E+00
         d1 = w(i)
         j1 = ipivot(i)
         adi = sqrtak*d(j1)

         if ( abs(wl) <= adi ) go to 150

 130     continue

         a = adi / wl
         b = d2 * a / d1
         t = a * b + 1.0E+00

         if ( t <= 2.5E+00 ) then

           w(i) = d1 / t
           d2 = d2 / t
           w(l) = t * wl
           a = -a
           do j1 = i, p
              l = l + j1
              step(j1) = a * w(l)
           end do
           go to 170

         end if

 150     continue

         b = wl / adi
         a = d1 * b / d2
         t = a * b + 1.0E+00

         if (t > 2.5E+00 ) go to 130

         w(i) = d2 / t
         d2 = d1 / t
         w(l) = t * adi
         do j1 = i, p
              l = l + j1
              wl = w(l)
              step(j1) = -wl
              w(l) = a * wl
         end do

 170     continue

         if ( i == p ) then
           exit
         end if
!
!  Now use Givens transformations to zero elements of temporary row.
!
         ip1 = i + 1
         do i1 = i + 1, p
              l = ( i1 * ( i1 + 1 ) ) / 2 + rmat0
              wl = w(l)
              si = step(i1-1)
              d1 = w(i1)
!
!  Rescale row I1 if necessary.
!
              if ( d1 < dtol ) then
                d1 = d1 * dfacsq
                wl = wl / dfac
                k = l
                do j1 = i1, p
                  k = k + j1
                  w(k) = w(k) / dfac
                end do
              end if
!
!  Use Givens transformations to zero next element of temporary row.
!
              if (abs(si) > abs(wl)) go to 220

              if (si == 0.0E+00) go to 260

 200          continue

              a = si / wl
              b = d2 * a / d1
              t = a * b + 1.0E+00

              if ( t <= 2.5E+00 ) then

                w(l) = t * wl
                w(i1) = d1 / t
                d2 = d2 / t
                do j1 = i1, p
                   l = l + j1
                   wl = w(l)
                   sj = step(j1)
                   w(l) = wl + b * sj
                   step(j1) = sj - a*wl
                end do

                go to 240

              end if

 220          b = wl / si
              a = d1 * b / d2
              t = a * b + 1.0E+00

              if (t > 2.5E+00 ) go to 200

              w(i1) = d2 / t
              d2 = d1 / t
              w(l) = t * si
              do j1 = i1, p
                   l = l + j1
                   wl = w(l)
                   sj = step(j1)
                   w(l) = a * wl + sj
                   step(j1) = b * sj - wl
              end do
!
!  Rescale temporary row if necessary.
!
 240          continue

              if ( d2 < dtol ) then
                   d2 = d2*dfacsq
                   step(i1:p) = step(i1:p) / dfac
              end if

 260          continue

        end do
      end do
!
!  Compute step.
!
 280  continue

      call litvmu ( p, w(res), w(rmat), w(res) )
!
!  Recover STEP and store permuted -D * STEP at W(RES).
!
      do i = 1, p
         j1 = ipivot(i)
         k = res0 + i
         t = w(k)
         step(j1) = -t
         w(k) = t * d(j1)
      end do

      dst = v2norm(p, w(res))
      phi = dst - rad
      if (phi <= phimax .and. phi >= phimin) go to 430
      if (oldphi == phi) go to 430
      oldphi = phi
!
!  Check for and handle special case.
!
      if ( phi <= 0.0E+00 ) then

        if ( kalim <= ka ) then
          go to 430
        end if

        twopsi = alphak * dst * dst - dotprd ( p, step, g )

        if ( alphak < twopsi * psifac ) then
          v(stppar) = -alphak
          go to 440
        end if

      end if

      if ( phi < 0.0E+00 ) then
        uk = alphak
      end if

 320  continue

      do i = 1, p
         j1 = ipivot(i)
         k = res0 + i
         step(i) = d(j1) * ( w(k) / dst )
      end do

      call livmul(p, step, w(rmat), step)
      step(1:p) = step(1:p) / sqrt ( w(1:p) )
      t = 1.0E+00 / v2norm(p, step)
      alphak = alphak + t * phi * t / rad
      lk = max ( lk, alphak )
      go to 110
!
!  Restart.
!
 370  continue

      lk = w(lk0)
      uk = w(uk0)

      if (v(dst0) > 0.0E+00 .and. v(dst0) - rad <= phimax) then
        go to 20
      end if

      alphak = abs ( v(stppar) )
      dst = w(dstsav)
      phi = dst - rad
      t = v(dgnorm) / rad
!
!  Smaller radius.
!
      if ( rad <= v(rad0) ) then
         uk = t
         if ( alphak <= 0.0E+00 ) then
           lk = 0.0E+00
         end if
         if (v(dst0) > 0.0E+00) lk = max ( lk, (v(dst0)-rad)*w(phipin) )
         if ( phi < 0.0E+00 ) then
           uk = min ( uk, alphak )
         end if
         go to 320
      end if
!
!  Bigger radius.
!
      if (alphak <= 0.0E+00 .or. uk > t) then
        uk = t
      end if

      if (v(dst0) > 0.0E+00) then
        lk = max ( lk, (v(dst0)-rad)*w(phipin) )
      else
        lk = 0.0E+00
      end if

      if ( phi < 0.0E+00 ) then
        uk = min ( uk, alphak )
      end if

      go to 320
!
!  Acceptable Gauss-Newton step.  Recover step from W.
!
 410  continue

      alphak = 0.0E+00
      do i = 1, p
         j1 = ipivot(i)
         step(j1) = -w(i)
      end do
!
!  Save values for use in a possible restart.
!
 430  continue

  v(stppar) = alphak

 440  continue

  v(gtstep) = dotprd ( p, step, g )
  v(preduc) = 0.5E+00 * (alphak*dst*dst - v(gtstep))
  v(dstnrm) = dst
  w(dstsav) = dst
  w(lk0) = lk
  w(uk0) = uk
  v(rad0) = rad

  return
end
subroutine lsqrt ( n1, n, l, a, irc )

!*****************************************************************************80
!
!! LSQRT computes the Cholesky factor of a lower triangular matrix.
!
!  Discussion:
!
!    Compute rows N1 through N of the Cholesky factor L of
!    A = L * L', where L and the lower triangle of A are both
!    stored compactly by rows, and may occupy the same storage.
!
!    IRC = 0 means all went well.  IRC = J means the leading
!    principal J x J submatrix of A is not positive definite,
!    and L(J*(J+1)/2) contains the nonpositive reduced J-th diagonal.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input, integer N1, N, the first and last rows to be computed.
!
!    Output, real L((N*(N+1))/2), contains rows N1 through N of the
!    Cholesky factorization of A, stored compactly by rows as a lower
!    triangular matrix.
!
!    Input, real A((N*(N+1))/2), the matrix whose Cholesky factorization
!    is desired.
!
!    Output, integer IRC, an error flag.  If IRC = 0, then the factorization
!    was carried out successfully.  Otherwise, the principal J x J subminor
!    of A was not positive definite.
!
  implicit none

  integer n

  real a(n*(n+1)/2)
  integer i
  integer i0
  integer ij
  integer ik
  integer irc
  integer j
  integer j0
  integer jk
  integer k
  real l(n*(n+1)/2)
  integer n1
  real t
  real td

  i0 = ( n1 * ( n1 - 1 ) ) / 2

  do i = n1, n

    td = 0.0E+00
    j0 = 0

    do j = 1, i - 1

      t = 0.0E+00

      do k = 1, j - 1
        ik = i0 + k
        jk = j0 + k
        t = t + l(ik) * l(jk)
      end do

      ij = i0 + j
      j0 = j0 + j
      t = ( a(ij) - t ) / l(j0)
      l(ij) = t
      td = td + t * t

    end do

    i0 = i0 + i
    t = a(i0) - td

    if ( t <= 0.0E+00 ) then
      l(i0) = t
      irc = i
      return
    end if

    l(i0) = sqrt ( t )

  end do

  irc = 0

  return
end
function lsvmin ( p, l, x, y )

!*****************************************************************************80
!
!! LSVMIN estimates the smallest singular value of a lower triangular matrix.
!
!  Discussion:
!
!    This function returns a good over-estimate of the smallest
!    singular value of the packed lower triangular matrix L.
!
!    The matrix L is a lower triangular matrix, stored compactly by rows.
!
!    The algorithm is based on Cline, Moler, Stewart and Wilkinson,
!    with the additional provision that LSVMIN = 0 is returned if the
!    smallest diagonal element of L in magnitude is not more than the unit
!    roundoff times the largest.
!
!    The algorithm uses a random number generator proposed by Smith,
!    which passes the spectral test with flying colors; see Hoaglin and
!    Knuth.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Reference:
!
!    A Cline, Cleve Moler, Pete Stewart, James Wilkinson,
!    An Estimate of the Condition Number of a Matrix,
!    Report TM-310,
!    Applied Math Division,
!    Argonne National Laboratory, 1977.
!
!    D C Hoaglin,
!    Theoretical Properties of Congruential Random-Number Generators,
!    An Empirical View,
!    Memorandum NS-340,
!    Department of Statistics,
!    Harvard University, 1976.
!
!    D E Knuth,
!    The Art of Computer Programming,
!    Volume 2, Seminumerical Algorithms,
!    Addison Wesley, 1969.
!
!    C S Smith,
!    Multiplicative Pseudo-Random Number Generators with Prime Modulus,
!    Journal of the Association for Computing Machinery,
!    Volume 19, pages 586-593, 1971.
!
!  Parameters:
!
!    Input, integer P, the order of L.
!
!    Input, real L((P*(P+1))/2), the elements of the lower triangular
!    matrix in row order, that is, L(1,1), L(2,1), L(2,2), L(3,1), L(3,2),
!    L(3,3), and so on.
!
!    Output, real X(P).  If LSVMIN returns a positive value, then X
!    is a normalized approximate left singular vector corresponding to
!    the smallest singular value.  This approximation may be very
!    crude.  If LSVMIN returns zero, then some components of X are zero
!    and the rest retain their input values.
!
!    Output, real Y(P).  If LSVMIN returns a positive value, then
!    Y = inverse ( L ) * X is an unnormalized approximate right singular
!    vector corresponding to the smallest singular value.  This
!    approximation may be crude.  If LSVMIN returns zero, then Y
!    retains its input value.  The caller may pass the same vector for X
!    and Y, in which case Y overwrites X, for nonzero LSVMIN returns.
!
  implicit none

  integer p

  real b
  integer i
  integer ii
  integer, save :: ix = 2
  integer j
  integer j0
  integer ji
  integer jj
  integer jjj
  real l((p*(p+1))/2)
  real lsvmin
  integer pplus1
  real psj
  real sminus
  real splus
  real t
  real v2norm
  real x(p)
  real xminus
  real xplus
  real y(p)
!
!  First check whether to return LSVMIN = 0 and initialize X.
!
  ii = 0

  do i = 1, p

    x(i) = 0.0E+00
    ii = ii + i

    if ( l(ii) == 0.0E+00 ) then
      lsvmin = 0.0E+00
      return
    end if

  end do

  if ( mod ( ix, 9973 ) == 0 ) then
    ix = 2
  end if
!
!  Solve L' * X = B, where the components of B have randomly
!  chosen magnitudes in ( 0.5, 1 ) with signs chosen to make X large.
!
  do j = p, 1, -1
!
!  Determine X(J) in this iteration.  Note for I = 1, 2,..., J
!  that X(I) holds the current partial sum for row I.
!
    ix = mod ( 3432 * ix, 9973 )
    b = 0.5E+00 * ( 1.0E+00 + real ( ix ) / 9973.0E+00 )
    xplus = ( b - x(j) )
    xminus = ( -b - x(j) )
    splus = abs ( xplus )
    sminus = abs ( xminus )
    j0 = ( j * ( j - 1 ) ) / 2
    jj = j0 + j
    xplus = xplus / l(jj)
    xminus = xminus / l(jj)

    do i = 1, j - 1
      ji = j0 + i
      splus = splus + abs ( x(i) + l(ji) * xplus )
      sminus = sminus + abs ( x(i) + l(ji) * xminus )
    end do

    if ( splus < sminus ) then
      xplus = xminus
    end if

    x(j) = xplus
!
!  Update partial sums.
!
    do i = 1, j - 1
      ji = j0 + i
      x(i) = x(i) + l(ji) * xplus
    end do

  end do
!
!  Normalize X.
!
  t = 1.0E+00 / v2norm ( p, x )
  x(1:p) = t * x(1:p)
!
!  Solve L * Y = X.
!  return SVMIN = 1 / twonorm ( Y ).
!
  do j = 1, p

    psj = 0.0E+00
    j0 = ( j * ( j - 1 ) ) / 2

    do i = 1, j - 1
      ji = j0 + i
      psj = psj + l(ji) * y(i)
    end do

    jj = j0 + j
    y(j) = ( x(j) - psj ) / l(jj)

  end do

 lsvmin = 1.0E+00 / v2norm ( p, y )

  return
end
subroutine ltsqar ( n, a, l )

!*****************************************************************************80
!
!! LTSQAR sets A to the lower triangle of L' * L.
!
!  Discussion:
!
!    L is an N by N lower triangular matrix, stored by rows.
!
!    A is also stored by rows, and may share storage with L.
!
!  Modified:
!
!    03 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input, integer N, the order of L and A.
!
!    Output, real A((N*(N+1))/2), the lower triangle of L' * L,
!    stored by rows.
!
!    Input, real L((N*(N+1))/2), the lower triangular matrix,
!    stored by rows.
!
  implicit none

  integer n

  real a((n*(n+1))/2)
  integer i
  integer i1
  integer ii
  integer j
  integer k
  real l((n*(n+1))/2)
  integer m

  ii = 0

  do i = 1, n

    i1 = ii + 1
    ii = ii + i
    m = 1

    do j = i1, ii - 1
      do k = i1, j
        a(m) = a(m) + l(j) * l(k)
        m = m + 1
      end do
    end do

    do j = i1, ii
      a(j) = l(ii) * l(j)
    end do

  end do

  return
end
subroutine nl2sol ( n, p, x, calcr, calcj, iv, v, uiparm, urparm, ufparm )

!*****************************************************************************80
!
!! NL2SOL minimizes a nonlinear sum of squares using an analytic jacobian.
!
!  Purpose:
!
!    Given a P-vector X of parameters, CALCR computes an N-vector
!    R = R(X) of residuals corresponding to X.  R(X) probably arises
!    from a nonlinear model involving P parameters and N observations.
!
!    This routine interacts with NL2ITR to seek a parameter vector X
!    that minimizes the sum of the squares of the components of R(X),
!    i.e., that minimizes the sum-of-squares function
!    F(X) = R(X)' * R(X) / 2.  R(X) is assumed to be a twice
!    continuously differentiable function of X.
!
!    See reference 1 for a description of the algorithm used.
!    On problems which are naturally well scaled, better performance
!    may be obtained by setting V(D0INIT) = 1.0 and IV(DTYPE) = 0,
!    which will cause the scale vector D to be set to all ones.
!
!    After a return with IV(1) <= 11, it is possible to restart,
!    that is, to change some of the IV and V input values and continue
!    the algorithm from the point where it was interrupted.  IV(1)
!    should not be changed, nor should any entries of IV
!    and V other than the input values (those supplied by DFAULT).
!
!    Those who do not wish to write a CALCJ which computes the jacobian
!    matrix analytically should call NL2SNO rather than NL2SOL.
!    NL2SNO uses finite differences to compute an approximate jacobian.
!
!    Those who would prefer to provide R and J (the residual and
!    jacobian) by reverse communication rather than by writing subroutines
!    CALCR and CALCJ may call on NL2ITR directly.  See the comments at the
!    beginning of NL2ITR.
!
!    Those who use NL2SOL interactively may wish to supply their
!    own STOPX function, which should return TRUE if the break key
!    has been pressed since stopx was last invoked.  This makes it possible
!    to externally interrupt NL2SOL (which will return with
!    IV(1) = 11 if STOPX returns TRUE).
!
!    Storage for J is allocated at the end of V.  Thus the caller
!    may make V longer than specified above and may allow CALCJ to use
!    elements of J beyond the first N*P as scratch storage.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Reference:
!
!    John Dennis, David Gay, Roy Welsch,
!    An Adaptive Nonlinear Least Squares Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, 1981.
!
!  Parameters:
!
!    Input, integer N, the number of observations, that is, the number of
!    components in R(X).  P <= N.
!
!    Input, integer P, the number of parameters, or components in X.  P must
!    be positive.
!
!    Input/output, real X(P).  On input, X is an initial guess at the
!    desired parameter estimate.  On output, X contains the best parameter
!    estimate found.
!
!    Input, external CALCR, a subroutine which, given X, computes R(X).
!    CALCR must be declared external in the calling program.
!    It is invoked by
!      call calcr ( n, p, x, nf, r, uiparm, urparm, ufparm )
!    When CALCR is called, NF is the invocation count for CALCR.  It is
!    included for possible use with CALCJ.  If X is out of bounds, for
!    instance, if it would cause overflow in computing R(X), then CALCR
!    should set NF to 0.  This will cause a shorter step to be attempted.
!    The other parameters are as described above and below.  CALCR
!    should not change N, P, or X.
!
!    Input, external CALCJ, a subroutine which, given X, computes the
!    jacobian matrix J of R at X, that is, the N by P matrix whose
!    (I,K) entry is the partial derivative of the I-th component of R
!    with respect to X(K).  CALCJ must be declared external in the
!    calling program.  It is invoked by
!      call calcj(n,p,x,nf,j,uiparm,urparm,ufparm)
!    NF is the invocation count for CALCR at the time R(X) was evaluated.
!    The X passed to CALCJ is usually the one passed to CALC on either its
!    most recent invocation or the one prior to it.  If CALCR saves
!    intermediate results for use by CALCJ, then it is possible to tell
!    from NF whether they are valid for the current X (or which copy is
!    valid if two copies are kept).  If J cannot be computed at X,
!    then CALCJ should set NF to 0.  In this case, NL2SOL will return
!    with IV(1) = 15.  The other parameters to CALCJ are as described
!    above and below.  CALCJ should not change N, P, or X.
!
!    Input/output, integer IV(60+P), helps control the NL2SOL algorithm
!    and is used to store various intermediate quantities.  Of particular
!    interest are the initialization/return code IV(1) and the entries
!    in that control printing and limit the number of iterations and
!    function evaluations.  See the section on IV input values.
!
! v........ (input/output) a floating-point value array of length at
!                  least 93 + n*p + 3*n + p*(3*p+33)/2 that helps con-
!                  trol the nl2sol algorithm and that is used to store
!                  various intermediate quantities.  of particular in-
!                  terest are the entries in v that limit the length of
!                  the first step attempted (lmax0), specify conver-
!                  gence tolerances (afctol, rfctol, xctol, xftol),
!                  and help choose the step size used in computing the
!                  covariance matrix (delta0).  see the section on
!                  (selected) v input values below.
!
! uiparm... (input) user integer parameter array passed without change
!                  to calcr and calcj.
!
! urparm... (input) user floating-point parameter array passed without
!                  change to calcr and calcj.
!
! ufparm... (input) user external subroutine or function passed without
!                  change to calcr and calcj.
!
!  iv input values (from subroutine dfault)
!
! iv(1)...  on input, iv(1) should have a value between 0 and 12......
!             0 and 12 mean this is a fresh start.  0 means that
!             dfault(iv, v) is to be called to provide all default
!             values to iv and v.  12 (the value that dfault assigns to
!             iv(1)) means the caller has already called dfault(iv, v)
!             and has possibly changed some iv and/or v entries to non-
!             default values.  default = 12.
! iv(covprt)... iv(14) = 1 means print a covariance matrix at the solu-
!             tion.  (this matrix is computed just before a return with
!             iv(1) = 3, 4, 5, 6.)
!             iv(covprt) = 0 means skip this printing.  default = 1.
! iv(covreq)... iv(15) = nonzero means compute a covariance matrix
!             just before a return with iv(1) = 3, 4, 5, 6.  in
!             this case, an approximate covariance matrix is obtained
!             in one of several ways.  let k = abs(iv(covreq)) and let
!             scale = 2*f(x)/max(1,n-p),  where 2*f(x) is the residual
!             sum of squares.  if k = 1 or 2, then a finite-difference
!             hessian approximation h is obtained.  if h is positive
!             definite (or, for k = 3, if the jacobian matrix j at x
!             is nonsingular), then one of the following is computed...
!                  k = 1....  scale * h**-1 * (j**t * j) * h**-1.
!                  k = 2....  scale * h**-1.
!                  k = 3....  scale * (j**t * j)**-1.
!             (j**t is the transpose of j, while **-1 means inverse.)
!             if iv(covreq) is positive, then both function and grad-
!             ient values (calls on calcr and calcj) are used in com-
!             puting h (with step sizes determined using v(delta0) --
!             see below), while if iv(covreq) is negative, then only
!             function values (calls on calcr) are used (with step
!             sizes determined using v(dltfdc) -- see below).  if
!             iv(covreq) = 0, then no attempt is made to compute a co-
!             variance matrix (unless iv(covprt) = 1, in which case
!             iv(covreq) = 1 is assumed).  see iv(covmat) below.
!             default = 1.
! iv(dtype).... iv(16) tells how the scale vector D (see ref. 1) should
!             be chosen.  iv(dtype) >= 1 means choose d as described
!             below with v(dfac).  iv(dtype) <= 0 means the caller
!             has chosen d and has stored it in v starting at
!             v(94 + 2*n + p*(3*p + 31)/2).  default = 1.
! iv(inits).... iv(25) tells how the S matrix (see ref. 1) should be
!             initialized.  0 means initialize S to 0 (and start with
!             the Gauss-Newton model).  1 and 2 mean that the caller
!             has stored the lower triangle of the initial S rowwise in
!             v starting at v(87+2*p).  iv(inits) = 1 means start with
!             the Gauss-Newton model, while iv(inits) = 2 means start
!             with the augmented model (see ref. 1).  default = 0.
! iv(mxfcal)... iv(17) gives the maximum number of function evaluations
!             (calls on calcr, excluding those used to compute the co-
!             variance matrix) allowed.  if this number does not suf-
!             fice, then nl2sol returns with iv(1) = 9.  default = 200.
! iv(mxiter)... iv(18) gives the maximum number of iterations allowed.
!             it also indirectly limits the number of gradient evalua-
!             tions (calls on calcj, excluding those used to compute
!             the covariance matrix) to iv(mxiter) + 1.  if iv(mxiter)
!             iterations do not suffice, then nl2sol returns with
!             iv(1) = 10.  default = 150.
! iv(outlev)... iv(19) controls the number and length of iteration sum-
!             mary lines printed (by itsmry).  iv(outlev) = 0 means do
!             not print any summary lines.  otherwise, print a summary
!             line after each abs(iv(outlev)) iterations.  if iv(outlev)
!             is positive, then summary lines of length 117 (plus carri-
!             age control) are printed, including the following...  the
!             iteration and function evaluation counts, current func-
!             tion value (v(f) = half the sum of squares), relative
!             difference in function values achieved by the latest step
!             (i.e., reldf = (f0-v(f))/f0, where f0 is the function
!             value from the previous iteration), the relative function
!             reduction predicted for the step just taken (i.e.,
!             preldf = v(preduc) / f0, where v(preduc) is described
!             below), the scaled relative change in x (see v(reldx)
!             below), the models used in the current iteration (g =
!             Gauss-Newton, s=augmented), the Marquardt parameter
!             STPPAR used in computing the last step, the sizing factor
!             used in updating s, the 2-norm of the scale vector d
!             times the step just taken (see ref. 1), and npreldf, i.e.,
!             v(nreduc)/f0, where v(nreduc) is described below -- if
!             npreldf is positive, then it is the relative function
!             reduction predicted for a Newton step (one with
!             STPPAR = 0).  if npreldf is zero, either the gradient
!             vanishes (as does preldf) or else the augmented model
!             is being used and its hessian is indefinite (with preldf
!             positive).  if npreldf is negative, then it is the nega-
!             of the relative function reduction predicted for a step
!             computed with step bound v(lmax0) for use in testing for
!             singular convergence.
!                  if iv(outlev) is negative, then lines of maximum
!             length 79 (or 55 is iv(covprt) = 0) are printed, includ-
!             ing only the first 6 items listed above (through reldx).
!             default = 1.
! iv(parprt)... iv(20) = 1 means print any nondefault v values on a
!             fresh start or any changed v values on a restart.
!             iv(parprt) = 0 means skip this printing.  default = 1.
! iv(prunit)... iv(21) is the output unit number on which all printing
!             is done.  iv(prunit) = 0 means suppress all printing.
!             (setting iv(prunit) to 0 is the only way to suppress the
!             one-line termination reason message printed by itsmry.)
!             default = standard output unit (unit 6 on most systems).
! iv(solprt)... iv(22) = 1 means print out the value of x returned (as
!             well as the corresponding gradient and scale vector d).
!             iv(solprt) = 0 means skip this printing.  default = 1.
! iv(statpr)... iv(23) = 1 means print summary statistics upon return-
!             ing.  these consist of the function value (half the sum
!             of squares) at x, v(reldx) (see below), the number of
!             function and gradient evaluations (calls on calcr and
!             calcj respectively, excluding any calls used to compute
!             the covariance), the relative function reductions predict-
!             ed for the last step taken and for a Newton step (or per-
!             haps a step bounded by v(lmax0) -- see the descriptions
!             of preldf and npreldf under iv(outlev) above), and (if an
!             attempt was made to compute the covariance) the number of
!             calls on calcr and calcj used in trying to compute the
!             covariance.  iv(statpr) = 0 means skip this printing.
!             default = 1.
! iv(x0prt).... iv(24) = 1 means print the initial x and scale vector d
!             (on a fresh start only).  iv(x0prt) = 0 means skip this
!             printing.  default = 1.
!
!  (selected) iv output values
!
! iv(1)........ on output, iv(1) is a return code....
!             3 = x-convergence.  the scaled relative difference be-
!                  tween the current parameter vector x and a locally
!                  optimal parameter vector is very likely at most
!                  v(xctol).
!             4 = relative function convergence.  the relative differ-
!                  ence between the current function value and its lo-
!                  cally optimal value is very likely at most v(rfctol).
!             5 = both x- and relative function convergence (i.e., the
!                  conditions for iv(1) = 3 and iv(1) = 4 both hold).
!             6 = absolute function convergence.  the current function
!                  value is at most v(afctol) in absolute value.
!             7 = singular convergence.  the hessian near the current
!                  iterate appears to be singular or nearly so, and a
!                  step of length at most v(lmax0) is unlikely to yield
!                  a relative function decrease of more than v(rfctol).
!             8 = false convergence.  the iterates appear to be converg-
!                  ing to a noncritical point.  this may mean that the
!                  convergence tolerances (v(afctol), v(rfctol),
!                  v(xctol)) are too small for the accuracy to which
!                  the function and gradient are being computed, that
!                  there is an error in computing the gradient, or that
!                  the function or gradient is discontinuous near x.
!             9 = function evaluation limit reached without other con-
!                  vergence (see iv(mxfcal)).
!            10 = iteration limit reached without other convergence
!                  (see iv(mxiter)).
!            11 = stopx returned .true. (external interrupt).  see the
!                  usage notes below.
!            13 = f(x) cannot be computed at the initial x.
!            14 = bad parameters passed to assess (which should not
!                  occur).
!            15 = the jacobian could not be computed at x (see calcj
!                  above).
!            16 = n or p (or parameter nn to nl2itr) out of range --
!                  p <= 0 or n < p or nn < n.
!            17 = restart attempted with n or p (or par. nn to nl2itr)
!                  changed.
!            18 = iv(inits) is out of range.
!            19...45 = v(iv(1)) is out of range.
!            50 = iv(1) was out of range.
!            87...(86+p) = jtol(iv(1)-86) (i.e., v(iv(1)) is not
!                  positive (see v(dfac) below).
! iv(covmat)... iv(26) tells whether a covariance matrix was computed.
!             if (iv(covmat) is positive, then the lower triangle of
!             the covariance matrix is stored rowwise in v starting at
!             v(iv(covmat)).  if iv(covmat) = 0, then no attempt was
!             made to compute the covariance.  if iv(covmat) = -1,
!             then the finite-difference hessian was indefinite.  and
!             and if iv(covmat) = -2, then a successful finite-differ-
!             encing step could not be found for some component of x
!             (i.e., calcr set nf to 0 for each of two trial steps).
!             note that iv(covmat) is reset to 0 after each successful
!             step, so if such a step is taken after a restart, then
!             the covariance matrix will be recomputed.
! iv(d)........ iv(27) is the starting subscript in v of the current
!             scale vector d.
! iv(g)........ iv(28) is the starting subscript in v of the current
!             least-squares gradient vector (j**t)*r.
! iv(nfcall)... iv(6) is the number of calls so far made on calcr (i.e.,
!             function evaluations, including those used in computing
!             the covariance).
! iv(nfcov).... iv(40) is the number of calls made on calcr when
!             trying to compute covariance matrices.
! iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
!             calcj) so far done (including those used for computing
!             the covariance).
! iv(ngcov).... iv(41) is the number of calls made on calcj when
!             trying to compute covariance matrices.
! iv(niter).... iv(31) is the number of iterations performed.
! iv(r)........ iv(50) is the starting subscript in v of the residual
!             vector r corresponding to x.
!
! (selected) v input values (from subroutine dfault)
!
! v(afctol)... v(31) is the absolute function convergence tolerance.
!             if nl2sol finds a point where the function value (half
!             the sum of squares) is less than v(afctol), and if nl2sol
!             does not return with iv(1) = 3, 4, or 5, then it returns
!             with iv(1) = 6.  default = max(10**-20, machep**2), where
!             machep is the unit roundoff.
! v(delta0)... v(44) is a factor used in choosing the finite-difference
!             step size used in computing the covariance matrix when
!             iv(covreq) = 1 or 2.  for component i, step size
!                  v(delta0) * max(abs(x(i)), 1/d(i)) * sign(x(i))
!             is used, where d is the current scale vector (see ref. 1).
!             (if this step results in calcr setting nf to 0, then -0.5
!             times this step is also tried.)  default = machep**0.5,
!             where machep is the unit roundoff.
! v(dfac)..... v(41) and the d0 and jtol arrays (see v(d0init) and
!             v(jtinit)) are used in updating the scale vector d when
!             iv(dtype) > 0.  (d is initialized according to
!             v(dinit).)  let d1(i) =
!               max(sqrt(jcnorm(i)**2 + max(s(i,i),0)), v(dfac)*d(i)),
!             where jcnorm(i) is the 2-norm of the i-th column of the
!             current jacobian matrix and s is the s matrix of ref. 1.
!             if iv(dtype) = 1, then d(i) is set to d1(i) unless
!             d1(i) < jtol(i), in which case d(i) is set to
!                                max(d0(i), jtol(i)).
!             if iv(dtype) >= 2, then d is updated during the first
!             iteration as for iv(dtype) = 1 (after any initialization
!             due to v(dinit)) and is left unchanged thereafter.
!             default = 0.6.
! v(dinit).... v(38), if nonnegative, is the value to which the scale
!             vector d is initialized.  default = 0.
! v(dltfdc)... v(40) helps choose the step size used when computing the
!             covariance matrix when iv(covreq) = -1 or -2.  for
!             differences involving x(i), the step size first tried is
!                       v(dltfdc) * max(abs(x(i)), 1/d(i)),
!             where d is the current scale vector (see ref. 1).  (if
!             this step is too big the first time it is tried, i.e., if
!             calcr sets nf to 0, then -0.5 times this step is also
!             tried.)  default = machep**(1/3), where machep is the
!             unit roundoff.
! v(d0init)... v(37), if positive, is the value to which all components
!             of the d0 vector (see v(dfac)) are initialized.  if
!             v(dfac) = 0, then it is assumed that the caller has
!             stored d0 in v starting at v(p+87).  default = 1.0.
! v(jtinit)... v(39), if positive, is the value to which all components
!             of the jtol array (see v(dfac)) are initialized.  if
!             v(jtinit) = 0, then it is assumed that the caller has
!             stored jtol in v starting at v(87).  default = 10**-6.
! v(lmax0).... v(35) gives the maximum 2-norm allowed for d times the
!             very first step that nl2sol attempts.  it is also used
!             in testing for singular convergence -- if the function
!             reduction predicted for a step of length bounded by
!             v(lmax0) is at most v(rfctol) * abs(f0), where  f0  is
!             the function value at the start of the current iteration,
!             and if nl2sol does not return with iv(1) = 3, 4, 5, or 6,
!             then it returns with iv(1) = 7.    default = 100.
! v(rfctol)... v(32) is the relative function convergence tolerance.
!             if the current model predicts a maximum possible function
!             reduction (see v(nreduc)) of at most v(rfctol)*abs(f0) at
!             the start of the current iteration, where  f0  is the
!             then current function value, and if the last step attempt-
!             ed achieved no more than twice the predicted function
!             decrease, then nl2sol returns with iv(1) = 4 (or 5).
!             default = max(10**-10, machep**(2/3)), where machep is
!             the unit roundoff.
! v(tuner1)... v(26) helps decide when to check for false convergence
!             and to consider switching models.  this is done if the
!             actual function decrease from the current step is no more
!             than v(tuner1) times its predicted value.  default = 0.1.
! v(xctol).... v(33) is the x-convergence tolerance.  if a Newton step
!             (see v(nreduc)) is tried that has v(reldx) <= v(xctol)
!             and if this step yields at most twice the predicted func-
!             tion decrease, then nl2sol returns with iv(1) = 3 (or 5).
!             (see the description of v(reldx) below.)
!             default = machep**0.5, where machep is the unit roundoff.
! v(xftol).... v(34) is the false convergence tolerance.  if a step is
!             tried that gives no more than v(tuner1) times the predict-
!             ed function decrease and that has v(reldx) <= v(xftol),
!             and if nl2sol does not return with iv(1) = 3, 4, 5, 6, or
!             7, then it returns with iv(1) = 8.  (see the description
!             of v(reldx) below.)  default = 100*machep, where
!             machep is the unit roundoff.
! v(*)........ dfault supplies to v a number of tuning constants, with
!             which it should ordinarily be unnecessary to tinker.  see
!             version 2.2 of the nl2sol usage summary (which is an
!             appendix to ref. 1).
!
!  (selected) v output values
!
! v(dgnorm)... v(1) is the 2-norm of (d**-1)*g, where g is the most re-
!             cently computed gradient and d is the corresponding scale
!             vector.
! v(dstnrm)... v(2) is the 2-norm of d * step, where step is the most re-
!             cently computed step and d is the current scale vector.
! v(f)........ v(10) is the current function value (half the sum of
!             squares).
! v(f0)....... v(13) is the function value at the start of the current
!             iteration.
! v(nreduc)... v(6), if positive, is the maximum function reduction
!             possible according to the current model, i.e., the func-
!             tion reduction predicted for a Newton step (i.e.,
!             step = -h**-1 * g,  where  g = (j**t) * r  is the current
!             gradient and h is the current hessian approximation --
!             h = (j**t)*j  for the Gauss-Newton model and
!             h = (j**t)*j + s  for the augmented model).
!                  v(nreduc) = zero means h is not positive definite.
!                  if v(nreduc) is negative, then it is the negative of
!             the function reduction predicted for a step computed with
!             a step bound of v(lmax0) for use in testing for singular
!             convergence.
! v(preduc)... v(7) is the function reduction predicted (by the current
!             quadratic model) for the current step.  this (divided by
!             v(f0)) is used in testing for relative function
!             convergence.
! v(reldx).... v(17) is the scaled relative change in x caused by the
!             current step, computed as
!                  max(abs(d(i)*(x(i)-x0(i)), 1 <= i <= p) /
!                     max(d(i)*(abs(x(i))+abs(x0(i))), 1 <= i <= p),
!             where x = x0 + step.
!
  implicit none

  integer p

  external calcj
  external calcr
  integer, parameter :: d = 27
  integer d1
  integer iv(60+p)
  integer, parameter :: j = 33
  integer j1
  integer n
  integer nf
  integer, parameter :: nfcall = 6
  integer, parameter :: nfgcal = 7
  integer, parameter :: r = 50
  integer r1
  logical strted
  integer, parameter :: toobig = 2
  external ufparm
  integer uiparm(*)
  real urparm(*)
  real v(93 + n*p + 3*n + (p*(3*p+33))/2)
  real x(p)

  d1 = 94 + 2*n + ( p * ( 3 * p + 31 ) ) / 2
  iv(d) = d1
  r1 = d1 + p
  iv(r) = r1
  j1 = r1 + n
  iv(j) = j1
  strted = .true.

  if ( iv(1) /= 0 .and. iv(1) /= 12 ) go to 40

         strted = .false.
         iv(nfcall) = 1
         iv(nfgcal) = 1

 10   continue

      nf = iv(nfcall)

      call calcr(n, p, x, nf, v(r1), uiparm, urparm, ufparm)

      if ( strted ) then

        if ( nf <= 0 ) then
          iv(toobig) = 1
        end if

        go to 40

      end if

      if ( nf <= 0 ) then
        iv(1) = 13
        call itsmry ( v(d1), iv, p, v, x )
        return
      end if

30    continue

      call calcj(n, p, x, iv(nfgcal), v(j1), uiparm, urparm, ufparm)

      if ( iv(nfgcal) == 0 ) then
        iv(1) = 15
        call itsmry ( v(d1), iv, p, v, x )
        return
      end if

      strted = .true.

 40   continue

      call nl2itr ( v(d1), iv, v(j1), n, n, p, v(r1), v, x )

  if ( iv(1) == 2 ) then
    go to 30
  end if

  if ( iv(1) < 2 ) then
    go to 10
  end if

  return
end
subroutine nl2sno ( n, p, x, calcr, iv, v, uiparm, urparm, ufparm )

!*****************************************************************************80
!
!! NL2SNO is like NL2SOL, but uses a finite difference jacobian.
!
!  Discussion:
!
!    NL2SNO is like NL2SOL, but without calcj -- minimize nonlinear sum of
!    squares using finite-difference jacobian approximations
!
!    The parameters for NL2SNO are the same as those for NL2SOL
!    except that CALCJ is omitted.  Instead of calling
!    CALCJ to obtain the jacobian matrix of R at X, NL2SNO computes
!    an approximation to it by forward finite differences.  See
!    V(DLTFDJ) below.  NL2SNO uses function values only when comput-
!    the covariance matrix, rather than the functions and gradients
!    that NL2SOL may use.  To do so, NL2SNO sets IV(COVREQ) to -1 if
!    IV(COVPRT) = 1 with IV(COVREQ) = 0 and to minus its absolute
!    value otherwise.  Thus V(DELTA0) is never referenced and only
!    V(DLTFDC) matters.  See NL2SOL for a description of V(DLTFDC).
!
!    The number of extra calls on CALCR used in computing the jacobian
!    approximation are not included in the function evaluation
!    count IV(NFCALL) and are not otherwise reported.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Reference:
!
!    John Dennis, David Gay, Roy Welsch,
!    An Adaptive Nonlinear Least Squares Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, 1981.
!
!  Parameters:
!
!    V(DLTFDJ) helps choose the step size used when computing the
!    finite-difference jacobian matrix.  For differences involving X(I),
!    the step size first tried is
!      V(DLTFDJ) * max ( abs ( X(I) ), 1/D(I)),
!    where D is the current scale vector; see reference 1.  If this step is
!    too big, so that CALCR sets NF to 0, then smaller steps are tried
!    until the step size is shrunk below 1000 * MACHEP, where MACHEP
!    is the unit roundoff.  Default = sqrt ( MACHEP ).
!
  implicit none

  integer p

  external calcr
  integer, parameter :: covprt = 14
  integer, parameter :: covreq = 15
  integer, parameter :: d = 27
  integer d1
  integer, parameter :: dinit = 38
  integer dk
  integer, parameter :: dltfdj = 36
  integer, parameter :: dtype = 16
  real h
  real, parameter :: hfac = 1000.0E+00
  real, save :: hlim = 0.0E+00
  integer i
  integer iv(60+p)
  integer, parameter :: j = 33
  integer j1
  integer j1k
  integer k
  integer n
  integer nf
  integer, parameter :: nfcall = 6
  integer, parameter :: nfgcal = 7
  integer, parameter :: r = 50
  integer r1
  integer rn
  logical strted
  integer, parameter :: toobig = 2
  external ufparm
  integer uiparm(*)
  real urparm(*)
  real v(93 + n*p + 3*n + (p*(3*p+33))/2)
  real x(p)
  real xk

  d1 = 94 + 2 * n + ( p * ( 3 * p + 31 ) ) / 2
  iv(d) = d1
  r1 = d1 + p
  iv(r) = r1
  j1 = r1 + n
  iv(j) = j1
  rn = j1 - 1

  if ( iv(1) == 0 ) then
    call dfault ( iv, v )
  end if

  iv(covreq) = -abs ( iv(covreq) )
  if ( iv(covprt) /= 0 .and. iv(covreq) == 0 ) then
    iv(covreq) = -1
  end if

  strted = .true.

  if (iv(1) /= 12) go to 80

  strted = .false.
  iv(nfcall) = 1
  iv(nfgcal) = 1
!
!  Initialize scale vector D to ones for computing initial jacobian.
!
  if ( 0 < iv(dtype) ) then
    v(d1:d1+p-1) = 1.0E+00
  end if

10 continue

  nf = iv(nfcall)

  call calcr ( n, p, x, nf, v(r1), uiparm, urparm, ufparm )

  if ( strted ) then

    if ( nf <= 0 ) then
      iv(toobig) = 1
    end if

    go to 80

  end if

  if ( nf <= 0 ) then
    iv(1) = 13
    call itsmry ( v(d1), iv, p, v, x )
    return
  end if
!
!  Compute finite-difference jacobian.
!
30 continue

  j1k = j1
  dk = d1

  do k = 1, p

    xk = x(k)
    h = v(dltfdj) * max ( abs ( xk ), 1.0E+00 / v(dk) )
    dk = dk + 1

    do

      x(k) = xk + h
      nf = iv(nfgcal)
      call calcr ( n, p, x, nf, v(j1k), uiparm, urparm, ufparm )

      if ( 0 < nf ) then
        exit
      end if

      if ( hlim == 0.0E+00 ) then
        hlim = hfac * epsilon ( hlim )
      end if

      h = -0.5E+00 * h

      if ( abs ( h ) < hlim ) then
        iv(1) = 15
        call itsmry ( v(d1), iv, p, v, x )
        return
      end if

    end do

    x(k) = xk

    do i = r1, rn
      v(j1k) = ( v(j1k) - v(i) ) / h
      j1k = j1k + 1
    end do

  end do

  strted = .true.

80 continue

  call nl2itr ( v(d1), iv, v(j1), n, n, p, v(r1), v, x )

  if ( iv(1) < 2 ) then
    go to 10
  else if ( iv(1) == 2 ) then
    go to 30
  end if

  return
end
subroutine nl2itr ( d, iv, j, n, nn, p, r, v, x )

!*****************************************************************************80
!
!! NL2ITR carries out iterations for NL2SOL.
!
!  Discussion:
!
!    Parameters IV, N, P, V, and X are the same as the corresponding
!    ones to NL2SOL, except that V can be shorter, since the part of V
!    that NL2SOL uses for storing D, J, and R is not needed.
!
!    Moreover, compared with NL2SOL, IV(1) may have the
!    two additional output values 1 and 2, which are explained below,
!    as is the use of IV(TOOBIG) and IV(NFGCAL).  The values IV(D),
!    IV(J), and IV(R), which are output values from NL2SOL (and
!    NL2SNO), are not referenced by NL2ITR or the subroutines it calls.
!
!    On a fresh start, that is, a call on NL2ITR with IV(1) = 0 or 12,
!    NL2ITR assumes that R = R(X), the residual at X, and J = J(X),
!    the corresponding jacobian matrix of R at X.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
! iv(1) = 1 means the caller should set r to r(x), the residual at x,
!             and call nl2itr again, having changed none of the other
!             parameters.  an exception occurs if r cannot be evaluated
!             at x (e.g. if r would overflow), which may happen because
!             of an oversized step.  in this case the caller should set
!             iv(toobig) = iv(2) to 1, which will cause nl2itr to ig-
!             nore r and try a smaller step.  the parameter nf that
!             nl2sol passes to CALCR (for possible use by calcj) is a
!             copy of iv(nfcall) = iv(6).
! iv(1) = 2 means the caller should set j to j(x), the jacobian matrix
!             of r at x, and call nl2itr again.  the caller may change
!             d at this time, but should not change any of the other
!             parameters.  the parameter nf that nl2sol passes to
!             calcj is iv(nfgcal) = iv(7).  if j cannot be evaluated
!             at x, then the caller may set iv(nfgcal) to 0, in which
!             case nl2itr will return with iv(1) = 15.
!
!  Parameters:
!
!    Input, real D(N), the scale vector.
!
!    Input/output, integer IV(*), the NL2SOL integer parameter array.
!
!    j   n by p jacobian matrix (lead dimension nn).
!
!    n   number of observations (components in r).
!
!    nn  lead dimension of j.
!
!    p   number of parameters (components in x).
!
!    r   residual vector.
!
!    v   floating-point value array.
!
!    x   parameter vector.
!
  implicit none

  integer n
  integer nn
  integer p

  real d(p)
  integer dummy
  integer dig1
  real dotprd
  real e
  integer g01
  integer g1
  integer h0
  integer h1
  integer i
  integer ipiv1
  integer ipivi
  integer ipivk
  integer ipk
  integer iv(60+p)
  real j(nn,p)
  integer k
  integer km1
  integer l
  integer lky1
  integer lmat1
  integer lstgst
  integer m
  integer pp1o2
  integer qtr1
  real r(n)
  integer rd0
  integer rd1
  integer rdk
  real rdof1
  integer rsave1
  integer s1
  integer smh
  integer sstep
  integer step1
  logical stopx
  integer stpmod
  integer, parameter :: stppar = 5
  real sttsst
  real t
  real t1
  integer temp1
  integer temp2
  real v(93 + 2*n + (p*(3*p+31))/2)
  real v2norm
  integer w1
  real x(p)
  integer x01
!
!  external functions and subroutines
!
      external assess, covclc, dotprd, dupdat, gqtstp, itsmry, lmstep, &
               parchk, qapply, rptmul, slupdt, slvmul, stopx, &
               v2norm!
! subscripts for iv and v
!
      integer cnvcod, cosmin, covmat, covprt, covreq, dgnorm, dig, &
              dinit, dstnrm, dtype, d0init, f, fdif, fuzz, &
              f0, g, gtstep, h, ierr, incfac, inits, ipivot, ipiv0, irc, &
              jtinit, jtol1, kagqt, kalm, lky, lmat, lmax0, mode, model, &
              mxfcal, mxiter, nfcall, nfgcal, nfcov, ngcov, ngcall, &
              niter, nvsave, phmxfc, preduc, qtr, radfac, radinc, &
              radius, rad0, rd, restor, rlimit, rsave, s, size, step, &
              stglim, stlstg, sused, switch, toobig, tuner4, &
              tuner5, vsave1, w, wscale, xirc, x0
!
! iv subscript values
!
parameter ( cnvcod=34, covmat=26, covprt=14 )
parameter (covreq=15, dig=43, dtype=16, g=28, h=44 )
parameter (ierr=32, inits=25, ipivot=61, ipiv0=60 )
parameter (irc=3, kagqt=35, kalm=36, lky=37, lmat=58 )
parameter (mode=38, model=5, mxfcal=17, mxiter=18 )
parameter (nfcall=6, nfgcal=7, nfcov=40, ngcov=41 )
parameter (ngcall=30, niter=31, qtr=49 )
parameter (radinc=8, rd=51, restor=9, rsave=52, s=53 )
parameter (step=55, stglim=11, stlstg=56, sused=57 )
parameter (switch=12, toobig=2, w=59, xirc=13, x0=60)
!
! v subscript values.
!
     parameter (cosmin=43, dgnorm=1, dinit=38, dstnrm=2)
     parameter ( d0init=37, f=10, fdif=11, fuzz=45 )
    parameter ( f0=13, gtstep=4, incfac=23 )
    parameter ( jtinit=39, jtol1=87, lmax0=35 )
    parameter ( nvsave=9, phmxfc=21, preduc=7 )
    parameter ( radfac=16, radius=8, rad0=9, rlimit=42 )
    parameter ( size=47, tuner4=29, tuner5=30 )
    parameter ( vsave1=78, wscale=48)

      i = iv(1)
      if (i == 1) go to 20
      if (i == 2) go to 50
!
!  Check validity of iv and v input values.
!
!  If iv(1) = 0, then PARCHK calls dfault(iv, v).
!
      call parchk ( iv, n, nn, p, v )
      i = iv(1) - 2

      if ( 10 < i ) then
        return
      end if

      go to (350, 350, 350, 350, 350, 350, 195, 160, 195, 10), i
!
!  Initialization and storage allocation.
!
 10   continue

      iv(niter) = 0
      iv(nfcall) = 1
      iv(ngcall) = 1
      iv(nfgcal) = 1
      iv(mode) = -1
      iv(stglim) = 2
      iv(toobig) = 0
      iv(cnvcod) = 0
      iv(covmat) = 0
      iv(nfcov) = 0
      iv(ngcov) = 0
      iv(kalm) = -1
      iv(radinc) = 0
      iv(s) = jtol1 + 2*p
      pp1o2 = p * (p + 1) / 2
      iv(x0) = iv(s) + pp1o2
      iv(step) = iv(x0) + p
      iv(stlstg) = iv(step) + p
      iv(dig) = iv(stlstg) + p
      iv(g) = iv(dig) + p
      iv(lky) = iv(g) + p
      iv(rd) = iv(lky) + p
      iv(rsave) = iv(rd) + p
      iv(qtr) = iv(rsave) + n
      iv(h) = iv(qtr) + n
      iv(w) = iv(h) + pp1o2
      iv(lmat) = iv(w) + 4*p + 7
!
!  Length of w = p*(p+9)/2 + 7.  lmat is contained in w.
!
      if (v(dinit) >= 0.0E+00 ) then
        d(1:p) = v(dinit)
      end if

      if (v(jtinit) > 0.0E+00 ) then
        v(jtol1:jtol1+p-1) = v(jtinit)
      end if

      i = jtol1 + p

      if (v(d0init) > 0.0E+00 ) then
        v(i:i+p-1) = v(d0init)
      end if

      v(rad0) = 0.0E+00
      v(stppar) = 0.0E+00
      v(radius) = v(lmax0) / ( 1.0E+00 + v(phmxfc) )
!
!  Set initial model and S matrix.
!
      iv(model) = 1
      if (iv(inits) == 2) iv(model) = 2
      s1 = iv(s)
      if (iv(inits) == 0) then
        v(s1:s1+pp1o2-1) = 0.0E+00
      end if
!
!  Compute function value (half the sum of squares).
!
 20   continue

      t = v2norm(n, r)

      if ( v(rlimit) < t ) then
        iv(toobig) = 1
      end if

      if ( iv(toobig) == 0 ) then
        v(f) = 0.5E+00 * t**2
      end if

 30   continue

      if ( iv(mode) == 0 ) then
        go to 350
      end if

      if ( 0 < iv(mode) ) then
        go to 730
      end if

 40   continue

      if ( iv(toobig) /= 0 ) then
         iv(1) = 13
         call itsmry ( d, iv, p, v, x )
         return
      end if

      go to 60
!
!  Make sure jacobian could be computed.
!
 50   continue

      if ( iv(nfgcal) == 0 ) then
         iv(1) = 15
         call itsmry ( d, iv, p, v, x )
         return
      end if
!
!  Compute gradient.
!
 60   continue

      iv(kalm) = -1
      g1 = iv(g)
      do i = 1, p
         v(g1) = dot_product ( r(1:n), j(1:n,i) )
         g1 = g1 + 1
      end do

      if ( 0 < iv(mode) ) then
        go to 710
      end if
!
!  Update D and make copies of R for possible use later.
!
      if ( 0 < iv(dtype) ) then
        call dupdat(d, iv, j, n, nn, p, v)
      end if

      rsave1 = iv(rsave)
      v(rsave1:rsave1+n-1) = r(1:n)
      qtr1 = iv(qtr)
      v(qtr1:qtr1+n-1) = r(1:n)
!
!  Compute inverse ( D ) * gradient.
!
      g1 = iv(g)
      dig1 = iv(dig)
      k = dig1

      do i = 1, p
         v(k) = v(g1) / d(i)
         k = k + 1
         g1 = g1 + 1
      end do

      v(dgnorm) = v2norm(p, v(dig1))

      if (iv(cnvcod) /= 0) go to 700
      if (iv(mode) == 0) go to 570
      iv(mode) = 0
!
!  Main loop.
!
!  Print iteration summary, check iteration limit.
!
 150  continue

      call itsmry(d, iv, p, v, x)

 160  k = iv(niter)

      if ( iv(mxiter) <= k ) then
         iv(1) = 10
         call itsmry ( d, iv, p, v, x )
         return
      end if

170   continue

      iv(niter) = k + 1
!
!  Update radius.
!
      if ( k /= 0 ) then

        step1 = iv(step)
        do i = 1, p
          v(step1) = d(i) * v(step1)
          step1 = step1 + 1
        end do
        step1 = iv(step)
        v(radius) = v(radfac) * v2norm(p, v(step1))

      end if
!
!  Initialize for start of next iteration.
!
      x01 = iv(x0)
      v(f0) = v(f)
      iv(kagqt) = -1
      iv(irc) = 4
      iv(h) = -abs ( iv(h) )
      iv(sused) = iv(model)
!
!  Copy X to X0.
!
      v(x01:x01+p-1) = x(1:p)
!
!  Check STOPX and function evaluation limit.
!
 190  if ( .not. stopx ( dummy ) ) go to 200
         iv(1) = 11
         go to 205
!
!  Come here when restarting after function evaluation limit or STOPX.
!
 195  continue

      if ( v(f) < v(f0) ) then
         v(radfac) = 1.0E+00
         k = iv(niter)
         go to 170
      end if

 200  continue

      if (iv(nfcall) < iv(mxfcal) + iv(nfcov)) go to 210
         iv(1) = 9
 205     continue

         if (v(f) >= v(f0)) then
           call itsmry ( d, iv, p, v, x )
           return
         end if
!
!  In case of STOPX or function evaluation limit with
!  improved V(F), evaluate the gradient at X.
!
         iv(cnvcod) = iv(1)
         go to 560
!
!  Compute candidate step.
!
 210  continue

      step1 = iv(step)
      w1 = iv(w)
!
!  Compute Levenberg-Marquardt step.
!
      if (iv(model) /= 2) then

         qtr1 = iv(qtr)

         if ( iv(kalm) < 0 ) then

           rd1 = iv(rd)

           if (-1 == iv(kalm)) then
             call qrfact ( nn, n, p, j, v(rd1), &
             iv(ipivot), iv(ierr), 0, v(w1) )
           end if

           call qapply ( nn, n, p, j, v(qtr1), iv(ierr) )

         end if

         h1 = iv(h)
!
!  Copy R matrix to H.
!
         if ( h1 <= 0 ) then

              h1 = -h1
              iv(h) = h1
              k = h1
              rd1 = iv(rd)
              v(k) = v(rd1)

              do i = 2, p
                   call vcopy(i-1, v(k+1), j(1,i))
                   k = k + i
                   rd1 = rd1 + 1
                   v(k) = v(rd1)
              end do
         end if

         g1 = iv(g)
         call lmstep(d, v(g1), iv(ierr), iv(ipivot), iv(kalm), p, &
                     v(qtr1), v(h1), v(step1), v, v(w1))
!
!  Compute Goldfeld-Quandt-Trotter step (augmented model).
!
      else

      if ( iv(h) <= 0 ) then
!
!  Set H to inverse ( D ) * ( J' * J + s) ) * inverse ( D ).
!
        h1 = -iv(h)
        iv(h) = h1
        s1 = iv(s)
!
!  J is in its original form.
!
        if ( iv(kalm) == -1 ) then

          do i = 1, p
            t = 1.0E+00 / d(i)
            do k = 1, i
              v(h1) = t * (dotprd(n,j(1,i),j(1,k))+v(s1)) / d(k)
              h1 = h1 + 1
              s1 = s1 + 1
            end do
          end do
!
!  LMSTEP has applied QRFACT to J.
!
        else

          smh = s1 - h1
          h0 = h1 - 1
          ipiv1 = iv(ipivot)
          t1 = 1.0E+00 / d(ipiv1)
          rd0 = iv(rd) - 1
          rdof1 = v(rd0 + 1)

          do i = 1, p

            l = ipiv0 + i
            ipivi = iv(l)
            h1 = h0 + ( ipivi*(ipivi-1) ) / 2
            l = h1 + ipivi
            m = l + smh
!
!  v(l) = h(ipivot(i), ipivot(i))
!  v(m) = s(ipivot(i), ipivot(i))
!
            t = 1.0E+00 / d(ipivi)
            rdk = rd0 + i
            e = v(rdk)**2
            if ( 1 < i ) then
              e = e + dotprd(i-1, j(1,i), j(1,i))
            end if
            v(l) = (e + v(m)) * t**2

            if ( i /= 1 ) then

              l = h1 + ipiv1
              if (ipivi < ipiv1) then
                l = l + ((ipiv1-ipivi)*(ipiv1+ipivi-3)) / 2
              end if
              m = l + smh
!
!  v(l) = h(ipivot(i), ipivot(1))
!  v(m) = s(ipivot(i), ipivot(1))
!
              v(l) = t * (rdof1 * j(1,i)  +  v(m)) * t1

              do k = 2, i - 1
                ipk = ipiv0 + k
                ipivk = iv(ipk)
                l = h1 + ipivk
                if (ipivi < ipivk) then
                  l = l + ((ipivk-ipivi)*(ipivk+ipivi-3)) / 2
                end if
                m = l + smh
!
!  v(l) = h(ipivot(i), ipivot(k))
!  v(m) = s(ipivot(i), ipivot(k))
!
                km1 = k - 1
                rdk = rd0 + k
                v(l) = t * (dotprd(km1, j(1,i), j(1,k)) + &
                  v(rdk)*j(k,i) + v(m)) / d(ipivk)
              end do

            end if

          end do

        end if

      end if
!
!  Compute actual Goldfeld-Quandt-Trotter step.
!
      h1 = iv(h)
      dig1 = iv(dig)
      lmat1 = iv(lmat)
      call gqtstp(d, v(dig1), v(h1), iv(kagqt), v(lmat1), p, v(step1), &
                  v, v(w1))
    end if
!
!  Compute R(X0 + STEP).
!
 310  continue

      if ( iv(irc) /= 6 ) then
        x01 = iv(x0)
        step1 = iv(step)
        x(1:p) = v(step1:step1+p-1) + v(x01:x01+p-1)
        iv(nfcall) = iv(nfcall) + 1
        iv(1) = 1
        iv(toobig) = 0
        return
      end if
!
!  Assess candidate step.
!
350   continue

      step1 = iv(step)
      lstgst = iv(stlstg)
      x01 = iv(x0)
      call assess(d, iv, p, v(step1), v(lstgst), v, x, v(x01))
!
!  If necessary, switch models and/or restore R.
!
      if ( iv(switch) /= 0 ) then
        iv(h) = -abs ( iv(h) )
        iv(sused) = iv(sused) + 2
        v(1:nvsave) = v(vsave1:vsave1+nvsave-1)
      end if

 360  continue

      if ( iv(restor) /= 0 ) then
         rsave1 = iv(rsave)
         r(1:n) = v(rsave1:rsave1+n-1)
      end if

 390  continue

      l = iv(irc) - 4
      stpmod = iv(model)

      if (l > 0) go to (410,440,450,450,450,450,450,450,640,570), l
!
!  Decide whether to change models.
!
      e = v(preduc) - v(fdif)
      sstep = iv(lky)
      s1 = iv(s)
      call slvmul ( p, v(sstep), v(s1), v(step1) )
      sttsst = 0.5E+00 * dotprd(p, v(step1), v(sstep))

      if ( iv(model) == 1 ) then
        sttsst = -sttsst
      end if
!
!  Switch models.
!
      if (abs(e + sttsst) * v(fuzz) >= abs(e)) go to 400

         iv(model) = 3 - iv(model)
         if (iv(model) == 1) iv(kagqt) = -1
         if (iv(model) == 2 .and. iv(kalm) > 0) iv(kalm) = 0
         if (-2 < l) go to 480
              iv(h) = -abs ( iv(h) )
              iv(sused) = iv(sused) + 2
              v(vsave1:vsave1+nvsave-1) = v(1:nvsave)
              go to 420

 400  if (-3 < l) go to 480
!
!  Recompute STEP with decreased radius.
!
         v(radius) = v(radfac) * v(dstnrm)
         go to 190
!
!  Recompute STEP, saving V values and R if necessary.
!
 410  continue

      v(radius) = v(radfac) * v(dstnrm)

 420  continue

      if ( v(f) < v(f0) ) then
        rsave1 = iv(rsave)
        v(rsave1:rsave1+n-1) = r(1:n)
      end if

      go to 190
!
!  Compute step of length V(LMAX0) for singular convergence test.
!
 440  continue

      v(radius) = v(lmax0)
      go to 210
!
!  Convergence or false convergence.
!
 450  continue

      iv(cnvcod) = l
      if (v(f) >= v(f0)) go to 700
         if (iv(xirc) == 14) go to 700
              iv(xirc) = 14
!
!  Process acceptable step.
!
 480  iv(covmat) = 0
!
!  Set LKY = J(X0)' * R(X).
!
      lky1 = iv(lky)
!
!  Jacobian has not been modified.
!
      if ( iv(kalm) < 0 ) then

         do i = 1, p
           v(lky1) = dotprd(n, j(1,i), r)
           lky1 = lky1 + 1
         end do
!
!  QRFACT has been applied to J.  Store copy of R in QTR and
!  apply Q to it.
!
      else

        qtr1 = iv(qtr)
        v(qtr1:qtr1+n-1) = r(1:n)
        call qapply(nn, n, p, j, v(qtr1), iv(ierr))
!
!  Multiply top P-vector in QTR by permuted upper triangle
!  stored by QRFACT in J and RD.
!
        rd1 = iv(rd)
        temp1 = iv(stlstg)
        call rptmul(3, iv(ipivot), j, nn, p, v(rd1), v(qtr1), v(lky1), &
                  v(temp1))

      end if
!
!  See whether to set V(RADFAC) by gradient tests.
!
 510  continue

      if (iv(irc) == 3 ) then

        step1 = iv(step)
        temp1 = iv(stlstg)
        temp2 = iv(x0)
!
!  Set TEMP1 = hessian * STEP for use in gradient tests
!
!  STEP computed using Gauss-Newton model.
!  QRFACT has been applied to J.
!
        if ( stpmod /= 2 ) then

          rd1 = iv(rd)
          call rptmul(2, iv(ipivot), j, nn, p, v(rd1), &
            v(step1), v(temp1), v(temp2))
!
!  STEP computed using augmented model.
!
        else

          h1 = iv(h)
          k = temp2

          do i = 1, p
            v(k) = d(i) * v(step1)
            k = k + 1
            step1 = step1 + 1
          end do

          call slvmul(p, v(temp1), v(h1), v(temp2))

          do i = 1, p
            v(temp1) = d(i) * v(temp1)
            temp1 = temp1 + 1
          end do

        end if

      end if
!
!  Save old gradient and compute new one.
!
 560  continue

      iv(ngcall) = iv(ngcall) + 1
      g1 = iv(g)
      g01 = iv(w)
      v(g01:g01+p-1) = v(g1:g1+p-1)
      iv(1) = 2
      return
!
!  Initializations -- g0 = g - g0, etc.
!
 570  continue

      g01 = iv(w)
      g1 = iv(g)
      v(g01:g01+p-1) = - v(g01:g01+p-1) + v(g1:g1+p-1)
      step1 = iv(step)
      temp1 = iv(stlstg)
      temp2 = iv(x0)
!
!  Set V(RADFAC) by gradient tests.
!
!  Set TEMP1 = d**-1 * (hessian * STEP  +  ( G(x0) - G(x) ) ).
!
      if ( iv(irc) == 3 ) then

         k = temp1
         l = g01
         do i = 1, p
           v(k) = (v(k) - v(l)) / d(i)
           k = k + 1
           l = l + 1
         end do
!
!  Do gradient tests.
!
         if ( v2norm(p, v(temp1)) <= v(dgnorm) * v(tuner4) .or. &
           dotprd(p, v(g1), v(step1)) < v(gtstep) * v(tuner5) ) then
           v(radfac) = v(incfac)
         end if

      end if
!
!  Finish computing LKY = ( J(X) - J(X0) )' * R.
!
!  Currently LKY = J(X0)' * R.
!
      lky1 = iv(lky)
      v(lky1:lky1+p-1) = - v(lky1:lky1+p-1) + v(g1:g1+p-1)
!
!  Determine sizing factor V(SIZE).
!
!  Set TEMP1 = S * STEP.
!
      s1 = iv(s)
      call slvmul(p, v(temp1), v(s1), v(step1))

      t1 = abs(dotprd(p, v(step1), v(temp1)))
      t = abs(dotprd(p, v(step1), v(lky1)))
      v(size) = 1.0E+00

      if ( t < t1 ) then
        v(size) = t / t1
      end if
!
!  Update S.
!
      call slupdt(v(s1), v(cosmin), p, v(size), v(step1), v(temp1), &
                  v(temp2), v(g01), v(wscale), v(lky1))
      iv(1) = 2
      go to 150
!
!  Bad parameters to ASSESS.
!
 640  iv(1) = 14
      call itsmry ( d, iv, p, v, x )
      return
!
!  Convergence obtained.  Compute covariance matrix if desired.
!
 700  continue

      if ( ( iv(covreq) == 0 .and. iv(covprt) == 0 ) .or. &
        iv(covmat) /= 0 .or. &
        iv(cnvcod) >= 7 ) then
        iv(1) = iv(cnvcod)
        iv(cnvcod) = 0
        call itsmry(d, iv, p, v, x)
        return
      end if

      iv(mode) = 0

 710  continue

      call covclc(i, d, iv, j, n, nn, p, r, v, x)

      if ( i == 3 ) then

        iv(ngcov) = iv(ngcov) + 1
        iv(ngcall) = iv(ngcall) + 1
        iv(1) = 2

      else if ( i == 4 ) then

        if ( iv(niter) == 0 ) then
          iv(mode) = -1
        else
          iv(mode) = 0
        end if

        iv(1) = iv(cnvcod)
        iv(cnvcod) = 0

        call itsmry(d, iv, p, v, x)

      else

        iv(nfcov) = iv(nfcov) + 1
        iv(nfcall) = iv(nfcall) + 1
        iv(restor) = i
        iv(1) = 1

      end if

      return

 730  continue

  if ( iv(restor) == 1 .or. iv(toobig) /= 0 ) then
    go to 710
  end if

  iv(nfgcal) = iv(nfcall)
  iv(ngcov) = iv(ngcov) + 1
  iv(ngcall) = iv(ngcall) + 1
  iv(1) = 2

  return
end
subroutine parchk ( iv, n, nn, p, v )

!*****************************************************************************80
!
!! PARCHK checks the NL2SOL parameters.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
  implicit none

  real, save :: big = 0.0E+00
  integer i
  integer iv(*)
  integer iv1
  integer jtolp
  integer k
  integer l
  integer m
  real machep
  integer n
  integer nn
  integer, parameter :: nvdflt = 27
  integer p
  integer, parameter :: parsv1 = 51
  integer, parameter :: prunit = 21
  integer pu
  real, save :: teensy = 1.0E+00
  real v(*)
  real vk
  real vm(27)
  character*8 vn(27)
  real vx(27)

     character*4 cngd(3), dflt(3), which(3)
!
! iv and v subscripts
!
      integer dtype, dtype0, d0init, epslon, inits, jtinit, jtol0, &
              jtol1, oldn, oldnn, oldp, parprt

     parameter (dtype=16, dtype0=29, d0init=37, epslon=19 )
     parameter ( inits=25, jtinit=39, jtol0=86, jtol1=87 )
     parameter ( oldn=45, oldnn=46, oldp=47, parprt=20 )

  data vn / &
    'epslon..', 'phmnfc..', 'phmxfc..', 'decfac..', 'incfac..', &
    'rdfcmn..', 'rdfcmx..', 'tuner1..', 'tuner2..', 'tuner3..', &
    'tuner4..', 'tuner5..', 'afctol..', 'rfctol..', 'xctol...', &
    'xftol...', 'lmax0...', 'dltfdj..', 'd0init..', 'dinit...', &
    'jtinit..', 'dltfdc..', 'dfac....', 'rlimit..', 'cosmin..', &
    'delta0..', 'fuzz....' /

      data vm(1)/1.0e-3/, vm(2)/-0.99e+0/, vm(3)/1.0e-3/, vm(4)/1.0e-2/, &
           vm(5)/1.2e+0/, vm(6)/1.e-2/, vm(7)/1.2e+0/, vm(8)/0.e+0/, &
           vm(9)/0.e+0/, vm(10)/1.e-3/, vm(11)/-1.e+0/, vm(15)/0.e+0/, &
           vm(16)/0.e+0/, vm(19)/0.e+0/, vm(20)/-10.e+0/, vm(21)/0.e+0/, &
           vm(23)/0.e+0/, vm(24)/1.e+10/, vm(27)/1.01e+0/
      data vx(1)/0.9e+0/, vx(2)/-1.e-3/, vx(3)/1.e+1/, vx(4)/0.8e+0/, &
           vx(5)/1.e+2/, vx(6)/0.8e+0/, vx(7)/1.e+2/, vx(8)/0.5e+0/, &
           vx(9)/0.5e+0/, vx(10)/1.e+0/, vx(11)/1.e+0/, vx(14)/0.1e+0/, &
           vx(15)/1.e+0/, vx(16)/1.e+0/, vx(18)/1.e+0/, vx(22)/1.e+0/, &
           vx(23)/1.e+0/, vx(25)/1.e+0/, vx(26)/1.e+0/, vx(27)/1.e+2/

     data cngd(1),cngd(2),cngd(3)/'---c','hang','ed v'/
     data dflt(1),dflt(2),dflt(3)/'nond','efau','lt v'/

  if ( iv(1) == 0 ) then
    call dfault ( iv, v )
  end if

  pu = iv(prunit)
  iv1 = iv(1)

  if ( iv1 == 12) then

    if ( nn < n .or. n < p .or. p < 1) then
      iv(1) = 16
      if ( pu /= 0 ) then
        write ( pu, '(a)' ) ' '
        write ( pu, '(a)' ) '  Bad NN, N or P:'
        write ( pu, '(a)' ) ' '
        write ( pu, '(a,i5)' ) '  NN = ', nn
        write ( pu, '(a,i5)' ) '  N =  ', n
        write ( pu, '(a,i5)' ) '  P =  ', p
      end if
      return
    end if

    k = iv(21)
    call dfault(iv(21), v(33))
    iv(21) = k
    iv(dtype0) = iv(dtype+20)
    iv(oldn) = n
    iv(oldnn) = nn
    iv(oldp) = p
    which(1) = dflt(1)
    which(2) = dflt(2)
    which(3) = dflt(3)

  else

    if ( n  /= iv(oldn)  .or. &
         nn /= iv(oldnn) .or. &
         p  /= iv(oldp) ) then

      iv(1) = 17

      if ( pu /= 0 ) then
        write ( pu, '(a)' ) ' '
        write ( pu, '(a)' ) '(NN,N,P) changed from:'
        write ( pu, '(a,i8)' ) '  NN = ', iv(oldnn)
        write ( pu, '(a,i8)' ) '  N =  ', iv(oldn)
        write ( pu, '(a,i8)' ) '  P =  ', iv(oldp)
        write ( pu, '(a)' ) ' to:'
        write ( pu, '(a,i8)' ) '  NN = ', nn
        write ( pu, '(a,i8)' ) '  N =  ', n
        write ( pu, '(a,i8)' ) '  P =  ', p
      end if

      return

    end if

    if ( iv1 < 1 .or. 11 < iv1 ) then
      iv(1) = 50
      if (pu /= 0) then
        write(pu,60) iv1
      end if
 60   format('  iv(1) =', i5, ' should be between 0 and 12.')
      return
    end if

    which(1) = cngd(1)
    which(2) = cngd(2)
    which(3) = cngd(3)

  end if

  if ( big <= teensy ) then
    teensy = tiny ( teensy )
    machep = epsilon ( machep )
    big = huge ( big )
    vm(12) = machep
    vx(12) = big
    vm(13) = teensy
    vx(13) = big
    vm(14) = machep
    vm(17) = teensy
    vx(17) = big
    vm(18) = machep
    vx(19) = big
    vx(20) = big
    vx(21) = big
    vm(22) = machep
    vx(24) = sqrt ( 0.999E+00 * huge ( vx(24) ) )
    vm(25) = machep
    vm(26) = machep
  end if

  m = 0

  if (iv(inits) >= 0 .and. iv(inits) <= 2) go to 110
         m = 18
         if (pu /= 0) write(pu,100) iv(inits)
 100     format('inits... iv(25) =',i4,' should be between 0 and 2.')
 110 continue

  k = epslon

  do i = 1, nvdflt
    vk = v(k)
    if (vk >= vm(i) .and. vk <= vx(i)) go to 130
      m = k
      if (pu /= 0) write(pu,120) vn(i), k, vk, vm(i), vx(i)
120   format( a8, '.. v(',i2, ') =', e11.3, &
      ' should be between ',e11.3, ' and', d11.3 )
130 continue
    k = k + 1

  end do
!
!  Check JTOL values.
!
  if ( iv1 /= 12 .or. v(jtinit) <= 0.0E+00) then

    jtolp = jtol0 + p
    do i = jtol1, jtolp
      if ( v(i) <= 0.0E+00 ) then
        k = i - jtol0
        if (pu /= 0) write(pu,150) k, i, v(i)
 150    format( 'jtol(', i3, ') = v(', i3, ') =', e11.3, &
                ' should be positive.' )
        m = i
      end if
    end do

  end if

  if ( m /= 0 ) then
    iv(1) = m
    return
  end if

 180  continue

  if ( pu == 0 .or. iv(parprt) == 0 ) then
    return
  end if

  if ( iv1 == 12 .and. iv(inits) /= 0) then
    m = 1
    write(pu,190) iv(inits)
190 format( 'nondefault values....inits..... iv(25) =', i3)
  end if

  if ( iv(dtype) /= iv(dtype0) ) then
    if (m == 0) write(pu,215) which
    m = 1
    write ( pu, '(a,i3)' ) 'DTYPE..... IV(16) = ', iv(dtype)
  end if

  k = epslon
  l = parsv1

  do i = 1, nvdflt

    if ( v(k) /= v(l) ) then
      if (m == 0) write(pu,215) which
 215  format(3a4,'alues....')
      m = 1
      write(pu,220) vn(i), k, v(k)
 220  format(1x,a8,'.. v(',i2,') =',e15.7)
    end if

    k = k + 1
    l = l + 1

  end do

  iv(dtype0) = iv(dtype)
  v(parsv1:parsv1+nvdflt-1) = v(epslon:epslon+nvdflt-1)

  if ( iv1 /= 12 ) then
    return
  end if

  if ( v(jtinit) <= 0.0E+00 ) then
    write ( pu, '(a)' ) '(Initial) JTOL array'
    write ( pu, '(6e12.3)' ) v(jtol1:jtol0+p)
  end if

  if ( v(d0init) <= 0.0E+00 ) then
    k = jtol1 + p
    write ( pu, '(a)' ) '(Initial) D0 array'
    write ( pu, '(6e12.3)' ) v(k:k+p-1)
  end if

  return
end
subroutine qapply ( nn, n, p, j, r, ierr )

!*****************************************************************************80
!
!! QAPPLY applies orthogonal transformation to the residual R.
!
!  Discussion:
!
!    This subroutine applies to R the orthogonal transformations
!    stored in J by QRFACT.
!
!    The vectors U which determine the Householder transformations
!    are normalized so that their 2-norm squared is 2.  The use of
!    these transformations here is in the spirit of Businger and Golub.
!
!  Modified:
!
!    06 April 2006
!
!  Author:
!
!    David Gay
!
!  Reference:
!
!    P A Businger and Gene Golub,
!    Linear Least Squares Solutions by Householder Transformations,
!    Numerische Mathematik,
!    Volume 7, pages 269-276, 1965.
!
!  Parameters:
!
!    Input, integer NN, the row dimension of the matrix J as declared in
!    the calling program dimension statement
!
!    Input, integer N, the number of rows of J and the size of the R.
!
!    Input, integer P, the number of columns of J and the size of SIGMA.
!
!    Input, real J(NN,P), an N by P matrix.  It contains on its diagonal
!    and below its diagonal the column vectors U which determine the
!    Householder transformations (identity - U*U').
!
!    Input/output, real R(N).  On input, the right hand side vector to
!    which the orthogonal transformations will be applied.  On output,
!    R has been transformed.
!
!    Input, integer IERR, if non-zero, indicates that not all the
!    transformations were successfully determined and only the first
!    abs(IERR) - 1 transformations will be used.
!
  implicit none

  integer n
  integer nn
  integer p

  real dotprd
  integer i
  integer ierr
  real j(nn,p)
  integer k
  integer l
  integer nl1
  real r(n)
  real t

  if ( ierr /= 0 ) then
    k = abs(ierr) - 1
  else
    k = p
  end if

  do l = 1, k

    nl1 = n - l + 1
    t = -dotprd ( nl1, j(l,l), r(l) )
    r(l:n) = r(l:n) + t * j(l:n,l)

  end do

  return
end
subroutine qrfact ( nm, m, n, qr, alpha, ipivot, ierr, nopivk, sum )

!*****************************************************************************80
!
!! QRFACT computes the QR decomposition of a matrix.
!
!  Discussion:
!
!    This subroutine does a QR decomposition on the M x N matrix QR,
!    with an optionally modified column pivoting, and returns the
!    upper triangular R matrix, as well as the orthogonal vectors
!    used in the transformations.
!
!    This may be used when solving linear least-squares problems.
!    See subroutine QR1 of ROSEPACK.  It is called for this purpose
!    in the NL2SOL package.
!
!    This version of QRFACT tries to eliminate the occurrence of
!    underflows during the accumulation of inner products.  RKTOL1
!    is chosen below so as to insure that discarded terms have no
!    effect on the computed two-norms.
!
!    This routine was adapted from Businger and Golub's ALGOL
!    routine "SOLVE".
!
!    This routine is equivalent to the subroutine QR1 of ROSEPACK
!    with RKTOL1 used in place of RKTOL below, with V2NORM used
!    to initialize (and sometimes update) the sum array, and
!    with calls on DOTPRD in place of some loops.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Reference:
!
!    P Businger and Gene Golub,
!    Linear Least Squares Solutions by Householder Transformations,
!    in Handbook for Automatic Computation,
!    Volume II, Linear Algebra,
!    edited by James Wilkinson and C Reinsch,
!    Springer Verlag, pages 111-118, 1971;
!    prepublished in Numerische Mathematik,
!    Volume 7, pages 269-276, 1965.
!
!  Parameters:
!
!    Input, integer NM, the row dimension of the two dimensional
!    array parameters as declared in the calling program dimension statement.
!
!    Input, integer M, the number of rows in the matrix.
!
!    Input, integer N, the number of columns in the matrix.
!
!    Input/output, real QR(NM,N), on input, the M by N rectangular matrix
!    to be decomposed.  On output, contains the non-diagonal elements of
!    the R matrix in the strict upper triangle.  The vectors U, which
!    define the Householder transformations (Identity - U*U'), are in the
!    columns of the lower triangle.  These vectors U are scaled so that
!    the square of their 2-norm is 2.0.
!
!    Output, real ALPHA(N), the diagonal elements of R.
!
!    Output, integer IPIVOT(N), reflects the column pivoting performed
!    on the input matrix to accomplish the decomposition.  The J-th
!    element of IPIVOT gives the column of the original matrix which was
!    pivoted into column J during the decomposition.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    K, if no non-zero pivot could be found for the K-th transformation,
!    -K, for an error exit on the K-th transformation.
!    If an error exit was taken, the first (K - 1) transformations are correct.
!
!    Input, integer NOPIVK, controls pivoting.  Columns 1 through NOPIVK
!    will remain fixed in position.
!
!    Workspace, real SUM(N).
!
!  Local Parameters:
!
!    Local, real UFETA, the smallest positive floating point number
!    such that UFETA and -UFETA can both be represented.
!
!    Local, real RKTOL, the square root of the relative precision
!    of floating point arithmetic (MACHEP).
!
  implicit none

  integer n
  integer nm

  real alpha(n)
  real alphak
  real beta
  real dotprd
  integer i
  integer ierr
  integer ipivot(n)
  integer j
  integer jbar
  integer k
  integer m
  integer minum
  integer mk1
  integer nopivk
  real qr(nm,n)
  real qrkk
  real qrkmax
  real, save :: rktol = 0.0E+00
  real rktol1
  real sigma
  real sum(n)
  real sumj
  real temp
  real, save :: ufeta = 0.0E+00
  real v2norm

  if ( ufeta <= 0.0E+00 ) then
    ufeta = tiny ( ufeta )
    rktol = sqrt ( 0.999E+00 * epsilon ( rktol ) )
  end if

  ierr = 0
  rktol1 = 0.01E+00 * rktol

  do j = 1, n
    sum(j) = v2norm ( m, qr(1,j) )
    ipivot(j) = j
  end do

  minum = min ( m, n )

  do k = 1, minum

    mk1 = m - k + 1
!
!  K-th Householder transformation.
!
    sigma = 0.0E+00
    jbar = 0
!
!  Find largest column sum.
!
    if ( nopivk < k ) then

      do j = k, n
        if ( sigma < sum(j) ) then
          sigma = sum(j)
          jbar = j
        end if
      end do

      if ( jbar == 0 ) then
        ierr = k
        do i = k, n
          alpha(i) = 0.0E+00
          if ( k < i ) then
            qr(k:i-1,i) = 0.0E+00
          end if
        end do
        return
      end if
!
!  Column interchange.
!
      i = ipivot(k)
      ipivot(k) = ipivot(jbar)
      ipivot(jbar) = i

      sum(jbar) = sum(k)
      sum(k) = sigma

      do i = 1, m
        sigma = qr(i,k)
        qr(i,k) = qr(i,jbar)
        qr(i,jbar) = sigma
      end do

    end if
!
!  Second inner product.
!
    qrkmax = maxval ( abs ( qr(k:m,k) ) )

    if ( qrkmax < ufeta ) then
      ierr = -k
      do i = k, n
        alpha(i) = 0.0E+00
        if ( k < i ) then
          qr(k:i-1,i) = 0.0E+00
        end if
      end do
      return
    end if

    alphak = v2norm ( mk1, qr(k,k) ) / qrkmax
    sigma = alphak**2
!
!  End second inner product.
!
    qrkk = qr(k,k)
    if ( 0.0E+00 <= qrkk ) then
      alphak = -alphak
    end if

    alpha(k) = alphak * qrkmax
    beta = qrkmax * sqrt ( sigma - ( qrkk * alphak / qrkmax ) )
    qr(k,k) = qrkk - alpha(k)
    qr(k:m,k) =  qr(k:m,k) / beta

    do j = k + 1, n

      temp = -dotprd ( mk1, qr(k,k), qr(k,j) )

      qr(k:m,j) = qr(k:m,j) + temp * qr(k:m,k)

      if ( k + 1 <= m ) then
        sumj = sum(j)
        if ( ufeta <= sumj ) then
          temp = abs ( qr(k,j) / sumj )
          if ( rktol1 <= temp ) then
            if ( 0.99E+00 <= temp ) then
              sum(j) = v2norm ( m-k, qr(k+1,j) )
            else
              sum(j) = sumj * sqrt ( 1.0E+00 - temp**2 )
            end if

          end if
        end if
      end if

    end do

  end do

  return
end
function reldst ( p, d, x, x0 )

!*****************************************************************************80
!
!! RELDST computes the relative difference between two real values.
!
!  Modified:
!
!    03 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input, integer P, the length of the vectors.
!
!    Input, real D(P), a scaling vector.
!
!    Input, real X(P), X0(P), two vectors whose relative difference
!    is desired.
!
!    Output, real RELDST, the relative difference between X and X0.
!
  implicit none

  integer p

  real d(p)
  real emax
  integer i
  real reldst
  real x(p)
  real x0(p)
  real xmax

  emax = 0.0E+00
  xmax = 0.0E+00
  do i = 1, p
    emax = max ( emax, abs ( d(i) * ( x(i) - x0(i) ) ) )
    xmax = max ( xmax, d(i) * ( abs ( x(i) ) + abs ( x0(i) ) ) )
  end do

  if ( 0.0E+00 < xmax ) then
    reldst = emax / xmax
  else
    reldst = 0.0E+00
  end if

  return
end
subroutine rptmul ( func, ipivot, j, nn, p, rd, x, y, z )

!*****************************************************************************80
!
!! RPTMUL multiplies the R factor times a vector X.
!
!  Discussion:
!
!    This routine computes one of:
!
!      Y = R * P' * X
!      Y = P * R' * R * P' * X
!      Y = P * R' * X.
!
!    where P is a permutation matrix represented by a permutation
!    vector, and R is an upper triangular matrix, the R factor of
!    a QR factorization.
!
!    The strict upper triangle of R is stored in the strict upper triangle
!    of the array J, and the diagonal of R is stored in the vector RD.
!
!    X and Y may share storage.
!
!  Modified:
!
!    10 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input, integer FUNC, determines which product to compute:
!    1, Y =                RMAT * PERM' * X.
!    2, Y = PERM * PERM' * RMAT * PERM' * X.
!    3, Y = PERM *                PERM' * X.
!
!    Input, integer IPIVOT(P), the permutation vector.
!
!    Input, real J(NN,P), contains the strict upper triangle of the
!    matrix RMAT.
!
!    Input, integer NN, the leading dimension of J.
!
!    Input, integer P, the length of X and Y, and the order of RMAT.
!
!    Input, real RD(P), the diagonal elements of the matrix RMAT.
!
!    Input, real X(P), the input vector.
!
!    Output, real Y(P), the output vector.
!
!    Workspace, real Z(P).
!
  implicit none

  integer nn
  integer p

  real dotprd
  integer func
  integer i
  integer ipivot(p)
  real j(nn,p)
  integer k
  integer km1
  real rd(p)
  real x(p)
  real y(p)
  real z(p)
  real zk

  if ( func <= 2 ) then
!
!  Set Z = PERM' * X.
!
    do i = 1, p
      k = ipivot(i)
      z(i) = x(k)
    end do
!
!  Set Y = RMAT * Z.
!
    y(1) = z(1) * rd(1)

    do k = 2, p
      zk = z(k)
      do i = 1, k-1
        y(i) = y(i) + j(i,k) * zk
      end do
      y(k) = zk * rd(k)
    end do

    if ( func <= 1 ) then
      return
    end if

  else

    y(1:p) = x(1:p)

  end if
!
!  Set Z = RMAT' * Y.
!
  z(1) = y(1) * rd(1)

  do i = 2, p
    z(i) = y(i) * rd(i) + dotprd ( i-1, j(1,i), y )
  end do
!
!  Set Y = PERM * Z.
!
  do i = 1, p
    k = ipivot(i)
    y(k) = z(i)
  end do

  return
end
subroutine slupdt ( a, cosmin, p, size, step, u, w, wchmtd, wscale, y )

!*****************************************************************************80
!
!! SLUPDT updates a symmetric matrix A so that A * STEP = Y.
!
!  Discussion:
!
!    Update the symmetric matrix A so that A * STEP = Y.  Only the lower
!    triangle of A is stored, by rows.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
  implicit none

  integer p

  real a((p*(p+1))/2)
  real cosmin
  real denmin
  real dotprd
  integer i
  integer j
  integer k
  real sdotwm
  real size
  real step(p)
  real t
  real u(p)
  real v2norm
  real w(p)
  real wchmtd(p)
  real wscale
  real y(p)

  sdotwm = dot_product ( step(1:p), wchmtd(1:p) )

  denmin = cosmin * v2norm ( p, step ) * v2norm ( p, wchmtd )

  if ( denmin /= 0.0E+00 ) then
    wscale = min ( 1.0E+00, abs ( sdotwm / denmin ) )
  else
    wscale = 1.0E+00
  end if

  if ( sdotwm /= 0.0E+00 ) then
    t = wscale / sdotwm
  else
    t = 0.0E+00
  end if

  w(1:p) = t * wchmtd(1:p)

  call slvmul ( p, u, a, step )

  t = 0.5E+00 * ( size * dotprd ( p, step, u ) - dotprd ( p, step, y ) )

  u(1:p) = t * w(1:p) + y(1:p) - size * u(1:p)
!
!  Set A = A + U * W' + W * U'.
!
  k = 1
  do i = 1, p
    do j = 1, i
      a(k) = size * a(k) + u(i) * w(j) + w(i) * u(j)
      k = k + 1
    end do
  end do

  return
end
subroutine slvmul ( p, y, s, x )

!*****************************************************************************80
!
!! SLVMUL sets Y = S * X, where S is a P by P symmetric matrix.
!
!  Discussion:
!
!    This routine sets Y = S * X,  where X is a given vector and
!    S is a P by P symmetric matrix.  The lower triangle of S is
!    stored by rows.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input, integer P, the order of S, X and Y.
!
!    Output, real Y(P), the product S * X.
!
!    Input, real S((P*(P+1))/2), the P by P symmetric matrix.  Only the
!    lower triangle is stored, by rows.
!
!    Input, real X(P), the vector to be multiplied by S.
!
  implicit none

  integer p

  integer i
  integer j
  integer k
  real s((p*(p+1))/2)
  real x(p)
  real xi
  real y(p)
!
!  Compute the lower triangle of S times X.
!
  j = 1
  do i = 1, p
    y(i) = dot_product ( s(j:j+i-1), x(1:i) )
    j = j + i
  end do
!
!  Compute the strict upper triangle of S times X.
!
  j = 1
  do i = 2, p
    j = j + 1
    do k = 1, i - 1
      y(k) = y(k) + s(j) * x(i)
      j = j + 1
    end do
  end do

  return
end
function stopx ( idummy )

!*****************************************************************************80
!
!! STOPX is called to stop execution.
!
!  Discussion:
!
!    This function may serve as the STOPX (asynchronous interruption)
!    function for the NL2SOL package at those installations which do not
!    wish to implement a dynamic STOPX.
!
!    At installations where the NL2SOL system is used
!    interactively, this dummy STOPX should be replaced by a
!    function that returns TRUE if and only if the interrupt
!    (break) key has been pressed since the last call on STOPX.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
  implicit none

  integer idummy
  logical stopx

  stopx = .false.

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
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
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
subroutine vcopy ( p, y, x )

!*****************************************************************************80
!
!! VCOPY copies a vector.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input, integer P, the length of the vector.
!
!    Output, real Y(P), the copy of the input vector.
!
!    Input, real X(P), the vector to be copied.
!
  implicit none

  integer p

  real x(p)
  real y(p)

  y(1:p) = x(1:p)

  return
end
function v2norm ( p, x )

!*****************************************************************************80
!
!! V2NORM computes the L2 norm of a vector.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    David Gay
!
!  Parameters:
!
!    Input, integer P, the length of the vector.
!
!    Input, real X(P), the vector.
!
!    Output, real V2NORM, the Euclidean norm of the vector.
!
!  Local Parameters:
!
!    SQTETA is (slightly larger than) the square root of the
!    smallest positive floating point number on the machine.
!    The tests involving SQTETA are done to prevent underflows.
!
  implicit none

  integer p

  integer i
  integer j
  real r
  real scale
  real, save :: sqteta = 0.0E+00
  real t
  real v2norm
  real x(p)
  real xi

  if ( p <= 0 ) then
    v2norm = 0.0E+00
    return
  end if

  i = 0

  do j = 1, p
    if ( x(j) /= 0.0E+00 ) then
      i = j
      exit
    end if
  end do

  if ( i == 0 ) then
    v2norm = 0.0E+00
    return
  end if

  scale = abs ( x(i) )

  t = 1.0E+00

  if ( sqteta == 0.0E+00 ) then
    sqteta = sqrt ( 1.001E+00 * tiny ( sqteta ) )
  end if

  j = i + 1

  do i = j, p

    xi = abs ( x(i) )

    if ( xi <= scale ) then

      r = xi / scale

      if ( sqteta < r ) then
        t = t + r * r
      end if

    else

      r = scale / xi

      if ( sqteta < r ) then
        t = 1.0E+00 + t * r * r
      else
        t = 1.0E+00
      end if

      scale = xi

    end if

  end do

  v2norm = scale * sqrt ( t )

  return
end
