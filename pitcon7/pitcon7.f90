subroutine pitcon ( df, fpar, fx, ierror, ipar, iwork, liw, nvar, rwork, &
  lrw, xr, slname )

!*****************************************************************************80
!
!! PITCON is the user-interface routine for the continuation code.
!
!  A) Introduction:
!
!  PITCON solves nonlinear systems with one degree of freedom.
!
!  PITCON is given an N dimensional starting point X, and N-1 nonlinear
!  functions F, with F(X) = 0.  Generally, there will be a connected
!  curve of points Y emanating from X and satisfying F(Y) = 0.  PITCON
!  produces successive points along this curve.
!
!  The program can be used to study many sorts of parameterized problems,
!  including structures under a varying load, or the equilibrium
!  behavior of a physical system as some quantity is varied.
!
!  PITCON is a revised version of ACM TOMS algorithm 596.
!
!  Both versions are available via NETLIB, the electronic software
!  distribution service.  NETLIB has the original version in its TOMS
!  directory, and the current version in its CONTIN directory.
!  For more information, send the message "send index from contin"
!  to "netlib@research.att.com".
!
!
!  B) Acknowledgements:
!
!  PITCON was written by
!
!    Professor Werner C Rheinboldt and John Burkardt,
!    Department of Mathematics and Statistics
!    University of Pittsburgh,
!    Pittsburgh, Pennsylvania, 15260, USA.
!
!    E-Mail: wcrhein@vms.cis.pitt.edu
!            burkardt@psc.edu
!
!  The original work on this package was partially supported by the National
!  Science Foundation under grants MCS-78-05299 and MCS-83-09926.
!
!
!  C) Overview:
!
!  PITCON computes a sequence of solution points along a one dimensional
!  manifold of a system of nonlinear equations F(X) = 0 involving NVAR-1
!  equations and an NVAR dimensional unknown vector X.
!
!  The operation of PITCON is somewhat analogous to that of an initial value
!  ODE solver.  In particular, the user must begin the computation by
!  specifying an approximate initial solution, and subsequent points returned
!  by PITCON lie on the curve which passes through this initial point and is
!  implicitly defined by F(X) = 0.  The extra degree of freedom in the system is
!  analogous to the role of the independent variable in a differential
!  equations.
!
!  However, PITCON does not try to solve the algebraic problem by turning it
!  into a differential equation system.  Unlike differential equations, the
!  solution curve may bend and switch back in any direction, and there may be
!  many solutions for a fixed value of one of the variables.  Accordingly,
!  PITCON is not required to parametrize the implicitly defined curve with a
!  fixed parameter.  Instead, at each step, PITCON selects a suitable variable
!  as the current parameter and then determines the other variables as
!  functions of it.  This allows PITCON to go around relatively sharp bends.
!  Moreover, if the equations were actually differentiated - that is, replaced
!  by some linearization - this would introduce an inevitable "drift" away from
!  the true solution curve.  Errors at previous steps would be compounded in a
!  way that would make later solution points much less reliable than earlier
!  ones.  Instead, PITCON solves the algebraic equations explicitly and each
!  solution has to pass an acceptance test in an iterative solution process
!  with tolerances provided by the user.
!
!  PITCON is only designed for systems with one degree of freedom.  However,
!  it may be used on systems with more degrees of freedom if the user reduces
!  the available degrees of freedom by the introduction of suitable constraints
!  that are added to the set of nonlinear equations.  In this sense, PITCON may
!  be used to investigate the equilibrium behavior of physical systems with
!  several degrees of freedom.
!
!  Program options include the ability to search for solutions for which a
!  given component has a specified value.  Another option is a search for a
!  limit or turning point with respect to a given component; that is, of a
!  point where this particular solution component has a local extremum.
!
!  Another feature of the program is the use of two work arrays, IWORK and
!  RWORK.  All information required for continuing any interrupted computation
!  is saved in these two arrays.
!
!
!  D) PITCON Calling Sequence:
!
!  subroutine PITCON(DF,FPAR,FX,IERROR,IPAR,IWORK,LIW,NVAR,RWORK,LRW,XR,SLVNAM)
!
!  On the first call, PITCON expects a point XR and a routine FX defining a
!  nonlinear function F.  Together, XR and FX specify a curve of points Y
!  with the property that F(Y) = 0.
!
!  On the first call, PITCON simply verifies that F(XR) = 0.  If this is not
!  the case, the program attempts to correct XR to a new value satisfying
!  the equation.
!
!  On subsequent calls, PITCON assumes that the input vector XR contains the
!  point which had been computed on the previous call.  It also assumes that
!  the work arrays IWORK and RWORK contain the results of the prior
!  calculations.  PITCON estimates an appropriate stepsize, computes the
!  tangent direction to the curve at the given input point, and calculates a
!  predicted new point on the curve.  A form of Newton's method is used to
!  correct this point so that it lies on the curve.  If the iteration is
!  successful, the code returns with a new point XR.  Otherwise, the stepsize
!  may be reduced, and the calculation retried.
!
!  Aside from its ability to produce successive points on the solution curve,
!  PITCON may be asked to search for "target points" or "limit points".
!  Target points are solution vectors for which a certain component has a
!  specified value.  One might ask for all solutions for which XR(17) = 4.0, for
!  instance.  Limit points occur when the curve turns back in a given
!  direction, and have the property that the corresponding component of the
!  tangent vector vanishes there.
!
!  If the user has asked for the computation of target or limit points, then
!  PITCON will usually proceed as described earlier, producing a new
!  continuation point on each return.  But if a target or limit point is
!  found, such a point is returned as the value of XR, temporarily interrupting
!  the usual form of the computation.
!
!
!  E) Overview of PITCON parameters:
!
!  Names of routines:
!
!    DF     Input,        external DF, evaluates the Jacobian of F.
!    FX     Input,        external FX, evaluates the function F.
!    SLVNAM Input,        external SLVNAM, solves the linear systems.
!
!  Information about the solution point:
!
!    NVAR   Input,        integer NVAR, number of variables, dimension of XR.
!    XR     Input/output, real XR(NVAR), the current solution point.
!
!  Workspace:
!
!    LIW    Input,        integer LIW, the dimension of IWORK.
!    IWORK  Input/output, integer IWORK(LIW), work array.
!
!    LRW    Input,        integer LRW, the dimension of RWORK.
!    RWORK  Input/output, real RWORK(LRW), work array.
!
!    FPAR   "Throughput", real FPAR(*), user defined parameter array.
!    IPAR   "Throughput", integer IPAR(*), user defined parameter array.
!
!  Error indicator:
!
!    IERROR Output,       integer IERROR, error return flag.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  DF     Input, external DF, the name of the Jacobian evaluation routine.
!         This name must be declared external in the calling program.
!
!         DF is not needed if the finite difference option is used
!         (IWORK(9) = 1 or 2). In that case, only a dummy name is needed for DF.
!
!         Otherwise, the user must write a routine which evaluates the
!         Jacobian matrix of the function FX at a given point X and stores it
!         in the FJAC array in accordance with the format used by the solver
!         specified in SLVNAM.
!
!         In the simplest case, when the full matrix solverDENSLV solver
!         provided with the package is used, DF must store  D F(I)/D X(J) into
!         FJAC(I,J).
!
!         The array to contain the Jacobian will be zeroed out before DF is
!         called, so that only nonzero elements need to be stored.  DF must
!         have the form:
!
!       subroutine DF(NVAR,FPAR,IPAR,X,FJAC,IERROR)
!
!           NVAR   Input, integer NVAR, number of variables.
!
!           FPAR   Input, real FPAR(*), vector for passing parameters.
!                  This vector is not used by the program, and is only provided
!                  for the user's convenience.
!
!           IPAR   Input, integer IPAR(*), vector for passing integer
!                  parameters.  This vector is not used by the program, and is
!                  only provided for the user's convenience.
!
!           X      Input, real X(NVAR), the point at which the
!                  Jacobian is desired.
!
!           FJAC   Output, real FJAC(*), array containing Jacobian.
!
!                  If DENSLV is the solver:  FJAC must be dimensioned
!                  FJAC(NVAR,NVAR) as shown above, and DF sets
!                  FJAC(I,J) = D F(I)/DX(J).
!
!                  If BANSLV is the solver:  the main portion of the Jacobian,
!                  rows and columns 1 through NVAR-1, is assumed to be a banded
!                  matrix in the standard LINPACK form with lower bandwidth ML
!                  and upper bandwidth MU.  However, the final column of the
!                  Jacobian is allowed to be full.
!
!                  BANSLV will pass to DF the beginning of the storage for
!                  FJAC, but it is probably best not to doubly dimension FJAC
!                  inside of DF, since it is a "hybrid" object.  The first
!                  portion of it is a (2*ML+MU+1, NEQN) array, followed by a
!                  single column of length NEQN (the last column of the
!                  Jacobian).  Thus the simplest approach is to declare FJAC to
!                  be a vector, and then then to store values as follows:
!
!                    If J is less than NVAR, then
!                      if I-J .LE. ML and J-I .LE. MU,
!                        set K = (2*ML+MU)*J + I - ML
!                        set FJAC(K) = D F(I)/DX(J).
!                      else
!                        do nothing, index is outside the band
!                    endif
!                    Else if J equals NVAR, then
!                      set K = (2*ML+MU+1)*(NVAR-1)+I,
!                      set FJAC(K) = D F(I)/DX(J).
!                  endif.
!
!           IERROR Output, integer IERROR, error return flag.  DF should set
!                  this to 0 for normal return, nonzero if trouble.
!
!  FPAR   Input/output, real FPAR(*), a user defined parameter array.
!
!         FPAR is not used in any way by PITCON.  It is provided for the user's
!         convenience.  It is passed to DF and FX, and hence may be used to
!         transmit information between the user calling program and these user
!         subprograms. The dimension of FPAR and its contents are up to the
!         user.  Internally, the program declares DIMENSION FPAR(*) but never
!         references its contents.
!
!  FX     Input, external FX, the name of the routine to evaluate the function.
!         FX computes the value of the nonlinear function.  This name must be
!         declared external in the calling program.  FX should evaluate the
!         NVAR-1 function components at the input point X, and store the result
!         in the vector FVEC.  An augmenting equation will be stored in entry
!         NVAR of FVEC by the PITCON program.  FX should have the form:
!
!           subroutine FX ( NVAR, FPAR, IPAR, X, FVEC, IERROR )
!
!           Input, integer NVAR, number of variables.
!           Input/output, real FPAR(*), user parameters.
!           Input/output, integer IPAR(*), array of user parameters.
!           Input, real X(NVAR), the point of evaluation.
!           Output, real FVEC(NVAR-1), the value of the function at X.
!           Output, integer IERROR, 0 for no errors, nonzero for an error.
!
!  IERROR Output, integer IERROR, error return flag.
!
!         On return from PITCON, a nonzero value of IERROR is a warning of some
!         problem which may be minor, serious, or fatal.
!
!         0, No errors occurred.
!
!         1, Insufficient storage provided in RWORK and IWORK, or NVAR is less
!            than 2.  This is a fatal error, which occurs on the first call to
!            PITCON.
!
!         2, A user defined error condition occurred in the FX or DF
!          subroutines.  PITCON treats this as a fatal error.
!
!         3, A numerically singular matrix was encountered.  Continuation
!            cannot proceed without some redefinition of the problem.  This is
!            a fatal error.
!
!         4, Unsuccessful corrector iteration.  Loosening the tolerances
!            RWORK(1) and RWORK(2), or decreasing the maximum stepsize RWORK(4)
!            might help.  This is a fatal error.
!
!         5, Too many corrector steps.  The corrector iteration was proceeding
!            properly, but too slowly.  Increase number of Newton steps
!            IWORK(17), increase the error tolerances RWORK(1) or RWORK(2), or
!            decrease RWORK(4).  This is a fatal error.
!
!         6, Null tangent vector.  A serious error which indicates a data
!            problem or singularity in the nonlinear system.  This is a fatal
!            error.
!
!         7, Root finder failed while searching for a limit point.
!            This is a warning.  It means that the limit point
!            computation has failed, but the main computation (computing the
!            continuation curve itself) may continue.
!
!         8, Limit point iteration took too many steps.  This is a warning
!            error.  It means that the limit point computation has failed, but
!            the main computation (computing the continuation curve itself) may
!            continue.
!
!         9, Target point calculation failed.  This generally means that
!            the program detected the existence of a target point, set up
!            an initial estimate of its value, and applied the corrector
!            iteration, but that this corrector iteration failed.
!            This is a only a warning message.  PITCON can proceed to compute
!            new points on the curve.  However, if the target point was
!            really desired, PITCON has no automatic recovery method to
!            retry the calculation.  The best prescription in that case
!            is to try to guarantee that PITCON is taking smaller steps
!            when it detects the target point, which you may do by reducing
!            HMAX, stored as RWORK(4).
!
!         10, Undiagnosed error condition.  This is a fatal error.
!
!  IPAR   Input/output, integer IPAR(*), user defined parameter array.
!
!         IPAR is not used in any way by PITCON.  It is provided for the user's
!         convenience in transmitting parameters between the calling program
!         and the user routines FX and DF.  IPAR is declared in the PITCON
!         program and passed through it to FX and DF, but otherwise ignored.
!         Note, however, that if BANSLV is used for the solver routine, then
!         IPAR(1) must contain the lower bandwidth, and IPAR(2) the upper
!         bandwidth of the Jacobian matrix.
!
!  IWORK  Input/output, integer IWORK(LIW).  Communication and workspace array.
!
!         The specific allocation of IWORK is described in the section devoted
!         to the work arrays.  Some elements of IWORK MUST be set by the user,
!         others may be set to change the way PITCON works.
!
!  LIW    Input, integer LIW, the dimension of IWORK.
!
!         The minimum acceptable value of LIW depends on the solver chosen,
!         but for either DENSLV or BANSLV, setting LIW = 29+NVAR is sufficient.
!
!  NVAR   Input, integer NVAR, the number of variables, the dimension of X.
!
!         This is, of course, one greater than the number of equations or
!         functions.  NVAR must be at least 2.
!
!  RWORK  Input/output, real RWORK(LRW), work array.
!
!         The specific allocation of RWORK is described in the section
!         devoted to the work arrays.
!
!  LRW    Input, integer LRW, the dimension of RWORK.
!
!         The minimum acceptable value depends heavily on the solver options.
!         There is storage required for scalars, vectors, and the Jacobian
!         array.  The minimum acceptable value of LRW is the sum of three
!         corresponding numbers.
!
!         For DENSLV with user-supplied Jacobian,
!
!           LRW = 29 + 4*NVAR + NVAR*NVAR.
!
!         For DENSLV with internally approximated Jacobian,
!
!           LRW = 29 + 6*NVAR + NVAR*NVAR.
!
!         For BANSLV, with a Jacobian matrix with upper bandwidth MU and lower
!         bandwidth ML, and NBAND = 2*ML+MU+1, with user supplied Jacobian,
!
!           LRW = 29 + 6*NVAR + (NVAR-1)*NBAND.
!
!         For BANSLV with internally approximated Jacobian,
!
!           LRW = 29 + 9*NVAR + (NVAR-1)*NBAND.
!
!  XR     Input/output, real XR(NVAR), the current solution point.
!
!         On the first call, the user should set XR to a starting point which
!         at least approximately satisfies F(XR) = 0.  The user need never
!         update XR again.
!
!         Thereafter, on each return from the program with IERROR = 0, XR will
!         hold the most recently computed point, whether a continuation, target
!         or limit point.
!
!  SLVNAM Input, external SLVNAM, the name of the solver to use on linear
!         systems.
!
!         The linear systems have the form A*x = b, where A is the augmented
!         Jacobian matrix.  A will be square, and presumably nonsingular.
!         The routine must return the value of the solution x.
!
!         Two possible choices for SLVNAM are "DENSLV" and "BANSLV", which are
!         the names of routines provided with the package.  DENSLV is
!         appropriate for a full storage jacobian, and BANSLV for a jacobian
!         which is banded except for the last column.
!
!         The advanced user may study the source code for these two routines
!         and write an equivalent solver more suited to a given problem.
!
!
!  G) The Integer Work Array IWORK:
!
!  Input to PITCON includes the setting of some of the entries in IWORK.
!  Some of this input is optional.  The user input section of IWORK involves
!  entries 1 through 9, and, possibly also 17.
!
!  IWORK(1) must be set by the user.  All other entries have default values.
!
!
!  IWORK(1)        On first call only, the user must set IWORK(1) = 0.
!                  Thereafter, the program sets IWORK(1) before return to
!                  explain what kind of point is being returned.  This return
!                  code is:
!
!                      1 return of corrected starting point.
!                      2 return of continuation point.
!                      3 return of target point.
!                      4 return of limit point.
!
!                  NOTE:  At any time, PITCON may be called with a negative
!                  value of IWORK(1). This requests a check of the
!                  jacobian routine against a finite difference approximation,
!                  or a printout of the jacobian or its approximation.
!
!                      -1, compare Jacobian against forward difference,
!                          print maximum discrepancy only.
!                      -2, compare Jacobian against central difference,
!                          print maximum discrepancy only.
!                      -3, compare Jacobian against forward difference,
!                          print out the entire discrepancy matrix.
!                      -4, compare Jacobian against central difference,
!                          print out the entire discrepancy matrix.
!                      -5, compute forward difference Jacobian,
!                          print maximum entry only.
!                      -6, compute central difference Jacobian,
!                          print maximum entry only.
!                      -7, compute forward difference Jacobian,
!                          print entire matrix.
!                      -8, compute central difference Jacobian,
!                          print entire matrix.
!                      -9, request user supplied Jacobian,
!                          print maximum entry only.
!                     -10, request user supplied Jacobian,
!                          print entire matrix.
!
!                  Before a call with negative IWORK(1), the current value of
!                  IWORK(1) should be saved, and then restored to the previous
!                  value after the call, in order to resume calculation.
!
!                  IWORK(1) does not have a default value.  The user MUST set
!                  it.
!
!  IWORK(2)        IPC, the component of the current continuation point XR
!                  which is to be used as the continuation parameter.  On first
!                  call, the program is willing to use the index NVAR as a
!                  default, but the user should set this value if better
!                  information is available.
!
!                  After the first call, the program sets this value for each
!                  step automatically unless the user prevents this by setting
!                  the parameterization option IWORK(3) to a non-zero valus.
!                  Note that a poor choice of IWORK(2) may cause the algorithm
!                  to fail.  IWORK(2) defaults to NVAR on the first step.
!
!  IWORK(3)        Parameterization option.  The program would prefer to be
!                  free to choose a new local parameter from step to step.
!                  The value of IWORK(3) allows or prohibits this action.
!                  IWORK(3) = 0 allows the program to vary the parameter,
!                  IWORK(3) = 1 forces the program to use whatever the contents
!                  of IWORK(2) are, which will not be changed from the user's
!                  input or the default.  The default is IWORK(3) = 0.
!
!  IWORK(4)        Newton iteration Jacobian update option.
!                  0, the Jacobian is reevaluated at every step of the
!                     Newton iteration.  This is costly, but may result in
!                     fewer Newton steps and fewer Newton iteration rejections.
!                  1, the Jacobian is evaluated only on the first and
!                     IWORK(17)-th steps of the Newton process.
!                  2, the Jacobian is evaluated only when absolutely
!                     necessary, namely, at the first step, and when the
!                     process fails. This option is most suitable for problems
!                     with mild nonlinearities.
!
!                  The default is IWORK(4) = 0.
!
!  IWORK(5)        IT, target point index.  If IWORK(5) is not zero, it is
!                  presumed to be the component index between 1 and NVAR for
!                  which target points are sought.  In this case, the value of
!                  RWORK(7) is assumed to be the target value.  The program
!                  will monitor every new continuation point, and if it finds
!                  that a target point may lie between the new point and the
!                  previous point, it will compute this target point and
!                  return.  This target point is defined by the property that
!                  its component with the index prescribed in IWORK(5) will
!                  have the value given in RWORK(7).  For a given problem there
!                  may be zero, one, or many target points requested.
!                  The default of IWORK(5) is 0.
!
!  IWORK(6)        LIM, the limit point index.  If IWORK(6) is nonzero, then
!                  the program will search for limit points with respect to
!                  the component with index IWORK(6); that is, of points for
!                  which the IWORK(6)-th variable has a local extremum, or
!                  equivalently where the IWORK(6)-th component of the tangent
!                  vector is zero.  The default of IWORK(6) is zero.
!
!  IWORK(7)        IWRITE, which controls the amount of output produced by the
!                  program. IWORK(7) may have a value between 0 and 3.
!                  For IWORK(7) = 0 there is almost no output while for
!                  IWORK(7) = 3 the most output is produced.
!                  The default is 1.
!
!  IWORK(9)        Control of the Jacobian option specifying whether the user
!                  has supplied a Jacobian routine, or wants the program
!                  to approximate the Jacobian.
!                  0, the user has supplied the Jacobian.
!                  1, program is to use forward difference approximation.
!                  2, program is to use central difference approximation.
!                  IWORK(9) defaults to 0.
!
!  IWORK(10)       State indicator of the progress of the program.
!                  The values are:
!                  0, start up with unchecked starting point.
!                  1, first step.  Corrected starting point available.
!                  2, two successive continuation points available, as well
!                     as the tangent vector at the oldest of them.
!                  3, two successive continuation points available, as well
!                     as the tangent vector at the newest of them.
!
!  IWORK(11)       Index of the last computed target point. This is used to
!                  avoid repeated computation of a target point.  If a target
!                  point has been found, then the target index IWORK(5) is
!                  copied into IWORK(11).
!
!  IWORK(12)       Second best choice for the local parameterization index.
!                  This index may be tried if the first choice causes poor
!                  performance in the Newton corrector.
!
!  IWORK(13)       Beginning location in IWORK of unused integer work space
!                  available for use by the solver.
!
!  IWORK(14)       LIW, the user declared dimension of the array IWORK.
!
!  IWORK(15)       Beginning location in RWORK of unused real work space
!                  available for use by the solver.
!
!  IWORK(16)       LRW, the user declared dimension of RWORK.
!
!  IWORK(17)       Maximum number of corrector steps allowed during one run
!                  of the Newton process in which the Jacobian is updated at
!                  every step.  If the Jacobian is only evaluated at
!                  the beginning of the Newton iteration then 2*IWORK(17) steps
!                  are allowed.
!                  IWORK(17) must be greater than 0.  It defaults to 10.
!
!  IWORK(18)       Number of stepsize reductions that were needed for
!                  producing the last continuation point.
!
!  IWORK(19)       Total number of calls to the user Jacobian routine DF.
!
!  IWORK(20)       Total number of calls to the matrix factorization routine.
!                  If DENSLV is the chose solver then factorization is done by
!                  the LINPACK routine SGEFA.  If BANSLV is the solver, the
!                  LINPACK routine SGBFA will be used.
!
!  IWORK(21)       Total number of calls to the back-substitution routine.
!                  If DENSLV is the chosen solver, then back substitution is
!                  done by the LINPACK routine SGESL.  If BANSLV is used, then
!                  the LINPACK routine SGBSL will be used.
!
!  IWORK(22)       Total number of calls to the user function routine FX.
!
!  IWORK(23)       Total number of steps taken in limit point iterations.
!                  Each step involves determining an approximate limit point
!                  and applying a Newton iteration to correct it.
!
!  IWORK(24)       Total number of Newton corrector steps used during the
!                  computation of target points.
!
!  IWORK(25)       Total number of Newton steps taken during the correction
!                  of a starting point or the continuation points.
!
!  IWORK(26)       Total number of predictor stepsize-reductions needed
!                  since the start of the continuation procesds.
!
!  IWORK(27)       Total number of calls to the program.  This also
!                  corresponds to the number of points computed.
!
!  IWORK(28)       Total number of Newton steps taken during current iteration.
!
!  IWORK(30)       and on are reserved for use by the linear equation solver,
!                  and typically are used for pivoting.
!
!
!  H) The Real Work Array RWORK:
!
!  Input to PITCON includes the setting of some of the entries in RWORK.
!  Some of this input is optional.  The user input section of RWORK involves
!  entries 1 through 7 and possibly 18 and 20.
!
!  All entries of RWORK have default values.
!
!
!  RWORK(1)        ABSERR, absolute error tolerance.   This value is used
!                  mainly during the Newton iteration.  RWORK(1) defaults to
!                  SQRT(EPMACH) where EPMACH is the machine relative precision
!                  stored in RWORK(8).
!
!  RWORK(2)        RELERR, relative error tolerance.  This value is used mainly
!                  during the Newton iteration.  RWORK(2) defaults to
!                  SQRT(EPMACH) where EPMACH is the machine relative precision
!                  stored in RWORK(8).
!
!  RWORK(3)        HMIN, minimum allowable predictor stepsize.  If failures of
!                  the Newton correction force the stepsize down to this level,
!                  then the program will give up.  The default value is
!                  SQRT(EPMACH).
!
!  RWORK(4)        HMAX, maximum allowable predictor step.  Too generous a value
!                  may cause erratic behavior of the program.  The default
!                  value is SQRT(NVAR).
!
!  RWORK(5)        HTAN,  the predictor stepsize.  On first call, it should be
!                  set by the user.  Thereafter it is set by the program.
!                  RWORK(5) should be positive.  In order to travel in the
!                  negative direction, see RWORK(6).
!                  The default initial value equals 0.5*(RWORK(3)+RWORK(4)).
!
!  RWORK(6)        The local continuation direction, which is either +1.0
!                  or -1.0 .  This asserts that the program is moving in the
!                  direction of increasing or decreasing values of the local
!                  continuation variable, whose index is in IWORK(2).  On first
!                  call, the user must choose IWORK(2).  Therefore, by setting
!                  RWORK(6), the user may also specify whether the program is
!                  to move initially to increase or decrease the variable whose
!                  index is IWORK(2).
!                  RWORK(6) defaults to +1.
!
!  RWORK(7)        A target value.  It is only used if a target index
!                  has been specified through IWORK(5).  In that case, solution
!                  points with the IWORK(5) component equal to RWORK(7) are
!                  to be computed. The code will return each time it finds such
!                  a point.  RWORK(7) does not have a default value.  The
!                  program does not set it, and it is not referenced unless
!                  IWORK(5) has been set.
!
!  RWORK(8)        EPMACH, the value of the machine precision.  The computer
!                  can distinguish 1.0+EPMACH from 1.0, but it cannot
!                  distinguish 1.0+(EPMACH/2) from 1.0. This number is used
!                  when estimating a reasonable accuracy request on a given
!                  computer.  PITCON computes a value for EPMACH internally.
!
!  RWORK(9)        STEPX, the size, using the maximum-norm, of the last
!                  step of Newton correction used on the most recently
!                  computed point, whether a starting point, continuation
!                  point, limit point or target point.
!
!  RWORK(10)       A minimum angle used in the steplength computation,
!                  equal to 2.0*ARCCOS(1-EPMACH).
!
!  RWORK(11)       Estimate of the angle between the tangent vectors at the
!                  last two continuation points.
!
!  RWORK(12)       The pseudo-arclength coordinate of the previous continuation
!                  pointl; that is, the sum of the Euclidean distances between
!                  all computed continuation points beginning with the start
!                  point.  Thus each new point should have a larger coordinate,
!                  except for target and limit points which lie between the two
!                  most recent continuation points.
!
!  RWORK(13)       Estimate of the pseudo-arclength coordinate of the current
!                  continuation point.
!
!  RWORK(14)       Estimate of the pseudoarclength coordinate of the current
!                  limit or target point, if any.
!
!  RWORK(15)       Size of the correction of the most recent continuation
!                  point; that is, the maximum norm of the distance between the
!                  predicted point and the accepted corrected point.
!
!  RWORK(16)       Estimate of the curvature between the last two
!                  continuation points.
!
!  RWORK(17)       Sign of the determinant of the augmented matrix at the
!                  last continuation point whose tangent vector has been
!                  calculated.
!
!  RWORK(18)       This quantity is only used if the jacobian matrix is to
!                  be estimated using finite differences.  In that case,
!                  this value determines the size of the increments and
!                  decrements made to the current solution values, according
!                  to the following formula:
!
!                    Delta X(J) = RWORK(18) * (1.0 + ABS(X(J))).
!
!                  The value of every entry of the approximate jacobian could
!                  be extremely sensitive to RWORK(18).  Obviously, no value
!                  is perfect.  Values too small will surely cause singular
!                  jacobians, and values too large will surely cause inaccuracy.
!                  Little more is certain.  However, for many purposes, it
!                  is suitable to set RWORK(18) to the square root of the
!                  machine epsilon, or perhaps to the third or fourth root,
!                  if singularity seems to be occuring.
!
!                  RWORK(18) defaults to SQRT(SQRT(RWORK(8))).
!
!                  The user may set RWORK(18).  If it has a nonzero value on
!                  input, that value will be used.  Otherwise, the default
!                  value will be used.
!
!  RWORK(19)       Not currently used.
!
!  RWORK(20)       Maximum growth factor for the predictor stepsize based
!                  on the previous secant stepsize.  The stepsize algorithm
!                  will produce a suggested step that is no less that the
!                  previous secant step divided by this factor, and no greater
!                  than the previous secant step multiplied by that factor.
!                  RWORK(20) defaults to 3.
!
!  RWORK(21)       The (Euclidean) secant distance between the last two
!                  computed continuation points.
!
!  RWORK(22)       The previous value of RWORK(21).
!
!  RWORK(23)       A number judging the quality of the Newton corrector
!                  convergence at the last continuation point.
!
!  RWORK(24)       Value of the component of the current tangent vector
!                  corresponding to the current continuation index.
!
!  RWORK(25)       Value of the component of the previous tangent vector
!                  corresponding to the current continuation index.
!
!  RWORK(26)       Value of the component of the current tangent vector
!                  corresponding to the limit index in IWORK(6).
!
!  RWORK(27)       Value of the component of the previous tangent vector
!                  corresponding to the limit index in IWORK(6).
!
!  RWORK(28)       Value of RWORK(7) when the last target point was
!                  computed.
!
!  RWORK(29)       Sign of the determinant of the augmented matrix at the
!                  previous continuation point whose tangent vector has been
!                  calculated.
!
!  RWORK(30)       through RWORK(30+4*NVAR-1) are used by the program to hold
!                  an old and new continuation point, a tangent vector and a
!                  work vector.  Subsequent entries of RWORK are used by the
!                  linear solver.
!
!
!  I) Programming Notes:
!
!  The minimal input and user routines required to apply the program are:
!
!    Write a function routine FX of the form described above.
!    Use DENSLV as the linear equation solver by setting SLVNAM to DENSLV.
!    Skip writing a Jacobian routine by using the finite difference option.
!    Pass the name of FX as the Jacobian name as well.
!    Declare the name of the function FX as external.
!    Set NVAR in accordance with your problem.
!
!  Then:
!
!    Dimension the vector IWORK to the size LIW = 29+NVAR.
!    Dimension the vector RWORK to the size LRW = 29+NVAR*(NVAR+6).
!    Dimension IPAR(1) and FPAR(1) as required in the function routine.
!    Dimension XR(NVAR) and set it to an approximate solution of F(XR) = 0.
!
!  Set the work arrays as follows:
!
!    Initialize IWORK to 0 and RWORK to 0.0.
!
!    Set IWORK(1) = 0 (Problem startup)
!    Set IWORK(7) = 3 (Maximum internally generated output)
!    Set IWORK(9) = 1 (Forward difference Jacobian)
!
!  Now call the program repeatedly, and never change any of its arguments.
!  Check IERROR to decide whether the code is working satisfactorily.
!  Print out the vector XR to see the current solution point.
!
!
!  The most obvious input to try to set appropriately after some practice
!  would be the error tolerances RWORK(1) and RWORK(2), the minimum, maximum
!  and initial stepsizes RWORK(3), RWORK(4) and RWORK(5), and the initial
!  continuation index IWORK(2).
!
!  For speed and efficiency, a Jacobian routine should be written. It can be
!  checked by comparing its results with the finite difference Jacobian.
!
!  For a particular problem, the target and limit point input can be very
!  useful.  For instance, in the case of a discretized boundary value problem
!  with a real parameter it may be desirable to compare the computed solutions
!  for different discretization-dimensions and equal values of the parameter.
!  For this the target option can be used with the prescribed values of the
!  parameter. Limit points usually are of importance in connection with
!  stability considerations.
!
!
!  The routine REPS attempts to compute the machine precision, which in practice
!  is simply the smallest power of two that can be added to 1 to produce a
!  result greater than 1.  If the REPS routine misbehaves, you can replace
!  it by code that assigns a constant precomputed value, or by a call to
!  the PORT/SLATEC routine R1MACH.  REPS is called just once, in the
!  PITCON routine, and the value returned is stored into RWORK(8).
!
!  In subroutines DENSLV and BANSLV, the parameter "EPS" is set to RWORK(18)
!  and used in estimating jacobian matrices via finite differences.  If EPS is
!  too large, the jacobian will be very inaccurate.  Unfortunately, if EPS is
!  too small, the "finite" differences may become "infinitesmal" differences.
!  That is, the difference of two function values at close points may be
!  zero.  This is a very common problem, and occurs even with a function
!  like F(X) = X*X.  A singular jacobian is much worse than an inaccurate one,
!  so we have tried setting the default value of RWORK(18) to
!  SQRT(SQRT(RWORK(8)).  Such a value attempts to ward off singularity at the
!  expense of accuracy.  You may find for a particular problem or machine that
!  this value is too large and should be adjusted.  It is an utterly arbitrary
!  value.
!
!
!  J) Reference:
!
!  1.
!  Werner Rheinboldt,
!  Solution Field of Nonlinear Equations and Continuation Methods,
!  SIAM Journal of Numerical Analysis,
!  Volume 17, 1980, pages 221-237.
!
!  2.
!  Cor den Heijer and Werner Rheinboldt,
!  On Steplength Algorithms for a Class of Continuation Methods,
!  SIAM Journal of Numerical Analysis,
!  Volume 18, 1981, pages 925-947.
!
!  3.
!  Werner Rheinboldt,
!  Numerical Analysis of Parametrized Nonlinear Equations
!  John Wiley and Sons, New York, 1986
!
!  4.
!  Werner Rheinboldt and John Burkardt,
!  A Locally Parameterized Continuation Process,
!  ACM Transactions on Mathematical Software,
!  Volume 9, Number 2, June 1983, pages 215-235.
!
!  5.
!  Werner Rheinboldt and John Burkardt,
!  Algorithm 596, A Program for a Locally Parameterized Continuation Process,
!  ACM Transactions on Mathematical Software,
!  Volume 9, Number 2, June 1983, Pages 236-241.
!
!  6.
!  J J Dongarra, J R Bunch, C B Moler and G W Stewart,
!  LINPACK User's Guide,
!  Society for Industrial and Applied Mathematics,
!  Philadelphia, 1979.
!
!  7.
!  Richard Brent,
!  Algorithms for Minimization without Derivatives,
!  Prentice Hall, 1973.
!
!  8.
!  Tony Chan,
!  Deflated Decomposition of Solutions of Nearly Singular Systems,
!  Technical Report 225,
!  Computer Science Department,
!  Yale University,
!  New Haven, Connecticut, 06520,
!  1982.
!
!
!  L)  Sample programs:
!
!
!  A number of sample problems are included with PITCON.  There are several
!  examples of the Freudenstein-Roth function, which is a nice example
!  with only a few variables, and nice whole number starting and stopping
!  points.  Other sample problems demonstrate or test various options in
!  PITCON.
!
!
!  1:
!  pcprb1.f
!  The Freudenstein-Roth function.  3 variables.
!
!  This is a simple problem, whose solution curve has some severe bends.
!  This file solves the problem with a minimum of fuss.
!
!
!  2:
!  pcprb2.f
!  The Aircraft Stability problem.  8 variables.
!
!  This is a mildly nonlinear problem, whose solution curve has some
!  limit points that are difficult to track.
!
!
!  3:
!  pcprb3.f
!  A boundary value problem.  22 variables.
!
!  This problem has a limit point in the LAMBDA parameter, which we seek.  We
!  solve this problem 6 times, illustrating the use of full and banded
!  jacobians, and of user-generated, or forward or central difference
!  approximated jacobian matrices.
!
!  The first 21 variables represent the values of a function along a grid
!  of 21 points between 0 and 1.  By increasing the number of grid points,
!  the problem can be set up with arbitrarily many variables.  This change
!  requires changing a single parameter in the main program.
!
!
!  4:
!  pcprb4.f
!  The Freudenstein Roth function.  3 variables.
!
!  This version of the problem tests the parameterization fixing option.
!  The user may demand that PITCON always use a given index for continuation,
!  rather than varying the index from step to step.  The Freudenstein-Roth
!  curve may be parameterized by the variable of index 2, although this
!  may increase the number of steps required to traverse the curve.
!
!  This file carries out 6 runs, with and without the parameterization fixed,
!  and with no limit point search, or with limit points sought for index 1
!  or for index 3.
!
!
!  5:
!  pcprb5.f
!  A boundary value problem.  21 variables.
!
!  This problem has a limit point in the LAMBDA parameter, which we do not
!  seek.  We do seek target points where LAMBDA = 0.8, which occurs twice,
!  before and after LAMBDA "goes around the bend".
!
!  We run this test 3 times, comparing the behavior of the full Newton
!  iteration, the Newton iteration that fixes the Jacobian during the
!  correction of each point, and the Newton iteration that fixes the
!  Jacobian as long as possible.
!
!
!  6:
!  pcprb6.f
!  Freudenstein-Roth function.  3 variables.
!
!  This version of the problem demonstrates the jacobian checking option.
!  Two runs are made.  Each is allowed only five steps.  The first
!  run is with the correct jacobian.  The second run uses a defective
!  jacobian, and demonstrates not only the jacobian checker, but also
!  shows that "slightly" bad jacobians can cause the Newton convergence
!  to become linear rather than quadratic.
!
!
!  7:
!  pcprb7.f
!  The materially nonlinear problem, up to 71 variables.
!
!  This is a two point boundary value problem, representing the behavior
!  of a rod under a load.  The shape of the rod is represented by piecewise
!  polynomials, whose form and continuity are specifiable as options.
!
!  The problem as programmed can allow up to 71 variables, but a simple
!  change of a few parameters will allow the problem to be arbitrarily
!  large.
!
!  This is a complicated program.  It is intended to demonstrate the
!  advanced kinds of problems that can be solved with PITCON, as opposed
!  to the "toy" problems like the Freudenstein-Roth function.
!
!
!  8:
!  pcprb8.f
!  The Freudenstein-Roth function.  3 variables.
!
!  This version of the problem tests the jacobian approximation options.
!  We run the problem 3 times, first with the user jacobian, second with
!  the forward difference approximation, and third with the central
!  difference approximation.
!
!
!  Recent Changes:
!
!
!  Changes made for version 7.0:
!
!    Inserted simple versions of LAPACK routines.
!
!  Changes made for version 6.6:
!
!    Some calculations in COQUAL, for ITOP and IBOT, could involve integers
!    to negative powers.  These were replaced by real calculations for TOP
!    and BOT.
!
!    There is a stringent convergence test in CORRECTOR which severely slows
!    down the traversal of the Freudenstein Roth curve, because it forces
!    a last very small step, which causes the computation of the stepsize
!    to be skewed.  I temporarily turned this off.
!
!  Changes made for version 6.5:
!
!    Spun off "STAJAC" from START.
!
!    Code written in lower case.
!
!    Replaced all labeled DO statements with DO/ENDDO loops.
!
!  Changes made for version 6.4:
!
!    Calls to LINPACK replaced by calls to LAPACK
!    Added routines SGEDET and SGBDET to compute determinants.
!
!    Call to SGBTRF, SGBTRS, SGBDET had incorrect value of LDA,
!    which was corrected 21 January 1994.
!
!    Dropped LOUNIT.  All WRITE's are to unit * now.
!
!
!  Changes made for version 6.3:
!
!    27 July 1992: Printout of "Possible bifurcation" message will
!    now occur if IWRITE.GE.1.
!
!    HMIN was reset to MAX(HMIN,SQRT(EPS)) before.  Now it is only
!    set to SQRT(EPS) if HMIN is nonpositive.
!
!    30 September 1991: SKALE was not declared in LIMIT, causing
!    problems for the double precision code.
!
!    26 September 1991: corrected the formulation of the limit
!    point calculation.
!
!    12 September 1991: Jyun-Ming Chen reported that the program did not
!    actually fix the parameter to a user specified value when requested to do
!    so.  That is, when IWORK(3) = 1, the program is never supposed to alter the
!    input value of IWORK(2).  However, the program was doing so.  This error
!    has been corrected now.  The only exception occurs when the initial
!    input value of IWORK(2) is less than 1 or greater than NVAR, in which
!    case the program sets IWORK(2) to NVAR, no matter what the value
!    of IWORK(3).
!
!    The target and limit point calculations, and many other operations,
!    were "spun off" into subroutines.
!
!
!  Changes made for version 6.2:
!
!    The internal documentation was corrected.  It originally stated that
!    the user input portions of IWORK were entries 1 through 8, 17 and 29.
!    This was corrected to 1 through 9 and 17.
!
!    The entry RWORK(18) was previously unused.  It is now used to allow
!    the user to set the value of the finite difference increment used
!    for the approximation of the jacobian.
!
!    Modified DENSLV and BANSLV to allow more user control over the
!    Jacobian comparison or printout.  IWORK(1) may now have the values -1
!    through -10.  The user thus chooses to use forward or central
!    differences, and minimal or maximal printout.
!
!    Modified "Programming notes" section to explain the role of the
!    REPS routine, the new use of IWORK(1), and the new parameter RWORK(18).
!
!
!  Changes made for version 6.0:
!
!    1) When computing a starting point, the Newton convergence
!    criterion was relaxed to require only that the function
!    norm decrease, but not to require that the step norm
!    decrease.
!
!    2) In the Newton correction routine, the MAX norm was replaced
!    by the L2 norm when computing the size of the step and the
!    residual.  This was an attempt to control the 'jerkiness'
!    of the convergence for poorly scaled problems.
!
!    3) Added option to check Jacobian.
!
!    4) Apparently, a programming error in the previous version meant
!    that IWORK(29) was set to zero on startup.  This meant that the
!    program would not realize that the user was using a band storage
!    scheme.  This has been corrected.
!
  external df
  external fx
  external slname

  double precision, parameter :: one = 1.0D+00

  integer liw
  integer lrw
  integer nvar

  double precision dets
  double precision fpar(*)
  integer i
  integer ierror
  integer ifound
  integer ipar(*)
  integer iwork(liw)
  integer iwrite
  integer job
  integer ltc
  integer lwk
  integer lxc
  integer lxf
  integer modnew
  double precision qual
  double precision rwork(lrw)
  double precision xr(nvar)
!
!
!  1.  Set a few parameters.
!
!
  ierror = 0
  iwrite = iwork(7)

  lxc = 29            + 1
  lxf = 29 +     nvar + 1
  ltc = 29 + 2 * nvar + 1
  lwk = 29 + 3 * nvar + 1

  if ( iwork(1) < -10 ) then
    iwork(1) = - 2
  else if ( iwork(1) > 5 ) then
    iwork(1) = 0
  end if
!
!  2.  Preparations.
!
!  Check entries of IWORK and RWORK, compute some constants.
!
  if ( iwork(1) == 0 ) then

    if ( iwrite >= 1 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PITCON 7.0'
      write ( *, * ) '  University of Pittsburgh continuation code'
      write ( *, * ) ' '
      write ( *, * ) '  Modified:        12 November 1999'
      write ( *, * ) '  Linear algebra:  LAPACK'
      write ( *, * ) '  Precision:       double'
      write ( *, * ) ' '
    end if

    call checkw ( ierror, iwork, liw, lrw, nvar, rwork )

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PITCON - Fatal error!'
      write ( *, * ) '  An error was detected in the user input.'
      write ( *, * ) '  PITCON can not begin.'
      stop
    end if

  end if
!
!  3.  Check the user Jacobian routine against a finite difference
!  approximation.
!
  if ( iwork(1) < 0 ) then

    job = 3

    call slname ( dets, fx, df, fpar, ierror, iwork(2), ipar, iwork, &
      liw, job, nvar, rwork, lrw, xr, rwork(lwk) )

    if ( ierror /= 0 ) then
      if ( iwrite > 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'PITCON - Warning!'
        write ( *, * ) '  An error occurred during the jacobian '
        write ( *, * ) '  check.'
      end if
    end if

    return

  end if
!
!  4.  Starting point check
!
!  On first call for a given problem, check that F(XR) is small
!  enough so that the starting point may be considered to lie on
!  the curve.
!
!  If this is not the case, call the corrector to try to enforce it.
!
  if ( iwork(1) == 0 ) then
!
!  See if the initial jacobian should be set up.
!
    if ( iwork(4) == 2 ) then

      call stajac ( df,fpar,fx,ierror,ipar,iwork(2),iwork,liw, &
        lrw,nvar,rwork,rwork(lwk),xr,slname)

      if ( ierror /= 0 ) then
        if ( iwrite > 0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'PITCON - Fatal error!'
          write ( *, * ) '  An error occurred during the'
          write ( *, * ) '  initial jacobian setup.'
          write ( *, * ) ' '
          write ( *, * ) '  The program can not continue!'
          stop
        end if
        return
      end if

    end if
!
!  Check that the starting point satisfies the equations.
!
    call start ( df,fpar,fx,ierror,ipar,iwork(2),iwork,liw, &
      lrw,nvar,rwork,rwork(lwk),rwork(lxc),rwork(lxf),xr,slname)

    if ( ierror /= 0 ) then
      if ( iwrite > 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'PITCON - Fatal error!'
        write ( *, * ) '  An error occurred during the starting'
        write ( *, * ) '  point check.'
        write ( *, * ) ' '
        write ( *, * ) '  The starting point does not satisfy the'
        write ( *, * ) '  accuracy requirements, and PITCON could'
        write ( *, * ) '  not correct it.'
        write ( *, * ) ' '
        write ( *, * ) '  The program can not continue!'
      end if
    end if

    do i = 1, nvar
      rwork(ltc+i-1) = 0.0D+00
    end do

    rwork(ltc+iwork(2)-1) = 1.0D+00
    return

  end if
!
!  5.  Target point
!
!  If IWORK(5) is nonzero, target points are sought.  Check to see if
!  target component IWORK(5), also called "IT", has value lying between
!  XC(IT) and XF(IT).  If so, get linearly interpolated starting point,
!  and use Newton's method to get target point.
!
  call target ( df, fpar, fx, ierror, ifound, ipar, iwork, liw, lrw, &
    nvar, rwork, slname, rwork(lwk), rwork(lxc), rwork(lxf), xr )

  if ( ifound == 1 ) then
    if ( ierror /= 0 ) then
      ierror = 9
      if ( iwrite > 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'PITCON - Warning!'
        write ( *, * ) '  An error occurred during the'
        write ( *, * ) '  target point computation.'
        write ( *, * ) ' '
        write ( *, * ) '  The target point returned does'
        write ( *, * ) '  not satisfy the accuracy requirements.'
        write ( *, * ) ' '
        write ( *, * ) '  However, the code can continue.'
      end if
    end if
    return
  end if
!
!  6.  Tangent and local continuation parameter calculation.
!
!  Unless the tangent and limit point calculations were already
!  performed, (because the loop was interrupted for a limit point
!  calculation), set up and solve the equation for the tangent vector.
!
!  Force the tangent vector to be of unit length, and try to preserve
!  the "sense" or "direction" of the curve by forcing the IPL-th
!  component of the tangent vector to agree in sign with the IPL-th
!  component of the previous secant vector.  (On the first step, however,
!  we have to use the user's input direction to choose the sign).
!
!  Set the local continuation parameter IPC.
!
!  If IWORK(3) is 0, the program is free to vary IPC from step to step.
!  In that case, IPC is normally set to the index of the component of
!  greatest magnitude in the tangent vector.  However, if a limit point
!  appears to be coming in that direction, the index of the second
!  greatest magnitude component might be chosen instead.
!
  if ( iwork(10) /= 3 ) then

    call tanpar ( df,fpar,fx,ierror,ipar,iwork,liw,lrw,nvar,rwork,slname, &
      rwork(ltc),rwork(lwk),rwork(lxc),rwork(lxf),xr)

    if ( ierror/= 0 ) then
      if ( iwrite>0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'PITCON - Fatal error.'
        write ( *, * ) '  The computation failed while computing'
        write ( *, * ) '  the parameter and the tangent vector.'
        write ( *, * ) ' '
        write ( *, * ) '  The program can not proceed!'
      end if
      return
    end if
!
!  7.  Limit point check.
!
!  Skip this section if IWORK(6) = 0.
!
!  Otherwise, user has requested a search for limit points in a given
!  index by setting IWORK(6), also called "LIM", to a nonzero value.
!
!  Compare LIM-th components of previous and current tangent vectors.
!  If a sign difference occurs, we assume a limit point has been passed.
!  Attempt to compute a point XR between the previous and current points,
!  for which the LIM-th component of the tangent is zero.
!
!  This search will be guided by a rootfinder.  The scalar function
!  to be zeroed out is the LIM-th tangent vector component.
!
    if ( (iwork(6)/= 0) .and. (iwork(1)/=4).and. (iwork(10) == 3).and. &
       (sign(one,rwork(26)) /= sign(one,rwork(27))) ) then

      call limit(df,fpar,fx,ierror,ipar,iwork,liw,lrw,nvar,rwork, &
        slname,rwork(ltc),rwork(lwk),rwork(lxc),rwork(lxf),xr)

      if ( ierror/= 0 ) then
        if ( iwrite>0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'PITCON - Warning!'
          write ( *, * ) '  An error occurred during the'
          write ( *, * ) '  limit point computation.'
          write ( *, * ) ' '
          write ( *, * ) '  The computed limit point does not'
          write ( *, * ) '  satisfy the accuracy requirements.'
          write ( *, * ) ' '
          write ( *, * ) '  However, the code can continue.'
        end if
      end if

      return

    end if

  end if
!
!  8.  Compute next predictor step length, HTAN.
!
  if ( iwork(10) > 1 ) then
    call setstp ( iwork, liw, lrw, rwork )
  end if
!
!  9.  Continuation step
!
!  Our current data is the current point XC, its tangent vector TC, and
!  a steplength HTAN.  We predict the location of the next point on the
!  curve using the Euler approximation XR = XC+HTAN*TC.
!
!  Newton iteration is applied to this point, to force it to lie on the
!  curve.  In order to make the system square, an augmenting equation
!  is added to the system, specifying that XR(IPC) = XC(IPC)+HTAN*TC(IPC).
!  (The right hand side is a constant.)
!
!  If the Newton correction process fails, the stepsize is reduced and
!  prediction and correction retried.  Failure will most likely be
!  signaled by repeated step reductions, until the minimum allowable
!  stepsize is reached.  If this occurs, PITCON has failed, and cannot
!  proceed along the curve any more.
!
  call trystp ( df, fpar, fx, ierror, ipar, iwork, liw, lrw, nvar, rwork, &
    slname, rwork(ltc), rwork(lwk), rwork(lxf), xr )

  if ( ierror /= 0 ) then
    if ( iwrite > 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PITCON - Fatal error.'
      write ( *, * ) 'The computation failed while trying'
      write ( *, * ) 'to compute the next point.'
      write ( *, * ) ' '
      write ( *, * ) 'The program can not proceed!'
    end if
    return
  end if
!
!  10.  Successful step.  Update information.
!
  call update ( iwork, liw, lrw, nvar, rwork, rwork(ltc), rwork(lxc), &
    rwork(lxf), xr )
!
!  Compute the convergence "quality", a factor between 1/8 and 8, which
!  tries to estimate how deeply inside the Newton attraction region we
!  are, and hence how much bolder or more timid we could be on the next
!  prediction.
!
  modnew = iwork(4)

  call coqual ( modnew, qual, iwork, liw, rwork, lrw )

  if ( 3 <= iwrite ) then
    write ( *, * ) 'PITCON - Corrector convergence quality factor = ', qual
  end if

  rwork(23) = qual

  return
end
subroutine checkw ( ierror, iwork, liw, lrw, nvar, rwork )

!*****************************************************************************80
!
!! CHECKW checks the entries of IWORK and RWORK on the first call.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
  integer liw
  integer lrw

  integer i
  integer ierror
  integer iwork(liw)
  integer iwrite
  integer nvar
  double precision rwork(lrw)
  double precision tcos
  double precision temp
  double precision tsin

  iwrite = iwork(7)

  rwork(11:17) = 0.0D+00
  rwork(21:29) = 0.0D+00

  iwork(10) = 0
  iwork(11) = 0
  iwork(12) = 0
  do i = 18,28
    iwork(i) = 0
  end do

  rwork(8) = epsilon ( rwork(8) )
  if ( iwrite>= 2 ) then
    write ( *, * ) 'CHECKW - Machine epsilon = ',rwork(8)
  end if
!
!  Set the value of the parameter used in estimating the jacobian
!  via finite differences.
!
  if ( rwork(18) == 0.0D+00 ) then
    rwork(18)=sqrt(sqrt(rwork(8)))
  end if

  if ( iwrite>= 3) then
    write ( *, * ) 'CHECKW - Jacobian finite difference increment = ',rwork(18)
  end if
!
!  Set the value of an angle which is "almost" a zero angle.
!
  tcos = sqrt ( 1.0D+00 - rwork(8) )
  tsin = sqrt ( rwork(8) )
  rwork(10) = 2.0D+00 * atan2 ( tsin, tcos )
  if ( nvar <= 1 ) then
    ierror = 1
    if ( iwrite >= 1 ) then
      write ( *, * ) 'CHECKW - Fatal error!'
      write ( *, * ) '  The number of variables, NVAR, must be at least 2.'
      write ( *, * ) '  The input value is ', nvar
    end if
    return
  end if
!
!  Set entries of IWORK which point to the next free location, and
!  total available space in IWORK and RWORK.
!
  iwork(13) = 30
  iwork(14) = liw
  iwork(15) = 29+4*nvar+1
  iwork(16) = lrw
!
!  Check allocated sizes of IWORK and RWORK versus minimal requirements.
!  The actual need for RWORK will generally be even greater than is checked
!  for here, depending on the solver and storage method chosen.
!
  if ( liw<iwork(13) ) then
    if ( iwrite>= 1 ) then
      write ( *, * ) 'CHECKW - Fatal error!'
      write ( *, * ) '  LIW is too small!  Input value of LIW = ',liw
      write ( *, * ) '  The minimum acceptable value = ',iwork(13)
    end if
    ierror = 1
    return
  end if

  if ( lrw<iwork(15) ) then
    if ( iwrite>= 1 ) then
      write ( *, * ) 'CHECKW - Fatal error!'
      write ( *, * ) '  LRW is too small!  Input value of LRW = ',lrw
      write ( *, * ) '  The minimum acceptable value = ',iwork(15)
    end if
    ierror = 1
    return
  end if
!
!  Check entries of IWORK
!
  if ( iwork(2)<1.or.iwork(2)>nvar ) then
    iwork(2) = nvar
    if ( iwrite>= 1 ) then
      write ( *, * ) 'CHECKW - Note:'
      write ( *, * ) '  Start continuation component IWORK(2) at ', iwork(2)
    end if
  end if

  if ( iwork(3) /= 1 ) then
    iwork(3) = 0
  end if

  if ( iwork(4)<0.or.iwork(4)>2)iwork(4) = 0
  if ( iwork(5)<1.or.iwork(5)>nvar)iwork(5) = 0
  if ( iwork(6)<1.or.iwork(6)>nvar)iwork(6) = 0
  if ( iwork(9)<0.or.iwork(9)>2)iwork(9) = 0
  if ( iwork(17)<1)iwork(17) = 10
!
!  Check entries of RWORK.
!
  if ( rwork(1) <= 0.0D+00)rwork(1)=sqrt(rwork(8))
  if ( rwork(2) <= 0.0D+00)rwork(2)=sqrt(rwork(8))

  if ( rwork(3) <= 0.0D+00 ) then
    rwork(3) = sqrt(rwork(8))
    if ( iwork(7)>= 1 ) then
      write ( *, * ) 'CHECKW - The minimum stepsize is too small.'
      write ( *, * ) '  Increasing RWORK(3) to ',rwork(3)
    end if
  end if

  if ( rwork(4)<rwork(3) ) then
    temp = nvar
    rwork(4) = max(rwork(3),sqrt(temp))
  end if

  if ( rwork(5)<0.0D+00 ) then
    rwork(5) = -rwork(5)
    rwork(6) = -rwork(6)
  end if

  if ( rwork(5)<rwork(3).or.rwork(5)>rwork(4) ) then
    rwork(5) = (rwork(3)+rwork(4))/2.0D+00
  end if

  if ( rwork(6) /= (-1.0D+00) ) then
    rwork(6)=1.0D+00
  end if

  if ( rwork(20)<1.0D+00.or.rwork(20)>100.0D+00)rwork(20) = 3.0D+00

  return
end
subroutine coqual ( modnew, qual, iwork, liw, rwork, lrw )

!*****************************************************************************80
!
!! COQUAL computes the correction rate quality.
!
!  Discussion:
!
!    Considerations used include the number of steps taken versus the
!    maximum allowed and the ratio of the last corrector step to the total
!    correction.
!
!    The quality factor, locally called QUAL, is stored in RWORK(23)
!    on return.  Its value is between 1/8 and 8, with 1/8 signifying
!    a poor correction process, 1 an average value, and 8 superior.
!    The value of QUAL is a factor in the size of the next step.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cor den Heijer, Werner Rheinboldt,
!    On Steplength Algorithms for a Class of Continuation Methods,
!    SIAM Journal on Numerical Analysis,
!    Volume 18, Number 5, October 1981, pages 925-947.
!
!  Parameters:
!
!    Input, integer MODNEW, Newton iteration Jacobian update option.
!    * 0, the Jacobian is reevaluated at every step of the Newton iteration.
!    This is costly, but may result in fewer Newton steps and fewer Newton
!    iteration rejections.
!    * 1, the Jacobian is evaluated only on the first and IWORK(17)-th steps
!    of the Newton process.
!    * 2, the Jacobian is evaluated only when absolutely necessary, namely,
!    at the first step, and when the process fails. This option is most suitable
!    for problems with mild nonlinearities.
!
!    Output, double precision QUAL, the correction rate quality.
!
!    Input, integer IWORK(LIW), the integer work array.
!
!    Input, integer LIW, the size of IWORK.
!
!    Input/output, double precision RWORK(LRW), the real work array.
!
!    Input, integer LRW, the size of RWORK.
!
  integer liw
  integer lrw

  double precision base
  double precision bot
  double precision cordis
  double precision esab
  double precision expo
  integer iwork(liw)
  integer maxcor
  integer modnew
  integer nave
  integer nmax
  double precision qual
  double precision rwork(lrw)
  double precision stepx
  double precision term
  double precision test
  double precision top

  stepx = rwork(9)
!
!  RWORK(15) is the Euclidean distance between the first and last
!  points in the Newton correction process.
!
!  RWORK(9) is the size of the last correction.
!
!  IWORK(17) is the maximum number of Newton corrections allowed.
!
!  IWORK(4) defines the type of Newton corrector process used.
!
!  IWORK(28) contains the actual number of Newton correction steps used.
!
  cordis = rwork(15)
  maxcor = iwork(17)
!
!  Was the correction very fast?
!
  if ( iwork(28) <= 1 ) then
    qual = 8.0D+00
    return
  end if
!
!  Was the correction very small?
!
  if ( cordis <= 8.0D+00 * rwork(8) ) then
    qual = 8.0D+00
    return
  end if
!
!  Were an "average" number of steps take?
!
  if ( modnew == 0 ) then
    nave = ( maxcor - 1 ) / 2
  else
    nave = maxcor
  end if

  if ( iwork(28) == nave ) then
    qual = 1.0D+00
    return
  end if
!
!  Were many steps taken?
!
  if ( modnew == 0 ) then
    nmax = maxcor
  else
    nmax = 2 * maxcor
  end if

  if ( nmax <= iwork(28) ) then
    qual = 0.125D+00
    return
  end if
!
!  For MODNEW = 0,
!    W = (STEPX/CORDIS),
!    IEXP = 1/(2**(NCOR-1)-1)
!    U = W**IEXP
!    JEXP = 2**(NCOR-NAVE)
!    QUAL = (U+1+(1/U)) / (U**JEXP+1+(1/U**JEXP))
!
  if ( modnew == 0 ) then
!
!  ?? This TOP isn't used...
!
    top = 2.0D+00**( iwork(28) - nave )
    bot = 2.0D+00**( iwork(28) - 1 ) - 1
    expo = 1.0D+00 / bot
    base = ( stepx / cordis )**expo

    if ( base <= epsilon ( base ) ) then
      qual = 8.0D+00
    else
      top = base + 1.0D+00 + ( 1.0D+00 / base )
      term = base**top
      if ( term <= epsilon ( term ) ) then
        qual = 8.0D+00
      else
        bot = term + 1.0D+00 + ( 1.0D+00 / term )
        qual = top / bot
      end if
    end if
!
!  For MODNEW nonzero,
!    EXP = (NCOR-NAVE)/(NCOR-1)
!    W = (STEPX/CORDIS)
!    QUAL = W**EXP
!
  else

    top = dble ( iwork(28) - nave )
    bot = dble ( iwork(28) - 1 )
    expo = bot / top
    test = 8.0D+00**expo
    base = stepx / cordis
    esab = cordis / stepx

    if ( ( iwork(28) < nave .and. test > base ) .or. &
         ( iwork(28) > nave .and. test < base ) ) then
      qual = 8.0D+00
      return
    end if

    if ( ( iwork(28) < nave .and. test > esab ) .or. &
         ( iwork(28) > nave .and. test < esab ) ) then
      qual = 0.125D+00
      return
    end if

    expo = top / bot
    qual = base**expo

  end if
!
!  Store QUAL in RWORK
!
  qual = min ( qual, 8.0D+00 )
  qual = max ( qual, 0.125D+00 )

  return
end
subroutine corrector ( df, fpar, fx, ierror, ihold, ipar, iwork, nvar, rwork, &
  wk, xr, lrw, liw, icrit, slname )

!*****************************************************************************80
!
!! CORRECTOR carries out Newton correction of an approximate solution.
!
!  Discussion:
!
!  It is required that the output value of X satisfy this equation to
!  within a certain tolerance.
!
!  Either Newton's method, or the chord Newton method may be used, depending
!  on a user chosen parameter.  In the latter case, the jacobian is only
!  evaluated at the starting point.
!
!  The system of NVAR-1 equations is temporarily augmented by an NVAR-th
!  equation which makes the system square and (presumably) nonsingular.
!  The equation is particularly simple:
!
!    X(IHOLD) = B
!
!  were B is some fixed value. In fact, the value B is set to the input value
!  of X(IHOLD).  This corresponds to simply holding the IHOLD-th entry of
!  X fixed during the iteration, which in turn, corresponds to treating
!  the IHOLD-th entry of X as a "parameter" in terms of which the other
!  entries of X may be solved for.
!
!  The linear system to be solved has the form
!
!    DFA(X,IHOLD) * (-DELX) = FA(X)
!
!  where the Jacobian DFA is augmented with an NVAR-th row containing a
!  1 in the IHOLD-th column, and the vector FA is augmented with an NVAR-th
!  value of 0.
!
!
!  After each Newton step, a decision is made whether to continue the iteration,
!  or to accept the current point, or to reject the entire Newton correction
!  process.  The criteria use the following parameters:
!
!    ABSERR is the user's absolute error tolerance.
!    RWORK(8) is the machine epsilon, also called EPMACH.
!    FMP is an adjustment factor, 2.0 if NCOR = 1, 1.05 otherwise.
!    FNRM is the maximum-norm of the function value of X.
!    FNRML is the maximum-norm of the function value of the previous iterate.
!    MAXCOR is the maximum number of Newton iterations allowed.
!    NCOR is the iteration counter, with NCOR = 0 for the starting point.
!    RELERR is the user's relative error tolerance.
!    STEPX is the maximum-norm of the difference of X and the previous iterate.
!    STEPXL is the previous value of STEPX.
!    X is the current Newton iterate.
!    XNRM is the maximum-norm of X.
!
!  At each step, the current point X will be accepted as the solution to
!  the nonlinear system if any of the following conditions hold:
!
!  Strong acceptance criterion:
!
!  1.  (FNRM.LE.ABSERR) and (STEPX.LE.(ABSERR+RELERR*XNRM))
!
!  Weak acceptance criteria:
!
!  2.  (NCOR.EQ.0) and (FNRM.LE.(.5*ABSERR))
!  3.  (FNRM.LE.8*EPMACH) or (STEPX.LE.8*EPMACH)
!  4.  (NCOR.GE.2) and
!      (FNRM+FNRML).LE.ABSERR and STEPX.LE.8*(ABSERR+RELERR*XNRM)
!  5.  (NCOR.GE.2) and
!      (FNRM.LE.8.0*ABSERR) and (STEPX+STEPXL).LE.(ABSERR+RELERR*XNRM)
!
!  The Newton iteration process is to be aborted if any of the following
!  criteria hold:
!
!  1.  FNRM.GT.(FMP*FNRML+ABSERR)
!  2.  (NCOR.GE.2) AND (ICRIT.EQ.0) and (STEPX.GT.(FMP*STEPXL+ABSERR))
!  3.  NCOR.GE.MAXCOR
!
!  Error conditions returned in IERROR:
!
!  0, no errors detected.
!  1, data or storage error.
!  2, error condition returned from user function or jacobian routine.
!  3, error condition returned from the solver routine.
!  4, a correction step was rejected (function norm increased or
!     size of correction increased)
!  5, too many correction steps were taken without convergence.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
  external df
  external fx
  external slname

  integer liw
  integer lrw
  integer nvar

  double precision abserr
  double precision dets
  double precision fmp
  double precision fnrm
  double precision fnrml
  double precision fpar(*)
  integer i
  integer icrit
  integer ierror
  integer ifmax
  integer ihold
  integer ipar(*)
  integer idamax
  integer iwork(liw)
  integer iwrite
  integer ixmax
  integer j
  integer job
  integer ksmax
  integer maxcor
  integer maxnew
  integer modnew
  double precision relerr
  double precision rwork(lrw)
  double precision dnrm2
  double precision stepxl
  double precision tlstep
  double precision wk(nvar)
  double precision xnrm
  double precision xr(nvar)
  double precision xvalue
!
!  Initialize.
!
  abserr = rwork(1)
  relerr = rwork(2)
  modnew = iwork(4)
  iwrite = iwork(7)
  maxcor = iwork(17)
  ierror = 0

  iwork(28) = 0

  if ( modnew == 0 ) then
    maxnew = maxcor
  else
    maxnew = 2 * maxcor
  end if

  fmp = 2.0D+00
  rwork(9) = 0.0D+00
  xvalue = xr(ihold)
!
!  Get the function value at the starting point.
!
  call fx ( nvar, fpar, ipar, xr, wk )
  iwork(22) = iwork(22)+1

  wk(nvar) = xr(ihold)-xvalue
  ifmax = idamax(nvar,wk,1)
  fnrm = dnrm2(nvar,wk,1)
  ixmax = idamax(nvar,xr,1)
  xnrm = dnrm2(nvar,xr,1)

  if ( iwrite>= 2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Step      FX             X             DX'
    write ( *, * ) ' '
    write(*,'(1x,i4,2g14.6)')iwork(28),fnrm,xnrm
    write(*,'(12x,i3,11x,i3)')ifmax,ixmax
  end if

  if ( 2.0D+00 * fnrm <= abserr ) then
    return
  end if
!
!  Carry out up to MAXNEW Newton corrections.
!
  do i = 1, maxnew

    iwork(28) = i
!
!  Decide whether it is desired to evaluate the jacobian on this step.
!
    if ( ( i /= 1 .and. i /= maxcor .and. modnew == 1 ) .or. modnew == 2 ) then
      job = 1
    else
      job = 0
    end if
!
!  Solve FPRIME * WK = FX.
!
    call slname ( dets, fx, df, fpar, ierror, ihold, ipar, iwork, liw, &
      job, nvar, rwork, lrw, xr, wk )

    if ( ierror/= 0 ) then
      if ( iwrite>= 1 ) then
        write ( *, * ) 'CORRECTOR - Fatal error!'
        write ( *, * ) '  Solver returned IERROR = ',ierror
      end if
      return
    end if
!
!  Subtract WK from XR to get the next iterate.
!
    xr(1:nvar) = xr(1:nvar) - wk(1:nvar)

    stepxl = rwork(9)
    ksmax = idamax ( nvar, wk, 1 )
    rwork(9) = abs ( wk(ksmax) )
    ixmax = idamax ( nvar, xr, 1 )
    xnrm = dnrm2 ( nvar, xr, 1 )
!
!  Compute function value at new iterate and take its norm.
!
    call fx ( nvar, fpar, ipar, xr, wk )
    iwork(22) = iwork(22)+1

    wk(nvar) = xr(ihold) - xvalue
    fnrml = fnrm
    ifmax = idamax ( nvar, wk, 1 )
    fnrm = dnrm2 ( nvar, wk, 1 )

    if ( iwrite>= 2) then
      write(*,'(1x,4x,28x,g14.6)')rwork(9)
      write(*,'(1x,4x,28x,7x,i3)')ksmax
      write(*,'(1x,i4,2g14.6)')iwork(28),fnrm,xnrm
      write(*,'(12x,i3,11x,i3)')ifmax,ixmax
    end if
!
!  Check for strong acceptance of function and stepsize.
!
    tlstep = abserr+relerr*xnrm

    if ( fnrm <= abserr .and. rwork(9) <= tlstep ) then
      return
    end if
!
!  Check for weak acceptance of function and stepsize.
!
    if ( fnrm <= 8.0D+00*rwork(8).or.rwork(9)<=8.0D+00*rwork(8)) then
      return
    end if

    if ( iwork(28)>1 ) then
      if ( (fnrm+fnrml) <= abserr.and.rwork(9)<=8.0D+00*tlstep)return
      if ( fnrm <= 8.0D+00*abserr.and.(rwork(9)+stepxl)<=tlstep)return
    end if
!
!  Decide if iteration should be aborted
!
    if ( iwork(28)>1 ) then
      if ( icrit<1.and.rwork(9)>(fmp*stepxl+abserr) ) then
        ierror = 4
        if ( iwrite>= 2 ) then
          write ( *, * ) 'CORRECTOR - Warning!'
          write ( *, * ) '  The correction DX is not decreasing.'
        end if
        return
      end if
    end if

    if ( icrit<2.and.fnrm>(fmp*fnrml+abserr)) then

       ierror = 4

       if ( iwrite>= 2 ) then
         write ( *, * ) 'CORRECTOR - Warning!'
         write ( *, * ) '  The residual FX is not decreasing.'
       end if

       return

     end if

    fmp = 1.05D+00

  end do
!
!  Reached maximum number of steps without acceptance or rejection.
!
  ierror = 5

  if ( iwrite>= 2 ) then
    write ( *, * ) 'CORRECTOR - Warning!'
    write ( *, * ) '  Convergence is too slow.'
  end if

  return
end
subroutine dgb_det ( a, lda, n, ml, mu, ipivot, det )

!*****************************************************************************80
!
!! DGB_DET computes the determinant of a band matrix factored by SGB_FA.
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision A(LDA,N), the band matrix, as factored by DGB_FA.
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, integer IPIVOT(N), the pivot vector, as computed by DGB_FA.
!
!    Output, double precision DET, the determinant of the matrix.
!
  integer lda
  integer n

  double precision a(lda,n)
  double precision det
  integer i
  integer ipivot(n)
  integer ml
  integer mu
  integer mband

  mband = ml + mu + 1

  det = 1.0D+00

  do i = 1, n
    if ( ipivot(i) /= i ) then
      det = - det
    end if
  end do

  det = det * product ( a(mband,1:n) )

  return
end
subroutine dgb_jac ( eps, fcol, fpar, fprime, frow, fx, ipar, ipc, iwork, &
  jac, liw, nband, nvar, x )

!*****************************************************************************80
!
!! DGB_JAC approximates a banded jacobian matrix.
!
!  Discussion:
!
!    DGB_JAC estimates the jacobian matrix FPRIME of the function FX,
!    using forward or central finite differences.  DGB_JAC is called by
!    DGB_SLV when the user has specified the jacobian option as 1 or 2.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real EPS, a tolerance used for shifting the X values.
!    A value of the square root of the machine precision is
!    usually appropriate.
!
!    Output, real FCOL(NEQN), the last column of the approximate
!    jacobian, which is allowed to be "full".  This comprises
!    matrix entries FPRIME(1,NVAR) through FPRIME(NEQN,NVAR).
!
!    Input/output, real FPAR(*), user parameter vector, not
!    touched by this routine, but passed on to user routines.
!
!    Output, real FPRIME(NBAND,NEQN), is the array into which the
!    the banded portion of the computed jacobian will be stored.
!    The LAPACK general band format is used, assigning entry (I,J)
!    to FPRIME(I-J+ML+MU+1,J), where ML and MU are the lower and
!    upper half bandwidths respectively.
!
!    Output, real FROW(NVAR), storage for the last (augmenting) row
!    of the jacobian, which will be all zero except for a 1 in
!    location IPC.
!
!    Input, external FX, the name of the routine which evaluates the
!    function.
!
!    FX computes the value of the nonlinear function.  This name
!    must be declared external in the calling program.  FX should
!    evaluate the NVAR-1 function components at the input point X,
!    and store the result in the vector FVEC.  An augmenting
!    equation will be stored in entry NVAR of FVEC by the PITCON
!    program.
!
!    FX should have the following form:
!
!    subroutine fx ( nvar, fpar, ipar, x, fvec )
!
!      Input, integer NVAR, number of variables.
!
!      Input/output, real FPAR(*), array of user parameters.
!
!      Input/output, integer IPAR(*), array of user parameters.
!
!      Input, real X(NVAR), the point at which function evaluation is required.
!
!      Output, real FVEC(NVAR), the value of the function at point
!      X.  Only the first NVAR-1 entries of FVEC are to be set by
!      the routine.  PITCON sets the final value itself.
!
!    Input, integer IPAR(*), a user parameter vector passed to FX.
!    However, because this is a problem with a banded jacobian, entries
!    IPAR(1) and IPAR(2) are read by this routine.  IPAR(1) contains
!    ML, the lower half bandwidth of the jacobian, and IPAR(2) contains
!    MU, the upper half bandwidth of the jacobian.
!
!    Input, integer IPC, the index of the current continuation parameter,
!    which is needed to determine the form of FROW.
!
!    Input, integer IWORK(LIW), work and statistics vector.  Only
!    required here so that we can count the number of function
!    evaluations.
!
!    Input, integer JAC, the user requested jacobian option.  For
!    our purposes, the only two values of interest are:
!
!      1 = estimate jacobian with forward differences,
!      2 = estimate jacobian with central differences (twice the work)
!
!    Input, integer LIW, the dimension of IWORK.
!
!    Input, integer NBAND, the first dimension of the jacobian matrix
!    FPRIME, NBAND = ML+MU+1.
!
!    Input, integer NVAR, the number of variables.
!
!    Input, real X(NVAR), the point at which the jacobian is desired.
!
  external fx

  integer liw
  integer nband
  integer nvar

  double precision eps
  double precision fcol(nvar-1)
  double precision fpar(*)
  double precision fprime(nband,nvar-1)
  double precision frow(nvar)
  integer i
  integer ihi
  integer ilo
  integer ipar(2)
  integer ipc
  integer iwork(liw)
  integer j
  integer jac
  integer kcall
  integer mband
  integer ml
  integer mu
  double precision skale
  double precision x(nvar)
  double precision xjac
  double precision xtemp(nvar)
  double precision work1(nvar)
  double precision work2(nvar)

  ml = ipar(1)
  mu = ipar(2)
  mband = ml + mu + 1

  if ( jac == 1 ) then
    call fx ( nvar, fpar, ipar, x, work2 )
    iwork(22) = iwork(22)+1
  end if

  if ( jac == 2 ) then
    xjac = 2.0D+00
  else
    xjac = 1.0D+00
  end if

  do kcall = 1, mband

    xtemp(1:nvar) = x(1:nvar)

    do j = kcall, nvar - 1, mband
      xtemp(j) = x(j) + eps * ( 1.0D+00 + abs ( x(j) ) )
    end do

    call fx ( nvar, fpar, ipar, xtemp, work1 )
    iwork(22) = iwork(22)+1

    if ( jac == 2 ) then

      xtemp(1:nvar) = x(1:nvar)

      do j = kcall, nvar - 1, mband
        xtemp(j) = x(j) - eps * ( 1.0D+00 + abs ( x(j) ) )
      end do

      call fx ( nvar, fpar, ipar, xtemp, work2 )
      iwork(22) = iwork(22)+1

    end if

    do j = kcall, nvar - 1, mband

      ilo = max ( 1, j - mu )
      ihi = min ( nvar - 1, j + ml )

      skale = 1.0D+00 / (xjac*eps*(1.0D+00+abs(x(j))))

      do i = ilo, ihi
        fprime(i-j+ml+mu+1,j) = fprime(i-j+ml+mu+1,j) &
          + ( work1(i) - work2(i) ) * skale
      end do

    end do

  end do
!
!  Compute last column of jacobian, rows 1 to NEQN
!
  xtemp(1:nvar) = x(1:nvar)

  xtemp(nvar) = x(nvar) + eps*(1.0D+00+abs(x(nvar)))

  call fx ( nvar, fpar, ipar, xtemp, work1 )
  iwork(22) = iwork(22)+1

  if ( jac == 2 ) then

    xtemp(nvar) = x(nvar) - eps*(1.0D+00+abs(x(nvar)))
    call fx ( nvar, fpar, ipar, xtemp, work2 )
    iwork(22) = iwork(22)+1

  end if

  do i = 1, nvar - 1
    fcol(i) = fcol(i) + ( work1(i) - work2(i) ) &
      / ( xjac * eps * ( 1.0D+00 + abs ( x(nvar) ) ) )
  end do
!
!  Do the last row, J = 1,NVAR
!
  frow(ipc) = frow(ipc) + 1.0D+00

  return
end
subroutine dgb_slv(dets,fx,df,fpar,ierror,ipc,ipar,iwork,liw,job, &
  nvar,rwork,lrw,x,y)

!*****************************************************************************80
!
!! DGB_SLV solves a dense banded linear system.
!
!  Discussion:
!
!    The linear system has the form
!
!            ( DFDY | DFDZ )
!      A*X = (-------------) * X = B
!            (   E(IPC)    )
!
!    where
!
!      B is a given vector of length NVAR,
!      DFDY is "logically" an NVAR-1 by NVAR-1 matrix, with band structure,
!      DFDZ is an NVAR-1 vector,
!      E(IPC) is an NVAR vector whose only nonzero entry is a 1 in position IPC.
!
!    DFDY is actually stored compactly, using LAPACK general band storage,
!
!    DFDY and DFDZ represent the jacobian of an NVAR-1 dimensional function
!    of NVAR variables, and E(IPC) is the augmenting row, determined by
!    the choice of "continuation parameter" IPC.
!
!
!    DGB_SLV factors and solves the linear system A*x = b, taking advantage
!    of the bandedness of the DFDY subsystem, which would be lost if the
!    full system was factored and solved directly.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real DETS, the sign of the determinant of the full matrix A.
!
!    Input, external FX, the name of the routine which evaluates the function.
!
!    FX computes the value of the nonlinear function.  This name must be
!    declared external in the calling program.  FX should evaluate the
!    NVAR-1 function components at the input point X, and store the result
!    in the vector FVEC.  An augmenting equation will be stored in entry
!    NVAR of FVEC by the PITCON program.
!
!    FX should have the following form:
!
!    subroutine fx ( nvar, fpar, ipar, x, fvec )
!
!      Input, integer NVAR, number of variables.
!
!      Input/output, real FPAR(*), user parameters.
!
!      Input/output, integer IPAR(*), user parameters.
!
!      Input, real X(NVAR), the point at which function evaluation is required.
!
!      Output, real FVEC(NVAR), value of the function at X.  Only the first
!      NVAR-1 entries of FVEC are to be set by the routine.  PITCON sets the
!      final value itself.
!
!    Input, external FP, the name of the user supplied routine which
!    computes the jacobian (DFDY, DFDZ) of the NVAR-1 nonlinear functions.
!    Jacobian entries for the NVAR-th, augmenting, function are inserted by
!    this routine.  Because of the special banded storage arrangement, care
!    must be taken to store the information properly.
!
!    Let MU and ML be the upper and lower bandwidths of DFDY.
!    Set NBAND = 2*ML+MU+1.
!
!    Then the index K of RWORK into which we want to store entry
!    (I,J) of DFDY is determined as follows:
!
!      K = I+J*(NBAND-1)-ML
!
!    Here we are assuming that I-ML <= J <= I+MU.  For any other values
!    of I and J, the corresponding entry of DFDY is assumed to be zero.
!
!    There are NVAR-1 entries in DFDZ, and the I-th entry of this
!    vector is stored in RWORK(K), where
!
!      K = (NVAR-1)*NBAND+I
!
!    DF must be a subroutine of the form:
!
!    subroutine df ( nvar, fpar, ipar, x, rwork )
!
!      integer nvar
!
!      real fpar(*)
!      integer ipar(*)
!      real rwork(*)
!      real x(nvar)
!
!      ml = ipar(1)
!      mu = ipar(2)
!      nband = 2*ml+mu+1
!
!      do i = 1,nvar-1
!        do j = max(i-ml,1), min(i+mu,nvar-1)
!          k = i+j*(nband-1)-ml
!          rwork(k) = d f(i)/d x(j)
!        end do
!      end do
!
!      do i = 1,nvar-1
!        k = (nvar-1)*nband+i
!        rwork(k) = d f(i)/d x(nvar)
!      end do
!
!      return
!    end
!
!    Input/output, real FPAR(*), a user defined parameter array.
!
!    Output, integer IERROR, the error return flag.
!    0, No errors were detected.
!    1, data or storage error, including illegal values for NVAR,
!      IPC, MU, ML, or insufficient storage in RWORK or IWORK.
!    2, The user set a nonzero error return in DF.
!    3, The matrix DFDY or A is numerically singular.
!      In some cases, a different choice of IPC could rectify this problem.
!
!    Input, integer IPC, the continuation parameter.
!    IPC determines the form of the final, augmenting row of the
!    jacobian matrix.
!
!    If IPC is equal to NVAR, then the linear system can easily
!    be solved using a standard band matrix solver, once we
!    have solved for X(NVAR).
!
!    But in the general case when IPC is not NVAR, we have to
!    do some work to modify the system, and still be able to
!    use a band solver.
!
!    Input, integer IPAR(*).
!    IPAR(1) = ML, the lower bandwidth of DFDY,
!    IPAR(2) = MU, the upper bandwidth of DFDY.
!    The other entries of IPAR are not referenced by PITCON, and
!    may be used to pass information between the calling program
!    and the user routines FX and FP.
!
!    Workspace, IWORK(*), a work array used by the continuation code
!    to store statistics, pointers, and the pivot vector for the
!    linear equation solver.
!
!    IWORK(13) stores the address of the first entry in IWORK
!    used for pivoting.
!
!    IWORK(15) stores the address of the first entry in RWORK
!    used to store the Jacobian.
!
!    IWORK(20) counts the number of matrix factorizations.
!
!    IWORK(21) counts the number of linear system back-solves.
!
!    Input, integer LIW, the dimension of IWORK.
!
!    Input, integer JOB, controls the action of the routine:
!
!    0, evaluate jacobian, decompose jacobian, compute determinant,
!    and solve linear system.
!
!    1, solve a system with a new right hand side, and a previously
!    factored jacobian.
!
!    2, evaluate jacobian, decompose jacobian, and compute determinant.
!
!    3, Check jacobian matrix.  Call user jacobian routine,
!    multiply by -1.0, add finite difference jacobian,
!    print largest entry.
!
!    Input, integer NVAR, the number of variables.  The dimension of X.
!
!    Workspace, real RWORK(LRW).  RWORK contains workspace used
!    for storing the jacobian, as well as various vectors and
!    scalars.  The IWORK array contains pointers to the beginning
!    locations of some of these objects.
!
!    Input, integer LRW, the dimension of RWORK.
!
!    Input, real X(NVAR), the point at which the jacobian is to be evaluated.
!
!    Input/output, real Y(NVAR).  Right hand side/solution vector.
!    On input, Y contains the right hand side.
!    On output, Y contains the solution to that same linear system.
!
  external df
  external fx

  integer liw
  integer lrw
  integer nvar

  double precision ak
  double precision det
  double precision dets
  double precision fpar(*)
  integer i
  integer ierror
  integer info
  integer ipar(*)
  integer ipc
  integer irl
  integer iru
  integer idamax
  integer itemp
  integer iwork(liw)
  integer iwrite
  integer j
  integer jac
  integer jack
  integer job
  integer jtemp
  integer k
  integer lda
  integer ldfl
  integer ldfx
  integer ldx
  integer lfxm
  integer lfxp
  integer lilst
  integer lpiv
  integer lrlst
  integer lrowip
  integer mband
  integer ml
  integer mu
  integer nband
  integer ndim
  integer neqn
  integer nswap
  double precision rwork(lrw)
  double precision temp
  double precision x(nvar)
  double precision y(nvar)

  ierror = 0
  iwrite = iwork(7)
  ml = ipar(1)
  mu = ipar(2)
  neqn = nvar-1

  if ( ml<0 ) then
    ierror = 1
    write ( *, * ) 'DGB_SLV - Fatal error!'
    write ( *, * ) '  Illegal lower bandwidth ML = ',ml
    write ( *, * ) '  ML must be at least 0.'
    stop
  else if ( ml>nvar - 2  ) then
    ierror = 1
    write ( *, * ) 'DGB_SLV - Fatal error!'
    write ( *, * ) '  Illegal lower bandwidth ML = ',ml
    write ( *, * ) '  ML must be no more than ',nvar-2
    stop
  else if ( mu<0 ) then
    ierror = 1
    write ( *, * ) 'DGB_SLV - Fatal error!'
    write ( *, * ) '  Illegal upper bandwidth MU = ',mu
    write ( *, * ) '  MU must be at least 0.'
    stop
  else if ( mu>nvar-2 ) then
    ierror = 1
    write ( *, * ) 'DGB_SLV - Fatal error!'
    write ( *, * ) '  Illegal upper bandwidth MU = ',mu
    write ( *, * ) '  MU must be no more than ',nvar-2
    stop
  end if

  lda = 2*ml+mu+1
  mband = ml+mu+1
  nband = mband+ml
  lpiv = iwork(13)
  lilst = lpiv+nvar-2
  ldfx = iwork(15)
  ldfl = ldfx+nband*( nvar - 1 )
  lrowip = ldfl+nvar
  jac = iwork(9)
!
!  Make sure that the jacobian routine is available if the user
!  has specified an option that requires it!
!
  if ( jac/= 0.and. &
       (iwork(1) == (-1).or.iwork(1)==(-2).or. &
        iwork(1) == (-3).or.iwork(1)==(-4).or. &
        iwork(1) == (-9).or.iwork(1)==(-10)) ) then
    ierror = 4
    if ( iwrite>= 1 ) then
      write ( *, * ) 'DGB_SLV - Fatal error!'
      write ( *, * ) '  An option was selected requiring'
      write ( *, * ) '  a user jacobian, but the value of'
      write ( *, * ) '  IWORK(9) indicates no such routine'
      write ( *, * ) '  is available!'
    end if
    return
  end if
!
!  Check that enough storage is available.
!
  if ( jac == 0.and.job/=3 ) then
    lrlst = lrowip+nvar-1
  else
    lfxp = lrowip+nvar
    lfxm = lfxp+nvar
    ldx = lfxm+nvar
    lrlst = ldx+nvar-1
  end if
  ndim = ( nvar - 1 ) * nband+nvar+nvar-1

  if ( lilst>liw ) then
    ierror = 1
    if ( iwrite>= 1 ) then
      write ( *, * ) 'DGB_SLV - Fatal error!'
      write ( *, * ) '  Insufficient integer workspace in IWORK!'
      write ( *, * ) '  Need workspace size    = ',lilst
      write ( *, * ) '  Available workspace LIW = ',liw
    end if
    return
  end if

  if ( lrlst>lrw ) then
    ierror = 1
    if ( iwrite>= 1 ) then
      write ( *, * ) 'DGB_SLV - Fatal error!'
      write ( *, * ) '  Insufficient real workspace in RWORK!'
      write ( *, * ) '  Needed workspace size    = ',lrlst
      write ( *, * ) '  Available workspace LRW =  ',lrw
    end if
    return
  end if
!
!  If JOB = 1, we need to solve a linear system.
!  Two entirely different methods are used to solve a linear system,
!  depending on whether the parameter is the last variable or not.
!
  if ( job == 1 ) then
    if ( ipc == nvar ) then
      go to 70
    else
      go to 40
    end if
  end if
!
!  JOB = 0 OR 2 means we must
!    evaluate the jacobian,
!    factor it,
!    get the sign of the determinant.
!
  do i = 1,ndim
    rwork(ldfx+i-1) = 0.0D+00
  end do
!
!  If a user jacobian routine is available, invoke it.
!
  if ( jac == 0 ) then
    if ( iwork(1)>(-5).or.iwork(1)<(-8) ) then
      call df ( nvar, fpar, ipar, x, rwork(ldfx) )
      iwork(19) = iwork(19)+1
      rwork(lrowip-1+ipc) = 1.0D+00
    end if
  end if
!
!  If we're going to compare the user jacobian to an approximate
!  one, negate the user jacobian.
!
  if ( job == 3 ) then
    if ( iwork(1) <= (-1).and.iwork(1)>=(-4) ) then

      do i = ldfx, ldfx + ndim - 1
        rwork(i) = - rwork(i)
      end do

    end if
  end if
!
!  If we don't have a user jacobian, or we need to compare the
!  user jacobian to an approximate one, get the approximation.
!
  if ( (iwork(1)>= 0.and.(jac == 1.or.jac==2)).or. &
       (iwork(1) <= (-1).and.iwork(1)>=(-8)) ) then

    jack = jac

    if ( job == 3 ) then
      if ( mod(iabs(iwork(1)),2) == 1 ) then
        jack = 1
      else
        jack = 2
      end if
    end if

    call dgb_jac(rwork(18),rwork(ldfl),fpar,rwork(ldfx),rwork(lrowip), &
      fx,ipar,ipc,iwork,jack,liw,nband,nvar,x )

  end if
!
!  If JOB = 3, print out information:
!    IWORK(1) = -1, -2, -5, -6, -9, print largest entry of matrix.
!    IWORK(1) = -3, -4, -7, -8, -10, also print out entire matrix.
!
  if ( job == 3 ) then

    k = idamax(ndim,rwork(ldfx),1)
    ak = rwork(ldfx+k-1)

    if ( k <= (nvar-1)*nband ) then
      j = ((k-1)/nband)+1
      i = k-(j-1)*nband+j-ml-mu-1
    else if ( k <= (nvar-1)*nband+nvar - 1 ) then
      i = k-(nvar-1)*nband
      j = nvar
    else
      i = nvar
      j = k-(nvar-1)*nband-neqn
    end if

    if ( iwork(1) <= (-1).and.iwork(1)>=(-4) ) then
      write ( *, * ) ' '
      write ( *, * ) 'DGB_SLV - Maximum value of FP_Approx(I,J)-FP_User(I,J)'
    else if ( iwork(1) <= (-5).and.iwork(1)>=(-8) ) then
      write ( *, * ) ' '
      write ( *, * ) 'DGB_SLV - Maximum value of finite difference jacobian:'
    else if ( iwork(1) <= (-9).and.iwork(1)>=(-10) ) then
      write ( *, * ) ' '
      write ( *, * ) 'DGB_SLV - Maximum value of user supplied jacobian:'
    end if

    write ( *, * ) ak,' I, J = ',i,j
    write ( *, * ) ' '

    if ( iwork(1) == (-3).or.iwork(1)==(-4).or. &
         iwork(1) == (-7).or.iwork(1)==(-8).or. &
         iwork(1) == (-10) ) then
      if ( iwork(1) == (-3).or.iwork(1)==(-4) ) then
        write ( *, * ) ' '
        write ( *, * ) 'DGB_SLV - Entire difference matrix:'
        write ( *, * ) 'FP_Approx(I,J)-FP_User(I,J)'
        write ( *, * ) ' '
      else if ( iwork(1) == (-7).or.iwork(1)==(-8) ) then
        write ( *, * ) ' '
        write ( *, * ) 'DGB_SLV - Finite difference jacobian:'
        write ( *, * ) ' '
      else if ( iwork(1) == (-10) ) then
        write ( *, * ) ' '
        write ( *, * ) 'DGB_SLV - User supplied jacobian:'
        write ( *, * ) ' '
      end if

      do i = 1,nvar
        do j = 1,nvar
          if ( j == nvar ) then
            k = ldfx-1+(nvar-1)*nband+i
            write ( *, * ) rwork(k),' I, J = ',i,j
          else if ( (j-i <= ml).and.(i-j<=mu) ) then
            k = ldfx-1+i+j*(nband-1)-ml
            write ( *, * ) rwork(k),' I, J = ',i,j
          end if
        end do
        write ( *, * ) ' '
      end do

    end if
    return
  end if

  if ( ipc == nvar)go to 60
!
!  Switch the NVAR-th and IPC-th rows.
!
  irl = max(1,ipc-ml)
  iru = min(nvar-1,ipc+mu)
  nswap = iru+1-irl
  itemp = ldfx-ml-1+ipc+irl*(nband-1)
  jtemp = nband-1

  do i = 0, nswap-1
    temp                     = rwork(lrowip-1+irl+i)
    rwork(lrowip-1+irl+i)    = rwork(itemp+i*(nband-1))
    rwork(itemp+i*(nband-1)) = temp
  end do

!  call dswap(nswap,rwork(lrowip-1+irl),1,rwork(itemp),jtemp)
!
!  Decompose the submatrix.
!
  call dgb_trf(neqn,neqn,ml,mu,rwork(ldfx),lda,iwork(lpiv),info)
  iwork(20) = iwork(20)+1

  if ( info/= 0 ) then
    if ( iwrite>= 1 ) then
      write ( *, * ) 'DGB_SLV - Fatal error!'
      write ( *, * ) '  LU factor routine DGBTRF returns INFO = ',info
      write ( *, * ) '  Current index IPC = ',ipc
    end if
    ierror = 3
    return
  end if
!
!  Compute the determinant.
!
  call dgb_det ( rwork(ldfx), lda, neqn, ml, mu, iwork(lpiv), det )

  dets = 0.0D+00
  if ( det>0.0D+00 ) then
    dets = -1.0D+00
  else if ( det<0.0D+00 ) then
    dets = 1.0D+00
  end if
!
!  Set the right hand side of the auxilliary system to the last
!  column of the jacobian, minus E(IPC).
!
!  Shuffle the IPC-th and NVAR-th entries of this right hand
!  side, to reflect the pivoting of the equations.
!
!  Then solve the system.
!
  temp = rwork(ldfl+ipc-1)-1.0D+00
  rwork(ldfl+ipc-1) = rwork(ldfl+nvar-1)
  rwork(ldfl+nvar-1) = temp
  call dgb_trs ( 'n', neqn, ml, mu, 1, rwork(ldfx), lda, iwork(lpiv), &
    rwork(ldfl), neqn, info )
  iwork(21) = iwork(21)+1

  if ( info/= 0 ) then
    if ( iwrite>= 1 ) then
      write ( *, * ) 'DGB_SLV - Fatal error!'
      write ( *, * ) '  LU backsolve routine DGBTRS returns INFO = ',info
    end if
    ierror = 1
    return
  end if
!
!  Solve for last entry of auxilliary solution.
!
  temp = 0
  do i = 0, nvar - 2
    temp = temp + rwork(lrowip+i) * rwork(ldfl+i)
  end do

  rwork(ldfl+nvar-1) = rwork(ldfl+nvar-1) - temp
!
!  Adjust the sign of the determinant.
!
  if ( 1.0D+00+rwork(ldfl+nvar-1) == 0.0D+00 ) then
    ierror = 3
    if ( iwrite>= 1 ) then
      write ( *, * ) 'DGB_SLV - Algorithm fails, DENOM = 0.0'
    end if
    return
  end if

  if ( 1.0D+00+rwork(ldfl+nvar-1)<0.0D+00)dets = -dets
  if ( job == 2)return
!
!  Solve the system.
!
40    continue

  if ( ipc == nvar)go to 70
!
!  Modify right hand side of main system.
!
  temp    = y(ipc)
  y(ipc)  = y(nvar)
  y(nvar) = temp
!
!  Solve subsystem.
!
  call dgb_trs ( 'n', neqn, ml, mu, 1, rwork(ldfx), lda, iwork(lpiv), y, &
    neqn, info )

  iwork(21) = iwork(21)+1

  if ( info/= 0 ) then
    if ( iwrite>= 1 ) then
      write ( *, * ) 'DGB_SLV - Fatal error!'
      write ( *, * ) '  DGBTRS returns INFO = ',info
    end if
    ierror = 1
    return
  end if
!
!  Solve for last entry of main solution.
!
  temp = 0.0D+00
  do i = 1, nvar - 1
    temp = temp + rwork(lrowip+i-1) * y(i)
  end do

  y(nvar) = y(nvar) - temp
!
!  Correct the main solution, using a multiple of the "subsolution".
!
  do i = 1, nvar
    y(i) = y(i) - y(nvar) * rwork(ldfl+i-1) / ( 1.0D+00 + rwork(ldfl+nvar-1) )
  end do

  return
!
!  Factor the matrix for the special case of IPC = NVAR.
!
   60 continue

  call dgb_trf(neqn,neqn,ml,mu,rwork(ldfx),lda,iwork(lpiv),info)

  iwork(20) = iwork(20)+1

  if ( info/= 0 ) then
    if ( iwrite>= 1 ) then
      write ( *, * ) 'DGB_SLV - Fatal error!'
      write ( *, * ) '  DGBTRF returns INFO = ',info
      write ( *, * ) '  Current index IPC = ',ipc
    end if
    ierror = 3
    return
  end if
!
!  Compute the determinant.
!
  call dgb_det ( rwork(ldfx), lda, neqn, ml, mu, iwork(lpiv), det )

  dets = 0.0D+00
  if ( det<0.0D+00 ) then
    dets = -1.0D+00
  else if ( det>0.0D+00 ) then
    dets = 1.0D+00
  end if

  if ( job == 2)return
!
!  Solve the linear system for the special case where IPC = NVAR.
!
   70 continue

  do i = 1, nvar - 1
    y(i) = y(i) - y(nvar) * rwork(ldfl+i-1)
  end do

  call dgb_trs('n',neqn,ml,mu,1,rwork(ldfx),lda,iwork(lpiv),y,neqn,info)
  iwork(21) = iwork(21)+1

  return
end
subroutine dge_det ( a, lda, n, ipivot, det )

!*****************************************************************************80
!
!! DGE_DET computes the determinant of a matrix factored by DGE_TRF.
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision A(LDA,N), the LU factors computed by DGE_TRF.
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer IPIVOT(N), as computed by DGE_TRF.
!
!    Output, double precision DET, the determinant of the matrix.
!
  integer lda
  integer n

  double precision a(lda,n)
  double precision det
  integer i
  integer ipivot(n)

  det = 1.0D+00

  do i = 1, n
    det = det * a(i,i)
  end do

  do i = 1, n
    if ( ipivot(i) /= i ) then
      det = - det
    end if
  end do

  return
end
subroutine dge_jac ( eps, fpar, fprime, fx, ipar, ipc, iwork, jac, liw, nvar, &
  x, work1, work2 )

!*****************************************************************************80
!
!! DGE_JAC approximates a dense jacobian matrix.
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real EPS, a tolerance to be used for shifting the X
!    values during the finite differencing.  No single value of EPS
!    will be reliable for all vectors X and functions FX.  Values of
!    EPS have typically been chosen between SQRT(EPSMCH) and
!    SQRT(SQRT(EPSMCH)) where EPSMCH is the machine tolerance.
!
!    Input, real FPAR(*), a real vector, available for the user
!    to communicate with the FX routine.
!
!    Output, real FPRIME(NVAR,NVAR), the array into which the
!    jacobian will be stored, including the augmenting last row.
!
!    Input, external FX, the name of the user supplied routine
!    that defines the nonlinear function, and which has the form:
!
!    subroutine fx ( nvar, fpar, ipar, x, f )
!
!      integer nvar
!
!      real f(nvar)
!      real fpar(*)
!      integer ipar(*)
!      real x(nvar)
!
!      do i = 1,nvar-1
!        f(i) = function(i)(x)
!      end do
!
!      return
!    end
!
!    Input, integer IPAR(*), an integer vector, available for the user
!    to communicate with the FX routine.
!
!    Input, integer IPC, the index of the current continuation parameter,
!    which determines the location of the "1" entry in the last row
!    of the jacobian.
!
!    Input/output, integer IWORK(LIW), an array containing pointers
!    and counters used by the continuation code.  In particular,
!    IWORK(22) contains a count of the number of function evaluations.
!
!    Input, integer JAC, the jacobian option.
!    0, the user has supplied the jacobian routine DF.
!    1, the program must estimate the jacobian using forward differences.
!    2, the program must estimate the jacobian using central differences.
!
!    Input, integer LIW, the dimension of IWORK.
!
!    Input, integer NVAR, the number of variables.
!
!    Input, real X(NVAR), the point at which the jacobian is to be estimated.
!
!    Workspace, real WORK1(NVAR), WORK2(NVAR).
!
  external fx

  integer liw
  integer nvar

  double precision delm
  double precision delp
  double precision eps
  double precision fpar(*)
  double precision fprime(nvar,nvar)
  integer i
  integer ipar(*)
  integer ipc
  integer iwork(liw)
  integer j
  integer jac
  double precision x(nvar)
  double precision xsave
  double precision work1(nvar)
  double precision work2(nvar)
!
!  If we are using forward differences, then evaluate the function F
!  at the base point X, and save the value for use in the difference
!  quotient.
!
  if ( jac == 1 ) then
    call fx ( nvar, fpar, ipar, x, work2 )
    iwork(22) = iwork(22) + 1
    delm = 0.0D+00
  end if
!
!  Increment each variable X(J) by a small amount, and evaluate the
!  function there.
!
  do j = 1,nvar
    xsave = x(j)
    delp = eps*(1.0D+00+abs(x(j)))
    x(j) = x(j)+delp
    call fx ( nvar, fpar, ipar, x, work1 )
    iwork(22) = iwork(22)+1
!
!  For central difference approximations, decrement each variable X(J)
!  by a small amount, and evaluate the function there.
!
    if ( jac == 2 ) then

      delm = -delp
      x(j) = xsave + delm

      call fx ( nvar, fpar, ipar, x, work2 )

      iwork(22) = iwork(22) + 1

    end if

    x(j) = xsave
!
!  Compute DFDX(*,J) = (F(X+)-F(X-))/DELX.
!
!  Note that we contrive to ADD this quantity to DFDX, without
!  overwriting anything that may already be in DFDX.
!
!  This makes it possible to use this routine also to check
!  the accuracy of the user's jacobian routine.
!
    do i = 1, nvar - 1
      fprime(i,j) = fprime(i,j) + ( work1(i) - work2(i) ) / ( delp - delm )
    end do

  end do

  fprime(nvar,ipc) = fprime(nvar,ipc) + 1.0D+00

  return
end
subroutine dge_slv ( dets, fx, df, fpar, ierror, ipc, ipar, iwork, liw, job, &
  nvar, rwork, lrw, x, y )

!*****************************************************************************80
!
!! DGE_SLV solves the NVAR by NVAR dense linear system
!
!         (  DF(X)  )
!  (1)    (---------) * Y  = Y
!         (  E(IPC) )
!
!  where the NVAR-1 by NVAR submatrix DF(X) is the jacobian of the
!  nonlinear function F(X), or an approximation thereto, and E(IPC)
!  is a row vector of NVAR entries, consisting of all zeroes except
!  for a 1 in column IPC.
!
!  As a special application, DGE_SLV may also be called to compare
!  the user supplied jacobian to a finite difference approximation.
!
!  DGE_SLV is used when the jacobian DF(X) is assumed to be dense.
!  In this case, the corresponding matrix is stored in NVAR*NVAR
!  consecutive entries in the RWORK array, RWORK(LRBEG) through
!  RWORK(LRBEG+NVAR*NVAR-1), with the value of LRBEG stored in
!  IWORK(15).
!
!  If DGE_SLV is required to factor the matrix, the factored version will
!  overwrite the original matrix, using the same storage area.  As part of
!  the factorization, DGE_SLV will also need a vector of length NVAR
!  to store the pivot information.  This is done in the tail end of the
!  IWORK array.  The first entry of IWORK devoted to this is recorded
!  in IWORK(13).  In other words, the pivot information starts in
!  IWORK(IWORK(13)).
!
!  If the jacobian of the nonlinear function is dense, the user should
!  specify DGE_SLV as the external solver argument SLNAME.  The user
!  should either request a finite difference jacobian be calculated,
!  or provide a jacobian routine of the following form:
!
!      subroutine df ( nvar, fpar, ipar, x, a )
!
!        integer nvar
!
!        real a(nvar,nvar)
!        real fpar(*)
!        integer ipar(*)
!        real x(nvar)
!
!        do i = 1,nvar-1
!          do j = 1,nvar
!            a(i,j) = df(i)/dx(j)
!          end do
!        end do
!
!        return
!      end
!
!  where "DF(I)/DX(J)" denotes the derivative of the I-th component of the
!  nonlinear function F(X) with respect to the J-th component of X.
!  The user need not supply the NVAR-th, augmenting, row of the matrix, since
!  this is taken care of automatically.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real DETS, the sign of the determinant of the matrix
!    in equation (1).
!
!    Input, external FX, the name of the user supplied routine
!    which evaluates the nonlinear function.
!
!    Input, external DF, the name of the user supplied routine
!    which evaluates the jacobian matrix.
!
!    Input, real FPAR(*), a real array which allows the user to pass
!    double precision parameters "through" PITCON to FX and DF.
!
!    Output, integer IERROR, error flag.
!    0, normal return.
!    1, data or storage error.
!    2, error returned by the derivative routine DF.
!    3, the matrix containing the augmented jacobian is singular.
!
!    Input, integer IPC, the index of the continuation parameter,
!    which determines the location of the "1" entry in the final
!    row of the jacobian.
!
!    Input, integer IPAR(*), an integer array available for the user
!    to transmit information to DF and FX.
!
!    Workspace, integer IWORK(LIW), an array containing pointers,
!    counters, and space for a pivot array.
!
!    Input, integer LIW, the dimension of IWORK.
!
!    Input, integer JOB, the action switch.
!
!    0, Evaluate jacobian, factor jacobian, compute determinant,
!    and solve linear system.
!
!    1, solve linear system with previously factored jacobian.
!
!    2, Evaluate jacobian, factor jacobian, compute determinant.
!
!    3, Check jacobian matrix.  Call user jacobian routine,
!    multiply by -1.0, add finite difference jacobian,
!    print largest entry.
!
!    Input, integer NVAR, the number of variables.
!
!    Input/output, real RWORK(*).  RWORK contains various scalars,
!    vectors, and the jacobian matrix.  On input, RWORK may contain
!    a previously factored jacobian matrix.  On output, RWORK
!    will generally contain the factored jacobian matrix.
!
!    Input, integer LRW, the dimension of RWORK.
!
!    Input, real X(NVAR), the point at which the jacobian is to be evaluated.
!
!    Input/output, real Y(NVAR).  On input, Y may contain the right
!    hand side of a linear system to be solved, in which case Y
!    will contain the solution of that linear system on output.
!
  external df
  external fx

  integer liw
  integer lrw
  integer nvar

  double precision det
  double precision dets
  double precision fpar(*)
  integer i
  integer ierror
  integer info
  integer ipar(*)
  integer ipc
  integer idamax
  integer iwork(liw)
  integer iwrite
  integer j
  integer jac
  integer jack
  integer job
  integer k
  integer ldf
  integer lfxm
  integer lfxp
  integer lilst
  integer lpiv
  integer lrlst
  integer ndim
  double precision rwork(lrw)
  double precision x(nvar)
  double precision y(nvar)

  ierror = 0
  iwrite = iwork(7)
  lpiv = iwork(13)
  ldf = iwork(15)
  jac = iwork(9)
!
!  Make sure that the jacobian routine is available if the user
!  has specified an option that requires it!
!
  if ( jac/= 0.and. &
       (iwork(1) == (-1).or.iwork(1)==(-2).or. &
        iwork(1) == (-3).or.iwork(1)==(-4).or. &
        iwork(1) == (-9).or.iwork(1)==(-10)) ) then
    ierror = 4
    if ( iwrite>= 1 ) then
      write ( *, * ) 'DGE_SLV - Fatal error!'
      write ( *, * ) 'An option was selected requiring'
      write ( *, * ) 'a user jacobian, but the value of'
      write ( *, * ) 'IWORK(9) indicates no such routine'
      write ( *, * ) 'is available!'
    end if
    return
  end if
!
!  Check the amount of storage available
!
  lilst = lpiv+nvar-1

  if ( jac == 0.and.job/=3 ) then
    lrlst = ldf-1+nvar*nvar
  else
    lrlst = ldf+nvar*nvar+2*nvar-1
  end if

  if ( lilst>liw ) then
    ierror = 1
    if ( iwrite>= 1 ) then
      write ( *, * ) 'DGE_SLV - Fatal error!'
      write ( *, * ) 'Insufficient integer workspace in IWORK!'
      write ( *, * ) 'Need workspace size    = ',lilst
      write ( *, * ) 'Available workspace LIW = ',liw
    end if
    return
  end if

  if ( lrlst>lrw ) then
    ierror = 1
    if ( iwrite >= 1 ) then
      write ( *, * ) 'DGE_SLV - Fatal error!'
      write ( *, * ) '  Insufficient real workspace in RWORK!'
      write ( *, * ) 'Need workspace size    = ',lrlst
      write ( *, * ) 'Available workspace LRW = ',lrw
    end if
    return
  end if
!
!  JOB = 1.
!  A linear system is to be solved.  The matrix has already been factored,
!  and the right hand side is supplied.  Return the solution.
!
  if ( job == 1 ) then

    call dge_trs ( 'n', nvar, 1, rwork(ldf), nvar, iwork(lpiv), y, nvar, info )

    iwork(21) = iwork(21) + 1

    if ( info /= 0 ) then
      if ( iwrite >= 1 ) then
        write ( *, * ) 'DGE_SLV - DGE_TRS returns INFO = ', info
      end if
      ierror = 3
    end if

    return

  end if

  ndim = nvar * nvar
  do i = 1, ndim
    rwork(ldf+i-1) = 0.0D+00
  end do
!
!  If a user jacobian routine is available, invoke it.
!
  if ( jac == 0 ) then
    if ( iwork(1)>(-5).or.iwork(1)<(-8) ) then
      call df ( nvar, fpar, ipar, x, rwork(ldf) )
      iwork(19) = iwork(19)+1
      rwork(ldf+ipc*nvar-1) = 1.0D+00
    end if
  end if
!
!  If we're going to compare the user jacobian to an approximate one,
!  negate the user jacobian now.
!
  if ( job == 3 ) then
    if ( iwork(1) <= (-1).and.iwork(1)>=(-4) ) then

      do i = ldf, ldf + ndim - 1
        rwork(i) = - rwork(i)
      end do

    end if
  end if
!
!  If we don't have a user jacobian, or we need to compare the
!  user jacobian to an approximate one, get the approximation.
!
  if ( (iwork(1)>= 0.and.(jac == 1.or.jac==2)).or. &
       (iwork(1) <= (-1).and.iwork(1)>=(-8)) ) then

    jack = jac

    if ( job == 3 ) then
      if ( mod(iabs(iwork(1)),2) == 1 ) then
        jack = 1
      else
        jack = 2
      end if
    end if

    lfxp = ldf+nvar*nvar
    lfxm = ldf+nvar*nvar+nvar

    call dge_jac ( rwork(18), fpar, rwork(ldf), fx, ipar, &
      ipc, iwork, jack, liw, nvar, x, rwork(lfxp), rwork(lfxm) )

  end if
!
!  If JOB = 3, print out information:
!    IWORK(1) = -1, -2, -5, -6, -9, print largest entry of matrix.
!    IWORK(1) = -3, -4, -7, -8, -10, also print out entire matrix.
!
  if ( job == 3 ) then

    k = idamax(ndim,rwork(ldf),1)
    i = mod(k-1,nvar)+1
    j = (k-i)/nvar+1

    if ( iwork(1) <= (-1).and.iwork(1)>=(-4) ) then
      write ( *, * ) ' '
      write ( *, * ) 'DGE_SLV - Maximum value of FP_Approx(I,J)-FP_User(I,J)'
    else if ( iwork(1) <= (-5).and.iwork(1)>=(-8) ) then
      write ( *, * ) ' '
      write ( *, * ) 'DGE_SLV - Maximum value of finite difference jacobian:'
    else if ( iwork(1) <= (-9).and.iwork(1)>=(-10) ) then
      write ( *, * ) ' '
      write ( *, * ) 'DGE_SLV - Maximum value of user supplied jacobian:'
    end if

    write ( *, * ) rwork(ldf+k-1),' I, J = ',i,j
    write ( *, * ) ' '

    if ( iwork(1) == (-3).or.iwork(1)==(-4).or. &
         iwork(1) == (-7).or.iwork(1)==(-8).or. &
         iwork(1) == (-10) ) then

      if ( iwork(1) == (-3).or.iwork(1)==(-4) ) then
        write ( *, * ) ' '
        write ( *, * ) 'DGE_SLV - Entire difference matrix:'
        write ( *, * ) '         FP_Approx(I,J)-FP_User(I,J)'
        write ( *, * ) ' '
      else if ( iwork(1) == (-7).or.iwork(1)==(-8) ) then
        write ( *, * ) ' '
        write ( *, * ) 'DGE_SLV - Finite difference jacobian:'
        write ( *, * ) ' '
      else if ( iwork(1) == (-10) ) then
        write ( *, * ) ' '
        write ( *, * ) 'DGE_SLV - User supplied jacobian:'
        write ( *, * ) ' '
      end if

      do i = 1,nvar
        do j = 1,nvar
          k = ldf+(j-1)*nvar+i-1
          write ( *, * ) rwork(k),' i, j = ',i,j
        end do
        write ( *, * ) ' '
      end do

    end if

    return

  end if
!
!  Decompose matrix into LU factors.
!
  call dge_trf ( nvar, nvar, rwork(ldf), nvar, iwork(lpiv), info )

  iwork(20) = iwork(20) + 1

  if ( info /= 0 ) then
    if ( iwrite >= 1 ) then
      write ( *, * ) 'DGE_SLV - Fatal error!'
      write ( *, * ) '  DGE_TRF returns INFO = ', info
    end if
    ierror = 3
    return
  end if
!
!  Compute matrix determinant, and record its sign.
!
  call dge_det ( rwork(ldf), nvar, nvar, iwork(lpiv), det )

  dets = 0.0D+00

  if ( det>0.0D+00 ) then
    dets = 1.0D+00
  else if ( det<0.0D+00 ) then
    dets = -1.0D+00
  end if

  if ( job == 2)return
!
!  Solve the linear system.
!
  call dge_trs('n',nvar,1,rwork(ldf),nvar,iwork(lpiv),y,nvar,info)

  iwork(21) = iwork(21)+1
  if ( info/= 0 ) then
    if ( iwrite>= 1 ) then
      write ( *, * ) 'DGE_SLV - Fatal error!'
      write ( *, * ) '  DGE_TRS returns INFO = ',info
    end if
    ierror = 3
    return
  end if

  return
end
subroutine dge_trf ( m, n, a, lda, ipiv, info )

!*****************************************************************************80
!
!! DGE_TRF computes the PLU factorization of a general M by N matrix.
!
!  Discussion:
!
!    DGE_TRF is a standalone version of the LAPACK routine DGETRF.
!
!    The factorization uses partial pivoting with row interchanges,
!    and has the form
!      A = P * L * U
!    where P is a permutation matrix, L is lower triangular with unit
!    diagonal elements (lower trapezoidal if M > N), and U is upper
!    triangular (upper trapezoidal if M < N).
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix A.  M >= 0.
!
!    Input, integer N, the number of columns of the matrix A.  N >= 0.
!
!    Input/output, double precision A(LDA,N).
!    On entry, the M by N matrix to be factored.
!    On exit, the factors L and U from the factorization
!    A = P*L*U; the unit diagonal elements of L are not stored.
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA >= max(1,M).
!
!    Output, integer IPIV(min(M,N)), the pivot indices;
!    for 1 <= I <= min(M,N), row i of the matrix was interchanged with
!    row IPIV(I).
!
!    Output, integer INFO.
!   = 0: successful exit
!    < 0: if INFO = -K, the K-th argument had an illegal value
!    > 0: if INFO = K, U(K,K) is exactly zero. The factorization
!         has been completed, but the factor U is exactly
!         singular, and division by zero will occur if it is used
!         to solve a system of equations.
!
  integer lda
  integer n

  double precision a(lda,n)
  integer i
  integer ii
  integer info
  integer ipiv(*)
  integer j
  integer jj
  integer jp
  integer m
  double precision temp
!
!  Test the input parameters.
!
  info = 0

  if ( m < 0 ) then
    info = - 1
    return
  else if (  n < 0 ) then
    info = - 2
    return
  else if ( lda < max ( 1, m ) ) then
    info = - 4
    return
  end if

  if ( m == 0 .or. n == 0 ) then
    return
  end if

  do j = 1, min ( m, n )
!
!  Find the pivot.
!
    temp = abs ( a(j,j) )
    jp = j
    do i = j+1, m
      if ( abs ( a(i,j) ) > temp ) then
        temp = abs ( a(i,j) )
        jp = i
      end if
    end do

    ipiv(j) = jp
!
!  Apply the interchange to columns 1:N.
!  Compute elements J+1:M of the J-th column.
!
    if ( a(jp,j) /= 0.0D+00 ) then

      if ( jp /= j ) then
        do jj = 1, n
          temp = a(j,jj)
          a(j,jj) = a(jp,jj)
          a(jp,jj) = temp
        end do
      end if

      if ( j < m ) then
        do ii = j + 1, m
          a(ii,j) = a(ii,j) / a(j,j)
        end do
      end if

     else if ( info == 0 ) then

        info = j

     end if
!
!  Update the trailing submatrix.
!
     if ( j < min ( m, n ) ) then

        do ii = j+1, m
          do jj = j+1, n
            a(ii,jj) = a(ii,jj) - a(ii,j) * a(j,jj)
          end do
        end do

     end if

  end do

  return
end
subroutine dge_trs ( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

!*****************************************************************************80
!
!! DGE_TRS solves a system of linear equations factored by DGE_TRF.
!
!  Discussion:
!
!    DGE_TRS is a standalone version of the LAPACK routine DGETRS.
!
!    DGE_TRS solves a system of linear equations
!      A * X = B  or  A' * X = B
!    with a general N by N matrix A using the LU factorization computed
!    by DGE_TRF.
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*1 TRANS, pecifies the form of the system of equations:
!    'N':  A * X = B  (No transpose)
!    'T':  A'* X = B  (Transpose)
!    'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!    Input, integer N, the order of the matrix A.  N >= 0.
!
!    Input, integer NRHS, the number of right hand sides.  NRHS >= 0.
!
!    Input, double precision A(LDA,N), the factors L and U from the
!    factorization A = P*L*U as computed by DGE_TRF.
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA >= max(1,N).
!
!    Input, integer IPIV(N), the pivot indices from DGE_TRF;
!    for 1 <= i<=N, row i of the matrix was interchanged with row IPIV(I).
!
!    Input/output, double precision B(LDB,NRHS).
!    On entry, the right hand side matrix B.
!    On exit, the solution matrix X.
!
!    Input, integer LDB, the leading dimension of the array B.
!    LDB >= max(1,N).
!
!    Output, integer INFO
!   = 0:  successful exit
!    < 0:  if INFO = -I, the I-th argument had an illegal value.
!
  integer lda
  integer ldb
  integer n
  integer nrhs

  double precision a(lda,n)
  double precision b(ldb,nrhs)
  integer i
  integer info
  integer ipiv(n)
  integer j
  integer k
  double precision temp
  character trans

  info = 0

  if ( trans /= 'n' .and. trans /= 'N' .and. &
       trans /= 't' .and. trans /= 'T' .and. &
       trans /= 'c' .and. trans /= 'C' ) then
    info = - 1
    return
  else if ( n < 0 ) then
    info = - 2
    return
  else if ( nrhs < 0 ) then
    info = - 3
    return
  else if ( lda < max ( 1, n ) ) then
    info = - 5
    return
  else if ( ldb < max ( 1, n ) ) then
    info = - 8
    return
  end if

  if ( n == 0 .or. nrhs == 0 ) then
    return
  end if

  if ( trans == 'n' .or. trans == 'N' ) then
!
!  Apply row interchanges to the right hand sides.
!
    do i = 1, n
      if ( ipiv(i) /= i ) then
        do k = 1, nrhs
          temp = b(i,k)
          b(i,k) = b(ipiv(i),k)
          b(ipiv(i),k) = temp
        end do
      end if
    end do
!
!  Solve L*X = B, overwriting B with X.
!
    do k = 1, nrhs
      do j = 1, n - 1
        do i = j + 1, n
          b(i,k) = b(i,k) - a(i,j) * b(j,k)
        end do
      end do
    end do
!
!  Solve U*X = B, overwriting B with X.
!
    do k = 1, nrhs
      do j = n, 1, -1
        b(j,k) = b(j,k) / a(j,j)
        do i = 1, j - 1
          b(i,k) = b(i,k) - a(i,j) * b(j,k)
        end do
      end do
    end do

  else
!
!  Solve U'*X = B, overwriting B with X.
!
    do k = 1, nrhs
      do j = 1, n
        b(j,k) = b(j,k) / a(j,j)
        do i = j + 1, n
          b(i,k) = b(i,k) - a(j,i) * b(j,k)
        end do
      end do
    end do
!
!  Solve L'*X = B, overwriting B with X.
!
    do k = 1, nrhs
      do j = n, 2, -1
        do i = 1, j - 1
          b(i,k) = b(i,k) - a(j,i) * b(j,k)
        end do
      end do
    end do
!
!  Apply row interchanges to the solution vectors.
!
    do i = n, 1, -1
      if ( ipiv(i) /= i ) then
        do k = 1, nrhs
          temp         = b(i,k)
          b(i,k)       = b(ipiv(i),k)
          b(ipiv(i),k) = temp
        end do
      end if
    end do

  end if

  return
end
subroutine limit ( df, fpar, fx, ierror, ipar, iwork, liw, lrw, nvar, rwork, &
  slname, tc, wk, xc, xf, xr )

!*****************************************************************************80
!
!! LIMIT seeks a limit point between two continuation points.
!
!  Discussion:
!
!    The continuation points are XF and XC, and there is also a tangent
!    vector at XF in TC, and a tangent vector at XC if WK.  It is assumed
!    that the LIM-th components of these tangent vectors differ in sign,
!    indicating the presence of a limit point having a LIM-th component
!    which is exactly zero.
!
!    We solve this problem using a one-dimensional zero finder.  We
!    set up a variable SN, which measures the proportional length along
!    the secant between XF and XC.  XF can now be thought of as the
!    point on the curve corresponding to SN = 0, and XC as the point
!    corresponding to SN = 1.  For other values of SN, we compute a linear
!    combination of XF and XC, and use Newton iteration.
!
!    During this Newton iteration, we must fix a component of the solution.
!    The component is chosen as the index of the entry of largest magnitude
!    in the secant.  However, should that index be LIM, then the second
!    largest magnitude will be chosen instead.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
  external df
  external fx
  external slname

  integer, parameter :: maxit = 25
  double precision, parameter :: one = 1.0D+00

  integer liw
  integer lrw
  integer nvar

  double precision a
  double precision b
  double precision dirlpc
  double precision fa
  double precision fb
  double precision fpar(*)
  integer i
  integer icrit
  integer ierror
  integer iflag
  integer imitl
  integer ipar(*)
  integer idamax
  integer iwork(liw)
  integer iwrite
  integer lim
  integer lpc
  integer modsav
  double precision rwork(lrw)
  double precision skale
  double precision sn
  double precision snl
  double precision tc(nvar)
  double precision temp
  double precision tsn
  double precision wk(nvar)
  double precision xabs
  double precision xc(nvar)
  double precision xdif
  double precision xf(nvar)
  double precision xr(nvar)

  iwrite = iwork(7)
  lim = iwork(6)
  if ( 2 <= iwrite ) then
    write ( *, * )  'LIMIT  - Attempt correction of approximate limit point.'
  end if

  modsav = iwork(4)
!
!  The limit point occurs somewhere between the current and previous
!  points.  See if either already computed point can be accepted as
!  the limit point.  If the current point XC is the limit point, then
!  its tangent is already in WK.
!
  if ( abs ( rwork(27) ) <= 0.5D+00 * rwork(1) ) then

    xr(1:nvar) = xc(1:nvar)

    go to 30

  else if ( abs ( rwork(26) ) <= 0.5D+00 * rwork(1) ) then

    xr(1:nvar) = xf(1:nvar)
    wk(1:nvar) = tc(1:nvar)

    go to 30
  end if
!
!  If interval is extremely small, simply assign one
!  endpoint of the interval as the answer.
!
  xdif = abs ( xc(lim) - xf(lim) )
  xabs = max ( abs ( xc(lim) ), abs ( xf(lim) ) )

  if ( xdif <= 8.0D+00 * rwork(8) * ( 1.0D+00 + xabs ) ) then

    if ( abs ( rwork(27) ) > abs ( rwork(26) ) ) then

      xr(1:nvar) = xf(1:nvar)
      wk(1:nvar) = tc(1:nvar)

    else

      xr(1:nvar) = xc(1:nvar)

    end if

    go to 30

  end if
!
!  Begin root-finding iteration on interval (0,1), with function values
!  TLLIM and TCLIM.
!
  a = 0.0D+00
  fa = rwork(27)
  b = 1.0D+00
  fb = rwork(26)
!
!  Find LPC, the index of the entry of maximum absolute value in the
!  secant between the two continuation points.  However, we will not
!  allow LPC to equal LIM.  We save the sign of the LPC-th entry of
!  the secant so that new tangents may be properly signed.
!
!  Note: If the value of IWORK(3) is 1, then the user has requested that
!  the parameterization index always be held fixed at the value in IWORK(2).
!  The user's choice overrides all considerations, even in this routine.
!  Certain doom would occur if the user chooses unwisely, but that is
!  not our concern!
!
  if ( iwork(3) /= 1 ) then
    temp = xc(lim)
    xc(lim) = xf(lim)
    xr(1:nvar) = xf(1:nvar) - xc(1:nvar)
    lpc = idamax ( nvar, xr, 1 )
    xc(lim) = temp
  else
    lpc = iwork(2)
  end if

  dirlpc = sign ( one, xf(lpc) - xc(lpc) )
!
!  The first approximation to the limit point will be whichever endpoint
!  has the smallest LIM-th component of the tangent vector.
!
  if ( abs ( rwork(26) ) >= abs ( rwork(27) ) ) then
    sn = 0.0D+00
    tsn = rwork(27)
    xr(1:nvar) = xc(1:nvar)
  else
    sn = 1.0D+00
    tsn = rwork(26)
    xr(1:nvar) = xf(1:nvar)
  end if

  imitl = 0
  if ( iwrite >= 2 ) then
    write ( *, * ) 'LIMIT - For S = ',0.0D+00,' Tan(X(S))(Lim)=',rwork(27)
    write ( *, * ) 'LIMIT - For S = ',1.0D+00,' Tan(X(S))(Lim)=',rwork(26)
  end if
!
!  Call the rootfinder repeatedly for the approximate root SN.
!  Use a linear combination of the points X(SNL) and X(0.0) or X(1.0)
!  to get a starting point for correction to the curve, returning
!  to the curve along the line X(LPC) = constant.  Compute the
!  tangent vector there, and return the LIM-th component of the
!  tangent as the function value whose zero we seek.
!
10    continue

  snl = sn
  call root_finder ( a, fa, b, fb, sn, tsn, imitl, iflag, ierror )
  iwork(23) = iwork(23)+1

  if ( ierror/= 0 ) then
    if ( iwrite>= 1 ) then
      write ( *, * ) 'LIMIT - Fatal error!'
      write ( *, * ) '  ROOT_FINDER returns IERROR = ',ierror
    end if
    ierror = 7
    return
  end if

  if ( iflag == (-1) .or. iflag == 0 ) go to 30
!
!  Find whether SN lies in (0.0,SNL) or (SNL,1.0).  This will determine
!  how we construct our linear combination of the current points to
!  get a starting point for X(SN).  This somewhat cumbersome procedure
!  occurs because we limit the number of vectors we use.
!
  if ( sn <= snl ) then
!
!  If SN lies between 0.0 and SNL, then set the approximation to the point
!  on the curve parameterized by the value SN as:
!
!  X(SN) = (SNL-SN)/(SNL-0.0) * X(0.0) + (SN-0.0)/(SNL-0.0) * X(SNL)
!
    if ( snl <= 0.0D+00 ) then
      skale = 0.0D+00
    else
      skale = sn / snl
      skale = max ( skale, 0.0D+00 )
      skale = min ( skale, 1.0D+00 )
    end if

    xr(1:nvar) = skale * xr(1:nvar) + ( 1.0D+00 - skale ) * xc(1:nvar)
!
!  Otherwise, if SN lies between SNL and 1.0, set
!
!  X(SN) = (SN-SNL)/(1.0-SNL)*X(1.0)+(1.0-SN)/(1.0-SNL)*X(SNL)
!
  else
    if ( snl>= 1.0D+00 ) then
      skale = 0.0D+00
    else
      skale = (1.0D+00-sn) / (1.0D+00-snl)
      skale = max ( skale, 0.0D+00 )
      skale = min ( skale, 1.0D+00 )
    end if

    xr(1:nvar) = skale * xr(1:nvar) + ( 1.0D+00 - skale ) * xf(1:nvar)

  end if
!
!  Try to correct the approximate point so that it lies on the curve.
!  If the user is trying to economize, by using a nonzero value of
!  IWORK(4), then we may try to recover from a failed correction by
!  retrying it with a smaller value of IWORK(4).
!
20    continue

  icrit = 0
  call corrector ( df,fpar,fx,ierror,lpc,ipar,iwork,nvar,rwork,wk,xr, &
    lrw,liw,icrit,slname)
!
!  If the correction fails, see if we can retry.
!
  if ( ierror/= 0 .and. iwork(4)>0 ) then
    ierror = 0
    iwork(4) = iwork(4)-1
    xr(1:nvar) = ( 1.0D+00 - sn ) * xc(1:nvar) + sn * xf(1:nvar)

    if ( iwrite>= 1) then
      write ( *, * ) 'LIMIT - Retry limit computation with IWORK(4) = ',iwork(4)
    end if

    go to 20
  end if

  iwork(4) = modsav
  if ( ierror /= 0 ) then
    if ( iwrite >= 1 ) then
      write ( *, * ) 'LIMIT - Fatal error!'
      write ( *, * ) '  Corrector returned IERROR = ',ierror
    end if
    ierror = 7
    return
  end if
!
!  Compute the tangent at the new point.
!
  call tangnt ( temp, fx, df, fpar, ierror, lpc, ipar, iwork, nvar, &
    rwork, wk, xr, liw, lrw, slname )

  if ( ierror /= 0 ) then
    if ( iwrite >= 1 ) then
      write ( *, * ) 'LIMIT - Fatal error!'
      write ( *, * ) '  TANGNT returned IERROR = ',ierror
    end if
    ierror = 7
    return
  end if
!
!  Adjust the sign of the tangent vector so the LPC-th component
!  has the same sign as the LPC-th component of the secant.
!
  if ( dirlpc /= sign(one,wk(lpc)) ) then
    wk(1:nvar) = - wk(1:nvar)
  end if
!
!  See if we can accept the new point as the actual limit point, because
!  the LIM-th component of the tangent is acceptably small.
!
  tsn = wk(lim)
  if ( iwrite >= 2 ) then
    write ( *, * ) 'LIMIT - For S = ',sn,' Tan(X(S))(Lim)=',tsn
  end if

  if ( abs(tsn) <= rwork(1) ) go to 30
!
!  See if we have reached MAXIT iterations
!
  if ( imitl < maxit ) go to 10
!
!  The limit point iteration has not produced a point which satisfies
!  our requirement.  We set an error flag, but we return the partial results
!  of the abortive computation anyway.
!
  ierror = 8
  if ( iwrite >= 1 ) then
    write ( *, * ) 'LIMIT - Warning!'
    write ( *, * ) '  Iteration did not reach limit point after '
    write ( *, * ) '  taking ',maxit,' steps.'
  end if
!
!  The limit point iteration is over.
!  Compute and store information.
!
30    continue

  temp = 0.0D+00
  do i = 1, nvar
    wk(i) = xr(i) - xc(i)
    temp = temp + wk(i)**2
  end do
  temp = sqrt ( temp )

  rwork(14) = rwork(12) + temp
  iwork(27) = iwork(27) + 1
  iwork(1) = 4

  return
end
subroutine root_finder ( a, fa, b, fb, u, fu, kount, iflag, ierror )

!*****************************************************************************80
!
!! ROOT_FINDER seeks a root of the scalar equation F(X) = 0.0.
!
!  Discussion:
!
!    ROOT must be called repeatedly to find the root.  On the first call,
!    ROOT is given a starting interval (A,B), on which F changes sign,
!    and the function values FA and FB.  ROOT returns a new point U
!    at which the value of F is to be computed.
!
!    The user may accept U as the root, or more likely, return FU
!    to ROOT, allowing it to make a new and better guess for the root.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization without Derivatives,
!    Prentice Hall, 1973.
!
!  Parameters:
!
!    Input/output, real A.  A is one endpoint of the current
!    interval in which the root is sought.  The user must set this
!    before the first call only.  Thereafter, ROOT adjusts A as
!    the interval shrinks.
!
!    Input/output, real FA.  On the very first call only, the
!    user must set FA to the value of F(A).  Thereafter, the
!    program resets FA as the value of A changes, and the user
!    should not alter the value of FA.
!
!    Input/output, real B.  B is one endpoint of the current
!    interval in which the root is sought.  The user must set
!    this before the first call only.  Thereafter, ROOT adjusts
!    B as the interval shrinks.
!
!    Input/output, real FB.  On the very first call only, the
!    user must set FB to the value of F(B).  Thereafter, the
!    program resets FB as the value of B changes, and the user
!    should not alter the value of FB.
!
!    Output, real U, the current approximation to the root.  After every
!    call to ROOT, U will contain the routine's best approximation to
!    the location of the root.
!
!    Input, real FU.  On the first call, the user should not set FU.
!    Before the second call, the user should evaluate the function at
!    the point U, returned on the first call, and send that value as FU.
!    This behavior should then be repeated on each subsequent call.
!    The output value of U should be used to evaluate the function, with
!    the result sent back as FU on the next call.
!
!    Input/output, integer KOUNT.  On the first call, set KOUNT to zero.
!    Thereafter, the program will update KOUNT, which counts the number
!    of calls made to ROOT.
!
!    Output, integer IFLAG, reports the status of the search.
!
!    -1, the current bracketing interval (A,B) or (B,A) is smaller than
!    4*EPMACH*ABS(U)+EPMACH, and so U should be accepted as the root.
!
!    0, the input value FU is exactly zero, so U should be accepted
!    as the root.
!
!    Positive values of IFLAG indicated that the program has found a
!    new approximation U to the root.  If a better approximation is
!    desired, return the value of the function in FU, and the program
!    will proceed with the search.  The actual value of IFLAG tells
!    what method was used to produce the current U:
!
!    1, bisection
!    2, linear interpolation
!    3, inverse quadratic interpolation.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    7, FA*FB is greater than 0, and so the given interval is unacceptable.
!
  double precision, parameter :: one = 1.0D+00

  double precision a
  double precision b
  double precision epmach
  double precision fa
  double precision fb
  double precision fu
  double precision halfub
  integer ierror
  integer iflag
  integer kount
  double precision p
  double precision q
  double precision r
  double precision s
  double precision sdel1
  double precision sdel2
  double precision sdel3
  double precision sdel4
  double precision step
  double precision toler
  double precision u

  epmach = epsilon ( epmach )
!
!  Segment 1.  The first call is handled specially.
!
  if ( kount <= 0 ) then
    if ( ( fa > 0.0D+00 .and. fb > 0.0D+00 ) .or. &
         ( fa < 0.0D+00 .and. fb < 0.0D+00 ) ) then
      ierror = 7
      kount = 0
      return
    end if
    kount = 1
    sdel1 = 2.0D+00 * abs ( b - a )
    sdel2 = 2.0D+00 * sdel1
    sdel3 = 2.0D+00 * sdel2
    u = b
    b = a
    fu = fb
    fb = fa
  else
!
!  On calls after the first call, increment the counter, and check
!  whether F(U) is zero.
!
    kount = kount + 1
    if ( fu == 0.0D+00 ) then
      iflag = 0
      return
    end if
!
!  If FU has the same sign as FB, then store the value of A in B.
!
    if ( sign ( one, fu ) == sign ( one, fb ) ) then
      b = a
      fb = fa
    end if
  end if
!
!  Segment 2.  Rearrange points if necessary to ensure that ABS(FU)<ABS(FB).
!
  if ( abs ( fb ) < abs ( fu ) ) then
    a = u
    u = b
    b = a
    fa = fu
    fu = fb
    fb = fa
  end if
!
!  Segment 3.  Check to see if we can accept the current estimate because
!  the current change-in-sign interval (B,U) or (U,B) is very small.
!
  toler = 2.0D+00 * epmach * abs ( u ) + epmach
  halfub = ( b - u ) / 2.0D+00
  sdel4 = sdel3
  sdel3 = sdel2
  sdel2 = sdel1
  sdel1 = abs ( b - u )

  if ( abs ( halfub ) <= toler ) then
    iflag = -1
    a = u
    fa = fu
    return
  end if
!
!  Segment 4.  Compute a new approximate root, of the form U(new) = U(old)+STEP.
!  Methods availabe are linear interpolation, inverse quadratic interpolation,
!  and bisection.
!
  if ( abs ( fu ) >= abs ( fa ) ) then
    iflag = 1
    step = halfub
    go to 10
  end if
!
!  If only two points are available, use linear interpolation.
!
  if ( a == b ) then
    iflag = 2
    s = fu / fa
    p = 2.0D+00 * halfub * s
    q = 1.0D+00 - s
!
!  If three points are available, try inverse quadratic
!  interpolation.
!
  else
    iflag = 3
    s = fu / fa
    q = fa / fb
    r = fu / fb
    p = s * ( 2.0D+00 * halfub * q * ( q - r ) - ( u - a ) * ( r - 1.0D+00 ) )
    q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )
  end if
!
!  Correct the signs of P and Q.
!
  if ( p > 0.0D+00 ) then
    q = -q
  else
    p = -p
  end if
!
!  If P/Q is too large, use bisection instead.
!
  if ( (8.0D+00*sdel1>sdel4).or.(p>= 1.5D+00*abs(halfub*q)-abs(toler*q)) ) then
    iflag = 1
    step = halfub
    go to 10
  end if
  step = p / q
!
!  Segment 5.  The value of STEP is known.  Update information.
!  The change in sign intevarl is now (A,B) or (B,A).
!
10    continue

  a = u
  fa = fu
  if ( abs ( step ) <= toler ) then
    step = sign ( toler, halfub )
  end if
  u = u + step

  return
end
subroutine setstp ( iwork, liw, lrw, rwork )

!*****************************************************************************80
!
!! SETSTP computes the stepsize to be used by the Euler prediction step.
!
!  Discussion:
!
!  The formulas underlying the algorithm are:
!
!    ALFMIN = A minimal angle, ALFMIN=2*ARCCOS(1-EPMACH).
!
!    ALPHLC = Angle between last two tangents, value of ARCCOS(TL dot TC),
!             except that ALPHLC must be at least equal to ALFMIN.
!
!    HSEC  = Euclidean norm of the secant step, NORM2(XC-XF).
!
!    HSECL = Euclidean norm of previous secant step, NORM2(XL-XC).
!
!    ABSNLC = ABS(SIN(.5*ALPHLC))
!
!    CURV  = Previous value of the curvature.
!
!    CURVN = 2*ABSNLC/HSEC
!
!    CORDIS = Distance between predicted and corrected points.
!             Adjust CORDIS to lie between 0.01*HSEC and HSEC,
!
!             But if CORDIS = 0, meaning the predicted point was accepted,
!             set HTAN = HFACT*HSEC instead of using the first estimate
!             for HTAN.
!
!  Then
!
!    CURVX = CURVN+HSEC*(CURVN-CURV)/(HSEC+HSECL)
!
!  A simpler formula is used if we do not have data at two old points.
!
!  If (IWORK(10).GE.2) CURVX must be at least as large as the maximum
!  of 0.001 and 0.01/HSEC.
!
!  First estimate for stepsize (unless CORDIS = 0):
!
!    HTAN = SQRT(2*QUAL*CORDIS/CURVX)
!
!  Adjusted value:
!
!    HTAN = HTAN*(1.0+HTAN*(TC(IPC)-TL(IPC))/(2*HSEC*TC(IPC)))
!
!  Readjustments to the calculated value:
!
!  If stepsize reduction occurred during the correction of the previous
!    continuation point, HTAN is forced to be less than (HFACT-1)*HSEC/2.
!  The calculated HTAN is forced to lie between (HSEC/HFACT) and (HSEC*HFACT).
!  The calculated HTAN is also forced to lie between HMIN and HMAX.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
  integer liw
  integer lrw

  double precision curv
  double precision curvx
  double precision htan
  integer iwork(liw)
  double precision rwork(lrw)
  double precision temp

  rwork(11) = max ( rwork(11), rwork(10) )
!
!  Update estimates of curvature.
!
  curv = rwork(16)
  rwork(16) = 2.0D+00 * abs ( sin ( 0.5D+00 * rwork(11) ) ) / rwork(21)

  if ( curv == 0.0D+00) then
    curv = rwork(16)
  end if

  if ( rwork(22) == 0.0D+00 ) then
    curvx = rwork(16)
  else
    curvx = rwork(16) + rwork(21)*(rwork(16)-curv)/(rwork(21)+rwork(22))
  end if

  curvx = max ( curvx, 0.001D+00 )
  curvx = max ( curvx, 0.01D+00 / rwork(21) )
!
!  If the convergence distance was zero, then set the tentative next
!  stepsize to the maximum growth factor times the size of the secant step.
!
!  Otherwise, use the curvature estimate and other information to
!  compute an optimal step.
!
  if ( rwork(15) == 0.0D+00 ) then
    htan = rwork(20)*rwork(21)
  else
    temp = rwork(23) * rwork(15)
    temp = max ( temp, rwork(21) / 100.0D+00 )
    temp = min ( temp, rwork(21) )
    htan = sqrt ( 2.0D+00 * temp / curvx )
  end if
!
!  Adjust the step to account for estimated curvature in the direction
!  of the parameter.
!
  if ( iwork(18)>0) then
    htan = min(htan,(rwork(20)-1.0D+00)*rwork(21)/2.0D+00)
  end if

  if ( iwork(3)/= 1 ) then
    temp = 1.0D+00+(1.0D+00-rwork(25)/rwork(24))*(htan/2.0D+00)/rwork(21)
    htan = htan*temp
  end if
!
!  Enforce growth rate restrictions.
!
  htan = max ( htan, rwork(21) / rwork(20) )
  htan = min ( htan, rwork(21) * rwork(20) )
!
!  Enforce size restrictions.
!
  htan = max ( htan, rwork(3) )
  htan = min ( htan, rwork(4) )

  rwork(5) = htan
  if ( iwork(7)>= 2 ) then
    write ( *, * ) 'SETSTP - Next stepsize HTAN = ',htan
  end if

  return
end
subroutine stajac(df,fpar,fx,ierror,ipar,ipc,iwork,liw,lrw,nvar,rwork,wk,xr, &
  slname)

!*****************************************************************************80
!
!! STAJAC generates and factors the jacobian matrix.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
  external df
  external fx
  external slname

  integer liw
  integer lrw
  integer nvar

  double precision dets
  double precision fpar(*)
  integer ierror
  integer ipar(*)
  integer ipc
  integer iwork(liw)
  integer iwrite
  integer job
  double precision rwork(lrw)
  double precision wk(nvar)
  double precision xr(nvar)

  iwrite = iwork(7)
!
!  If user is requesting that Jacobian be used as long as possible,
!  then go ahead, generate and factor the first one now.
!
  if ( iwork(4) == 2 ) then

    if ( iwrite>= 2 ) then
      write ( *, * ) 'STAJAC - Generating initial jacobian.'
    end if

    job = 2
    call slname(dets,fx,df,fpar,ierror,ipc,ipar,iwork,liw,job,nvar,rwork, &
      lrw,xr,wk)
    rwork(17) = dets

    if ( ierror/= 0 ) then

      if ( iwrite>= 1 ) then
        write ( *, * ) 'STAJAC - Serious error!'
        write ( *, * ) '  Could not factor initial jacobian.'
      end if

      return

    end if

  end if

  return
end
subroutine start(df,fpar,fx,ierror,ipar,ipc,iwork,liw,lrw,nvar,rwork,wk, &
  xc,xf,xr,slname)

!*****************************************************************************80
!
!! START forces the starting point to satisfy the nonlinear system.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
  external df
  external fx
  external slname

  integer liw
  integer lrw
  integer nvar
  double precision fpar(*)
  integer i
  integer icrit
  integer ierror
  integer imax
  integer ipar(*)
  integer ipc
  integer idamax
  integer iwork(liw)
  integer iwrite
  integer modnew
  integer modsav
  double precision qual
  double precision rwork(lrw)
  double precision wk(nvar)
  double precision xc(nvar)
  double precision xf(nvar)
  double precision xr(nvar)

  iwrite = iwork(7)

  if ( iwrite>= 2 ) then
    write ( *, * ) 'START  - Checking the initial point.'
    write ( *, * ) '  Fixing variable number',ipc
  end if

  xc(1:nvar) = xr(1:nvar)
  modsav = iwork(4)
  icrit = 1

20    continue

  xr(1:nvar) = xc(1:nvar)

  call corrector ( df,fpar,fx,ierror,ipc,ipar,iwork,nvar,rwork,wk,xr, &
    lrw,liw,icrit,slname)

  iwork(25) = iwork(25)+iwork(28)
!
!  If an error occurred, then, if possible, retry with ICRIT = 2.
!
  if ( ierror /= 0 ) then

    if ( icrit == 1 ) then

      if ( iwrite >= 1 ) then
        write ( *, * ) 'START - Warning!'
        write ( *, * ) '  The starting point needs to be corrected.'
        write ( *, * ) '  The first correction attempt failed.'
        write ( *, * ) '  Correction will be retried once.'
      end if

      icrit = 2
      go to 20

    else

      icrit = 1

    end if

    if ( iwork(4)>0 ) then

      iwork(4) = iwork(4)-1
      ierror = 0

      if ( iwrite>= 1 ) then
        write ( *, * ) 'START  - Retry starting point correction.'
        write ( *, * ) '  Set Newton option IWORK(4) = ',iwork(4)
      end if

      go to 20

    end if

  end if

  iwork(4) = modsav

  if ( ierror/= 0 ) then

    if ( iwrite>= 1 ) then
      write ( *, * ) 'START - Fatal error!'
      write ( *, * ) '  Correction of the starting point failed.'
    end if

    return

  end if
!
!  If we were able to correct the starting point, then
!  record necessary data.
!
  xc(1:nvar) = xc(1:nvar) - xr(1:nvar)
  imax = idamax ( nvar, xc, 1 )
  rwork(15) = abs ( xc(imax) )

  xc(1:nvar) = xr(1:nvar)
  xf(1:nvar) = xr(1:nvar)

  modnew = iwork(4)
  call coqual ( modnew, qual, iwork, liw, rwork, lrw )

  rwork(23) = qual
  rwork(14) = rwork(13)
  iwork(27) = iwork(27)+1
  iwork(10) = 1
  iwork(1) = 1

  return
end
subroutine tangnt ( detsn, fx, df, fpar, ierror, ip, ipar, iwork, nvar, &
  rwork, tan, xr, liw, lrw, slname )

!*****************************************************************************80
!
!! TANGNT computes a tangent vector to the solution curve.
!
!  Discussion:
!
!    The length of the tangent vector in the Euclidean norm will be 1.
!    There are two such tangent vectors at each point, one the negative
!    of the other.  This routine produces one such vector.  The choice of
!    which is appropriate must be made outside this routine.
!
!    The tangent vector TAN is the solution of the linear system
!
!      DFA(X,IP)*TAN = E(NVAR)
!
!    Here E(I) denotes the I-th basis vector, that is, the vector with
!    a 1 in the I-th position and 0 elsewhere.  DFA(X,IP) is the NVAR by
!    NVAR matrix whose first NVAR-1 rows are the jacobian of FX evaluated
!    at X, and whose last row is the transpose of E(IP).
!
!    After computation, the tangent vector is normalized so that it has
!    a Euclidean norm of 1.  The adjustment of sign is performed outside
!    this routine.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, data or storage error.
!    2, error occurred in user jacobian routine.
!    3, an error occurred in the solver.
!    6, the computed tangent vector is entirely zero.
!
  external fx
  external df
  external slname

  integer liw
  integer lrw
  integer nvar

  double precision detsn
  double precision fpar(*)
  integer i
  integer ierror
  integer ip
  integer ipar(*)
  integer iwork(liw)
  integer iwrite
  integer job
  double precision rwork(lrw)
  double precision tan(nvar)
  double precision tnorm
  double precision xr(nvar)

  iwrite = iwork(7)
!
!  Set the right hand side of the linear system.
!
  tan(1:nvar-1) = 0.0D+00
  tan(nvar) = 1.0D+00
!
!  Call the user-specified solver.
!
  if ( iwork(4) == 2) then
    job = 1
  else
    job = 0
  end if

  call slname(detsn,fx,df,fpar,ierror,ip,ipar,iwork,liw,job,nvar,rwork, &
    lrw,xr,tan)

  if ( ierror/= 0 ) then

    if ( iwrite>= 1 ) then
      write ( *, * ) 'TANGNT - Warning!'
      write ( *, * ) 'The linear solver returned error flag IERROR = ',ierror
    end if

    return

  end if
!
!  Normalize the tangent vector.
!
  tnorm = 0.0D+00
  do i = 1, nvar
    tnorm = tnorm + tan(i)**2
  end do
  tnorm = sqrt ( tnorm )

  if ( tnorm == 0.0D+00 ) then

    ierror = 6

    if ( iwrite >= 1 ) then
      write ( *, * ) 'TANGNT - Warning!'
      write ( *, * ) 'The computed tangent has zero norm!'
      return
    end if

  end if

  tan(1:nvar) = tan(1:nvar) / tnorm

  return
end
subroutine tanpar ( df, fpar, fx, ierror, ipar, iwork, liw, lrw, nvar, &
  rwork, slname, tc, wk, xc, xf, xr )

!*****************************************************************************80
!
!! TANPAR computes the tangent vector and the next continuation index.
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external DF, the name of the Jacobian evaluation routine.
!
!    Input/output, real FPAR(*), a user defined parameter array.
!
  double precision, parameter :: one = 1.0D+00

  external df
  external fx
  external slname

  integer liw
  integer lrw
  integer nvar

  double precision atcipc
  double precision atcjpc
  double precision fpar(*)
  integer i
  integer idamax
  integer ierror
  integer ipar(*)
  integer ipc
  integer ipl
  integer jpc
  integer itemp
  integer iwork(liw)
  integer iwrite
  double precision rwork(lrw)
  double precision tc(nvar)
  double precision tcipc
  double precision tcos
  double precision temp
  double precision tlipc
  double precision tsin
  double precision wk(nvar)
  double precision xc(nvar)
  double precision xf(nvar)
  double precision xr(nvar)

  iwrite = iwork(7)
!
!  Save data that's about to become outdated:
!
!    IPL       <= the value of the continuation parameter.
!    WK        <= the tangent.
!    RWORK(29) <= the sign of the determinant.
!
  ipl = iwork(2)

  wk(1:nvar) = tc(1:nvar)

  rwork(29) = rwork(17)
!
!  Compute the tangent for the current point and store in TC.
!
  call tangnt ( rwork(17),fx,df,fpar,ierror,iwork(2),ipar,iwork,nvar, &
    rwork,tc,xf,liw,lrw,slname)

  if ( ierror/= 0 ) then

    if ( iwrite>= 1 ) then
      write ( *, * ) 'TANPAR - Serious error!'
      write ( *, * ) '  The tangent calculation failed.'
    end if

    return

  end if
!
!  Find the two largest components of tangent vector, which are our
!  candidates for continuation directions, unless the user demands
!  a fixed parameterization.
!
  if ( iwork(3)/= 1 ) then
    ipc = idamax(nvar,tc,1)
    temp = tc(ipc)
    tc(ipc) = 0.0D+00
    jpc = idamax(nvar,tc,1)
    if ( jpc == 0) then
      jpc = ipc
    end if
    tc(ipc) = temp
  else
    ipc = iwork(2)
    jpc = ipc
  end if
!
!  Adjust the sign of the tangent vector TC.  To do so, compare the sign
!  of the IPL-th component with the sign of the IPL-th component of
!  the secant vector XF-XC.  If a target or limit point, XR, has been
!  computed, use XF-XR instead of XF-XC.  And on the first step, we'll
!  have to use the user's input direction for our comparison.
!
  if ( iwork(10) <= 1 ) then
    temp = rwork(6)
  else if ( iwork(1) == 3.or.iwork(1)==4 ) then
    temp = xf(ipl)-xr(ipl)
  else
    temp = xf(ipl)-xc(ipl)
  end if

  if ( sign(one,tc(ipl))/= sign(one,temp) ) then

    tc(1:nvar) = - tc(1:nvar)
    rwork(17) = -rwork(17)

  end if
!
!  Unless we are computing a starting point, record the new state.
!
  if ( iwork(10)>1) then
    iwork(10) = 3
  end if
!
!  The index of the entry of largest magnitude in the tangent vector will
!  be used for the local parameterization, unless the user has ordered
!  us to stick to a given parameter, or if a limit point in the first
!  choice appears to be approaching.
!
!  To check this, we compare TC(IPC) and the second largest component, TC(JPC).
!  If TC(JPC) is no less than 0.1 of TC(IPC), and
!     TC(JPC) is larger than TL(JPC), and
!     TC(IPC) is smaller than TL(IPC), then
!  we suspect a limit point may be coming in the IPC direction and we switch
!  our choice to JPC.
!
  if ( iwork(3)/= 1 ) then

    atcipc = abs(tc(ipc))
    atcjpc = abs(tc(jpc))

    if ( jpc/= ipc.and.iwork(10)>1 ) then

      tlipc = wk(ipc)
      tcipc = tc(ipc)
      temp = abs(wk(jpc))

      if ( (sign(one,tcipc) == sign(one,tlipc)) &
           .and.(atcipc<abs(tlipc)) &
           .and.(atcjpc>= max(0.1*atcipc,temp)) ) then

        if ( iwrite>= 3 ) then
          write ( *, * ) 'TANPAR - A limit point may be coming'
          write ( *, * ) 'in the preferable index ',ipc
          write ( *, * ) 'We''ll try second-best index.'
        end if

        itemp = ipc
        ipc = jpc
        jpc = itemp

      end if

    end if

    if ( iwrite>= 3 ) then
      write ( *, * ) 'TANPAR - Continuation index: First choice = ',ipc
      if ( jpc/= ipc) then
        write ( *, * ) '                             Second choice = ',jpc
      end if
    end if

  end if
!
!  Record the values of the IPC-th component of the new tangent vector TC,
!  as well as the IPC-th component of the old tangent vector TL.
!  Set the sign of the determinant.
!  Record the value of the LIM-th compoment of the new tangent vector,
!  if a limit point check is being done.
!
  iwork(2) = ipc
  iwork(12) = jpc
  rwork(24) = tc(ipc)
  rwork(25) = wk(ipc)
  rwork(6) = rwork(17)

  if ( iwork(6)>0 ) then
    rwork(27) = rwork(26)
    rwork(26) = tc(iwork(6))
    if ( iwrite>= 2) then
      write ( *, * ) 'TANPAR - Tangent vector has limit component = ',rwork(26)
    end if
  end if
!
!  Compute the angle between the old and new tangents.
!
  if ( iwork(10) <= 1 ) then
    return
  end if

  tcos = dot_product ( wk(1:nvar), tc(1:nvar) )
  tcos = min ( tcos,   1.0D+00 )
  tcos = max ( tcos, - 1.0D+00 )

  tsin = sqrt ( 1.0D+00 - tcos**2 )
  rwork(11) = atan2 ( tsin, tcos )
!
!  Note possible bifurcation point, if determinant has changed sign.
!
!     if ( rwork(17)/= rwork(29) ) then
!       if ( iwrite>= 1)write ( *, * ) 'TANPAR - Possible bifurcation point.'
!     end if

  return
end
subroutine target(df,fpar,fx,ierror,ifound,ipar,iwork,liw,lrw,nvar,rwork, &
  slname,wk,xc,xf,xr)

!*****************************************************************************80
!
!! TARGET controls the computation of a target point.
!
!  Discussion:
!
!    The calling code has detected that between the points XC and XF must lie a
!    point XR whose IT-th component has the value XIT desired.  TARGET
!    must compute that point.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IFOUND.  Reports whether a target point was detected.
!     0, no target point was detected.
!    +1, A target point was detected.
!
  external df
  external fx
  external slname

  integer liw
  integer lrw
  integer nvar

  double precision fpar(*)
  integer i
  integer icrit
  integer ierror
  integer ifound
  integer ipar(*)
  integer it
  integer iwork(liw)
  integer iwrite
  integer modsav
  double precision rwork(lrw)
  double precision skale
  double precision temp
  double precision wk(nvar)
  double precision xlow
  double precision xup
  double precision xc(nvar)
  double precision xf(nvar)
  double precision xr(nvar)

  it = iwork(5)
  if ( it < 1 .or. it > nvar ) then
    return
  end if

  iwrite = iwork(7)
  ifound = 0

  if ( iwork(10) <= 1 ) then
    return
  end if

  if (  iwork(1) == 3 .and. rwork(7) == rwork(28) .and. it == iwork(11) ) then
    return
  end if

  xlow = xc(it)
  xup = xf(it)

  if ( ( rwork(7) < xlow .and. rwork(7) < xup ) .or. &
       ( rwork(7) > xlow .and. rwork(7) > xup ) ) then
    return
  end if

  ifound = 1

  if ( iwrite >= 2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'TARGET - A target point has been detected.'
  end if
!
!  Approximate the target point using the bracketing solutions.
!
  modsav = iwork(4)

10    continue

  if ( xlow /= xup ) then
    skale = ( rwork(7) - xlow ) / ( xup - xlow )
  else
    skale = 1.0D+00
  end if

  xr(1:nvar) = skale * xf(1:nvar) + ( 1.0D+00 - skale ) * xc(1:nvar)

  xr(it) = rwork(7)
!
!  Call CORRECTOR to compute the exact target point, holding index IT fixed.
!
  icrit = 0

  call corrector ( df,fpar,fx,ierror,it,ipar,iwork,nvar,rwork,wk,xr, &
    lrw,liw,icrit,slname)

  iwork(24) = iwork(24) + iwork(28)

  if ( ierror /= 0 .and. iwork(4) > 0 ) then
    ierror = 0
    iwork(4) = iwork(4)-1
    if ( iwrite >= 1 ) then
      write ( *, * ) 'TARGET - Retry computation with IWORK(4) = ',iwork(4)
    end if
    go to 10
  end if

  iwork(4) = modsav

  if ( ierror /= 0 ) then
    if ( iwrite >= 1 ) then
      write ( *, * ) 'TARGET - Target point calculation failed.'
    end if
    return
  end if
!
!  Record the values of IT and XIT, and compute the arclength to the target
!  point.
!
  iwork(1) = 3
  iwork(11) = it
  rwork(28) = rwork(7)

  temp = 0.0D+00
  do i = 1, nvar
    wk(i) = xr(i) - xc(i)
    temp = temp + wk(i)**2
  end do
  temp = sqrt ( temp )

  rwork(14) = rwork(12) + temp
  iwork(27) = iwork(27)+1

  return
end
subroutine trystp ( df, fpar, fx, ierror, ipar, iwork, liw, lrw, nvar, &
  rwork, slname, tc, wk, xf, xr )

!*****************************************************************************80
!
!! TRYSTP tries to carry out a continuation step.
!
!  Discussion:
!
!    Given a point X, a tangent vector TAN, and a stepsize H, TRYSTP
!    estimates the value of the next point on the curve as X + H*TAN,
!    and then uses Newton iteration to correct this estimate.
!
!    TRYSTP tries various fixes if the Newton iteration fails to converge,
!    or does not converge rapidly enough.
!
!  Modified:
!
!    30 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  integer liw
  integer lrw
  integer nvar

  double precision fpar(*)
  integer i
  integer icrit
  integer ierror
  integer ihold
  integer ipar(*)
  integer ipc
  integer iwork(liw)
  integer iwrite
  integer jpc
  integer modsav
  double precision rwork(lrw)
  double precision tc(nvar)
  double precision temp
  double precision wk(nvar)
  double precision xf(nvar)
  double precision xr(nvar)

  external df
  external fx
  external slname

  iwork(18) = 0
  ihold = iwork(2)
  ipc = iwork(2)
  iwrite = iwork(7)
  jpc = iwork(12)
  modsav = iwork(4)
!
!  Set initial guess for next point to XR = XF + RWORK(5) * TC
!
10    continue

  if ( iwrite>= 2 ) then
    write ( *, * ) 'TRYSTP - Predictor using stepsize ',rwork(5)
  end if

  xr(1:nvar) = xf(1:nvar) + rwork(5) * tc(1:nvar)
!
!  If you return to this statement from a later statement, we will
!  NOT be changing the current iterate, but rather, altering the
!  parameter index, trying to "save" the iteration.
!
20    continue

  if ( iwrite >= 2 .and. iwork(3) == 0 ) then
    write ( *, * ) 'TRYSTP - Corrector is fixing index ',ihold
  end if

  icrit = 0

  call corrector ( df, fpar, fx, ierror, ihold, ipar, iwork, nvar, rwork, &
    wk, xr, lrw, liw, icrit, slname )

  iwork(25) = iwork(25) + iwork(28)
!
!  If the corrector succeeds, we're done.  However, before we
!  return, we may have to restore the default values of the
!  continuation parameter and the Newton correction algorithm,
!  because these may have been temporarily altered in a desperate
!  situation!
!
  if ( ierror == 0 ) then
    if ( iwork(3) /= 1 ) then
      iwork(2) = ihold
    end if
    iwork(4) = modsav
    return
  end if
!
!  Only VERY fatal errors should abort the correction process.
!  Abort this process only after kicking and screaming.
!
  if ( ierror == 2 ) then
    if ( iwrite>= 1 ) then
      write ( *, * ) 'TRYSTP - Fatal error during Newton correction'
      write ( *, * ) 'of a continutation point!'
    end if
    return
  end if
!
!  We reach this point if the corrector has failed to converge, but
!  not in a disastrous way.  We have two options: vary the parameter
!  held fixed, or reduce the stepsize.  We prefer to vary the parameter,
!  since that would allow us to take the same "healthy" step.
!
!  We will try varying the parameter if:
!
!    There is a nonzero second choice parameter (JPC); and
!    The corrector was not already using the second choice parameter; and
!    The user did not request via IWORK(3) that we always use a particular
!      parameter.
!
  if ( iwork(3)/= 1.and.ihold/=jpc.and.jpc/=0 ) then
    ihold = jpc
    if ( ierror == 5)go to 20
    go to 10
  end if
!
!  If JPC fails as the index held fixed during the corrector iteration,
!  return to using IPC.
!
  ihold = ipc
!
!  If we are using some modified form of Newton's method, and
!  we have reached the minimum stepsize, or
!  we have had two failures in a row on this step, then
!  we will try a better form of Newton's method.
!
  if ( iwork(4)>0.and. (rwork(5)<rwork(20)*rwork(3).or.iwork(18)>= 2)  ) then
    iwork(4) = iwork(4)-1
    if ( iwrite>= 1 ) then
      write ( *, * ) 'TRYSTP - Retrying step with IWORK(4) = ',iwork(4)
    end if
    if ( ierror/= 5)go to 10
    go to 20
  end if
!
!  No convergence, so reduce the stepsize.  At this point, if the
!  stepsize falls below the user-specified minimum, we have to quit.
!
  if ( rwork(5) < rwork(20) * rwork(3) ) then
    ierror = 4
    if ( 1 <= iwrite ) then
      temp = rwork(5)/rwork(20)
      write ( *, * ) 'TRYSTP - Warning!'
      write ( *, * ) '  The predictor stepsize fell below minimum.'
      write ( *, * ) '  Current step is now ',temp
      write ( *, * ) '  Minimum step is     ',rwork(3)
    end if
    return
  end if

  rwork(5) = rwork(5) / rwork(20)
  iwork(18) = iwork(18) + 1

  if ( ierror /= 5 ) then
    go to 10
  end if
!
!  We're reducing the stepsize, but the corrector iteration was converging,
!  though slowly.  We'll try to salvage the work we had done by using
!  as our new predicted point the linear interpolant between our current
!  accepted point and the corrector iterate that we gave up on.
!
  xr(1:nvar) = xf(1:nvar) + ( xr(1:nvar) - xf(1:nvar) ) / rwork(20)

  go to 20
end
subroutine update ( iwork, liw, lrw, nvar, rwork, tc, xc, xf, xr )

!*****************************************************************************80
!
!! UPDATE updates information after a successful continuation step.
!
!  Modified:
!
!    16 May 2008
!
!  Author:
!
!    John Burkardt
!
  integer liw
  integer lrw
  integer nvar

  integer i
  integer iwork(liw)
  integer iwrite
  double precision rwork(lrw)
  double precision tc(nvar)
  double precision temp
  double precision xc(nvar)
  double precision xf(nvar)
  double precision xr(nvar)

  iwrite = iwork(7)
  iwork(1) = 2
!
!  Note that IWORK(2) may be set to IHOLD here.  This simply reflects
!  the fact that the corrector iteration may have failed with the
!  preferable index, and had to try the "second-best" parameter.
!  Hence, our initial expectation that the parameter would be
!  determined by the maximum entry of the tangent is tempered by our
!  experience with the corrector iteration.
!
  iwork(10) = 2
  iwork(26) = iwork(26) + iwork(18)
  iwork(27) = iwork(27) + 1
  rwork(22) = rwork(21)

  if ( 1 <= iwrite .and. 0 < iwork(18) ) then
    write ( *, * ) ' '
    write ( *, * ) 'UPDATE:'
    write ( *, * ) '  Predictor stepsize reductions: ', iwork(18)
  end if

  temp = 0.0D+00
  do i = 1, nvar
    temp = temp + ( xr(i) - xf(i) )**2
  end do
  temp = sqrt ( temp )
  rwork(21) = temp

  rwork(12) = rwork(13)
  rwork(13) = rwork(12) + rwork(21)
  rwork(14) = rwork(13)
!
!  Compute and store "CORDIS", the maximum-norm of the total
!  correction, MaxNorm(XF+HTAN*TC - XR).
!
  if ( iwork(28) == 0 ) then

   rwork(15) = 0.0D+00

  else

    temp = 0.0D+00
    do i = 1, nvar
      temp = max ( temp, abs ( xr(i) - ( xf(i) + rwork(5) * tc(i) ) ) )
    end do
    rwork(15) = temp

  end if
!
!  Update XC and XF, so that XC has the older continuation point, and
!  XF has the most recent one.
!
  xc(1:nvar) = xf(1:nvar)
  xf(1:nvar) = xr(1:nvar)

  return
end
subroutine dgb_trf ( m, n, ml, mu, a, lda, ipivot, info )

!*****************************************************************************80
!
!! DGB_TRF performs a PLU factorization of an M by N band matrix.
!
!  Discussion:
!
!    DGB_TRF is a simplified, standalone version of the LAPACK
!    routine DGBTRF.
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix A.  M >= 0.
!
!    Input, integer N, the number of columns of the matrix A.  N >= 0.
!
!    Input, integer ML, the number of subdiagonals within the band of A.
!    ML >= 0.
!
!    Input, integer MU, the number of superdiagonals within the band of A.
!    MU >= 0.
!
!    Input/output, double precision A(LDA,N).
!
!    On input, the matrix A in band storage, in rows ML+1 to
!    2*ML+MU+1; rows 1 to ML of the array need not be set.
!    The j-th column of A is stored in the j-th column of the
!    array A as follows:
!    A(ml+mu+1+i-j,j) = A(i,j) for max(1,j-mu) <= i<=min(m,j+ml)
!
!    On exit, details of the factorization: U is stored as an
!    upper triangular band matrix with ML+MU superdiagonals in
!    rows 1 to ML+MU+1, and the multipliers used during the
!    factorization are stored in rows ML+MU+2 to 2*ML+MU+1.
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA >= 2*ML+MU+1.
!
!    Output, integer IPIVOT(min(M,N)), the pivot indices;
!    for 1 <= i <= min(M,N), row i of the matrix was interchanged with
!    row IPIV(i).
!
!    Output, integer INFO, error flag.
!    = 0: successful exit;
!    < 0: an input argument was illegal;
!    > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!         has been completed, but the factor U is exactly
!         singular, and division by zero will occur if it is used
!         to solve a system of equations.
!
  integer lda
  integer m
  integer n

  double precision a(lda,n)
  integer i
  integer info
  integer ipivot(*)
  integer j
  integer jp
  integer ju
  integer k
  integer ml
  integer km
  integer mu
  integer kv
  double precision piv
  double precision temp

  info = 0
!
!  KV is the number of superdiagonals in the factor U, allowing for fill-in.
!
  kv = mu + ml
!
!  Set fill-in elements in columns MU+2 to KV to zero.
!
  do j = mu + 2, min ( kv, n )
    do i = kv - j + 2, ml
      a(i,j) = 0.0D+00
    end do
  end do
!
!  JU is the index of the last column affected by the current stage
!  of the factorization.
!
  ju = 1

  do j = 1, min ( m, n )
!
!  Set the fill-in elements in column J+KV to zero.
!
    if ( j + kv <= n ) then

      do i = 1, ml
        a(i,j+kv ) = 0.0D+00
      end do

    end if
!
!  Find pivot and test for singularity.
!  KM is the number of subdiagonal elements in the current column.
!
    km = min ( ml, m-j )

    piv = abs ( a(kv+1,j) )
    jp = kv+1

    do i = kv + 2, kv + km + 1
      if ( abs ( a(i,j) ) > piv ) then
        piv = abs ( a(i,j ) )
        jp = i
      end if
    end do

    jp = jp - kv

    ipivot(j) = jp + j - 1

    if (  a(kv+jp,j) /= 0.0D+00 ) then

      ju = max ( ju, min ( j+mu+jp-1, n ) )
!
!  Apply interchange to columns J to JU.
!
      if ( jp /= 1 ) then

        do i = 0, ju - j
          temp = a(kv+jp-i,j+i)
          a(kv+jp-i,j+i) = a(kv+1-i,j+i)
          a(kv+1-i,j+i) = temp
        end do

      end if
!
!  Compute multipliers.
!
      if ( km > 0 ) then

        do i = kv+2, kv+km+1
          a(i,j) = a(i,j) / a(kv+1,j)
        end do
!
!  Update the trailing submatrix within the band.
!
        if ( ju > j ) then

          do k = 1, ju-j

            if ( a(kv+1-k,j+k) /= 0.0D+00 ) then

              do i = 1, km

                a(kv+i+1-k,j+k) = a(kv+i+1-k,j+k) - a(kv+i+1,j) * a(kv+1-k,j+k)

              end do
            end if
          end do

        end if

      end if

    else
!
!  If pivot is zero, set INFO to the index of the pivot
!  unless a zero pivot has already been found.
!
      if ( info == 0 ) then
        info = j
      end if

    end if

  end do

  return
end
subroutine dgb_trs ( trans, n, ml, mu, nrhs, a, lda, ipivot, b, ldb, info )

!*****************************************************************************80
!
!! DGB_TRS solves a linear system factored by DGB_TRF.
!
!  Discussion:
!
!    DGB_TRS is a simplified, standalone version of the LAPACK
!    routine DGBTRS.
!
!  Modified:
!
!    19 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*1 TRANS, specifies the form of the system.
!    'N':  A * X = B  (No transpose)
!    'T':  A'* X = B  (Transpose)
!    'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!    Input, integer N, the order of the matrix A.
!    N must be positive.
!
!    Input, integer ML, the number of subdiagonals within the band of A.
!    ML must be at least 0, and no greater than N - 1.
!
!    Input. integer MU, the number of superdiagonals within the band of A.
!    MU must be at least 0, and no greater than N - 1.
!
!    Input, integer NRHS, the number of right hand sides and the number of
!    columns of the matrix B.  NRHS must be positive.
!
!    Input, double precision A(LDA,N), contains the LU factorization of
!    the band matrix A, computed by DGB_TRF.  U is stored as an upper
!    triangular band matrix with ML+MU superdiagonals in rows 1 to ML+MU+1,
!    and the multipliers used during the factorization are stored in
!    rows ML+MU+2 to 2*ML+MU+1.
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer IPIVOT(N), the pivot indices; for 1 <= I <= N, row I
!    of the matrix was interchanged with row IPIVOT(I).
!
!    Input/output, double precision B(LDB,NRHS),
!    On entry, the right hand side vectors B for the system of linear equations.
!    On exit, the solution vectors, X.
!
!    Input, integer LDB, the leading dimension of the array B.
!    LDB must be at least N.
!
!    Output, integer INFO, error flag.
!   = 0:  successful exit
!    < 0: if INFO = -K, the K-th argument had an illegal value
!
  integer lda
  integer ldb
  integer n
  integer nrhs

  double precision a(lda,n)
  double precision b(ldb,nrhs)
  integer i
  integer info
  integer ipivot(*)
  integer j
  integer k
  integer kd
  integer l
  integer lm
  integer ml
  integer mu
  double precision temp
  character trans
!
!  Test the input parameters.
!
  info = 0

  if ( trans /= 'N' .and. trans /= 'n' .and. &
       trans /= 'T' .and. trans /= 't' .and. &
       trans /= 'C' .and. trans /= 'c' ) then
    info = -1
  else if ( n <= 0 ) then
    info = -2
  else if ( ml < 0 ) then
    info = -3
  else if ( mu < 0 ) then
    info = -4
  else if ( nrhs <= 0 ) then
    info = -5
  else if ( lda < ( 2*ml+mu+1 ) ) then
    info = -7
  else if ( ldb < max ( 1, n ) ) then
    info = -10
  end if

  if ( info /= 0 ) then
    return
  end if

  kd = mu + ml + 1
!
!  Solve  A*X = B.
!
!  Solve L*X = B, overwriting B with X.
!
!  L is represented as a product of permutations and unit lower
!  triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!  where each transformation L(i) is a rank-one modification of
!  the identity matrix.
!
  if ( trans == 'N' .or. trans == 'n' ) then

    if ( ml > 0 ) then

      do j = 1, n - 1

        lm = min ( ml, n-j )
        l = ipivot(j)

        do i = 1, nrhs
          temp = b(l,i)
          b(l,i) = b(j,i)
          b(j,i) = temp
        end do

        do k = 1, nrhs
          if ( b(j,k) /= 0.0D+00 ) then
            do i = 1, lm
              b(j+i,k) = b(j+i,k) - a(kd+i,j) * b(j,k)
            end do
          end if
        end do

      end do

    end if
!
!  Solve U*X = B, overwriting B with X.
!
    do i = 1, nrhs

      do j = n, 1, -1
        if ( b(j,i) /= 0.0D+00 ) then
          l = ml + mu + 1 - j
          b(j,i) = b(j,i) / a(ml+mu+1,j)
          do k = j - 1, max ( 1, j - ml - mu ), -1
            b(k,i) = b(k,i) - a(l+k,j) * b(j,i)
          end do
        end if
      end do

    end do

  else
!
!  Solve A'*X = B.
!
!  Solve U'*X = B, overwriting B with X.
!
    do i = 1, nrhs

      do j = 1, n
        temp = b(j,i)
        l = ml + mu + 1 - j
        do k = max ( 1, j - ml - mu ), j - 1
          temp = temp - a(l+k,j) * b(k,i)
        end do
        temp = temp / a(ml+mu+1,j)
        b(j,i) = temp
      end do

    end do
!
!  Solve L'*X = B, overwriting B with X.
!
    if ( ml > 0 ) then

      do j = n - 1, 1, -1

        lm = min ( ml, n-j )

        do k = 1, nrhs
          temp = 0.0D+00
          do i = 1, lm
            temp = temp + b(j+i,k) * a(kd+i,j)
          end do
          b(j,k) = b(j,k) - temp
        end do

        l = ipivot(j)

        do i = 1, nrhs
          temp = b(l,i)
          b(l,i) = b(j,i)
          b(j,i) = temp
        end do

      end do

    end if

  end if

  return
end
function dnrm2 ( n, x, incx )

!*****************************************************************************80
!
!! DNRM2 computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    The original DNRM2 algorithm is accurate but written in a bizarre,
!    unreadable and obsolete format.  This version goes for clarity.
!
!  Modified:
!
!    01 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, double precision X(*), the vector whose norm is to be computed.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, double precision DNRM2, the Euclidean norm of X.
!
  double precision damax
  double precision dnrm2
  integer i
  integer incx
  integer ix
  integer n
  double precision temp
  double precision x(*)
  double precision xmax

  if ( n <= 0 ) then

    dnrm2 = 0.0D+00

  else

    xmax = damax ( n, x, incx )

    if ( xmax == 0.0D+00 ) then

      dnrm2 = 0.0D+00

    else

      if ( incx >= 0 ) then
        ix = 1
      else
        ix = ( - n + 1 ) * incx + 1
      end if

      temp = 0.0D+00
      do i = 1, n
        temp = temp + ( x(ix) / xmax )**2
        ix = ix + incx
      end do

      dnrm2 = xmax * sqrt ( temp )

    end if

  end if

  return
end
function damax ( n, x, incx )

!*****************************************************************************80
!
!! DAMAX returns the maximum absolute value of the entries in a vector.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, double precision X(*), the vector to be examined.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, double precision DAMAX, the maximum absolute value of an
!    element of X.
!
  double precision damax
  integer i
  integer incx
  integer ix
  integer n
  double precision x(*)

  if ( n <= 0 ) then

    damax = 0.0D+00

  else if ( n == 1 ) then

    damax = abs ( x(1) )

  else if ( incx == 1 ) then

    damax = abs ( x(1) )

    do i = 2, n
      if ( abs ( x(i) ) > damax ) then
        damax = abs ( x(i) )
      end if
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    damax = abs ( x(ix) )
    ix = ix + incx

    do i = 2, n
      if ( abs ( x(ix) ) > damax ) then
        damax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
function idamax ( n, x, incx )

!*****************************************************************************80
!
!! IDAMAX finds the index of the vector element of maximum absolute value.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, double precision X(*), the vector to be examined.
!
!    Input, integer INCX, the increment between successive entries of SX.
!
!    Output, integer IDAMAX, the index of the element of SX of maximum
!    absolute value.
!
  integer i
  integer incx
  integer idamax
  integer ix
  integer n
  double precision damax
  double precision x(*)

  if ( n <= 0 ) then

    idamax = 0

  else if ( n == 1 ) then

    idamax = 1

  else if ( incx == 1 ) then

    idamax = 1
    damax = abs ( x(1) )

    do i = 2, n

      if ( abs ( x(i) ) > damax ) then
        idamax = i
        damax = abs ( x(i) )
      end if

    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    idamax = 1
    damax = abs ( x(ix) )

    ix = ix + incx

    do i = 2, n
      if ( abs ( x(ix) ) > damax ) then
        idamax = i
        damax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

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
!    31 May 2001   9:45:54.872 AM
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
