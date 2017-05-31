!###############################################################################
!#                                                                             #
!# aed2_vode.F90                                                               #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created June 2012                                                           #
!#                                                                             #
!###############################################################################

!------------------------------------------------------------------------------!
MODULE aed2_vode

   USE aed2_gctypes,only :SEDP,aed2_sed_candi_t

   IMPLICIT NONE

   PRIVATE
!  PUBLIC :: DVODE, dgefa, dgedia
   PUBLIC :: DVODE, dgefa


CONTAINS


!------------------------------------------------------------------------------!
!    D VODE Variable-coefficient Ordinary Differential Equation solver,        !
!          with fixed-leading coefficient implementation.                      !
!          This is a REAL(SEDP) version of VODE.                               !
!                                                                              !
!         DVODE solves the initial value problem for stiff or nonstiff         !
!         systems of first order ODEs,                                         !
!                                                                              !
!         dy/dt = f(t,y) ,  or, in component form,                             !
!         dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).       !
!                                                                              !
!         DVODE is a package based on the EPISODE and EPISODEB packages,       !
!         and on the ODEPACK user interface standard, with minor               !
!         modifications.                                                       !
!-----------------------------------------------------------------------       !
! Revision History (YYMMDD)                                                    !
!   890615  Date Written                                                       !
!   890922  Added interrupt/restart ability, minor changes throughout.         !
!   910228  Minor revisions in line format,  prologue, etc.                    !
!   920227  Modifications by D. Pang:                                          !
!           (1) Applied subgennam to get generic intrinsic names.              !
!           (2) Changed intrinsic names to generic in comments.                !
!           (3) Added *DECK lines before each routine.                         !
!   920721  Names of routines and labeled Common blocks changed, so as         !
!           to be unique in combined single/REAL(SEDP) code (ACH).             !
!   920722  Minor revisions to prologue (ACH).                                 !
!   920831  Conversion to REAL  done (ACH).                                    !
!-----------------------------------------------------------------------       !
! References..                                                                 !
                                                                               !
! 1. P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, VODE: A Variable           !
!    Coefficient ODE Solver, SIAM J. Sci. Stat. Comput., 10 (1989),            !
!    pp. 1038-1051.  Also, LLNL Report UCRL-98412, June 1988.                  !
! 2. G. D. Byrne and A. C. Hindmarsh, A Polyalgorithm for the                  !
!    Numerical Solution of Ordinary Differential Equations,                    !
!    ACM Trans. Math. Software, 1 (1975), pp. 71-96.                           !
! 3. A. C. Hindmarsh and G. D. Byrne, EPISODE: An Effective Package            !
!    for the Integration of Systems of Ordinary Differential                   !
!    Equations, LLNL Report UCID-30112, Rev. 1, April 1977.                    !
! 4. G. D. Byrne and A. C. Hindmarsh, EPISODEB: An Experimental                !
!    Package for the Integration of Systems of Ordinary Differential           !
!    Equations with Banded Jacobians, LLNL Report UCID-30132, April            !
!    1976.                                                                     !
! 5. A. C. Hindmarsh, ODEPACK, a Systematized Collection of ODE                !
!    Solvers, in Scientific Computing, R. S. Stepleman et al., eds.,           !
!    North-Holland, Amsterdam, 1983, pp. 55-64.                                !
! 6. K. R. Jackson and R. Sacks-Davis, An Alternative Implementation           !
!    of Variable Step-Size Multistep Formulas for Stiff ODEs, ACM              !
!    Trans. Math. Software, 6 (1980), pp. 295-318.                             !
!-----------------------------------------------------------------------       !
! Authors..                                                                    !
                                                                               !
!               Peter N. Brown and Alan C. Hindmarsh                           !
!               Computing and Mathematics Research Division, L-316             !
!               Lawrence Livermore National Laboratory                         !
!               Livermore, CA 94550                                            !
! and                                                                          !
!               George D. Byrne                                                !
!               Exxon Research and Engineering Co.                             !
!               Clinton Township                                               !
!               Route 22 East                                                  !
!               Annandale, NJ 08801                                            !
!-----------------------------------------------------------------------       !
! Summary of usage.                                                            !
                                                                               !
! Communication between the user and the DVODE package, for normal             !
! situations, is summarized here.  This summary describes only a subset        !
! of the full set of options available.  See the full description for          !
! details, including optional communication, nonstandard options,              !
! and instructions for special situations.  See also the example               !
! problem (with program and output) following this summary.                    !
                                                                               !
! A. First provide a subroutine of the form..                                  !
                                                                               !
!           SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)                         !
!           REAL(SEDP) T, Y, YDOT, RPAR                                        !
!           DIMENSION Y(NEQ), YDOT(NEQ)                                        !
                                                                               !
! which supplies the vector function f by loading YDOT(i) with f(i).           !
                                                                               !
! B. Next determine (or guess) whether or not the problem is stiff.            !
! Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue            !
! whose REAL(SEDP) part is negative and large in magnitude, compared to the    !
! reciprocal of the t span of interest.  If the problem is nonstiff,           !
! use a method flag MF = 10.  If it is stiff, there are four standard          !
! choices for MF (21, 22, 24, 25), and DVODE requires the Jacobian             !
! matrix in some form.  In these cases (MF .gt. 0), DVODE will use a           !
! saved copy of the Jacobian matrix.  If this is undesirable because of        !
! storage limitations, set MF to the corresponding negative value              !
! (-21, -22, -24, -25).  (See full description of MF below.)                   !
! The Jacobian matrix is regarded either as full (MF = 21 or 22),              !
! or banded (MF = 24 or 25).  In the banded case, DVODE requires two           !
! half-bandwidth parameters ML and MU.  These are, respectively, the           !
! widths of the lower and upper parts of the band, excluding the main          !
! diagonal.  Thus the band consists of the locations (i,j) with                !
! i-ML .le. j .le. i+MU, and the full bandwidth is ML+MU+1.                    !
                                                                               !
! C. If the problem is stiff, you are encouraged to supply the Jacobian        !
! directly (MF = 21 or 24), but if this is not feasible, DVODE will            !
! compute it internally by difference quotients (MF = 22 or 25).               !
! If you are supplying the Jacobian, provide a subroutine of the form..        !
                                                                               !
!           SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)         !
!           REAL(SEDP) T, Y, PD, RPAR                                          !
!           DIMENSION Y(NEQ), PD(NROWPD,NEQ)                                   !
                                                                               !
! which supplies df/dy by loading PD as follows..                              !
!     For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),            !
! the partial derivative of f(i) with respect to y(j).  (Ignore the            !
! ML and MU arguments in this case.)                                           !
!     For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with                !
! df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of          !
! PD from the top down.                                                        !
!     In either case, only nonzero elements need be loaded.                    !
                                                                               !
! D. Write a main program which calls subroutine DVODE once for                !
! each point at which answers are desired.  This should also provide           !
! for possible use of logical unit 6 for output of error messages              !
! by DVODE.  On the first CALL to DVODE, supply arguments as follows..         !
! F      = Name of subroutine for right-hand side vector f.                    !
!          This name must be declared external in calling program.             !
! NEQ    = Number of first order ODE-s.                                        !
! Y      = Array of initial values, of length NEQ.                             !
! T      = The initial value of the independent variable.                      !
! TOUT   = First point where output is desired (.ne. T).                       !
! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.              !
! RTOL   = Relative tolerance parameter (scalar).                              !
! ATOL   = Absolute tolerance parameter (scalar or array).                     !
!          The estimated local error in Y(i) will be controlled so as          !
!          to be roughly less (in magnitude) than                              !
!             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or               !
!             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2.                  !
!          Thus the local error test passes if, in each component,             !
!          either the absolute error is less than ATOL (or ATOL(i)),           !
!          or the relative error is less than RTOL.                            !
!          Use RTOL = 0.0 for pure absolute error control, and                 !
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error           !
!          control.  Caution.. Actual (global) errors may exceed these         !
!          local tolerances, so choose them conservatively.                    !
! ITASK  = 1 for normal computation of output values of Y at t = TOUT.         !
! ISTATE = Integer flag (input and output).  Set ISTATE = 1.                   !
! IOPT   = 0 to indicate no optional input used.                               !
! RWORK  = FLOAT work array of length at least..                               !
!             20 + 16*NEQ                      for MF = 10,                    !
!             22 +  9*NEQ + 2*NEQ**2           for MF = 21 or 22,              !
!             22 + 11*NEQ + (3*ML + 2*MU)*NEQ  for MF = 24 or 25.              !
! LRW    = Declared length of RWORK (in users DIMENSION statement).            !
! IWORK  = Integer work array of length at least..                             !
!             30        for MF = 10,                                           !
!             30 + NEQ  for MF = 21, 22, 24, or 25.                            !
!          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower              !
!          and upper half-bandwidths ML,MU.                                    !
! LIW    = Declared length of IWORK (in users DIMENSION).                      !
! JAC    = Name of subroutine for Jacobian matrix (MF = 21 or 24).             !
!          If used, this name must be declared external in calling             !
!          program.  If not used, pass a dummy name.                           !
! MF     = Method flag.  Standard values are..                                 !
!          10 for nonstiff (Adams) method, no Jacobian used.                   !
!          21 for stiff (BDF) method, user-supplied full Jacobian.             !
!          22 for stiff method, internally generated full Jacobian.            !
!          24 for stiff method, user-supplied banded Jacobian.                 !
!          25 for stiff method, internally generated banded Jacobian.          !
! RPAR,IPAR = user-defined REAL(SEDP) and integer arrays passed to F and JAC.  !
! Note that the main program must declare arrays Y, RWORK, IWORK,              !
! and possibly ATOL, RPAR, and IPAR.                                           !
                                                                               !
! E. The output from the first CALL (or any call) is..                         !
!      Y = Array of computed values of y(t) vector.                            !
!      T = Corresponding value of independent variable (normally TOUT).        !
! ISTATE = 2  if DVODE was successful, negative otherwise.                     !
!          -1 means excess work done on this CALL. (Perhaps wrong MF.)         !
!          -2 means excess accuracy requested. (Tolerances too small.)         !
!          -3 means illegal input detected. (See printed message.)             !
!          -4 means repeated error test failures. (Check all input.)           !
!          -5 means repeated convergence failures. (Perhaps bad                !
!             Jacobian supplied or wrong choice of MF or tolerances.)          !
!          -6 means error weight became zero during problem. (Solution         !
!             component i vanished, and ATOL or ATOL(i) = 0.)                  !
                                                                               !
! F. To continue the integration after a successful return, simply             !
! reset TOUT and CALL DVODE again.  No other parameters need be reset.         !
                                                                               !
!-----------------------------------------------------------------------       !
! Other Routines in the DVODE Package.                                         !
                                                                               !
! In addition to subroutine DVODE, the DVODE package includes the              !
! following subroutines and function routines..                                !
!  DVHIN     computes an approximate step size for the initial step.           !
!  DVINDY    computes an interpolated value of the y vector at t = TOUT.       !
!  DVSTEP    is the core integrator, which does one step of the                !
!            integration and the associated error control.                     !
!  DVSET     sets all method coefficients and test constants.                  !
!  DVNLSD    solves the underlying nonlinear system -- the corrector.          !
!  DVJAC     computes and preprocesses the Jacobian matrix J = df/dy           !
!            and the Newton iteration matrix P = I - (h/l1)*J.                 !
!  DVSOL     manages solution of linear system in chord iteration.             !
!  DVJUST    adjusts the history array on a change of order.                   !
!  DEWSET    sets the error weight vector EWT before each step.                !
!  DVNORM    computes the weighted r.m.s. norm of a vector.                    !
!  DVSRCO    is a user-callable routines to save and restore                   !
!            the contents of the internal COMMON blocks.                       !
!  DACOPY    is a routine to copy one two-dimensional array to another.        !
!  DGEFA and DGESL   are routines from LINPACK for solving full                !
!            systems of linear algebraic equations.                            !
!  DGBFA and DGBSL   are routines from LINPACK for solving banded              !
!            linear systems.                                                   !
!  DAXPY, DSCAL, and DCOPY are basic linear algebra modules (BLAS).            !
!  D1MACH    sets the unit roundoff of the machine.                            !
!  XERRWD, XSETUN, XSETF, LUNSAV, and MFLGSV handle the printing of all        !
!            error messages and warnings.  XERRWD is machine-dependent.        !
! Note..  DVNORM, D1MACH, LUNSAV, and MFLGSV are function routines.            !
! All the others are subroutines.                                              !
                                                                               !
! The intrinsic and external routines used by the DVODE package are..          !
! ABS, MAX, MIN, FLOAT, SIGN, SQRT, and WRITE.                                 !
                                                                               !
!------------------------------------------------------------------------------!
 SUBROUTINE dvode (candi, f, neq, y, t, tout, itol, rtol, atol, itask,         &
                   istate, iopt, rwork, lrw, iwork, liw, jac, mf,              &
                   rpar, ipar)

   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout):: candi
   INTEGER,    INTENT(IN)                      :: neq
   INTEGER,    INTENT(IN)                      :: lrw
   INTEGER,    INTENT(IN)                      :: liw
   REAL(SEDP), INTENT(INOUT)                   :: y(neq)
   REAL(SEDP), INTENT(INOUT)                   :: t
   REAL(SEDP), INTENT(IN)                      :: tout
   INTEGER,    INTENT(IN)                      :: itol
   REAL(SEDP), INTENT(IN)                      :: rtol(lrw)
   REAL(SEDP), INTENT(IN)                      :: atol(liw)
   INTEGER,    INTENT(IN)                      :: itask
   INTEGER,    INTENT(INOUT)                   :: istate
   INTEGER,    INTENT(IN)                      :: iopt
   REAL(SEDP), INTENT(INOUT)                   :: rwork(lrw)
   INTEGER,    INTENT(INOUT)                   :: iwork(liw)
   INTEGER,    INTENT(INOUT)                   :: mf
   INTEGER,    INTENT(INOUT)                   :: ipar(*)
   REAL(SEDP), DIMENSION(:), INTENT(INOUT)     :: rpar
   EXTERNAL f, jac


   ! Type declarations for labeled COMMON block DVOD01 --------------------

   REAL(SEDP) :: acnrm, ccmxj, conp, crate, drc, el,                           &
                 eta, etamax, h, hmin, hmxi, hnew, hscal, prl1, rc, rl1,       &
                 tau, tq, tn, uround
   INTEGER  :: icf, init, ipup, jcur, jstart, jsv, kflag, kuth,              &
                 l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                  &
                 locjs, maxord, meth, miter, msbj, mxhnil, mxstep,             &
                 n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh

   ! Type declarations for labeled COMMON block DVOD02 --------------------

   REAL(SEDP) :: hu
   INTEGER  :: ncfn, netf, nfe, nje, nlu, nni, nqu, nst

   ! Type declarations for local variables --------------------------------

   LOGICAL    :: ihit
   REAL(SEDP) :: atoli, big, ewti, four, h0, hmax, hmx, hun, one,              &
                 pt2, rh, rtoli, size, tcrit, tnext, tolsf, tp, two, zero
   INTEGER  :: i, ier, iflag, imxer, jco, kgo, leniw, lenj, lenp, lenrw,     &
                 lenwm, lf0, mband, ml, mord, mu, mxhnl0, mxstp0, niter, nslast
   CHARACTER (LEN=80) :: msg

   ! Type declaration for function subroutines called ---------------------

   !REAL(SEDP) :: d1mach, dvnorm

   DIMENSION mord(2)
   SAVE mord, mxhnl0, mxstp0
   SAVE zero, one, two, four, pt2, hun
   COMMON /dvod01/ acnrm, ccmxj, conp, crate, drc, el(13),                     &
                   eta, etamax, h, hmin, hmxi, hnew, hscal, prl1,              &
                   rc, rl1, tau(13), tq(5), tn, uround,                        &
                   icf, init, ipup, jcur, jstart, jsv, kflag, kuth,            &
                   l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                &
                   locjs, maxord, meth, miter, msbj, mxhnil, mxstep,           &
                   n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh
   COMMON /dvod02/ hu, ncfn, netf, nfe, nje, nlu, nni, nqu, nst

   DATA  mord(1) /12/, mord(2) /5/, mxstp0 /500/, mxhnl0 /10/
   DATA zero /0.0D+00/, one /1.0D+00/, two /2.0D+00/, four /4.0D+00/,          &
        pt2 /0.2D+00/, hun /100.0D+00/

   !-----------------------------------------------------------------------
   ! Block A.
   ! This code block is executed on every CALL.
   ! It tests ISTATE and ITASK for legality and branches appropriately.
   ! If ISTATE .gt. 1 but the flag INIT shows that initialization has
   ! not yet been done, an error return occurs.
   ! If ISTATE = 1 and TOUT = T, return immediately.
   !-----------------------------------------------------------------------
   IF (istate < 1 .OR. istate > 3) GO TO 601
   IF (itask < 1 .OR. itask > 5) GO TO 602
   IF (istate == 1) GO TO 10
   IF (init /= 1) GO TO 603
   IF (istate == 2) GO TO 200
   GO TO 20
   10 init = 0
   IF (tout == t) RETURN
   !-----------------------------------------------------------------------
   ! Block B.
   ! The next code block is executed for the initial CALL (ISTATE = 1),
   ! or for a continuation CALL with parameter changes (ISTATE = 3).
   ! It contains checking of all input and various initializations.

   ! First check legality of the non-optional input NEQ, ITOL, IOPT,
   ! MF, ML, and MU.
   !-----------------------------------------------------------------------
   20 IF (neq <= 0) GO TO 604
   IF (istate == 1) GO TO 25
   IF (neq > n) GO TO 605
   25 n = neq
   IF (itol < 1 .OR. itol > 4) GO TO 606
   IF (iopt < 0 .OR. iopt > 1) GO TO 607
   jsv = SIGN(1,mf)
   mf = ABS(mf)
   meth = mf/10
   miter = mf - 10*meth
   IF (meth < 1 .OR. meth > 2) GO TO 608
   IF (miter < 0 .OR. miter > 5) GO TO 608
   IF (miter <= 3) GO TO 30
   ml = iwork(1)
   mu = iwork(2)
   IF (ml < 0 .OR. ml >= n) GO TO 609
   IF (mu < 0 .OR. mu >= n) GO TO 610
   30 CONTINUE
   ! Next process and check the optional input. ---------------------------
   IF (iopt == 1) GO TO 40
   maxord = mord(meth)
   mxstep = mxstp0
   mxhnil = mxhnl0
   IF (istate == 1) h0 = zero
   hmxi = zero
   hmin = zero
   GO TO 60
   40 maxord = iwork(5)
   IF (maxord < 0) GO TO 611
   IF (maxord == 0) maxord = 100
   maxord = MIN(maxord,mord(meth))
   mxstep = iwork(6)
   IF (mxstep < 0) GO TO 612
   IF (mxstep == 0) mxstep = mxstp0
   mxhnil = iwork(7)
   IF (mxhnil < 0) GO TO 613
   IF (mxhnil == 0) mxhnil = mxhnl0
   IF (istate /= 1) GO TO 50
   h0 = rwork(5)
   IF ((tout - t)*h0 < zero) GO TO 614
   50 hmax = rwork(6)
   IF (hmax < zero) GO TO 615
   hmxi = zero
   IF (hmax > zero) hmxi = one/hmax
   hmin = rwork(7)
   IF (hmin < zero) GO TO 616
   !-----------------------------------------------------------------------
   ! Set work array pointers and check lengths LRW and LIW.
   ! Pointers to segments of RWORK and IWORK are named by prefixing L to
   ! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
   ! Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
   ! Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
   !-----------------------------------------------------------------------
   60 lyh = 21
   IF (istate == 1) nyh = n
   lwm = lyh + (maxord + 1)*nyh
   jco = MAX(0,jsv)
   IF (miter == 0) lenwm = 0
   IF (miter == 1 .OR. miter == 2) THEN
     lenwm = 2 + (1 + jco)*n*n
     locjs = n*n + 3
   END IF
   IF (miter == 3) lenwm = 2 + n
   IF (miter == 4 .OR. miter == 5) THEN
     mband = ml + mu + 1
     lenp = (mband + ml)*n
     lenj = mband*n
     lenwm = 2 + lenp + jco*lenj
     locjs = lenp + 3
   END IF
   lewt = lwm + lenwm
   lsavf = lewt + n
   lacor = lsavf + n
   lenrw = lacor + n - 1
   iwork(17) = lenrw
   liwm = 1
   leniw = 30 + n
   IF (miter == 0 .OR. miter == 3) leniw = 30
   iwork(18) = leniw
   IF (lenrw > lrw) GO TO 617
   IF (leniw > liw) GO TO 618
   ! Check RTOL and ATOL for legality. ------------------------------------
   rtoli = rtol(1)
   atoli = atol(1)
   DO  i = 1,n
     IF (itol >= 3) rtoli = rtol(i)
     IF (itol == 2 .OR. itol == 4) atoli = atol(i)
     IF (rtoli < zero) GO TO 619
     IF (atoli < zero) GO TO 620
   END DO
   IF (istate == 1) GO TO 100
   ! If ISTATE = 3, set flag to signal parameter changes to DVSTEP. -------
   jstart = -1
   IF (nq <= maxord) GO TO 90
   ! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
   !CALL dcopy (n, rwork(lwm), 1, rwork(lsavf), 1)
   do i=0,n-1
     rwork(lsavf+i)=rwork(lwm+i)
   enddo
   ! Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
   90 IF (miter > 0) rwork(lwm) = SQRT(uround)
   !-----------------------------------------------------------------------
   ! Block C.
   ! The next block is for the initial CALL only (ISTATE = 1).
   ! It contains all remaining initializations, the initial CALL to F,
   ! and the calculation of the initial step size.
   ! The error weights in EWT are inverted after being loaded.
   !-----------------------------------------------------------------------
   100 uround = d1mach(4)
   tn = t
   IF (itask /= 4 .AND. itask /= 5) GO TO 110
   tcrit = rwork(1)
   IF ((tcrit - tout)*(tout - t) < zero) GO TO 625
   IF (h0 /= zero .AND. (t + h0 - tcrit)*h0 > zero) h0 = tcrit - t
   110 jstart = 0
   IF (miter > 0) rwork(lwm) = SQRT(uround)
   ccmxj = pt2
   msbj = 50
   nhnil = 0
   nst = 0
   nje = 0
   nni = 0
   ncfn = 0
   netf = 0
   nlu = 0
   nslj = 0
   nslast = 0
   hu = zero
   nqu = 0
   ! Initial CALL to F.  (LF0 points to YH(*,2).) -------------------------
   lf0 = lyh + nyh

   CALL f (candi, n, t, y, rwork(lf0), rpar, ipar)  ! Here is where FEX is called
   ! Wow how on earth do you know that? - Dan

   nfe = 1
   ! Load the initial value vector in YH. ---------------------------------
   CALL dcopy (n, y, 1, rwork(lyh), 1)
   ! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
   nq = 1
   h = one
   CALL dewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
   DO  i = 1,n
     IF (rwork(i+lewt-1) <= zero) GO TO 621
     rwork(i+lewt-1) = one/rwork(i+lewt-1)
   END DO
   IF (h0 /= zero) GO TO 180
   ! Call DVHIN to set initial step size H0 to be attempted. --------------
   CALL dvhin (candi, n, t, rwork(lyh), rwork(lf0), f, rpar, ipar, tout,  &
       uround, rwork(lewt), itol, atol, y, rwork(lacor), h0, niter, ier)
   nfe = nfe + niter
   IF (ier /= 0) GO TO 622
   ! Adjust H0 if necessary to meet HMAX bound. ---------------------------
   180 rh = ABS(h0)*hmxi
   IF (rh > one) h0 = h0/rh
   ! Load H with H0 and scale YH(*,2) by H0. ------------------------------
   h = h0
   !CALL dscal (n, h0, rwork(lf0), 1)
   do i=0,n-1
     rwork(lf0+i)=rwork(lf0+i)*h0
   enddo
   GO TO 270
   !-----------------------------------------------------------------------
   ! Block D.
   ! The next code block is for continuation calls only (ISTATE = 2 or 3)
   ! and is to check stop conditions before taking a step.
   !-----------------------------------------------------------------------
   200 nslast = nst
   kuth = 0
   SELECT CASE ( itask )
     CASE (    1)
       GO TO 210
     CASE (    2)
       GO TO  250
     CASE (    3)
       GO TO  220
     CASE (    4)
       GO TO  230
     CASE (    5)
       GO TO  240
   END SELECT
   210 IF ((tn - tout)*h < zero) GO TO 250
   CALL dvindy (tout, 0, rwork(lyh), nyh, y, iflag)
   IF (iflag /= 0) GO TO 627
   t = tout
   GO TO 420
   220 tp = tn - hu*(one + hun*uround)
   IF ((tp - tout)*h > zero) GO TO 623
   IF ((tn - tout)*h < zero) GO TO 250
   GO TO 400
   230 tcrit = rwork(1)
   IF ((tn - tcrit)*h > zero) GO TO 624
   IF ((tcrit - tout)*h < zero) GO TO 625
   IF ((tn - tout)*h < zero) GO TO 245
   CALL dvindy (tout, 0, rwork(lyh), nyh, y, iflag)
   IF (iflag /= 0) GO TO 627
   t = tout
   GO TO 420
   240 tcrit = rwork(1)
   IF ((tn - tcrit)*h > zero) GO TO 624
   245 hmx = ABS(tn) + ABS(h)
   ihit = ABS(tn - tcrit) <= hun*uround*hmx
   IF (ihit) GO TO 400
   tnext = tn + hnew*(one + four*uround)
   IF ((tnext - tcrit)*h <= zero) GO TO 250
   h = (tcrit - tn)*(one - four*uround)
   kuth = 1
   !-----------------------------------------------------------------------
   ! Block E.
   ! The next block is normally executed for all calls and contains
   ! the CALL to the one-step core integrator DVSTEP.

   ! This is a looping point for the integration steps.

   ! First check for too many steps being taken, update EWT (if not at
   ! start of problem), check for too much accuracy being requested, and
   ! check for H below the roundoff level in T.
   !-----------------------------------------------------------------------
   250 CONTINUE
   IF ((nst-nslast) >= mxstep) GO TO 500
   CALL dewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
   DO  i = 1,n
     IF (rwork(i+lewt-1) <= zero) GO TO 510
     rwork(i+lewt-1) = one/rwork(i+lewt-1)
   END DO
   270 tolsf = uround*dvnorm (n, rwork(lyh), rwork(lewt))
   IF (tolsf <= one) GO TO 280
   tolsf = tolsf*two
   IF (nst == 0) GO TO 626
   GO TO 520
   280 IF ((tn + h) /= tn) GO TO 290
   nhnil = nhnil + 1
   IF (nhnil > mxhnil) GO TO 290
   msg = 'DVODE--  Warning..internal T (=R1) and H (=R2) are'
   CALL xerrwd (msg, 50, 101, 1, 0, 0, 0, 0, zero, zero)
   msg='      such that in the machine, T + H = T on the next step  '
   CALL xerrwd (msg, 60, 101, 1, 0, 0, 0, 0, zero, zero)
   msg = '      (H = step size). solver will continue anyway'
   CALL xerrwd (msg, 50, 101, 1, 0, 0, 0, 2, tn, h)
   IF (nhnil < mxhnil) GO TO 290
   msg = 'DVODE--  Above warning has been issued I1 times.  '
   CALL xerrwd (msg, 50, 102, 1, 0, 0, 0, 0, zero, zero)
   msg = '      it will not be issued again for this problem'
   CALL xerrwd (msg, 50, 102, 1, 1, mxhnil, 0, 0, zero, zero)
   290 CONTINUE
   !-----------------------------------------------------------------------
   ! CALL DVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR,
   !              WM, IWM, F, JAC, F, DVNLSD, RPAR, IPAR)
   !-----------------------------------------------------------------------
   CALL dvstep (candi, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),  &
       rwork(lsavf), y, rwork(lacor), rwork(lwm), iwork(liwm),  &
       f, jac, f, dvnlsd, rpar, ipar)
   kgo = 1 - kflag
   ! Branch on KFLAG.  Note..In this version, KFLAG can not be set to -3.
   !  KFLAG .eq. 0,   -1,  -2
   SELECT CASE ( kgo )
     CASE (    1)
       GO TO 300
     CASE (    2)
       GO TO  530
     CASE (    3)
       GO TO  540
   END SELECT
   !-----------------------------------------------------------------------
   ! Block F.
   ! The following block handles the case of a successful return from the
   ! core integrator (KFLAG = 0).  Test for stop conditions.
   !-----------------------------------------------------------------------
   300 init = 1
   kuth = 0
   SELECT CASE ( itask )
     CASE (    1)
       GO TO 310
     CASE (    2)
       GO TO  400
     CASE (    3)
       GO TO  330
     CASE (    4)
       GO TO  340
     CASE (    5)
       GO TO  350
   END SELECT
   ! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
   310 IF ((tn - tout)*h < zero) GO TO 250
   CALL dvindy (tout, 0, rwork(lyh), nyh, y, iflag)
   t = tout
   GO TO 420
   ! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
   330 IF ((tn - tout)*h >= zero) GO TO 400
   GO TO 250
   ! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
   340 IF ((tn - tout)*h < zero) GO TO 345
   CALL dvindy (tout, 0, rwork(lyh), nyh, y, iflag)
   t = tout
   GO TO 420
   345 hmx = ABS(tn) + ABS(h)
   ihit = ABS(tn - tcrit) <= hun*uround*hmx
   IF (ihit) GO TO 400
   tnext = tn + hnew*(one + four*uround)
   IF ((tnext - tcrit)*h <= zero) GO TO 250
   h = (tcrit - tn)*(one - four*uround)
   kuth = 1
   GO TO 250
   ! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
   350 hmx = ABS(tn) + ABS(h)
   ihit = ABS(tn - tcrit) <= hun*uround*hmx
   !-----------------------------------------------------------------------
   ! Block G.
   ! The following block handles all successful returns from DVODE.
   ! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
   ! ISTATE is set to 2, and the optional output is loaded into the work
   ! arrays before returning.
   !-----------------------------------------------------------------------
   400 CONTINUE
   CALL dcopy (n, rwork(lyh), 1, y, 1)
   t = tn
   IF (itask /= 4 .AND. itask /= 5) GO TO 420
   IF (ihit) t = tcrit
   420 istate = 2
   rwork(11) = hu
   rwork(12) = hnew
   rwork(13) = tn
   iwork(11) = nst
   iwork(12) = nfe
   iwork(13) = nje
   iwork(14) = nqu
   iwork(15) = newq
   iwork(19) = nlu
   iwork(20) = nni
   iwork(21) = ncfn
   iwork(22) = netf
   RETURN
   !-----------------------------------------------------------------------
   ! Block H.
   ! The following block handles all unsuccessful returns other than
   ! those for illegal input.  First the error message routine is called.
   ! if there was an error test or convergence test failure, IMXER is set.
   ! Then Y is loaded from YH, T is set to TN, and the illegal input
   ! The optional output is loaded into the work arrays before returning.
   !-----------------------------------------------------------------------
   ! The maximum number of steps was taken before reaching TOUT. ----------
   500  msg = 'DVODE--  At current T (=R1), MXSTEP (=I1) steps   '
   CALL xerrwd (msg, 50, 201, 1, 0, 0, 0, 0, zero, zero)
   msg = '      taken on this CALL before reaching TOUT     '
   CALL xerrwd (msg, 50, 201, 1, 1, mxstep, 0, 1, tn, zero)
   istate = -1
   GO TO 580
   ! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
   510 ewti = rwork(lewt+i-1)
   msg = 'DVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.'
   CALL xerrwd (msg, 50, 202, 1, 1, i, 0, 2, tn, ewti)
   istate = -6
   GO TO 580
   ! Too much accuracy requested for machine precision. -------------------
   520 msg = 'DVODE--  At T (=R1), too much accuracy requested  '
   CALL xerrwd (msg, 50, 203, 1, 0, 0, 0, 0, zero, zero)
   msg = '      for precision of machine..  see TOLSF (=R2) '
   CALL xerrwd (msg, 50, 203, 1, 0, 0, 0, 2, tn, tolsf)
   rwork(14) = tolsf
   istate = -2
   GO TO 580
   ! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
   530 msg = 'DVODE--  At T(=R1) and step size H(=R2), the error'
   CALL xerrwd (msg, 50, 204, 1, 0, 0, 0, 0, zero, zero)
   msg = '      test failed repeatedly or with abs(H) = HMIN'
   CALL xerrwd (msg, 50, 204, 1, 0, 0, 0, 2, tn, h)
   istate = -4
   GO TO 560
   ! KFLAG = -2.  Convergence failed repeatedly or with abs(H) = HMIN. ----
   540  msg = 'DVODE--  At T (=R1) and step size H (=R2), the    '
   CALL xerrwd (msg, 50, 205, 1, 0, 0, 0, 0, zero, zero)
   msg = '      corrector convergence failed repeatedly     '
   CALL xerrwd (msg, 50, 205, 1, 0, 0, 0, 0, zero, zero)
   msg = '      or with abs(H) = HMIN   '
   CALL xerrwd (msg, 30, 205, 1, 0, 0, 0, 2, tn, h)
   istate = -5
   ! Compute IMXER if relevant. -------------------------------------------
   560 big = zero
   imxer = 1
   DO  i = 1,n
     size = ABS(rwork(i+lacor-1)*rwork(i+lewt-1))
     IF (big >= size) CYCLE
     big = size
     imxer = i
   END DO
   iwork(16) = imxer
   ! Set Y vector, T, and optional output. --------------------------------
   580 CONTINUE
   CALL dcopy (n, rwork(lyh), 1, y, 1)
   t = tn
   rwork(11) = hu
   rwork(12) = h
   rwork(13) = tn
   iwork(11) = nst
   iwork(12) = nfe
   iwork(13) = nje
   iwork(14) = nqu
   iwork(15) = nq
   iwork(19) = nlu
   iwork(20) = nni
   iwork(21) = ncfn
   iwork(22) = netf
   RETURN
   !-----------------------------------------------------------------------
   ! Block I.
   ! The following block handles all error returns due to illegal input
   ! (ISTATE = -3), as detected before calling the core integrator.
   ! First the error message routine is called.   If the illegal input
   ! is a negative ISTATE, the run is aborted (apparent infinite loop).
   !-----------------------------------------------------------------------
   601 msg = 'DVODE--  ISTATE (=I1) illegal '
   CALL xerrwd (msg, 30, 1, 1, 1, istate, 0, 0, zero, zero)
   IF (istate < 0) GO TO 800
   GO TO 700
   602 msg = 'DVODE--  ITASK (=I1) illegal  '
   CALL xerrwd (msg, 30, 2, 1, 1, itask, 0, 0, zero, zero)
   GO TO 700
   603 msg='DVODE--  ISTATE (=I1) .gt. 1 but DVODE not initialized      '
   CALL xerrwd (msg, 60, 3, 1, 1, istate, 0, 0, zero, zero)
   GO TO 700
   604 msg = 'DVODE--  NEQ (=I1) .lt. 1     '
   CALL xerrwd (msg, 30, 4, 1, 1, neq, 0, 0, zero, zero)
   GO TO 700
   605 msg = 'DVODE--  ISTATE = 3 and NEQ increased (I1 to I2)  '
   CALL xerrwd (msg, 50, 5, 1, 2, n, neq, 0, zero, zero)
   GO TO 700
   606 msg = 'DVODE--  ITOL (=I1) illegal   '
   CALL xerrwd (msg, 30, 6, 1, 1, itol, 0, 0, zero, zero)
   GO TO 700
   607 msg = 'DVODE--  IOPT (=I1) illegal   '
   CALL xerrwd (msg, 30, 7, 1, 1, iopt, 0, 0, zero, zero)
   GO TO 700
   608 msg = 'DVODE--  MF (=I1) illegal     '
   CALL xerrwd (msg, 30, 8, 1, 1, mf, 0, 0, zero, zero)
   GO TO 700
   609 msg = 'DVODE--  ML (=I1) illegal.. .lt.0 or .ge.NEQ (=I2)'
   CALL xerrwd (msg, 50, 9, 1, 2, ml, neq, 0, zero, zero)
   GO TO 700
   610 msg = 'DVODE--  MU (=I1) illegal.. .lt.0 or .ge.NEQ (=I2)'
   CALL xerrwd (msg, 50, 10, 1, 2, mu, neq, 0, zero, zero)
   GO TO 700
   611 msg = 'DVODE--  MAXORD (=I1) .lt. 0  '
   CALL xerrwd (msg, 30, 11, 1, 1, maxord, 0, 0, zero, zero)
   GO TO 700
   612 msg = 'DVODE--  MXSTEP (=I1) .lt. 0  '
   CALL xerrwd (msg, 30, 12, 1, 1, mxstep, 0, 0, zero, zero)
   GO TO 700
   613 msg = 'DVODE--  MXHNIL (=I1) .lt. 0  '
   CALL xerrwd (msg, 30, 13, 1, 1, mxhnil, 0, 0, zero, zero)
   GO TO 700
   614 msg = 'DVODE--  TOUT (=R1) behind T (=R2)      '
   CALL xerrwd (msg, 40, 14, 1, 0, 0, 0, 2, tout, t)
   msg = '      integration direction is given by H0 (=R1)  '
   CALL xerrwd (msg, 50, 14, 1, 0, 0, 0, 1, h0, zero)
   GO TO 700
   615 msg = 'DVODE--  HMAX (=R1) .lt. 0.0  '
   CALL xerrwd (msg, 30, 15, 1, 0, 0, 0, 1, hmax, zero)
   GO TO 700
   616 msg = 'DVODE--  HMIN (=R1) .lt. 0.0  '
   CALL xerrwd (msg, 30, 16, 1, 0, 0, 0, 1, hmin, zero)
   GO TO 700
   617 CONTINUE
   msg='DVODE--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
   CALL xerrwd (msg, 60, 17, 1, 2, lenrw, lrw, 0, zero, zero)
   GO TO 700
   618 CONTINUE
   msg='DVODE--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
   CALL xerrwd (msg, 60, 18, 1, 2, leniw, liw, 0, zero, zero)
   GO TO 700
   619 msg = 'DVODE--  RTOL(I1) is R1 .lt. 0.0        '
   CALL xerrwd (msg, 40, 19, 1, 1, i, 0, 1, rtoli, zero)
   GO TO 700
   620 msg = 'DVODE--  ATOL(I1) is R1 .lt. 0.0        '
   CALL xerrwd (msg, 40, 20, 1, 1, i, 0, 1, atoli, zero)
   GO TO 700
   621 ewti = rwork(lewt+i-1)
   msg = 'DVODE--  EWT(I1) is R1 .le. 0.0         '
   CALL xerrwd (msg, 40, 21, 1, 1, i, 0, 1, ewti, zero)
   GO TO 700
   622 CONTINUE
   msg='DVODE--  TOUT (=R1) too close to T(=R2) to start integration'
   CALL xerrwd (msg, 60, 22, 1, 0, 0, 0, 2, tout, t)
   GO TO 700
   623 CONTINUE
   msg='DVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
   CALL xerrwd (msg, 60, 23, 1, 1, itask, 0, 2, tout, tp)
   GO TO 700
   624 CONTINUE
   msg='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
   CALL xerrwd (msg, 60, 24, 1, 0, 0, 0, 2, tcrit, tn)
   GO TO 700
   625 CONTINUE
   msg='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
   CALL xerrwd (msg, 60, 25, 1, 0, 0, 0, 2, tcrit, tout)
   GO TO 700
   626 msg = 'DVODE--  At start of problem, too much accuracy   '
   CALL xerrwd (msg, 50, 26, 1, 0, 0, 0, 0, zero, zero)
   msg='      requested for precision of machine..  see TOLSF (=R1) '
   CALL xerrwd (msg, 60, 26, 1, 0, 0, 0, 1, tolsf, zero)
   rwork(14) = tolsf
   GO TO 700
   627 msg='DVODE--  Trouble from DVINDY.  ITASK = I1, TOUT = R1.       '
   CALL xerrwd (msg, 60, 27, 1, 1, itask, 0, 1, tout, zero)

   700  CONTINUE
   istate = -3
   RETURN

   800 msg = 'DVODE--  Run aborted.. apparent infinite loop     '
   CALL xerrwd (msg, 50, 303, 2, 0, 0, 0, 0, zero, zero)
   RETURN
   !----------------------- End of Subroutine DVODE -----------------------
 END SUBROUTINE dvode
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! dvhin                                                                        !
! Call sequence input -- N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,         !
!                        EWT, ITOL, ATOL, Y, TEMP                              !
! Call sequence output -- H0, NITER, IER                                       !
! COMMON block variables accessed -- None                                      !
                                                                               !
! Subroutines called by DVHIN.. F                                              !
! Function routines called by DVHIN.. DVNORM                                   !
!-----------------------------------------------------------------------       !
! This routine computes the step size, H0, to be attempted on the              !
! first step, when the user has not supplied a value for this.                 !
                                                                               !
! First we check that TOUT - T0 differs significantly from zero.  Then         !
! an iteration is done to approximate the initial second derivative            !
! and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.        !
! A bias factor of 1/2 is applied to the resulting h.                          !
! The sign of H0 is inferred from the initial values of TOUT and T0.           !
                                                                               !
! Communication with DVHIN is done with the following variables..              !
                                                                               !
! N      = Size of ODE system, input.                                          !
! T0     = Initial value of independent variable, input.                       !
! Y0     = Vector of initial conditions, input.                                !
! YDOT   = Vector of initial first derivatives, input.                         !
! F      = Name of subroutine for right-hand side f(t,y), input.               !
! RPAR, IPAR = Dummy names for users REAL(SEDP) and integer work arrays.       !
! TOUT   = First output value of independent variable                          !
! UROUND = Machine unit roundoff                                               !
! EWT, ITOL, ATOL = Error weights and tolerance parameters                     !
!                   as described in the driver routine, input.                 !
! Y, TEMP = Work arrays of length N.                                           !
! H0     = Step size to be attempted, output.                                  !
! NITER  = Number of iterations (and of f evaluations) to compute H0,          !
!          output.                                                             !
! IER    = The error flag, returned with the value                             !
!          IER = 0  if no trouble occurred, or                                 !
!          IER = -1 if TOUT and T0 are considered too close to proceed.        !
!------------------------------------------------------------------------------!
 SUBROUTINE dvhin (candi, n, t0, y0, ydot, f, rpar, ipar, tout, uround,        &
                   ewt, itol, atol, y, temp, h0, niter, ier)

   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi

   INTEGER,    INTENT(IN)             :: n
   REAL(SEDP), INTENT(IN)             :: t0
   REAL(SEDP), INTENT(IN)             :: y0(*)
   REAL(SEDP), INTENT(IN)             :: ydot(*)
   REAL(SEDP), INTENT(INOUT)          :: rpar(*)
   INTEGER,    INTENT(INOUT)          :: ipar(*)
   REAL(SEDP), INTENT(IN)             :: tout
   REAL(SEDP), INTENT(IN)             :: uround
   REAL(SEDP), INTENT(INOUT)          :: ewt(*)
   INTEGER,    INTENT(IN)             :: itol
   REAL(SEDP), INTENT(IN)             :: atol(*)
   REAL(SEDP), INTENT(OUT)            :: y(*)
   REAL(SEDP), INTENT(OUT)            :: temp(*)
   REAL(SEDP), INTENT(OUT)            :: h0
   INTEGER,    INTENT(OUT)            :: niter
   INTEGER,    INTENT(OUT)            :: ier
   EXTERNAL f

   !-----------------------------------------------------------------------
   ! Type declarations for local variables --------------------------------

   REAL(SEDP) :: afi, atoli, delyi, hg, hlb, hnew, hrat,                       &
                 hub, t1, tdist, tround, yddnrm

   !REAL(SEDP) :: half,hun,pt1,two
   INTEGER  :: i
   INTEGER  :: iter

   ! Type declaration for function subroutines called ---------------------
   !REAL(SEDP) :: dvnorm
   !-----------------------------------------------------------------------
   ! The following Fortran-77 declaration is to cause the values of the
   ! listed (local) variables to be saved between calls to this integrator.
   !-----------------------------------------------------------------------
   !SAVE half, hun, pt1, two
   !DATA half /0.5D+00/, hun /100.0D+00/, pt1 /0.1D+00/, two /2.0D+00/

   niter = 0
   tdist = ABS(tout - t0)
   tround = uround*MAX(ABS(t0),ABS(tout))
   IF (tdist < 2.0*tround) GO TO 100

   ! Set a lower bound on h based on the roundoff level in T0 and TOUT. ---
   hlb = 100.0*tround
   ! Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT. -
   hub = 0.5*tdist
   atoli = atol(1)
   !print *, "REIN",n,atoli,itol,hub
   DO  i = 1, n
   ! Leave the print inside to aviod problems on a NEC-Bullshit mainfraim
     if (i .eq. 0) print *, n,atoli,itol,hub,hlb,tround
     IF (itol == 2 .OR. itol == 4) atoli = atol(i)
     delyi = 0.5*ABS(y0(i)) + atoli
     afi = ABS(ydot(i))
     IF (afi*hub > delyi) hub = delyi/afi
   END DO
   ! Set initial guess for h as geometric mean of upper and lower bounds. -
   iter = 0
   hg = SQRT(hlb*hub)
   ! If the bounds have crossed, exit with the mean value. ----------------
   IF (hub < hlb) THEN
     h0 = hg
     GO TO 90
   END IF

   ! Looping point for iteration. -----------------------------------------
   50  CONTINUE
   ! Estimate the second derivative as a difference quotient in f. --------
   t1 = t0 + hg
   DO  i = 1, n
     y(i) = y0(i) + hg*ydot(i)
   END DO
   CALL f (candi, n, t1, y, temp, rpar, ipar)
   DO  i = 1, n
     temp(i) = (temp(i) - ydot(i))/hg
   END DO
   yddnrm = dvnorm (n, temp, ewt)
   ! Get the corresponding new value of h. --------------------------------
   IF (yddnrm*hub*hub > 2.0) THEN
     hnew = SQRT(2.0/yddnrm)
   ELSE
     hnew = SQRT(hg*hub)
   END IF
   iter = iter + 1
   !-----------------------------------------------------------------------
   ! Test the stopping conditions.
   ! Stop if the new and previous h values differ by a factor of .lt. 2.
   ! Stop if four iterations have been done.  Also, stop with previous h
   ! if HNEW/HG .gt. 2 after first iteration, as this probably means that
   ! the second derivative value is bad because of cancellation error.
   !-----------------------------------------------------------------------
   IF (iter >= 4) GO TO 80
   hrat = hnew/hg
   IF ( (hrat > 0.5) .AND. (hrat < 2.0) ) GO TO 80
   IF ( (iter >= 2) .AND. (hnew > 2.0*hg) ) THEN
     hnew = hg
     GO TO 80
   END IF
   hg = hnew
   GO TO 50

   ! Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ----
   80 h0 = hnew*0.5
   IF (h0 < hlb) h0 = hlb
   IF (h0 > hub) h0 = hub
   90 h0 = SIGN(h0, tout - t0)
   niter = iter
   ier = 0
   RETURN
   ! Error return for TOUT - T0 too small. --------------------------------
   100 ier = -1
   RETURN
   !----------------------- End of Subroutine DVHIN -----------------------
 END SUBROUTINE dvhin
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! dvindy                                                                       !
! Call sequence input -- T, K, YH, LDYH                                        !
! Call sequence output -- DKY, IFLAG                                           !
! COMMON block variables accessed..                                            !
!     /DVOD01/ --  H, TN, UROUND, L, N, NQ                                     !
!     /DVOD02/ --  HU                                                          !
                                                                               !
! Subroutines called by DVINDY.. DSCAL, XERRWD                                 !
! Function routines called by DVINDY.. None                                    !
!-----------------------------------------------------------------------       !
! DVINDY computes interpolated values of the K-th derivative of the            !
! dependent variable vector y, and stores it in DKY.  This routine             !
! is called within the package with K = 0 and T = TOUT, but may                !
! also be called by the user for any K up to the current order.                !
! (See detailed instructions in the usage documentation.)                      !
!-----------------------------------------------------------------------       !
! The computed values in DKY are gotten by interpolation using the             !
! Nordsieck history array YH.  This array corresponds uniquely to a            !
! vector-valued polynomial of degree NQCUR or less, and DKY is set             !
! to the K-th derivative of this polynomial at T.                              !
! The formula for DKY is..                                                     !
!              q                                                               !
!  DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)              !
!             j=K                                                              !
! where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.         !
! The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are                       !
! communicated by COMMON.  The above sum is done in reverse order.             !
! IFLAG is returned negative if either K or T is out of bounds.                !
                                                                               !
! Discussion above and comments in driver explain all variables.               !
!------------------------------------------------------------------------------!
 SUBROUTINE dvindy (t, k, yh, ldyh, dky, iflag)

   INTEGER,    INTENT(IN)                      :: k
   INTEGER,    INTENT(INOUT)                   :: ldyh
   REAL(SEDP), INTENT(IN)                      :: yh(ldyh,*)
   REAL(SEDP), INTENT(IN)                      :: t
   REAL(SEDP), INTENT(OUT)                     :: dky(*)
   INTEGER,    INTENT(OUT)                     :: iflag

   !-----------------------------------------------------------------------
   ! Type declarations for labeled COMMON block DVOD01 --------------------

   REAL(SEDP) :: acnrm, ccmxj, conp, crate, drc, el,                           &
                 eta, etamax, h, hmin, hmxi, hnew, hscal, prl1, rc, rl1,       &
                 tau, tq, tn, uround

   INTEGER  :: icf, init, ipup, jcur, jstart, jsv, kflag, kuth,              &
                 l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                  &
                 locjs, maxord, meth, miter, msbj, mxhnil, mxstep,             &
                 n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh

   ! Type declarations for labeled COMMON block DVOD02 --------------------

   REAL(SEDP) :: hu
   INTEGER  :: ncfn, netf, nfe, nje, nlu, nni, nqu, nst

   ! Type declarations for local variables --------------------------------

   REAL(SEDP) :: c, hun, r, s, tfuzz, tn1, tp, zero
   INTEGER  :: i, ic, j, jb, jb2, jj, jj1, jp1
   CHARACTER (LEN=80) :: msg

   !-----------------------------------------------------------------------
   ! The following Fortran-77 declaration is to cause the values of the
   ! listed (local) variables to be saved between calls to this integrator.
   !-----------------------------------------------------------------------
   SAVE hun, zero

   COMMON /dvod01/ acnrm, ccmxj, conp, crate, drc, el(13),                     &
       eta, etamax, h, hmin, hmxi, hnew, hscal, prl1,                          &
       rc, rl1, tau(13), tq(5), tn, uround,                                    &
       icf, init, ipup, jcur, jstart, jsv, kflag, kuth,                        &
       l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                            &
       locjs, maxord, meth, miter, msbj, mxhnil, mxstep,                       &
       n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh
   COMMON /dvod02/ hu, ncfn, netf, nfe, nje, nlu, nni, nqu, nst

   DATA hun /100.0D+00/, zero /0.0D+00/

   iflag = 0
   IF (k < 0 .OR. k > nq) GO TO 80
   tfuzz = hun*uround*(tn + hu)
   tp = tn - hu - tfuzz
   tn1 = tn + tfuzz
   IF ((t-tp)*(t-tn1) > zero) GO TO 90

   s = (t - tn)/h
   ic = 1
   IF (k == 0) GO TO 15
   jj1 = l - k
   DO  jj = jj1, nq
     ic = ic*jj
   END DO
   15 c = FLOAT(ic)
   DO  i = 1, n
     dky(i) = c*yh(i,l)
   END DO
   IF (k == nq) GO TO 55
   jb2 = nq - k
   DO  jb = 1, jb2
     j = nq - jb
     jp1 = j + 1
     ic = 1
     IF (k == 0) GO TO 35
     jj1 = jp1 - k
     DO  jj = jj1, j
       ic = ic*jj
     END DO
     35   c = FLOAT(ic)
     DO  i = 1, n
       dky(i) = c*yh(i,jp1) + s*dky(i)
     END DO
   END DO
   IF (k == 0) RETURN
   55 r = h**(-k)
   !CALL dscal (n, r, dky, 1)
   do i=1,n
     dky(i)=dky(i)*r
   enddo
   RETURN

   80 msg = 'DVINDY-- K (=I1) illegal      '
   CALL xerrwd (msg, 30, 51, 1, 1, k, 0, 0, zero, zero)
   iflag = -1
   RETURN
   90 msg = 'DVINDY-- T (=R1) illegal      '
   CALL xerrwd (msg, 30, 52, 1, 0, 0, 0, 1, t, zero)
   msg='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
   CALL xerrwd (msg, 60, 52, 1, 0, 0, 0, 2, tp, tn)
   iflag = -2
   RETURN
   !----------------------- End of Subroutine DVINDY ----------------------
 END SUBROUTINE dvindy
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! dvstep                                                                       !
! Call sequence input -- Y, YH, LDYH, YH1, EWT, SAVF, VSAV,                    !
!                        ACOR, WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR         !
! Call sequence output -- YH, ACOR, WM, IWM                                    !
! COMMON block variables accessed..                                            !
!     /DVOD01/  ACNRM, EL(13), H, HMIN, HMXI, HNEW, HSCAL, RC, TAU(13),        !
!               TQ(5), TN, JCUR, JSTART, KFLAG, KUTH,                          !
!               L, LMAX, MAXORD, MITER, N, NEWQ, NQ, NQWAIT                    !
!     /DVOD02/  HU, NCFN, NETF, NFE, NQU, NST                                  !
                                                                               !
! Subroutines called by DVSTEP.. F, DAXPY, DCOPY, DSCAL,                       !
!                               DVJUST, VNLS, DVSET                            !
! Function routines called by DVSTEP.. DVNORM                                  !
!-----------------------------------------------------------------------       !
! DVSTEP performs one step of the integration of an initial value              !
! problem for a system of ordinary differential equations.                     !
! DVSTEP calls subroutine VNLS for the solution of the nonlinear system        !
! arising in the time step.  Thus it is independent of the problem             !
! Jacobian structure and the type of nonlinear system solution method.         !
! DVSTEP returns a completion flag KFLAG (in COMMON).                          !
! A return with KFLAG = -1 or -2 means either ABS(H) = HMIN or 10              !
! consecutive failures occurred.  On a return with KFLAG negative,             !
! the values of TN and the YH array are as of the beginning of the last        !
! step, and H is the last step size attempted.                                 !
                                                                               !
! Communication with DVSTEP is done with the following variables..             !
                                                                               !
! Y      = An array of length N used for the dependent variable vector.        !
! YH     = An LDYH by LMAX array containing the dependent variables            !
!          and their approximate scaled derivatives, where                     !
!          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate              !
!          j-th derivative of y(i), scaled by H**j/factorial(j)                !
!          (j = 0,1,...,NQ).  On entry for the first step, the first           !
!          two columns of YH must be set from the initial values.              !
! LDYH   = A constant integer .ge. N, the first dimension of YH.               !
!          N is the number of ODEs in the system.                              !
! YH1    = A one-dimensional array occupying the same space as YH.             !
! EWT    = An array of length N containing multiplicative weights              !
!          for local error measurements.  Local errors in y(i) are             !
!          compared to 1.0/EWT(i) in various error tests.                      !
! SAVF   = An array of working storage, of length N.                           !
!          also used for input of YH(*,MAXORD+2) when JSTART = -1              !
!          and MAXORD .lt. the current order NQ.                               !
! VSAV   = A work array of length N passed to subroutine VNLS.                 !
! ACOR   = A work array of length N, used for the accumulated                  !
!          corrections.  On a successful return, ACOR(i) contains              !
!          the estimated one-step local error in y(i).                         !
! WM,IWM = FLOAT and integer work arrays associated with matrix                !
!          operations in VNLS.                                                 !
! F      = Dummy name for the user supplied subroutine for f.                  !
! JAC    = Dummy name for the user supplied Jacobian subroutine.               !
! PSOL   = Dummy name for the subroutine passed to VNLS, for                   !
!          possible use there.                                                 !
! VNLS   = Dummy name for the nonlinear system solving subroutine,             !
!          whose REAL(SEDP) name is dependent on the method used.              !
! RPAR, IPAR = Dummy names for users REAL(SEDP) and integer work arrays.       !
!------------------------------------------------------------------------------!
 SUBROUTINE dvstep (candi, y, yh, ldyh, yh1, ewt, savf, vsav, acor,                   &
                    wm, iwm, f, jac, psol, vnls, rpar, ipar)

   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   REAL(SEDP), INTENT(INOUT)         :: y(*)
   INTEGER,    INTENT(IN)            :: ldyh
   REAL(SEDP), INTENT(OUT)           :: yh(ldyh,*)
   REAL(SEDP), INTENT(OUT)           :: yh1(*)
   REAL(SEDP), INTENT(INOUT)         :: ewt(*)
   REAL(SEDP), INTENT(INOUT)         :: savf(*)
   REAL(SEDP), INTENT(INOUT)         :: vsav(*)
   REAL(SEDP), INTENT(INOUT)         :: acor(*)
   REAL(SEDP), INTENT(INOUT)         :: wm(*)
   INTEGER,    INTENT(INOUT)         :: iwm(*)
   REAL(SEDP), INTENT(INOUT)         :: rpar(*)
   INTEGER,    INTENT(INOUT)         :: ipar(*)
   EXTERNAL f, jac, psol, vnls

   !-----------------------------------------------------------------------
   ! Type declarations for labeled COMMON block DVOD01 --------------------

   REAL(SEDP) :: acnrm, ccmxj, conp, crate, drc, el,                           &
                 eta, etamax, h, hmin, hmxi, hnew, hscal, prl1, rc, rl1,       &
                 tau, tq, tn, uround
   INTEGER  :: icf, init, ipup, jcur, jstart, jsv, kflag, kuth,              &
                 l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                  &
                 locjs, maxord, meth, miter, msbj, mxhnil, mxstep,             &
                 n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh

   ! Type declarations for labeled COMMON block DVOD02 --------------------

   REAL(SEDP) :: hu
   INTEGER  :: ncfn, netf, nfe, nje, nlu, nni, nqu, nst

   ! Type declarations for local variables --------------------------------

   REAL(SEDP) :: addon, bias1,bias2,bias3, cnquot, ddn, dsm, dup,              &
                 etacf, etamin, etamx1, etamx2, etamx3, etamxf,                &
                 etaq, etaqm1, etaqp1, flotl, one, onepsm, r, thresh, told, zero
   INTEGER  :: i, i1, i2, iback, j, jb, kfc, kfh, mxncf, ncf, nflag

   ! Type declaration for function subroutines called ---------------------

   !REAL(SEDP) :: dvnorm

   !-----------------------------------------------------------------------
   ! The following Fortran-77 declaration is to cause the values of the
   ! listed (local) variables to be saved between calls to this integrator.
   !-----------------------------------------------------------------------
   SAVE addon, bias1, bias2, bias3,  &
       etacf, etamin, etamx1, etamx2, etamx3, etamxf,  &
       kfc, kfh, mxncf, onepsm, thresh, one, zero
   !-----------------------------------------------------------------------
   COMMON /dvod01/ acnrm, ccmxj, conp, crate, drc, el(13),  &
       eta, etamax, h, hmin, hmxi, hnew, hscal, prl1,  &
       rc, rl1, tau(13), tq(5), tn, uround,  &
       icf, init, ipup, jcur, jstart, jsv, kflag, kuth,  &
       l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,  &
       locjs, maxord, meth, miter, msbj, mxhnil, mxstep,  &
       n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh
   COMMON /dvod02/ hu, ncfn, netf, nfe, nje, nlu, nni, nqu, nst

   DATA kfc/-3/, kfh/-7/, mxncf/10/
   DATA addon  /1.0D-06/,   bias1  /6.0D+00/,     bias2  /6.0D+00/,  &
       bias3  /10.0D+00/,  etacf  /0.25D+00/,    etamin /0.1D+00/,  &
       etamxf /0.2D+00/,   etamx1 /1.0D+04/,     etamx2 /10.0D+00/,  &
       etamx3 /10.0D+00/,  onepsm /1.00001D+00/, thresh /1.5D+00/
   DATA one/1.0D+00/, zero/0.0D+00/

   kflag = 0
   told = tn
   ncf = 0
   jcur = 0
   nflag = 0
   IF (jstart > 0) GO TO 20
   IF (jstart == -1) GO TO 100
   !-----------------------------------------------------------------------
   ! On the first CALL, the order is set to 1, and other variables are
   ! initialized.  ETAMAX is the maximum ratio by which H can be increased
   ! in a single step.  It is normally 1.5, but is larger during the
   ! first 10 steps to compensate for the small initial H.  If a failure
   ! occurs (in corrector convergence or error test), ETAMAX is set to 1
   ! for the next increase.
   !-----------------------------------------------------------------------
   lmax = maxord + 1
   nq = 1
   l = 2
   nqnyh = nq*ldyh
   tau(1) = h
   prl1 = one
   rc = zero
   etamax = etamx1
   nqwait = 2
   hscal = h
   GO TO 200
   !-----------------------------------------------------------------------
   ! Take preliminary actions on a normal continuation step (JSTART.GT.0).
   ! If the driver changed H, then ETA must be reset and NEWH set to 1.
   ! If a change of order was dictated on the previous step, then
   ! it is done here and appropriate adjustments in the history are made.
   ! On an order decrease, the history array is adjusted by DVJUST.
   ! On an order increase, the history array is augmented by a column.
   ! On a change of step size H, the history array YH is rescaled.
   !-----------------------------------------------------------------------
   20 CONTINUE
   IF (kuth == 1) THEN
     eta = MIN(eta,h/hscal)
     newh = 1
   END IF
   50 IF (newh == 0) GO TO 200
   IF (newq == nq) GO TO 150
   IF (newq < nq) THEN
     CALL dvjust (yh, ldyh, -1)
     nq = newq
     l = nq + 1
     nqwait = l
     GO TO 150
   END IF
   IF (newq > nq) THEN
     CALL dvjust (yh, ldyh, 1)
     nq = newq
     l = nq + 1
     nqwait = l
     GO TO 150
   END IF
   !-----------------------------------------------------------------------
   ! The following block handles preliminaries needed when JSTART = -1.
   ! If N was reduced, zero out part of YH to avoid undefined references.
   ! If MAXORD was reduced to a value less than the tentative order NEWQ,
   ! then NQ is set to MAXORD, and a new H ratio ETA is chosen.
   ! Otherwise, we take the same preliminary actions as for JSTART .gt. 0.
   ! In any case, NQWAIT is reset to L = NQ + 1 to prevent further
   ! changes in order for that many steps.
   ! The new H ratio ETA is limited by the input H if KUTH = 1,
   ! by HMIN if KUTH = 0, and by HMXI in any case.
   ! Finally, the history array YH is rescaled.
   !-----------------------------------------------------------------------
   100 CONTINUE
   lmax = maxord + 1
   IF (n == ldyh) GO TO 120
   i1 = 1 + (newq + 1)*ldyh
   i2 = (maxord + 1)*ldyh
   IF (i1 > i2) GO TO 120
   DO  i = i1, i2
     yh1(i) = zero
   END DO
   120 IF (newq <= maxord) GO TO 140
   flotl = FLOAT(lmax)
   IF (maxord < nq-1) THEN
     ddn = dvnorm (n, savf, ewt)/tq(1)
     eta = one/((bias1*ddn)**(one/flotl) + addon)
   END IF
   IF (maxord == nq .AND. newq == nq+1) eta = etaq
   IF (maxord == nq-1 .AND. newq == nq+1) THEN
     eta = etaqm1
     CALL dvjust (yh, ldyh, -1)
   END IF
   IF (maxord == nq-1 .AND. newq == nq) THEN
     ddn = dvnorm (n, savf, ewt)/tq(1)
     eta = one/((bias1*ddn)**(one/flotl) + addon)
     CALL dvjust (yh, ldyh, -1)
   END IF
   eta = MIN(eta,one)
   nq = maxord
   l = lmax
   140 IF (kuth == 1) eta = MIN(eta,ABS(h/hscal))
   IF (kuth == 0) eta = MAX(eta,hmin/ABS(hscal))
   eta = eta/MAX(one,ABS(hscal)*hmxi*eta)
   newh = 1
   nqwait = l
   IF (newq <= maxord) GO TO 50
   ! Rescale the history array for a change in H by a factor of ETA. ------
   150 r = one
   DO  j = 2, l
     r = r*eta
   !  CALL dscal (n, r, yh(1,j), 1 )
     do i=1,n
       yh(i,j)=yh(i,j)*r
     enddo
   END DO
   h = hscal*eta
   hscal = h
   rc = rc*eta
   nqnyh = nq*ldyh
   !-----------------------------------------------------------------------
   ! This section computes the predicted values by effectively
   ! multiplying the YH array by the Pascal triangle matrix.
   ! DVSET is called to calculate all integration coefficients.
   ! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
   !-----------------------------------------------------------------------
   200 tn = tn + h
   i1 = nqnyh + 1
   DO  jb = 1, nq
     i1 = i1 - ldyh
     DO  i = i1, nqnyh
       yh1(i) = yh1(i) + yh1(i+ldyh)
     END DO
   END DO
   CALL dvset
   rl1 = one/el(2)
   rc = rc*(rl1/prl1)
   prl1 = rl1

   ! Call the nonlinear system solver. ------------------------------------

   CALL vnls (candi, y, yh, ldyh, vsav, savf, ewt, acor, iwm, wm,  &
       f, jac, psol, nflag, rpar, ipar)

   IF (nflag == 0) GO TO 450
   !-----------------------------------------------------------------------
   ! The VNLS routine failed to achieve convergence (NFLAG .NE. 0).
   ! The YH array is retracted to its values before prediction.
   ! The step size H is reduced and the step is retried, if possible.
   ! Otherwise, an error exit is taken.
   !-----------------------------------------------------------------------
   ncf = ncf + 1
   ncfn = ncfn + 1
   etamax = one
   tn = told
   i1 = nqnyh + 1
   DO  jb = 1, nq
     i1 = i1 - ldyh
     DO  i = i1, nqnyh
       yh1(i) = yh1(i) - yh1(i+ldyh)
     END DO
   END DO
   IF (nflag < -1) GO TO 680
   IF (ABS(h) <= hmin*onepsm) GO TO 670
   IF (ncf == mxncf) GO TO 670
   eta = etacf
   eta = MAX(eta,hmin/ABS(h))
   nflag = -1
   GO TO 150
   !-----------------------------------------------------------------------
   ! The corrector has converged (NFLAG = 0).  The local error test is
   ! made and control passes to statement 500 if it fails.
   !-----------------------------------------------------------------------
   450 CONTINUE
   dsm = acnrm/tq(2)
   IF (dsm > one) GO TO 500
   !-----------------------------------------------------------------------
   ! After a successful step, update the YH and TAU arrays and decrement
   ! NQWAIT.  If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved
   ! for use in a possible order increase on the next step.
   ! If ETAMAX = 1 (a failure occurred this step), keep NQWAIT .ge. 2.
   !-----------------------------------------------------------------------
   kflag = 0
   nst = nst + 1
   hu = h
   nqu = nq
   DO  iback = 1, nq
     i = l - iback
     tau(i+1) = tau(i)
   END DO
   tau(1) = h
   DO  j = 1, l
   !  CALL daxpy (n, el(j), acor, yh(1,j))
   DO  i = 1,n
     yh(i,j)=yh(i,j)+el(j)*acor(i)
   END DO
   END DO
   nqwait = nqwait - 1
   IF ((l == lmax) .OR. (nqwait /= 1)) GO TO 490
   !CALL dcopy (n, acor, 1, yh(1,lmax), 1 )
   do i=1,n
     yh(i,lmax)=acor(i)
   enddo
   conp = tq(5)
   490 IF (etamax /= one) GO TO 560
   IF (nqwait < 2) nqwait = 2
   newq = nq
   newh = 0
   eta = one
   hnew = h
   GO TO 690
   !-----------------------------------------------------------------------
   ! The error test failed.  KFLAG keeps track of multiple failures.
   ! Restore TN and the YH array to their previous values, and prepare
   ! to try the step again.  Compute the optimum step size for the
   ! same order.  After repeated failures, H is forced to decrease
   ! more rapidly.
   !-----------------------------------------------------------------------
   500 kflag = kflag - 1
   netf = netf + 1
   nflag = -2
   tn = told
   i1 = nqnyh + 1
   DO  jb = 1, nq
     i1 = i1 - ldyh
     DO  i = i1, nqnyh
       yh1(i) = yh1(i) - yh1(i+ldyh)
     END DO
   END DO
   IF (ABS(h) <= hmin*onepsm) GO TO 660
   etamax = one
   IF (kflag <= kfc) GO TO 530
   ! Compute ratio of new H to current H at the current order. ------------
   flotl = FLOAT(l)
   eta = one/((bias2*dsm)**(one/flotl) + addon)
   eta = MAX(eta,hmin/ABS(h),etamin)
   IF ((kflag <= -2) .AND. (eta > etamxf)) eta = etamxf
   GO TO 150
   !-----------------------------------------------------------------------
   ! Control reaches this section if 3 or more consecutive failures
   ! have occurred.  It is assumed that the elements of the YH array
   ! have accumulated errors of the wrong order.  The order is reduced
   ! by one, if possible.  Then H is reduced by a factor of 0.1 and
   ! the step is retried.  After a total of 7 consecutive failures,
   ! an exit is taken with KFLAG = -1.
   !-----------------------------------------------------------------------
   530 IF (kflag == kfh) GO TO 660
   IF (nq == 1) GO TO 540
   eta = MAX(etamin,hmin/ABS(h))
   CALL dvjust (yh, ldyh, -1)
   l = nq
   nq = nq - 1
   nqwait = l
   GO TO 150
   540 eta = MAX(etamin,hmin/ABS(h))
   h = h*eta
   hscal = h
   tau(1) = h
   CALL f (candi, n, tn, y, savf, rpar, ipar)
   nfe = nfe + 1
   DO  i = 1, n
     yh(i,2) = h*savf(i)
   END DO
   nqwait = 10
   GO TO 200
   !-----------------------------------------------------------------------
   ! If NQWAIT = 0, an increase or decrease in order by one is considered.
   ! Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
   ! be multiplied at order q, q-1, or q+1, respectively.
   ! The largest of these is determined, and the new order and
   ! step size set accordingly.
   ! A change of H or NQ is made only if H increases by at least a
   ! factor of THRESH.  If an order change is considered and rejected,
   ! then NQWAIT is set to 2 (reconsider it after 2 steps).
   !-----------------------------------------------------------------------
   ! Compute ratio of new H to current H at the current order. ------------
   560 flotl = FLOAT(l)
   etaq = one/((bias2*dsm)**(one/flotl) + addon)
   IF (nqwait /= 0) GO TO 600
   nqwait = 2
   etaqm1 = zero
   IF (nq == 1) GO TO 570
   ! Compute ratio of new H to current H at the current order less one. ---
   ddn = dvnorm (n, yh(1,l), ewt)/tq(1)
   etaqm1 = one/((bias1*ddn)**(one/(flotl - one)) + addon)
   570 etaqp1 = zero
   IF (l == lmax) GO TO 580
   ! Compute ratio of new H to current H at current order plus one. -------
   cnquot = (tq(5)/conp)*(h/tau(2))**l
   DO  i = 1, n
     savf(i) = acor(i) - cnquot*yh(i,lmax)
   END DO
   dup = dvnorm (n, savf, ewt)/tq(3)
   etaqp1 = one/((bias3*dup)**(one/(flotl + one)) + addon)
   580  IF (etaq >= etaqp1) GO TO 590
   IF (etaqp1 > etaqm1) GO TO 620
   GO TO 610
   590  IF (etaq < etaqm1) GO TO 610
   600  eta = etaq
   newq = nq
   GO TO 630
   610  eta = etaqm1
   newq = nq - 1
   GO TO 630
   620  eta = etaqp1
   newq = nq + 1
   !CALL dcopy (n, acor, 1, yh(1,lmax), 1)
   do i=1,n
     yh(i,lmax)=acor(i)
   enddo
   ! Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ----
   630  IF (eta < thresh .OR. etamax == one) GO TO 640
   eta = MIN(eta,etamax)
   eta = eta/MAX(one,ABS(h)*hmxi*eta)
   newh = 1
   hnew = h*eta
   GO TO 690
   640  newq = nq
   newh = 0
   eta = one
   hnew = h
   GO TO 690
   !-----------------------------------------------------------------------
   ! All returns are made through this section.
   ! On a successful return, ETAMAX is reset and ACOR is scaled.
   !-----------------------------------------------------------------------
   660  kflag = -1
   GO TO 720
   670  kflag = -2
   GO TO 720
   680  IF (nflag == -2) kflag = -3
   IF (nflag == -3) kflag = -4
   GO TO 720
   690  etamax = etamx3
   IF (nst <= 10) etamax = etamx2
   r = one/tq(2)
   !CALL dscal (n, r, acor, 1)
   do i=1,n
     acor(i)=acor(i)*r
   enddo
   720  jstart = 1
   RETURN
   !----------------------- End of Subroutine DVSTEP ----------------------
 END SUBROUTINE dvstep
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! DVSET                                                                        !
! Call sequence communication.. None                                           !
! COMMON block variables accessed..                                            !
!     /DVOD01/ -- EL(13), H, TAU(13), TQ(5), L(= NQ + 1),                      !
!                 METH, NQ, NQWAIT                                             !
                                                                               !
! Subroutines called by DVSET.. None                                           !
! Function routines called by DVSET.. None                                     !
!-----------------------------------------------------------------------       !
! DVSET is called by DVSTEP and sets coefficients for use there.               !
                                                                               !
! For each order NQ, the coefficients in EL are calculated by use of           !
!  the generating polynomial lambda(x), with coefficients EL(i).               !
!      lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).                   !
! For the backward differentiation formulas,                                   !
!                                     NQ-1                                     !
!      lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) .                  !
!                                     i = 1                                    !
! For the Adams formulas,                                                      !
!                              NQ-1                                            !
!      (d/dx) lambda(x) = c * product (1 + x/xi(i) ) ,                         !
!                              i = 1                                           !
!      lambda(-1) = 0,    lambda(0) = 1,                                       !
! where c is a normalization constant.                                         !
! In both cases, xi(i) is defined by                                           !
!      H*xi(i) = t sub n  -  t sub (n-i)                                       !
!              = H + TAU(1) + TAU(2) + ... TAU(i-1).                           !
                                                                               !
                                                                               !
! In addition to variables described previously, communication                 !
! with DVSET uses the following..                                              !
!   TAU    = A vector of length 13 containing the past NQ values               !
!            of H.                                                             !
!   EL     = A vector of length 13 in which vset stores the                    !
!            coefficients for the corrector formula.                           !
!   TQ     = A vector of length 5 in which vset stores constants               !
!            used for the convergence test, the error test, and the            !
!            selection of H at a new order.                                    !
!   METH   = The basic method indicator.                                       !
!   NQ     = The current order.                                                !
!   L      = NQ + 1, the length of the vector stored in EL, and                !
!            the number of columns of the YH array being used.                 !
!   NQWAIT = A counter controlling the frequency of order changes.             !
!            An order change is about to be considered if NQWAIT = 1.          !
!------------------------------------------------------------------------------!
 SUBROUTINE dvset()

   !-----------------------------------------------------------------------
   ! Type declarations for labeled COMMON block DVOD01 --------------------

   REAL(SEDP) :: acnrm, ccmxj, conp, crate, drc, el,                           &
                 eta, etamax, h, hmin, hmxi, hnew, hscal, prl1, rc, rl1,       &
                 tau, tq, tn, uround
   INTEGER  :: icf, init, ipup, jcur, jstart, jsv, kflag, kuth,              &
                 l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                  &
                 locjs, maxord, meth, miter, msbj, mxhnil, mxstep,             &
                 n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh

   ! Type declarations for local variables --------------------------------

   REAL(SEDP) :: ahatn0, alph0, cnqm1, cortes, csum, elp, em,  &
       em0, floti, flotl, flotnq, hsum, one, rxi, rxis, s, six,  &
       t1, t2, t3, t4, t5, t6, two, xi, zero
   INTEGER  :: i, iback, j, jp1, nqm1, nqm2

   DIMENSION em(13)
   !-----------------------------------------------------------------------
   ! The following Fortran-77 declaration is to cause the values of the
   ! listed (local) variables to be saved between calls to this integrator.
   !-----------------------------------------------------------------------
   SAVE cortes, one, six, two, zero

   COMMON /dvod01/ acnrm, ccmxj, conp, crate, drc, el(13),                     &
                   eta, etamax, h, hmin, hmxi, hnew, hscal, prl1,              &
                   rc, rl1, tau(13), tq(5), tn, uround,                        &
                   icf, init, ipup, jcur, jstart, jsv, kflag, kuth,            &
                   l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                &
                   locjs, maxord, meth, miter, msbj, mxhnil, mxstep,           &
                   n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh

   DATA cortes /0.1D+00/
   DATA one  /1.0D+00/, six /6.0D+00/, two /2.0D+00/, zero /0.0D+00/

   flotl = FLOAT(l)
   nqm1 = nq - 1
   nqm2 = nq - 2
   SELECT CASE ( meth )
     CASE (    1)
       GO TO 100
     CASE (    2)
       GO TO  200
   END SELECT

   ! Set coefficients for Adams methods. ----------------------------------
   100  IF (nq /= 1) GO TO 110
   el(1) = one
   el(2) = one
   tq(1) = one
   tq(2) = two
   tq(3) = six*tq(2)
   tq(5) = one
   GO TO 300
   110  hsum = h
   em(1) = one
   flotnq = flotl - one
   DO  i = 2, l
     em(i) = zero
   END DO
   DO  j = 1, nqm1
     IF ((j /= nqm1) .OR. (nqwait /= 1)) GO TO 130
     s = one
     csum = zero
     DO  i = 1, nqm1
       csum = csum + s*em(i)/FLOAT(i+1)
       s = -s
     END DO
     tq(1) = em(nqm1)/(flotnq*csum)
     130    rxi = h/hsum
     DO  iback = 1, j
       i = (j + 2) - iback
       em(i) = em(i) + em(i-1)*rxi
     END DO
     hsum = hsum + tau(j)
   END DO
   ! Compute integral from -1 to 0 of polynomial and of x times it. -------
   s = one
   em0 = zero
   csum = zero
   DO  i = 1, nq
     floti = FLOAT(i)
     em0 = em0 + s*em(i)/floti
     csum = csum + s*em(i)/(floti+one)
     s = -s
   END DO
   ! In EL, form coefficients of normalized integrated polynomial. --------
   s = one/em0
   el(1) = one
   DO  i = 1, nq
     el(i+1) = s*em(i)/FLOAT(i)
   END DO
   xi = hsum/h
   tq(2) = xi*em0/csum
   tq(5) = xi/el(l)
   IF (nqwait /= 1) GO TO 300
   ! For higher order control constant, multiply polynomial by 1+x/xi(q). -
   rxi = one/xi
   DO  iback = 1, nq
     i = (l + 1) - iback
     em(i) = em(i) + em(i-1)*rxi
   END DO
   ! Compute integral of polynomial. --------------------------------------
   s = one
   csum = zero
   DO  i = 1, l
     csum = csum + s*em(i)/FLOAT(i+1)
     s = -s
   END DO
   tq(3) = flotl*em0/csum
   GO TO 300

   ! Set coefficients for BDF methods. ------------------------------------
   200 DO  i = 3, l
     el(i) = zero
   END DO
   el(1) = one
   el(2) = one
   alph0 = -one
   ahatn0 = -one
   hsum = h
   rxi = one
   rxis = one
   IF (nq == 1) GO TO 240
   DO  j = 1, nqm2
   ! In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
     hsum = hsum + tau(j)
     rxi = h/hsum
     jp1 = j + 1
     alph0 = alph0 - one/FLOAT(jp1)
     DO  iback = 1, jp1
       i = (j + 3) - iback
       el(i) = el(i) + el(i-1)*rxi
     END DO
   END DO
   alph0 = alph0 - one/FLOAT(nq)
   rxis = -el(2) - alph0
   hsum = hsum + tau(nqm1)
   rxi = h/hsum
   ahatn0 = -el(2) - rxi
   DO  iback = 1, nq
     i = (nq + 2) - iback
     el(i) = el(i) + el(i-1)*rxis
   END DO
   240 t1 = one - ahatn0 + alph0
   t2 = one + FLOAT(nq)*t1
   tq(2) = ABS(alph0*t2/t1)
   tq(5) = ABS(t2/(el(l)*rxi/rxis))
   IF (nqwait /= 1) GO TO 300
   cnqm1 = rxis/el(l)
   t3 = alph0 + one/FLOAT(nq)
   t4 = ahatn0 + rxi
   elp = t3/(one - t4 + t3)
   tq(1) = ABS(elp/cnqm1)
   hsum = hsum + tau(nq)
   rxi = h/hsum
   t5 = alph0 - one/FLOAT(nq+1)
   t6 = ahatn0 - rxi
   elp = t2/(one - t6 + t5)
   tq(3) = ABS(elp*rxi*(flotl + one)*t5)
   300 tq(4) = cortes*tq(2)
   RETURN
   !----------------------- End of Subroutine DVSET -----------------------
 END SUBROUTINE dvset
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! dvnlsd
! Call sequence input -- Y, YH, LDYH, SAVF, EWT, ACOR, IWM, WM,
!                        F, JAC, NFLAG, RPAR, IPAR
! Call sequence output -- YH, ACOR, WM, IWM, NFLAG
! COMMON block variables accessed..
!     /DVOD01/ ACNRM, CRATE, DRC, H, RC, RL1, TQ(5), TN, ICF,
!                JCUR, METH, MITER, N, NSLP
!     /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST

! Subroutines called by DVNLSD.. F, DAXPY, DCOPY, DSCAL, DVJAC, DVSOL
! Function routines called by DVNLSD.. DVNORM
!-----------------------------------------------------------------------
! Subroutine DVNLSD is a nonlinear system solver, which uses functional
! iteration or a chord (modified Newton) method.  For the chord method
! direct linear algebraic system solvers are used.  Subroutine DVNLSD
! then handles the corrector phase of this integration package.

! Communication with DVNLSD is done with the following variables. (For
! more details, please see the comments in the driver subroutine.)

! Y          = The dependent variable, a vector of length N, input.
! YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input
!              and output.  On input, it contains predicted values.
! LDYH       = A constant .ge. N, the first dimension of YH, input.
! VSAV       = Unused work array.
! SAVF       = A work array of length N.
! EWT        = An error weight vector of length N, input.
! ACOR       = A work array of length N, used for the accumulated
!              corrections to the predicted y vector.
! WM,IWM     = FLOAT and integer work arrays associated with matrix
!              operations in chord iteration (MITER .ne. 0).
! F          = Dummy name for user supplied routine for f.
! JAC        = Dummy name for user supplied Jacobian routine.
! PDUM       = Unused dummy subroutine name.  Included for uniformity
!              over collection of integrators.
! NFLAG      = Input/output flag, with values and meanings as follows..
!              INPUT
!                  0 first CALL for this time step.
!                 -1 convergence failure in previous CALL to DVNLSD.
!                 -2 error test failure in DVSTEP.
!              OUTPUT
!                  0 successful completion of nonlinear solver.
!                 -1 convergence failure or singular matrix.
!                 -2 unrecoverable error in matrix preprocessing
!                    (cannot occur here).
!                 -3 unrecoverable error in solution (cannot occur
!                    here).
! RPAR, IPAR = Dummy names for users REAL(SEDP) and integer work arrays.

! IPUP       = Own variable flag with values and meanings as follows..
!              0,            do not update the Newton matrix.
!              MITER .ne. 0, update Newton matrix, because it is the
!                            initial step, order was changed, the error
!                            test failed, or an update is indicated by
!                            the scalar RC or step counter NST.

! For more details, see comments in driver subroutine.
!------------------------------------------------------------------------------!
 SUBROUTINE dvnlsd (candi_, y, yh, ldyh, vsav, savf, ewt, acor, iwm, wm,        &
                    f, jac, pdum, nflag, rpar, ipar)

   TYPE(aed2_sed_candi_t),TARGET,INTENT(inout) :: candi_
   INTEGER,    INTENT(INOUT)          :: ldyh
   REAL(SEDP), INTENT(OUT)            :: y(*)
   REAL(SEDP), INTENT(IN)             :: yh(ldyh,*)
   REAL(SEDP), INTENT(INOUT)          :: vsav(*)
   REAL(SEDP), INTENT(OUT)            :: savf(*)
   REAL(SEDP), INTENT(INOUT)          :: ewt(*)
   REAL(SEDP), INTENT(OUT)            :: acor(*)
   INTEGER,    INTENT(INOUT)          :: iwm(*)
   REAL(SEDP), INTENT(INOUT)          :: wm(*)
   INTEGER,    INTENT(INOUT)          :: nflag
   REAL(SEDP), INTENT(INOUT)          :: rpar(*)
   INTEGER,    INTENT(INOUT)          :: ipar(*)
   EXTERNAL f, jac, pdum

   !-----------------------------------------------------------------------

   ! Type declarations for labeled COMMON block DVOD01 --------------------

   REAL(SEDP) :: acnrm, ccmxj, conp, crate, drc, el,                           &
                 eta, etamax, h, hmin, hmxi, hnew, hscal, prl1, rc, rl1,       &
                 tau, tq, tn, uround
   INTEGER  :: icf, init, ipup, jcur, jstart, jsv, kflag, kuth,              &
                 l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                  &
                 locjs, maxord, meth, miter, msbj, mxhnil, mxstep,             &
                 n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh

   ! Type declarations for labeled COMMON block DVOD02 --------------------

   REAL(SEDP) :: hu
   INTEGER  :: ncfn, netf, nfe, nje, nlu, nni, nqu, nst

   ! Type declarations for local variables --------------------------------

   REAL(SEDP) :: ccmax, crdown, cscale, dcon, del, delp, one,                  &
                 rdiv, two, zero
   INTEGER  :: i, ierpj, iersl, m, maxcor, msbp
   TYPE(aed2_sed_candi_t),POINTER :: candi

   ! Type declaration for function subroutines called ---------------------

   !REAL(SEDP) :: dvnorm

   !-----------------------------------------------------------------------
   ! The following Fortran-77 declaration is to cause the values of the
   ! listed (local) variables to be saved between calls to this integrator.
   !-----------------------------------------------------------------------
   SAVE ccmax, crdown, maxcor, msbp, rdiv, one, two, zero

   COMMON /dvod01/ acnrm, ccmxj, conp, crate, drc, el(13),                     &
                   eta, etamax, h, hmin, hmxi, hnew, hscal, prl1,              &
                   rc, rl1, tau(13), tq(5), tn, uround,                        &
                   icf, init, ipup, jcur, jstart, jsv, kflag, kuth,            &
                   l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                &
                   locjs, maxord, meth, miter, msbj, mxhnil, mxstep,           &
                   n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh
   COMMON /dvod02/ hu, ncfn, netf, nfe, nje, nlu, nni, nqu, nst

   DATA ccmax /0.3D+00/, crdown /0.3D+00/, maxcor /3/, msbp /20/, rdiv  /2.0D+00/
   DATA one /1.0D+00/, two /2.0D+00/, zero /0.0D+00/

   candi => candi_

   !-----------------------------------------------------------------------
   ! On the first step, on a change of method order, or after a
   ! nonlinear convergence failure with NFLAG = -2, set IPUP = MITER
   ! to force a Jacobian update when MITER .ne. 0.
   !-----------------------------------------------------------------------
   IF (jstart == 0) nslp = 0
   IF (nflag == 0) icf = 0
   IF (nflag == -2) ipup = miter
   IF ( (jstart == 0) .OR. (jstart == -1) ) ipup = miter
   ! If this is functional iteration, set CRATE .eq. 1 and drop to 220
   IF (miter == 0) THEN
     crate = one
     GO TO 220
   END IF
   !-----------------------------------------------------------------------
   ! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
   ! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
   ! to force DVJAC to be called, if a Jacobian is involved.
   ! In any case, DVJAC is called at least every MSBP steps.
   !-----------------------------------------------------------------------
   drc = ABS(rc-one)
   IF (drc > ccmax .OR. nst >= nslp+msbp) ipup = miter
   !-----------------------------------------------------------------------
   ! Up to MAXCOR corrector iterations are taken.  A convergence test is
   ! made on the r.m.s. norm of each correction, weighted by the error
   ! weight vector EWT.  The sum of the corrections is accumulated in the
   ! vector ACOR(i).  The YH array is not altered in the corrector loop.
   !-----------------------------------------------------------------------
   220  m = 0
   delp = zero
   !CALL dcopy (n, yh(1,1), 1, y, 1 )
   do i=1,n
     y(i)=yh(i,1)
   enddo
   CALL f (candi, n, tn, y, savf, rpar, ipar)
   nfe = nfe + 1
   IF (ipup <= 0) GO TO 250
   !-----------------------------------------------------------------------
   ! If indicated, the matrix P = I - h*rl1*J is reevaluated and
   ! preprocessed before starting the corrector iteration.  IPUP is set
   ! to 0 as an indicator that this has been done.
   !-----------------------------------------------------------------------
   CALL dvjac (candi, y, yh, ldyh, ewt, acor, savf, wm, iwm, f, jac, ierpj, rpar, ipar)
   ipup = 0
   rc = one
   drc = zero
   crate = one
   nslp = nst
   ! If matrix is singular, take error return to force cut in step size. --
   IF (ierpj /= 0) GO TO 430
   250 DO  i = 1,n
     acor(i) = zero
   END DO
   ! This is a looping point for the corrector iteration. -----------------
   270  IF (miter /= 0) GO TO 350
   !-----------------------------------------------------------------------
   ! In the case of functional iteration, update Y directly from
   ! the result of the last function evaluation.
   !-----------------------------------------------------------------------
   DO  i = 1,n
     savf(i) = rl1*(h*savf(i) - yh(i,2))
   END DO
   DO  i = 1,n
     y(i) = savf(i) - acor(i)
   END DO
   del = dvnorm (n, y, ewt)
   DO  i = 1,n
     y(i) = yh(i,1) + savf(i)
   END DO
   !CALL dcopy (n, savf, 1, acor, 1)
   do i=1,n
     acor(i)=savf(i)
   enddo
   GO TO 400
   !-----------------------------------------------------------------------
   ! In the case of the chord method, compute the corrector error,
   ! and solve the linear system with that as right-hand side and
   ! P as coefficient matrix.  The correction is scaled by the factor
   ! 2/(1+RC) to account for changes in h*rl1 since the last DVJAC CALL.
   !-----------------------------------------------------------------------
   350  DO  i = 1,n
     y(i) = (rl1*h)*savf(i) - (rl1*yh(i,2) + acor(i))
   END DO
   CALL dvsol (wm, iwm, y, iersl)
   nni = nni + 1
   IF (iersl > 0) GO TO 410
   IF (meth == 2 .AND. rc /= one) THEN
     cscale = two/(one + rc)
   !  CALL dscal (n, cscale, y, 1)
     do i=1,n
       y(i)=y(i)*cscale
     enddo
   END IF
   del = dvnorm (n, y, ewt)
   !CALL daxpy (n, one, y, acor)
   !print *, "ROGER1" ,n
   DO  i = 1,n
     acor(i)=acor(i)+y(i)
   END DO

   DO  i = 1,n
     y(i) = yh(i,1) + acor(i)
   END DO
   !-----------------------------------------------------------------------
   ! Test for convergence.  If M .gt. 0, an estimate of the convergence
   ! rate constant is stored in CRATE, and this is used in the test.
   !-----------------------------------------------------------------------
   400  IF (m /= 0) crate = MAX(crdown*crate,del/delp)
   dcon = del*MIN(one,crate)/tq(4)
   IF (dcon <= one) GO TO 450
   m = m + 1
   IF (m == maxcor) GO TO 410
   IF (m >= 2 .AND. del > rdiv*delp) GO TO 410
   delp = del
   CALL f (candi, n, tn, y, savf, rpar, ipar)
   nfe = nfe + 1
   GO TO 270

   410  IF (miter == 0 .OR. jcur == 1) GO TO 430
   icf = 1
   ipup = miter
   GO TO 220

   430  CONTINUE
   nflag = -1
   icf = 2
   ipup = miter
   RETURN

   ! Return for successful step. ------------------------------------------
   450  nflag = 0
   jcur = 0
   icf = 0
   IF (m == 0) acnrm = del
   IF (m > 0) acnrm = dvnorm (n, acor, ewt)
   RETURN
   !----------------------- End of Subroutine DVNLSD ----------------------
 END SUBROUTINE dvnlsd
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! DVJAC                                                                        !
! Call sequence input -- Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM,                !
!                        F, JAC, RPAR, IPAR                                    !
! Call sequence output -- WM, IWM, IERPJ                                       !
! COMMON block variables accessed..                                            !
!     /DVOD01/  CCMXJ, DRC, H, RL1, TN, UROUND, ICF, JCUR, LOCJS,              !
!               MSBJ, NSLJ                                                     !
!     /DVOD02/  NFE, NST, NJE, NLU                                             !
                                                                               !
! Subroutines called by DVJAC.. F, JAC, DACOPY, DCOPY, DGBFA, DGEFA,           !
!                              DSCAL                                           !
! Function routines called by DVJAC.. DVNORM                                   !
!-----------------------------------------------------------------------       !
! DVJAC is called by DVSTEP to compute and process the matrix                  !
! P = I - h*rl1*J , where J is an approximation to the Jacobian.               !
! Here J is computed by the user-supplied routine JAC if                       !
! MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.             !
! If MITER = 3, a diagonal approximation to J is used.                         !
! If JSV = -1, J is computed from scratch in all cases.                        !
! If JSV = 1 and MITER = 1, 2, 4, or 5, and if the saved value of J is         !
! considered acceptable, then P is constructed from the saved J.               !
! J is stored in wm and replaced by P.  If MITER .ne. 3, P is then             !
! subjected to LU decomposition in preparation for later solution              !
! of linear systems with P as coefficient matrix. This is done                 !
! by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.                  !
                                                                               !
! Communication with DVJAC is done with the following variables.  (For         !
! more details, please see the comments in the driver subroutine.)             !
! Y          = Vector containing predicted values on entry.                    !
! YH         = The Nordsieck array, an LDYH by LMAX array, input.              !
! LDYH       = A constant .ge. N, the first dimension of YH, input.            !
! EWT        = An error weight vector of length N.                             !
! SAVF       = Array containing f evaluated at predicted y, input.             !
! WM         = FLOAT work space for matrices.  In the output, it containS      !
!              the inverse diagonal matrix if MITER = 3 and the LU             !
!              decomposition of P if MITER is 1, 2 , 4, or 5.                  !
!              Storage of matrix elements starts at WM(3).                     !
!              Storage of the saved Jacobian starts at WM(LOCJS).              !
!              WM also contains the following matrix-related data..            !
!              WM(1) = SQRT(UROUND), used in numerical Jacobian step.          !
!              WM(2) = H*RL1, saved for later use if MITER = 3.                !
! IWM        = Integer work space containing pivot information,                !
!              starting at IWM(31), if MITER is 1, 2, 4, or 5.                 !
!              IWM also contains band parameters ML = IWM(1) and               !
!              MU = IWM(2) if MITER is 4 or 5.                                 !
! F          = Dummy name for the user supplied subroutine for f.              !
! JAC        = Dummy name for the user supplied Jacobian subroutine.           !
! RPAR, IPAR = Dummy names for users REAL(SEDP) and integer work arrays.       !
! RL1        = 1/EL(2) (input).                                                !
! IERPJ      = Output error flag,  = 0 if no trouble, 1 if the P               !
!              matrix is found to be singular.                                 !
! JCUR       = Output flag to indicate whether the Jacobian matrix             !
!              (or approximation) is now current.                              !
!              JCUR = 0 means J is not current.                                !
!              JCUR = 1 means J is current.                                    !
!------------------------------------------------------------------------------!
 SUBROUTINE dvjac (candi, y, yh, ldyh, ewt, ftem, savf, wm, iwm, f, jac,       &
                   ierpj, rpar, ipar)

   TYPE(aed2_sed_candi_t),POINTER,INTENT(inout) :: candi
   REAL(SEDP), INTENT(INOUT)         :: y(*)
   INTEGER,    INTENT(IN)            :: ldyh
   REAL(SEDP), INTENT(IN)            :: yh(ldyh,*)
   REAL(SEDP), INTENT(IN)            :: ewt(*)
   REAL(SEDP), INTENT(INOUT)         :: ftem(*)
   REAL(SEDP), INTENT(IN)            :: savf(*)
   REAL(SEDP), INTENT(OUT)           :: wm(*)
   INTEGER,    INTENT(OUT)           :: iwm(*)
   INTEGER,    INTENT(OUT)           :: ierpj
   REAL(SEDP), INTENT(INOUT)         :: rpar(*)
   INTEGER,    INTENT(INOUT)         :: ipar(*)
   EXTERNAL f, jac


   !-----------------------------------------------------------------------
   ! Type declarations for labeled COMMON block DVOD01 --------------------

   REAL(SEDP) :: acnrm, ccmxj, conp, crate, drc, el,                         &
                 eta, etamax, h, hmin, hmxi, hnew, hscal, prl1, rc, rl1,     &
                 tau, tq, tn, uround
   INTEGER  :: icf, init, ipup, jcur, jstart, jsv, kflag, kuth,              &
                 l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                &
                 locjs, maxord, meth, miter, msbj, mxhnil, mxstep,           &
                 n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh

   ! Type declarations for labeled COMMON block DVOD02 --------------------

   REAL(SEDP) :: hu
   INTEGER  :: ncfn, netf, nfe, nje, nlu, nni, nqu, nst

   ! Type declarations for local variables --------------------------------

   REAL(SEDP) :: con, di, fac, hrl1, one, pt1, r, r0, srur, thou,            &
                 yi, yj, yjj, zero
   INTEGER  :: i, i1, i2, ier, ii, j, j1, jj, jok, lenp, mba, mband,         &
                 meb1, meband, ml, ml3, mu, np1

   ! Type declaration for function subroutines called ---------------------

   !REAL(SEDP) :: dvnorm

   !-----------------------------------------------------------------------
   ! The following Fortran-77 declaration is to cause the values of the
   ! listed (local) variables to be saved between calls to this subroutine.
   !-----------------------------------------------------------------------
   SAVE one, pt1, thou, zero
   !-----------------------------------------------------------------------
   COMMON /dvod01/ acnrm, ccmxj, conp, crate, drc, el(13),                     &
                   eta, etamax, h, hmin, hmxi, hnew, hscal, prl1,              &
                   rc, rl1, tau(13), tq(5), tn, uround,                        &
                   icf, init, ipup, jcur, jstart, jsv, kflag, kuth,            &
                   l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                &
                   locjs, maxord, meth, miter, msbj, mxhnil, mxstep,           &
                   n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh
   COMMON /dvod02/ hu, ncfn, netf, nfe, nje, nlu, nni, nqu, nst

   DATA one /1.0D+00/,thou /1000.0D+00/,zero /0.0D+00/,pt1 /0.1D+00/

   ierpj = 0
   hrl1 = h*rl1
   ! See whether J should be evaluated (JOK = -1) or not (JOK = 1). -------
   jok = jsv
   IF (jsv == 1) THEN
     IF (nst == 0 .OR. nst > nslj+msbj) jok = -1
     IF (icf == 1 .AND. drc < ccmxj) jok = -1
     IF (icf == 2) jok = -1
   END IF
   ! End of setting JOK. --------------------------------------------------

   IF (jok == -1 .AND. miter == 1) THEN
   ! If JOK = -1 and MITER = 1, CALL JAC to evaluate Jacobian. ------------
     nje = nje + 1
     nslj = nst
     jcur = 1
     lenp = n*n
     DO  i = 1,lenp
       wm(i+2) = zero
     END DO
     CALL jac (candi, n, tn, y, 0, 0, wm(3), n, rpar, ipar)
     IF (jsv == 1) CALL dcopy (lenp, wm(3), 1, wm(locjs),1)
   END IF

   IF (jok == -1 .AND. miter == 2) THEN
   ! If MITER = 2, make N calls to F to approximate the Jacobian. ---------
     nje = nje + 1
     nslj = nst
     jcur = 1
     fac = dvnorm (n, savf, ewt)
     r0 = thou*ABS(h)*uround*FLOAT(n)*fac
     IF (r0 == zero) r0 = one
     srur = wm(1)
     j1 = 2
     DO  j = 1,n
       yj = y(j)
       r = MAX(srur*ABS(yj),r0/ewt(j))
       y(j) = y(j) + r
       fac = one/r
       CALL f (candi, n, tn, y, ftem, rpar, ipar)
       DO  i = 1,n
         wm(i+j1) = (ftem(i) - savf(i))*fac
       END DO
       y(j) = yj
       j1 = j1 + n
     END DO
     nfe = nfe + n
     lenp = n*n
     IF (jsv == 1) CALL dcopy (lenp, wm(3), 1, wm(locjs),1)
   END IF

   IF (jok == 1 .AND. (miter == 1 .OR. miter == 2)) THEN
     jcur = 0
     lenp = n*n
     CALL dcopy (lenp, wm(locjs), 1, wm(3), 1)
   END IF

   IF (miter == 1 .OR. miter == 2) THEN
   ! Multiply Jacobian by scalar, add identity, and do LU decomposition. --
     con = -hrl1
     CALL dscal (lenp, con, wm(3), 1)
     j = 3
     np1 = n + 1
     DO  i = 1,n
       wm(j) = wm(j) + one
       j = j + np1
     END DO
     nlu = nlu + 1
     CALL dgefa (wm(3), n, n, iwm(31), ier)
     IF (ier /= 0) ierpj = 1
     RETURN
   END IF
   ! End of code block for MITER = 1 or 2. --------------------------------

   IF (miter == 3) THEN
   ! If MITER = 3, construct a diagonal approximation to J and P. ---------
     nje = nje + 1
     jcur = 1
     wm(2) = hrl1
     r = rl1*pt1
     DO  i = 1,n
       y(i) = y(i) + r*(h*savf(i) - yh(i,2))
     END DO
     CALL f (candi,n, tn, y, wm(3), rpar, ipar)
     nfe = nfe + 1
     DO  i = 1,n
       r0 = h*savf(i) - yh(i,2)
       di = pt1*r0 - h*(wm(i+2) - savf(i))
       wm(i+2) = one
       IF (ABS(r0) < uround/ewt(i)) CYCLE
       IF (ABS(di) == zero) GO TO 330
       wm(i+2) = pt1*r0/di
     END DO
     RETURN
     330  ierpj = 1
     RETURN
   END IF
   ! End of code block for MITER = 3. -------------------------------------

   ! Set constants for MITER = 4 or 5. ------------------------------------
   ml = iwm(1)
   mu = iwm(2)
   ml3 = ml + 3
   mband = ml + mu + 1
   meband = mband + ml
   lenp = meband*n

   IF (jok == -1 .AND. miter == 4) THEN
   ! If JOK = -1 and MITER = 4, CALL JAC to evaluate Jacobian. ------------
     nje = nje + 1
     nslj = nst
     jcur = 1
     DO  i = 1,lenp
       wm(i+2) = zero
     END DO
     CALL jac (candi, n, tn, y, ml, mu, wm(ml3), meband, rpar, ipar)
     IF (jsv == 1) CALL dacopy (mband, n, wm(ml3), meband, wm(locjs), mband)
   END IF

   IF (jok == -1 .AND. miter == 5) THEN
   ! If MITER = 5, make N calls to F to approximate the Jacobian. ---------
     nje = nje + 1
     nslj = nst
     jcur = 1
     mba = MIN(mband,n)
     meb1 = meband - 1
     srur = wm(1)
     fac = dvnorm (n, savf, ewt)
     r0 = thou*ABS(h)*uround*FLOAT(n)*fac
     IF (r0 == zero) r0 = one
     DO  j = 1,mba
       DO  i = j,n,mband
         yi = y(i)
         r = MAX(srur*ABS(yi),r0/ewt(i))
         y(i) = y(i) + r
       END DO
       CALL f (candi,n, tn, y, ftem, rpar, ipar)
       DO  jj = j,n,mband
         y(jj) = yh(jj,1)
         yjj = y(jj)
         r = MAX(srur*ABS(yjj),r0/ewt(jj))
         fac = one/r
         i1 = MAX(jj-mu,1)
         i2 = MIN(jj+ml,n)
         ii = jj*meb1 - ml + 2
         DO  i = i1,i2
           wm(ii+i) = (ftem(i) - savf(i))*fac
         END DO
       END DO
     END DO
     nfe = nfe + mba
     IF (jsv == 1) THEN
       CALL dacopy (mband, n, wm(ml3), meband, wm(locjs), mband)
     END IF
   END IF

   IF (jok == 1) THEN
     jcur = 0
     CALL dacopy (mband, n, wm(locjs), mband, wm(ml3), meband)
   END IF

   ! Multiply Jacobian by scalar, add identity, and do LU decomposition.
   con = -hrl1
   CALL dscal (lenp, con, wm(3), 1 )
   ii = mband + 2
   DO  i = 1,n
     wm(ii) = wm(ii) + one
     ii = ii + meband
   END DO
   nlu = nlu + 1
   CALL dgbfa (wm(3), meband, n, ml, mu, iwm(31), ier)
   IF (ier /= 0) ierpj = 1
   RETURN
   ! End of code block for MITER = 4 or 5. --------------------------------

   !----------------------- End of Subroutine DVJAC -----------------------
 END SUBROUTINE dvjac
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! DVSOL                                                                        !
! Call sequence input -- WM, IWM, X                                            !
! Call sequence output -- X, IERSL                                             !
! COMMON block variables accessed..                                            !
!     /DVOD01/ -- H, RL1, MITER, N                                             !
                                                                               !
! Subroutines called by DVSOL.. DGESL, DGBSL                                   !
! Function routines called by DVSOL.. None                                     !
!-----------------------------------------------------------------------       !
! This routine manages the solution of the linear system arising from          !
! a chord iteration.  It is called if MITER .ne. 0.                            !
! If MITER is 1 or 2, it calls DGESL to accomplish this.                       !
! If MITER = 3 it updates the coefficient H*RL1 in the diagonal                !
! matrix, and then computes the solution.                                      !
! If MITER is 4 or 5, it calls DGBSL.                                          !
! Communication with DVSOL uses the following variables..                      !
! WM    = FLOAT work space containing the inverse diagonal matrix if           !
!         MITER = 3 and the LU decomposition of the matrix otherwise.          !
!         Storage of matrix elements starts at WM(3).                          !
!         WM also contains the following matrix-related data..                 !
!         WM(1) = SQRT(UROUND) (not used here),                                !
!         WM(2) = HRL1, the previous value of H*RL1, used if MITER = 3.        !
! IWM   = Integer work space containing pivot information, starting at         !
!         IWM(31), if MITER is 1, 2, 4, or 5.  IWM also contains band          !
!         parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.           !
! X     = The right-hand side vector on input, and the solution vector         !
!         on output, of length N.                                              !
! IERSL = Output flag.  IERSL = 0 if no trouble occurred.                      !
!         IERSL = 1 if a singular matrix arose with MITER = 3.                 !
!------------------------------------------------------------------------------!
 SUBROUTINE dvsol (wm, iwm, x, iersl)

   REAL(SEDP), INTENT(INOUT)          :: wm(*)
   INTEGER,    INTENT(IN)             :: iwm(*)
   REAL(SEDP), INTENT(OUT)            :: x(*)
   INTEGER,    INTENT(OUT)            :: iersl

   ! Type declarations for labeled COMMON block DVOD01 --------------------

   REAL(SEDP) :: acnrm, ccmxj, conp, crate, drc, el,                           &
                 eta, etamax, h, hmin, hmxi, hnew, hscal, prl1, rc, rl1,       &
                 tau, tq, tn, uround
   INTEGER  :: icf, init, ipup, jcur, jstart, jsv, kflag, kuth,              &
                 l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                  &
                 locjs, maxord, meth, miter, msbj, mxhnil, mxstep,             &
                 n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh

   ! Type declarations for local variables --------------------------------

   INTEGER  :: i, meband, ml, mu
   REAL(SEDP) :: di, hrl1, one, phrl1, r, zero

   !-----------------------------------------------------------------------
   ! The following Fortran-77 declaration is to cause the values of the
   ! listed (local) variables to be saved between calls to this integrator.
   !-----------------------------------------------------------------------
   SAVE one, zero

   COMMON /dvod01/ acnrm, ccmxj, conp, crate, drc, el(13),                     &
                   eta, etamax, h, hmin, hmxi, hnew, hscal, prl1,              &
                   rc, rl1, tau(13), tq(5), tn, uround,                        &
                   icf, init, ipup, jcur, jstart, jsv, kflag, kuth,            &
                   l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                &
                   locjs, maxord, meth, miter, msbj, mxhnil, mxstep,           &
                   n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh

   DATA one /1.0D+00/, zero /0.0D+00/

   iersl = 0
   SELECT CASE ( miter )
     CASE (    1)
       GO TO 100
     CASE (    2)
       GO TO  100
     CASE (    3)
       GO TO  300
     CASE (    4)
       GO TO  400
     CASE (    5)
       GO TO  400
   END SELECT
   100  CALL dgesl (wm(3), n, n, iwm(31), x, 0)
   RETURN

   300  phrl1 = wm(2)
   hrl1 = h*rl1
   wm(2) = hrl1
   IF (hrl1 == phrl1) GO TO 330
   r = hrl1/phrl1
   DO  i = 1,n
     di = one - r*(one - one/wm(i+2))
     IF (ABS(di) == zero) GO TO 390
     wm(i+2) = one/di
   END DO

   330  DO  i = 1,n
     x(i) = wm(i+2)*x(i)
   END DO
   RETURN
   390 iersl = 1
   RETURN

   400  ml = iwm(1)
   mu = iwm(2)
   meband = 2*ml + mu + 1
   CALL dgbsl (wm(3), meband, n, ml, mu, iwm(31), x, 0)
   RETURN
   !----------------------- End of Subroutine DVSOL -----------------------
 END SUBROUTINE dvsol
!------------------------------------------------------------------------------!



!------------------------------------------------------------------------------!
! DVJUST                                                                       !
! Call sequence input -- YH, LDYH, IORD                                        !
! Call sequence output -- YH                                                   !
! COMMON block input -- NQ, METH, LMAX, HSCAL, TAU(13), N                      !
! COMMON block variables accessed..                                            !
!     /DVOD01/ -- HSCAL, TAU(13), LMAX, METH, N, NQ,                           !
                                                                               !
! Subroutines called by DVJUST.. DAXPY                                         !
! Function routines called by DVJUST.. None                                    !
!-----------------------------------------------------------------------       !
! This subroutine adjusts the YH array on reduction of order,                  !
! and also when the order is increased for the stiff option (METH = 2).        !
! Communication with DVJUST uses the following..                               !
! IORD  = An integer flag used when METH = 2 to indicate an order              !
!         increase (IORD = +1) or an order decrease (IORD = -1).               !
! HSCAL = Step size H used in scaling of Nordsieck array YH.                   !
!         (If IORD = +1, DVJUST assumes that HSCAL = TAU(1).)                  !
! See References 1 and 2 for details.                                          !
!------------------------------------------------------------------------------!
 SUBROUTINE dvjust (yh, ldyh, iord)

   INTEGER,    INTENT(IN)             :: ldyh
   INTEGER,    INTENT(IN)             :: iord
   REAL(SEDP), INTENT(OUT)            :: yh(ldyh,*)

   ! Type declarations for labeled COMMON block DVOD01 --------------------

   REAL(SEDP) :: acnrm, ccmxj, conp, crate, drc, el,                           &
                 eta, etamax, h, hmin, hmxi, hnew, hscal, prl1, rc, rl1,       &
                 tau, tq, tn, uround
   INTEGER  :: icf, init, ipup, jcur, jstart, jsv, kflag, kuth,                 &
              l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                     &
              locjs, maxord, meth, miter, msbj, mxhnil, mxstep,                &
              n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh

   ! Type declarations for local variables --------------------------------

   REAL(SEDP) :: alph0, alph1, hsum, one, prod, t1, xi,xiold, zero
   INTEGER  :: i, iback, j, jp1, lp1, nqm1, nqm2, nqp1

   !-----------------------------------------------------------------------
   ! The following Fortran-77 declaration is to cause the values of the
   ! listed (local) variables to be saved between calls to this integrator.
   !-----------------------------------------------------------------------
   SAVE one, zero

   COMMON /dvod01/ acnrm, ccmxj, conp, crate, drc, el(13),                     &
                   eta, etamax, h, hmin, hmxi, hnew, hscal, prl1,              &
                   rc, rl1, tau(13), tq(5), tn, uround,                        &
                   icf, init, ipup, jcur, jstart, jsv, kflag, kuth,            &
                   l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,                &
                   locjs, maxord, meth, miter, msbj, mxhnil, mxstep,           &
                   n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh

   DATA one /1.0D+00/, zero /0.0D+00/

   IF ((nq == 2) .AND. (iord /= 1)) RETURN
   nqm1 = nq - 1
   nqm2 = nq - 2
   SELECT CASE ( meth )
     CASE (    1)
       GO TO 100
     CASE (    2)
       GO TO  200
   END SELECT
   !-----------------------------------------------------------------------
   ! Nonstiff option...
   ! Check to see if the order is being increased or decreased.
   !-----------------------------------------------------------------------
   100  CONTINUE
   IF (iord == 1) GO TO 180
   ! Order decrease. ------------------------------------------------------
   DO  j = 1, lmax
     el(j) = zero
   END DO
   el(2) = one
   hsum = zero
   DO  j = 1, nqm2
   ! Construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
     hsum = hsum + tau(j)
     xi = hsum/hscal
     jp1 = j + 1
     DO  iback = 1, jp1
       i = (j + 3) - iback
       el(i) = el(i)*xi + el(i-1)
     END DO
   END DO
   ! Construct coefficients of integrated polynomial. ---------------------
   DO  j = 2, nqm1
     el(j+1) = FLOAT(nq)*el(j)/FLOAT(j)
   END DO
   ! Subtract correction terms from YH array. -----------------------------
   DO  j = 3, nq
     DO  i = 1, n
       yh(i,j) = yh(i,j) - yh(i,l)*el(j)
     END DO
   END DO
   RETURN
   ! Order increase. ------------------------------------------------------
   ! Zero out next column in YH array. ------------------------------------
   180  CONTINUE
   lp1 = l + 1
   DO  i = 1, n
     yh(i,lp1) = zero
   END DO
   RETURN
   !-----------------------------------------------------------------------
   ! Stiff option...
   ! Check to see if the order is being increased or decreased.
   !-----------------------------------------------------------------------
   200  CONTINUE
   IF (iord == 1) GO TO 300
   ! Order decrease. ------------------------------------------------------
   DO  j = 1, lmax
     el(j) = zero
   END DO
   el(3) = one
   hsum = zero
   DO  j = 1,nqm2
   ! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
     hsum = hsum + tau(j)
     xi = hsum/hscal
     jp1 = j + 1
     DO  iback = 1, jp1
       i = (j + 4) - iback
       el(i) = el(i)*xi + el(i-1)
     END DO
   END DO
   ! Subtract correction terms from YH array. -----------------------------
   DO  j = 3,nq
     DO  i = 1, n
       yh(i,j) = yh(i,j) - yh(i,l)*el(j)
     END DO
   END DO
   RETURN
   ! Order increase. ------------------------------------------------------
   300  DO  j = 1, lmax
     el(j) = zero
   END DO
   el(3) = one
   alph0 = -one
   alph1 = one
   prod = one
   xiold = one
   hsum = hscal
   IF (nq == 1) GO TO 340
   DO  j = 1, nqm1
   ! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
     jp1 = j + 1
     hsum = hsum + tau(jp1)
     xi = hsum/hscal
     prod = prod*xi
     alph0 = alph0 - one/FLOAT(jp1)
     alph1 = alph1 + one/xi
     DO  iback = 1, jp1
       i = (j + 4) - iback
       el(i) = el(i)*xiold + el(i-1)
     END DO
     xiold = xi
   END DO
   340 CONTINUE
   t1 = (-alph0 - alph1)/prod
   ! Load column L + 1 in YH array. ---------------------------------------
   lp1 = l + 1
   DO  i = 1, n
     yh(i,lp1) = t1*yh(i,lmax)
   END DO
   ! Add correction terms to YH array. ------------------------------------
   nqp1 = nq + 1
   DO  j = 3, nqp1
   !  CALL daxpy (n, el(j), yh(1,lp1), yh(1,j) )
   DO  i = 1,n
     yh(i,j)=yh(i,j)+el(j)*yh(i,lp1)
   END DO
   END DO
   RETURN
   !----------------------- End of Subroutine DVJUST ----------------------
 END SUBROUTINE dvjust
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! dewset                                                                       !
! Call sequence input -- N, ITOL, RTOL, ATOL, YCUR                             !
! Call sequence output -- EWT                                                  !
! COMMON block variables accessed -- None                                      !
                                                                               !
! Subroutines/functions called by DEWSET.. None                                !
!-----------------------------------------------------------------------       !
! This subroutine sets the error weight vector EWT according to                !
!     EWT(i) = RTOL(i)*abs(YCUR(i)) + ATOL(i),  i = 1,...,N,                   !
! with the subscript on RTOL and/or ATOL possibly replaced by 1 above,         !
! depending on the value of ITOL.                                              !
!------------------------------------------------------------------------------!
 SUBROUTINE dewset (n, itol, rtol, atol, ycur, ewt)

   INTEGER,    INTENT(IN)             :: n
   INTEGER,    INTENT(IN)             :: itol
   REAL(SEDP), INTENT(IN)             :: rtol(*)
   REAL(SEDP), INTENT(IN)             :: atol(*)
   REAL(SEDP), INTENT(INOUT)          :: ycur(n)
   REAL(SEDP), INTENT(OUT)            :: ewt(n)

   !-----------------------------------------------------------------------
   INTEGER  :: i

   SELECT CASE ( itol )
     CASE (    1)
       GO TO 10
     CASE (    2)
       GO TO  20
     CASE (    3)
       GO TO  30
     CASE (    4)
       GO TO  40
   END SELECT
   10 CONTINUE

   DO  i = 1, n
     ewt(i) = rtol(1)*ABS(ycur(i)) + atol(1)
   END DO
   RETURN
   20 CONTINUE
   DO  i = 1, n
     ewt(i) = rtol(1)*ABS(ycur(i)) + atol(i)
   END DO
   RETURN
   30 CONTINUE
   DO  i = 1, n
     ewt(i) = rtol(i)*ABS(ycur(i)) + atol(1)
   END DO
   RETURN
   40 CONTINUE
   DO  i = 1, n
     ewt(i) = rtol(i)*ABS(ycur(i)) + atol(i)
   END DO
   RETURN
   !----------------------- End of Subroutine DEWSET ----------------------

 END SUBROUTINE dewset
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! DVNORM                                                                       !
! Call sequence input -- N, V, W                                               !
! Call sequence output -- None                                                 !
! COMMON block variables accessed -- None                                      !
                                                                               !
! Subroutines/functions called by DVNORM.. None                                !
!-----------------------------------------------------------------------       !
! This function routine computes the weighted root-mean-square norm            !
! of the vector of length N contained in the array V, with weights             !
! contained in the array W of length N..                                       !
!   DVNORM = sqrt( (1/N) * sum( V(i)*W(i) )**2 )                               !
!------------------------------------------------------------------------------!
 REAL(SEDP) FUNCTION dvnorm (n, v, w)

   INTEGER,    INTENT(IN)             :: n
   REAL(SEDP), INTENT(IN)             :: v(n)
   REAL(SEDP), INTENT(IN)             :: w(n)

   !-----------------------------------------------------------------------
   REAL(SEDP) :: sum
   INTEGER  :: i

   sum = 0.0D+00
   DO  i = 1, n
     sum = sum + (v(i)*w(i))**2
   END DO
   dvnorm = SQRT(sum/FLOAT(n))
   RETURN
   !----------------------- End of Function DVNORM ------------------------

 END FUNCTION dvnorm
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! DVSRCO                                                                       !
! Call sequence input -- RSAV, ISAV, JOB                                       !
! Call sequence output -- RSAV, ISAV                                           !
! COMMON block variables accessed -- All of /DVOD01/ and /DVOD02/              !
                                                                               !
! Subroutines/functions called by DVSRCO.. None                                !
!-----------------------------------------------------------------------       !
! This routine saves or restores (depending on JOB) the contents of the        !
! COMMON blocks DVOD01 and DVOD02, which are used internally by DVODE.         !
                                                                               !
! RSAV = REAL(SEDP) array of length 49 or more.                                !
! ISAV = integer array of length 41 or more.                                   !
! JOB  = flag indicating to save or restore the COMMON blocks..                !
!        JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV).             !
!        JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV).           !
!        A CALL with JOB = 2 presumes a prior call with JOB = 1.               !
!------------------------------------------------------------------------------!
 SUBROUTINE dvsrco (rsav, isav, job)

   REAL(SEDP), INTENT(OUT)                     :: rsav(*)
   INTEGER,    INTENT(OUT)                     :: isav(*)
   INTEGER,    INTENT(IN)                      :: job

   !-----------------------------------------------------------------------
   REAL(SEDP) :: rvod1, rvod2
   INTEGER  :: ivod1, ivod2
   INTEGER  :: i
   !-----------------------------------------------------------------------
   ! The following Fortran-77 declaration is to cause the values of the
   ! listed (local) variables to be saved between calls to this integrator.
   !-----------------------------------------------------------------------

   COMMON /dvod01/ rvod1(48), ivod1(33)
   COMMON /dvod02/ rvod2(1), ivod2(8)

   IF (job == 1) THEN
     DO i = 1,48
       rsav(i) = rvod1(i)
     END DO
     rsav(49) = rvod2(1)
     DO i = 1,33
       isav(i) = ivod1(i)
     END DO
     DO i = 1,8
       isav(33+i) = ivod2(i)
     END DO
   ELSE
     DO i = 1,48
       rvod1(i) = rsav(i)
     END DO
     rvod2(1) = rsav(49)
     DO i = 1,33
       ivod1(i) = isav(i)
     END DO
     DO i = 1,8
       ivod2(i) = isav(33+i)
     END DO
   END IF
   RETURN
   !----------------------- End of Subroutine DVSRCO ----------------------
 END SUBROUTINE dvsrco
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! DACOPY                                                                       !
! Call sequence input -- NROW, NCOL, A, NROWA, NROWB                           !
! Call sequence output -- B                                                    !
! COMMON block variables accessed -- None                                      !
                                                                               !
! Subroutines called by DACOPY.. DCOPY                                         !
! Function routines called by DACOPY.. None                                    !
!-----------------------------------------------------------------------       !
! This routine copies one rectangular array, A, to another, B,                 !
! where A and B may have different row dimensions, NROWA and NROWB.            !
! The data copied consists of NROW rows and NCOL columns.                      !
!------------------------------------------------------------------------------!
 SUBROUTINE dacopy (nrow, ncol, a, nrowa, b, nrowb)

   INTEGER,    INTENT(IN)                     :: nrow
   INTEGER,    INTENT(IN)                     :: ncol
   INTEGER,    INTENT(IN)                     :: nrowa
   REAL(SEDP), INTENT(IN)                     :: a(nrowa,ncol)
   INTEGER,    INTENT(IN)                     :: nrowb
   REAL(SEDP), INTENT(INOUT)                  :: b(nrowb,ncol)

   !-----------------------------------------------------------------------
   INTEGER  :: ic,i

   DO  ic = 1,ncol
   !  CALL dcopy (nrow, a(1,ic), 1, b(1,ic), 1)
     do i = 1,nrow
       b(i,ic) = a(i,ic)
     enddo
   END DO

   RETURN
   !----------------------- End of Subroutine DACOPY ----------------------
 END SUBROUTINE dacopy
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! DGEFA                                                                        !
!     DGEFA FACTORS A REAL(SEDP) MATRIX BY GAUSSIAN ELIMINATION.               !
!------------------------------------------------------------------------------!
 SUBROUTINE dgefa(a,lda,n,ipvt,info)

   INTEGER,    INTENT(IN)                      :: lda
   INTEGER,    INTENT(IN)                      :: n
   REAL(SEDP), INTENT(IN OUT)                  :: a(lda,*)
   INTEGER,    INTENT(OUT)                     :: ipvt(*)
   INTEGER,    INTENT(OUT)                     :: info

   REAL(SEDP) :: t
   !INTEGER  :: idamax1,i,j,k,kp1,l,nm1
   INTEGER  :: i,j,k,kp1,l,nm1

   !     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING

   info = 0
   nm1 = n - 1
   IF (nm1 < 1) GO TO 70
   DO  k = 1, nm1
     kp1 = k + 1

   !        FIND L = PIVOT INDEX

     l = idamax1(n-k+1,a(k,k)) + k - 1
     ipvt(k) = l

   !        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED

     IF (a(l,k) == 0.0D+00) GO TO 40

   !           INTERCHANGE IF NECESSARY

   !           IF (L .EQ. K) GO TO 10
     t = a(l,k)
     a(l,k) = a(k,k)
     a(k,k) = t
   !  10       CONTINUE

   !           COMPUTE MULTIPLIERS

     t = -1.0D+00/a(k,k)
   !  CALL dscal(n-k,t,a(k+1,k),1)
     DO i=1,n-k
       a(k+i,k)=a(k+i,k)*t
     END DO

   !           ROW ELIMINATION WITH COLUMN INDEXING

     DO  j = kp1, n
       t = a(l,j)
   !              IF (L .EQ. K) GO TO 20
       a(l,j) = a(k,j)
       a(k,j) = t
       !  20          CONTINUE
       !    CALL daxpy(n-k,t,a(k+1,k),a(k+1,j))
       DO  i = 1,n-k
         a(i+k,j)=a(i+k,j)+t*a(i+k,k)
       END DO
     END DO
     GO TO 50
     40    CONTINUE
     info = k
     50    CONTINUE
   END DO
   70 CONTINUE
   ipvt(n) = n
   IF (a(n,n) == 0.0D+00) info = n
   RETURN
 END SUBROUTINE dgefa
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
!     DGESL SOLVES THE REAL(SEDP) SYSTEM                                       !
!     A * X = B  OR  TRANS(A) * X = B                                          !
!     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.                            !
!------------------------------------------------------------------------------!
 SUBROUTINE dgesl(a,lda,n,ipvt,b,job)

   INTEGER,    INTENT(INOUT)                   :: lda
   INTEGER,    INTENT(IN)                      :: n
   INTEGER,    INTENT(IN)                      :: ipvt(*)
   REAL(SEDP), INTENT(IN)                      :: a(lda,*)
   REAL(SEDP), INTENT(INOUT)                   :: b(*)
   INTEGER,    INTENT(IN)                      :: job

   REAL(SEDP) :: t
   !REAL(SEDP) :: ddot,t
   INTEGER  :: k,kb,l,nm1,i

   nm1 = n - 1
   IF (job /= 0) GO TO 50

   !        JOB = 0 , SOLVE  A * X = B
   !        FIRST SOLVE  L*Y = B

   IF (nm1 < 1) GO TO 30
   DO  k = 1, nm1
     l = ipvt(k)
     t = b(l)
     IF (l == k) GO TO 10
     b(l) = b(k)
     b(k) = t
     10       CONTINUE
   !  CALL daxpy(n-k,t,a(k+1,k),b(k+1))
   DO  i = 1,n-k
     b(i+k)=b(i+k)+t*a(i+k,k)
   END DO
   END DO
   30    CONTINUE

   !        NOW SOLVE  U*X = Y

   DO  kb = 1, n
     k = n + 1 - kb
     b(k) = b(k)/a(k,k)
     t = -b(k)
     !  CALL daxpy(k-1,t,a(1,k),b(1))
     DO  i = 1,k-1
       b(i)=b(i)+t*a(i,k)
     END DO
   END DO
   GO TO 100
   50 CONTINUE

   !        JOB = NONZERO, SOLVE  TRANS(A) * X = B
   !        FIRST SOLVE  TRANS(U)*Y = B

   DO  k = 1, n
     t = ddot(k-1,a(1,k),b(1))
     b(k) = (b(k) - t)/a(k,k)
   END DO

   !        NOW SOLVE TRANS(L)*X = Y

   IF (nm1 < 1) GO TO 90
   DO  kb = 1, nm1
     k = n - kb
     b(k) = b(k) + ddot(n-k,a(k+1,k),b(k+1))
     l = ipvt(k)
     IF (l == k) GO TO 70
     t = b(l)
     b(l) = b(k)
     b(k) = t
     70       CONTINUE
   END DO
   90    CONTINUE
   100 CONTINUE
   RETURN
 END SUBROUTINE dgesl
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! This routine computes the unit roundoff of the machine.                      !
! This is defined as the smallest positive machine number                      !
! u such that  1.0 + u .ne. 1.0                                                !
!                                                                              !
! Subroutines/functions called by D1MACH.. None                                !
!------------------------------------------------------------------------------!
 REAL(SEDP) FUNCTION d1mach (idum)

   INTEGER, INTENT(IN)                  :: idum

   REAL(SEDP) :: u, comp

   u = 1.0D+00
   10   u = u*0.5D+00
   comp = 1.0D+00 + u
   IF (comp /= 1.0D+00) GO TO 10
   d1mach = u*2.0D+00
   RETURN
   !----------------------- End of Function D1MACH ------------------------
 END FUNCTION d1mach
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! XERRWD                                                                       !
! Subroutines XERRWD, XSETF, XSETUN, and the two function routines             !
! MFLGSV and LUNSAV, as given here, constitute a simplified version of         !
! the SLATEC error handling package.                                           !
! Written by A. C. Hindmarsh and P. N. Brown at LLNL.                          !
! Version of 13 April, 1989.                                                   !
! This version is in REAL(SEDP).                                               !
                                                                               !
! All arguments are input arguments.                                           !
                                                                               !
! MSG    = The message (character array).                                      !
! NMES   = The length of MSG (number of characters).                           !
! NERR   = The error number (not used).                                        !
! LEVEL  = The error level..                                                   !
!          0 or 1 means recoverable (control returns to caller).               !
!          2 means fatal (run is aborted--see note below).                     !
! NI     = Number of integers (0, 1, or 2) to be printed with message.         !
! I1,I2  = Integers to be printed, depending on NI.                            !
! NR     = Number of REAL(SEDP)s (0, 1, or 2) to be printed with message.      !
! R1,R2  = FLOATs to be printed, depending on NR.                              !
                                                                               !
! Note..  this routine is machine-dependent and specialized for use            !
! in limited context, in the following ways..                                  !
! 1. The argument MSG is assumed to be of type CHARACTER, and                  !
!    the message is printed with a format of (1X,80A1).                        !
! 2. The message is assumed to take only one line.                             !
!    Multi-line messages are generated by repeated calls.                      !
! 3. If LEVEL = 2, control passes to the statement   STOP                      !
!    to abort the run.  This statement may be machine-dependent.               !
! 4. R1 and R2 are assumed to be in REAL(SEDP) and are printed                 !
!    in D21.13 format.                                                         !
                                                                               !
! For a different default logical unit number, change the data                 !
! statement in function routine LUNSAV.                                        !
! For a different run-abort command, change the statement following            !
! statement 100 at the end.                                                    !
!-----------------------------------------------------------------------       !
! Subroutines called by XERRWD.. None                                          !
! Function routines called by XERRWD.. MFLGSV, LUNSAV                          !
!------------------------------------------------------------------------------!
 SUBROUTINE xerrwd (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)

   INTEGER,    INTENT(IN)                  :: nmes
   INTEGER,    INTENT(IN)                  :: nerr
   INTEGER,    INTENT(IN)                  :: level
   INTEGER,    INTENT(IN)                  :: ni
   INTEGER,    INTENT(IN)                  :: i1
   INTEGER,    INTENT(IN)                  :: i2
   INTEGER,    INTENT(IN)                  :: nr
   REAL(SEDP), INTENT(IN)                  :: r1
   REAL(SEDP), INTENT(IN)                  :: r2
   CHARACTER (LEN=*), INTENT(IN)           :: msg !(nmes)
   !##MATT CHARACTER (LEN=*), INTENT(IN)        :: msg(nmes)

   !-----------------------------------------------------------------------

   !INTEGER  :: i, lunit, lunsav, mesflg, mflgsv
   INTEGER  :: i, lunit, mesflg

   ! Get message print flag and logical unit number. ----------------------
   mesflg = mflgsv (0,.false.)
   lunit = lunsav (0,.false.)
   IF (mesflg == 0) GO TO 100
   ! Write the message. ---------------------------------------------------
   !WRITE (lunit,10) (msg(i),i=1,nmes)
   WRITE (lunit,10) msg
   10   FORMAT(1X,80A1)
   IF (ni == 1) WRITE (lunit, 20) i1
   20   FORMAT(6X,'In above message,  I1 =',i10)
   IF (ni == 2) WRITE (lunit, 30) i1,i2
   30   FORMAT(6X,'In above message,  I1 =',i10,3X,'I2 =',i10)
   IF (nr == 1) WRITE (lunit, 40) r1
   40   FORMAT(6X,'In above message,  R1 =',e21.13)
   IF (nr == 2) WRITE (lunit, 50) r1,r2
   50   FORMAT(6X,'In above,  R1 =',e21.13,3X,'R2 =',e21.13)
   ! Abort the run if LEVEL = 2. ------------------------------------------
   100  IF (level /= 2) RETURN
   STOP
   !----------------------- End of Subroutine XERRWD ----------------------
 END SUBROUTINE xerrwd
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! This routine resets the logical unit number for messages.                    !
                                                                               !
! Subroutines called by XSETUN.. None                                          !
! Function routines called by XSETUN.. LUNSAV                                  !
!------------------------------------------------------------------------------!
 SUBROUTINE xsetun (lun)

   INTEGER, INTENT(IN)                      :: lun
   !INTEGER  :: junk, lunsav
   INTEGER  :: junk

   IF (lun > 0) junk = lunsav (lun,.true.)
   RETURN
   !----------------------- End of Subroutine XSETUN ----------------------
 END SUBROUTINE xsetun
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! This routine resets the print control flag MFLAG.

! Subroutines called by XSETF.. None
! Function routines called by XSETF.. MFLGSV
!------------------------------------------------------------------------------!
 SUBROUTINE xsetf (mflag)

   INTEGER, INTENT(IN)                      :: mflag
   !INTEGER  :: junk, mflgsv
   INTEGER  :: junk

   IF (mflag == 0 .OR. mflag == 1) junk = mflgsv (mflag,.true.)
   RETURN
   !----------------------- End of Subroutine XSETF -----------------------
 END SUBROUTINE xsetf
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! LUNSAV                                                                       !
! LUNSAV saves and recalls the parameter LUNIT which is the logical            !
! unit number to which error messages are printed.                             !
                                                                               !
! Saved local variable..                                                       !
                                                                               !
!  LUNIT   = Logical unit number for messages.                                 !
!            The default is 6 (machine-dependent).                             !
                                                                               !
! On input..                                                                   !
                                                                               !
!   IVALUE = The value to be set for the LUNIT parameter,                      !
!            if ISET is .TRUE. .                                               !
                                                                               !
!   ISET   = Logical flag to indicate whether to read or write.                !
!            If ISET=.TRUE., the LUNIT parameter will be given                 !
!            the value IVALUE.  If ISET=.FALSE., the LUNIT                     !
!            parameter will be unchanged, and IVALUE is a dummy                !
!            parameter.                                                        !
                                                                               !
! On return..                                                                  !
                                                                               !
!   The (old) value of the LUNIT parameter will be returned                    !
!   in the function value, LUNSAV.                                             !
                                                                               !
! This is a modification of the SLATEC library routine J4SAVE.                 !
                                                                               !
! Subroutines/functions called by LUNSAV.. None                                !
!------------------------------------------------------------------------------!
 INTEGER FUNCTION lunsav (ivalue, iset)

   INTEGER, INTENT(IN)                      :: ivalue
   LOGICAL, INTENT(IN)                      :: iset

   !-----------------------------------------------------------------------
   INTEGER  :: lunit

   !-----------------------------------------------------------------------
   ! The following Fortran-77 declaration is to cause the values of the
   ! listed (local) variables to be saved between calls to this integrator.
   !-----------------------------------------------------------------------
   SAVE lunit
   DATA lunit/6/

   lunsav = lunit
   IF (iset) lunit = ivalue
   RETURN
   !----------------------- End of Function LUNSAV ------------------------
 END FUNCTION lunsav
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! MFLGSV                                                                       !
! MFLGSV saves and recalls the parameter MESFLG which controls the             !
! printing of the error messages.                                              !
                                                                               !
! Saved local variable..                                                       !
                                                                               !
!   MESFLG = Print control flag..                                              !
!            1 means print all messages (the default).                         !
!            0 means no printing.                                              !
                                                                               !
! On input..                                                                   !
                                                                               !
!   IVALUE = The value to be set for the MESFLG parameter,                     !
!            if ISET is .TRUE. .                                               !
                                                                               !
!   ISET   = Logical flag to indicate whether to read or write.                !
!            If ISET=.TRUE., the MESFLG parameter will be given                !
!            the value IVALUE.  If ISET=.FALSE., the MESFLG                    !
!            parameter will be unchanged, and IVALUE is a dummy                !
!            parameter.                                                        !
                                                                               !
! On return..                                                                  !
                                                                               !
!   The (old) value of the MESFLG parameter will be returned                   !
!   in the function value, MFLGSV.                                             !
                                                                               !
! This is a modification of the SLATEC library routine J4SAVE.                 !
                                                                               !
! Subroutines/functions called by MFLGSV.. None                                !
!------------------------------------------------------------------------------!
 INTEGER FUNCTION mflgsv (ivalue, iset)

   INTEGER, INTENT(IN)                      :: ivalue
   LOGICAL, INTENT(IN)                      :: iset

   !-----------------------------------------------------------------------
   INTEGER  :: mesflg
   !-----------------------------------------------------------------------
   ! The following Fortran-77 declaration is to cause the values of the
   ! listed (local) variables to be saved between calls to this integrator.
   !-----------------------------------------------------------------------
   SAVE mesflg
   DATA mesflg/1/

   mflgsv = mesflg
   IF (iset) mesflg = ivalue
   RETURN
   !----------------------- End of Function MFLGSV ------------------------
 END FUNCTION mflgsv
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
!     DGBSL SOLVES THE REAL(SEDP) BAND SYSTEM
!     A * X = B  OR  TRANS(A) * X = B
!     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.
!------------------------------------------------------------------------------!
 SUBROUTINE dgbsl(abd,lda,n,ml,mu,ipvt,b,job)

   INTEGER,    INTENT(IN)             :: job
   INTEGER,    INTENT(INOUT)          :: lda
   INTEGER,    INTENT(IN)             :: n
   INTEGER,    INTENT(IN)             :: ml
   INTEGER,    INTENT(IN)             :: mu
   INTEGER,    INTENT(IN)             :: ipvt(*)
   REAL(SEDP), INTENT(IN)             :: abd(lda,*)
   REAL(SEDP), INTENT(INOUT)          :: b(*)

   !REAL(SEDP) :: ddot,t
   REAL(SEDP) :: t
   INTEGER  :: k,kb,l,la,lb,lm,m,nm1,i

   m = mu + ml + 1
   nm1 = n - 1
   IF (job /= 0) GO TO 50

   !        JOB = 0 , SOLVE  A * X = B
   !        FIRST SOLVE L*Y = B

   IF (ml == 0) GO TO 30
   !        IF (NM1 .LT. 1) GO TO 30
   DO  k = 1, nm1
     lm = MIN(ml,n-k)
     l = ipvt(k)
     t = b(l)
     IF (l == k) GO TO 10
     b(l) = b(k)
     b(k) = t
     10          CONTINUE
   !  CALL daxpy(lm,t,abd(m+1,k),b(k+1))
   !print *, "ROGER1" ,n
   DO  i = 1,lm
     b(k+i)=b(k+i)+t*abd(m+i,k)
   END DO

   END DO
   30    CONTINUE

   !        NOW SOLVE  U*X = Y

   DO  kb = 1, n
     k = n + 1 - kb
     b(k) = b(k)/abd(m,k)
     lm = MIN(k,m) - 1
     la = m - lm
     lb = k - lm
     t = -b(k)
   !  CALL daxpy(lm,t,abd(la,k),b(lb))
   DO  i = 0,lm-1
     b(lb+i)=b(lb+i)+t*abd(la+i,k)
   END DO

   END DO
   RETURN

   !        JOB = NONZERO, SOLVE  TRANS(A) * X = B
   !        FIRST SOLVE  TRANS(U)*Y = B

   50 CONTINUE
   DO  k = 1, n
     lm = MIN(k,m) - 1
     la = m - lm
     lb = k - lm
     t = ddot(lm,abd(la,k),b(lb))
     b(k) = (b(k) - t)/abd(m,k)
   END DO

   !        NOW SOLVE TRANS(L)*X = Y

   IF (ml == 0) GO TO 90
   IF (nm1 < 1) GO TO 90
   DO  kb = 1, nm1
     k = n - kb
     lm = MIN(ml,n-k)
     b(k) = b(k) + ddot(lm,abd(m+1,k),b(k+1))
     l = ipvt(k)
     IF (l == k) GO TO 70
     t = b(l)
     b(l) = b(k)
     b(k) = t
     70          CONTINUE
   END DO
   90    CONTINUE
   RETURN
END SUBROUTINE dgbsl
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
!     DGBFA FACTORS A REAL(SEDP) BAND MATRIX BY ELIMINATION.                   !
!------------------------------------------------------------------------------!
 SUBROUTINE dgbfa(abd,lda,n,ml,mu,ipvt,info)

   INTEGER,    INTENT(IN)             :: lda
   REAL(SEDP), INTENT(OUT)            :: abd(lda,*)
   INTEGER,    INTENT(IN)             :: n
   INTEGER,    INTENT(IN)             :: ml
   INTEGER,    INTENT(IN)             :: mu
   INTEGER,    INTENT(OUT)            :: ipvt(*)
   INTEGER,    INTENT(OUT)            :: info

   REAL(SEDP) :: t
   !INTEGER  :: i,idamax1,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
   INTEGER  :: i,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1

   m = ml + mu + 1
   info = 0

   !     ZERO INITIAL FILL-IN COLUMNS

   j0 = mu + 2
   j1 = MIN(n,m) - 1
   !     IF (J1 .LT. J0) GO TO 30
   DO  jz = j0, j1
     i0 = m + 1 - jz
     DO  i = i0, ml
       abd(i,jz) = 0.0D+00
     END DO
   END DO
   jz = j1
   ju = 0

   !     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING

   nm1 = n - 1
   !     IF (NM1 .LT. 1) GO TO 130
   DO  k = 1, nm1
     kp1 = k + 1

   !        ZERO NEXT FILL-IN COLUMN

     jz = jz + 1
     IF (jz > n) GO TO 50
   !        IF (ML .LT. 1) GO TO 50
     DO  i = 1, ml
       abd(i,jz) = 0.0D+00
     END DO
     50    CONTINUE

   !        FIND L = PIVOT INDEX

     lm = MIN(ml,n-k)
     l = idamax1(lm+1,abd(m,k)) + m - 1
     ipvt(k) = l + k - m

   !        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED

     IF (abd(l,k) == 0.0D+00) GO TO 100

   !           INTERCHANGE IF NECESSARY

     IF (l == m) GO TO 60
     t = abd(l,k)
     abd(l,k) = abd(m,k)
     abd(m,k) = t
     60       CONTINUE

   !           COMPUTE MULTIPLIERS

     t = -1.0D+00/abd(m,k)
   !  CALL dscal(lm,t,abd(m+1,k),1)
     do i=1,lm
       abd(m+i,k)=abd(m+i,k)*t
     enddo

   !           ROW ELIMINATION WITH COLUMN INDEXING

     ju = MIN(MAX(ju,mu+ipvt(k)),n)
     mm = m
   !           IF (JU .LT. KP1) GO TO 90
     DO  j = kp1, ju
       l = l - 1
       mm = mm - 1
       t = abd(l,j)
       IF (l == mm) GO TO 70
       abd(l,j) = abd(mm,j)
       abd(mm,j) = t
       70          CONTINUE
   !    CALL daxpy(lm,t,abd(m+1,k),abd(mm+1,j))
   do i=1,lm
     abd(mm+i,j)=abd(mm+i,j)+t*abd(m+i,k)
   enddo
     END DO
     GO TO 110
     100    CONTINUE
     info = k
     110    CONTINUE
   END DO
   ipvt(n) = n
   IF (abd(m,n) == 0.0D+00) info = n
   RETURN
 END SUBROUTINE dgbfa
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
!     FORMS THE DOT PRODUCT OF TWO VECTORS.                                    !
!        CODE FOR BOTH INCREMENTS EQUAL TO 1                                   !
!     JACK DONGARRA, LINPACK, 3/11/78.                                         !
!------------------------------------------------------------------------------!
 REAL(SEDP) FUNCTION ddot(n,dx,dy)

   INTEGER,    INTENT(IN)             :: n
   REAL(SEDP), INTENT(IN)             :: dx(*)
   REAL(SEDP), INTENT(IN)             :: dy(*)

   REAL(SEDP) :: dtemp
   INTEGER  :: i

   ddot = 0.0D+00
   dtemp = 0.0D+00
   IF(n <= 0)RETURN
   DO  i = 1,n
     dtemp = dtemp + dx(i)*dy(i)
   END DO
   ddot = dtemp
   RETURN
 END FUNCTION ddot
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
!     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.                   !
!     JACK DONGARRA, LINPACK, 3/11/78.                                         !
!     modified for increment 1 only R. Luff                                    !
!------------------------------------------------------------------------------!
 INTEGER FUNCTION idamax1(n,dx)

   INTEGER,    INTENT(IN)               :: n
   REAL(SEDP), INTENT(IN)               :: dx(*)

   REAL(SEDP) :: dmax,tmax
   INTEGER  :: i

   idamax1 = 0
   dmax = 0.0D+00
   DO  i = 1,n
     tmax=ABS(dx(i))
     IF (tmax > dmax) THEN
       idamax1 = i
       dmax = tmax
     END IF
   END DO
   RETURN
 END FUNCTION idamax1
!------------------------------------------------------------------------------!






 !     program testdet
!     INTEGER IPVT(5),INFO
!     DOUBLE PRECISION A(5,5), work(5),det(2)
!     A(1,1) = 2.0
!     A(1,2) = 2.0
!     A(1,3) = 4.0
!     A(1,4) = 1.0
!     A(2,1) = 9.0
!     A(2,2) = -3.0
!     A(2,3) = 8.0
!     A(2,4) = 2.0
!     A(3,1) = 9.0
!     A(3,2) = 12.0
!     A(3,3) = 3.0
!     A(3,4) = 6.0
!     A(4,1) = 4.0
!     A(4,2) = 8.0
!     A(4,3) = -5.0
!     A(4,4) = 4.0
!     do i=1,4
!       print *, (a(j,i),j=1,4)
!     enddo
!     CALL DGEFA(A,5,4,ipvt,info)
!     do i=1,4
!       print *, (a(j,i),j=1,4)
!     enddo
!     CALL dgedi(a,5,4,ipvt,det,work,11)
!     print *, det(1) * 10.0**det(2)
!     end


!------------------------------------------------------------------------------!
! DGEDI                                                                        !
!     dgedi computes the determinant and inverse of a matrix                   !
!     using the factors computed by dgeco or dgefa.                            !
                                                                               !
!     on entry                                                                 !
                                                                               !
!        a       DOUBLE PRECISION(lda, n)                                      !
!                the output from dgeco or dgefa.                               !
                                                                               !
!        lda     integer                                                       !
!                the leading dimension of the array  a .                       !
!                groesse des feldes a                                          !
                                                                               !
!        n       integer                                                       !
!                the order of the matrix  a .                                  !
!                symetrische matrix spalten=zeilen=a                           !
                                                                               !
!        ipvt    integer(n)                                                    !
!                the pivot vector from dgeco or dgefa.                         !
                                                                               !
!        work    DOUBLE PRECISION(n)                                           !
!                work vector.  contents destroyed.                             !
                                                                               !
!        job     integer                                                       !
!                = 11   both determinant and inverse.                          !
!                = 01   inverse only.                                          !
!                = 10   determinant only.                                      !
                                                                               !
!     on return                                                                !
                                                                               !
!        a       inverse of original matrix if requested.                      !
!                otherwise unchanged.                                          !
                                                                               !
!        det     DOUBLE PRECISION(2)                                           !
!                determinant of original matrix if requested.                  !
!                otherwise not referenced.                                     !
!                determinant = det(1) * 10.0**det(2)                           !
!                with  1.0 .le. abs(det(1)) .lt. 10.0                          !
!                or  det(1) .eq. 0.0 .                                         !
                                                                               !
!     error condition                                                          !
                                                                               !
!        a division by zero will occur if the input factor contains            !
!        a zero on the diagonal and the inverse is requested.                  !
!        it will not occur if the subroutines are called correctly             !
!        and if dgeco has set rcond .gt. 0.0 or dgefa has set                  !
!        info .eq. 0 .                                                         !
                                                                               !
!     linpack. this version dated 08/14/78 .                                   !
!     cleve moler, university of new mexico, argonne national lab.             !
                                                                               !
!     subroutines and functions                                                !
                                                                               !
!     blas daxpy,dscal,dswap                                                   !
!     fortran abs,mod                                                          !
!------------------------------------------------------------------------------!
 SUBROUTINE dgedi(a,lda,n,ipvt,det,work,job)

   INTEGER,    INTENT(IN)                  :: lda
   INTEGER,    INTENT(IN)                  :: n
   INTEGER,    INTENT(IN)                  :: ipvt(*)
   INTEGER,    INTENT(IN)                  :: job
   REAL(SEDP), INTENT(INOUT)               :: a(lda,*)
   REAL(SEDP), INTENT(OUT)                 :: det(2)
   REAL(SEDP), INTENT(OUT)                 :: work(*)

   !     internal variables

   REAL(SEDP) :: t
   REAL(SEDP) :: ten
   INTEGER  :: i,j,k,kb,kp1,l,nm1,ii

   !     compute determinant

   IF (job/10 == 0) GO TO 70
   det(1) = 1.0D+00
   det(2) = 0.0D+00
   ten = 10.0D+00
   DO  i = 1, n
     IF (ipvt(i) /= i) det(1) = -det(1)
     det(1) = a(i,i)*det(1)
   !        ...exit
     IF (det(1) == 0.0D+00) EXIT
     10       IF (ABS(det(1)) >= 1.0D+00) GO TO 20
     det(1) = ten*det(1)
     det(2) = det(2) - 1.0D+00
     GO TO 10
     20       CONTINUE
     30       IF (ABS(det(1)) < ten) GO TO 40
     det(1) = det(1)/ten
     det(2) = det(2) + 1.0D+00
     GO TO 30
     40       CONTINUE
   END DO
   70 CONTINUE

   !     compute inverse(u)

   IF (MOD(job,10) == 0) GO TO 150
   DO  k = 1, n
     a(k,k) = 1.0D+00/a(k,k)
     t = -a(k,k)
     DO i = 1,k-1
       a(i,k) = t*a(i,k)
     END DO
     kp1 = k + 1
     DO  j = kp1, n
       t = a(k,j)
       a(k,j) = 0.0D+00
       CALL daxpy(k,t,a(1,k),a(1,j))
     END DO
   END DO

   !        form inverse(u)*inverse(l)

   nm1 = n - 1
   DO  kb = 1, nm1
     k = n - kb
     kp1 = k + 1
     DO  i = kp1, n
       work(i) = a(i,k)
       a(i,k) = 0.0D+00
     END DO
     DO  j = kp1, n
       t = work(j)
       CALL daxpy(n,t,a(1,j),a(1,k))
     END DO
     l = ipvt(k)
     IF (l /= k) CALL dswap(n,a(1,k),1,a(1,l),1)
   END DO
   150 CONTINUE
   RETURN
 END SUBROUTINE dgedi
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
!     interchanges two vectors.                                                !
!     uses unrolled loops for increments equal one.                            !
!     jack dongarra, linpack, 3/11/78.                                         !
!------------------------------------------------------------------------------!
 SUBROUTINE  dswap (n,dx,incx,dy,incy)

   INTEGER,    INTENT(IN)                      :: incy
   INTEGER,    INTENT(IN)                      :: n
   INTEGER,    INTENT(IN)                      :: incx
   REAL(SEDP), INTENT(INOUT)                   :: dx(n)
   REAL(SEDP), INTENT(INOUT)                   :: dy(n)
   REAL(SEDP)                                  :: dtemp1(n)
   REAL(SEDP)                                  :: dtemp2
   INTEGER  :: i, ix,iy

   IF(n <= 0)RETURN
   IF(incx == 1.AND.incy == 1) THEN
   ! code for both increments equal to 1
     dtemp1(:) = dx(:)
     dx(:)     = dy(:)
     dy(:)     = dtemp1(:)
   ELSE
   ! code for unequal increments or equal increments not equal to 1
     ix = 1
     iy = 1
     IF(incx < 0)ix = (-n+1)*incx + 1
     IF(incy < 0)iy = (-n+1)*incy + 1
     DO  i = 1,n
       dtemp2 = dx(ix)
       dx(ix) = dy(iy)
       dy(iy) = dtemp2
       ix = ix + incx
       iy = iy + incy
     END DO
   ENDIF
   RETURN
 END SUBROUTINE  dswap
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
!     COPIES A VECTOR, X, TO A VECTOR, Y.                                      !
!     JACK DONGARRA, LINPACK, 3/11/78.                                         !
!------------------------------------------------------------------------------!
 SUBROUTINE dcopy(n,dx,incx,dy,incy)

   INTEGER,    INTENT(IN)             :: n
   INTEGER,    INTENT(IN)             :: incx
   INTEGER,    INTENT(IN)             :: incy
   REAL(SEDP), INTENT(IN)             :: dx(n)
   REAL(SEDP), INTENT(OUT)            :: dy(n)

   INTEGER  :: i, ix,iy,m,mp1

   IF(n <= 0)RETURN
   IF(incx == 1.AND.incy == 1) THEN
   ! CODE FOR BOTH INCREMENTS EQUAL TO 1
     dy(:) = dx(:)
   else
   ! CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1
     ix = 1
     iy = 1
     IF(incx < 0)ix = (-n+1)*incx + 1
     IF(incy < 0)iy = (-n+1)*incy + 1
     DO  i = 1,n
       dy(iy) = dx(ix)
       ix = ix + incx
       iy = iy + incy
     END DO
   ENDIF

   RETURN
 END SUBROUTINE dcopy
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
!     SCALES A VECTOR BY A CONSTANT.                                           !
!     JACK DONGARRA, LINPACK, 3/11/78.                                         !
!------------------------------------------------------------------------------!
 SUBROUTINE dscal(n,da,dx,incx)

   INTEGER,    INTENT(IN)             :: n
   REAL(SEDP), INTENT(IN)             :: da
   REAL(SEDP), INTENT(INOUT)          :: dx(n)
   INTEGER,    INTENT(IN)             :: incx

   INTEGER  :: i, m,mp1, nincx

   IF(n <= 0)RETURN
   IF(incx == 1) THEN
   ! CODE FOR INCREMENT EQUAL TO 1
     dx(:) = da*dx(:)
   ELSE
   ! CODE FOR INCREMENT NOT EQUAL TO 1
     nincx = n*incx
     DO  i = 1,nincx,incx
       dx(i) = da*dx(i)
     END DO
   ENDIF

   RETURN
 END SUBROUTINE dscal
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
!     CONSTANT TIMES A VECTOR PLUS A VECTOR.                                   !
!     JACK DONGARRA, LINPACK, 3/11/78.                                         !
!------------------------------------------------------------------------------!
 SUBROUTINE daxpy(n,da,dx,dy)

   INTEGER,    INTENT(IN)             :: n
   REAL(SEDP), INTENT(IN)             :: da
   REAL(SEDP), INTENT(IN)             :: dx(n)
   REAL(SEDP), INTENT(INOUT)          :: dy(n)

   INTEGER  :: i

   IF (n <= 0)RETURN
   IF (da == 0.0D+00) RETURN
   dy(:) = dy(:) + da*dx(:)

   RETURN
 END SUBROUTINE daxpy
!------------------------------------------------------------------------------!


END MODULE aed2_vode
!------------------------------------------------------------------------------!
