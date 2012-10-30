C *** THIS IS WITH NEW ERROR CONTROL NORMS, 
C     NEW CHANGE WITH IPAR(10),
C     AND CHANGE IN CASE OF SUCCESSFUL STEP  *** LD, April 2004
C *** THIS ALSO WITH CHANGES IN NONLINEAR DISCRETE QR
C     IF NO ERROR CONTROL ON LEs *** LD, May 15 2004
C *** ALSO WITH NO "HAT" COMPUTATION IN FIXED STEPSIZE *** LD, May 17 2004
C *** ALSO WITH NO UNNEEDED COMPUTATION FOR TRAP RULE *** LD, May 17 2004
C *** ALSO WITH MODIFICATIONS FOR IPAR(3)=1: NEW CALLS 
C     STRUCTURED WITH A SUBROUTINE CALL  *** LD, May 26, 27 2004
C *** ALSO WITH DIFFERENT TOLS FOR TRAJECTORY, Q, LEs *** LD, May 27 2004 


      SUBROUTINE LESNLS(GETF,GETDF,M,N,APPLES,T0,TE,DT,Y0,X0,TOLT,
     *                  TOLQ,TOLL,IPAR,FWORK,IFLAG,INARR,REARR)

C
C ***** WARNING: before calling LESNLS for the very first time,
C                the user must call INIT, which will be called only
C                once, regardless of integration options
C                In INIT, the code checks for valid Input, sets
C                up workspace and defaults.  The user must call
C                INIT after having set the Input data, in particular
C                the array IPAR.
C                See the sample driver "DRIV" for how this is
C                done *****
C
C	LESNLS approximates N Lyapunov exponents of a 
C       nonlinear system y'=f(y) on the interval [T0,TE]
C
C	LESNLS refers to: Lyapunov ExponentS NonLinear Small
C       in that it approximates Lyapunov exponents of a 
C       nonlinear system whose Jacobian matrix Df(y(t)) is
C       known and it is small enough that can be stored.
C       In case of large systems, for which the Jacobian
C       matrix cannot be stored and/or is not known, the user
C       should consider using LESNLL (Lyapunov Exponents
C       NonLinear Large).  For linear systems, see the
C       documentations of LESLIS and LESLIL. 
C
C	LESNLS is an interface to subroutines which approximate
C       the Lyapunov exponents by using either the discrete or 
C       the continuous QR method with a number of integration 
C       options.  See IPAR(*), which is the array allowing the user 
C       to specify which options are desired.
C
C       In a nutshell, one can proceed in fixed or variable 
C       step size.  In the latter case, different options of
C       local error control are allowed: error control can be 
C       performed on the trajectory, or on the Lyapunov exponents,
C       or on the Q-factor, or on combinations of these.  See 
C       IPAR(10), but see also IPAR(6), IPAR(8) and IPAR(9) for
C       compatibilities.
C
C       Solution trajectory is always approximated, and
C       LESNLS can be used also just to approximate the solution
C       of y'=f(y) (see setting of IPAR(6))
C
C       The details of the algorithms used can be found in 
C
C       [1] L. Dieci, M. Jolly, E. Van Vleck
C           "LESNLS and LESNLL: Suites of codes for approximating Lyapunov 
C            exponents of nonlinear systems", Technical Report .... 
C
C       [2] L. Dieci, E. Van Vleck
C           "LESLIS and LESLIL: Suites of codes for approximating Lyapunov 
C            exponents of linear systems", Technical Report .... 
C
C	Authors:
C
C		Luca Dieci
C		School of Mathematics
C 		Georgia Institute of Technology
C		Atlanta, GA, 30332
C		dieci@math.gatech.edu
C
C		Michael S. Jolly
C		Department of Mathematics
C 		Indiana University
C		Bloomington, IN, 47405
C		msjolly@indiana.edu
C
C		Erik S. Van Vleck
C		Department of Mathematics
C		University of Kansas
C		Lawrence, KS, 66045
C		evanvleck@math.ukansas.edu
C
C	Inquiries, Bugs, Malfunctioning, Fixes, HowTo:
C
C		Luca Dieci
C		School of Mathematics
C 		Georgia Institute of Technology
C		Atlanta, GA, 30332
C		dieci@math.gatech.edu
C

C   ***** User Interface *****
C
C	Variables:
C
C       GETF - subroutine that returns the vector field; must
C           be declared EXTERNAL and have parameter list of the form
C                        GETF(M,Y,YDOT,INARR,REARR)
C           where M is the dimension of the system,
C           Y are the unknown variables, and YDOT is the vector field.
C           INARR and REARR are integer and double precision arrays
C              that the user can use for communication between the calling
C              program and the routine GETDF.  Treat them as dummy
C              arguments if they are not needed.
C
C       GETDF - subroutine that returns the Jacobian matrix Df; must
C           be declared EXTERNAL and have parameter list of the form
C                       GETDF(M,Y,A,INARR,REARR)
C           where M is the dimension of the matrix,
C           Y are the M dependent variables and 
C           A is the (M,M) Jacobian matrix  
C	
C       M - integer. Number of rows, dimension of the problem
C           Must have 1.LE.M
C	
C       N - integer. Number of Lyapunov exponents desired
C           Must have 1.LE.N.LE.M 
C       --> IMPORTANT: must set N.GE.1 even if set IPAR(6)=1
C
C       APPLES - double precision vector of dimension N. 
C           These are the approximate Lyapunov exponents.  For random 
C           initial conditions of the fundamental matrix (see X0), 
C           they are typically ordered from largest to smallest, at 
C           the required output points.
C           Content irrelevant on INPUT, but DO NOT CHANGE between calls.
C
C	T0 - double precision
C            Initial time on INPUT, final time reached on OUTPUT
C            Thus, on OUTPUT should have T0=TE when IPAR(2)=0.
C
C	TE - Ending time, double precision.  
C            Must be different from T0.  Can have TE<T0 or TE>T0
C            But cannot change direction of integration between calls.
C
C	DT - Chosen stepsize in fixed stepsize mode (IPAR(1)=1), 
C            stores the current stepsize in variable stepsize mode. 
C
C	Y0 - double precision, dimension M.  On INPUT, this is the
C            IC on the trajectory.  On OUTPUT, Y0 stores the 
C            approximate solution at the point reached (typically, TE).
C
C	X0 - double precision, dimension (M,N).  On INPUT,
C            these are ICs on the fundamental matrix (see IPAR(4)) to be
C            used when computing LEs.
C            On OUTPUT, it contains the Q-factor of the "transition matrix" 
C            at the point we reached (typically, TE)
C            X0 must be declared even if user does not want to approximate
C            the Lyapunov exponents (i.e., even with IPAR(6)=1) 
C
C	TOLT - double precision scalar
C             This is local error tolerance on the solution (sup-norm)
C             Must have 1.E-2>=TOLT>=1.E-13 
C
C	TOLQ - double precision scalar.
C             This is local error tolerance on the sup-norm error 
C             for the columns of Q (see IPAR(8) and IPAR(9)).  
C             Must have 1.E-2>=TOLQ>=1.E-13 
C
C	TOLL - double precision vector of dimension N.
C             These are local error tolerance on the Lyapunov exponents 
C             (see IPAR(8) and IPAR(9)).  
C             Must have 1.E-2>=TOLL(I)>=1.E-13 
C             Given the fact that it is harder to approximate the smaller
C             exponents, user may want to ask for tighter tolerances for
C             largest exponents (the first few ones)
C
C    --->     RECOMMENDED values of TOLT, TOLQ, TOLL are in the range 
C             [1.E-10, 1.E-5] if use the 5th order scheme (IPAR(8)=0,2,4)
C             and in the range [1.E-7, 1.E-3] is use the 4th order scheme
C             (IPAR(8)=1,2,5)
C
C    --->     In constant stepsize mode, IPAR(1)=1, treat TOLT,TOLQ,TOLL as 
C             dummy arguments
C
C	IPAR - integer array of dimension 13.  This serves to communicate
C             to the code several integration options.
C             Some entries can be defaulted (indicated by DEFAULT below).
C	      IPAR(1) - INPUT: set IPAR(1)=0 for variable stepsize,
C                              set IPAR(1)=1 for fixed stepsize (in DT)
C                       [if IPAR(1).NE.1, code sets DEFAULT IPAR(1)=0]
C	      IPAR(2) - INPUT: from T0 to TE --> 0, in one step mode --> 1
C                       [if IPAR(2).NE.1, code sets DEFAULT IPAR(2)=0]
C                       One-step-mode means that the code returns control to 
C                       the user after a successful step in direction of TE
C	      IPAR(3) - INPUT: must be set to 1 before calling INIT, and also
C                       for every "new" call to the code
C                       (in one step mode or otherwise, and in constant 
C                       or variable stepsize mode)
C                       A "new" call to the code is one for
C                       which values of any of IPAR(1,2,6,9,10) have been
C                       changed.   Moreover, changes in the TOLerances may have
C                       taken place as well.  Of course, T0 and TE must be
C                       specified afresh for a new call as well.
C                       Setting IPAR(3)=1 forces to zero the counters
C                       in IPAR(12) and IPAR(13)
C                 -->   N.B.:  Every new call carries a certain amount 
C                       of overhead, such as determining a new stepsize.
C             IPAR(4) - INPUT: yes/no ICs on X. 
C                       If IPAR(4) = 0, the DEFAULT, the code sets 
C                          ICs given by X0 = [I_n 0]
C                       If IPAR(4) = 1, you must specify ICs in X0
C	      IPAR(5) - INPUT: dimension of FWORK; must have
C                       IPAR(5).GE.M*M+11*M*N+13*M+8*N+63
C	      IPAR(6) - INPUT: set to 1 if want to integrate only for solution
C                       trajectory, set to 0 if want to approximate also the LEs.
C                       If IPAR(6).NE.1, the code sets the DEFAULT IPAR(6)=0.
C                       [Use of IPAR(6)=1 is justified if the user wants to make sure
C                       that the trajectory has passed its transient behavior].
C                 -->   IMPORTANT: If user then wants to approximate the LEs, then
C                       must set IPAR(6)=0, and call LESNLS again providing values 
C                       of T0, TE, ICs for the solution trajectory, setting 
C                       IPAR(3)=1, choosing IPAR(9), and possibly changing 
C                       IPAR(1,2,10) and changing values of TOLerances
C                       No other changes are allowed
C             IPAR(7) - 
C             IPAR(8) - INPUT: this specifies which method and RK formula must be
C                       used.  This refers to integration for: the trajectory, 
C                       the transitions matrix or the Q-factor, and the LEs.
C                       Of course, if IPAR(6)=1, only the solution trajectory is
C                       approximated.  As far as the approximation of LEs, 
C                       recall that we have two different approaches and three 
C                       different methods.  Projection and Hybrid are different 
C                       implementations of the continuous QR approach, while 
C                       Disc QR refers to the discrete QR method 
C                       (see the Technical Reports quoted above).
C                       DP5 refers to the Dormand-Prince embedded 
C                       pair of order 5
C                       RK38 refers to an embedded pair of order 4 
C                       (the RK 3/8 rule)
C                       IPAR(8)=0 (DEFAULT) is Projection/DP5
C                       IPAR(8)=1 is Hybrid/DP5
C                       IPAR(8)=2 is Projection/RK38
C                       IPAR(8)=3 is Hybrid/RK38
C                       IPAR(8)=4 is Disc QR/DP5
C                       IPAR(8)=5 is Disc QR/RK38
C	      IPAR(9) - INPUT: In case in which IPAR(6)=0 and 
C                       IPAR(8)=0,1,2,3, IPAR(9) specifies if LEs 
C                       are to be found by integration with DP5 or 
C                       RK38 on the nu-variables, or by the composite
C                       trapezoidal rule. 
C                       [Content of IPAR(9) is ignored if IPAR(6)=1]
C                       IPAR(9)=0, the DEFAULT, the LEs are to be found
C                          using DP5 or RK38 on the nu-variables
C                       IPAR(9)=1, the LEs are found by the composite 
C                          trapezoidal rule.
C                       In case IPAR(9)=1, and IPAR(1)=0, error control
C                       can only be performed on trajectory and Q; thus,
C                       you cannot set IPAR(10)=1,10,21,210 
C                -->    N.B.: When you choose IPAR(9)=1, then 
C                       the schemes used for the Q-integration are not the 
C                       complete projection or hybrid method, but the 
C                       simple projection/hybrid methods.  
C             IPAR(10) - INPUT: if IPAR(1)=0   
C                       then IPAR(10) specifies if user wants to do error 
C                       control on the trajectory, or
C                       on the Lyapunov exponents or on Q, or on any
C                       combination of these, subject to below restrictions:
C                       IPAR(10)=0 --> error control on trajectory.  
C                       IPAR(10)=1 --> error control on exponents 
C                       IPAR(10)=2 --> error control on Q (only if IPAR(8)=0,1,2,3)
C                       IPAR(10)=10 --> error control on trajectory and exponents 
C                       IPAR(10)=20 --> error control on trajectory and Q 
C                                       (only when IPAR(8)=0,1,2,3)
C                       IPAR(10)=21 --> error control on exponents and Q
C                                       (only when IPAR(8)=0,1,2,3)
C                       IPAR(10)=210 --> error control on trajectory, exponents, Q 
C                                       (only when IPAR(8)=0,1,2,3)
C                       if IPAR(6)=1, can choose only IPAR(10)=0
C                       if IPAR(9)=0, can choose above values for IPAR(10)
C                                     subject to stated restrictions on IPAR(8)
C                       if IPAR(9)=1, can choose IPAR(10)=0, 2, 20
C                                     subject to stated restrictions on IPAR(8)
C                       The DEFAULT set by the code is IPAR(10)=0.
C             IPAR(11) - OUTPUT: set to 0 if code reached TE,
C                               otherwise it is set to 1 
C             IPAR(12) - OUTPUT: NREJ (number of rejections)
C             IPAR(13) - OUTPUT: NST (number of steps taken)
C
C	FWORK - double precision array of dimension IPAR(5)
C               Do not touch between calls
C
C	IFLAG - OUTPUT flag
C             IFLAG=0 successful return; either one step was completed
C                     in one step mode, or solution reached TE
C             IFLAG=1 unsuccessful return; wrong dimension of FWORK
C             IFLAG=2 unsuccessful return: calling INIT with IPAR(3).NE.1
C             IFLAG=3 unsuccessful return; wrong value for TOLerances
C             IFLAG=4 mildly unsuccessful return; TE=T0, nothing done
C             IFLAG=5 unsuccessful return; abs-value stepsize
C                     smaller than 100*eps (about 2.2E-14).
C                     No further integration allowed.
C             IFLAG=6 unsuccessful return; wrong sign of DT in fixed 
C                     stepsize mode, or changed direction of integration
C                     for a "new" call (as flagged by IPAR(3)=1)
C             IFLAG=7 unsuccessful return; wrong value of N
C             IFLAG=8 unsuccessful return; you gave ICs on X0 not of full rank
C             IFLAG=9 unsuccessful return; wrong input: M.LT.1
C             IFLAG=10 unsuccessful return; something unexpected happened,
C                      most likely attempted QR factorization of a matrix
C                      which is not full rank to working precision.  This
C                      is unlikely to happen in variable stepsize: are you
C                      proceeding in constant stepsize with stepsize too large?
C             IFLAG=11 unsuccessful return.  You want to integrate only for 
C                      trajectory, but required error control on the LEs and/or Q.
C             IFLAG=12 unsuccessful return.  You want to proceed in variable 
C                      step-size, to approximate the LEs by the composite   
C                      trap-rule, but still required error control on the LEs.  
C                      Check syntax of IPAR(9) and IPAR(10).
C             IFLAG=13 unsuccessful return.  Wrong request for variable stepsize
C                      discrete QR method.
C                      Either you want to approximate the LEs by the composite 
C                      trap-rule, or want to control the error on Q.
C                      Check syntax of IPAR(9) and IPAR(10).
C
C       INARR, REARR: integer and double precision arrays
C           that the user can use for communication between the calling
C           program and the routine GETF.  Treat them as dummy
C           arguments if they are not needed.
C
C*****Author: L Dieci 
C*****Date: May 27, 2004
C
CCCCCCCC DECLARE VARIABLES
      IMPLICIT NONE
CCCCCCC INPUT/OUTPUT
      INTEGER M,N
      DOUBLE PRECISION APPLES(*),T0,TE,DT,TOLT,TOLQ,TOLL(*)
      INTEGER IPAR(*), INARR(*)
      DOUBLE PRECISION FWORK(*), Y0(M),X0(M,N), REARR(*)
      INTEGER IFLAG
      EXTERNAL GETF, GETDF

CCCCCCC LOCAL
      DOUBLE PRECISION T,HOLD
      INTEGER NREJ 

CCCCCCC POINTERS, AND COMMON VARIABLES
      LOGICAL LAST
      DOUBLE PRECISION HMIN, HNEXT, TFIRST
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      INTEGER INTDIR, IFIRST, IWHAT
      INTEGER IAPTR, IXPTR, IXHPTR, IQPTR, IBPTR
      INTEGER K1PTR, K2PTR, K3PTR, K4PTR, K5PTR, K6PTR
      INTEGER ISTAGE
      INTEGER INUPTR, INUHTP
      INTEGER K1NUPT, K2NUPT, K3NUPT, K4NUPT, K5NUPT, K6NUPT
      INTEGER IARKPT,ICRKPT,IBHTPT
      INTEGER IYPTR, IYHTPT
      INTEGER K1YPTR, K2YPTR, K3YPTR, K4YPTR, K5YPTR, K6YPTR
      INTEGER I1YPTR, IS1YPT, IS2YPT, IS3YPT, IS4YPT
      COMMON /STEP/ HMIN, HNEXT, TFIRST, INTDIR, IFIRST, IWHAT, LAST
      SAVE /STEP/
      COMMON /NSQUAX/ IAPTR, IXPTR, IQPTR, IXHPTR, IBPTR
      SAVE /NSQUAX/
      COMMON /RKPTRS/ K1PTR, K2PTR, K3PTR, K4PTR, K5PTR, K6PTR,
     *                ISTAGE
      SAVE /RKPTRS/
      COMMON /RNUPTR/INUPTR, INUHTP, K1NUPT, K2NUPT, K3NUPT,
     *              K4NUPT, K5NUPT, K6NUPT 
      SAVE /RNUPTR/
      COMMON /RKCF/ IARKPT, ICRKPT, IBHTPT
      SAVE /RKCF/
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/
      COMMON /YRKPTR/ IYPTR, IYHTPT,  
     *                K1YPTR, K2YPTR, K3YPTR, K4YPTR, K5YPTR, K6YPTR, 
     *                I1YPTR, IS1YPT, IS2YPT, IS3YPT, IS4YPT
      SAVE /YRKPTR/

      IFLAG = 0
      IF (IPAR(3).NE.1) GOTO 10
C IT IS A "NEW" CALL
      LAST=.FALSE.
C IN NEWYES WE CHECK/DEFAULT ALL THINGS TO BE DONE ON A NEW CALL
      CALL NEWYES(GETF,GETDF,M,N,APPLES,T0,TE,DT,Y0,X0,TOLT,TOLQ,
     *                     TOLL,IPAR,FWORK,IFLAG,INARR,REARR)
      IF (IFLAG.NE.0) RETURN
 10   CONTINUE
CCCCCCC LOOP
      T = T0
 20   CONTINUE  
      NREJ=0
C MONITOR IF THIS IS (PRESUMABLY) LAST STEP
C ADJUST STEPSIZE AS NEEDED, AND STORE 
C IN HNEXT THE STEP WHICH WOULD HAVE BEEN USED
      IF (LAST) THEN
         LAST=.FALSE.
         IF (IPAR(1).EQ.0) DT=HNEXT
      ENDIF
      IF (DABS(DT).GE.DABS(TE-T)) THEN
         LAST=.TRUE.
         HNEXT=DT
         IF (IPAR(1).EQ.0) DT=INTDIR*DABS(TE-T)
        ELSE
         LAST=.FALSE.
      ENDIF
      IF ((DABS(DT).LE.HMIN).AND.(.NOT.LAST)) THEN
         IFLAG=5
         RETURN
      ENDIF
C CALL TO ONE OF THE DIFFERENT METHODS.
      IF (IPAR(8).EQ.0) THEN 
CCCCCCC PQDPNS
         CALL PQDPNS(GETF,GETDF,M,N,T,DT,HOLD,Y0,X0,APPLES,TOLT,
     *    TOLQ,TOLL,FWORK(IAPTR), FWORK(IQPTR),FWORK(IXHPTR),
     *    FWORK(IBPTR),FWORK(K1PTR),
     *    FWORK(K2PTR), FWORK(K3PTR),FWORK(K4PTR),FWORK(K5PTR),
     *    FWORK(K6PTR), FWORK(ISTAGE),
     *    FWORK(INUPTR),FWORK(INUHTP),FWORK(K1NUPT),FWORK(K2NUPT), 
     *    FWORK(K3NUPT), FWORK(K4NUPT), FWORK(K5NUPT), FWORK(K6NUPT),
     *    FWORK(IARKPT),FWORK(IBHTPT),
     *    FWORK(IYPTR), FWORK(IYHTPT), FWORK(K1YPTR), FWORK(K2YPTR),
     *    FWORK(K3YPTR), FWORK(K4YPTR), FWORK(K5YPTR), FWORK(K6YPTR), 
     *    FWORK(I1YPTR), FWORK(IS1YPT), FWORK(IS2YPT), FWORK(IS3YPT), 
     *    FWORK(IS4YPT), IPAR,NREJ,IFLAG,INARR,REARR)
       ELSEIF (IPAR(8).EQ.1) THEN
CCCCCCC HQDPNS
         CALL HQDPNS(GETF,GETDF,M,N,T,DT,HOLD,Y0,X0,APPLES,TOLT,
     *    TOLQ,TOLL,FWORK(IAPTR), FWORK(IXPTR),FWORK(IXHPTR), 
     *    FWORK(IQPTR),FWORK(IBPTR),FWORK(K1PTR),
     *    FWORK(K2PTR), FWORK(K3PTR),FWORK(K4PTR),FWORK(K5PTR),
     *    FWORK(K6PTR), FWORK(ISTAGE),
     *    FWORK(INUPTR),FWORK(INUHTP),FWORK(K1NUPT),FWORK(K2NUPT), 
     *    FWORK(K3NUPT), FWORK(K4NUPT), FWORK(K5NUPT), FWORK(K6NUPT),
     *    FWORK(IARKPT),FWORK(IBHTPT),
     *    FWORK(IYPTR), FWORK(IYHTPT), FWORK(K1YPTR), FWORK(K2YPTR),
     *    FWORK(K3YPTR), FWORK(K4YPTR), FWORK(K5YPTR), FWORK(K6YPTR), 
     *    FWORK(I1YPTR), FWORK(IS1YPT), FWORK(IS2YPT), FWORK(IS3YPT), 
     *    FWORK(IS4YPT), IPAR,NREJ,IFLAG,INARR,REARR)
       ELSEIF (IPAR(8).EQ.2) THEN
CCCCCCC PQ38NS
         CALL PQ38NS(GETF,GETDF,M,N,T,DT,HOLD,Y0,X0,APPLES,TOLT,
     *    TOLQ,TOLL,FWORK(IAPTR), FWORK(IQPTR),FWORK(IXHPTR), 
     *    FWORK(IBPTR),FWORK(K1PTR),
     *    FWORK(K2PTR), FWORK(K3PTR),FWORK(K4PTR),FWORK(ISTAGE),
     *    FWORK(INUPTR),FWORK(INUHTP),FWORK(K1NUPT),FWORK(K2NUPT), 
     *    FWORK(K3NUPT), FWORK(K4NUPT),  
     *    FWORK(IARKPT),FWORK(IBHTPT),
     *    FWORK(IYPTR), FWORK(IYHTPT), FWORK(K1YPTR), FWORK(K2YPTR),
     *    FWORK(K3YPTR), FWORK(K4YPTR), FWORK(I1YPTR), FWORK(IS1YPT),
     *    FWORK(IS2YPT), IPAR,NREJ,IFLAG,INARR,REARR)
       ELSEIF (IPAR(8).EQ.3) THEN
CCCCCCC HQ38NS
         CALL HQ38NS(GETF,GETDF,M,N,T,DT,HOLD,Y0,X0,APPLES,TOLT,
     *    TOLQ,TOLL,FWORK(IAPTR), FWORK(IXPTR),FWORK(IXHPTR),  
     *    FWORK(IQPTR),FWORK(IBPTR),FWORK(K1PTR),
     *    FWORK(K2PTR), FWORK(K3PTR),FWORK(K4PTR),FWORK(ISTAGE),
     *    FWORK(INUPTR),FWORK(INUHTP),FWORK(K1NUPT),FWORK(K2NUPT), 
     *    FWORK(K3NUPT), FWORK(K4NUPT),  
     *    FWORK(IARKPT),FWORK(IBHTPT),
     *    FWORK(IYPTR), FWORK(IYHTPT), FWORK(K1YPTR), FWORK(K2YPTR),
     *    FWORK(K3YPTR), FWORK(K4YPTR), FWORK(I1YPTR), FWORK(IS1YPT),
     *    FWORK(IS2YPT), IPAR, NREJ,IFLAG,INARR,REARR)
       ELSEIF (IPAR(8).EQ.4) THEN
CCCCCCC XDPNS
         CALL XDPNS(GETF,GETDF,M,N,T,DT,HOLD,Y0,X0,APPLES,TOLT,
     *    TOLL,FWORK(IAPTR), FWORK(IXPTR),FWORK(IXHPTR),
     *    FWORK(K1PTR),FWORK(K2PTR), FWORK(K3PTR),FWORK(K4PTR),
     *    FWORK(K5PTR),FWORK(K6PTR), FWORK(ISTAGE),
     *    FWORK(IARKPT),FWORK(IBHTPT),
     *    FWORK(IYPTR), FWORK(IYHTPT), FWORK(K1YPTR), FWORK(K2YPTR),
     *    FWORK(K3YPTR), FWORK(K4YPTR), FWORK(K5YPTR), FWORK(K6YPTR), 
     *    FWORK(I1YPTR), FWORK(IS1YPT), FWORK(IS2YPT), FWORK(IS3YPT), 
     *    FWORK(IS4YPT), IPAR(1),IPAR(10),NREJ,IFLAG,INARR,REARR)
       ELSEIF (IPAR(8).EQ.5) THEN
CCCCCCC X38NS
         CALL X38NS(GETF,GETDF,M,N,T,DT,HOLD,Y0,X0,APPLES,TOLT,
     *    TOLL,FWORK(IAPTR), FWORK(IXPTR),FWORK(IXHPTR),
     *    FWORK(K1PTR),FWORK(K2PTR), FWORK(K3PTR),FWORK(K4PTR),
     *    FWORK(ISTAGE),FWORK(IARKPT),FWORK(IBHTPT),
     *    FWORK(IYPTR), FWORK(IYHTPT), FWORK(K1YPTR), FWORK(K2YPTR),
     *    FWORK(K3YPTR), FWORK(K4YPTR), FWORK(I1YPTR), FWORK(IS1YPT),
     *    FWORK(IS2YPT), IPAR(1),IPAR(10),NREJ,IFLAG,INARR,REARR)
      ENDIF
C CHECK IF ALL WENT WELL, UPDATE COUNTERS, RETURN CONTROL IF NEEDED
      IF (IFLAG.NE.0) RETURN  
      IPAR(12) = IPAR(12)+NREJ
      IPAR(13) = IPAR(13)+1
CCCCCCCCCC UPDATE FOR NEXT STEP
      T = T + HOLD
      T0=T
C CHECK IF REACHED DESIRED END POINT TE 
      IF (DABS(T).LT.DABS(TE)) THEN
         IPAR(11)=1
        ELSE
         IPAR(11)=0
      END IF
C CHECK IF FAILED ON WHAT HAD TO BE LAST STEP 
      IF((LAST).AND.(IPAR(11).EQ.1)) LAST=.FALSE.
CCCCCCC RETURN IF IN ONESTEP MODE 
      IF ((IPAR(2).EQ.0).AND.(IPAR(11).EQ.1)) GOTO 20 

      RETURN
      END
CCCCC END OF LESNLS CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: INIT
C
C*****Purpose
C
C checks for correctness of INPUT,
C sets up pointers, RK coefficients, and other quantities
C which need to be computed only once
C
C input: 
C     M, dimension of coefficient matrix
C     N, number of exponents to approximate
C     IPAR, integer vector of length 13 with information
C        on what to do
C     T0, TE : initial and final times 
C     FWORK : work array with memory to be split
C 
C output:
C     indices and quantities to retain in common blocks
C     IFLAG = 0  all went well.  If IFLAG is not 0, then
C        some mishap occurred
C
C*****Authors: L Dieci 
C*****Date: May 26, 2004
C
CCCCCC INIT CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INIT(M,N,IPAR,T0,TE,FWORK,IFLAG)
      IMPLICIT NONE
CCCCCCC INPUT/OUTPUT
      INTEGER M,N
      INTEGER IPAR(*)
      DOUBLE PRECISION T0, TE, FWORK(*)
      INTEGER IFLAG

      DOUBLE PRECISION ACC
      
CCCCCCC POINTERS, AND COMMON VARIABLES
      LOGICAL LAST
      DOUBLE PRECISION HMIN, HNEXT, TFIRST
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      INTEGER INTDIR, IFIRST, IWHAT
      INTEGER IAPTR, IXPTR, IXHPTR, IQPTR, IBPTR
      INTEGER K1PTR, K2PTR, K3PTR, K4PTR, K5PTR, K6PTR
      INTEGER ISTAGE
      INTEGER INUPTR, INUHTP
      INTEGER K1NUPT, K2NUPT, K3NUPT, K4NUPT, K5NUPT, K6NUPT
      INTEGER IARKPT,ICRKPT,IBHTPT
      INTEGER IYPTR, IYHTPT
      INTEGER K1YPTR, K2YPTR, K3YPTR, K4YPTR, K5YPTR, K6YPTR
      INTEGER I1YPTR, IS1YPT, IS2YPT, IS3YPT, IS4YPT
      COMMON /STEP/ HMIN, HNEXT, TFIRST, INTDIR, IFIRST, IWHAT, LAST
      SAVE /STEP/
      COMMON /NSQUAX/ IAPTR, IXPTR, IQPTR, IXHPTR, IBPTR
      SAVE /NSQUAX/
      COMMON /RKPTRS/ K1PTR, K2PTR, K3PTR, K4PTR, K5PTR, K6PTR,
     *                ISTAGE
      SAVE /RKPTRS/
      COMMON /RNUPTR/INUPTR, INUHTP, K1NUPT, K2NUPT, K3NUPT,
     *              K4NUPT, K5NUPT, K6NUPT 
      SAVE /RNUPTR/
      COMMON /RKCF/ IARKPT, ICRKPT, IBHTPT
      SAVE /RKCF/
      COMMON /CONSTS/ ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/
      COMMON /YRKPTR/ IYPTR, IYHTPT,  
     *                K1YPTR, K2YPTR, K3YPTR, K4YPTR, K5YPTR, K6YPTR, 
     *                I1YPTR, IS1YPT, IS2YPT, IS3YPT, IS4YPT
      SAVE /YRKPTR/

      IFLAG = 0
      LAST=.FALSE.
C CHECK CORRECTENESS OF INPUT ON M AND N
      IF (M.LT.1) THEN
         IFLAG=9
         RETURN
      END IF
      IF ((N.GT.M).OR.(N.LT.1)) THEN
         IFLAG=7
         RETURN
      END IF
C CHECK VALUES OF IPAR(4,5,8) AND SET DEFAULTS IF NEEDED
      IF (IPAR(3).NE.1) THEN
         IFLAG=2
         RETURN
      END IF
      IF (IPAR(4).NE.1) IPAR(4)=0
CCCCCCC CHECK SIZE OF FWORK 
CC       CURRENTLY NEED DIM(FWORK).GE.M*M+11*M*N+13*M+8*N+63
      IF (IPAR(5).LT.(M*M+11*M*N+13*M+8*N+63)) THEN
         IFLAG=1
         RETURN
      END IF
      IF ((IPAR(8).GT.0).AND.(IPAR(8).LT.6)) THEN
         IPAR(8) = IPAR(8)
        ELSE
         IPAR(8) = 0
      END IF
C SET UP VALUES OF CONSTANTS NEEDED IN COMMON BLOCK CONSTS
      ONE=1.0D0
      ZERO=0.0D0
C THE VALUES BELOW MAY BE CHANGED, BUT THIS IS DISCOURAGED
      EXP1=0.20D0
      EXP2=0.25D0
      HINCR=5.0D0
      HSFTY=0.8D0
      HDRSTC=0.2D0
C CHECK TO.NE.TE 
      IF (T0.EQ.TE) THEN
         IFLAG=4
         RETURN
        ELSE
C SET INTEGRATION DIRECTION 
         IF(TE.GT.T0)THEN
            INTDIR=1
           ELSE
            INTDIR=-1   
         END IF
      END IF
C COMPUTE MINIMAL ALLOWED STEPSIZE AND EPS
      HMIN=ONE
 5    HMIN=0.5D0*HMIN
      ACC=ONE+HMIN
      IF (ACC.GT.ONE) GOTO 5
      TEPS=2.0D0*HMIN
      HMIN=200.0D0*HMIN
CCCCCCC GENERATE POINTERS TO BREAK UP FWORK 
C  JACOBIAN A, TRANSITIONS X AND XHAT,  WORKSPACE QMAT AND BMAT IN /NSQUAX/  
      IAPTR=1        
      IXPTR = IAPTR+M*M
      IQPTR=IXPTR+M*N
      IXHPTR = IQPTR+M*N
      IBPTR=IXHPTR+M*N
C  QUANTITIES NEEDED TO THE RK INTEGRATORS IN /RKPTRS/  
      K1PTR = IBPTR+M*N
      K2PTR = K1PTR+M*N
      K3PTR = K2PTR+M*N
      K4PTR = K3PTR+M*N
      K5PTR = K4PTR+M*N
      K6PTR = K5PTR+M*N
      ISTAGE = K6PTR+M*N
C  QUANTITIES NEEDED TO THE RK INTEGRATORS FOR NU IN /RNUPTR/
      INUPTR = ISTAGE+M*N
      INUHTP = INUPTR+N
      K1NUPT = INUHTP+N
      K2NUPT = K1NUPT+N
      K3NUPT = K2NUPT+N
      K4NUPT = K3NUPT+N
      K5NUPT = K4NUPT+N
      K6NUPT = K5NUPT+N
C COEFFICIENTS NEEDED TO THE RK INTEGRATORS IN /RKCF/  
      IARKPT = K6NUPT+N
      ICRKPT = IARKPT+49
      IBHTPT = ICRKPT+7
C COMPUTE RK COEFFICIENTS ONCE AND FOR ALL
      CALL CFRK(IPAR(8),FWORK(IARKPT),FWORK(ICRKPT),
     *           FWORK(IBHTPT))
C  QUANTITIES NEEDED TO THE RK TRAJECTORY INTEGRATOR IN /YRKPTR/  
      IYPTR = IBHTPT+7
      IYHTPT = IYPTR+M
      K1YPTR = IYHTPT+M
      K2YPTR = K1YPTR+M
      K3YPTR = K2YPTR+M
      K4YPTR = K3YPTR+M
      K5YPTR = K4YPTR+M
      K6YPTR = K5YPTR+M
      I1YPTR = K6YPTR+M
      IS1YPT = I1YPTR+M
      IS2YPT = IS1YPT+M
      IS3YPT = IS2YPT+M
      IS4YPT = IS3YPT+M

      RETURN 
      END 
CCCCC END OF INIT CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: NEWYES
C
C*****Purpose
C
C checks for correctness of INPUT, determine Defaults,
C decide on initial stepsize, for every "new call" to LESNLS
C 
C output:
C     indices and quantities to retain in common blocks
C     IFLAG = 0  all went well.  If IFLAG is not 0, then
C        some mishap occurred
C
C*****Authors: L Dieci 
C*****Date: May 27, 2004
C
CCCCCC NEWYES CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NEWYES(GETF,GETDF,M,N,APPLES,T0,TE,DT,Y0,X0,TOLT,
     *                  TOLQ,TOLL,IPAR,FWORK,IFLAG,INARR,REARR)

CCCCCCCC DECLARE VARIABLES
      IMPLICIT NONE
CCCCCCC INPUT/OUTPUT
      INTEGER M,N
      DOUBLE PRECISION APPLES(*),T0,TE,DT,TOLT,TOLQ,TOLL(*)
      INTEGER IPAR(*), INARR(*)
      DOUBLE PRECISION FWORK(*), Y0(M),X0(M,N), REARR(*)
      INTEGER IFLAG
      EXTERNAL GETF, GETDF

CCCCCCC LOCAL
      DOUBLE PRECISION YDNORM, TOLM, TOLMIN, TOLMAX
      PARAMETER (TOLMIN=1.0D-13, TOLMAX=1.0D-2)
      INTEGER INTDIR, IFIRST, IWHAT
      INTEGER I, J

CCCCCCC POINTERS, AND COMMON VARIABLES
      LOGICAL LAST
      DOUBLE PRECISION HMIN, HNEXT, TFIRST
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      INTEGER IAPTR, IXPTR, IXHPTR, IQPTR, IBPTR
      INTEGER K1PTR, K2PTR, K3PTR, K4PTR, K5PTR, K6PTR
      INTEGER ISTAGE
      INTEGER IYPTR, IYHTPT
      INTEGER K1YPTR, K2YPTR, K3YPTR, K4YPTR, K5YPTR, K6YPTR
      INTEGER I1YPTR, IS1YPT, IS2YPT, IS3YPT, IS4YPT

      COMMON /STEP/ HMIN, HNEXT, TFIRST, INTDIR, IFIRST, IWHAT, LAST
      SAVE /STEP/
      COMMON /NSQUAX/ IAPTR, IXPTR, IQPTR, IXHPTR, IBPTR
      SAVE /NSQUAX/
      COMMON /RKPTRS/ K1PTR, K2PTR, K3PTR, K4PTR, K5PTR, K6PTR,
     *                ISTAGE
      SAVE /RKPTRS/
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/
      COMMON /YRKPTR/ IYPTR, IYHTPT,  
     *                K1YPTR, K2YPTR, K3YPTR, K4YPTR, K5YPTR, K6YPTR, 
     *                I1YPTR, IS1YPT, IS2YPT, IS3YPT, IS4YPT
      SAVE /YRKPTR/


C CHECK CORRECTNESS OF INPUT ON T0 AND TE
      IF (T0.EQ.TE) THEN
         IFLAG=4
         RETURN
        ELSE
C CHECK CONSISTENCY OF INTEGRATION DIRECTION 
         IF((TE-T0)*INTDIR.LT.ZERO) THEN
            IFLAG=6
            RETURN
         END IF
      END IF
C CHECK IPAR(1) AND IPAR(2) AND SET DEFAULTS
      IF (IPAR(1).NE.1) IPAR(1)=0
      IF (IPAR(2).NE.1) IPAR(2)=0
C CHECK CONSISTENCY OF INPUT IF FIXED STEPSIZE
      IF((IPAR(1).EQ.1) .AND. (INTDIR*DT.LT.ZERO)) THEN
         IFLAG=6
         RETURN
      ENDIF
      IF ((IPAR(1).EQ.1).AND.(DABS(DT).LE.HMIN)) THEN
         IFLAG=5
         RETURN
      END IF
C ZERO COUNTERS OF STEPS AND REJECTIONS
      IPAR(12)=0
      IPAR(13)=0
C CHECK IPAR(6) 
      IF (IPAR(6).NE.1) IPAR(6)=0 
      IWHAT=IPAR(6)
C CHECK IPAR(9) AND IPAR(10)
      IF ((IPAR(10).NE.1).AND.(IPAR(10).NE.2).AND.
     *   (IPAR(10).NE.10).AND.(IPAR(10).NE.20).AND.
     *   (IPAR(10).NE.21).AND.(IPAR(10).NE.210)) IPAR(10)=0
      IF (IPAR(8).LT.4) THEN
         IF (IPAR(9).NE.1) IPAR(9)=0      
      END IF
      IF ((IPAR(8).GE.4).AND.((IPAR(10).EQ.2).OR.(IPAR(10).EQ.20)
     *    .OR.(IPAR(10).EQ.21).OR.(IPAR(10).EQ.210)).AND.
     *    (IPAR(1).EQ.0)) THEN
         IFLAG=13
         RETURN
      END IF
      IF (IPAR(9).EQ.1.AND.IPAR(8).GT.3.AND.IWHAT.NE.1) THEN
         IFLAG=13
         RETURN
      END IF
      IF ((IPAR(1).EQ.0).AND.(IPAR(9).EQ.1).AND.(IWHAT.NE.1)
     *     .AND.((IPAR(10).EQ.1).OR.(IPAR(10).EQ.10).OR.
     *     (IPAR(10).EQ.21).OR.(IPAR(10).EQ.210))) THEN
         IFLAG=12
         RETURN
      ENDIF
C CHECK TOLERANCES TOLT, TOLQ, TOLL .  THE VALUES OF TOLMIN AND TOLMAX 
C (SEE PARAMETER STATEMENT) MAY BE CHANGED, BUT THIS IS DISCOURAGED
      IF (IPAR(1).EQ.0) THEN
         IF (IPAR(10).EQ.0.OR.IPAR(10).EQ.10.OR.
     *       IPAR(10).EQ.20.OR.IPAR(10).EQ.210) THEN
            IF ((TOLT.LT.TOLMIN).OR.(TOLT.GT.TOLMAX)) THEN
               IFLAG=3
               RETURN
            END IF
         ENDIF
         IF (IPAR(10).EQ.2.OR.IPAR(10).EQ.20.OR.
     *       IPAR(10).EQ.21.OR.IPAR(10).EQ.210) THEN
            IF ((TOLQ.LT.TOLMIN).OR.(TOLQ.GT.TOLMAX)) THEN
               IFLAG=3
               RETURN
            END IF
         ENDIF
         IF (IPAR(10).EQ.1.OR.IPAR(10).EQ.10.OR.
     *       IPAR(10).EQ.21.OR.IPAR(10).EQ.210) THEN
            DO 2 I=1,N
               IF (TOLL(I).LT.TOLMIN.OR.TOLL(I).GT.TOLMAX) THEN
                  IFLAG=3
                  RETURN
               END IF
 2          CONTINUE
         ENDIF
      END IF
C IF ONLY WANT TRAJECTORY, SKIP INTIALIZATION FOR LES STUFF
      IF (IWHAT.EQ.1) GOTO 8 
C MONITOR FROM WHERE NEED TO APPROXIMATE LYAPUNOV EXPONENTS
      TFIRST=T0
C FIRST TIME, INITIALIZE LES FROM QR OF X0, AND SAVE Q FACTOR
C Q-FACTOR OF X0 IS IN X0 ITSELF
      IF (IPAR(4).EQ.1) THEN
         CALL MODGRS(M,N,X0,M,FWORK(IXHPTR),IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG = 8
            RETURN
         END IF
        ELSE
         DO 5 J=1,N
            DO 4 I=1,M
               X0(I,J)=ZERO
 4          CONTINUE
            X0(J,J)=ONE
 5       CONTINUE
      END IF
      DO 7 I=1,N
         APPLES(I)=ZERO
 7    CONTINUE
 8    CONTINUE
      IFIRST=1
      IPAR(3)=0
      IF (IPAR(1).EQ.0) THEN
C CHECK IF ERROR CONTROL IS COMPATIBLE WITH WHAT USER WANTS
         IF ((IWHAT.EQ.1).AND.(IPAR(10).NE.0)) THEN
            IFLAG=11
            RETURN
         END IF
CCCCCCC PREPARE INITIAL STEP SIZE IN VARIABLE STEP MODE
C TO CHOOSE INITIAL STEP, 
C IF IPAR(10)=0,10,20,210 USE Y-DOT (DERIVATIVE OF VECTOR FIELD)
C   OTHERWISE USE Q-DOT IF IPAR(8)=0,2, OR
C   USE X-DOT WHEN IPAR(8)=1,3,4,5 
         IF ((IPAR(10).EQ.0).OR.(IPAR(10).EQ.10).OR.
     1       (IPAR(10).EQ.20).OR.(IPAR(10).EQ.210)) THEN 
            CALL GETF(M,Y0,FWORK(K1YPTR),INARR,REARR)
C HAVE Y-DOT IN FWORK(K1YPTR).  FIND ITS MAX-NORM
            YDNORM=ZERO
            DO 9 I=1,M
 9            YDNORM=DMAX1(YDNORM, DABS(FWORK(K1YPTR+I-1)))
           ELSE
            IF (IPAR(8).EQ.0.OR.IPAR(8).EQ.2) THEN 
C QDOT RETURNS QDOT IN RK1
               CALL QDOT(GETDF,M,N,Y0,FWORK(IAPTR),FWORK(IBPTR),
     *             X0,FWORK(K1PTR),INARR,REARR)
              ELSE
C XDOT RETURNS XDOT IN RK1
               CALL XDOT(GETDF,M,N,Y0,FWORK(IAPTR),
     1                X0,FWORK(K1PTR),INARR,REARR)
            ENDIF
C HAVE X-DOT OR Q-DOT IN FWORK(K1PTR).  FIND ITS MAX-NORM
            YDNORM=ZERO
            DO 10 J=1,N
               DO 10 I=1,M
 10               YDNORM=DMAX1(YDNORM,DABS(FWORK(K1PTR+(I-1)*N+J-1)))
         END IF
C TAKE DT=(5-TH ROOT OF TOL) FOR RKDP OR
C      DT=(4-TH ROOT OF TOL) FOR RK38
C WEIGHTED AGAINST DERIVATIVE OF Y, OR OF Q, OR OF X
         IF ((IPAR(10).EQ.0).OR.(IPAR(10).EQ.10).OR.
     1       (IPAR(10).EQ.20).OR.(IPAR(10).EQ.210)) THEN 
            TOLM=TOLT
           ELSEIF (IPAR(10).EQ.2) THEN
            TOLM=TOLQ
           ELSE
            TOLM=TOLL(1)
            DO 11 I=2,N
 11            TOLM=DMIN1(TOLM,TOLL(I))           
         ENDIF
         IF (IPAR(8).EQ.0.OR.IPAR(8).EQ.1.OR.IPAR(8).EQ.4) THEN 
            DT=INTDIR*TOLM**EXP1/(DMAX1(ONE,YDNORM))
           ELSE
            DT=INTDIR*TOLM**EXP2/(DMAX1(ONE,YDNORM))
         END IF
         IF (INTDIR*(T0+DT).GT.INTDIR*TE) DT=TE-T0
         IF (DABS(DT).LE.HMIN) THEN
            IFLAG=5
            RETURN
         END IF
      END IF

      RETURN
      END
CCCCC END OF NEWYES CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: CFRK
C
C*****Purpose
C
C forms the coefficients for RKDP or RK38
C
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCCC CFRK CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CFRK(IFORM,ARK,CRK,BHAT)

      IMPLICIT NONE

      INTEGER IFORM
      DOUBLE PRECISION ARK(7,7), CRK(7), BHAT(7)
      INTEGER I, J
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/
       
C SET UP THE RKDP COEFFICIENTS

      DO 2 I=1,7
         BHAT(I)=ZERO
         CRK(I)=ZERO
         DO 2 J=1,7
            ARK(J,I)=ZERO
 2    CONTINUE
      IF ((IFORM.EQ.0).OR.(IFORM.EQ.1).OR.(IFORM.EQ.4)) THEN
         ARK(2,1)=0.2D0
         CRK(2)=ARK(2,1)
         ARK(3,1)=0.075D0
         ARK(3,2)=0.225D0
         CRK(3)=ARK(3,1)+ARK(3,2)
         ARK(4,1)=0.9777777777777777777D0
         ARK(4,2)=-3.733333333333333333D0
         ARK(4,3)=3.5555555555555555555D0
         CRK(4)=0.8D0
         ARK(5,1)=19372.D0/6561.D0
         ARK(5,2)=-25360.D0/2187.D0
         ARK(5,3)=64448.D0/6561.D0
         ARK(5,4)=-212.D0/729.D0
         CRK(5)=ARK(5,1)+ARK(5,2)+ARK(5,3)+ARK(5,4)
         ARK(6,1)=9017.D0/3168.D0
         ARK(6,2)=-355.D0/33.D0
         ARK(6,3)=46732.D0/5247.D0
         ARK(6,4)=49.D0/176.D0
         ARK(6,5)=-5103.D0/18656.D0
         CRK(6)=ARK(6,1)+ARK(6,2)+ARK(6,3)+ARK(6,4)+ARK(6,5)
         ARK(7,1)=0.0911458333333333D0
         ARK(7,3)=0.449236298292902D0
         ARK(7,4)=0.651041666666667D0
         ARK(7,5)=-0.322376179245283D0
         ARK(7,6)=0.130952380952381D0
         CRK(7)=ONE
         BHAT(1)=0.0899131944444444D0 
         BHAT(3)=0.4534890685834082D0 
         BHAT(4)=0.6140625D0 
         BHAT(5)=-0.2715123820754717D0 
         BHAT(6)=0.08904761904761904D0 
         BHAT(7)=0.025D0 
       ELSE 
         ARK(2,1)=ONE/3.D0
         CRK(2)=ARK(2,1)
         ARK(3,1)=-ONE/3.D0
         ARK(3,2)=ONE
         CRK(3)=ARK(3,1)+ARK(3,2)
         ARK(4,1)=ONE
         ARK(4,2)=-ONE
         ARK(4,3)=ONE
         CRK(4)=ONE
         ARK(5,1)=0.125D0
         ARK(5,2)=0.375D0
         ARK(5,3)=0.375D0
         ARK(5,4)=0.125D0
         CRK(5)=ONE
         BHAT(1)=0.08333333333333D0 
         BHAT(2)=0.5D0 
         BHAT(3)=0.25D0 
         BHAT(5)=0.16666666666667D0 
      ENDIF

      RETURN
      END
CCCCC END OF CFRK CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: XDOT
C
C*****Purpose
C
C forms the derivative of the X-equation 
C (this is also used to select the first stepsize)
C
C input: 
C     GETDF, declared EXTERNAL, returns coefficient matrix
C     N, dimension of the coefficient matrix
C     NQ, number of columns of X (or Q) (i.e., to be triangularized)
C     Y0, value where to get derivative (Jacobian)
C     A, to store the Jacobian
C     XMAT, the value of X for which we form A*X 
C          (at first, this is X0, i.e., Q0)
C output:
C     XP, the first derivative of X
C
C BLAS: Blas3: DGEMM
C
C*****Authors: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC XDOT CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE XDOT(GETDF,N,NQ,Y0,A,XMAT,XP,INARR,REARR)
      IMPLICIT NONE

      EXTERNAL GETDF
      INTEGER N, NQ, INARR(*)
      DOUBLE PRECISION Y0(N), A(N,N), XMAT(N,NQ), XP(N,NQ), REARR(*)

      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

      CALL GETDF(N,Y0,A,INARR,REARR)
C FORM A*X0 WITH CALL TO BLAS
      CALL DGEMM ('N', 'N', N, NQ, N, ONE, A, N, XMAT, N,
     $                   ZERO, XP, N )

      RETURN
      END
CCCCC END OF XDOT CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: QDOT
C
C*****Purpose
C
C forms the derivative of the Q-equation
C
C input: 
C     GETDF, declared EXTERNAL, returns Jacobian matrix
C     N, dimension of the coefficient matrix
C     NQ, number of columns of Q (i.e., to be triangularized)
C     Y0, value where to get derivative (Jacobian)
C     A, to store the Jacobian
C     B, workspace
C     QMAT, the matrix Q
C output:
C     A, the coefficient matrix at T0
C     B, on diagonal it has diag(Q^TAQ)
C     QP, the first derivative of Q
C
C*****Authors: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC QDOT CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE QDOT(GETDF,N,NQ,Y0,A,B,QMAT,QP,INARR,REARR)
      IMPLICIT NONE

      EXTERNAL GETDF
      INTEGER N, NQ, INARR(*)
      DOUBLE PRECISION Y0(N), A(N,N), B(N,NQ), QMAT(N,NQ), QP(N,NQ)
      DOUBLE PRECISION REARR(*)
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

      INTEGER I, J

      CALL GETDF(N,Y0,A,INARR,REARR)
C FORM A*Q
      CALL DGEMM ('N', 'N', N, NQ, N, ONE, A, N, QMAT, N,
     $                   ZERO, QP, N )
C FORM Q^T*A*Q-S (S IS SKEW MATRIX)
      CALL DGEMM ('T', 'N', NQ, NQ, N, ONE, QMAT, N, QP, N,
     $                   ZERO, B, N )
C FINISH OFF Q^T*A*Q-S
      DO 11 I=1,NQ
         DO 11 J=1,NQ
            IF (I.LT.J) THEN 
               B(I,J)=B(I,J)+B(J,I)
              ELSE IF (I.GT.J) THEN 
               B(I,J)=0.0D0
            END IF
 11   CONTINUE
C UPDATE QDOT
      CALL DGEMM ('N', 'N', N, NQ, NQ, -1.0D0*ONE, QMAT, N, B, N,
     $                   ONE, QP, N )

      RETURN
      END
CCCCC END OF QDOT CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: FORMAX
C
C*****Purpose
C
C forms the matrix product AX
C 
C input: 
C     LDA, integer (leading dimension of A)
C     A, double precision (coefficient matrix; it is (M,M))
C     m, integer (true number of rows/columns of A)
C     LDX, integer (leading dimension of X)
C     X, double precision (the matrix X; it is (M,N))
C     N, integer (number of columns of X)
C     LDR, integer (leading dimension of RHS)
C output:
C     RHS, double precision (LDR,N) matrix (holds A*X in (M,N) corner)
C
C BLAS: Blas3: DGEMM
C
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC FORMAX CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FORMAX(LDA,A,M,LDX,X,N,LDR,RHS)

      IMPLICIT NONE

      INTEGER LDA,M,LDX,N,LDR
      DOUBLE PRECISION A(LDA,M), X(LDX,N), RHS(LDR,N)
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

C FORM A*X IN RHS
      CALL DGEMM ('N', 'N', M, N, M, ONE, A, LDA, X, LDX,
     $                   ZERO, RHS, LDR )

      RETURN
      END
CCCCC END OF FORMAX CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: PQDPNS
C
C*****Purpose
C
C performs the one step integration of y'=f(y)
C and (if required, see IPAR(6)), of Q'=F(A,Q), Q0 given.
C Monitors errors on the one step approximations to the 
C trajectory, and (if required) on the LEs and/or Q.
C It updates stepsize and approximation to the solutions
C trajectory and (if required) to the LEs and Q.
C
C       Method used: Complete Projection continuous QR method with 
C          DP5 as integrator for Q(i.e., also
C          for stage values). As specified by
C          IPAR(9), two options are provided for 
C          approximating the LEs:
C          (i) a quadrature rule for the nu-variables of the 
C              same order as DP5 (the default), or
C          (ii) the composite trap-rule.  
C
C          Code can proceed in fixed or variable
C          step size.  In the latter case, error control is 
C          performed on the trajectory, and/or on the Lyapunov 
C          exponents, and/or on the Q-factor.
C
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC PQDPNS CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PQDPNS(GETF,GETDF,M,N,T,H,HOLD,Y0,X0,RLYAP,TOLT,
     *          TOLQ,TOLL,AMAT,Q,QHAT,BMAT,RK1,RK2,RK3,RK4,RK5,  
     *          RK6,STAGE,RNU,RNUHAT,RKNU1,RKNU2,RKNU3,RKNU4,
     *          RKNU5,RKNU6,ARK,BHAT,Y,YHAT,YRK1,YRK2,YRK3,YRK4,
     *          YRK5,YRK6,Y1,Y15,Y03,Y45,Y89,
     *          IPAR,NREJ,IFLAG,INARR,REARR)

      IMPLICIT NONE

      EXTERNAL GETF, GETDF
      INTEGER M,N,IPAR(*),NREJ,IFLAG, INARR(*)
      DOUBLE PRECISION T,H,HOLD, Y0(M), REARR(*)
      DOUBLE PRECISION X0(M,N), RLYAP(N), TOLT,TOLQ,TOLL(*)
      DOUBLE PRECISION AMAT(M,M), Q(M,N), QHAT(M,N)
      DOUBLE PRECISION BMAT(M,N)
      DOUBLE PRECISION RK1(M,N),RK2(M,N)
      DOUBLE PRECISION RK3(M,N),RK4(M,N)
      DOUBLE PRECISION RK5(M,N),RK6(M,N)
      DOUBLE PRECISION STAGE(M,N)     
      DOUBLE PRECISION RNU(N), RNUHAT(N)       
      DOUBLE PRECISION RKNU1(N), RKNU2(N), RKNU3(N), RKNU4(N)
      DOUBLE PRECISION RKNU5(N), RKNU6(N)
      DOUBLE PRECISION ARK(7,7), BHAT(7)
      DOUBLE PRECISION Y(M),YHAT(M),YRK1(M),YRK2(M),YRK3(M),YRK4(M)
      DOUBLE PRECISION YRK5(M),YRK6(M)
      DOUBLE PRECISION Y1(M),Y15(M),Y03(M),Y45(M),Y89(M)

      LOGICAL PREVIO, LAST
      INTEGER INTDIR, IFIRST, IFXDPF, IERRC, IHOWLE, IWHAT, I, J
      DOUBLE PRECISION HNEW, HMIN, HNEXT, TFIRST
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      DOUBLE PRECISION SCALE, RTOL, ERR
      COMMON /STEP/ HMIN, HNEXT, TFIRST, INTDIR, IFIRST, IWHAT, LAST
      SAVE /STEP/
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

      IFXDPF=IPAR(1)
      IERRC=IPAR(10)
      IHOWLE=IPAR(9)
      PREVIO=.FALSE.
C PREVIO IS .FALSE. IF PREVIOUS STEP WAS SUCCESSFUL
C IT IS .TRUE. IF PREVIOUS STEP WAS A FAILURE
 2    CONTINUE
CCCCCCCCCC INTEGRATE FOR X: RETURNS X AND XHAT
      CALL PQRKDP(GETF,GETDF,PREVIO,M,N,IFIRST,IWHAT,IFXDPF,
     *        IHOWLE,IERRC,H,Y0,X0,
     *        AMAT,Q,QHAT,BMAT,RK1,RK2,RK3,RK4,RK5,RK6,
     *        STAGE, RNU,RNUHAT,RKNU1,RKNU2,RKNU3,RKNU4,
     *        RKNU5,RKNU6,ARK,BHAT,
     *        Y,YHAT,YRK1,YRK2,YRK3,YRK4,YRK5,YRK6,
     *        Y1,Y15,Y03,Y45,Y89,IFLAG,INARR,REARR)
      IF (IFIRST.EQ.1) IFIRST=0
      IF (IFLAG.NE.0) THEN
         IF (IFXDPF.EQ.0) THEN
C TRY DRASTIC WAY TO RECOVER
            IF (DABS(H).GT.HMIN) THEN
               IFLAG=0
               H=INTDIR*HMIN
               PREVIO=.TRUE.
               NREJ=NREJ+1
               GOTO 2
              ELSE
               IFLAG=5
               RETURN
            ENDIF
           ELSE
            IFLAG=10
            RETURN
         END IF
      END IF
      IF (IFXDPF.EQ.0) THEN
C HERE WE PREDICT NEXT STEPSIZE HNEW, NEVER ALLOWING IT TO
C BE GREATER THAN FIVE TIMES CURRENT H 
C (THIS VALUE CAN BE CHANGED IN HINCR IN "INIT")
         HNEW=HINCR*H
         RTOL=ZERO
         IF (IERRC.EQ.0) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
           ELSE IF (IERRC.EQ.1) THEN
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
           ELSE IF (IERRC.EQ.2) THEN
            CALL NERRQ(M,N,Q,QHAT,X0,ERR,TOLQ,ONE,ZERO)
           ELSE IF (IERRC.EQ.10) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
           ELSE IF (IERRC.EQ.20) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRQ(M,N,Q,QHAT,X0,ERR,TOLQ,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
           ELSE IF (IERRC.EQ.21) THEN
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRQ(M,N,Q,QHAT,X0,ERR,TOLQ,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
           ELSE
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRQ(M,N,Q,QHAT,X0,ERR,TOLQ,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
         END IF
C PREDICT NEW STEPSIZE
         IF (ERR.GT.ZERO)         
     1      HNEW=DMIN1(DABS(HNEW),DABS(HSFTY*H*(ONE/ERR)**EXP1))
C GO BACK IF FAILED UNLESS H TOO SMALL
         IF (ERR.GT.ONE) THEN
C AVOID DRASTIC REDUCTION OF STEPSIZE (NEVER LESS THAN 1/5)
C (THIS VALUE CAN BE CHANGED IN HDRSTC IN "INIT")
            H=INTDIR*DMAX1(DABS(HNEW),HDRSTC*DABS(H))
            NREJ=NREJ+1
            IF (INTDIR*H.LE.HMIN) THEN
               PRINT*,' H = ', H
               IFLAG = 5
               RETURN
            ENDIF
            PREVIO=.TRUE.
            GOTO 2
         ENDIF            
      END IF
C HAVE COMPLETED THE STEP WITH STEPSIZE H
C SAVE LAST SUCCESSFUL STEPSIZE
      HOLD=H
      IF (IFXDPF.EQ.0) THEN
C         HNEW=INTDIR*DMAX1(DABS(HNEW),DABS(H))
         IF((PREVIO).AND.(DABS(HNEW).GT.DABS(H)))HNEW=H
         H=HNEW
      END IF
C FINALLY, UPDATE THE APPROXIMATION FOR THE LES, Q, AND TRAJECTORY
      IF (IWHAT.EQ.1) GOTO 29
      SCALE=ONE/(T+HOLD-TFIRST)      
      DO 25 I=1,N
         RLYAP(I)=((T-TFIRST)*SCALE)*RLYAP(I)+SCALE*RNU(I)
C LOWER ORDER APPROXIMATION
C         RLYAP(I)=((T-TFIRST)*SCALE)*RLYAP(I)+SCALE*RNUHAT(I)
 25   CONTINUE         
      DO 28 J=1,N
         DO 28 I=1,M
            X0(I,J)=Q(I,J)
 28   CONTINUE
 29   CONTINUE
      DO 32 I=1,M
         Y0(I)=Y(I)
 32   CONTINUE

      RETURN
      END
CCCCC END OF PQDPNS CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: PQRKDP
C
C*****Purpose
C
C integrates on one step the solution of x'=f(x).  If needed,
C also of Q'=F(Df,Q).  It returns two approximations, x and xhat,
C and Q, the Q-factor in the QR factorization of X.
C It also returns two approximations to the updates of
C the LEs, rnu and rnuhat. All of these can be used for error purposes.
C
C Integrator is the DP rule with embedded FSAL of order 5 (4)
C Scheme is total projection method
C 
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC PQRKDP CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PQRKDP(GETF,GETDF,PREVIO,M,N,IFIRST,IWHAT,IFXDPF,
     *        IHOWLE,IERRC,H,Y0,X0,
     *        AMAT,Q,QHAT,BMAT,RK1,RK2,RK3,RK4,RK5,RK6,
     *        STAGE, RNU,RNUHAT,RKNU1,RKNU2,RKNU3,RKNU4,
     *        RKNU5,RKNU6,ARK,BHAT,
     *        Y,YHAT,YRK1,YRK2,YRK3,YRK4,YRK5,YRK6,
     *        Y1,Y15,Y03,Y45,Y89,IFLAG,INARR,REARR)

      IMPLICIT NONE

      EXTERNAL GETF,GETDF
      LOGICAL PREVIO
      INTEGER M, N, IFIRST, IWHAT, IFXDPF, IHOWLE
      INTEGER IERRC, IFLAG,INARR(*)
      DOUBLE PRECISION H, Y0(M), X0(M,N), REARR(*)
      DOUBLE PRECISION AMAT(M,M),Q(M,N), QHAT(M,N)
      DOUBLE PRECISION BMAT(M,N)
      DOUBLE PRECISION RK1(M,N),RK2(M,N)
      DOUBLE PRECISION RK3(M,N),RK4(M,N)
      DOUBLE PRECISION RK5(M,N),RK6(M,N)
      DOUBLE PRECISION STAGE(M,N)
      DOUBLE PRECISION RNU(N),RNUHAT(N)
      DOUBLE PRECISION RKNU1(N),RKNU2(N),RKNU3(N)
      DOUBLE PRECISION RKNU4(N),RKNU5(N),RKNU6(N)
      DOUBLE PRECISION ARK(7,7), BHAT(7)
      DOUBLE PRECISION Y(M),YHAT(M),YRK1(M),YRK2(M),YRK3(M),YRK4(M)
      DOUBLE PRECISION YRK5(M),YRK6(M)
      DOUBLE PRECISION Y1(M),Y15(M),Y03(M),Y45(M),Y89(M)
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      INTEGER I, J
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

C INITIALIZE Y TO OLD Y0
      DO 2 I=1,M
         Y(I)=Y0(I)
 2       YHAT(I)=Y0(I)
C FORM RK1, RK2, RK3, RK4, RK5, RK6, RK7 
C FORM ALSO THE STAGE VALUES AT H/5, 3H/10, 4H/5, 8H/9 AND H 
C (THE LATTER IS NOT THE SOLUTION VALUE, BUT JUST THE STAGE VALUE)
C ... RK1 AND U15 ...
C VERY FIRST CALL
C      IF (IFIRST.EQ.1) THEN
      IF ((IFIRST.EQ.1).OR.(IFXDPF.EQ.1)) THEN
         CALL GETF(M,Y,YRK1,INARR,REARR) 
         GOTO 7
      END IF
C IF NOT FAILED PREVIOUS STEP, YRK IS IN YRK2
C OTHERWISE, IT IS IN YRK1
      IF (.NOT.PREVIO) THEN
          DO 5 I=1,M 
 5           YRK1(I)=YRK2(I)
      ENDIF
 7    CONTINUE
      DO 11 I=1,M
         Y15(I) = Y(I) + H*ARK(2,1)*YRK1(I)
 11   CONTINUE
C ... RK2 AND U03 ...
      CALL GETF(M,Y15,YRK2,INARR,REARR) 
      DO 14 I=1,M
         Y03(I) = Y(I) + H*(ARK(3,1)*YRK1(I)+ARK(3,2)*YRK2(I))
 14   CONTINUE
C ... RK3 AND U45 
      CALL GETF(M,Y03,YRK3,INARR,REARR) 
      DO 16 I=1,M
         Y45(I) = Y(I) + H*(ARK(4,1)*YRK1(I)+ARK(4,2)*YRK2(I)+
     1                ARK(4,3)*YRK3(I))
 16   CONTINUE
C ... RK4 AND U89 
      CALL GETF(M,Y45,YRK4,INARR,REARR) 
      DO 19 I=1,M
         Y89(I) = Y(I) + H*(ARK(5,1)*YRK1(I)+ARK(5,2)*YRK2(I)+
     1                ARK(5,3)*YRK3(I)+ARK(5,4)*YRK4(I))
 19   CONTINUE
C ... RK5 AND U1 (STAGE VALUE AT H) ...
      CALL GETF(M,Y89,YRK5,INARR,REARR) 
      DO 20 I=1,M
         Y1(I) = Y(I) + H*(ARK(6,1)*YRK1(I)+ARK(6,2)*YRK2(I)+
     1                ARK(6,3)*YRK3(I)+ARK(6,4)*YRK4(I)+
     2                ARK(6,5)*YRK5(I))
 20   CONTINUE
C ... RK6, SOLUTION AT H ...
      CALL GETF(M,Y1,YRK6,INARR,REARR) 
      DO 23 I=1,M
         Y(I) = Y(I) + H*(ARK(7,1)*YRK1(I)+
     1         ARK(7,3)*YRK3(I)+ARK(7,4)*YRK4(I)+
     2         ARK(7,5)*YRK5(I)+ARK(7,6)*YRK6(I))
 23   CONTINUE
      IF (IFXDPF.EQ.1) GOTO 25
C ... RK7 (STORED IN RK2), AND COMPARISON SOLUTION XHAT ...
C TO BE USED FOR ERROR MONITOR/STEP-SIZE SELECTION
      CALL GETF(M,Y,YRK2,INARR,REARR) 
      DO 24 I=1,M 
         YHAT(I)=YHAT(I)+H*(BHAT(1)*YRK1(I)+BHAT(3)*YRK3(I)+
     1       BHAT(4)*YRK4(I)+BHAT(5)*YRK5(I)+BHAT(6)*YRK6(I)+
     2       BHAT(7)*YRK2(I))
 24   CONTINUE
C NOW DO IT FOR THE TRANSITION MATRIX (IF IWHAT=0)
 25   IF (IWHAT.EQ.1) RETURN
C INTIALIZE Q TO OLD Q
C THE RNU'S ARE INITIALLY AT 0, SINCE WE TAKE ONE STEP OF INTEGRAL
      DO 26 J=1,N
         RNU(J)=ZERO
         RNUHAT(J)=ZERO
         DO 26 I=1,M
            Q(I,J)=X0(I,J)
            QHAT(I,J)=X0(I,J)
 26    CONTINUE
C FORM RK1, RK2, RK3, RK4, RK5, RK6, RK7 
C FORM ALSO THE STAGE VALUES AT H/5, 3H/10, 4H/5, 8H/9 AND H (NB THE LATTER IS
C NOT THE SOLUTION VALUE, BUT JUST THE STAGE VALUE)
C ... RK1 AND U15 ...
C VERY FIRST CALL NEED AMAT AND RK1
C IF FAILED PREVIOUS STEP, HAVE RK1 IN RK1
C IF NOT FAILED PREVIOUS STEP, NEED RK1.  
CC **** SEE IF CAN ECONOMIZE HERE ****
      CALL QDOT(GETDF,M,N,Y0,AMAT,BMAT,Q,RK1,INARR,REARR)
      DO 29 J=1,N
         DO 29 I=1,M
            STAGE(I,J) = Q(I,J) + H*ARK(2,1)*RK1(I,J)
 29   CONTINUE
      DO 30 I=1,N
 30      RKNU1(I)=BMAT(I,I)
      IF (IHOWLE.EQ.0) THEN
         CALL MODGRS(M,N,STAGE,M,BMAT,IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG=10
            RETURN
         ENDIF
      ENDIF
C ... RK2 AND U03 ...
      CALL QDOT(GETDF,M,N,Y15,AMAT,BMAT,STAGE,RK2,INARR,REARR)
      DO 34 J=1,N
         DO 34 I=1,M
            STAGE(I,J) = Q(I,J) + H*(ARK(3,1)*RK1(I,J)+
     1                 ARK(3,2)*RK2(I,J))
 34   CONTINUE
      IF (IHOWLE.EQ.0) THEN
         DO 40 I=1,N
 40         RKNU2(I)=BMAT(I,I)
         CALL MODGRS(M,N,STAGE,M,BMAT,IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG=10
            RETURN
         ENDIF
      ENDIF
C ... RK3 AND U45 ...
      CALL QDOT(GETDF,M,N,Y03,AMAT,BMAT,STAGE,RK3,INARR,REARR)
      DO 43 J=1,N
         DO 43 I=1,M
            STAGE(I,J) = Q(I,J) + H*(ARK(4,1)*RK1(I,J)+
     1                ARK(4,2)*RK2(I,J)+ARK(4,3)*RK3(I,J))
 43   CONTINUE
      IF (IHOWLE.EQ.0) THEN
         DO 50 I=1,N
 50         RKNU3(I)=BMAT(I,I)
         CALL MODGRS(M,N,STAGE,M,BMAT,IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG=10
            RETURN
         ENDIF
      ENDIF
C ... RK4 AND U89 
      CALL QDOT(GETDF,M,N,Y45,AMAT,BMAT,STAGE,RK4,INARR,REARR)
       DO 55 J=1,N
         DO 55 I=1,M
            STAGE(I,J) = Q(I,J) + H*(ARK(5,1)*RK1(I,J)+
     1                 ARK(5,2)*RK2(I,J)+
     2                 ARK(5,3)*RK3(I,J)+ARK(5,4)*RK4(I,J))
 55   CONTINUE
      IF (IHOWLE.EQ.0) THEN
         DO 60 I=1,N
 60         RKNU4(I)=BMAT(I,I)
         CALL MODGRS(M,N,STAGE,M,BMAT,IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG=10
            RETURN
         ENDIF
      ENDIF
C ... RK5 AND U1 (STAGE VALUE AT H) ...
      CALL QDOT(GETDF,M,N,Y89,AMAT,BMAT,STAGE,RK5,INARR,REARR)
      DO 65 J=1,N
         DO 65 I=1,M
            STAGE(I,J) = Q(I,J) + H*(ARK(6,1)*RK1(I,J)+
     1                ARK(6,2)*RK2(I,J)+ARK(6,3)*RK3(I,J)+
     2                ARK(6,4)*RK4(I,J)+ARK(6,5)*RK5(I,J))
 65   CONTINUE
      IF (IHOWLE.EQ.0) THEN
         DO 70 I=1,N
 70         RKNU5(I)=BMAT(I,I)
         CALL MODGRS(M,N,STAGE,M,BMAT,IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG=10
            RETURN
         ENDIF
      ENDIF
C ... RK6, SOLUTION AT H ...
      CALL QDOT(GETDF,M,N,Y1,AMAT,BMAT,STAGE,RK6,INARR,REARR)
      IF (IHOWLE.EQ.0) THEN
         DO 77 I=1,N
 77         RKNU6(I)=BMAT(I,I)
      ENDIF
      DO 80 J=1,N
         IF (IHOWLE.EQ.0) THEN 
            RNU(J) = RNU(J) + H*(ARK(7,1)*RKNU1(J)+
     1         ARK(7,3)*RKNU3(J)+ARK(7,4)*RKNU4(J)+
     2         ARK(7,5)*RKNU5(J)+ARK(7,6)*RKNU6(J))
         ENDIF
         DO 80 I=1,M
            Q(I,J) = Q(I,J) + H*(ARK(7,1)*RK1(I,J)+
     1         ARK(7,3)*RK3(I,J)+ARK(7,4)*RK4(I,J)+
     2         ARK(7,5)*RK5(I,J)+ARK(7,6)*RK6(I,J))
 80   CONTINUE
C ... RK7 (STORED IN RK2), AND COMPARISON SOLUTION QHAT AND RNUHAT ...
C TO BE USED FOR ERROR MONITOR/STEP-SIZE SELECTION
      CALL MODGRS(M,N,Q,M,BMAT,IFLAG) 
      IF (IFLAG.NE.0) THEN
         IFLAG=10
         RETURN
      ENDIF
      IF (((IERRC.EQ.2).OR.(IERRC.EQ.20).OR.(IERRC.EQ.21).OR.
     *    (IERRC.EQ.210)).AND.IFXDPF.EQ.0) THEN
         CALL QDOT(GETDF,M,N,Y,AMAT,BMAT,Q,RK2,INARR,REARR)
        ELSE
         CALL GETDF(M,Y,AMAT,INARR,REARR)
         CALL DGQTAQ(M,M,N,Q,AMAT,M,BMAT,M)
      ENDIF
      DO 83 I=1,N
 83      RKNU2(I)=BMAT(I,I)
      DO 90 J=1,N 
         IF (IHOWLE.EQ.0.AND.IFXDPF.EQ.0.AND.(IERRC.EQ.1.OR.
     *       IERRC.EQ.10.OR.IERRC.EQ.21.OR.IERRC.EQ.210)) THEN
            RNUHAT(J)=RNUHAT(J)+H*(BHAT(1)*RKNU1(J)+
     1       BHAT(3)*RKNU3(J)+BHAT(4)*RKNU4(J)+
     2       BHAT(5)*RKNU5(J)+BHAT(6)*RKNU6(J)+BHAT(7)*RKNU2(J))
           ELSEIF (IHOWLE.EQ.1) THEN
             CALL TRAPRL(RKNU1(J),RKNU2(J),H,RNU(J))
         ENDIF
         IF (IFXDPF.EQ.0.AND.(IERRC.EQ.2.OR.IERRC.EQ.20.OR.
     *       IERRC.EQ.21.OR.IERRC.EQ.210)) THEN
            DO 86 I=1,M
               QHAT(I,J)=QHAT(I,J)+H*(BHAT(1)*RK1(I,J)+
     1           BHAT(3)*RK3(I,J)+BHAT(4)*RK4(I,J)+
     2           BHAT(5)*RK5(I,J)+BHAT(6)*RK6(I,J)+BHAT(7)*RK2(I,J))
 86         CONTINUE
         ENDIF
 90   CONTINUE
      IF (IFXDPF.EQ.0.AND.(IERRC.EQ.2.OR.IERRC.EQ.20.OR.
     *    IERRC.EQ.21.OR.IERRC.EQ.210)) THEN
         CALL MODGRS(M,N,QHAT,M,BMAT,IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG=10
            RETURN
         ENDIF
      ENDIF
      RETURN
      END
CCCCC END OF PQRKDP CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: HQDPNS
C
C*****Purpose
C
C performs the one step integration of the system x'=f(x)
C and (depending on input options) of Q'=F(A,Q), Q0 given.
C The latter equation is integrated by integrating 
C the linear variational problem X'=DF(x(t))X, X0=Q0,
C and extracting Q.
C Monitors errors on the one step approximations to the 
C trajectory, and (if required) on the LEs and/or Q.
C It updates stepsize and approximation to the solutions
C trajectory and (if required) to the LEs and Q.
C 
C       Method used: Hybrid continuous QR method with DP5 as 
C          integrator for X,and  extraction of Q, As specified by
C          IPAR(9), two options are provided for 
C          approximating the LEs:
C          (i) a quadrature rule for the nu-variables of the 
C              same order as DP5 (the default), or
C          (ii) the composite trap-rule.  
C
C          Code can proceed in fixed or variable
C          step size.  In the latter case, error control is 
C          performed on the trajectory, and/or on the Lyapunov 
C          exponents, and/or on the Q-factor.
C
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC HQDPNS CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HQDPNS(GETF,GETDF,M,N,T,H,HOLD,Y0,X0,RLYAP,TOLT,
     *          TOLQ,TOLL,AMAT,X,XHAT,QMAT,BMAT,RK1,RK2,RK3,RK4,  
     *          RK5,RK6,STAGE,RNU,RNUHAT,RKNU1,RKNU2,RKNU3,RKNU4,
     *          RKNU5,RKNU6,ARK,BHAT,Y,YHAT,YRK1,YRK2,YRK3,YRK4,
     *          YRK5,YRK6,Y1,Y15,Y03,Y45,Y89,
     *          IPAR,NREJ,IFLAG,INARR,REARR)

      IMPLICIT NONE

      EXTERNAL GETF, GETDF
      INTEGER M,N,IPAR(*),NREJ,IFLAG,INARR(*)
      DOUBLE PRECISION T,H,HOLD, Y0(M),REARR(*)
      DOUBLE PRECISION X0(M,N), RLYAP(N), TOLT,TOLQ,TOLL(*)
      DOUBLE PRECISION AMAT(M,M), X(M,N), XHAT(M,N)
      DOUBLE PRECISION QMAT(M,N), BMAT(M,N)
      DOUBLE PRECISION RK1(M,N),RK2(M,N)
      DOUBLE PRECISION RK3(M,N),RK4(M,N)
      DOUBLE PRECISION RK5(M,N),RK6(M,N)
      DOUBLE PRECISION STAGE(M,N)     
      DOUBLE PRECISION RNU(N), RNUHAT(N)       
      DOUBLE PRECISION RKNU1(N), RKNU2(N), RKNU3(N), RKNU4(N)
      DOUBLE PRECISION RKNU5(N), RKNU6(N)
      DOUBLE PRECISION ARK(7,7), BHAT(7)
      DOUBLE PRECISION Y(M),YHAT(M),YRK1(M),YRK2(M),YRK3(M),YRK4(M)
      DOUBLE PRECISION YRK5(M),YRK6(M)
      DOUBLE PRECISION Y1(M),Y15(M),Y03(M),Y45(M),Y89(M)

      LOGICAL PREVIO, LAST
      INTEGER INTDIR, IFIRST, IFXDPF, IERRC, IHOWLE, IWHAT, I, J
      DOUBLE PRECISION HNEW, HMIN, HNEXT, TFIRST
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      DOUBLE PRECISION SCALE, RTOL, ERR
      COMMON /STEP/ HMIN, HNEXT, TFIRST, INTDIR, IFIRST, IWHAT, LAST
      SAVE /STEP/
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

      IFXDPF=IPAR(1)
      IERRC=IPAR(10)
      IHOWLE=IPAR(9)
      PREVIO=.FALSE.
C PREVIO IS .FALSE. IF PREVIOUS STEP WAS SUCCESSFUL
C IT IS .TRUE. IF PREVIOUS STEP WAS A FAILURE
 2    CONTINUE
CCCCCCCCCC INTEGRATOR FOR TRAJECTORY, AND Q-FACTOR
C OF TRANSITIONS.  RETURNS ALSO UPDATE FOR LES
      CALL HQRKDP(GETF,GETDF,PREVIO,M,N,IFIRST,IWHAT,IFXDPF,
     *        IHOWLE,IERRC,H,Y0,X0,
     *        AMAT,X,XHAT,QMAT,BMAT,RK1,RK2,RK3,RK4,RK5,RK6,
     *        STAGE, RNU,RNUHAT,RKNU1,RKNU2,RKNU3,RKNU4,
     *        RKNU5,RKNU6,ARK,BHAT,
     *        Y,YHAT,YRK1,YRK2,YRK3,YRK4,YRK5,YRK6,
     *        Y1,Y15,Y03,Y45,Y89,IFLAG,INARR,REARR)
      IF (IFIRST.EQ.1) IFIRST=0
      IF (IFLAG.NE.0) THEN
         IF (IFXDPF.EQ.0) THEN
C TRY DRASTIC WAY TO RECOVER
            IF (DABS(H).GT.HMIN) THEN
               IFLAG=0
               H=INTDIR*HMIN
               PREVIO=.TRUE.
               NREJ=NREJ+1
               GOTO 2
              ELSE
               IFLAG=5
               RETURN
            ENDIF
           ELSE
            IFLAG=10
            RETURN
         END IF
      END IF
      IF (IFXDPF.EQ.0) THEN
C HERE WE PREDICT NEXT STEPSIZE HNEW, NEVER ALLOWING IT TO
C BE GREATER THAN FIVE TIMES CURRENT H 
C (THIS VALUE CAN BE CHANGED IN HINCR IN "INIT")
         HNEW=HINCR*H
         RTOL=ZERO
         IF (IERRC.EQ.0) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
           ELSE IF (IERRC.EQ.1) THEN
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
           ELSE IF (IERRC.EQ.2) THEN
            CALL NERRQ(M,N,X,XHAT,X0,ERR,TOLQ,ONE,ZERO)
           ELSE IF (IERRC.EQ.10) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
           ELSE IF (IERRC.EQ.20) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRQ(M,N,X,XHAT,X0,ERR,TOLQ,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
           ELSE IF (IERRC.EQ.21) THEN
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRQ(M,N,X,XHAT,X0,ERR,TOLQ,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
           ELSE
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRQ(M,N,X,XHAT,X0,ERR,TOLQ,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
         END IF
C PREDICT NEW STEPSIZE; USE SAFETY FACTOR OF 0.8
C (THIS VALUE CAN BE CHANGED IN HSFTY IN "INIT")
         IF (ERR.GT.ZERO)         
     1      HNEW=DMIN1(DABS(HNEW),DABS(HSFTY*H*(ONE/ERR)**EXP1))
C GO BACK IF FAILED UNLESS H TOO SMALL
         IF (ERR.GT.ONE) THEN
C AVOID DRASTIC REDUCTION OF STEPSIZE (NEVER LESS THAN 1/5)
C (THIS VALUE CAN BE CHANGED IN HDRSTC IN "INIT")
            H=INTDIR*DMAX1(DABS(HNEW),HDRSTC*DABS(H))
            NREJ=NREJ+1
            IF (INTDIR*H.LE.HMIN) THEN
               PRINT*,' H = ', H
               IFLAG = 5
               RETURN
            ENDIF
            PREVIO=.TRUE.
            GOTO 2
         ENDIF            
      END IF
C HAVE COMPLETED THE STEP WITH STEPSIZE H
C SAVE LAST SUCCESSFUL STEPSIZE
      HOLD=H
      IF (IFXDPF.EQ.0) THEN
C         HNEW=INTDIR*DMAX1(DABS(HNEW),DABS(H))
         IF((PREVIO).AND.(DABS(HNEW).GT.DABS(H)))HNEW=H
         H=HNEW
      END IF
C FINALLY, UPDATE THE APPROXIMATION FOR THE LES, Q, AND TRAJECTORY
      IF (IWHAT.EQ.1) GOTO 29
      SCALE=ONE/(T+HOLD-TFIRST)      
      DO 25 I=1,N
         RLYAP(I)=((T-TFIRST)*SCALE)*RLYAP(I)+SCALE*RNU(I)
C LOWER ORDER APPROXIMATION
C         RLYAP(I)=((T-TFIRST)*SCALE)*RLYAP(I)+SCALE*RNUHAT(I)
 25   CONTINUE         
      DO 28 J=1,N
         DO 28 I=1,M
            X0(I,J)=X(I,J)
 28   CONTINUE
 29   CONTINUE
      DO 32 I=1,M
         Y0(I)=Y(I)
 32   CONTINUE

      RETURN
      END
CCCCC END OF HQDPNS CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: HQRKDP
C
C*****Purpose
C
C Integrates on one step the solution of x'=f(x).  If needed,
C approximates the Q factor in the QR factorization of X, 
C here X solves X'=A(t)X, X0=Q0, by using the hynrid scheme.
C It returns two approximations, x and xhat, to the solution, 
C two approximations X and XHAT to the Q-factor, and two 
C approximations rnu and rnuhat to the one-step contribution
C to the LEs.  All of these can be used for error purposes.
C
C Integrator is the DP rule with embedded FSAL of order 5 (4)
C 
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC HQRKDP CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HQRKDP(GETF,GETDF,PREVIO,M,N,IFIRST,IWHAT,IFXDPF,
     *        IHOWLE,IERRC,H,Y0,X0,
     *        AMAT,X,XHAT,QMAT,BMAT,RK1,RK2,RK3,RK4,RK5,RK6,
     *        STAGE, RNU,RNUHAT,RKNU1,RKNU2,RKNU3,RKNU4,
     *        RKNU5,RKNU6,ARK,BHAT,
     *        Y,YHAT,YRK1,YRK2,YRK3,YRK4,YRK5,YRK6,
     *        Y1,Y15,Y03,Y45,Y89,IFLAG,INARR,REARR)

      IMPLICIT NONE

      EXTERNAL GETF,GETDF
      LOGICAL PREVIO
      INTEGER M, N, IFIRST, IWHAT, IFXDPF, IHOWLE
      INTEGER IERRC, IFLAG,INARR(*)
      DOUBLE PRECISION H, Y0(M), X0(M,N),REARR(*)
      DOUBLE PRECISION X(M,N),AMAT(M,M), XHAT(M,N)
      DOUBLE PRECISION QMAT(M,N), BMAT(M,N)
      DOUBLE PRECISION RK1(M,N),RK2(M,N)
      DOUBLE PRECISION RK3(M,N),RK4(M,N)
      DOUBLE PRECISION RK5(M,N),RK6(M,N)
      DOUBLE PRECISION STAGE(M,N)
      DOUBLE PRECISION RNU(N),RNUHAT(N)
      DOUBLE PRECISION RKNU1(N),RKNU2(N),RKNU3(N)
      DOUBLE PRECISION RKNU4(N),RKNU5(N),RKNU6(N)
      DOUBLE PRECISION ARK(7,7), BHAT(7)
      DOUBLE PRECISION Y(M),YHAT(M),YRK1(M),YRK2(M),YRK3(M),YRK4(M)
      DOUBLE PRECISION YRK5(M),YRK6(M)
      DOUBLE PRECISION Y1(M),Y15(M),Y03(M),Y45(M),Y89(M)
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      INTEGER I, J
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

C INITIALIZE Y TO OLD Y0
      DO 2 I=1,M
         Y(I)=Y0(I)
 2       YHAT(I)=Y0(I)
C FORM RK1, RK2, RK3, RK4, RK5, RK6, RK7 
C FORM ALSO THE STAGE VALUES AT H/5, 3H/10, 4H/5, 8H/9 AND H (NB THE LATTER IS
C NOT THE SOLUTION VALUE, BUT JUST THE STAGE VALUE)
C ... RK1 AND U15 ...
C VERY FIRST CALL
C      IF (IFIRST.EQ.1) THEN
      IF ((IFIRST.EQ.1).OR.(IFXDPF.EQ.1)) THEN
         CALL GETF(M,Y,YRK1,INARR,REARR) 
         GOTO 7
      END IF
C IF NOT FAILED PREVIOUS STEP, YRK IS IN YRK2
C OTHERWISE, IT IS IN YRK1
      IF (.NOT.PREVIO) THEN
          DO 5 I=1,M 
 5           YRK1(I)=YRK2(I)
      ENDIF
 7    CONTINUE
      DO 11 I=1,M
         Y15(I) = Y(I) + H*ARK(2,1)*YRK1(I)
 11   CONTINUE
C ... RK2 AND U03 ...
      CALL GETF(M,Y15,YRK2,INARR,REARR) 
      DO 14 I=1,M
         Y03(I) = Y(I) + H*(ARK(3,1)*YRK1(I)+ARK(3,2)*YRK2(I))
 14   CONTINUE
C ... RK3 AND U45 
      CALL GETF(M,Y03,YRK3,INARR,REARR) 
      DO 16 I=1,M
         Y45(I) = Y(I) + H*(ARK(4,1)*YRK1(I)+ARK(4,2)*YRK2(I)+
     1                ARK(4,3)*YRK3(I))
 16   CONTINUE
C ... RK4 AND U89 
      CALL GETF(M,Y45,YRK4,INARR,REARR) 
      DO 19 I=1,M
         Y89(I) = Y(I) + H*(ARK(5,1)*YRK1(I)+ARK(5,2)*YRK2(I)+
     1                ARK(5,3)*YRK3(I)+ARK(5,4)*YRK4(I))
 19   CONTINUE
C ... RK5 AND U1 (STAGE VALUE AT H) ...
      CALL GETF(M,Y89,YRK5,INARR,REARR) 
      DO 20 I=1,M
         Y1(I) = Y(I) + H*(ARK(6,1)*YRK1(I)+ARK(6,2)*YRK2(I)+
     1                ARK(6,3)*YRK3(I)+ARK(6,4)*YRK4(I)+
     2                ARK(6,5)*YRK5(I))
 20   CONTINUE
C ... RK6, SOLUTION AT H ...
      CALL GETF(M,Y1,YRK6,INARR,REARR) 
      DO 23 I=1,M
         Y(I) = Y(I) + H*(ARK(7,1)*YRK1(I)+
     1         ARK(7,3)*YRK3(I)+ARK(7,4)*YRK4(I)+
     2         ARK(7,5)*YRK5(I)+ARK(7,6)*YRK6(I))
 23   CONTINUE
      IF (IFXDPF.EQ.1) GOTO 25
C ... RK7 (STORED IN RK2), AND COMPARISON SOLUTION XHAT ...
C TO BE USED FOR ERROR MONITOR/STEP-SIZE SELECTION
      CALL GETF(M,Y,YRK2,INARR,REARR) 
      DO 24 I=1,M 
         YHAT(I)=YHAT(I)+H*(BHAT(1)*YRK1(I)+BHAT(3)*YRK3(I)+
     1       BHAT(4)*YRK4(I)+BHAT(5)*YRK5(I)+BHAT(6)*YRK6(I)+
     2       BHAT(7)*YRK2(I))
 24   CONTINUE
C NOW DO IT FOR THE TRANSITION MATRIX (IF IWHAT=0)
 25   IF (IWHAT.EQ.1) RETURN
C INTIALIZE X TO OLD Q
C THE RNU'S ARE INITIALLY AT 0, SINCE WE TAKE ONE STEP OF INTEGRAL
      DO 26 J=1,N
         RNU(J)=ZERO
         RNUHAT(J)=ZERO
         DO 26 I=1,M
            X(I,J)=X0(I,J)
            XHAT(I,J)=X0(I,J)
 26    CONTINUE
C FORM RK1, RK2, RK3, RK4, RK5, RK6, RK7 
C FORM ALSO THE STAGE VALUES AT H/5, 3H/10, 4H/5, 8H/9 AND H (NB THE LATTER IS
C NOT THE SOLUTION VALUE, BUT JUST THE STAGE VALUE)
C ... RK1 AND U15 ...
C VERY FIRST CALL NEED AMAT AND RK1
      IF (IFIRST.EQ.1) THEN
         CALL XDOT(GETDF,M,N,Y0,AMAT,X,RK1,INARR,REARR) 
         GOTO 28
      END IF
C OTHERWISE IF NOT FAILED PREVIOUS STEP, HAVE AMAT IN AMAT
C BUT NEED RK1.  IF FAILED PREVIOUS STEP, HAVE RK1 IN RK1
      IF (.NOT.PREVIO) THEN
          CALL FORMAX(M,AMAT,M,M,X,N,M,RK1)
      ENDIF
 28   CONTINUE
      CALL DGQTB(M,M,N,X,RK1,M,BMAT,M)
      DO 29 I=1,N
 29      RKNU1(I)=BMAT(I,I)
      DO 30 J=1,N
         DO 30 I=1,M
            STAGE(I,J) = X(I,J) + H*ARK(2,1)*RK1(I,J)
 30   CONTINUE
C ... RK2 AND U03 ...
      CALL XDOT(GETDF,M,N,Y15,AMAT,STAGE,RK2,INARR,REARR) 
      IF (IHOWLE.EQ.1) GOTO 36
      DO 32 J=1,N
         DO 32 I=1,M
 32         QMAT(I,J)=STAGE(I,J)
      CALL MODGRS(M,N,QMAT,M,BMAT,IFLAG) 
      IF (IFLAG.NE.0) THEN
         IFLAG=10
         RETURN
      ENDIF
      CALL DGQTAQ(M,M,N,QMAT,AMAT,M,BMAT,M)
      DO 34 I=1,N
 34      RKNU2(I)=BMAT(I,I)
 36   CONTINUE
      DO 40 J=1,N
         DO 40 I=1,M
            STAGE(I,J) = X(I,J) + H*(ARK(3,1)*RK1(I,J)+
     1                 ARK(3,2)*RK2(I,J))
 40   CONTINUE
C ... RK3 AND U45 ...
      CALL XDOT(GETDF,M,N,Y03,AMAT,STAGE,RK3,INARR,REARR) 
      IF (IHOWLE.EQ.1) GOTO 47
      DO 41 J=1,N
         DO 41 I=1,M
 41         QMAT(I,J)=STAGE(I,J)
      CALL MODGRS(M,N,QMAT,M,BMAT,IFLAG) 
      IF (IFLAG.NE.0) THEN
         IFLAG=10
         RETURN
      ENDIF
      CALL DGQTAQ(M,M,N,QMAT,AMAT,M,BMAT,M)
      DO 43 I=1,N
 43      RKNU3(I)=BMAT(I,I)
 47   CONTINUE
      DO 50 J=1,N
         DO 50 I=1,M
            STAGE(I,J) = X(I,J) + H*(ARK(4,1)*RK1(I,J)+
     1                ARK(4,2)*RK2(I,J)+ARK(4,3)*RK3(I,J))
 50   CONTINUE
C ... RK4 AND U89 
      CALL XDOT(GETDF,M,N,Y45,AMAT,STAGE,RK4,INARR,REARR) 
      IF (IHOWLE.EQ.1) GOTO 56
      DO 52 J=1,N
         DO 52 I=1,M
 52         QMAT(I,J)=STAGE(I,J)
      CALL MODGRS(M,N,QMAT,M,BMAT,IFLAG) 
      IF (IFLAG.NE.0) THEN
         IFLAG=10
         RETURN
      ENDIF
      CALL DGQTAQ(M,M,N,QMAT,AMAT,M,BMAT,M)
      DO 55 I=1,N
 55      RKNU4(I)=BMAT(I,I)
 56   CONTINUE
       DO 60 J=1,N
         DO 60 I=1,M
            STAGE(I,J) = X(I,J) + H*(ARK(5,1)*RK1(I,J)+
     1                 ARK(5,2)*RK2(I,J)+
     2                 ARK(5,3)*RK3(I,J)+ARK(5,4)*RK4(I,J))
 60   CONTINUE
C ... RK5 AND U1 (STAGE VALUE AT H) ...
      CALL XDOT(GETDF,M,N,Y89,AMAT,STAGE,RK5,INARR,REARR) 
      IF (IHOWLE.EQ.1) GOTO 68
      DO 61 J=1,N
         DO 61 I=1,M
 61         QMAT(I,J)=STAGE(I,J)
      CALL MODGRS(M,N,QMAT,M,BMAT,IFLAG) 
      IF (IFLAG.NE.0) THEN
         IFLAG=10
         RETURN
      ENDIF
      CALL DGQTAQ(M,M,N,QMAT,AMAT,M,BMAT,M)
      DO 65 I=1,N
 65      RKNU5(I)=BMAT(I,I)
 68   CONTINUE
      DO 70 J=1,N
         DO 70 I=1,M
            STAGE(I,J) = X(I,J) + H*(ARK(6,1)*RK1(I,J)+
     1                ARK(6,2)*RK2(I,J)+ARK(6,3)*RK3(I,J)+
     2                ARK(6,4)*RK4(I,J)+ARK(6,5)*RK5(I,J))
 70   CONTINUE
C ... RK6, SOLUTION AT H ...
      CALL XDOT(GETDF,M,N,Y1,AMAT,STAGE,RK6,INARR,REARR) 
      IF (IHOWLE.EQ.1) GOTO 78
      DO 74 J=1,N
         DO 74 I=1,M
 74         QMAT(I,J)=STAGE(I,J)
      CALL MODGRS(M,N,QMAT,M,BMAT,IFLAG) 
      IF (IFLAG.NE.0) THEN
         IFLAG=10
         RETURN
      ENDIF
      CALL DGQTAQ(M,M,N,QMAT,AMAT,M,BMAT,M)
      DO 77 I=1,N
 77      RKNU6(I)=BMAT(I,I)
 78   CONTINUE
      DO 80 J=1,N
         IF (IHOWLE.EQ.0) THEN 
            RNU(J) = RNU(J) + H*(ARK(7,1)*RKNU1(J)+
     1         ARK(7,3)*RKNU3(J)+ARK(7,4)*RKNU4(J)+
     2         ARK(7,5)*RKNU5(J)+ARK(7,6)*RKNU6(J))
         ENDIF
         DO 80 I=1,M
            X(I,J) = X(I,J) + H*(ARK(7,1)*RK1(I,J)+
     1         ARK(7,3)*RK3(I,J)+ARK(7,4)*RK4(I,J)+
     2         ARK(7,5)*RK5(I,J)+ARK(7,6)*RK6(I,J))
 80   CONTINUE
C ... RK7 (STORED IN RK2), AND COMPARISON SOLUTION XHAT AND RNUHAT ...
C TO BE USED FOR ERROR MONITOR/STEP-SIZE SELECTION
      CALL XDOT(GETDF,M,N,Y,AMAT,X,RK2,INARR,REARR) 
      CALL MODGRS(M,N,X,M,BMAT,IFLAG) 
      IF (IFLAG.NE.0) THEN
         IFLAG=10
         RETURN
      ENDIF
      CALL DGQTAQ(M,M,N,X,AMAT,M,BMAT,M)
      DO 83 I=1,N
 83      RKNU2(I)=BMAT(I,I)
      DO 92 J=1,N 
         IF (IHOWLE.EQ.0.AND.IFXDPF.EQ.0.AND.(IERRC.EQ.1.OR.
     *       IERRC.EQ.10.OR.IERRC.EQ.21.OR.IERRC.EQ.210)) THEN
           RNUHAT(J)=RNUHAT(J)+H*(BHAT(1)*RKNU1(J)+
     1        BHAT(3)*RKNU3(J)+BHAT(4)*RKNU4(J)+
     2        BHAT(5)*RKNU5(J)+BHAT(6)*RKNU6(J)+BHAT(7)*RKNU2(J))
          ELSEIF (IHOWLE.EQ.1) THEN
            CALL TRAPRL(RKNU1(J),RKNU2(J),H,RNU(J))
         ENDIF
         IF (IFXDPF.EQ.0.AND.(IERRC.EQ.2.OR.IERRC.EQ.20.OR.
     *       IERRC.EQ.21.OR.IERRC.EQ.210)) THEN
            DO 90 I=1,M
               XHAT(I,J)=XHAT(I,J)+H*(BHAT(1)*RK1(I,J)+
     1          BHAT(3)*RK3(I,J)+BHAT(4)*RK4(I,J)+
     2          BHAT(5)*RK5(I,J)+BHAT(6)*RK6(I,J)+BHAT(7)*RK2(I,J))
 90         CONTINUE
         ENDIF
 92   CONTINUE
      IF (IFXDPF.EQ.0.AND.(IERRC.EQ.2.OR.IERRC.EQ.20.OR.
     *    IERRC.EQ.21.OR.IERRC.EQ.210)) THEN
         CALL MODGRS(M,N,XHAT,M,BMAT,IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG=10
            RETURN
         ENDIF
      ENDIF

      RETURN
      END
CCCCC END OF HQRKDP CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: PQ38NS
C
C*****Purpose
C
C performs the one step integration of y'=f(y)
C and (if required, see IPAR(6)), of Q'=F(A,Q), Q0 given.
C Monitors errors on the one step approximations to the 
C trajectory, and (if required) on the LEs and/or Q.
C It updates stepsize and approximation to the solutions
C trajectory and (if required) to the LEs and Q.
C
C       Method used: Complete Projection continuous QR method with 
C          RK38 as integrator for Q(i.e., also
C          for stage values). As specified by
C          IPAR(9), two options are provided for 
C          approximating the LEs:
C          (i) a quadrature rule for the nu-variables of the 
C              same order as RK38 (the default), or
C          (ii) the composite trap-rule.  
C
C          Code can proceed in fixed or variable
C          step size.  In the latter case, error control is 
C          performed on the trajectory, and/or on the Lyapunov 
C          exponents, and/or on the Q-factor.
C
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC PQ38NS CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PQ38NS(GETF,GETDF,M,N,T,H,HOLD,Y0,X0,RLYAP,TOLT,
     *          TOLQ,TOLL,AMAT,Q,QHAT,BMAT,RK1,RK2,RK3,RK4,STAGE,
     *          RNU,RNUHAT,RKNU1,RKNU2,RKNU3,RKNU4,
     *          ARK,BHAT,
     *          Y,YHAT,YRK1,YRK2,YRK3,YRK4,Y1,Y13,Y23,
     *          IPAR,NREJ,IFLAG,INARR,REARR)

      IMPLICIT NONE

      EXTERNAL GETF, GETDF
      INTEGER M,N,IPAR(*),NREJ,IFLAG, INARR(*)
      DOUBLE PRECISION T,H,HOLD, Y0(M), REARR(*)
      DOUBLE PRECISION X0(M,N), RLYAP(N), TOLT, TOLQ, TOLL(*)
      DOUBLE PRECISION AMAT(M,M), Q(M,N), QHAT(M,N)
      DOUBLE PRECISION BMAT(M,N)
      DOUBLE PRECISION RK1(M,N),RK2(M,N)
      DOUBLE PRECISION RK3(M,N),RK4(M,N)
      DOUBLE PRECISION STAGE(M,N)
      DOUBLE PRECISION RNU(N), RNUHAT(N)       
      DOUBLE PRECISION RKNU1(N), RKNU2(N), RKNU3(N), RKNU4(N)       
      DOUBLE PRECISION ARK(7,7), BHAT(7)
      DOUBLE PRECISION Y(M),YHAT(M),YRK1(M),YRK2(M),YRK3(M),YRK4(M)
      DOUBLE PRECISION Y1(M),Y13(M),Y23(M)

      LOGICAL PREVIO, LAST
      INTEGER INTDIR, IFIRST, IFXDPF, IERRC, IHOWLE, IWHAT, I, J
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      DOUBLE PRECISION HNEW, HMIN, HNEXT, TFIRST
      DOUBLE PRECISION SCALE, RTOL, ERR
      COMMON /STEP/ HMIN, HNEXT, TFIRST, INTDIR, IFIRST, IWHAT, LAST
      SAVE /STEP/
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

      IFXDPF=IPAR(1)
      IERRC=IPAR(10)
      IHOWLE=IPAR(9)
      PREVIO=.FALSE.
C PREVIO IS .FALSE. IF PREVIOUS STEP WAS SUCCESSFUL
C IT IS .TRUE. IF PREVIOUS STEP WAS A FAILURE
 2    CONTINUE
CCCCCCCCCC INTEGRATOR FOR TRAJECTORY, AND Q-FACTOR
C OF TRANSITIONS.  RETURNS ALSO UPDATE FOR LES
      CALL PQRK38(GETF,GETDF,PREVIO,M,N,IFIRST,IWHAT,IFXDPF,
     *        IHOWLE,IERRC,H,Y0,X0,
     *        AMAT,Q,QHAT,BMAT,RK1,RK2,RK3,RK4,STAGE,
     *        RNU,RNUHAT,RKNU1,RKNU2,RKNU3,RKNU4,ARK,BHAT,
     *        Y,YHAT,YRK1,YRK2,YRK3,YRK4,Y1,Y13,Y23,IFLAG,
     *        INARR,REARR)
      IF (IFIRST.EQ.1) IFIRST=0
      IF (IFLAG.NE.0) THEN
         IF (IFXDPF.EQ.0) THEN
C TRY DRASTIC WAY TO RECOVER
            IF (DABS(H).GT.HMIN) THEN
               IFLAG=0
               H=INTDIR*HMIN
               PREVIO=.TRUE.
               NREJ=NREJ+1
               GOTO 2
              ELSE
               IFLAG=5
               RETURN
            ENDIF
           ELSE
            IFLAG=10
            RETURN
         END IF
      END IF
      IF (IFXDPF.EQ.0) THEN
C HERE WE PREDICT NEXT STEPSIZE HNEW, NEVER ALLOWING IT TO
C BE GREATER THAN FIVE TIMES CURRENT H 
C (THIS VALUE CAN BE CHANGED IN HINCR IN "INIT")
         HNEW=HINCR*H
         RTOL=ZERO
         IF (IERRC.EQ.0) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
           ELSE IF (IERRC.EQ.1) THEN
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
           ELSE IF (IERRC.EQ.2) THEN
            CALL NERRQ(M,N,Q,QHAT,X0,ERR,TOLQ,ONE,ZERO)
           ELSE IF (IERRC.EQ.10) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
           ELSE IF (IERRC.EQ.20) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRQ(M,N,Q,QHAT,X0,ERR,TOLQ,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
           ELSE IF (IERRC.EQ.21) THEN
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRQ(M,N,Q,QHAT,X0,ERR,TOLQ,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
           ELSE
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRQ(M,N,Q,QHAT,X0,ERR,TOLQ,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
         END IF
C PREDICT NEW STEPSIZE; USE SAFETY FACTOR OF 0.8
C (THIS VALUE CAN BE CHANGED IN HSFTY IN "INIT")
         IF (ERR.GT.ZERO)         
     1         HNEW=DMIN1(DABS(HNEW),DABS(HSFTY*H*(ONE/ERR)**EXP2))
C GO BACK IF FAILED UNLESS H TOO SMALL
         IF (ERR.GT.ONE) THEN
C AVOID DRASTIC REDUCTION OF STEPSIZE (NEVER LESS THAN 1/5)
C (THIS VALUE CAN BE CHANGED IN HDRSTC IN "INIT")
            H=INTDIR*DMAX1(DABS(HNEW),HDRSTC*DABS(H))
            NREJ=NREJ+1
            IF (INTDIR*H.LE.HMIN) THEN
               PRINT*,' H = ', H
               IFLAG = 5
               RETURN
            ENDIF
            PREVIO=.TRUE.
            GOTO 2
         ENDIF            
      END IF
C HAVE COMPLETED THE STEP WITH STEPSIZE H
C SAVE LAST SUCCESSFUL STEPSIZE
      HOLD=H
      IF (IFXDPF.EQ.0) THEN
C         HNEW=INTDIR*DMAX1(DABS(HNEW),DABS(H))
         IF((PREVIO).AND.(DABS(HNEW).GT.DABS(H)))HNEW=H
         H=HNEW
      END IF
C FINALLY, UPDATE THE APPROXIMATION FOR THE LES, Q, AND TRAJECTORY
      IF (IWHAT.EQ.1) GOTO 29
      SCALE=ONE/(T+HOLD-TFIRST)      
      DO 25 I=1,N
         RLYAP(I)=((T-TFIRST)*SCALE)*RLYAP(I)+SCALE*RNU(I)
C LOWER ORDER APPROXIMATION
C         RLYAP(I)=((T-TFIRST)*SCALE)*RLYAP(I)+SCALE*RNUHAT(I)
 25   CONTINUE         
      DO 28 J=1,N
         DO 28 I=1,M
            X0(I,J)=Q(I,J)
 28   CONTINUE
 29   CONTINUE
      DO 32 I=1,M
         Y0(I)=Y(I)
 32   CONTINUE

      RETURN
      END
CCCCC END OF PQ38NS CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: PQRK38
C
C*****Purpose
C
C integrates on one step the solution of x'=f(x).  If needed,
C also of Q'=F(Df,Q).  It returns two approximations, x and xhat,
C and Q, the Q-factor in the QR factorization of X.
C It also returns two approximations to the updates of
C the LEs, rnu and rnuhat. All of these can be used for error purposes.
C
C Integrator is the 3/8 rule with embedded FSAL of order 4 (3)
C Scheme is total projection method
C 
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC PQRK38 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PQRK38(GETF,GETDF,PREVIO,M,N,IFIRST,IWHAT,IFXDPF,
     *        IHOWLE,IERRC,H,Y0,X0,
     *        AMAT,Q,QHAT,BMAT,RK1,RK2,RK3,RK4,STAGE,
     *        RNU,RNUHAT,RKNU1,RKNU2,RKNU3,RKNU4,ARK,BHAT,
     *        Y,YHAT,YRK1,YRK2,YRK3,YRK4,Y1,Y13,Y23,IFLAG,
     *        INARR,REARR)

      IMPLICIT NONE

      EXTERNAL GETF,GETDF
      LOGICAL PREVIO
      INTEGER M, N, IFIRST, IWHAT, IFXDPF, IHOWLE
      INTEGER IERRC, IFLAG,INARR(*)
      DOUBLE PRECISION H, Y0(M), X0(M,N), REARR(*)
      DOUBLE PRECISION AMAT(M,M), Q(M,N), QHAT(M,N)
      DOUBLE PRECISION BMAT(M,N)
      DOUBLE PRECISION RK1(M,N),RK2(M,N),RK3(M,N),RK4(M,N)
      DOUBLE PRECISION STAGE(M,N)
      DOUBLE PRECISION RNU(N),RNUHAT(N)
      DOUBLE PRECISION RKNU1(N),RKNU2(N),RKNU3(N),RKNU4(N)
      DOUBLE PRECISION ARK(7,7), BHAT(7)
      DOUBLE PRECISION Y(M),YHAT(M),YRK1(M),YRK2(M),YRK3(M),YRK4(M)
      DOUBLE PRECISION Y1(M),Y13(M),Y23(M)
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      INTEGER I, J
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

C INITIALIZE Y TO OLD Y0
      DO 2 I=1,M
         Y(I)=Y0(I)
 2       YHAT(I)=Y0(I)
C FORM RK1, RK2, RK3, RK4 AND RK5 
C FORM ALSO THE STAGE VALUES AT H/3, 2H/3 AND H (NB THE LATTER IS
C NOT THE SOLUTION VALUE, BUT JUST THE STAGE VALUE)
C ... RK1 AND U13 ...
C VERY FIRST CALL
C      IF (IFIRST.EQ.1) THEN
      IF ((IFIRST.EQ.1).OR.(IFXDPF.EQ.1)) THEN
         CALL GETF(M,Y,YRK1,INARR,REARR) 
         GOTO 7
      END IF
C IF NOT FAILED PREVIOUS STEP, YRK1 IS IN YRK4
C OTHERWISE, IT IS IN YRK1
      IF (.NOT.PREVIO) THEN
          DO 5 I=1,M 
 5           YRK1(I)=YRK4(I)
      ENDIF
 7    CONTINUE
      DO 11 I=1,M
         Y13(I) = Y(I) + H*ARK(2,1)*YRK1(I)
 11   CONTINUE
C ... RK2 AND U23 ...
      CALL GETF(M,Y13,YRK2,INARR,REARR) 
      DO 14 I=1,M
         Y23(I) = Y(I) + H*(ARK(3,1)*YRK1(I)+ARK(3,2)*YRK2(I))
 14   CONTINUE
C ... RK3 AND U1 (STAGE VALUE AT H) ...
      CALL GETF(M,Y23,YRK3,INARR,REARR) 
      DO 16 I=1,M
         Y1(I) = Y(I) + H*(ARK(4,1)*YRK1(I)+ARK(4,2)*YRK2(I)+
     1                ARK(4,3)*YRK3(I))
 16   CONTINUE
C ... RK4, SOLUTION AT H ...
      CALL GETF(M,Y1,YRK4,INARR,REARR) 
      DO 19 I=1,M
         Y(I) = Y(I) + H*(ARK(5,1)*YRK1(I)+
     1         ARK(5,2)*YRK2(I)+ARK(5,3)*YRK3(I)+ARK(5,4)*YRK4(I))
 19   CONTINUE
      IF (IFXDPF.EQ.1) GOTO 22
C ... RK5 (STORED IN RK4), AND COMPARISON SOLUTION YHAT ...
C TO BE USED FOR ERROR MONITOR/STEP-SIZE SELECTION
      CALL GETF(M,Y,YRK4,INARR,REARR) 
      DO 20 I=1,M 
         YHAT(I)=YHAT(I)+H*(BHAT(1)*YRK1(I)+
     1       BHAT(2)*YRK2(I)+BHAT(3)*YRK3(I)+BHAT(5)*YRK4(I))
 20   CONTINUE
C NOW DO IT FOR THE TRANSITION MATRIX (IF IWHAT=0)
 22   IF (IWHAT.EQ.1) RETURN
C INTIALIZE Q TO OLD Q
C THE RNU'S ARE INITIALLY AT 0, SINCE WE TAKE ONE STEP OF INTEGRAL
      DO 23 J=1,N
         RNU(J)=ZERO
         RNUHAT(J)=ZERO
         DO 23 I=1,M
            Q(I,J)=X0(I,J)
            QHAT(I,J)=X0(I,J)
 23   CONTINUE
C FORM RK1, RK2, RK3, RK4 AND RK5 
C FORM ALSO THE STAGE VALUES AT H/3, 2H/3 AND H (NB THE LATTER IS
C NOT THE SOLUTION VALUE, BUT JUST THE STAGE VALUE)
C ... RK1 AND U13 ...
C VERY FIRST CALL NEED RK1
C IF FAILED PREVIOUS STEP, HAVE RK1 IN RK1
C IF NOT FAILED PREVIOUS STEP, NEED RK1.  
CC **** SEE IF CAN ECONOMIZE HERE ****
      CALL QDOT(GETDF,M,N,Y0,AMAT,BMAT,Q,RK1,INARR,REARR)
      DO 25 J=1,N
         DO 25 I=1,M
            STAGE(I,J) = Q(I,J) + H*ARK(2,1)*RK1(I,J)
 25   CONTINUE
      DO 29 I=1,N
 29      RKNU1(I)=BMAT(I,I)
      IF (IHOWLE.EQ.0) THEN
         CALL MODGRS(M,N,STAGE,M,BMAT,IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG=10
            RETURN
         ENDIF
      ENDIF
C ... RK2 AND U23 ...
      CALL QDOT(GETDF,M,N,Y13,AMAT,BMAT,STAGE,RK2,INARR,REARR)
      DO 32 J=1,N
         DO 32 I=1,M
            STAGE(I,J) = Q(I,J) + H*(ARK(3,1)*RK1(I,J)+
     1                 ARK(3,2)*RK2(I,J))
 32   CONTINUE
      IF (IHOWLE.EQ.0) THEN
         DO 34 I=1,N
 34         RKNU2(I)=BMAT(I,I)
         CALL MODGRS(M,N,STAGE,M,BMAT,IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG=10
            RETURN
         ENDIF
      ENDIF
C ... RK3 AND U1 (STAGE VALUE AT H) ...
      CALL QDOT(GETDF,M,N,Y23,AMAT,BMAT,STAGE,RK3,INARR,REARR)
      DO 38 J=1,N
         DO 38 I=1,M
            STAGE(I,J) = Q(I,J) + H*(ARK(4,1)*RK1(I,J)+
     1                ARK(4,2)*RK2(I,J)+ARK(4,3)*RK3(I,J))
 38   CONTINUE
      IF (IHOWLE.EQ.0) THEN
         DO 41 I=1,N
 41         RKNU3(I)=BMAT(I,I)
         CALL MODGRS(M,N,STAGE,M,BMAT,IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG=10
            RETURN
         ENDIF
      ENDIF
C ... RK4, SOLUTION AT H ...
      CALL QDOT(GETDF,M,N,Y1,AMAT,BMAT,STAGE,RK4,INARR,REARR)
      IF (IHOWLE.EQ.0) THEN
         DO 43 I=1,N
 43         RKNU4(I)=BMAT(I,I)
      ENDIF
      DO 47 J=1,N
         IF (IHOWLE.EQ.0) THEN 
            RNU(J)=RNU(J)+H*(ARK(5,1)*RKNU1(J)+ARK(5,2)*RKNU2(J)+
     1                ARK(5,3)*RKNU3(J)+ARK(5,4)*RKNU4(J))
         ENDIF
         DO 47 I=1,M
            Q(I,J) = Q(I,J) + H*(ARK(5,1)*RK1(I,J)+
     1         ARK(5,2)*RK2(I,J)+ARK(5,3)*RK3(I,J)+ARK(5,4)*RK4(I,J))
 47   CONTINUE
C ... RK5 (STORED IN RK4), AND COMPARISON SOLUTION QHAT AND RNUHAT...
C TO BE USED FOR ERROR MONITOR/STEP-SIZE SELECTION
      CALL MODGRS(M,N,Q,M,BMAT,IFLAG) 
      IF (IFLAG.NE.0) THEN
         IFLAG=10
         RETURN
      ENDIF
      IF (((IERRC.EQ.2).OR.(IERRC.EQ.20).OR.(IERRC.EQ.21).OR.
     *    (IERRC.EQ.210)).AND.IFXDPF.EQ.0) THEN
         CALL QDOT(GETDF,M,N,Y,AMAT,BMAT,Q,RK4,INARR,REARR)
        ELSE
         CALL GETDF(M,Y,AMAT,INARR,REARR)
         CALL DGQTAQ(M,M,N,Q,AMAT,M,BMAT,M)
      ENDIF
      DO 50 I=1,N
 50      RKNU4(I)=BMAT(I,I)
      DO 55 J=1,N 
         IF (IHOWLE.EQ.0.AND.IFXDPF.EQ.0.AND.(IERRC.EQ.1.OR.
     *       IERRC.EQ.10.OR.IERRC.EQ.21.OR.IERRC.EQ.210)) THEN
            RNUHAT(J)=RNUHAT(J)+H*(BHAT(1)*RKNU1(J)+BHAT(2)*RKNU2(J)
     1                +BHAT(3)*RKNU3(J)+BHAT(5)*RKNU4(J))
           ELSEIF (IHOWLE.EQ.1) THEN
             CALL TRAPRL(RKNU1(J),RKNU4(J),H,RNU(J))
         ENDIF
         IF (IFXDPF.EQ.0.AND.(IERRC.EQ.2.OR.IERRC.EQ.20.OR.
     *       IERRC.EQ.21.OR.IERRC.EQ.210)) THEN
            DO 52 I=1,M
               QHAT(I,J)=QHAT(I,J)+H*(BHAT(1)*RK1(I,J)+
     1           BHAT(2)*RK2(I,J)+BHAT(3)*RK3(I,J)+BHAT(5)*RK4(I,J))
 52         CONTINUE
         ENDIF
 55   CONTINUE
      IF (IFXDPF.EQ.0.AND.(IERRC.EQ.2.OR.IERRC.EQ.20.OR.
     *    IERRC.EQ.21.OR.IERRC.EQ.210)) THEN
         CALL MODGRS(M,N,QHAT,M,BMAT,IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG=10
            RETURN
         ENDIF
      ENDIF
      RETURN
      END
CCCCC END OF PQRK38 CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: HQ38NS
C
C*****Purpose
C
C performs the one step integration of the system x'=f(x)
C and (depending on input options) of Q'=F(A,Q), Q0 given.
C The latter equation is integrated by integrating 
C the linear variational problem X'=DF(x(t))X, X0=Q0,
C and extracting Q.
C Monitors errors on the one step approximations to the 
C trajectory, and (if required) on the LEs and/or Q.
C It updates stepsize and approximation to the solutions
C trajectory and (if required) to the LEs and Q.
C 
C       Method used: Hybrid continuous QR method with DRK38 as 
C          integrator for X,and  extraction of Q, As specified by
C          IPAR(9) (see below), two options are provided for 
C          approximating the LEs:
C          (i) a quadrature rule for the nu-variables of the 
C              same order as DP5 (the default), or
C          (ii) the composite trap-rule.  
C
C          Code can proceed in fixed or variable
C          step size.  In the latter case, error control is 
C          performed on the trajectory, and/or on the Lyapunov 
C          exponents, and/or on the Q-factor.
C
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC HQ38NS CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HQ38NS(GETF,GETDF,M,N,T,H,HOLD,Y0,X0,RLYAP,TOLT,
     *          TOLQ,TOLL,AMAT,X,XHAT,QMAT,BMAT,RK1,RK2,RK3,RK4,STAGE,
     *          RNU,RNUHAT,RKNU1,RKNU2,RKNU3,RKNU4,
     *          ARK,BHAT,
     *          Y,YHAT,YRK1,YRK2,YRK3,YRK4,Y1,Y13,Y23,
     *          IPAR,NREJ,IFLAG,INARR,REARR)

      IMPLICIT NONE

      EXTERNAL GETF, GETDF
      INTEGER M,N,IPAR(*),NREJ,IFLAG, INARR(*)
      DOUBLE PRECISION T,H,HOLD, Y0(M), REARR(*)
      DOUBLE PRECISION X0(M,N), RLYAP(N), TOLT,TOLQ,TOLL(*)
      DOUBLE PRECISION AMAT(M,M), X(M,N), XHAT(M,N)
      DOUBLE PRECISION QMAT(M,N), BMAT(M,N)
      DOUBLE PRECISION RK1(M,N),RK2(M,N)
      DOUBLE PRECISION RK3(M,N),RK4(M,N)
      DOUBLE PRECISION STAGE(M,N)
      DOUBLE PRECISION RNU(N), RNUHAT(N)       
      DOUBLE PRECISION RKNU1(N), RKNU2(N), RKNU3(N), RKNU4(N)       
      DOUBLE PRECISION ARK(7,7), BHAT(7)
      DOUBLE PRECISION Y(M),YHAT(M),YRK1(M),YRK2(M),YRK3(M),YRK4(M)
      DOUBLE PRECISION Y1(M),Y13(M),Y23(M)

      LOGICAL PREVIO, LAST
      INTEGER INTDIR, IFIRST, IFXDPF, IERRC, IHOWLE, IWHAT, I, J
      DOUBLE PRECISION HNEW, HMIN, HNEXT, TFIRST
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      DOUBLE PRECISION SCALE, RTOL, ERR
      COMMON /STEP/ HMIN, HNEXT, TFIRST, INTDIR, IFIRST, IWHAT, LAST
      SAVE /STEP/
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

      IFXDPF=IPAR(1)
      IERRC=IPAR(10)
      IHOWLE=IPAR(9)
      PREVIO=.FALSE.
C PREVIO IS .FALSE. IF PREVIOUS STEP WAS SUCCESSFUL
C IT IS .TRUE. IF PREVIOUS STEP WAS A FAILURE
 2    CONTINUE
CCCCCCCCCC INTEGRATOR FOR TRAJECTORY, AND Q-FACTOR
C OF TRANSITIONS.  RETURNS ALSO UPDATE FOR LES
      CALL HQRK38(GETF,GETDF,PREVIO,M,N,IFIRST,IWHAT,IFXDPF,
     *        IHOWLE,IERRC,H,Y0,X0,
     *        AMAT,X,XHAT,QMAT,BMAT,RK1,RK2,RK3,RK4,STAGE,
     *        RNU,RNUHAT,RKNU1,RKNU2,RKNU3,RKNU4,ARK,BHAT,
     *        Y,YHAT,YRK1,YRK2,YRK3,YRK4,Y1,Y13,Y23,IFLAG,
     *        INARR,REARR)
      IF (IFIRST.EQ.1) IFIRST=0
      IF (IFLAG.NE.0) THEN
         IF (IFXDPF.EQ.0) THEN
C TRY DRASTIC WAY TO RECOVER
            IF (DABS(H).GT.HMIN) THEN
               IFLAG=0
               H=INTDIR*HMIN
               PREVIO=.TRUE.
               NREJ=NREJ+1
               GOTO 2
              ELSE
               IFLAG=5
               RETURN
            ENDIF
           ELSE
            IFLAG=10
            RETURN
         END IF
      END IF
      IF (IFXDPF.EQ.0) THEN
C HERE WE PREDICT NEXT STEPSIZE HNEW, NEVER ALLOWING IT TO
C BE GREATER THAN FIVE TIMES CURRENT H 
C (THIS VALUE CAN BE CHANGED IN HINCR IN "INIT")
         HNEW=HINCR*H
         RTOL=ZERO
         IF (IERRC.EQ.0) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
           ELSE IF (IERRC.EQ.1) THEN
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
           ELSE IF (IERRC.EQ.2) THEN
            CALL NERRQ(M,N,X,XHAT,X0,ERR,TOLQ,ONE,ZERO)
           ELSE IF (IERRC.EQ.10) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
           ELSE IF (IERRC.EQ.20) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRQ(M,N,X,XHAT,X0,ERR,TOLQ,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
           ELSE IF (IERRC.EQ.21) THEN
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRQ(M,N,X,XHAT,X0,ERR,TOLQ,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
           ELSE
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRNU(N,RNU,RNUHAT,ERR,TOLL,ONE,ZERO)
            RTOL=DMAX1(RTOL,ERR)
            CALL NERRQ(M,N,X,XHAT,X0,ERR,TOLQ,ONE,ZERO)
            ERR=DMAX1(ERR,RTOL)
         END IF
C PREDICT NEW STEPSIZE; USE SAFETY FACTOR OF 0.8
C (THIS VALUE CAN BE CHANGED IN HSFTY IN "INIT")
         IF (ERR.GT.ZERO)         
     1         HNEW=DMIN1(DABS(HNEW),DABS(HSFTY*H*(ONE/ERR)**EXP2))
C GO BACK IF FAILED UNLESS H TOO SMALL
         IF (ERR.GT.ONE) THEN
C AVOID DRASTIC REDUCTION OF STEPSIZE (NEVER LESS THAN 1/5)
C (THIS VALUE CAN BE CHANGED IN HDRSTC IN "INIT")
            H=INTDIR*DMAX1(DABS(HNEW),HDRSTC*DABS(H))
            NREJ=NREJ+1
            IF (INTDIR*H.LE.HMIN) THEN
               PRINT*,' H = ', H
               IFLAG = 5
               RETURN
            ENDIF
            PREVIO=.TRUE.
            GOTO 2
         ENDIF            
      END IF
C HAVE COMPLETED THE STEP WITH STEPSIZE H
C SAVE LAST SUCCESSFUL STEPSIZE
      HOLD=H
      IF (IFXDPF.EQ.0) THEN
C         HNEW=INTDIR*DMAX1(DABS(HNEW),DABS(H))
         IF((PREVIO).AND.(DABS(HNEW).GT.DABS(H)))HNEW=H
         H=HNEW
      END IF
C FINALLY, UPDATE THE APPROXIMATION FOR THE LES, Q, AND TRAJECTORY
      IF (IWHAT.EQ.1) GOTO 29
      SCALE=ONE/(T+HOLD-TFIRST)      
      DO 25 I=1,N
         RLYAP(I)=((T-TFIRST)*SCALE)*RLYAP(I)+SCALE*RNU(I)
C LOWER ORDER APPROXIMATION
C         RLYAP(I)=((T-TFIRST)*SCALE)*RLYAP(I)+SCALE*RNUHAT(I)
 25   CONTINUE         
      DO 28 J=1,N
         DO 28 I=1,M
            X0(I,J)=X(I,J)
 28   CONTINUE
 29   CONTINUE
      DO 32 I=1,M
         Y0(I)=Y(I)
 32   CONTINUE

      RETURN
      END
CCCCC END OF HQ38NS CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: HQRK38
C
C*****Purpose
C
C Integrates on one step the solution of x'=f(x).  If needed,
C approximates the Q factor in the QR factorization of X, 
C here X solves X'=A(t)X, X0=Q0, by using the hynrid scheme.
C It returns two approximations, x and xhat, to the solution, 
C two approximations X and XHAT to the Q-factor, and two 
C approximations rnu and rnuhat to the one-step contribution
C to the LEs.  All of thesecan be used for error purposes.
C
C Integrator is the 3/8 rule with embedded FSAL of order 4 (3)
C 
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC HQRK38 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HQRK38(GETF,GETDF,PREVIO,M,N,IFIRST,IWHAT,IFXDPF,
     *        IHOWLE,IERRC,H,Y0,X0,
     *        AMAT,X,XHAT,QMAT,BMAT,RK1,RK2,RK3,RK4,STAGE,
     *        RNU,RNUHAT,RKNU1,RKNU2,RKNU3,RKNU4,ARK,BHAT,
     *        Y,YHAT,YRK1,YRK2,YRK3,YRK4,Y1,Y13,Y23,IFLAG,
     *        INARR,REARR)

      IMPLICIT NONE

      EXTERNAL GETF,GETDF
      LOGICAL PREVIO
      INTEGER M, N, IFIRST, IWHAT, IFXDPF, IHOWLE
      INTEGER IERRC, IFLAG,INARR(*)
      DOUBLE PRECISION H, Y0(M), X0(M,N), REARR(*)
      DOUBLE PRECISION X(M,N),AMAT(M,M), XHAT(M,N)
      DOUBLE PRECISION QMAT(M,N), BMAT(M,N)
      DOUBLE PRECISION RK1(M,N),RK2(M,N),RK3(M,N),RK4(M,N)
      DOUBLE PRECISION STAGE(M,N)
      DOUBLE PRECISION RNU(N),RNUHAT(N)
      DOUBLE PRECISION RKNU1(N),RKNU2(N),RKNU3(N),RKNU4(N)
      DOUBLE PRECISION ARK(7,7), BHAT(7)
      DOUBLE PRECISION Y(M),YHAT(M),YRK1(M),YRK2(M),YRK3(M),YRK4(M)
      DOUBLE PRECISION Y1(M),Y13(M),Y23(M)
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      INTEGER I, J
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

C INITIALIZE Y TO OLD Y0
      DO 2 I=1,M
         Y(I)=Y0(I)
 2       YHAT(I)=Y0(I)
C FORM RK1, RK2, RK3, RK4 AND RK5 
C FORM ALSO THE STAGE VALUES AT H/3, 2H/3 AND H (NB THE LATTER IS
C NOT THE SOLUTION VALUE, BUT JUST THE STAGE VALUE)
C ... RK1 AND U13 ...
C VERY FIRST CALL
C      IF (IFIRST.EQ.1) THEN
      IF ((IFIRST.EQ.1).OR.(IFXDPF.EQ.1)) THEN
         CALL GETF(M,Y,YRK1,INARR,REARR) 
         GOTO 7
      END IF
C IF NOT FAILED PREVIOUS STEP, YRK1 IS IN YRK4
C OTHERWISE, IT IS IN YRK1
      IF (.NOT.PREVIO) THEN
          DO 5 I=1,M 
 5           YRK1(I)=YRK4(I)
      ENDIF
 7    CONTINUE
      DO 11 I=1,M
         Y13(I) = Y(I) + H*ARK(2,1)*YRK1(I)
 11   CONTINUE
C ... RK2 AND U23 ...
      CALL GETF(M,Y13,YRK2,INARR,REARR) 
      DO 14 I=1,M
         Y23(I) = Y(I) + H*(ARK(3,1)*YRK1(I)+ARK(3,2)*YRK2(I))
 14   CONTINUE
C ... RK3 AND U1 (STAGE VALUE AT H) ...
      CALL GETF(M,Y23,YRK3,INARR,REARR) 
      DO 16 I=1,M
         Y1(I) = Y(I) + H*(ARK(4,1)*YRK1(I)+ARK(4,2)*YRK2(I)+
     1                ARK(4,3)*YRK3(I))
 16   CONTINUE
C ... RK4, SOLUTION AT H ...
      CALL GETF(M,Y1,YRK4,INARR,REARR) 
      DO 19 I=1,M
         Y(I) = Y(I) + H*(ARK(5,1)*YRK1(I)+
     1         ARK(5,2)*YRK2(I)+ARK(5,3)*YRK3(I)+ARK(5,4)*YRK4(I))
 19   CONTINUE
      IF (IFXDPF.EQ.1) GOTO 22
C ... RK5 (STORED IN RK4), AND COMPARISON SOLUTION XHAT ...
C TO BE USED FOR ERROR MONITOR/STEP-SIZE SELECTION
      CALL GETF(M,Y,YRK4,INARR,REARR) 
      DO 20 I=1,M 
         YHAT(I)=YHAT(I)+H*(BHAT(1)*YRK1(I)+
     1       BHAT(2)*YRK2(I)+BHAT(3)*YRK3(I)+BHAT(5)*YRK4(I))
 20   CONTINUE
C NOW DO IT FOR THE TRANSITION MATRIX (IF IWHAT=0)
 22   IF (IWHAT.EQ.1) RETURN
C INTIALIZE X TO OLD Q
C THE RNU'S ARE INITIALLY AT 0, SINCE WE TAKE ONE STEP OF INTEGRAL
      DO 23 J=1,N
         RNU(J)=ZERO
         RNUHAT(J)=ZERO
         DO 23 I=1,M
            X(I,J)=X0(I,J)
            XHAT(I,J)=X0(I,J)
 23   CONTINUE
C FORM RK1, RK2, RK3, RK4 AND RK5 
C FORM ALSO THE STAGE VALUES AT H/3, 2H/3 AND H (NB THE LATTER IS
C NOT THE SOLUTION VALUE, BUT JUST THE STAGE VALUE)
C ... RK1 AND U13 ...
C VERY FIRST CALL NEED AMAT AND RK1
      IF (IFIRST.EQ.1) THEN
         CALL XDOT(GETDF,M,N,Y0,AMAT,X,RK1,INARR,REARR) 
         GOTO 24
      END IF
C OTHERWISE IF NOT FAILED PREVIOUS STEP, HAVE AMAT IN AMAT
C BUT NEED RK1.  IF FAILED PREVIOUS STEP, HAVE RK1 IN RK1
      IF (.NOT.PREVIO) THEN
         CALL FORMAX(M,AMAT,M,M,X,N,M,RK1)
      ENDIF
 24   CONTINUE
      CALL DGQTB(M,M,N,X,RK1,M,BMAT,M)
      DO 25 I=1,N
 25      RKNU1(I)=BMAT(I,I)
      DO 30 J=1,N
         DO 30 I=1,M
            STAGE(I,J) = X(I,J) + H*ARK(2,1)*RK1(I,J)
 30   CONTINUE
C ... RK2 AND U23 ...
      CALL XDOT(GETDF,M,N,Y13,AMAT,STAGE,RK2,INARR,REARR) 
      IF (IHOWLE.EQ.1) GOTO 36
      DO 32 J=1,N
         DO 32 I=1,M
 32         QMAT(I,J)=STAGE(I,J)
      CALL MODGRS(M,N,QMAT,M,BMAT,IFLAG) 
      IF (IFLAG.NE.0) THEN
         IFLAG=10
         RETURN
      ENDIF
      CALL DGQTAQ(M,M,N,QMAT,AMAT,M,BMAT,M)
      DO 34 I=1,N
 34      RKNU2(I)=BMAT(I,I)
 36   CONTINUE
      DO 40 J=1,N
         DO 40 I=1,M
            STAGE(I,J) = X(I,J) + H*(ARK(3,1)*RK1(I,J)+
     1                 ARK(3,2)*RK2(I,J))
 40   CONTINUE
C ... RK3 AND U1 (STAGE VALUE AT H) ...
      CALL XDOT(GETDF,M,N,Y23,AMAT,STAGE,RK3,INARR,REARR) 
      IF (IHOWLE.EQ.1) GOTO 47
      DO 41 J=1,N
         DO 41 I=1,M
 41         QMAT(I,J)=STAGE(I,J)
      CALL MODGRS(M,N,QMAT,M,BMAT,IFLAG) 
      IF (IFLAG.NE.0) THEN
         IFLAG=10
         RETURN
      ENDIF
      CALL DGQTAQ(M,M,N,QMAT,AMAT,M,BMAT,M)
      DO 43 I=1,N
 43      RKNU3(I)=BMAT(I,I)
 47   CONTINUE
      DO 50 J=1,N
         DO 50 I=1,M
            STAGE(I,J) = X(I,J) + H*(ARK(4,1)*RK1(I,J)+
     1                ARK(4,2)*RK2(I,J)+ARK(4,3)*RK3(I,J))
 50   CONTINUE
C ... RK4, SOLUTION AT H ...
      CALL XDOT(GETDF,M,N,Y1,AMAT,STAGE,RK4,INARR,REARR) 
      IF (IHOWLE.EQ.1) GOTO 56
      DO 52 J=1,N
         DO 52 I=1,M
 52         QMAT(I,J)=STAGE(I,J)
      CALL MODGRS(M,N,QMAT,M,BMAT,IFLAG) 
      IF (IFLAG.NE.0) THEN
         IFLAG=10
         RETURN
      ENDIF
      CALL DGQTAQ(M,M,N,QMAT,AMAT,M,BMAT,M)
      DO 55 I=1,N
 55      RKNU4(I)=BMAT(I,I)
 56   CONTINUE
      DO 60 J=1,N
         IF (IHOWLE.EQ.0) THEN 
            RNU(J)=RNU(J)+H*(ARK(5,1)*RKNU1(J)+ARK(5,2)*RKNU2(J)+
     1                ARK(5,3)*RKNU3(J)+ARK(5,4)*RKNU4(J))
         ENDIF
         DO 60 I=1,M
            X(I,J) = X(I,J) + H*(ARK(5,1)*RK1(I,J)+
     1         ARK(5,2)*RK2(I,J)+ARK(5,3)*RK3(I,J)+ARK(5,4)*RK4(I,J))
 60   CONTINUE
C ... RK5 (STORED IN RK4), AND COMPARISON SOLUTION XHAT AND RNUHAT...
C TO BE USED FOR ERROR MONITOR/STEP-SIZE SELECTION
      CALL XDOT(GETDF,M,N,Y,AMAT,X,RK4,INARR,REARR) 
      CALL MODGRS(M,N,X,M,BMAT,IFLAG) 
      IF (IFLAG.NE.0) THEN
         IFLAG=10
         RETURN
      ENDIF
      CALL DGQTAQ(M,M,N,X,AMAT,M,BMAT,M)
      DO 65 I=1,N
 65      RKNU4(I)=BMAT(I,I)
      DO 74 J=1,N 
         IF (IHOWLE.EQ.0.AND.IFXDPF.EQ.0.AND.(IERRC.EQ.1.OR.
     *       IERRC.EQ.10.OR.IERRC.EQ.21.OR.IERRC.EQ.210)) THEN
            RNUHAT(J)=RNUHAT(J)+H*(BHAT(1)*RKNU1(J)+BHAT(2)*RKNU2(J)
     1                +BHAT(3)*RKNU3(J)+BHAT(5)*RKNU4(J))
          ELSEIF (IHOWLE.EQ.1) THEN
            CALL TRAPRL(RKNU1(J),RKNU4(J),H,RNU(J))
         ENDIF
         IF (IFXDPF.EQ.0.AND.(IERRC.EQ.2.OR.IERRC.EQ.20.OR.
     *       IERRC.EQ.21.OR.IERRC.EQ.210)) THEN
            DO 70 I=1,M
               XHAT(I,J)=XHAT(I,J)+H*(BHAT(1)*RK1(I,J)+
     1          BHAT(2)*RK2(I,J)+BHAT(3)*RK3(I,J)+BHAT(5)*RK4(I,J))
 70        CONTINUE
         ENDIF
 74   CONTINUE
      IF (IFXDPF.EQ.0.AND.(IERRC.EQ.2.OR.IERRC.EQ.20.OR.
     *    IERRC.EQ.21.OR.IERRC.EQ.210)) THEN
         CALL MODGRS(M,N,XHAT,M,BMAT,IFLAG) 
         IF (IFLAG.NE.0) THEN
            IFLAG=10
            RETURN
         ENDIF
      ENDIF

      RETURN
      END
CCCCC END OF HQRK38 CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: XDPNS
C
C*****Purpose
C
C Performs the one step integration of x'=f(x) and at the
C same time (if required) of X'=A(t)X, X0=[I_n 0]' in order
C to approximate the LEs.
C Monitors errors on the one step approximations to the 
C trajectory and/or the LEs and updates approximation to
C the trajectory, and to the LEs, and updates stepsize.
C
C Method used: discrete QR method with the Dormand-Prince 5(4) pair as 
C    integrator for transition matrix X.  Can proceed in fixed or 
C    variable step size.  In the latter case, local error control is 
C    performed on the trajectory or on the Lyapunov exponents
C 
C*****Author: L Dieci 
C*****Date: May 15, 2004
C
CCCCCC XDPNS CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE XDPNS(GETF,GETDF,M,N,T,H,HOLD,Y0,X0,RLYAP,TOLT,
     *          TOLL,AMAT,X,XHAT,RK1,RK2,RK3,RK4,RK5,RK6,STAGE,  
     *          ARK,BHAT,Y,YHAT,YRK1,YRK2,YRK3,YRK4,YRK5,YRK6,
     *          Y1,Y15,Y03,Y45,Y89,
     *          IFXDPF,IERRC,NREJ,IFLAG,INARR,REARR)

      IMPLICIT NONE

      EXTERNAL GETF, GETDF
      INTEGER M,N,IFXDPF,IERRC,NREJ,IFLAG, INARR(*)
      DOUBLE PRECISION T,H,HOLD, Y0(M), REARR(*)
      DOUBLE PRECISION X0(M,N), RLYAP(N), TOLT,TOLL(*)
      DOUBLE PRECISION AMAT(M,M), X(M,N), XHAT(M,N)
      DOUBLE PRECISION RK1(M,N),RK2(M,N)
      DOUBLE PRECISION RK3(M,N),RK4(M,N)
      DOUBLE PRECISION RK5(M,N),RK6(M,N)
      DOUBLE PRECISION STAGE(M,N)     
      DOUBLE PRECISION ARK(7,7), BHAT(7)
      DOUBLE PRECISION Y(M),YHAT(M),YRK1(M),YRK2(M),YRK3(M),YRK4(M)
      DOUBLE PRECISION YRK5(M),YRK6(M)
      DOUBLE PRECISION Y1(M),Y15(M),Y03(M),Y45(M),Y89(M)

      LOGICAL PREVIO, LAST
      INTEGER INTDIR, IFIRST, IWHAT, I, J
      DOUBLE PRECISION HNEW, HMIN, HNEXT, TFIRST
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      DOUBLE PRECISION SCALE, RTOL, ERR
      COMMON /STEP/ HMIN, HNEXT, TFIRST, INTDIR, IFIRST, IWHAT, LAST
      SAVE /STEP/
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

      PREVIO=.FALSE.
C PREVIO IS .FALSE. IF PREVIOUS STEP WAS SUCCESSFUL
C IT IS .TRUE. IF PREVIOUS STEP WAS A FAILURE
 2    CONTINUE
CCCCCCCCCC INTEGRATE FOR X: RETURNS X AND XHAT
      CALL XRKDP(GETF,GETDF,PREVIO,M,N,IFIRST,IWHAT,H,Y0,X0,
     *        AMAT,X,XHAT,RK1,RK2,RK3,RK4,RK5,RK6,
     *        STAGE,ARK,BHAT,
     *        Y,YHAT,YRK1,YRK2,YRK3,YRK4,YRK5,YRK6,
     *        Y1,Y15,Y03,Y45,Y89,IFXDPF,IERRC,INARR,REARR)
      IF (IFIRST.EQ.1) IFIRST=0
C SKIP IF ONLY WANT SOLUTION TRAJECTORY
      IF (IWHAT.EQ.1) GOTO 11
C NOW DO QR FACTORIZATIONS, GET 1-STEP APPROXIMATIONS
C OF THE LES, CHECK ERRORS AND UPDATE APPROXIMATIONS
C DO QR OF X AND PUT 1-STEP APPROXIMATIONS OF LES IN RK2(I,I)
C Q-FACTOR OF X IS IN X ITSELF
      CALL MODGRS(M,N,X,M,RK2,IFLAG) 
C      CALL QLAPQR(M,N,X,M,RK2,STAGE,RK3,IFLAG) 
      IF (IFLAG.NE.0) GOTO 10
      DO 5 I=1,N
         RK2(I,I)=DLOG(RK2(I,I))
 5    CONTINUE
      IF ((IFXDPF.EQ.0.AND.IERRC.EQ.0).OR.(IFXDPF.EQ.1)) GOTO 10
C DO QR OF XHAT AND PUT APPROXIMATE HAT-LES IN STAGE(I,I)
C Q-FACTOR OF XHAT IS IN XHAT
      CALL MODGRS(M,N,XHAT,M,STAGE,IFLAG) 
C      CALL QLAPQR(M,N,XHAT,M,STAGE,RK3,RK4,IFLAG) 
      IF (IFLAG.NE.0) GOTO 10
      DO 7 I=1,N
         STAGE(I,I)=DLOG(STAGE(I,I))
 7    CONTINUE
 10   IF (IFLAG.NE.0) THEN
         IF (IFXDPF.EQ.0) THEN
C TRY DRASTIC WAY TO RECOVER
            IF (DABS(H).GT.HMIN) THEN
               IFLAG=0
               H=INTDIR*HMIN
               PREVIO=.TRUE.
               NREJ=NREJ+1
               GOTO 2
              ELSE
               IFLAG=5
               RETURN
            ENDIF
           ELSE
            IFLAG=10
            RETURN
         END IF
      END IF
 11   CONTINUE
      IF (IFXDPF.EQ.0) THEN
C HERE WE PREDICT NEXT STEPSIZE HNEW, NEVER ALLOWING IT TO
C BE GREATER THAN FIVE TIMES CURRENT H 
C (THIS VALUE CAN BE CHANGED IN HINCR IN "INIT")
         HNEW=HINCR*H
C (A): CHECK ERROR RELATIVE TO TRAJECTORY OR LES, OR BOTH
C (B): PREDICT STEPSIZE RESULTING FROM THIS ERROR ESTIMATE
C (C): IF UNSUCCESSFUL, GO BACK 
         RTOL=ZERO
         IF ((IERRC.EQ.0).OR.(IERRC.EQ.10)) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
C PREDICT STEPSIZE; USE SAFETY FACTOR OF 0.8
C (THIS VALUE CAN BE CHANGED IN HSFTY IN "INIT")
            IF (ERR.GT.ZERO)         
     1         HNEW=DMIN1(DABS(HNEW),DABS(HSFTY*H*(ONE/ERR)**EXP1))
            RTOL=DMAX1(RTOL,ERR)
         ENDIF
         IF (IERRC.NE.0) THEN
            DO 20 I=1,N
C CHECK ERROR I-TH LE: 
               SCALE=ONE+DABS(RK2(I,I))
               ERR=DABS(RK2(I,I)-STAGE(I,I))/SCALE
               ERR=ERR/TOLL(I)
               RTOL=DMAX1(RTOL,ERR)
 20         CONTINUE
            ERR=RTOL
C UPDATE PREDICTION OF STEPSIZE
            IF (ERR.GT.ZERO)         
     1         HNEW=DMIN1(DABS(HNEW),DABS(HSFTY*H*(ONE/ERR)**EXP1))
         ENDIF
C (C) GO BACK IF FAILED UNLESS H TOO SMALL
         IF (ERR.GT.ONE) THEN
C AVOID DRASTIC REDUCTION OF STEPSIZE (NEVER LESS THAN 1/5)
C (THIS VALUE CAN BE CHANGED IN HDRSTC IN "INIT")
            H=INTDIR*DMAX1(DABS(HNEW),HDRSTC*DABS(H))
            NREJ=NREJ+1
            IF (INTDIR*H.LE.HMIN) THEN
               PRINT*,' H = ', H
               IFLAG = 5
               RETURN
            ENDIF
            PREVIO=.TRUE.
            GOTO 2
         ENDIF            
      END IF
C HAVE COMPLETED THE STEP WITH STEPSIZE H
C SAVE LAST SUCCESSFUL STEPSIZE
      HOLD=H
      IF (IFXDPF.EQ.0) THEN
C         HNEW=INTDIR*DMAX1(DABS(HNEW),DABS(H))
         IF((PREVIO).AND.(DABS(HNEW).GT.DABS(H)))HNEW=H
         H=HNEW
      END IF
C FINALLY, UPDATE THE APPROXIMATION FOR THE LES, Q, AND TRAJECTORY
      IF (IWHAT.EQ.1) GOTO 29
      SCALE=ONE/(T+HOLD-TFIRST)      
      DO 25 I=1,N
         RLYAP(I)=((T-TFIRST)*SCALE)*RLYAP(I)+SCALE*RK2(I,I)
C         RLYAP(I)=((T-TFIRST)*SCALE)*RLYAP(I)+SCALE*STAGE(I,I)
 25   CONTINUE         
      DO 28 J=1,N
         DO 28 I=1,M
            X0(I,J)=X(I,J)
 28   CONTINUE
 29   CONTINUE
      DO 32 I=1,M
         Y0(I)=Y(I)
 32   CONTINUE

      RETURN
      END
CCCCC END OF XDPNS CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: XRKDP
C
C*****Purpose
C
C integrates on one step the solution of x'=f(x).
C If needed (see IPAR) it also approximates X'=A(t)X, X0=Q0.
C It returns two approximations, x and xhat, to the solution, 
C two approximations X and XHAT to the X-factor.  From these, 
C one can compute (in HDPNS) two one-step contribution
C to the LEs.  All of these can be used for error purposes.
C 
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC XRKDP CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE XRKDP(GETF,GETDF,PREVIO,M,N,IFIRST,IWHAT,H,Y0,X0,
     *        AMAT,X,XHAT,RK1,RK2,RK3,RK4,RK5,RK6,
     *        STAGE, ARK,BHAT,
     *        Y,YHAT,YRK1,YRK2,YRK3,YRK4,YRK5,YRK6,
     *        Y1,Y15,Y03,Y45,Y89,IFXDPF,IERRC,INARR,REARR)

      IMPLICIT NONE

      EXTERNAL GETF,GETDF
      LOGICAL PREVIO
      INTEGER M, N, IFIRST, IWHAT, IFXDPF,IERRC, INARR(*)
      DOUBLE PRECISION H, Y0(M), X0(M,N), REARR(*)
      DOUBLE PRECISION X(M,N),AMAT(M,M), XHAT(M,N)
      DOUBLE PRECISION RK1(M,N),RK2(M,N)
      DOUBLE PRECISION RK3(M,N),RK4(M,N)
      DOUBLE PRECISION RK5(M,N),RK6(M,N)
      DOUBLE PRECISION STAGE(M,N)
      DOUBLE PRECISION ARK(7,7), BHAT(7)
      DOUBLE PRECISION Y(M),YHAT(M),YRK1(M),YRK2(M),YRK3(M),YRK4(M)
      DOUBLE PRECISION YRK5(M),YRK6(M)
      DOUBLE PRECISION Y1(M),Y15(M),Y03(M),Y45(M),Y89(M)
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      INTEGER I, J
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

C INITIALIZE Y TO OLD Y0
      DO 7 I=1,M
         Y(I)=Y0(I)
 7       YHAT(I)=Y0(I)
C FORM RK1, RK2, RK3, RK4, RK5, RK6, RK7 
C FORM ALSO THE STAGE VALUES AT H/5, 3H/10, 4H/5, 8H/9 AND H (NB THE LATTER IS
C NOT THE SOLUTION VALUE, BUT JUST THE STAGE VALUE)
C ... RK1 AND U15 ...
C VERY FIRST CALL
C      IF (IFIRST.EQ.1) THEN
      IF ((IFIRST.EQ.1).OR.(IFXDPF.EQ.1)) THEN
         CALL GETF(M,Y,YRK1,INARR,REARR) 
         GOTO 8
      END IF
C IF NOT FAILED PREVIOUS STEP, YRK IS IN YRK2
C OTHERWISE, IT IS IN YRK1
      IF (.NOT.PREVIO) THEN
          DO 2 I=1,M 
 2           YRK1(I)=YRK2(I)
      ENDIF
 8    CONTINUE
      DO 11 I=1,M
         Y15(I) = Y(I) + H*ARK(2,1)*YRK1(I)
 11   CONTINUE
C ... RK2 AND U03 ...
      CALL GETF(M,Y15,YRK2,INARR,REARR) 
      DO 22 I=1,M
         Y03(I) = Y(I) + H*(ARK(3,1)*YRK1(I)+ARK(3,2)*YRK2(I))
 22   CONTINUE
C ... RK3 AND U45 
      CALL GETF(M,Y03,YRK3,INARR,REARR) 
      DO 24 I=1,M
         Y45(I) = Y(I) + H*(ARK(4,1)*YRK1(I)+ARK(4,2)*YRK2(I)+
     1                ARK(4,3)*YRK3(I))
 24   CONTINUE
C ... RK4 AND U89 
      CALL GETF(M,Y45,YRK4,INARR,REARR) 
      DO 25 I=1,M
         Y89(I) = Y(I) + H*(ARK(5,1)*YRK1(I)+ARK(5,2)*YRK2(I)+
     1                ARK(5,3)*YRK3(I)+ARK(5,4)*YRK4(I))
 25   CONTINUE
C ... RK5 AND U1 (STAGE VALUE AT H) ...
      CALL GETF(M,Y89,YRK5,INARR,REARR) 
      DO 26 I=1,M
         Y1(I) = Y(I) + H*(ARK(6,1)*YRK1(I)+ARK(6,2)*YRK2(I)+
     1                ARK(6,3)*YRK3(I)+ARK(6,4)*YRK4(I)+
     2                ARK(6,5)*YRK5(I))
 26   CONTINUE
C ... RK6, SOLUTION AT H ...
      CALL GETF(M,Y1,YRK6,INARR,REARR) 
      DO 28 I=1,M
         Y(I) = Y(I) + H*(ARK(7,1)*YRK1(I)+
     1         ARK(7,3)*YRK3(I)+ARK(7,4)*YRK4(I)+
     2         ARK(7,5)*YRK5(I)+ARK(7,6)*YRK6(I))
 28   CONTINUE
      IF (IFXDPF.EQ.1) GOTO 32
C ... RK7 (STORED IN RK2), AND COMPARISON SOLUTION XHAT ...
C TO BE USED FOR ERROR MONITOR/STEP-SIZE SELECTION
      CALL GETF(M,Y,YRK2,INARR,REARR) 
      DO 29 I=1,M 
         YHAT(I)=YHAT(I)+H*(BHAT(1)*YRK1(I)+BHAT(3)*YRK3(I)+
     1       BHAT(4)*YRK4(I)+BHAT(5)*YRK5(I)+BHAT(6)*YRK6(I)+
     2       BHAT(7)*YRK2(I))
 29   CONTINUE
C NOW DO IT FOR THE TRANSITION MATRIX (IF IWHAT=0)
 32   IF (IWHAT.EQ.1) RETURN
C INTIALIZE X TO OLD Q
      DO 5 J=1,N
         DO 5 I=1,M
            X(I,J)=X0(I,J)
            XHAT(I,J)=X0(I,J)
 5    CONTINUE
C FORM RK1, RK2, RK3, RK4, RK5, RK6, RK7 
C FORM ALSO THE STAGE VALUES AT H/5, 3H/10, 4H/5, 8H/9 AND H (NB THE LATTER IS
C NOT THE SOLUTION VALUE, BUT JUST THE STAGE VALUE)
C ... RK1 AND U15 ...
C VERY FIRST CALL NEED AMAT AND RK1
C      IF (IFIRST.EQ.1) THEN
      IF ((IFIRST.EQ.1).OR.(IFXDPF.EQ.0.AND.IERRC.EQ.0).OR.
     *    (IFXDPF.EQ.1)) THEN
         CALL XDOT(GETDF,M,N,Y0,AMAT,X,RK1,INARR,REARR) 
         GOTO 83
      END IF
C OTHERWISE IF NOT FAILED PREVIOUS STEP, HAVE AMAT IN AMAT
C BUT NEED RK1.  IF FAILED PREVIOUS STEP, HAVE RK1 IN RK1
      IF (.NOT.PREVIO) THEN
          CALL FORMAX(M,AMAT,M,M,X,N,M,RK1)
      ENDIF
 83   CONTINUE
      DO 10 J=1,N
         DO 10 I=1,M
            STAGE(I,J) = X(I,J) + H*ARK(2,1)*RK1(I,J)
 10   CONTINUE
C ... RK2 AND U03 ...
      CALL XDOT(GETDF,M,N,Y15,AMAT,STAGE,RK2,INARR,REARR) 
      DO 20 J=1,N
         DO 20 I=1,M
            STAGE(I,J) = X(I,J) + H*(ARK(3,1)*RK1(I,J)+
     1                 ARK(3,2)*RK2(I,J))
 20   CONTINUE
C ... RK3 AND U45 ...
      CALL XDOT(GETDF,M,N,Y03,AMAT,STAGE,RK3,INARR,REARR) 
      DO 30 J=1,N
         DO 30 I=1,M
            STAGE(I,J) = X(I,J) + H*(ARK(4,1)*RK1(I,J)+
     1                ARK(4,2)*RK2(I,J)+ARK(4,3)*RK3(I,J))
 30   CONTINUE
C ... RK4 AND U89 
      CALL XDOT(GETDF,M,N,Y45,AMAT,STAGE,RK4,INARR,REARR) 
      DO 40 J=1,N
         DO 40 I=1,M
            STAGE(I,J) = X(I,J) + H*(ARK(5,1)*RK1(I,J)+
     1                 ARK(5,2)*RK2(I,J)+
     2                 ARK(5,3)*RK3(I,J)+ARK(5,4)*RK4(I,J))
 40   CONTINUE
C ... RK5 AND U1 (STAGE VALUE AT H) ...
      CALL XDOT(GETDF,M,N,Y89,AMAT,STAGE,RK5,INARR,REARR) 
      DO 50 J=1,N
         DO 50 I=1,M
            STAGE(I,J) = X(I,J) + H*(ARK(6,1)*RK1(I,J)+
     1                ARK(6,2)*RK2(I,J)+ARK(6,3)*RK3(I,J)+
     2                ARK(6,4)*RK4(I,J)+ARK(6,5)*RK5(I,J))
 50   CONTINUE
C ... RK6, SOLUTION AT H ...
      CALL XDOT(GETDF,M,N,Y1,AMAT,STAGE,RK6,INARR,REARR) 
      DO 60 J=1,N
         DO 60 I=1,M
            X(I,J) = X(I,J) + H*(ARK(7,1)*RK1(I,J)+
     1         ARK(7,3)*RK3(I,J)+ARK(7,4)*RK4(I,J)+
     2         ARK(7,5)*RK5(I,J)+ARK(7,6)*RK6(I,J))
 60   CONTINUE
C DO NOT FORM HAT-VALUE IF ERR-CONTROL ON TRAJECTORY ONLY
      IF ((IFXDPF.EQ.0.AND.IERRC.EQ.0).OR.(IFXDPF.EQ.1)) RETURN
C ... RK7 (STORED IN RK2), AND COMPARISON SOLUTION XHAT ...
C TO BE USED FOR ERROR MONITOR/STEP-SIZE SELECTION
      CALL XDOT(GETDF,M,N,Y,AMAT,X,RK2,INARR,REARR) 
      DO 70 J=1,N 
         DO 70 I=1,M
            XHAT(I,J)=XHAT(I,J)+H*(BHAT(1)*RK1(I,J)+
     1       BHAT(3)*RK3(I,J)+BHAT(4)*RK4(I,J)+
     2       BHAT(5)*RK5(I,J)+BHAT(6)*RK6(I,J)+BHAT(7)*RK2(I,J))
 70   CONTINUE

      RETURN
      END
CCCCC END OF XRKDP CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: X38NS
C
C*****Purpose
C
C Performs the one step integration of x'=f(x) and at the
C same time (if required) of X'=A(t)X, X0=[I_n 0]' in order
C to approximate the LEs.
C Monitors errors on the one step approximations to the 
C trajectory and/or the LEs and updates approximation to
C the trajectory, and to the LEs, and updates stepsize.
C
C Method used: discrete QR method with the RK38 4(3) pair as 
C    integrator for transition matrix X.  Can proceed in fixed or 
C    variable step size.  In the latter case, local error control is 
C    performed on the trajectory or on the Lyapunov exponents
C 
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC X38NS CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE X38NS(GETF,GETDF,M,N,T,H,HOLD,Y0,X0,RLYAP,TOLT,
     *          TOLL,AMAT,X,XHAT,RK1,RK2,RK3,RK4,STAGE,ARK,BHAT,
     *          Y,YHAT,YRK1,YRK2,YRK3,YRK4,Y1,Y13,Y23,
     *          IFXDPF,IERRC,NREJ,IFLAG,INARR,REARR)

      IMPLICIT NONE

      EXTERNAL GETF, GETDF
      INTEGER M,N,IFXDPF,IERRC,NREJ,IFLAG, INARR(*)
      DOUBLE PRECISION T,H,HOLD, Y0(M), REARR(*)
      DOUBLE PRECISION X0(M,N), RLYAP(N), TOLT, TOLL(*)
      DOUBLE PRECISION AMAT(M,M), X(M,N), XHAT(M,N)
      DOUBLE PRECISION RK1(M,N),RK2(M,N)
      DOUBLE PRECISION RK3(M,N),RK4(M,N)
      DOUBLE PRECISION STAGE(M,N)
      DOUBLE PRECISION ARK(7,7), BHAT(7)
      DOUBLE PRECISION Y(M),YHAT(M),YRK1(M),YRK2(M),YRK3(M),YRK4(M)
      DOUBLE PRECISION Y1(M),Y13(M),Y23(M)

      LOGICAL PREVIO, LAST
      INTEGER INTDIR, IFIRST, IWHAT, I, J
      DOUBLE PRECISION HNEW, HMIN, HNEXT, TFIRST
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      DOUBLE PRECISION SCALE, RTOL, ERR
      COMMON /STEP/ HMIN, HNEXT, TFIRST, INTDIR, IFIRST, IWHAT, LAST
      SAVE /STEP/
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

      PREVIO=.FALSE.
C PREVIO IS .FALSE. IF PREVIOUS STEP WAS SUCCESSFUL
C IT IS .TRUE. IF PREVIOUS STEP WAS A FAILURE
 2    CONTINUE
CCCCCCCCCC INTEGRATE FOR X: RETURNS X AND XHAT
      CALL XRK38(GETF,GETDF,PREVIO,M,N,IFIRST,IWHAT,H,Y0,X0,
     *        AMAT,X,XHAT,RK1,RK2,RK3,RK4,STAGE,ARK,BHAT,
     *        Y,YHAT,YRK1,YRK2,YRK3,YRK4,Y1,Y13,Y23,
     *        IFXDPF,IERRC,INARR,REARR)
      IF (IFIRST.EQ.1) IFIRST=0
C SKIP IF ONLY WANT SOLUTION TRAJECTORY
      IF (IWHAT.EQ.1) GOTO 11
C NOW DO QR FACTORIZATIONS, GET 1-STEP APPROXIMATIONS
C OF THE LES, CHECK ERRORS AND UPDATE APPROXIMATIONS
C DO QR OF X AND PUT 1-STEP APPROXIMATIONS OF LES IN RK2(I,I)
C Q-FACTOR OF X IS IN X ITSELF
      CALL MODGRS(M,N,X,M,RK2,IFLAG) 
C      CALL QLAPQR(M,N,X,M,RK2,STAGE,RK3,IFLAG) 
      IF (IFLAG.NE.0) GOTO 10
      DO 5 I=1,N
         RK2(I,I)=DLOG(RK2(I,I))
 5    CONTINUE
      IF ((IFXDPF.EQ.0.AND.IERRC.EQ.0).OR.(IFXDPF.EQ.1)) GOTO 10
C DO QR OF XHAT AND PUT APPROXIMATE HAT-LES IN STAGE(I,I)
C Q-FACTOR OF XHAT IS IN XHAT
      CALL MODGRS(M,N,XHAT,M,STAGE,IFLAG) 
C      CALL QLAPQR(M,N,XHAT,M,STAGE,RK3,RK4,IFLAG) 
      IF (IFLAG.NE.0) GOTO 10
      DO 7 I=1,N
         STAGE(I,I)=DLOG(STAGE(I,I))
 7    CONTINUE
 10   IF (IFLAG.NE.0) THEN
         IF (IFXDPF.EQ.0) THEN
C TRY DRASTIC WAY TO RECOVER
            IF (DABS(H).GT.HMIN) THEN
               IFLAG=0
               H=INTDIR*HMIN
               PREVIO=.TRUE.
               NREJ=NREJ+1
               GOTO 2
              ELSE
               IFLAG=5
               RETURN
            ENDIF
           ELSE
            IFLAG=10
            RETURN
         END IF
      END IF
 11   CONTINUE
      IF (IFXDPF.EQ.0) THEN
C HERE WE PREDICT NEXT STEPSIZE HNEW, NEVER ALLOWING IT TO
C BE GREATER THAN FIVE TIMES CURRENT H 
C (THIS VALUE CAN BE CHANGED IN HINCR IN "INIT")
         HNEW=HINCR*H
C (A): CHECK ERROR RELATIVE TO TRAJECTORY OR LES, OR BOTH
C (B): PREDICT STEPSIZE RESULTING FROM THIS ERROR ESTIMATE
C (C): IF UNSUCCESSFUL, GO BACK 
         RTOL=ZERO
         IF ((IERRC.EQ.0).OR.(IERRC.EQ.10)) THEN
            CALL NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)
C PREDICT STEPSIZE; USE SAFETY FACTOR OF 0.8
C (THIS VALUE CAN BE CHANGED IN HSFTY IN "INIT")
            IF (ERR.GT.ZERO)         
     1         HNEW=DMIN1(DABS(HNEW),DABS(HSFTY*H*(ONE/ERR)**EXP2))
            RTOL=DMAX1(RTOL,ERR)
         ENDIF
         IF (IERRC.NE.0) THEN
            DO 20 I=1,N
C CHECK ERROR I-TH LE: 
               SCALE=ONE+DABS(RK2(I,I))
               ERR=DABS(RK2(I,I)-STAGE(I,I))/SCALE
               ERR=ERR/TOLL(I)
               RTOL=DMAX1(RTOL,ERR)
 20         CONTINUE
            ERR=RTOL
C UPDATE PREDICTION OF STEPSIZE
            IF (ERR.GT.ZERO)         
     1         HNEW=DMIN1(DABS(HNEW),DABS(HSFTY*H*(ONE/ERR)**EXP2))
         ENDIF
C (C) GO BACK IF FAILED UNLESS H TOO SMALL
         IF (ERR.GT.ONE) THEN
C AVOID DRASTIC REDUCTION OF STEPSIZE (NEVER LESS THAN 1/5)
C (THIS VALUE CAN BE CHANGED IN HDRSTC IN "INIT")
            H=INTDIR*DMAX1(DABS(HNEW),HDRSTC*DABS(H))
            NREJ=NREJ+1
            IF (INTDIR*H.LE.HMIN) THEN
               PRINT*,' H = ', H
               IFLAG = 5
               RETURN
            ENDIF
            PREVIO=.TRUE.
            GOTO 2
         ENDIF            
      END IF
C HAVE COMPLETED THE STEP WITH STEPSIZE H
C SAVE LAST SUCCESSFUL STEPSIZE
      HOLD=H
      IF (IFXDPF.EQ.0) THEN
C         HNEW=INTDIR*DMAX1(DABS(HNEW),DABS(H))
         IF((PREVIO).AND.(DABS(HNEW).GT.DABS(H)))HNEW=H
         H=HNEW
      END IF
C FINALLY, UPDATE THE APPROXIMATION FOR THE LES, Q, AND TRAJECTORY
      IF (IWHAT.EQ.1) GOTO 29
      SCALE=ONE/(T+HOLD-TFIRST)      
      DO 25 I=1,N
         RLYAP(I)=((T-TFIRST)*SCALE)*RLYAP(I)+SCALE*RK2(I,I)
C         RLYAP(I)=((T-TFIRST)*SCALE)*RLYAP(I)+SCALE*STAGE(I,I)
 25   CONTINUE         
      DO 28 J=1,N
         DO 28 I=1,M
            X0(I,J)=X(I,J)
 28   CONTINUE
 29   CONTINUE
      DO 32 I=1,M
         Y0(I)=Y(I)
 32   CONTINUE

      RETURN
      END
CCCCC END OF X38NS CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: XRK38
C
C*****Purpose
C
C integrates on one step the solution of x'=f(x).
C If needed (see IPAR) it also approximates X'=A(t)X, X0=Q0.
C It returns two approximations, x and xhat, to the solution, 
C two approximations X and XHAT to the X-factor.  From these, 
C one can compute (in HDPNS) two one-step contribution
C to the LEs.  All of these can be used for error purposes.
C 
C*****Author: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC XRK38 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE XRK38(GETF,GETDF,PREVIO,M,N,IFIRST,IWHAT,H,Y0,X0,
     *        AMAT,X,XHAT,RK1,RK2,RK3,RK4,STAGE,ARK,BHAT,
     *        Y,YHAT,YRK1,YRK2,YRK3,YRK4,Y1,Y13,Y23,
     *        IFXDPF,IERRC,INARR,REARR)

      IMPLICIT NONE

      EXTERNAL GETF,GETDF
      LOGICAL PREVIO
      INTEGER M, N, IFIRST, IWHAT, IFXDPF,IERRC, INARR(*)
      DOUBLE PRECISION H, Y0(M), X0(M,N), REARR(*)
      DOUBLE PRECISION X(M,N),AMAT(M,M), XHAT(M,N)
      DOUBLE PRECISION RK1(M,N),RK2(M,N),RK3(M,N),RK4(M,N)
      DOUBLE PRECISION STAGE(M,N)
      DOUBLE PRECISION ARK(7,7), BHAT(7)
      DOUBLE PRECISION Y(M),YHAT(M),YRK1(M),YRK2(M),YRK3(M),YRK4(M)
      DOUBLE PRECISION Y1(M),Y13(M),Y23(M)
      DOUBLE PRECISION ONE, ZERO, EXP1,EXP2
      DOUBLE PRECISION HINCR, HSFTY, HDRSTC,TEPS
      INTEGER I, J
      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

C INITIALIZE Y TO OLD Y0
      DO 7 I=1,M
         Y(I)=Y0(I)
 7       YHAT(I)=Y0(I)
C FORM RK1, RK2, RK3, RK4 AND RK5 
C FORM ALSO THE STAGE VALUES AT H/3, 2H/3 AND H (NB THE LATTER IS
C NOT THE SOLUTION VALUE, BUT JUST THE STAGE VALUE)
C ... RK1 AND U13 ...
C VERY FIRST CALL
C      IF (IFIRST.EQ.1) THEN
      IF ((IFIRST.EQ.1).OR.(IFXDPF.EQ.1)) THEN
         CALL GETF(M,Y,YRK1,INARR,REARR) 
         GOTO 8
      END IF
C IF NOT FAILED PREVIOUS STEP, YRK IS IN YRK4
C OTHERWISE, IT IS IN YRK1
      IF (.NOT.PREVIO) THEN
          DO 2 I=1,M 
 2           YRK1(I)=YRK4(I)
      ENDIF
 8    CONTINUE
      DO 11 I=1,M
         Y13(I) = Y(I) + H*ARK(2,1)*YRK1(I)
 11   CONTINUE
C ... RK2 AND U23 ...
      CALL GETF(M,Y13,YRK2,INARR,REARR) 
      DO 22 I=1,M
         Y23(I) = Y(I) + H*(ARK(3,1)*YRK1(I)+ARK(3,2)*YRK2(I))
 22   CONTINUE
C ... RK3 AND U1 (STAGE VALUE AT H) ...
      CALL GETF(M,Y23,YRK3,INARR,REARR) 
      DO 24 I=1,M
         Y1(I) = Y(I) + H*(ARK(4,1)*YRK1(I)+ARK(4,2)*YRK2(I)+
     1                ARK(4,3)*YRK3(I))
 24   CONTINUE
C ... RK4, SOLUTION AT H ...
      CALL GETF(M,Y1,YRK4,INARR,REARR) 
      DO 26 I=1,M
         Y(I) = Y(I) + H*(ARK(5,1)*YRK1(I)+
     1         ARK(5,2)*YRK2(I)+ARK(5,3)*YRK3(I)+ARK(5,4)*YRK4(I))
 26   CONTINUE
      IF (IFXDPF.EQ.1) GOTO 32
C ... RK5 (STORED IN RK4), AND COMPARISON SOLUTION XHAT ...
C TO BE USED FOR ERROR MONITOR/STEP-SIZE SELECTION
      CALL GETF(M,Y,YRK4,INARR,REARR) 
      DO 29 I=1,M 
         YHAT(I)=YHAT(I)+H*(BHAT(1)*YRK1(I)+
     1       BHAT(2)*YRK2(I)+BHAT(3)*YRK3(I)+BHAT(5)*YRK4(I))
 29   CONTINUE
C NOW DO IT FOR THE TRANSITION MATRIX (IF IWHAT=0)
 32   IF (IWHAT.EQ.1) RETURN
C INTIALIZE X TO OLD Q
      DO 5 J=1,N
         DO 5 I=1,M
            X(I,J)=X0(I,J)
            XHAT(I,J)=X0(I,J)
 5    CONTINUE
C FORM RK1, RK2, RK3, RK4 AND RK5 
C FORM ALSO THE STAGE VALUES AT H/3, 2H/3 AND H (NB THE LATTER IS
C NOT THE SOLUTION VALUE, BUT JUST THE STAGE VALUE)
C ... RK1 AND U13 ...
C VERY FIRST CALL NEED AMAT AND RK1
C      IF (IFIRST.EQ.1) THEN
      IF ((IFIRST.EQ.1).OR.(IFXDPF.EQ.0.AND.IERRC.EQ.0).OR.
     *    (IFXDPF.EQ.1)) THEN
         CALL XDOT(GETDF,M,N,Y0,AMAT,X,RK1,INARR,REARR) 
         GOTO 83
      END IF
C OTHERWISE IF NOT FAILED PREVIOUS STEP, HAVE AMAT IN AMAT
C BUT NEED RK1.  IF FAILED PREVIOUS STEP, HAVE RK1 IN RK1
      IF (.NOT.PREVIO) THEN
          CALL FORMAX(M,AMAT,M,M,X,N,M,RK1)
      ENDIF
 83   CONTINUE
      DO 10 J=1,N
         DO 10 I=1,M
            STAGE(I,J) = X(I,J) + H*ARK(2,1)*RK1(I,J)
 10   CONTINUE
C ... RK2 AND U23 ...
      CALL XDOT(GETDF,M,N,Y13,AMAT,STAGE,RK2,INARR,REARR) 
      DO 20 J=1,N
         DO 20 I=1,M
            STAGE(I,J) = X(I,J) + H*(ARK(3,1)*RK1(I,J)+
     1                 ARK(3,2)*RK2(I,J))
 20   CONTINUE
C ... RK3 AND U1 (STAGE VALUE AT H) ...
      CALL XDOT(GETDF,M,N,Y23,AMAT,STAGE,RK3,INARR,REARR) 
      DO 23 J=1,N
         DO 23 I=1,M
            STAGE(I,J) = X(I,J) + H*(ARK(4,1)*RK1(I,J)+
     1                ARK(4,2)*RK2(I,J)+ARK(4,3)*RK3(I,J))
 23   CONTINUE
C ... RK4, SOLUTION AT H ...
      CALL XDOT(GETDF,M,N,Y1,AMAT,STAGE,RK4,INARR,REARR) 
      DO 25 J=1,N
         DO 25 I=1,M
            X(I,J) = X(I,J) + H*(ARK(5,1)*RK1(I,J)+
     1         ARK(5,2)*RK2(I,J)+ARK(5,3)*RK3(I,J)+ARK(5,4)*RK4(I,J))
 25   CONTINUE
      IF ((IFXDPF.EQ.0.AND.IERRC.EQ.0).OR.(IFXDPF.EQ.1)) RETURN
C ... RK5 (STORED IN RK4), AND COMPARISON SOLUTION XHAT ...
C TO BE USED FOR ERROR MONITOR/STEP-SIZE SELECTION
      CALL XDOT(GETDF,M,N,Y,AMAT,X,RK4,INARR,REARR) 
      DO 28 J=1,N 
         DO 28 I=1,M
            XHAT(I,J)=XHAT(I,J)+H*(BHAT(1)*RK1(I,J)+
     1       BHAT(2)*RK2(I,J)+BHAT(3)*RK3(I,J)+BHAT(5)*RK4(I,J))
 28   CONTINUE

      RETURN
      END
CCCCC END OF XRK38 CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: NERRQ
C
C*****Purpose
C     Computes max-element norm of errors of Y by columns
C     using a mixed absolute/relative error criterion
C
C*****Authors: L Dieci 
C*****Date: April 23, 2004
C
CCCCCCC NERRQ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NERRQ(M,N,Q,QHAT,Y0,ERR,TOLQ,ONE,ZERO)

      IMPLICIT NONE

      INTEGER M,N
      DOUBLE PRECISION Q(M,N), QHAT(M,N), Y0(M,N), TOLQ, ERR
      DOUBLE PRECISION ONE, ZERO

      DOUBLE PRECISION RTOL, SCALE
      INTEGER I,J

C Computes the weighted norm of the error for step-size selection
C Uses the 2-norm for each column and then the max of these errors

      RTOL=ZERO
      DO 20 J=1,N
         ERR=ZERO
C WEIGH AGAINST LARGEST OF INITIAL AND FINAL VALUES
         DO 10 I=1,M
            SCALE=ONE+DMAX1(DABS(Q(I,J)),DABS(Y0(I,J)))
            ERR=DMAX1(ERR,DABS(Q(I,J)-QHAT(I,J))/SCALE)
 10      CONTINUE
         ERR=ERR/TOLQ
         RTOL=DMAX1(RTOL,ERR)
 20   CONTINUE
      ERR=RTOL

      RETURN 
      END
CCCCC END OF NERRQ CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: NERRNU
C
C*****Purpose
C     Computes supnorm of errors of nu (integrals)
C     using a mixed absolute/relative error criterion
C
C*****Authors: L Dieci 
C*****Date: April 23, 2004
C
CCCCCCC NERRNU CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NERRNU(N,RNU,RNUHAT,ERR,TOL,ONE,ZERO)

      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION RNU(N), RNUHAT(N), TOL(N), ERR
      DOUBLE PRECISION ONE, ZERO

      DOUBLE PRECISION RTOL, SCALE
      INTEGER I

C Computes the weighted sup-norm of the error for step-size selection

      RTOL=ZERO
      DO 10 I=1,N
C RECALL THAT INITIAL VALUE IS 0
         SCALE=ONE+DABS(RNU(I))
         ERR=DABS(RNU(I)-RNUHAT(I))/SCALE
         ERR=ERR/TOL(I)
         RTOL=DMAX1(RTOL,ERR)
 10   CONTINUE
      ERR=RTOL

      RETURN 
      END
CCCCC END OF NERRNU CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: NERRY
C
C*****Purpose
C     Computes weighted supnorm of errors of solution
C     using a mixed absolute/relative error criterion
C
C*****Authors: L Dieci 
C*****Date: April 23, 2004
C
CCCCCCC NERRNU CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NERRY(M,Y,YHAT,Y0,ERR,TOLT,ONE,ZERO)

      IMPLICIT NONE

      INTEGER M
      DOUBLE PRECISION Y(M), YHAT(M), Y0(M), TOLT, ERR
      DOUBLE PRECISION ONE, ZERO

      DOUBLE PRECISION RTOL, SCALE
      INTEGER I

C Computes the weighted sup-norm of the error for step-size selection

      RTOL=ZERO
      DO 10 I=1,M
C SCALE WITH RESPECT TO MAX AT BEGINNING AND END OF STEP
         SCALE=ONE+DMAX1(DABS(Y(I)),DABS(Y0(I)))
         ERR=DABS(Y(I)-YHAT(I))/SCALE
         ERR=ERR/TOLT
         RTOL=DMAX1(RTOL,ERR)
 10   CONTINUE
      ERR=RTOL

      RETURN 
      END
CCCCC END OF NERRY CCCCC
