CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: TNORM
C
C*****Purpose
C
C computes the 2-norm of the vector X
C input: 
C     N, integer (dimension of X)
C     X, vector whose norm is sought
C     N and X are untouched
C output:
C     RNORM, the 2-norm of X
C        should monitor the value of RNORM on return
C        if RNORM.le.2*eps vector has all 0's to working precision
C
C*****Authors: L Dieci & ES Van Vleck
C*****Date: May 29, 2003
C
CCCCCC TNORM CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TNORM(N,X,RNORM)

      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION X(N), RNORM

      INTEGER I
      DOUBLE PRECISION XMAX
      DOUBLE PRECISION ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS

      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

      IF (N.LE.0) RETURN
      RNORM=ZERO
      XMAX=DABS(X(1))
      DO 2 I=2,N
         XMAX=DMAX1(DABS(X(I)),XMAX)
 2    CONTINUE
      IF (XMAX.EQ.ZERO) RETURN
      DO 10 I=1,N    
         RNORM=RNORM+(DABS(X(I))/XMAX)**2
 10   CONTINUE
      RNORM=DSQRT(RNORM)
      RNORM=RNORM*XMAX
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: MODGRS 
C
C*****Purpose
C
C computes the (modified) Gram-Schmidt factorization 
C of the full rank matrix A
C input: 
C     M, integer (number of rows of A)
C     N, integer (number of columns of A): MUST HAVE M.ge.N
C     A, double precision (matrix to be orthogonalized)
C     LDA, integer (leading dimension of A): MUST HAVE LDA.GE.M
C     R, double precision array of size M*N (work space)
C output:
C     A, the orthogonal factor in the QR of A
C     IFLAG, integer (it is = 0 on successful exit,
C                it is = k if k-th column had norm 0
C                during orthogonalization process)
C
C*****Authors: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC MODGRS CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MODGRS(M,N,A,LDA,R,IFLAG)

      IMPLICIT NONE

      INTEGER M, N, LDA, IFLAG
      DOUBLE PRECISION A(LDA,N), R(M,N)

      INTEGER I, J, K
      DOUBLE PRECISION RTNORM

      DOUBLE PRECISION ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS

      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

      IFLAG=0
      DO 10 K=1,N
         CALL TNORM(M,A(1,K),RTNORM)
         IF(RTNORM.LE.TEPS)THEN
            IFLAG=K
            RETURN
           ELSE 
            R(K,K)=RTNORM
         ENDIF
         DO 30 I=1,M
 30         A(I,K) = A(I,K)/R(K,K)
         DO 40 J=K+1,N
            R(K,J) = ZERO
            DO 50 I=1,M
 50            R(K,J) = R(K,J) + A(I,K)*A(I,J)
            DO 60 I=1,M
 60            A(I,J) = A(I,J) - R(K,J)*A(I,K)
 40      CONTINUE
 10   CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: QLAPQR 
C
C*****Purpose
C
C computes the Q factor in QR factorization of full 
C rank matrix A, using routines from LAPACK and signs' monitoring
C input: 
C     M, integer (number of rows of A)
C     N, integer (number of columns of A): MUST HAVE M.ge.N
C     A, double precision (matrix to be orthogonalized)
C     LDA, integer (leading dimension of A): MUST HAVE LDA.GE.M
C     R, double precision array of size M*N (work space)
C     TAU, double precision array of size N (work space)
C     WORK, double precision array of size N (work space)
C output:
C     A, the orthogonal factor in the QR of A
C     A, on the diagonal it has the (positive) diagonal of
C        the triangular factor in the QR of A
C     IFLAG, integer (it is = 0 on successful exit,
C                it is = k if k-th column had norm 0
C                during orthogonalization process
C
C LAPACK: DGEQRF, DORGQR
C
C*****Authors: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC QLAPQR CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE QLAPQR(M,N,A,LDA,R,TAU,WORK,IFLAG)

      IMPLICIT NONE

      INTEGER M, N, LDA, IFLAG
      DOUBLE PRECISION A(LDA,N), R(M,N), TAU(*), WORK(*)

      INTEGER I, J, NN
      DOUBLE PRECISION ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS

      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

      NN=N
      IF (N.EQ.M) NN=NN-1
      IFLAG=0
      CALL DGEQRF(M,N,A,LDA,TAU,WORK,M*N,IFLAG)
      IF (IFLAG.NE.0) RETURN
C STORE UPPER PART IN R FOR LATER SIGNS' MONITORING
      DO 7 J=1,N
         DO 7 I=J,N
 7          R(I,J)=A(I,J)
C BUILD ORTHOGONAL Q AS PRODUCT OF REFLECTORS INTO A
      CALL DORGQR(M,N,NN,A,LDA,TAU,WORK,M*N,IFLAG)
      IF (IFLAG.NE.0) RETURN
C ADJUST THE SIGNS OF COLUMNS OF Q AS NEEDED
      DO 11 J=1,N
         IF (R(J,J).LT.ZERO) THEN
            R(J,J)=-ONE*R(J,J)
            DO 10 I=1,M
 10            A(I,J)=-ONE*A(I,J)
         ENDIF
 11   CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: TRAPRL 
C
C*****Purpose
C
C takes one step of trapezoidal rule approximation for  
C approximating an integral
C input: 
C     FA and FB, double precision 
C                function values at end points [a,b]
C     H, double precision (it is (b-a))
C output:
C     APP, approximation by quadrature rule
C
C*****Authors: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC TRAPRL CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TRAPRL(FA,FB,H,APP)

      IMPLICIT NONE

      DOUBLE PRECISION FA,FB,H,APP, HALF
      PARAMETER (HALF=0.5D0)

C TRAP-RULE APPROXIMATION 
      APP=FA+FB
      APP=HALF*H*APP

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: DGQTAQ 
C
C*****Purpose
C
C computes the diagonal of Q^(tr)AQ using BLAS routines 
C input: 
C     LDQ, integer (leading dimension of Q)
C     M, integer (number of rows of Q)
C     N, integer (number of columns of Q): MUST HAVE M.ge.N
C     Q, double precision orthogonal (M,N) matrix, unchanged on output
C     A, double precision (square (M,M) matrix), unchanged on output
C     LDA, integer (leading dimension of A): MUST HAVE LDA.ge.M
C     C, double precision matrix of size at least (M,N) (work array)
C     LDC, integer (leading dimension of C): MUST HAVE LDC.ge.M
C output:
C     C, in its N diagonal elements contains the diagonal of the product 
C
C Blas: Blas1: DDOT, Blas3: DGEMM
C
C*****Authors: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC DGQTAQ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DGQTAQ(LDQ,M,N,Q,A,LDA,C,LDC)

      IMPLICIT NONE

      INTEGER M, N, LDQ, LDA, LDC
      DOUBLE PRECISION Q(LDQ,N), A(LDA,M), C(LDC,N)
      DOUBLE PRECISION DDOT

      INTEGER K, INC
      PARAMETER (INC = 1 )
      DOUBLE PRECISION ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS

      COMMON /CONSTS/ONE,ZERO,EXP1,EXP2,HINCR,HSFTY,HDRSTC,TEPS
      SAVE /CONSTS/

C FORM A*Q IN C
      CALL DGEMM ('N', 'N', M, N, M, ONE, A, LDA, Q, LDQ,
     $                   ZERO, C, LDC )
C FORM DIAG(Q'*C) IN DIAGONAL OF C
      DO 5 K=1,N
         C(K,K)=DDOT(M,Q(1,K),INC,C(1,K),INC)
 5    CONTINUE        

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****Subroutine name: DGQTB 
C
C*****Purpose
C
C computes the diagonal of Q^(tr)B using BLAS routines 
C input: 
C     LDQ, integer (leading dimension of Q)
C     M, integer (number of rows of Q)
C     N, integer (number of columns of Q): MUST HAVE M.ge.N
C     Q, double precision orthogonal (M,N) matrix, unchanged on output
C     B, double precision (M,N) matrix), unchanged on output
C     LDB, integer (leading dimension of B): MUST HAVE LDB.ge.M
C     C, double precision matrix of size at least (N,N) (work array)
C     LDC, integer (leading dimension of C): MUST HAVE LDC.ge.N
C output:
C     C, in its N diagonal elements contains the diagonal of the product 
C
C Blas: Blas1: DDOT
C
C*****Authors: L Dieci 
C*****Date: May 29, 2003
C
CCCCCC DGQTB CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DGQTB(LDQ,M,N,Q,B,LDB,C,LDC)

      IMPLICIT NONE

      INTEGER M, N, LDQ, LDB, LDC
      DOUBLE PRECISION Q(LDQ,N), B(LDB,N), C(LDC,N)
      DOUBLE PRECISION DDOT

      INTEGER K, INC
      PARAMETER (INC = 1 )

C FORM DIAG(Q'*B) IN DIAGONAL OF C
      DO 5 K=1,N
         C(K,K)=DDOT(M,Q(1,K),INC,B(1,K),INC)
 5    CONTINUE        

      RETURN
      END
