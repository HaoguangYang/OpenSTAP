SUBROUTINE JACOBI (A,B,X,EIGV,D,N,NWA,RTOL,NSMAX,IFPR,IOUT)
! ....................................................................
! .
! . P R O G R A M
! . TO SOLVE THE GENERALIZED EIGENPROBLEM USING THE
! . GENERALIZED JACOBI ITERATION
! ....................................................................
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: N,NWA,NSMAX,IFPR,IOUT
REAL(8), INTENT(INOUT) :: A(NWA),B(NWA),X(N,N),EIGV(N),D(N)
REAL(8), INTENT(IN) :: RTOL
!
INTEGER :: N1,I,J,K,II,IK,NSWEEP,NR,JP1,JM1,KP1,KM1,JK,KK,LJK,JJ
INTEGER :: IM1,IJ,JI,KI,LKI,LJI
REAL(8) :: AK,AKK,AJ,BJ,AJJ,AB,BB,BK,EPS,EPSA,EPSB,EPTOLA,EPTOLB
REAL(8) :: ABCH,AKKCH,AJJCH,CHECK,SQCH,D1,D2,DEN,CA,CG,DIF,TOL
REAL(8) :: SCALE,XK,XJ
!
! INITIALIZE EIGENVALUE AND EIGENVECTOR MATRICES
!
    N1=N + 1
    II=1
!
    DO I=1,N
      IF (A(II).LE.0. .OR. B(II).LE.0.) THEN
        WRITE (IOUT,2020) II,A(II),B(II)
        STOP
      END IF
      D(I)=A(II)/B(II)
      EIGV(I)=D(I)
      II=II + N1 - I
    END DO
!
    DO I=1,N
      DO J=1,N
        X(I,J)=0.
      END DO
      X(I,I)=1.
    END DO
!
  IF (N.EQ.1) RETURN
!
! INITIALIZE SWEEP COUNTER AND BEGIN ITERATION
!
    NSWEEP=0
    NR=N - 1
!
    DO WHILE (NSWEEP.LT.NSMAX)
      NSWEEP=NSWEEP + 1
      IF (IFPR.EQ.1) WRITE (IOUT,2000) NSWEEP
!
! CHECK IF PRESENT OFF-DIAGONAL ELEMENT IS LARGE ENOUGH TO REQUIRE
! ZEROING
!
      EPS=(.01)**(NSWEEP*2)
      DO J=1,NR
        JP1=J + 1
        JM1=J - 1
        LJK=JM1*N - JM1*J/2
        JJ=LJK + J
        DO K=JP1,N
          KP1=K + 1
          KM1=K - 1
          JK=LJK + K
          KK=KM1*N - KM1*K/2 + K
          EPTOLA=(A(JK)/A(JJ))*(A(JK)/A(KK))
          EPTOLB=(B(JK)/B(JJ))*(B(JK)/B(KK))
          IF (EPTOLA.LT.EPS .AND. EPTOLB.LT.EPS) CYCLE
!
! IF ZEROING IS REQUIRED, CALCULATE THE ROTATION MATRIX ELEMENTS CA
! AND CG
!
          AKK=A(KK)*B(JK) - B(KK)*A(JK)
          AJJ=A(JJ)*B(JK) - B(JJ)*A(JK)
          AB=A(JJ)*B(KK) - A(KK)*B(JJ)
          SCALE=A(KK)*B(KK)
          ABCH=AB/SCALE
          AKKCH=AKK/SCALE
          AJJCH=AJJ/SCALE
          CHECK=(ABCH*ABCH+4.0*AKKCH*AJJCH)/4.0
          IF (CHECK.LT.0) THEN
            WRITE (IOUT,2020) JJ,A(JJ),B(JJ)
            STOP
          END IF
          SQCH=SCALE*SQRT(CHECK)
          D1=AB/2. + SQCH
          D2=AB/2. - SQCH
          DEN=D1
          IF (ABS(D2).GT.ABS(D1)) DEN=D2
          IF (DEN.EQ.0) THEN
            CA=0.
            CG=-A(JK)/A(KK)
          ELSE
            CA=AKK/DEN
            CG=-AJJ/DEN
          END IF
!
! PERFORM THE GENERALIZED ROTATION TO ZERO THE PRESENT OFF-DIAGONAL
! ELEMENT
!
          IF (N-2.NE.0) THEN
            IF (JM1-1.GE.0) THEN
              DO I=1,JM1
                IM1=I - 1
                IJ=IM1*N - IM1*I/2 + J
                IK=IM1*N - IM1*I/2 + K
                AJ=A(IJ)
                BJ=B(IJ)
                AK=A(IK)
                BK=B(IK)
                A(IJ)=AJ + CG*AK
                B(IJ)=BJ + CG*BK
                A(IK)=AK + CA*AJ
                B(IK)=BK + CA*BJ
              END DO
            END IF
            IF (KP1-N.LE.0) THEN
              LJI=JM1*N - JM1*J/2
              LKI=KM1*N - KM1*K/2
              DO I=KP1,N
                JI=LJI + I
                KI=LKI + I
                AJ=A(JI)
                BJ=B(JI)
                AK=A(KI)
                BK=B(KI)
                A(JI)=AJ + CG*AK
                B(JI)=BJ + CG*BK
                A(KI)=AK + CA*AJ
                B(KI)=BK + CA*BJ
              END DO
            END IF
            IF (JP1-KM1.LE.0) THEN
              LJI=JM1*N - JM1*J/2
              DO I=JP1,KM1
                JI=LJI + I
                IM1=I - 1
                IK=IM1*N - IM1*I/2 + K
                AJ=A(JI)
                BJ=B(JI)
                AK=A(IK)
                BK=B(IK)
                A(JI)=AJ + CG*AK
                B(JI)=BJ + CG*BK
                A(IK)=AK + CA*AJ
                B(IK)=BK + CA*BJ
              END DO
            END IF
          END IF
          AK=A(KK)
          BK=B(KK)
          A(KK)=AK + 2.*CA*A(JK) + CA*CA*A(JJ)
          B(KK)=BK + 2.*CA*B(JK) + CA*CA*B(JJ)
          A(JJ)=A(JJ) + 2.*CG*A(JK) + CG*CG*AK
          B(JJ)=B(JJ) + 2.*CG*B(JK) + CG*CG*BK
          A(JK)=0.
          B(JK)=0.
!
! UPDATE THE EIGENVECTOR MATRIX AFTER EACH ROTATION
!
          DO I=1,N
            XJ=X(I,J)
            XK=X(I,K)
            X(I,J)=XJ + CG*XK
            X(I,K)=XK + CA*XJ
          END DO
        END DO
      END DO
!
! UPDATE THE EIGENVALUES AFTER EACH SWEEP
!
      II=1
      DO I=1,N
        IF (A(II).LE.0. .OR. B(II).LE.0.) THEN
          WRITE (IOUT,2020) II,A(II),B(II)
          STOP
        END IF
        EIGV(I)=A(II)/B(II)
        II=II + N1 - I
      END DO
      IF (IFPR.NE.0) THEN
      WRITE (IOUT,2030)
      WRITE (IOUT,2010) (EIGV(I),I=1,N)
      END IF
!
! CHECK FOR CONVERGENCE
!
      DO I=1,N
        TOL=RTOL*D(I)
        DIF=ABS(EIGV(I)-D(I))
        IF (DIF.GT.TOL) GO TO 280
      END DO
!
! CHECK ALL OFF-DIAGONAL ELEMENTS TO SEE IF ANOTHER SWEEP IS
! REQUIRED
!
      EPS=RTOL**2
      DO J=1,NR
        JM1=J - 1
        JP1=J + 1
        LJK=JM1*N - JM1*J/2
        JJ=LJK + J
        DO K=JP1,N
          KM1=K - 1
          JK=LJK + K
          KK=KM1*N - KM1*K/2 + K
          EPSA=(A(JK)/A(JJ))*(A(JK)/A(KK))
          EPSB=(B(JK)/B(JJ))*(B(JK)/B(KK))
          IF (EPSA.GE.EPS .OR. EPSB.GE.EPS) GO TO 280
        END DO
      END DO
!
! CONVERGENCE REACHED
!
      EXIT
!
! UPDATE D MATRIX AND START NEW SWEEP, IF ALLOWED
!
280   DO I=1,N
        D(I)=EIGV(I)
      END DO
    END DO
!
! SCALE EIGENVECTORS
!
    II=1
    DO I=1,N
      BB=SQRT(B(II))
      DO K=1,N
        X(K,I)=X(K,I)/BB
      END DO
      II=II + N1 - I
    END DO
!
    RETURN
!
2000 FORMAT (//,' SWEEP NUMBER IN *JACOBI* = ',I8)
2010 FORMAT (' ',6E20.12)
2020 FORMAT (' *** ERROR *** SOLUTION STOP',/, &
' MATRICES NOT POSITIVE DEFINITE',/, &
' II = ',I8,' A(II) = ',E20.12,' B(II) = ',E20.12)
2030 FORMAT (/,' CURRENT EIGENVALUES IN *JACOBI* ARE',/)
END