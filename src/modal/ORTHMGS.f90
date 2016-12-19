    SUBROUTINE ORTHMGS(Q,NQ,R,NR,NC,NN,B,NWM,MAXA,X,TT)
!
!-----[--.----+----.----+----.-----------------------------------------]
!      PURPOSE: USE GRAM-SCHMIDT ORTHOGONALIzATION Q

!      INPUTS:
!       Q(NN,NQ)    - ORIGIN ORTHOGONALIzATION MATRIx
!        NQ            - NUMBER OF COLUMNS TO ORTHOGONALIzE IN Q
!       R(NN,NR)    - ORIGIN ORTHOGONALIzATION MATRIx
!       NR          - NUMBER OF COLUMNS OF MATRIX R
!        NC            - NUMBER OF COLUMNS TO ORTHOGONALIzE IN Q
!                   - NC.GT.0, IF HAVE THE EIGENVALUE
!       NN            - NUMBER OF ROwS IN Q, R
!       B(*)        - MASS MATRIx
!        NWM         - NUMBER FOR MASS STORE

!      SCRATCH:
!        TT(NN)  - WORKING vECTOR

!      OUTPUTS:
!         X  - ORTHOGONALIzATION MATRIx
!-----[--.----+----.----+----.-----------------------------------------]
!
    IMPLICIT NONE
    INTEGER, INTENT(IN):: NN,NQ,NR,NC,NWM,MAXA(*)
    REAL(8), INTENT(IN):: Q(NN,*),R(NN,*),B(*)
    REAL*8   ALPHA,TT(NN),X(NN),Y(NN)
    INTEGER I,J,K

    Y = X

    DO K=1,NQ
      CALL MULT(TT,B(1),Q(1,K),MAXA,NN,NWM)
      ALPHA=DOT_PRODUCT(X,TT)
      DO I=1,NN
        Y(I)=Y(I)-ALPHA*Q(I,K)
      END DO
    END DO
    DO K=1,NC
      CALL MULT(TT,B(1),R(1,K),MAXA,NN,NWM)
      ALPHA=DOT_PRODUCT(X,TT)
      DO I=1,NN
        Y(I)=Y(I)-ALPHA*R(I,K)
      END DO
    END DO

    X = Y

    END