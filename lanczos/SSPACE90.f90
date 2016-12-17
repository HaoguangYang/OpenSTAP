      SUBROUTINE SSPACE90(A,B,MAXA,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO, &
      BUPC,NN,NNM,NWK,NWM,NROOT,RTOL,NC,NNC,NITEM,IFSS,IFPR,NSTIF,IOUT)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .
! .   P R O G R A M
! .        TO SOLVE FOR THE SMALLEST EIGENVALUES-- ASSUMED .GT. 0 --
! .        AND CORRESPONDING EIGENVECTORS IN THE GENERALIZED
! .        EIGENPROBLEM USING THE SUBSPACE ITERATION METHOD
! .
! .  - - INPUT VARIABLES - -
! .        A(NWK)    = STIFFNESS MATRIX IN COMPACTED FORM (ASSUMED
! .                    POSITIVE DEFINITE)
! .        B(NWM)    = MASS MATRIX IN COMPACTED FORM
! .        MAXA(NNM) = VECTOR CONTAINING ADDRESSES OF DIAGONAL
! .                    ELEMENTS OF STIFFNESS MATRIX A
! .        R(NN,NC)  = STORAGE FOR EIGENVECTORS
! .        EIGV(NC)  = STORAGE FOR EIGENVALUES
! .        TT(NN)    = WORKING VECTOR
! .        W(NN)     = WORKING VECTOR
! .        AR(NNC)   = WORKING MATRIX STORING PROJECTION OF K
! .        BR(NNC)   = WORKING MATRIX STORING PROJECTION OF M
! .        VEC(NC,NC)= WORKING MATRIX
! .        D(NC)     = WORKING VECTOR
! .        RTOLV(NC) = WORKING VECTOR
! .        BUP(NC)   = WORKING VECTOR
! .        BLO(NC)   = WORKING VECTOR
! .        BUPC(NC)  = WORKING VECTOR
! .        NN        = ORDER OF STIFFNESS AND MASS MATRICES
! .        NNM       = NN + 1
! .        NWK       = NUMBER OF ELEMENTS BELOW SKYLINE OF
! .                    STIFFNESS MATRIX
! .        NWM       = NUMBER OF ELEMENTS BELOW SKYLINE OF
! .                    MASS MATRIX
! .                      I. E. NWM=NWK FOR CONSISTENT MASS MATRIX
! .                            NWM=NN  FOR LUMPED MASS MATRIX
! .        NROOT     = NUMBER OF REQUIRED EIGENVALUES AND EIGENVECTORS.
! .        RTOL      = CONVERGENCE TOLERANCE ON EIGENVALUES
! .                    ( 1.E-06 OR SMALLER )
! .        NC        = NUMBER OF ITERATION VECTORS USED
! .                    (USUALLY SET TO MIN(2*NROOT, NROOT+8), BUT NC
! .                    CANNOT BE LARGER THAN THE NUMBER OF MASS
! .                    DEGREES OF FREEDOM)
! .        NNC       = NC*(NC+1)/2 DIMENSION OF STORAGE VECTORS AR,BR
! .        NITEM     = MAXIMUM NUMBER OF SUBSPACE ITERATIONS PERMITTED.
! .                    (USUALLY SET TO 16)
! .                    THE PARAMETERS NC AND/OR NITEM MUST BE
! .                    INCREASED IF A SOLUTION HAS NOT CONVERGED
! .        IFSS      = FLAG FOR STURM SEQUENCE CHECK
! .                      EQ.0  NO CHECK
! .                      EQ.1  CHECK
! .        IFPR      = FLAG FOR PRINTING DURING ITERATION
! .                      EQ.0  NO PRINTING
! .                      EQ.1  PRINT
! .        NSTIF     = SCRATCH FILE
! .        IOUT      = UNIT USED FOR OUTPUT
! .
! .  - - OUTPUT - -
! .        EIGV(NROOT) = EIGENVALUES
! .        R(NN,NROOT) = EIGENVECTORS
! .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!     
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NN,NNM,NWK,NWM,NROOT,NC,NNC,NITEM
      INTEGER, INTENT(IN) :: MAXA(NNM),IFPR,NSTIF,IOUT
      INTEGER, INTENT(INOUT) :: IFSS

      REAL(8), INTENT(IN) :: B(NWM)
      REAL(8), INTENT(IN) :: BUP(NC)
      REAL(8), INTENT(IN) :: BLO(NC),BUPC(NC),RTOL

      REAL(8), INTENT(INOUT) :: A(NWK),D(NC),W(NN),TT(NN),AR(NNC),BR(NNC)
      REAL(8), INTENT(INOUT) :: EIGV(NC),R(NN,NC),RTOLV(NC),VEC(NC,NC)

      REAL(8) :: ART,BRT,TOLJ,TOLJ2,RT,BT,PI,XX,SHIFT
      REAL(8) :: VDOT,VNORM,WNORM,EIGV2,DIF,RDIF,EIGVT

      INTEGER :: ICONV,NSCH,NSMAX,N1,NC1,I,J,K,IJ,L,II,ND,IX
      INTEGER :: ITEMP,IS,ISH,NMIS,NITE,NEI,LRT
!     
!     SET TOLERANCE FOR JACOBI ITERATION
      TOLJ=1.0D-12
!     
!     INITIALIZATION
!     
      ICONV=0
      NSCH=0
      NSMAX=12
      N1=NC + 1
      NC1=NC - 1
!     
      REWIND NSTIF
      WRITE (NSTIF) A
!     
      DO I=1,NC
         D(I)=0.
      END DO
!     
!     ESTABLISH STARTING ITERATION VECTORS
!     
      ND=NN/NC
      IF (NWM .LE. NN) THEN
         J=0
         DO I=1,NN
            II=MAXA(I)
            R(I,1)=B(I)
            IF (B(I).GT.0) J=J + 1
            W(I)=B(I)/A(II)
         END DO
         IF (NC.GT.J) THEN
            WRITE (IOUT,1007)
            STOP
         END IF
      ELSE
         DO I=1,NN
            II=MAXA(I)
            R(I,1)=B(II)
            W(I)=B(II)/A(II)
         END DO
      END IF
!     
      DO J=2,NC
         DO I=1,NN
            R(I,J)=0.
         END DO
      END DO
!     
      L=NN - ND
      DO J=2,NC
         RT=0.
         DO I=1,L
            IF (W(I).LT.RT) CYCLE
            RT=W(I)
            IJ=I
         END DO
!     
         DO I=L,NN
            IF (W(I).LE.RT) CYCLE
            RT=W(I)
            IJ=I
         END DO
!     
         TT(J)=FLOAT(IJ)
         W(IJ)=0.
         L=L - ND
         R(IJ,J)=1.
      END DO
!     
      WRITE (IOUT,1008)
      WRITE (IOUT,1002) (TT(J),J=2,NC)
!     
!     A RANDOM VECTOR IS ADDED TO THE LAST VECTOR
!     
      PI=3.141592654D0
      XX=0.5D0
      DO K=1,NN
         XX=(PI + XX)**5
         IX=INT(XX)
         XX=XX - FLOAT(IX)
         R(K,NC)=R(K,NC) + XX
      END DO
!     
!     FACTORIZE MATRIX A INTO (L)*(D)*(L(T))
!     
      ISH=0
      CALL DECOMP (A,MAXA,NN,ISH,IOUT)
!     
!     - - - S T A R T   O F   I T E R A T I O N   L O O P
!     
      NITE=0
      TOLJ2=1.0D-24
!     
      DO
         NITE=NITE + 1
         IF (IFPR.NE.0) WRITE (IOUT,1010) NITE
!     
!     CALCULATE THE PROJECTIONS OF A AND B
!     
         IJ=0
         DO  J=1,NC
            DO  K=1,NN
               TT(K)=R(K,J) 
            END DO
            CALL REDBAK (A,TT,MAXA,NN)
            DO I=J,NC
               ART=0.
               DO K=1,NN
                  ART=ART + R(K,I)*TT(K)
               END DO
               IJ=IJ + 1
               AR(IJ)=ART
            END DO
            DO K=1,NN
               R(K,J)=TT(K)
            END DO
         END DO
!     
         IJ=0
         DO J=1,NC
            CALL MULT (TT,B,R(1,J),MAXA,NN,NWM)
            DO I=J,NC
               BRT=0.
               DO K=1,NN
                  BRT=BRT + R(K,I)*TT(K)
               END DO
               IJ=IJ + 1
               BR(IJ)=BRT
            END DO
            IF (ICONV.GT.0) CYCLE
            DO K=1,NN
               R(K,J)=TT(K)
            END DO
         END DO
!     
!     SOLVE FOR EIGENSYSTEM OF SUBSPACE OPERATORS
!     
         IF (IFPR.NE.0) CALL DISPARBR(NC,NNC,N1,AR,BR,IOUT)
!     
         CALL JACOBI (AR,BR,VEC,EIGV,W,NC,NNC,TOLJ,NSMAX,IFPR,IOUT)
!     
         IF (IFPR.NE.0) CALL DISPARBR(NC,NNC,N1,AR,BR,IOUT)
!     
!     ARRANGE EIGENVALUES IN ASCENDING ORDER
!     
         DO
            IS=0
            II=1
            DO I=1,NC1
               ITEMP=II + N1 - I
               IF (EIGV(I+1).GE.EIGV(I)) CYCLE
               IS=IS + 1
               EIGVT=EIGV(I+1)
               EIGV(I+1)=EIGV(I)
               EIGV(I)=EIGVT
               BT=BR(ITEMP)
               BR(ITEMP)=BR(II)
               BR(II)=BT
               DO K=1,NC
                  RT=VEC(K,I+1)
                  VEC(K,I+1)=VEC(K,I)
                  VEC(K,I)=RT
               END DO
               II=ITEMP
            END DO
!     
            IF (IS.LE.0) EXIT
         END DO
!     
         IF (IFPR.NE.0) THEN
            WRITE (IOUT,1035)
            WRITE (IOUT,1006) (EIGV(I),I=1,NC)
         END IF
!     
!     CALCULATE B TIMES APPROXIMATE EIGENVECTORS (ICONV.EQ.0)
!     OR     FINAL EIGENVECTOR APPROXIMATIONS (ICONV.GT.0)
!     
         DO I=1,NN
            DO J=1,NC
               TT(J)=R(I,J)
            END DO
            DO K=1,NC
               RT=0.
               DO L=1,NC
                  RT=RT + TT(L)*VEC(L,K)
               END DO
               R(I,K)=RT
            END DO
         END DO
!     
!     CALCULATE ERROR BOUNDS AND CHECK FOR CONVERGENCE OF EIGENVALUES
!     
         DO I=1,NC
            VDOT=0.
            DO J=1,NC
               VDOT=VDOT + VEC(I,J)*VEC(I,J)
            END DO
            EIGV2=EIGV(I)*EIGV(I)
            DIF=VDOT - EIGV2
            RDIF=MAX(DIF,TOLJ2*EIGV2)/EIGV2
            RDIF=SQRT(RDIF)
            RTOLV(I)=RDIF
         END DO
!     
         IF (IFPR.NE.0 .OR. ICONV.NE.0) THEN
            WRITE (IOUT,1050)
            WRITE (IOUT,1005) (RTOLV(I),I=1,NC)
         END IF
!     
         IF (ICONV.GT.0) EXIT
!     
         LRT=0     
         DO I=1,NROOT
            IF (RTOLV(I).GT.RTOL) THEN
               LRT=1
               EXIT
            END IF
         END DO
!         
         IF (LRT.EQ.0) THEN
            WRITE (IOUT,1060) RTOL
            ICONV=1
         ELSE
            IF (NITE.GE.NITEM) THEN
               WRITE (IOUT,1070)
               ICONV=2
               IFSS=0
            END IF
         END IF
      END DO
!     
!     - - - E N D   O F   I T E R A T I O N   L O O P
!     
      WRITE (IOUT,1100)
      WRITE (IOUT,1006) (EIGV(I),I=1,NROOT)
      WRITE (IOUT,1110)
      DO J=1,NROOT
         WRITE (IOUT,1005) (R(K,J),K=1,NN)
      END DO
!     
!     CALCULATE AND PRINT ERROR MEASURES
!     
      REWIND NSTIF
      READ (NSTIF) A
!     
      DO L=1,NROOT
         RT=EIGV(L)
         CALL MULT(TT,A,R(1,L),MAXA,NN,NWK)
         VNORM=0.
         DO I=1,NN
            VNORM=VNORM + TT(I)*TT(I)
         END DO
         CALL MULT(W,B,R(1,L),MAXA,NN,NWM)
         WNORM=0.
         DO I=1,NN
            TT(I)=TT(I) - RT*W(I)
            WNORM=WNORM + TT(I)*TT(I)
         END DO
         VNORM=SQRT(VNORM)
         WNORM=SQRT(WNORM)
         D(L)=WNORM/VNORM
      END DO
!
      WRITE (IOUT,1115)
      WRITE (IOUT,1005) (D(I),I=1,NROOT)
!     
!     APPLY STURM SEQUENCE CHECK
!     
      IF (IFSS.EQ.0) RETURN

      CALL SCHECK(EIGV,RTOLV,BUP,BLO,BUPC,D,NC,NEI,RTOL,SHIFT,IOUT)
!     
      WRITE (IOUT,1120) SHIFT
!     
!     SHIFT MATRIX A
!     
      REWIND NSTIF
      READ (NSTIF) A
!
      IF (NWM.LE.NN) THEN
         DO I=1,NN
            II=MAXA(I)
            A(II)=A(II) - B(I)*SHIFT
         END DO
      ELSE
         DO I=1,NWK
            A(I)=A(I) - B(I)*SHIFT
         END DO
      END IF
!     
!     FACTORIZE SHIFTED MATRIX
!     
      ISH=1
      CALL DECOMP (A,MAXA,NN,ISH,IOUT)
!     
!     COUNT NUMBER OF NEGATIVE DIAGONAL ELEMENTS
!     
      NSCH=0
      DO I=1,NN
         II=MAXA(I)
         IF (A(II).LT.0.) NSCH=NSCH + 1
      END DO
!
      IF (NSCH.NE.NEI) THEN
         NMIS=NSCH - NEI
         WRITE (IOUT,1130) NMIS
      ELSE
         WRITE (IOUT,1140) NSCH
      END IF
!     
      RETURN
!     
 1002 FORMAT (' ',10F10.0)
 1005 FORMAT (' ',12E11.4)
 1006 FORMAT (' ',6E22.14)
 1007 FORMAT (///,' STOP, NC IS LARGER THAN THE NUMBER OF MASS ', &
                  'DEGREES OF FREEDOM')
 1008 FORMAT (///,' DEGREES OF FREEDOM EXCITED BY UNIT STARTING ', &
                  'ITERATION VECTORS')
 1010 FORMAT (//,' I T E R A T I O N   N U M B E R ',I8)
 1035 FORMAT (/,' EIGENVALUES OF AR-LAMBDA*BR')
 1040 FORMAT (//,' AR AND BR AFTER JACOBI DIAGONALIZATION')
 1050 FORMAT (/,' ERROR BOUNDS REACHED ON EIGENVALUES')
 1060 FORMAT (///,' CONVERGENCE REACHED FOR RTOL ',E10.4)
 1070 FORMAT (' *** NO CONVERGENCE IN MAXIMUM NUMBER OF ITERATIONS', &
              ' PERMITTED',/, &
              ' WE ACCEPT CURRENT ITERATION VALUES',/, &
              ' THE STURM SEQUENCE CHECK IS NOT PERFORMED')
 1100 FORMAT (///,' THE CALCULATED EIGENVALUES ARE')
 1115 FORMAT (//,' ERROR MEASURES ON THE EIGENVALUES')
 1110 FORMAT (//,' THE CALCULATED EIGENVECTORS ARE',/)
 1120 FORMAT (///,' CHECK APPLIED AT SHIFT ',E22.14)
 1130 FORMAT (//,' THERE ARE ',I8,' EIGENVALUES MISSING')
 1140 FORMAT (//,' WE FOUND THE LOWEST ',I8,' EIGENVALUES')
!     
      END

      SUBROUTINE DISPARBR(NC,NNC,N1,AR,BR,IOUT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NC,NNC,N1,IOUT
      REAL(8), INTENT(INOUT) :: AR(NNC),BR(NNC)
! 
      INTEGER I,J,II,ITEMP
!
      WRITE (IOUT,1020)
!
      II=1
      DO I=1,NC
         ITEMP=II + NC - I
         WRITE (IOUT,1005) (AR(J),J=II,ITEMP)
         II=II + N1 - I
      END DO
!
      WRITE (IOUT,1030)

      II=1
      DO I=1,NC
         ITEMP=II + NC - I
         WRITE (IOUT,1005) (BR(J),J=II,ITEMP)
         II=II + N1 - I
      END DO
!
      RETURN
 1020 FORMAT (/,' PROJECTION OF A (MATRIX AR)')
 1030 FORMAT (/,' PROJECTION OF B (MATRIX BR)')
 1005 FORMAT (' ',12E11.4)
      END
!
      SUBROUTINE DECOMP(A,MAXA,NN,ISH,IOUT)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .
! .   P R O G R A M
! .        TO CALCULATE (L)*(D)*(L)(T) FACTORIZATION OF
! .        STIFFNESS MATRIX
! .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!     
      IMPLICIT NONE
      REAL(8), INTENT(INOUT) :: A(1)
      INTEGER, INTENT(IN) :: MAXA(1),NN,ISH,IOUT
! 
      INTEGER :: J,N,K,L,IC,KI,KN,KL,KU,KH,KK,KLT,ND
      REAL(8) :: B,C
!
      IF (NN.EQ.1) RETURN
!     
      DO N=1,NN
         KN=MAXA(N)
         KL=KN + 1
         KU=MAXA(N+1) - 1
         KH=KU - KL
         IF (KH.GT.0) THEN
            K=N - KH
            IC=0
            KLT=KU
            DO J=1,KH
               IC=IC + 1
               KLT=KLT - 1
               KI=MAXA(K)
               ND=MAXA(K+1) - KI - 1
               IF (ND.GT.0) THEN
                  KK=MIN0(IC,ND)
                  C=0.
                  DO L=1,KK
                     C=C + A(KI+L)*A(KLT+L)
                  END DO
                  A(KLT)=A(KLT) - C
               END IF
               K=K + 1
            END DO
         END IF
!
         IF (KH.GE.0) THEN
            K=N
            B=0.
            DO KK=KL,KU
               K=K - 1
               KI=MAXA(K)
               C=A(KK)/A(KI)
               IF (ABS(C).GE.1.E07) THEN
                  WRITE (IOUT,2010) N,C
                  STOP
               END IF
               B=B + C*A(KK)
               A(KK)=C
            END DO
            A(KN)=A(KN) - B
         END IF
!
         IF (A(KN).LE.0) THEN
            IF (ISH.EQ.0) THEN
               WRITE (IOUT,2000) N,A(KN)
               STOP
            END IF
            IF (A(KN).EQ.0.) A(KN)=-1.E-16
         END IF
      END DO
!     
      RETURN
!     
 2000 FORMAT (//' STOP - STIFFNESS MATRIX NOT POSITIVE DEFINITE',//, &
                ' NONPOSITIVE PIVOT FOR EQUATION ',I8,//, &
                ' PIVOT = ',E20.12)
 2010 FORMAT (//' STOP - STURM SEQUENCE CHECK FAILED BECAUSE OF', &
                ' MULTIPLIER GROWTH FOR COLUMN NUMBER ',I8,//, &
                ' MULTIPLIER = ',E20.8)
      END
!
      SUBROUTINE REDBAK (A,V,MAXA,NN)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .
! .   P R O G R A M
! .        TO REDUCE AND BACK-SUBSTITUTE ITERATION VECTORS
! .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!     
      IMPLICIT NONE
!
      REAL(8), INTENT(INOUT) :: A(1),V(1)
      INTEGER, INTENT(IN) ::NN,MAXA(1)
!
      INTEGER :: L,K,N,KL,KU,KK
      REAL(8) :: C
!     
      DO N=1,NN
         KL=MAXA(N) + 1
         KU=MAXA(N+1) - 1
         IF (KU-KL.LT.0) CYCLE
         K=N
         C=0.
         DO KK=KL,KU
            K=K - 1
            C=C + A(KK)*V(K)
         END DO
         V(N)=V(N) - C
      END DO
!     
      DO N=1,NN
         K=MAXA(N)
         V(N)=V(N)/A(K)
      END DO
!     
      IF (NN.EQ.1) RETURN
!     
      N=NN
      DO L=2,NN
         KL=MAXA(N) + 1
         KU=MAXA(N+1) - 1
         IF (KU-KL.LT.0) CYCLE
         K=N
         DO KK=KL,KU
            K=K - 1
            V(K)=V(K) - A(KK)*V(N)
         END DO
         N=N - 1
      END DO
!     
      RETURN
      END

      SUBROUTINE MULT (TT,B,RR,MAXA,NN,NWM)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .
! .   P R O G R A M
! .        TO EVALUATE PRODUCT OF B TIMES RR AND STORE RESULT IN TT
! .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!     
      IMPLICIT NONE
      REAL(8), INTENT(INOUT) :: TT(1)
      REAL(8), INTENT(IN) :: B(1),RR(1)
      INTEGER, INTENT(IN) :: NN,NWM,MAXA(1)
!
      INTEGER :: I,II,KL,KU,KK
      REAL(8) :: AA,CC
!     
      IF (NWM.LE.NN) THEN
         DO I=1,NN
            TT(I)=B(I)*RR(I)
         END DO
         RETURN
      END IF
!     
      DO I=1,NN
         TT(I)=0.
      END DO
!     
      DO I=1,NN
         KL=MAXA(I)
         KU=MAXA(I+1) - 1
         II=I + 1
         CC=RR(I)
         DO KK=KL,KU
            II=II - 1
            TT(II)=TT(II) + B(KK)*CC
         END DO
      END DO
!
      IF (NN.EQ.1) RETURN
!
      DO I=2,NN
         KL=MAXA(I) + 1
         KU=MAXA(I+1) - 1
         IF (KU-KL.LT.0) CYCLE
         II=I
         AA=0.
         DO KK=KL,KU
            II=II - 1
            AA=AA + B(KK)*RR(II)
         END DO
         TT(I)=TT(I) + AA
      END DO
!     
      RETURN
      END

      SUBROUTINE SCHECK(EIGV,RTOLV,BUP,BLO,BUPC,NEIV,NC,NEI,RTOL,SHIFT,IOUT)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .
! .   P R O G R A M
! .        TO EVALUATE SHIFT FOR STURM SEQUENCE CHECK
! .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NC,IOUT
      REAL(8), INTENT(IN) :: EIGV(NC),RTOLV(NC),RTOL
      REAL(8), INTENT(INOUT) :: BUP(NC),BUPC(NC),BLO(NC),SHIFT
      INTEGER, INTENT(INOUT) :: NEI,NEIV(NC)
!
      INTEGER :: I,L,LM,LL,NROOT
      REAL(8) :: FTOL
!     
      FTOL=0.01
!     
      DO I=1,NC
         BUP(I)=EIGV(I)*(1.+FTOL)
         BLO(I)=EIGV(I)*(1.-FTOL)
      END DO
!
      NROOT=0
      DO I=1,NC
         IF (RTOLV(I).LT.RTOL) NROOT=NROOT + 1
      END DO
!
      IF (NROOT.LT.1) THEN
         WRITE (IOUT,1010)
         STOP
      END IF
!     
!     FIND UPPER BOUNDS ON EIGENVALUE CLUSTERS
!     
      DO I=1,NROOT
         NEIV(I)=1
      END DO
!
      IF (NROOT.EQ.1) THEN 
         BUPC(1)=BUP(1)
         LM=1
         L=1
         I=2
      ELSE 
         L=1
         I=2
         DO
            DO
               IF (BUP(I-1).LE.BLO(I)) EXIT
               NEIV(L)=NEIV(L) + 1
               I=I + 1
               IF (I.GT.NROOT) EXIT
            END DO
            BUPC(L)=BUP(I-1)
            IF (I.GT.NROOT) GO TO 290
            L=L + 1
            I=I + 1
            IF (I.GT.NROOT) EXIT
         END DO
         BUPC(L)=BUP(I-1)
 290     LM=L
         IF (NROOT.EQ.NC) GO TO 300
      END IF
!
      DO
         IF (BUP(I-1).LE.BLO(I)) EXIT
         IF (RTOLV(I).GT.RTOL) EXIT
         BUPC(L)=BUP(I)
         NEIV(L)=NEIV(L) + 1
         NROOT=NROOT + 1
         IF (NROOT.EQ.NC) EXIT
         I=I + 1
      END DO
!     
!     FIND SHIFT
!     
 300  WRITE (IOUT,1020)
      WRITE (IOUT,1005) (BUPC(I),I=1,LM)
      WRITE (IOUT,1030)
      WRITE (IOUT,1006) (NEIV(I),I=1,LM)
!
      LL=LM - 1
      IF (LM.NE.1) THEN
         DO
            DO I=1,LL
               NEIV(L)=NEIV(L) + NEIV(I)
            END DO
            L=L - 1
            LL=LL - 1
            IF (L.EQ.1) EXIT
         END DO
      END IF
!
      WRITE (IOUT,1040)
      WRITE (IOUT,1006) (NEIV(I),I=1,LM)
!
      L=0
      DO  I=1,LM
         L=L + 1
         IF (NEIV(I).GE.NROOT) EXIT
      END DO
      SHIFT=BUPC(L)
      NEI=NEIV(L)
!     
      RETURN
!     
 1005 FORMAT (' ',6E22.14)
 1006 FORMAT (' ',6I22)
 1010 FORMAT (' *** ERROR ***  SOLUTION STOP IN *SCHECK*',/, &
              ' NO EIGENVALUES FOUND',/)
 1020 FORMAT (///,' UPPER BOUNDS ON EIGENVALUE CLUSTERS')
 1030 FORMAT (//,' NO. OF EIGENVALUES IN EACH CLUSTER')
 1040 FORMAT (' NO. OF EIGENVALUES LESS THAN UPPER BOUNDS')
      END
!
      SUBROUTINE JACOBI (A,B,X,EIGV,D,N,NWA,RTOL,NSMAX,IFPR,IOUT)
! ....................................................................
! .
! .   P R O G R A M
! .        TO SOLVE THE GENERALIZED EIGENPROBLEM USING THE
! .        GENERALIZED JACOBI ITERATION
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
!     INITIALIZE EIGENVALUE AND EIGENVECTOR MATRICES
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
!     INITIALIZE SWEEP COUNTER AND BEGIN ITERATION
!     
      NSWEEP=0
      NR=N - 1
!
      DO WHILE (NSWEEP.LT.NSMAX)
         NSWEEP=NSWEEP + 1
         IF (IFPR.EQ.1) WRITE (IOUT,2000) NSWEEP
!     
!     CHECK IF PRESENT OFF-DIAGONAL ELEMENT IS LARGE ENOUGH TO REQUIRE
!     ZEROING
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
!     IF ZEROING IS REQUIRED, CALCULATE THE ROTATION MATRIX ELEMENTS CA
!     AND CG
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
!     PERFORM THE GENERALIZED ROTATION TO ZERO THE PRESENT OFF-DIAGONAL
!     ELEMENT
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
!     UPDATE THE EIGENVECTOR MATRIX AFTER EACH ROTATION
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
!     UPDATE THE EIGENVALUES AFTER EACH SWEEP
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
!     CHECK FOR CONVERGENCE
!     
         DO I=1,N
            TOL=RTOL*D(I)
            DIF=ABS(EIGV(I)-D(I))
            IF (DIF.GT.TOL) GO TO 280
         END DO
!     
!     CHECK ALL OFF-DIAGONAL ELEMENTS TO SEE IF ANOTHER SWEEP IS
!     REQUIRED
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
!     CONVERGENCE REACHED
!
         EXIT
!     
!     UPDATE  D  MATRIX AND START NEW SWEEP, IF ALLOWED
!     
 280     DO I=1,N
            D(I)=EIGV(I)
         END DO
      END DO
!     
!     SCALE EIGENVECTORS
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
