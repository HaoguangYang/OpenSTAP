! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       .
! .                            S T A P 9 0                                .
! .                                                                       .
! .     AN IN-CORE SOLUTION STATIC ANALYSIS PROGRAM IN FORTRAN 90         .
! .     Adapted from STAP (KJ Bath, FORTRAN IV) for teaching purpose      .
! .                                                                       .
! .     Xiong Zhang, (2013)                                               .
! .     Computational Dynamics Group, School of Aerospace                 .
! .     Tsinghua Univerity                                                .
! .                                                                       .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .

SUBROUTINE COLHT (MHT,ND,LM)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To calculate column heights                                     .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : NEQ
  IMPLICIT NONE
  INTEGER :: ND, LM(ND),MHT(NEQ)
  INTEGER :: I, LS, II, ME

  LS=HUGE(1)   ! The largest integer number

  DO I=1,ND
     IF (LM(I) .NE. 0) THEN
        IF (LM(I)-LS .LT. 0) LS=LM(I)
     END IF
  END DO

  DO I=1,ND
     II=LM(I)
     IF (II.NE.0) THEN
        ME=II - LS
        IF (ME.GT.MHT(II)) MHT(II)=ME
     END IF
  END DO

  RETURN
END SUBROUTINE COLHT


SUBROUTINE ADDRES (MAXA,MHT)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To calculate addresses of diagonal elements in banded           .
! .   matrix whose column heights are known                           .
! .                                                                   .
! .   MHT  = Active column heights                                    .
! .   MAXA = Addresses of diagonal elements                           .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : NEQ, MK, NWK

  IMPLICIT NONE
  INTEGER :: MAXA(NEQ+1),MHT(NEQ)
  INTEGER :: NN, I

! Clear array maxa

  NN=NEQ + 1
  DO I=1,NN
     MAXA(I)=0.0
  END DO

  MAXA(1)=1
  MAXA(2)=2
  MK=0
  IF (NEQ.GT.1) THEN
     DO I=2,NEQ
        IF (MHT(I).GT.MK) MK=MHT(I)
        MAXA(I+1)=MAXA(I) + MHT(I) + 1
     END DO
  END IF
  MK=MK + 1
  NWK=MAXA(NEQ+1) - MAXA(1)

  RETURN
END SUBROUTINE ADDRES


SUBROUTINE ASSEM (AA)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To call element subroutines for assemblage of the               .
! .   structure stiffness matrix                                      .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : IELMNT, NUMEG, NUMEST, NPAR

  IMPLICIT NONE
  REAL :: AA(*)
  INTEGER :: N, I

  REWIND IELMNT
  DO N=1,NUMEG
     READ (IELMNT) NUMEST,NPAR,(AA(I),I=1,NUMEST)
     CALL ELEMNT
  END DO

  RETURN
END SUBROUTINE ASSEM


SUBROUTINE ADDBAN (A,MAXA,S,LM,ND)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To assemble element stiffness into compacted global stiffness   .
! .                                                                   .
! .      A = GLOBAL STIFFNESS (1D skyline storage)                    .
! .      S = ELEMENT STIFFNESS                                        .
! .      ND = DEGREES OF FREEDOM IN ELEMENT STIFFNESS                 .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  USE GLOBALS, ONLY : NWK, NEQ
  IMPLICIT NONE
  REAL(8) :: A(NWK),S(ND,ND)
  INTEGER :: MAXA(NEQ+1),LM(ND)
  INTEGER :: ND, I, J, II, JJ, KK
  
  DO J=1,ND
     JJ=LM(J)
     IF (JJ .GT. 0) THEN
        DO I=1,J
           II=LM(I)
           IF (II .GT. 0) THEN
              IF (JJ .GE. II) THEN
                 KK= MAXA(JJ) + JJ - II
              ELSE
                 KK= MAXA(II) + II - JJ
              END IF              
              A(KK)=A(KK) + S(I,J)
           END IF
        END DO
     END IF
  END DO

  RETURN
END SUBROUTINE ADDBAN


SUBROUTINE COLSOL (A,V,MAXA,NN,NWK,NNM,KKK)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To solve finite element static equilibrium equations in         .
! .   core, using compacted storage and column reduction scheme       .
! .                                                                   .
! .  - - Input variables - -                                          .
! .        A(NWK)    = Stiffness matrix stored in compacted form      .
! .        V(NN)     = Right-hand-side load vector                    .
! .        MAXA(NNM) = Vector containing addresses of diagonal        .
! .                    elements of stiffness matrix in a              .
! .        NN        = Number of equations                            .
! .        NWK       = Number of elements below skyline of matrix     .
! .        NNM       = NN + 1                                         .
! .        KKK       = Input flag                                     .
! .            EQ. 1   Triangularization of stiffness matrix          .
! .            EQ. 2   Reduction and back-substitution of load vector .
! .        IOUT      = UNIT used for output                           .
! .                                                                   .
! .  - - OUTPUT - -                                                   .
! .        A(NWK)    = D and L - Factors of stiffness matrix          .
! .        V(NN)     = Displacement vector                            .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : IOUT

  IMPLICIT NONE
  INTEGER :: MAXA(NNM),NN,NWK,NNM,KKK
  REAL(8) :: A(NWK),V(NN),C,B
  INTEGER :: N,K,KN,KL,KU,KH,IC,KLT,KI,J,ND,KK,L
  INTEGER :: MIN0

! Perform L*D*L(T) factorization of stiffness matrix

  IF (KKK == 1) THEN

      DO N=1,NN
         KN=MAXA(N)
         KL=KN + 1
         KU=MAXA(N+1) - 1
         KH=KU - KL

         IF (KH > 0) THEN
             K=N - KH
             IC=0
             KLT=KU
             DO J=1,KH
                IC=IC + 1
                KLT=KLT - 1
                KI=MAXA(K)
                ND=MAXA(K+1) - KI - 1
                IF (ND .GT. 0) THEN
                   KK=MIN0(IC,ND)
                   C=0.
                   DO L=1,KK
                      C=C + A(KI+L)*A(KLT+L)
                   END DO
                   A(KLT)=A(KLT) - C
                END IF
                K=K + 1
             END DO
         ENDIF

         IF (KH >= 0) THEN
             K=N
             B=0.
             DO KK=KL,KU
                K=K - 1
                KI=MAXA(K)
                C=A(KK)/A(KI)
                B=B + C*A(KK)
                A(KK)=C
             END DO
             A(KN)=A(KN) - B
         ENDIF

         IF (A(KN) .LE. 0) THEN
            WRITE (IOUT,"(//' STOP - STIFFNESS MATRIX NOT POSITIVE DEFINITE',//,  &
                            ' NONPOSITIVE PIVOT FOR EQUATION ',I8,//,' PIVOT = ',E20.12 )") N,A(KN)
            STOP
         END IF
      END DO

  ELSE IF (KKK == 2) THEN

! REDUCE RIGHT-HAND-SIDE LOAD VECTOR

       DO N=1,NN
         KL=MAXA(N) + 1
         KU=MAXA(N+1) - 1
         IF (KU-KL .GE. 0) THEN
            K=N
            C=0.
            DO KK=KL,KU
               K=K - 1
               C=C + A(KK)*V(K)
            END DO
            V(N)=V(N) - C
         END IF
      END DO

! BACK-SUBSTITUTE

      DO N=1,NN
         K=MAXA(N)
         V(N)=V(N)/A(K)
      END DO

      IF (NN.EQ.1) RETURN

      N=NN
      DO L=2,NN
         KL=MAXA(N) + 1
         KU=MAXA(N+1) - 1
         IF (KU-KL .GE. 0) THEN
            K=N
            DO KK=KL,KU
               K=K - 1
               V(K)=V(K) - A(KK)*V(N)
            END DO
         END IF
         N=N - 1
      END DO

  END IF

END SUBROUTINE COLSOL
