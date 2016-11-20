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

SUBROUTINE COLHT (ColumnHeight,ElementDOF,ElementLocationMatrix)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To calculate column heights                                     .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : NumberOfEquations
  IMPLICIT NONE
  INTEGER :: ElementDOF, ElementLocationMatrix(ElementDOF),ColumnHeight(NumberOfEquations)
  INTEGER :: I, SmallestIndex, CurrentIndex, IndexDisparity

  SmallestIndex=HUGE(1)   ! The largest integer number

  DO I=1,ElementDOF
     IF (ElementLocationMatrix(I) .NE. 0) THEN
        IF (ElementLocationMatrix(I)-SmallestIndex .LT. 0) SmallestIndex=ElementLocationMatrix(I)
     END IF
  END DO
  
  !SmallestIndex=minval(ElementLocationMatrix(1:ElementDOF))

  DO I=1,ElementDOF
     CurrentIndex = ElementLocationMatrix(I)
     IF (CurrentIndex.NE.0) THEN
        IndexDisparity = CurrentIndex - SmallestIndex
        IF (IndexDisparity.GT.ColumnHeight(CurrentIndex)) ColumnHeight(CurrentIndex) = IndexDisparity
     END IF
  END DO

  RETURN
END SUBROUTINE COLHT


SUBROUTINE ADDRES (MAXA,ColumnHeight)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To calculate addresses of diagonal elements in banded           .
! .   matrix whose column heights are known                           .
! .                                                                   .
! .   ColumnHeight  = Active column heights                           .
! .   MAXA = Addresses of diagonal elements                           .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : NumberOfEquations, MaxHalfBandwidth, NumberOfMatrixElements

  IMPLICIT NONE
  INTEGER :: MAXA(NumberOfEquations+1),ColumnHeight(NumberOfEquations)
  INTEGER :: I

! Clear array maxa

  DO I=1,NumberOfEquations + 1
     MAXA(I)=0.0
  END DO

  MAXA(1)=1
  MAXA(2)=2
  MaxHalfBandwidth=0
  IF (NumberOfEquations.GT.1) THEN
     DO I=2,NumberOfEquations
        IF (ColumnHeight(I).GT.MaxHalfBandwidth) MaxHalfBandwidth = ColumnHeight(I)
        MAXA(I+1)=MAXA(I) + ColumnHeight(I) + 1
     END DO
  END IF
  MaxHalfBandwidth = MaxHalfBandwidth + 1       !+1 to correctly output max half bandwidth
  NumberOfMatrixElements = MAXA(NumberOfEquations+1) - MAXA(1)

  RETURN
END SUBROUTINE ADDRES


SUBROUTINE ASSEM (ElementStiffnessMatrix)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To call element subroutines for assemblage of the               .
! .   structure stiffness matrix                                      .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : ElementTmpFile, NumberOfElementGroups, ElementGroupArraySize, NPAR

  IMPLICIT NONE
  REAL :: ElementStiffnessMatrix(*)
  INTEGER :: N, I

  REWIND ElementTmpFile
  DO N=1,NumberOfElementGroups
     READ (ElementTmpFile) ElementGroupArraySize,NPAR, &
                           (ElementStiffnessMatrix(I),I=1,ElementGroupArraySize)                      !Read Tmp File In.
     CALL ELEMNT                                    !IND=2 Entrance
  END DO

  RETURN
END SUBROUTINE ASSEM


SUBROUTINE ADDBAN (SkylineK,MAXA,S,ElementLocationMatrix,ElementDOF)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To assemble element stiffness into compacted global stiffness   .
! .                                                                   .
! .      SkylineK = GLOBAL STIFFNESS (1D skyline storage)                    .
! .      S = ELEMENT STIFFNESS                                        .
! .      ElementDOF = DEGREES OF FREEDOM IN ELEMENT STIFFNESS                 .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  USE GLOBALS, ONLY : NumberOfMatrixElements, NumberOfEquations
  IMPLICIT NONE
  REAL(8) :: SkylineK(NumberOfMatrixElements),S(ElementDOF,ElementDOF)
  INTEGER :: MAXA(NumberOfEquations+1),ElementLocationMatrix(ElementDOF)
  INTEGER :: ElementDOF, I, J, IndexI, IndexJ, SkylineIndex
  
  DO J=1,ElementDOF
     IndexJ = ElementLocationMatrix(J)
     IF (IndexJ .GT. 0) THEN
        DO I=1,J
           IndexI = ElementLocationMatrix(I)
           IF (IndexI .GT. 0) THEN
              IF (IndexJ .GE. IndexI) THEN
                 SkylineIndex= MAXA(IndexJ) + IndexJ - IndexI
              ELSE
                 SkylineIndex= MAXA(IndexI) + IndexI - IndexJ
              END IF              
              SkylineK(SkylineIndex)=SkylineK(SkylineIndex) + S(I,J)
           END IF
        END DO
     END IF
  END DO

  RETURN
END SUBROUTINE ADDBAN


SUBROUTINE COLSOL (SkylineK, LoadToDisplacement, MAXA, NumberOfEquations, NumberOfMatrixElements , NNM, SolutionMode)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To solve finite element static equilibrium equations in         .
! .   core, using compacted storage and column reduction scheme       .
! .                                                                   .
! .  - - Input variables - -                                          .
! .        SkylineK(NumberOfMatrixElements)    = Stiffness matrix stored in compacted form.
! .        LoadToDisplacement(NumberOfEquations)     = Right-hand-side load vector                    .
! .        MAXA(NNM) = Vector containing addresses of diagonal        .
! .                    elements of stiffness matrix in a              .
! .        NumberOfEquations        = Number of equations                            .
! .        NumberOfMatrixElements       = Number of elements below skyline of matrix     .
! .        NNM       = NumberOfEquations + 1                                         .
! .        SolutionMode    = Input flag                                     .
! .            EQ. 1   Triangularization of stiffness matrix          .
! .            EQ. 2   Reduction and back-substitution of load vector .
! .        OutputFile      = UNIT used for output                           .
! .                                                                   .
! .  - - OUTPUT - -                                                   .
! .        SkylineK(NumberOfMatrixElements)    = D and L - Factors of stiffness matrix.
! .        LoadToDisplacement(NumberOfEquations)     = Displacement vector                            .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : OutputFile

  IMPLICIT NONE
  INTEGER :: MAXA(NNM), NumberOfEquations, NumberOfMatrixElements, NNM, SolutionMode
  REAL(8) :: SkylineK(NumberOfMatrixElements), LoadToDisplacement(NumberOfEquations), C,B
  INTEGER :: N,K,KN,KL,KU,KH,IC,KLT,KI,J, ElementDOF,KK,L
  INTEGER :: MIN0

! Perform L*D*L(T) factorization of stiffness matrix

  IF (SolutionMode == 1) THEN

      DO N=1,NumberOfEquations
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
                ElementDOF=MAXA(K+1) - KI - 1
                IF (ElementDOF .GT. 0) THEN
                   KK=MIN0(IC,ElementDOF)
                   C=0.
                   DO L=1,KK
                      C=C + SkylineK(KI+L) * SkylineK(KLT+L)
                   END DO
                   SkylineK(KLT)=SkylineK(KLT) - C
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
                C=SkylineK(KK)/SkylineK(KI)
                B=B + C*SkylineK(KK)
                SkylineK(KK)=C
             END DO
             SkylineK(KN)=SkylineK(KN) - B
         ENDIF

         IF (SkylineK(KN) .LE. 0) THEN
            WRITE (OutputFile,"(//' STOP - STIFFNESS MATRIX NOT POSITIVE DEFINITE',//,  &
                            ' NONPOSITIVE PIVOT FOR EQUATION ',I8,//,' PIVOT = ',E20.12 )") N,SkylineK(KN)
            STOP
         END IF
      END DO
      
      !write (*,*) "Skyline",SkylineK
      !write (*,*) "MAXA", MAXA

  ELSE IF (SolutionMode == 2) THEN

! REDUCE RIGHT-HAND-SIDE LOAD VECTOR

       DO N=1,NumberOfEquations
         KL=MAXA(N) + 1
         KU=MAXA(N+1) - 1
         
         !write (*,*) "Load(",N,")",LoadToDisplacement(N)
         
         IF (KU-KL .GE. 0) THEN
            !K=N
            !C=0.
            !DO KK=KL,KU
            !   K=K - 1
            !   C=C + SkylineK(KK)*LoadToDisplacement(K)
            !END DO
            
			C = dot_product(SkylineK(KL:KU),LoadToDisplacement(N-1:N-(KU-KL)-1:-1))
			
			!write (*,*) "C= dot_product(SkylineK(KL:KU),LoadToDisplacement(N:N-(KU-KL)))",C
            
            LoadToDisplacement(N) = LoadToDisplacement(N) - C
         END IF
      END DO

! BACK-SUBSTITUTE

      DO N=1,NumberOfEquations
         K=MAXA(N)
         LoadToDisplacement(N) = LoadToDisplacement(N) / SkylineK(K)
      END DO

      IF (NumberOfEquations.EQ.1) RETURN

      N=NumberOfEquations
      DO L=2,NumberOfEquations
         KL=MAXA(N) + 1
         KU=MAXA(N+1) - 1
         IF (KU-KL .GE. 0) THEN
            K=N
            DO KK=KL,KU
               K=K - 1
               LoadToDisplacement(K) = LoadToDisplacement(K) - SkylineK(KK) * LoadToDisplacement(N)
            END DO
         END IF
         N = N - 1
      END DO

  END IF

END SUBROUTINE COLSOL
