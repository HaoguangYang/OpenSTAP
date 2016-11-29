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

SUBROUTINE INPUT (NodeIdentifier,X,Y,Z,NumberOfNodalPoints,NumberOfEquations)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       .
! .   To read, generate, and print nodal point input data                 .
! .   To calculate equation numbers and store them in NodeIdentifier arrray           .
! .                                                                       .
! .      N = Element number                                               .
! .      NodeIdentifier = Boundary condition codes (0=free,1=deleted)                 .
! .      X,Y,Z = Coordinates                                              .
! .      KN = Generation code                                             .
! .           i.e. increment on nodal point number                        .
! .                                                                       .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : InputFile, OutputFile, NPAR

  IMPLICIT NONE
  INTEGER :: NumberOfNodalPoints, NumberOfEquations, NodeIdentifier(3,NumberOfNodalPoints)
  REAL(8) :: X(NumberOfNodalPoints), Y(NumberOfNodalPoints), Z(NumberOfNodalPoints)
  INTEGER :: I, J, N

! Read nodal point data

  N = 0
  NodeIdentifier(:,:) = 1
  X(:) = 0
  Y(:) = 0
  Z(:) = 0
  
  DO WHILE (N.NE.NumberOfNodalPoints)
      READ (InputFile,"(4I5,3F10.0,I5)") N,(NodeIdentifier(I,N),I=1,3),X(N),Y(N),Z(N)
  END DO

! Write complete nodal data

  WRITE (OutputFile,"(//,' N O D A L   P O I N T   D A T A',/)")

  WRITE (OutputFile,"('  NODE',10X,'BOUNDARY',25X,'NODAL POINT',/,  &
                ' NUMBER     CONDITION  CODES',21X,'COORDINATES', /,15X, &
                'X    Y    Z',15X,'X',12X,'Y',12X,'Z')")

  DO N=1,NumberOfNodalPoints
     WRITE (OutputFile,"(I5,6X,3I5,6X,3F13.3)") N,(NodeIdentifier(I,N),I=1,3),X(N),Y(N),Z(N)
  END DO

! Number unknowns

  NumberOfEquations=0
  DO N=1,NumberOfNodalPoints
     DO I=1,3
        IF (NodeIdentifier(I,N) .EQ. 0) THEN
           NumberOfEquations=NumberOfEquations + 1
           NodeIdentifier(I,N)=NumberOfEquations
        ELSE
           NodeIdentifier(I,N)=0
        END IF
     END DO
  END DO

! Write equation numbers
  WRITE (OutputFile,"(//,' EQUATION NUMBERS',//,'   NODE',9X,  &
                   'DEGREES OF FREEDOM',/,'  NUMBER',/,  &
                   '     N',13X,'X    Y    Z',/,(1X,I5,9X,3I5))") (N,(NodeIdentifier(I,N),I=1,3),N=1,NumberOfNodalPoints)

  RETURN

END SUBROUTINE INPUT


SUBROUTINE LOADS (R, NodeOfLoad, LoadDirection, FLOAD, NodeIdentifier, L, NumberOfConcentratedLoads, &
                  NumberOfEquations)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To read nodal load data                                         .
! .   To calculate the load vector r for each load case and           .
! .   write onto unit LoadTmpFile                                           .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  USE GLOBALS, ONLY : InputFile, OutputFile, LoadTmpFile, SolutionMode

  IMPLICIT NONE
  INTEGER :: NumberOfConcentratedLoads, NumberOfEquations, NodeIdentifier(3,*)
  INTEGER :: NodeOfLoad(NumberOfConcentratedLoads), LoadDirection(NumberOfConcentratedLoads)
  REAL(8) :: R(NumberOfEquations), FLOAD(NumberOfConcentratedLoads)
  INTEGER :: I,L,LI,LN,II

  WRITE (OutputFile,"(/,'    NODE       DIRECTION      LOAD',/, '   NUMBER',19X,'MAGNITUDE')")

  READ (InputFile,"(2I5,F10.0)") &
    (NodeOfLoad(I), LoadDirection(I), FLOAD(I), I=1, NumberOfConcentratedLoads)

  WRITE (OutputFile,"(' ',I6,9X,I4,7X,E12.5)") &
    (NodeOfLoad(I), LoadDirection(I), FLOAD(I),I=1, NumberOfConcentratedLoads)

  IF (SolutionMode.EQ.0) RETURN

  !DO I=1,NumberOfEquations
     R(:)=0.
  !END DO

  DO I = 1,NumberOfConcentratedLoads
     LN = NodeOfLoad(I)
     LI = LoadDirection(I)
     II = NodeIdentifier(LI,LN)
     IF (II > 0) R(II)=R(II) + FLOAD(I)
  END DO

  WRITE (LoadTmpFile, Rec=L) R
  !WRITE(*,*) "R",R
  RETURN
  
END SUBROUTINE LOADS


SUBROUTINE ObtainLoadVector (R, NumberOfEquations, CurrentLoadCase)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To obtain the load vector                                       .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
  USE GLOBALS, ONLY : LoadTmpFile

  IMPLICIT NONE
  INTEGER :: CurrentLoadCase, NumberOfEquations
  REAL(8) :: R(NumberOfEquations)

  READ (LoadTmpFile, Rec=CurrentLoadCase) R
  RETURN
END SUBROUTINE ObtainLoadVector
