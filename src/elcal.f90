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

SUBROUTINE ELCAL
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To loop over all element groups for reading,                    .
! .   generating and storing the element data                         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: N, I

  REWIND ElementTmpFile
  WRITE (OutputFile,"(//,' E L E M E N T   G R O U P   D A T A',//)")

! Loop over all element groups

  DO N=1,NumberOfElementGroups
     IF (N.NE.1) WRITE (OutputFile,'(1X)')

     READ (InputFile,'(10I5)') NPAR

     CALL ELEMNT

     IF (ElementGroupArraySize.GT.MaxElementGroupArraySize) MaxElementGroupArraySize = ElementGroupArraySize

     WRITE (ElementTmpFile) ElementGroupArraySize,NPAR,(A(I),I=NFIRST,NLAST)

  END DO

  RETURN

END SUBROUTINE ELCAL


SUBROUTINE ELEMNT
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To call the appropriate element subroutine                      .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS

  IMPLICIT NONE
  INTEGER :: ElementType

  ElementType = NPAR(1)

  IF      (ElementType == 1) THEN
     CALL TRUSS
<<<<<<< HEAD
  ELSE IF (ElementType == 2) THEN
=======
  ELSE IF (NPAR1 == 2) THEN
>>>>>>> test_h
     CALL ELEMENT_4Q
!    Other element types would be called here, identifying each
!    element type by a different NPAR(1) parameter
  ELSE IF (.TRUE.) THEN
     stop "ELEMENT TYPE STILL UNDER DEVELOPMENT..."
  END IF

  RETURN
END SUBROUTINE ELEMNT


SUBROUTINE STRESS (ElementGroupData)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To call the element subroutine for the calculation of stresses  .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : ElementTmpFile, CurrentElementGroup, ElementGroupArraySize, NPAR, NumberOfElementGroups

  IMPLICIT NONE
  REAL :: ElementGroupData(*)
  INTEGER:: N, I

! Loop over all element groups

  REWIND ElementTmpFile

  DO N = 1,NumberOfElementGroups
     CurrentElementGroup = N

     READ (ElementTmpFile) ElementGroupArraySize,NPAR,(ElementGroupData(I),I=1,ElementGroupArraySize)

     CALL ELEMNT
  END DO

  RETURN
END subroutine STRESS
