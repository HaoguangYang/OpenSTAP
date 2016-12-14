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

  REWIND IELMNT
  WRITE (IOUT,"(//,' E L E M E N T   G R O U P   D A T A',//)")

! Loop over all element groups
! every element group represent for a kind of element
  DO N=1,NUMEG
     IF (N.NE.1) WRITE (IOUT,'(1X)')

     READ (IIN,'(10I5)') NPAR

     CALL ELEMNT

     IF (MIDEST.GT.MAXEST) MAXEST=MIDEST

     WRITE (IELMNT) MIDEST,NPAR,(A(I),I=NFIRST,NLAST)

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
  INTEGER :: NPAR1

  NPAR1=NPAR(1)

  SELECT CASE (NPAR1)
!    Other element types would be called here, identifying each
!    element type by a different NPAR(1) parameter
  case (1)
     CALL TRUSS
  case (2)
     CALL ELEMENT_4Q
  case (3)
     CALL ELEMENT_9Q
  case (4)
     CALL HEXAHEDRAL
  case (5)
     CALL BEAM
  case (6)
     CALL PLATE
  case (7)
     CALL SHELL
  case (8)
     CALL PLATE8Q
  case (9)
     CALL SHELL8Q
  case (10)
     CALL ELEMENT_3T
!    Other element types would be called here, identifying each
!    element type by a different NPAR(1) parameter
  case (11:)
     stop "ELEMENT TYPE STILL UNDER DEVELOPMENT..."
  END SELECT

  RETURN
END SUBROUTINE ELEMNT


SUBROUTINE STRESS (AA)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To call the element subroutine for the calculation of stresses  .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : IELMNT, NG, MIDEST, NPAR, NUMEG

  IMPLICIT NONE
  REAL :: AA(*)
  INTEGER N, I

! Loop over all element groups

  REWIND IELMNT

  DO N=1,NUMEG
     NG=N

     READ (IELMNT) MIDEST,NPAR,(AA(I),I=1,MIDEST)

     CALL ELEMNT
  END DO

  RETURN
END subroutine STRESS
