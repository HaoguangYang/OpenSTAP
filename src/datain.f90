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

SUBROUTINE INPUT (ID,X,Y,Z,NUMNP,NEQ)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       .
! .   To read, generate, and print nodal point input data                 .
! .   To calculate equation numbers and store them in id arrray           .
! .                                                                       .
! .      N = Element number                                               .
! .      ID = Boundary condition codes (0=free,1=deleted)                 .
! .      X,Y,Z = Coordinates                                              .
! .      KN = Generation code                                             .
! .           i.e. increment on nodal point number                        .
! .                                                                       .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : IIN, IOUT

  IMPLICIT NONE
  INTEGER :: NUMNP,NEQ,ID(3,NUMNP)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP)
  INTEGER :: I, J, N

! Read and generate nodal point data

  N = 0
  DO WHILE (N.NE.NUMNP)
     READ (IIN,"(4I5,3F10.0,I5)") N,(ID(I,N),I=1,3),X(N),Y(N),Z(N)
  END DO

! Write complete nodal data

  WRITE (IOUT,"(//,' N O D A L   P O I N T   D A T A',/)")

  WRITE (IOUT,"('  NODE',10X,'BOUNDARY',25X,'NODAL POINT',/,  &
                ' NUMBER     CONDITION  CODES',21X,'COORDINATES', /,15X, &
                'X    Y    Z',15X,'X',12X,'Y',12X,'Z')")

  DO N=1,NUMNP
     WRITE (IOUT,"(I5,6X,3I5,6X,3F13.3)") N,(ID(I,N),I=1,3),X(N),Y(N),Z(N)
  END DO

! Number unknowns

  NEQ=0
  DO N=1,NUMNP
     DO I=1,3
        IF (ID(I,N) .EQ. 0) THEN
           NEQ=NEQ + 1
           ID(I,N)=NEQ
        ELSE
           ID(I,N)=0
        END IF
     END DO
  END DO

! Write equation numbers
  WRITE (IOUT,"(//,' EQUATION NUMBERS',//,'   NODE',9X,  &
                   'DEGREES OF FREEDOM',/,'  NUMBER',/,  &
                   '     N',13X,'X    Y    Z',/,(1X,I5,9X,3I5))") (N,(ID(I,N),I=1,3),N=1,NUMNP)

  RETURN

END SUBROUTINE INPUT


SUBROUTINE LOADS (R,NOD,IDIRN,FLOAD,ID,NLOAD,NEQ)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To read nodal load data                                         .
! .   To calculate the load vector r for each load case and           .
! .   write onto unit ILOAD                                           .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  USE GLOBALS, ONLY : IIN, IOUT, ILOAD, MODEX

  IMPLICIT NONE
  INTEGER :: NLOAD,NEQ,ID(3,*),NOD(NLOAD),IDIRN(NLOAD)
  REAL(8) :: R(NEQ),FLOAD(NLOAD)
  INTEGER :: I,L,LI,LN,II

  WRITE (IOUT,"(/,'    NODE       DIRECTION      LOAD',/, '   NUMBER',19X,'MAGNITUDE')")

  READ (IIN,"(2I5,F10.0)") (NOD(I),IDIRN(I),FLOAD(I),I=1,NLOAD)

  WRITE (IOUT,"(' ',I6,9X,I4,7X,E12.5)") (NOD(I),IDIRN(I),FLOAD(I),I=1,NLOAD)

  IF (MODEX.EQ.0) RETURN

  DO I=1,NEQ
     R(I)=0.
  END DO

  DO L=1,NLOAD
     LN=NOD(L)
     LI=IDIRN(L)
     II=ID(LI,LN)
     IF (II > 0) R(II)=R(II) + FLOAD(L)
  END DO

  WRITE (ILOAD) R

  RETURN
  
END SUBROUTINE LOADS


SUBROUTINE LOADV (R,NEQ)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To obtain the load vector                                       .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
  USE GLOBALS, ONLY : ILOAD

  IMPLICIT NONE
  INTEGER :: NEQ
  REAL(8) :: R(NEQ)

  READ (ILOAD) R
  
  RETURN
END SUBROUTINE LOADV
