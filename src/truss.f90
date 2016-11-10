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

SUBROUTINE TRUSS
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the truss element subroutine         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NumberOfElements, NumberOfMaterials, MM, N101, N102, N103, N104, N105, N106, N107, N108, N109

  NumberOfElements = NPAR(2)
  NumberOfMaterials = NPAR(3)

! Allocate storage for element group data
  IF (SolutionPhase == 1) THEN
      MM = 4 * NumberOfMaterials * ITWO + 9 * NumberOfElements + 6 * NumberOfElements * ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)
!pointer lists
! Calculate the pointer to the arrays in the element group data
! N101: E(NumberOfMaterials)
! N102: AREA(NumberOfMaterials)
! N103: Density(NumberOfMaterials)
! N104: Gravity(NumberOfMaterials)
! N105: LM(6,NumberOfElements)
! N106: PositionData(6,NumberOfElements)
! N107: MTAP(NumberOfElements)
! N108: Element Length(NumberOfElements)
  N101=NFIRST
  N102=N101+NumberOfMaterials*ITWO
  N103=N102+NumberOfMaterials*ITWO
  N104=N103+NumberOfMaterials*ITWO
  N105=N104+NumberOfMaterials*ITWO
  N106=N105+6*NumberOfElements
  N107=N106+6*NumberOfElements*ITWO
  N108=N107+NumberOfElements
  N109=N108+NumberOfElements*2
  NLAST=N109

  ElementGroupArraySize = NLAST - NFIRST

  CALL RUSS (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105),A(N106),A(N107),A(N108))

  RETURN

END SUBROUTINE TRUSS


SUBROUTINE RUSS (ID,X,Y,Z,U,MHT,E,Area ,Density, Gravity, LM, PositionData, MaterialData, LengthSquared)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(3,NumberOfNodalPoints), LM(6,NPAR(2)), MaterialData(NPAR(2)), MHT(NumberOfEquations)
  REAL(8) :: X(NumberOfNodalPoints),Y(NumberOfNodalPoints),Z(NumberOfNodalPoints),E(NPAR(3)),AREA(NPAR(3)),  &
             PositionData(6,NPAR(2)),U(NumberOfEquations), Density(NPAR(3)), Gravity(NPAR(3)), LengthSquared(NPAR(2))
  REAL(8) :: S(6,6),ST(6),D(3)

  INTEGER :: ElementType, NumberOfElements, NumberOfMaterials, ND, I, J, L, N, M
  INTEGER :: MaterialType, IPRINT
  REAL(8) :: XX, YY, Stress, Force

  ElementType       = NPAR(1)
  NumberOfElements  = NPAR(2)
  NumberOfMaterials = NPAR(3) 

  ND=6

! Read and generate element information
  IF (SolutionPhase .EQ. 1) THEN

     WRITE (OutputFile,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.1, TRUSS ELEMENTS',/,      &
                   '     EQ.2, ELEMENTS CURRENTLY',/,  &
                   '     EQ.3, NOT AVAILABLE',//,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") ElementType, NumberOfElements

     IF (NumberOfMaterials.EQ.0) NumberOfMaterials=1

     WRITE (OutputFile,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I5,/)") NumberOfMaterials

! Added Gravity Input.
     WRITE (OutputFile,"('  SET       YOUNG''S     CROSS-SECTIONAL    DENSITY    GRAVITY',/,  &
                   ' NUMBER     MODULUS',10X,'AREA',/,  &
                   15 X,'E',14X,'A')")

     DO I=1,NumberOfMaterials
        READ (InputFile,'(I5,4F10.0)') N, E(N), Area(N), Density(N), Gravity(N) ! Read material information
        WRITE (OutputFile,"(I5,4X,E12.5,2X,E14.6,2X,E12.5,2X,E12.5)") N, E(N), Area(N), Density(N), Gravity(N)
     END DO
     
     WRITE (OutputFile,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE       MATERIAL',/,   &
                      ' NUMBER-N      I        J       SET NUMBER')")

     N=0
     DO WHILE (N .NE. NumberOfElements)
        READ (InputFile,'(5I5)') N,I,J,MaterialType  ! Read in element information

!       Save element information
        PositionData(1,N)=X(I)  ! Coordinates of the element's left node
        PositionData(2,N)=Y(I)
        PositionData(3,N)=Z(I)

        PositionData(4,N)=X(J)  ! Coordinates of the element's right node
        PositionData(5,N)=Y(J)
        PositionData(6,N)=Z(J)

        MaterialData(N) = MaterialType  ! Material type
        
        DO L=1,6
           LM(L,N)=0
        END DO

        DO L=1,3
           LM(L,N)=ID(L,I)     ! Connectivity matrix
           LM(L+3,N)=ID(L,J)
           D(L)=PositionData(L,N) - PositionData(L+3,N)     !Length^2 of Element
        END DO
        
        LengthSquared(N)=dot_product(D,D)
        
!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (OutputFile,"(I5,6X,I5,4X,I5,7X,I5)") N,I,J,MaterialType

     END DO
     !write(*,*) LengthSquared
     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (SolutionPhase .EQ. 2) THEN

     DO N=1,NumberOfElements
        MaterialType = MaterialData(N)

     XX=E(MaterialType)*Area(MaterialType)*SQRT(LengthSquared(N))   !  E*A*l

     DO L=1,3
        !write(*,*) PositionData(L,N), PositionData(L+3,N), LengthSquared(N)
        ST(L)=(PositionData(L,N) - PositionData(L+3,N))/LengthSquared(N)
        !write(*,*) ST(L)
        ST(L+3)=-ST(L)
     END DO
        
     DO J=1,ND
        YY=ST(J)*XX
        DO I=1,J
           S(I,J)=ST(I)*YY
        END DO
     END DO

     CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)

     END DO
     
     !Reading Load Data, Adding Gravity, and Storing Data to Disk.
    !rewind LoadTmpFile
    do L=1,NumberOfLoadCases
        CALL ObtainLoadVector (DA(NP(4)),NumberOfEquations,L)
        !WRITE(*,*) DA(NP(4))
        DO N=1,NumberOfElements
            CALL LoadGravity (DA(NP(4)), N, LM, MaterialData(N), SQRT(LengthSquared(N)), Area, Density, Gravity)
        ENDDO
        !CALL MEMPRINT(4)
        write(LoadTmpFile, Rec=L) DA(NP(4):NP(4)+NumberOfEquations)       !???
     enddo

     RETURN

! Stress calculations
  ELSE IF (SolutionPhase .EQ. 3) THEN

     IPRINT=0
     DO N=1,NumberOfElements
        IPRINT=IPRINT + 1
        IF (IPRINT.GT.50) IPRINT=1
        IF (IPRINT.EQ.1) WRITE (OutputFile,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  ',  &
                                           'E L E M E N T  G R O U P',I4,//,   &
                                           '  ELEMENT',13X,'FORCE',12X,'STRESS',/,'  NUMBER')") CurrentElementGroup
        MaterialType = MaterialData(N)

        DO L=1,3
           ST(L)=((PositionData(L,N) - PositionData(L+3,N))/LengthSquared(N))*E(MaterialType)
           ST(L+3)=-ST(L)
        END DO

        Stress=0.0
        DO L=1,3
           I=LM(L,N)
           IF (I.GT.0) Stress=Stress + ST(L)*U(I)

           J=LM(L+3,N)
           IF (J.GT.0) Stress=Stress + ST(L+3)*U(J)
        END DO

        Force=Stress*AREA(MaterialType)

        WRITE (OutputFile,"(1X,I5,11X,E13.6,4X,E13.6)") N,Force,Stress
     END DO

  ELSE 
     STOP "*** ERROR *** Invalid SolutionPhase value."
  END IF

END SUBROUTINE RUSS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!             Add Gravity to Load                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LoadGravity (Load, ElementNumber, LocationMatrix, MaterialType, Length, Area, Density, Gravity)

use globals
USE MEMALLOCATE
implicit none

integer:: ID, MaterialType, LocationMatrix(6,NPAR(2)), ElementNumber
real(8):: Load(NumberOfEquations), Length, Area(NPAR(3)), Density(NPAR(3)), Gravity(NPAR(3))

    ID=LocationMatrix(3,ElementNumber)
    !WRITE(*,*) id
    IF (ID > 0) then
        Load(ID)=Load(ID) - Gravity(MaterialType) * Density(MaterialType) * Area(MaterialType) * Length /2
        !WRITE(*,*) LOAD(ID)
    endif
    ID=LocationMatrix(6,ElementNumber)
    !WRITE(*,*) id
    IF (ID > 0) then
        Load(ID)=Load(ID) - Gravity(MaterialType) * Density(MaterialType) * Area(MaterialType) * Length /2
        !WRITE(*,*) LOAD(ID)
    endif
end subroutine LoadGravity
