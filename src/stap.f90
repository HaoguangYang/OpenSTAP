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

PROGRAM STAP90

  USE GLOBALS
  USE memAllocate

  IMPLICIT NONE
  INTEGER :: NumberOfConcentratedLoads, MeanHalfBandwidth
  INTEGER :: CurrentLoadCase, LoadCaseIndex, I
  REAL :: TotalTime

! OPEN INPUT DATA FILE, RESULTS OUTPUT FILE AND TEMPORARY FILES
  CALL OpenFiles()

  MaxElementGroupArraySize=0

! * * * * * * * * * * * * * * * * * * * * * *
! *              INPUT PHASE                *
! * * * * * * * * * * * * * * * * * * * * * *

  WRITE(*,'("Input phase ... ")')

  CALL SECOND (Timing(1))

! Read control information

!   JobIdentifier    - The master heading informaiton for use in labeling the output
!   NumberOfNodalPoints  - Total number of nodal points
!            0 : program stop
!   NumberOfElementGroups  - Total number of element group (>0)
!   NumberOfLoadCases - Number of load case (>0)
!   SolutionMode  - Solution mode
!            0 : data check only;
!            1 : execution

  READ (InputFile,'(A80,/,4I5)') JobIdentifier, NumberOfNodalPoints, NumberOfElementGroups, NumberOfLoadCases, SolutionMode

  IF (NumberOfNodalPoints.EQ.0) STOP   ! Data check mode

  WRITE (OutputFile,"(/,' ',A80,//,  &
     ' C O N T R O L   I N F O R M A T I O N',//,  &
     '      NUMBER OF NODAL POINTS',10(' .'),' (NumberOfNodalPoints)  = ',I5,/,   &
     '      NUMBER OF ELEMENT GROUPS',9(' .'),' (NumberOfElementGroups)  = ',I5,/,  &
     '      NUMBER OF LOAD CASES',11(' .'),' (NumberOfLoadCases) = ',I5,/,     &
     '      SOLUTION MODE ',14(' .'),' (SolutionMode)  = ',I5,/,           &
     '         EQ.0, DATA CHECK',/,   &
     '         EQ.1, EXECUTION')") JobIdentifier, NumberOfNodalPoints, NumberOfElementGroups, NumberOfLoadCases, SolutionMode

! Read nodal point data

! ALLOCATE STORAGE
!   ID(3,NumberOfNodalPoints) : Boundary condition codes (0=free,1=deleted)
!   X(NumberOfNodalPoints)    : X coordinates
!   Y(NumberOfNodalPoints)    : Y coordinates
!   Z(NumberOfNodalPoints)    : Z coordinates

  CALL memAlloc(1,"ID   ", 3* NumberOfNodalPoints, 1)
  CALL memAlloc(2,"X    ", NumberOfNodalPoints, ITWO)
  CALL memAlloc(3,"Y    ", NumberOfNodalPoints, ITWO)
  CALL memAlloc(4,"Z    ", NumberOfNodalPoints, ITWO)

  CALL INPUT (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)), NumberOfNodalPoints, NumberOfEquations)
    CALL OpenTempFiles()
! Calculate and store load vectors
!   R(NumberOfEquations) : Load vector

  CALL memAlloc(5,"R    ",NumberOfEquations,ITWO)

  WRITE (OutputFile,"(//,' L O A D   C A S E   D A T A')")

  !REWIND LoadTmpFile

  DO CurrentLoadCase = 1,NumberOfLoadCases

!    LL    - Load case number
!    NumberOfConcentratedLoads - The number of concentrated loads applied in this load case

     READ (InputFile,'(2I5)') LoadCaseIndex , NumberOfConcentratedLoads

     WRITE (OutputFile,"(/,'     LOAD CASE NUMBER',7(' .'),' = ',I5,/, &
                     '     NUMBER OF CONCENTRATED LOADS . = ',I5)") LoadCaseIndex , NumberOfConcentratedLoads

     IF (LoadCaseIndex.NE.CurrentLoadCase) THEN
        WRITE (OutputFile,"(' *** ERROR *** LOAD CASES ARE NOT IN ORDER')")
        STOP
     ENDIF

!    Allocate storage
!       NOD(NumberOfConcentratedLoads)   : Node number to which this load is applied (1~NumberOfNodalPoints)
!       IDIRN(NumberOfConcentratedLoads) : Degree of freedom number for this load component
!                      1 : X-direction;
!                      2 : Y-direction;
!                      3 : Z-direction
!       FLOAD(NumberOfConcentratedLoads) : Magnitude of load

     CALL memAlloc(6,"NOD  ",NumberOfConcentratedLoads,1)
     CALL memAlloc(7,"IDIRN",NumberOfConcentratedLoads,1)
     CALL memAlloc(8,"FLOAD",NumberOfConcentratedLoads,ITWO)

     CALL LOADS (DA(NP(5)),IA(NP(6)),IA(NP(7)),DA(NP(8)),IA(NP(1)), CurrentLoadCase, &
                 NumberOfConcentratedLoads, NumberOfEquations)

  END DO

! Read, generate and store element data

! Clear storage
!   MHT(NumberOfEquations) - Vector of column heights

  CALL memFreeFrom(5)
  CALL memAlloc(5,"MHT  ",NumberOfEquations,1)

  SolutionPhase=1    ! Read and generate element information
  CALL ELCAL

  CALL SECOND (Timing(2))

! * * * * * * * * * * * * * * * * * * * * * *
! *               SOLUTION PHASE            *
! * * * * * * * * * * * * * * * * * * * * * *

  WRITE(*,'("Solution phase ... ")')

! Assemble stiffness matrix

! ALLOCATE STORAGE
!    MAXA(NumberOfEquations+1)
  CALL memFreeFrom(6)
  CALL memFreeFromTo(2,4)
  CALL memAlloc(2,"MAXA ",NumberOfEquations+1,1)

  CALL ADDRES (IA(NP(2)),IA(NP(5)))

! ALLOCATE STORAGE
!    A(NWK) - Global structure stiffness matrix K
!    R(NumberOfEquations) - Load vector R and then displacement solution U

  MeanHalfBandwidth = NumberOfMatrixElements / NumberOfEquations

  CALL memAlloc(3,"STFF ",NumberOfMatrixElements,ITWO)
  CALL memAlloc(4,"R    ",NumberOfEquations,ITWO)
  CALL memAlloc(11,"ELEGP",MaxElementGroupArraySize,1)

! Write total system data

  WRITE (OutputFile,"(//,' TOTAL SYSTEM DATA',//,   &
                   '     NUMBER OF EQUATIONS',14(' .'),'(NEQ) = ',I5,/,   &
                   '     NUMBER OF MATRIX ELEMENTS',11(' .'),'(NWK) = ',I5,/,   &
                   '     MAXIMUM HALF BANDWIDTH ',12(' .'),'(MK) = ',I5,/,     &
                   '     MEAN HALF BANDWIDTH',14(' .'),'(MM ) = ',I5)" )     &
                   NumberOfEquations, NumberOfMatrixElements, MaxHalfBandwidth, MeanHalfBandwidth

! In data check only mode we skip all further calculations

  IF (SolutionMode.LE.0) THEN
     CALL SECOND (Timing(3))
     CALL SECOND (Timing(4))
     CALL SECOND (Timing(5))
  ELSE
     SolutionPhase=2    ! Assemble structure stiffness matrix
     CALL ASSEM (A(NP(11)))                 !SolutionPhase=2 RUSS Entrance

     CALL SECOND (Timing(3))

!    Triangularize stiffness matrix
     CALL COLSOL (DA(NP(3)),DA(NP(4)),IA(NP(2)), NumberOfEquations, NumberOfMatrixElements, NumberOfEquations+1, 1)

     CALL SECOND (Timing(4))

     SolutionPhase=3    ! Stress calculations

     !REWIND LoadTmpFile
     DO CurrentLoadCase = 1,NumberOfLoadCases
        CALL ObtainLoadVector (DA(NP(4)),NumberOfEquations, CurrentLoadCase)   ! Read in the load vector
        
!       Solve the equilibrium equations to calculate the displacements
        CALL COLSOL (DA(NP(3)),DA(NP(4)),IA(NP(2)), NumberOfEquations, NumberOfMatrixElements, NumberOfEquations+1, 2)

        WRITE (OutputFile,"(//,' LOAD CASE ',I3)") CurrentLoadCase
        CALL WRITED (DA(NP(4)),IA(NP(1)), NumberOfEquations, NumberOfNodalPoints)  ! Print displacements

!       Calculation of stresses
        CALL STRESS (A(NP(11)))

     END DO

     CALL SECOND (Timing(5))
  END IF

! Print solution times

  TotalTime=0.
  DO I=1,4
     Timing(I) = Timing(I+1) - Timing(I)
     TotalTime = TotalTime + Timing(I)
  END DO

  WRITE (OutputFile,"(//,  &
     ' S O L U T I O N   T I M E   L O G   I N   S E C',//,   &
     '     TIME FOR INPUT PHASE ',14(' .'),' =',F12.2,/,     &
     '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . =',F12.2, /,   &
     '     TIME FOR FACTORIZATION OF STIFFNESS MATRIX  . . . =',F12.2, /,   &
     '     TIME FOR LOAD CASE SOLUTIONS ',10(' .'),' =',F12.2,//,   &
     '      T O T A L   S O L U T I O N   T I M E  . . . . . =',F12.2)") (Timing(I),I=1,4),TotalTime


  WRITE (*,"(//,  &
     ' S O L U T I O N   T I M E   L O G   I N   S E C',//,   &
     '     TIME FOR INPUT PHASE ',14(' .'),' =',F12.2,/,     &
     '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . =',F12.2, /,   &
     '     TIME FOR FACTORIZATION OF STIFFNESS MATRIX  . . . =',F12.2, /,   &
     '     TIME FOR LOAD CASE SOLUTIONS ',10(' .'),' =',F12.2,//,   &
     '      T O T A L   S O L U T I O N   T I M E  . . . . . =',F12.2)") (Timing(I),I=1,4),TotalTime
  STOP

END PROGRAM STAP90


SUBROUTINE SECOND (Timing)
! USE DFPORT   ! Only for Compaq Fortran
  IMPLICIT NONE
  REAL :: Timing

! This is a Fortran 95 intrinsic subroutine
! Returns the processor time in seconds

  CALL CPU_TIME(Timing)

  RETURN
END SUBROUTINE SECOND


SUBROUTINE WRITED (Displacement, Ident, NumberOfEquations, NumberOfNodalPoints)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To print displacements                                          .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : OutputFile

  IMPLICIT NONE
  INTEGER :: NumberOfEquations,NumberOfNodalPoints,Ident(3,NumberOfNodalPoints)
  REAL(8) :: Displacement(NumberOfEquations),DispToBePrinted(3)
  INTEGER :: RowBanner,NodeIndex,I,FixedOrNot

! Print displacements

  WRITE (OutputFile,"(//,' D I S P L A C E M E N T S',//,'  NODE ',10X,   &
                    'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')")

  RowBanner=4

  DO NodeIndex=1,NumberOfNodalPoints
     RowBanner = RowBanner + 1
     IF (RowBanner.GE.56) THEN
        WRITE (OutputFile,"(//,' D I S P L A C E M E N T S',//,'  NODE ',10X,   &
                          'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')")
        RowBanner = 4
     END IF

     DispToBePrinted(1:3)=0

     DO I=1,3
        FixedOrNot=Ident(I,NodeIndex)
        IF (FixedOrNot.NE.0) DispToBePrinted(I) = Displacement(FixedOrNot)
     END DO

     WRITE (OutputFile,'(1X,I3,8X,3E18.6)') NodeIndex,DispToBePrinted

  END DO

  RETURN

END SUBROUTINE WRITED


SUBROUTINE OpenFiles()
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   Open input data file, results output file and temporary files   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  use memAllocate
! use DFLIB ! for NARGS()  ! Only for Compaq Fortran

  IMPLICIT NONE
  LOGICAL :: IfExist
  CHARACTER*80 FileName

! Only for Compaq Fortran
! if(NARGS().ne.2) then
!    stop 'Usage: mpm3d InputFileName'
!  else
!    call GETARG(1,FileInp)
!  end if

  if(COMMAND_ARGUMENT_COUNT().ne.1) then
     stop 'Usage: STAP90 InputFileName'
  else
     call GET_COMMAND_ARGUMENT(1,FileName)
  end if

  INQUIRE(FILE = FileName, EXIST = IfExist)
  IF (.NOT. IfExist) THEN
     PRINT *, "*** STOP *** FILE", FileName, "DOES NOT EXIST !"
     STOP
  END IF

  OPEN(InputFile   , FILE = FileName,  STATUS = "OLD")
  OPEN(OutputFile  , FILE = "STAP90.OUT", STATUS = "REPLACE")
END SUBROUTINE OpenFiles


SUBROUTINE OpenTempFiles()
use GLOBALS
use memAllocate
implicit none
OPEN(LoadTmpFile , FILE = "LOAD.TMP",   FORM = "UNFORMATTED",  Access='Direct',  RecL=6*ITWO*NumberOfEquations)
OPEN(ElementTmpFile, FILE = "ELMNT.TMP",  FORM = "UNFORMATTED")
end subroutine OpenTempFiles


SUBROUTINE CLOSEFILES()
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   Close all data files                                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  IMPLICIT NONE
  CLOSE(InputFile)
  CLOSE(OutputFile)
  CLOSE(ElementTmpFile)
  CLOSE(LoadTmpFile)
END SUBROUTINE CLOSEFILES
