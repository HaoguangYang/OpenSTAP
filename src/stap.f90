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
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NLCASE, NEQ1, NLOAD, MM
  INTEGER :: L, LL, I
  REAL :: TT

! OPEN INPUT DATA FILE, RESULTS OUTPUT FILE AND TEMPORARY FILES
  CALL OPENFILES()

  NUMEST=0
  MAXEST=0

! * * * * * * * * * * * * * * * * * * * * * *
! *              INPUT PHASE                *
! * * * * * * * * * * * * * * * * * * * * * *

  WRITE(*,'("Input phase ... ")')
  
  CALL SECOND (TIM(1))

! Read control information
  
!   HED    - The master heading informaiton for use in labeling the output
!   NUMNP  - Total number of nodal points
!            0 : program stop
!   NUMEG  - Total number of element group (>0)
!   NLCASE - Number of load case (>0)
!   MODEX  - Solution mode
!            0 : data check only;
!            1 : execution
  
  READ (IIN,'(A80,/,4I5)') HED,NUMNP,NUMEG,NLCASE,MODEX

  IF (NUMNP.EQ.0) STOP   ! Data check mode

  WRITE (IOUT,"(/,' ',A80,//,  &
     ' C O N T R O L   I N F O R M A T I O N',//,  &
     '      NUMBER OF NODAL POINTS',10(' .'),' (NUMNP)  = ',I5,/,   &
     '      NUMBER OF ELEMENT GROUPS',9(' .'),' (NUMEG)  = ',I5,/,  &
     '      NUMBER OF LOAD CASES',11(' .'),' (NLCASE) = ',I5,/,     &
     '      SOLUTION MODE ',14(' .'),' (MODEX)  = ',I5,/,            &
     '         EQ.0, DATA CHECK',/,   &
     '         EQ.1, EXECUTION')") HED,NUMNP,NUMEG,NLCASE,MODEX

! Read nodal point data

! ALLOCATE STORAGE
!   ID(3,NUMNP) : Boundary condition codes (0=free,1=deleted)
!   X(NUMNP)    : X coordinates
!   Y(NUMNP)    : Y coordinates
!   Z(NUMNP)    : Z coordinates

  CALL MEMALLOC(1,"ID   ",3*NUMNP,1)
  CALL MEMALLOC(2,"X    ",NUMNP,ITWO)
  CALL MEMALLOC(3,"Y    ",NUMNP,ITWO)
  CALL MEMALLOC(4,"Z    ",NUMNP,ITWO)

  CALL INPUT (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),NUMNP,NEQ)

  NEQ1=NEQ + 1

! Calculate and store load vectors
!   R(NEQ) : Load vector

  CALL MEMALLOC(5,"R    ",NEQ,ITWO)

  WRITE (IOUT,"(//,' L O A D   C A S E   D A T A')")

  REWIND ILOAD

  DO L=1,NLCASE

!    LL    - Load case number
!    NLOAD - The number of concentrated loads applied in this load case

     READ (IIN,'(2I5)') LL,NLOAD

     WRITE (IOUT,"(/,'     LOAD CASE NUMBER',7(' .'),' = ',I5,/, &
                     '     NUMBER OF CONCENTRATED LOADS . = ',I5)") LL,NLOAD

     IF (LL.NE.L) THEN
        WRITE (IOUT,"(' *** ERROR *** LOAD CASES ARE NOT IN ORDER')")
        STOP
     ENDIF

!    Allocate storage
!       NOD(NLOAD)   : Node number to which this load is applied (1~NUMNP)
!       IDIRN(NLOAD) : Degree of freedom number for this load component
!                      1 : X-direction;
!                      2 : Y-direction;
!                      3 : Z-direction
!       FLOAD(NLOAD) : Magnitude of load

     CALL MEMALLOC(6,"NOD  ",NLOAD,1)
     CALL MEMALLOC(7,"IDIRN",NLOAD,1)
     CALL MEMALLOC(8,"FLOAD",NLOAD,ITWO)

     CALL LOADS (DA(NP(5)),IA(NP(6)),IA(NP(7)),DA(NP(8)),IA(NP(1)),NLOAD,NEQ)

  END DO

! Read, generate and store element data

! Clear storage
!   MHT(NEQ) - Vector of column heights

  CALL MEMFREEFROM(5)
  CALL MEMALLOC(5,"MHT  ",NEQ,1)

  IND=1    ! Read and generate element information
  CALL ELCAL

  CALL SECOND (TIM(2))

! * * * * * * * * * * * * * * * * * * * * * *
! *               SOLUTION PHASE            *
! * * * * * * * * * * * * * * * * * * * * * *

  WRITE(*,'("Solution phase ... ")')
  
! Assemble stiffness matrix

! ALLOCATE STORAGE
!    MAXA(NEQ+1)
  CALL MEMFREEFROM(6)
  CALL MEMFREEFROMTO(2,4)
  CALL MEMALLOC(2,"MAXA ",NEQ+1,1)

  CALL ADDRES (IA(NP(2)),IA(NP(5)))

! ALLOCATE STORAGE
!    A(NWK) - Global structure stiffness matrix K
!    R(NEQ) - Load vector R and then displacement solution U

  MM=NWK/NEQ

  CALL MEMALLOC(3,"STFF ",NWK,ITWO)
  CALL MEMALLOC(4,"R    ",NEQ,ITWO)
  CALL MEMALLOC(11,"ELEGP",MAXEST,1)

! Write total system data

  WRITE (IOUT,"(//,' TOTAL SYSTEM DATA',//,   &
                   '     NUMBER OF EQUATIONS',14(' .'),'(NEQ) = ',I5,/,   &
                   '     NUMBER OF MATRIX ELEMENTS',11(' .'),'(NWK) = ',I5,/,   &
                   '     MAXIMUM HALF BANDWIDTH ',12(' .'),'(MK ) = ',I5,/,     &
                   '     MEAN HALF BANDWIDTH',14(' .'),'(MM ) = ',I5)") NEQ,NWK,MK,MM

! In data check only mode we skip all further calculations

  IF (MODEX.LE.0) THEN
     CALL SECOND (TIM(3))
     CALL SECOND (TIM(4))
     CALL SECOND (TIM(5))
  ELSE
     IND=2    ! Assemble structure stiffness matrix
     CALL ASSEM (A(NP(11)))

     CALL SECOND (TIM(3))

!    Triangularize stiffness matrix
     CALL COLSOL (DA(NP(3)),DA(NP(4)),IA(NP(2)),NEQ,NWK,NEQ1,1)

     CALL SECOND (TIM(4))

     IND=3    ! Stress calculations

     REWIND ILOAD
     DO L=1,NLCASE
        CALL LOADV (DA(NP(4)),NEQ)   ! Read in the load vector

!       Solve the equilibrium equations to calculate the displacements
        CALL COLSOL (DA(NP(3)),DA(NP(4)),IA(NP(2)),NEQ,NWK,NEQ1,2)

        WRITE (IOUT,"(//,' LOAD CASE ',I3)") L
        CALL WRITED (DA(NP(4)),IA(NP(1)),NEQ,NUMNP)  ! Print displacements

!       Calculation of stresses
        CALL STRESS (A(NP(11)))

     END DO

     CALL SECOND (TIM(5))
  END IF

! Print solution times

  TT=0.
  DO I=1,4
     TIM(I)=TIM(I+1) - TIM(I)
     TT=TT + TIM(I)
  END DO

  WRITE (IOUT,"(//,  &
     ' S O L U T I O N   T I M E   L O G   I N   S E C',//,   &
     '     TIME FOR INPUT PHASE ',14(' .'),' =',F12.2,/,     &
     '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . =',F12.2, /,   &
     '     TIME FOR FACTORIZATION OF STIFFNESS MATRIX  . . . =',F12.2, /,   &
     '     TIME FOR LOAD CASE SOLUTIONS ',10(' .'),' =',F12.2,//,   &
     '      T O T A L   S O L U T I O N   T I M E  . . . . . =',F12.2)") (TIM(I),I=1,4),TT


  WRITE (*,"(//,  &
     ' S O L U T I O N   T I M E   L O G   I N   S E C',//,   &
     '     TIME FOR INPUT PHASE ',14(' .'),' =',F12.2,/,     &
     '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . =',F12.2, /,   &
     '     TIME FOR FACTORIZATION OF STIFFNESS MATRIX  . . . =',F12.2, /,   &
     '     TIME FOR LOAD CASE SOLUTIONS ',10(' .'),' =',F12.2,//,   &
     '      T O T A L   S O L U T I O N   T I M E  . . . . . =',F12.2)") (TIM(I),I=1,4),TT
  STOP

END PROGRAM STAP90


SUBROUTINE SECOND (TIM)
! USE DFPORT   ! Only for Compaq Fortran
  IMPLICIT NONE
  REAL :: TIM

! This is a Fortran 95 intrinsic subroutine
! Returns the processor time in seconds

  CALL CPU_TIME(TIM)

  RETURN
END SUBROUTINE SECOND


SUBROUTINE WRITED (DISP,ID,NEQ,NUMNP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To print displacements                                          .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : IOUT

  IMPLICIT NONE
  INTEGER :: NEQ,NUMNP,ID(3,NUMNP)
  REAL(8) :: DISP(NEQ),D(3)
  INTEGER :: IC,II,I,KK,IL

! Print displacements

  WRITE (IOUT,"(//,' D I S P L A C E M E N T S',//,'  NODE ',10X,   &
                    'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')")

  IC=4

  DO II=1,NUMNP
     IC=IC + 1
     IF (IC.GE.56) THEN
        WRITE (IOUT,"(//,' D I S P L A C E M E N T S',//,'  NODE ',10X,   &
                          'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')")
        IC=4
     END IF

     DO I=1,3
        D(I)=0.
     END DO

     DO I=1,3
        KK=ID(I,II)
        IL=I
        IF (KK.NE.0) D(IL)=DISP(KK)
     END DO

     WRITE (IOUT,'(1X,I3,8X,3E18.6)') II,D

  END DO

  RETURN

END SUBROUTINE WRITED


SUBROUTINE OPENFILES()
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   Open input data file, results output file and temporary files   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
! use DFLIB ! for NARGS()  ! Only for Compaq Fortran

  IMPLICIT NONE
  LOGICAL :: EX
  CHARACTER*80 FileInp

! Only for Compaq Fortran
! if(NARGS().ne.2) then
!    stop 'Usage: mpm3d InputFileName'
!  else
!    call GETARG(1,FileInp)
!  end if

  if(COMMAND_ARGUMENT_COUNT().ne.1) then
     stop 'Usage: STAP90 InputFileName'
  else
     call GET_COMMAND_ARGUMENT(1,FileInp)
  end if

  INQUIRE(FILE = FileInp, EXIST = EX)
  IF (.NOT. EX) THEN
     PRINT *, "*** STOP *** FILE STAP90.IN DOES NOT EXIST !"
     STOP
  END IF

  OPEN(IIN   , FILE = FileInp,  STATUS = "OLD")
  OPEN(IOUT  , FILE = "STAP90.OUT", STATUS = "REPLACE")

  OPEN(IELMNT, FILE = "ELMNT.TMP",  FORM = "UNFORMATTED")
  OPEN(ILOAD , FILE = "LOAD.TMP",   FORM = "UNFORMATTED")
END SUBROUTINE OPENFILES


SUBROUTINE CLOSEFILES()
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   Close all data files                                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  IMPLICIT NONE
  CLOSE(IIN)
  CLOSE(IOUT)
  CLOSE(IELMNT)
  CLOSE(ILOAD)
END SUBROUTINE CLOSEFILES
