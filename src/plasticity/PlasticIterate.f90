SUBROUTINE PLASTICITERATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        PLASTIC ITERATION                   !
!         LIU CHANGWU                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NEQ1, NLOAD, MM
  INTEGER :: LL, I
  INTEGER :: TT

  
! * * * * * * * * * * * * * * * * * * * * * *
! *               SOLUTION PHASE            *
! * * * * * * * * * * * * * * * * * * * * * *  

WRITE(*,'("Solution phase ... ")')
IND=1    ! Read and generate element information

CALL MEMFREEFROM(5)
CALL MEMALLOC(5,"MHT  ",NEQ,1)                  !if (.NOT. PARDISODOOR)
CALL ELCAL ! 到这里2,3,4才没用的
!CALL VTKgenerate (IND)        !Prepare Post-Processing Files.

! ********************************************************************8
! Read, generate and store element data
! 从这里开始，用不用pardiso会变得很不一样
! Clear storage
!   MHT(NEQ) - Vector of column heights

if(.not. pardisodoor) then
  CALL SECOND (TIM(2)) 
! ALLOCATE STORAGE
!    MAXA(NEQ+1)
  CALL MEMFREEFROM(7)
  CALL MEMFREEFROMTO(2,4)
  CALL MEMALLOC(2,"MAXA ",NEQ+1,1)
  CALL ADDRES (IA(NP(2)),IA(NP(5)))
  CALL SECOND (TIM(3))
! ALLOCATE STORAGE
!    A(NWK) - Global structure stiffness matrix K
!    R(NEQ) - Load vector R and then displacement solution U
 
  MM=NWK/NEQ
  CALL MEMALLOC(3,"STFF ",NWK,ITWO)
  CALL MEMALLOC(4,"R    ",NEQ,ITWO)
  CALL MEMALLOC(11,"ELEGP",MAXEST,1)

! Write total system data

  WRITE (IOUT,"(//,' TOTAL SYSTEM DATA',//,   &
                   '     NUMBER OF EQUATIONS',14(' .'),'(NEQ) = ',I10,/,   &
                   '     NUMBER OF MATRIX ELEMENTS',11(' .'),'(NWK) = ',I10,/,   &
                   '     MAXIMUM HALF BANDWIDTH ',12(' .'),'(MK ) = ',I10,/,     &
                   '     MEAN HALF BANDWIDTH',14(' .'),'(MM ) = ',I10)") NEQ,NWK,MK,MM
! ***************************************************************************************
else !如果使用pardiso
  CALL MEMFREEFROMTO(2,4)
  ! NP(2,3,4,5)均在这里被分配
  CALL pardiso_input(IA(NP(1)))
  CALL SECOND (TIM(3))
  CALL MEMALLOC(11,"ELEGP",MAXEST,1)
! Write total system data
end if

IF (DYNANALYSIS .EQV. .TRUE.) call prepare_MassMatrix

! In data check only mode we skip all further calculations
  IF (MODEX.LE.0) THEN
     CALL SECOND (TIM(4))
     CALL SECOND (TIM(5))
     CALL SECOND (TIM(6))
  ELSE
     IND=2    ! Assemble structure stiffness matrix
     CALL ASSEM (A(NP(11)))
     
     CALL SECOND (TIM(4))
     IF (DYNANALYSIS .EQV. .TRUE.) CALL EIGENVAL (DA(NP(3)), DA(NP(10)), IA(NP(2)), NEQ, NWK, NEQ1, 2)
     if(pardisodoor) then
        if (.not. DYNANALYSIS) call pardiso_crop(DA(NP(3)), IA(NP(2)), IA(NP(5)))          ! Condensing CSR format sparse matrix storage: deleting zeros
     else
        !    Triangularize stiffness matrix
        NEQ1=NEQ + 1
        CALL COLSOL (DA(NP(3)),DA(NP(4)),IA(NP(2)),NEQ,NWK,NEQ1,1)
     end if
     
      
     IND=3    ! Stress calculations

     REWIND ILOAD
     CALL SECOND (TIM(5))
     DO CURLCASE=1,NLCASE
        CALL LOADV (DA(NP(4)),NEQ)   ! Read in the load vector
        if(pardisodoor) then
            WRITE (IOUT,"(//,' TOTAL SYSTEM DATA',//,   &
                   '     NUMBER OF EQUATIONS',14(' .'),'(NEQ) = ',I10,/,   &
                   '     NUMBER OF MATRIX ELEMENTS',11(' .'),'(NWK) = ',I9)") NEQ,NWK
            call pardiso_solver(DA(NP(3)),DA(NP(4)),IA(NP(2)), IA(NP(5)))
        else
!       Solve the equilibrium equations to calculate the displacements
            IF (LOADANALYSIS .EQV. .TRUE.) CALL COLSOL (DA(NP(3)),DA(NP(4)),IA(NP(2)),NEQ,NWK,NEQ1,2)
        end if
        CALL SECOND (TIM(6))
        WRITE (IOUT,"(//,' LOAD CASE ',I3)") CURLCASE
        
        CALL WRITED (DA(NP(4)),IA(NP(1)),NEQ,NUMNP)  ! PRINT DISPLACEMENTS FOR OTHER SITUATIONS(THE FORMER ONE)
!           Calculation of stresses
            CALL STRESS (A(NP(11)))
            CALL SECOND (TIM(7))
     END DO
     

     
!     CALL VTKgenerate (IND)
     
  END IF

    
    
ENDSUBROUTINE PLASTICITERATE