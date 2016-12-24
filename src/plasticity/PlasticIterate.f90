!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               弹塑性杆分析模块                              !  
!                                 作者：刘畅武                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PLASTIC 
  USE GLOBALS
  USE memAllocate
  
  IMPLICIT NONE
  INTEGER :: I
  
  DO I=1,FLOOR(ITERATENUM)
      CALL PLASTICITERATE
  ENDDO
  
  CALL FINALWRITE(IA(NP(1)))
  
ENDSUBROUTINE PLASTIC

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

!WRITE(*,'("Solution phase ... ")')
!IND=1    ! Read and generate element information

!CALL MEMFREEFROM(5)
!CALL MEMALLOC(5,"MHT  ",NEQ,1)                  !if (.NOT. PARDISODOOR)
!CALL ELCAL ! 到这里2,3,4才没用的
!CALL VTKgenerate (IND)        !Prepare Post-Processing Files.

! ********************************************************************8
! Read, generate and store element data
! 从这里开始，用不用pardiso会变得很不一样
! Clear storage
!   MHT(NEQ) - Vector of column heights

!if(.not. pardisodoor) then
!  CALL SECOND (TIM(2)) 
! ALLOCATE STORAGE
!    MAXA(NEQ+1)
  CALL MEMFREEFROM(7)
  CALL MEMFREEFROMTO(2,4)
  CALL MEMALLOC(2,"MAXA ",NEQ+1,1)
  CALL ADDRES (IA(NP(2)),IA(NP(5)))
!  CALL SECOND (TIM(3))
! ALLOCATE STORAGE
!    A(NWK) - Global structure stiffness matrix K
!    R(NEQ) - Load vector R and then displacement solution U
 
  MM=NWK/NEQ
  CALL MEMALLOC(3,"STFF ",NWK,ITWO)
  CALL MEMALLOC(4,"R    ",NEQ,ITWO)
  CALL MEMALLOC(11,"ELEGP",MAXEST,1)

! Write total system data

!  WRITE (IOUT,"(//,' TOTAL SYSTEM DATA',//,   &
!                   '     NUMBER OF EQUATIONS',14(' .'),'(NEQ) = ',I10,/,   &
!                   '     NUMBER OF MATRIX ELEMENTS',11(' .'),'(NWK) = ',I10,/,   &
!                   '     MAXIMUM HALF BANDWIDTH ',12(' .'),'(MK ) = ',I10,/,     &
!                  '     MEAN HALF BANDWIDTH',14(' .'),'(MM ) = ',I10)") NEQ,NWK,MK,MM
! ***************************************************************************************
!else !如果使用pardiso
!  CALL MEMFREEFROMTO(2,4)
  ! NP(2,3,4,5)均在这里被分配
!  CALL pardiso_input(IA(NP(1)))
!  CALL SECOND (TIM(3))
!  CALL MEMALLOC(11,"ELEGP",MAXEST,1)
! Write total system data
!end if
!****************************************************************************************


! In data check only mode we skip all further calculations
  IF (MODEX.LE.0) THEN
     CALL SECOND (TIM(4))
     CALL SECOND (TIM(5))
     CALL SECOND (TIM(6))
  ELSE
     IND=2    ! Assemble structure stiffness matrix
     CALL ASSEM (A(NP(11)))
     
!     CALL SECOND (TIM(4))
!     IF (DYNANALYSIS .EQV. .TRUE.) CALL EIGENVAL (DA(NP(3)), DA(NP(10)), IA(NP(2)), NEQ, NWK, NEQ1, 2)
!     if(pardisodoor) then
!        if (.not. DYNANALYSIS) call pardiso_crop(DA(NP(3)), IA(NP(2)), IA(NP(5)))          ! Condensing CSR format sparse matrix storage: deleting zeros
!     else
        !    Triangularize stiffness matrix
        NEQ1=NEQ + 1
        CALL COLSOL (DA(NP(3)),DA(NP(4)),IA(NP(2)),NEQ,NWK,NEQ1,1)
!     end if
     
      
     IND=3    ! Stress calculations

     REWIND DELTALOAD
!     CALL SECOND (TIM(5))
     DO CURLCASE=1,NLCASE
        CALL DELTALOADV (DA(NP(4)),NEQ)   ! Read in the load vector
!        if(pardisodoor) then
!           WRITE (IOUT,"(//,' TOTAL SYSTEM DATA',//,   &
!                   '     NUMBER OF EQUATIONS',14(' .'),'(NEQ) = ',I10,/,   &
!                   '     NUMBER OF MATRIX ELEMENTS',11(' .'),'(NWK) = ',I9)") NEQ,NWK
!            call pardiso_solver(DA(NP(3)),DA(NP(4)),IA(NP(2)), IA(NP(5)))
!        else
!       Solve the equilibrium equations to calculate the displacements
!            IF (LOADANALYSIS .EQV. .TRUE.) CALL COLSOL (DA(NP(3)),DA(NP(4)),IA(NP(2)),NEQ,NWK,NEQ1,2)
!        end if
!        CALL SECOND (TIM(6))
!        WRITE (IOUT,"(//,' LOAD CASE ',I3)") CURLCASE

        CALL COLSOL (DA(NP(3)),DA(NP(4)),IA(NP(2)),NEQ,NWK,NEQ1,2)
        
        CALL WRITEDISPLACEMENT(DA(NP(4)))
        
        CALL WRITED (DA(NP(4)),IA(NP(1)),NEQ,NUMNP)  ! PRINT DISPLACEMENTS FOR OTHER SITUATIONS(THE FORMER ONE)
!           Calculation of stresses
        CALL STRESS (A(NP(11)))
!            CALL SECOND (TIM(7))
     END DO
     
  END IF
 
ENDSUBROUTINE PLASTICITERATE
    
SUBROUTINE DELTALOADV (R,NEQ)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To obtain the load vector                                       .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
  USE GLOBALS, ONLY : DELTALOAD

  IMPLICIT NONE
  INTEGER :: NEQ
  REAL(8) :: R(NEQ)

  READ (DELTALOAD) R
  
  RETURN
END SUBROUTINE DELTALOADV 
    
    
SUBROUTINE WRITEDISPLACEMENT (DISP)


  USE GLOBALS, ONLY : PRESENTDISPLACEMENT,NEQ

  IMPLICIT NONE
  REAL(8) :: DISP(NEQ),DISP0(NEQ),DISP2(NEQ)
  INTEGER :: I

  REWIND PRESENTDISPLACEMENT
  READ(PRESENTDISPLACEMENT) DISP0
  
  DISP2=DISP+DISP0
  
  REWIND PRESENTDISPLACEMENT
  WRITE(PRESENTDISPLACEMENT) DISP2
  
  RETURN

END SUBROUTINE WRITEDISPLACEMENT   
    
SUBROUTINE FINALWRITE(ID)
  USE GLOBALS
  USE memAllocate
  
  IMPLICIT NONE
  INTEGER:: ID(6,NUMNP),I,IC,II,KK
  REAL(8):: DISP(NEQ),ELESTRESS(NPAR(2)),D(6)
 
  !输出最终位移
  REWIND PRESENTDISPLACEMENT
  READ(PRESENTDISPLACEMENT) (DISP(I),I=1,NEQ)
  
  
  WRITE (IOUT,"(//,' D I S P L A C E M E N T S',//,'  NODE ',3X,   &
                    'X-DISPLACEMENT  Y-DISPLACEMENT  Z-DISPLACEMENT  X-ROTATION  Y-ROTATION  Z-ROTATION')")

  IC=4

 ! write(String, "('Displacement_Load_Case',I2.2)") CURLCASE
  !write (VTKTmpFile) String, 3, NUMNP
  
  DO II=1,NUMNP
     IC=IC + 1
     IF (IC.GE.56) THEN
        WRITE (IOUT,"(//,' D I S P L A C E M E N T S',//,'  NODE ',3X,   &
                          'X-DISPLACEMENT   Y-DISPLACEMENT  Z-DISPLACEMENT  X-ROTATION  Y-ROTATION  Z-ROTATION')")
        IC=4
     END IF

     DO I=1,6
        D(I)=0.
     END DO

     DO I=1,6
        KK=ID(I,II)
        IF (KK.NE.0) D(I)=DISP(KK)
     END DO

     WRITE (IOUT,'(1X,I10,5X,6E14.6)') II,D
  !   write (VTKTmpFile) D(1:3)                                    !Displacements

  END DO
  
  !输出最终应力
   REWIND PRESENTSTRESS
   READ(PRESENTSTRESS) (ELESTRESS(I),I=1,NPAR(2))
   
   WRITE (IOUT,"(//,' ELASTIC TRIAL SOLUTION  F O R  ',  &
                                           'E L E M E N T  G R O U P ',I4,//,   &
                                           '  ELEMENT',12X,'STRESS',/,'  NUMBER')") NG
   DO I=1,NPAR(2)
     WRITE(IOUT,"(1X,I5,11X,E13.6)") I,ELESTRESS(I)  
   ENDDO
  
ENDSUBROUTINE FINALWRITE