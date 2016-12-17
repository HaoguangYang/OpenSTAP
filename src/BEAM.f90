SUBROUTINE BEAM
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the BEAMELE element subroutine
! .   CALCULATED THE BEAM ELEMENT STIFFNESS MATRIX
!.    POSTPROCESSING FOR THE BEAM ELEMENT
!     DONE BY LIUCHANGWU                     .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N(14)

  NUME = NPAR(2)
  NUMMAT = NPAR(3)
  NPAR(5) = 2

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: G(NUMMAT)
! N103: AREA(NUMMAT)
! N104: I_X(NUMMAT)
! N105: I_Y(NUMMAT)
! N106: I_Z(NUMMAT)
! N107: J_X(NUMMAT)
! N108: J_Y(NUMMAT)
! N109: J_Z(NUMMAT)
! N110: LM(12,NUME)
! N111: XYZ(6,NUME)
! N112: MTAP(NUME)

  N(1)=0
  N(2)=N(1)+NUMMAT*ITWO
  N(3)=N(2)+NUMMAT*ITWO
  N(4)=N(3)+NUMMAT*ITWO
  N(5)=N(4)+NUMMAT*ITWO
  N(6)=N(5)+NUMMAT*ITWO
  N(7)=N(6)+NUMMAT*ITWO
  N(8)=N(7)+NUMMAT*ITWO
  N(9)=N(8)+NUMMAT*ITWO
  N(10)=N(9)+NUMMAT*ITWO
  N(11)=N(10)+12*NUME
  N(12)=N(11)+6*NUME*ITWO
  N(13)=N(12)+NUME
  N(14)=N(13)+NPAR(5)*NPAR(2)
  
  MIDEST=N(14)
  if (IND .EQ. 1) then
        ! Allocate storage for element group data
        call MemAlloc(11,"ELEGP",MIDEST,1)
  end if
  NFIRST = NP(11)   ! Pointer to the first entry in the element group data array in the unit of single precision (corresponding to A)
  N(:) = N(:) + NFIRST
  NLAST=N(14)

  CALL BEAMELE (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N(1)),A(N(2)),A(N(3)),A(N(4)),A(N(5)),A(N(6)),A(N(7)),A(N(8)),A(N(9)),A(N(10)),A(N(11)),A(N(12)),A(N(13)))

  RETURN

END SUBROUTINE BEAM


SUBROUTINE BEAMELE (ID,X,Y,Z,U,MHT,E,G,AREA,I_X,I_Y,I_Z,J_X,J_Y,J_Z,LM,XYZ,MATP,Node)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   BEAM ELEMENT subroutine                                        .
! .   
! .   OUTPUT: S(12,12) THE ELEMENT STIFFNESS MATRIX  
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(6,NUMNP),LM(12,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),G(NPAR(3)),AREA(NPAR(3)),  &
             I_X(NPAR(3)),I_Y(NPAR(3)),I_Z(NPAR(3)),J_X(NPAR(3)),J_Y(NPAR(3)),J_Z(NPAR(3)),XYZ(6,NPAR(2)),U(NEQ)
  REAL(8) :: S(12,12),D(3),UELE(12),FORCE(12),SIGMA   !UELE: THE DISPLACEMENT OF ELEMENT N; FORCE: THE INTERNAL FORCE AND MOMENT FOR ELEMENT N

  INTEGER :: NPAR1, NUME, NUMMAT, ND, I, J, L, N, Node(NPAR(2),NPAR(5))
  INTEGER :: MTYPE, IPRINT
  REAL(8) :: XL2, XL ,CXX,CXY,CXZ,CYX,CYY,CYZ,CZX,CZY,CZZ,T(12,12)!TRANSLATION MATRIX
  REAL(8) :: StressCollection(6,NPAR(2)), GaussianCollection(3,NPAR(2))

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=12

! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.1, TRUSS ELEMENTS',/,      &
                   '     EQ.2, ELEMENTS CURRENTLY',/,  &
                   '     EQ.5, 3D BEAM',//,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I5,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S     SHEARING     AREA     I_X     I_Y     I_Z     J_X     J_Y     J_Z',/,  &
                   ' NUMBER     MODULUS',10X,'MODULUS',/,  &
                   15 X,'E',14X,'A')")

     DO I=1,NUMMAT
        READ (IIN,'(I5,9F10.0)') N,E(N),G(N),AREA(N),I_X(N),I_Y(N),I_Z(N),J_X(N),J_Y(N),J_Z(N)  ! Read material information
        WRITE (IOUT,"(I5,1X,E12.5,1X,E12.5,1X,E12.5,1X,E12.5,1X,E12.5,1X, E12.5,1X,E12.5, 1X,E12.5,1X,E12.5)") &
        N,E(N),G(N),AREA(N),I_X(N),I_Y(N),I_Z(N),J_X(N),J_Y(N),J_Z(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE       MATERIAL',/,   &
                      ' NUMBER-N      I        J       SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(5I5)') N,Node(N,1:NPAR(5)),MTYPE  ! Read in element information

!       Save element information
        XYZ(1:NPAR(5)*3-1:3,N)=X(Node(N,:))  ! Coordinates of the element's nodes
        XYZ(2:NPAR(5)*3  :3,N)=Y(Node(N,:))
        XYZ(3:NPAR(5)*3+1:3,N)=Z(Node(N,:))

        MATP(N)=MTYPE  ! Material type

        DO L=1,12
           LM(L,N)=0
        END DO

        DO L=1,6
           LM(L,N)=ID(L,Node(N,1))     ! Connectivity matrix
           LM(L+6,N)=ID(L,Node(N,2))
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,4X,I5,7X,I5)") N,Node(N,1:NPAR(5)),MTYPE
        write (VTKNodeTmp) NPAR(5), Node(N,:)-1

     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     DO N=1,NUME
        MTYPE=MATP(N)

        !CALCULATE THE LENGTH OF THE ELEMENT
        
        XL2=0.
        DO L=1,3
           D(L)=XYZ(L,N) - XYZ(L+3,N)
           XL2=XL2 + D(L)*D(L)
        END DO
        XL=SQRT(XL2)  
        
        ! COORDINATE TRANSLATION MATRIX T
        
        DO I=1,12
            DO J=1,12
                T(I,J)=0
            ENDDO
        ENDDO
        
        !ASSUME THE DIRECTION FOR PRINCIPLE AXIS        
        CXX=(XYZ(4,N)-XYZ(1,N))/XL
        CXY=(XYZ(5,N)-XYZ(2,N))/XL
        CXZ=(XYZ(6,N)-XYZ(3,N))/XL
        CYX=-CXX*CXZ/SQRT(CXX*CXX+CXY*CXY)
        CYY=-CXY*CXZ/SQRT(CXX*CXX+CXY*CXY)
        CYZ=SQRT(CXX*CXX+CXY*CXY)
        CZX=CXY/SQRT(CXX*CXX+CXY*CXY)
        CZY=-CXX/SQRT(CXX*CXX+CXY*CXY)
        CZZ=0.0
        
        T(1,1)=CXX
        T(1,2)=CXY
        T(1,3)=CXZ
        T(2,1)=CYX
        T(2,2)=CYY
        T(2,3)=CYZ
        T(3,1)=CZX
        T(3,2)=CZY
        T(3,3)=CZZ
        
        DO I=1,3
           DO J=1,3
              T(I+3,J+3)=T(I,J)
              T(I+6,J+6)=T(I,J)
              T(I+9,J+9)=T(I,J)
           ENDDO
        ENDDO
        
        !CALCULATE ELEMENT STIFFNESS MATRIX
        
        DO I=1,12
           DO J=1,12
              S(I,J)=0
           ENDDO
        ENDDO
        
        S(1,1)=E(MTYPE)*AREA(MTYPE)/XL
        S(2,2)=12.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL*XL)
        S(3,3)=12.0*E(MTYPE)*I_Y(MTYPE)/(XL*XL*XL)
        S(4,4)=G(MTYPE)*J_Z(MTYPE)/XL
        S(5,5)=4.0*E(MTYPE)*I_Y(MTYPE)/XL
        S(6,6)=4.0*E(MTYPE)*I_Z(MTYPE)/XL
        S(7,7)=E(MTYPE)*AREA(MTYPE)/XL
        S(8,8)=12.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL*XL)
        S(9,9)=12.0*E(MTYPE)*I_Y(MTYPE)/(XL*XL*XL)
        S(10,10)=G(MTYPE)*J_Z(MTYPE)/XL
        S(11,11)=4.0*E(MTYPE)*I_Y(MTYPE)/XL
        S(12,12)=4.0*E(MTYPE)*I_Z(MTYPE)/XL
        S(5,3)=-6.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL)
        S(6,2)=6.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL)
        S(7,1)=-E(MTYPE)*AREA(MTYPE)/XL
        S(8,2)=-12.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL*XL)
        S(8,6)=-6.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL)
        S(9,3)=-12.0*E(MTYPE)*I_Y(MTYPE)/(XL*XL*XL)
        S(9,5)=6.0*E(MTYPE)*I_Y(MTYPE)/(XL*XL)
        S(10,4)=-G(MTYPE)*J_X(MTYPE)/XL
        S(11,3)=-6.0*E(MTYPE)*I_Y(MTYPE)/(XL*XL)
        S(11,5)=2.0*E(MTYPE)*I_Y(MTYPE)/XL
        S(11,9)=6.0*E(MTYPE)*I_Y(MTYPE)/(XL*XL)
        S(12,2)=6.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL)
        S(12,6)=2.0*E(MTYPE)*I_Z(MTYPE)/XL
        S(12,8)=-6.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL)
        
        DO I=1,11
           DO J=I+1,12
               S(I,J)=S(J,I)
           ENDDO
        ENDDO

        S=MATMUL(TRANSPOSE(T),MATMUL(S,T))
        
        CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)

     END DO

     RETURN

! Stress calculations
  ELSE IF (IND .EQ. 3) THEN

     IPRINT=0
     DO N=1,NUME
        IPRINT=IPRINT + 1
        IF (IPRINT.GT.50) IPRINT=1
        IF (IPRINT.EQ.1) WRITE (IOUT,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  ',  &
              'E L E M E N T  G R O U P',I4,//,   &
              '  ELEMENT',2X,'FORCE_X1',3X,'FORCE_Y1',3X,'FORCE_Z1',3X,'MOMENT_X1',3X,'MOMENT_Y1',3X, 'MOMENT_Z1',3X, &
              'FORCE_X2',3X,'FORCE_Y2',3X,'FORCE_Z2',3X,'MOMENT_X2',3X,'MOMENT_Y2',3X,'MOMENT_Z2',/,'  NUMBER')") NG
        MTYPE=MATP(N)
        
        !CALCULATE THE LENGTH OF THE ELEMENT
        
        XL2=0.
        DO L=1,3
           D(L) = XYZ(L,N) - XYZ(L+3,N)
           XL2=XL2 + D(L)*D(L)
        END DO
        
         XL=SQRT(XL2)  
        
        ! COORDINATE TRANSLATION MATRIX T
        
        DO I=1,12
            DO J=1,12
                T(I,J)=0
            ENDDO
        ENDDO
        
        !ASSUME THE DIRECTION FOR PRINCIPLE AXIS
        
        CXX=(XYZ(4,N)-XYZ(1,N))/XL
        CXY=(XYZ(5,N)-XYZ(2,N))/XL
        CXZ=(XYZ(6,N)-XYZ(3,N))/XL
        CYX=-CXX*CXZ/SQRT(CXX*CXX+CXY*CXY)
        CYY=-CXY*CXZ/SQRT(CXX*CXX+CXY*CXY)
        CYZ=SQRT(CXX*CXX+CXY*CXY)
        CZX=CXY/SQRT(CXX*CXX+CXY*CXY)
        CZY=-CXX/SQRT(CXX*CXX+CXY*CXY)
        CZZ=0.0
        
        T(1,1)=CXX
        T(1,2)=CXY
        T(1,3)=CXZ
        T(2,1)=CYX
        T(2,2)=CYY
        T(2,3)=CYZ
        T(3,1)=CZX
        T(3,2)=CZY
        T(3,3)=CZZ
        
        DO I=1,3
           DO J=1,3
              T(I+3,J+3)=T(I,J)
              T(I+6,J+6)=T(I,J)
              T(I+9,J+9)=T(I,J)
           ENDDO
        ENDDO
        
        !CALCULATE ELEMENT STIFFNESS MATRIX
        
        DO I=1,12
           DO J=1,12
              S(I,J)=0
           ENDDO
        ENDDO
        
        S(1,1)=E(MTYPE)*AREA(MTYPE)/XL
        S(2,2)=12.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL*XL)
        S(3,3)=12.0*E(MTYPE)*I_Y(MTYPE)/(XL*XL*XL)
        S(4,4)=G(MTYPE)*J_Z(MTYPE)/XL
        S(5,5)=4.0*E(MTYPE)*I_Y(MTYPE)/XL
        S(6,6)=4.0*E(MTYPE)*I_Z(MTYPE)/XL
        S(7,7)=E(MTYPE)*AREA(MTYPE)/XL
        S(8,8)=12.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL*XL)
        S(9,9)=12.0*E(MTYPE)*I_Y(MTYPE)/(XL*XL*XL)
        S(10,10)=G(MTYPE)*J_Z(MTYPE)/XL
        S(11,11)=4.0*E(MTYPE)*I_Y(MTYPE)/XL
        S(12,12)=4.0*E(MTYPE)*I_Z(MTYPE)/XL
        S(5,3)=-6.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL)
        S(6,2)=6.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL)
        S(7,1)=-E(MTYPE)*AREA(MTYPE)/XL
        S(8,2)=-12.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL*XL)
        S(8,6)=-6.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL)
        S(9,3)=-12.0*E(MTYPE)*I_Y(MTYPE)/(XL*XL*XL)
        S(9,5)=6.0*E(MTYPE)*I_Y(MTYPE)/(XL*XL)
        S(10,4)=-G(MTYPE)*J_X(MTYPE)/XL
        S(11,3)=-6.0*E(MTYPE)*I_Y(MTYPE)/(XL*XL)
        S(11,5)=2.0*E(MTYPE)*I_Y(MTYPE)/XL
        S(11,9)=6.0*E(MTYPE)*I_Y(MTYPE)/(XL*XL)
        S(12,2)=6.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL)
        S(12,6)=2.0*E(MTYPE)*I_Z(MTYPE)/XL
        S(12,8)=-6.0*E(MTYPE)*I_Z(MTYPE)/(XL*XL)
        
        DO I=1,11
           DO J=I+1,12
               S(I,J)=S(J,I)
           ENDDO
        ENDDO

        S=MATMUL(TRANSPOSE(T),MATMUL(S,T))

     !GET THE DISPLACEMENT FOR ELEMENT N   
        DO L=1,12
           UELE(L)=0.0
        ENDDO
        DO L=1,12
           I=LM(L,N)
           IF (I .GT. 0)THEN
               UELE(L)=U(I)
           ENDIF
        ENDDO
        
      !CALCULATE FORCE AND MOMENT OF THE ELEMENT
         
        DO I=1,12
           SIGMA=0.0
           DO J=1,12
               SIGMA=SIGMA+UELE(J)*S(I,J)
           ENDDO
           FORCE(I)=SIGMA
        ENDDO
        
       !PRINT THE RESULTS FOR POSTPROCESSING
       
       WRITE (IOUT,"(1X,I5,1X,12(E12.6,1X))") N,(FORCE(I),I=1,12)
       GaussianCollection(:,N) = 0.5*(XYZ(4:6,N)+XYZ(1:3,N))
       StressCollection((/1,2,3,5,6,4/),N) = 0.5*(FORCE(7:12)-FORCE(1:6))
     END DO
     !call PostProcessor(NPAR(1), 3, XYZ, &
     !                   Node, 1, GaussianCollection, StressCollection, U)
  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

 END SUBROUTINE BEAMELE
    

    

