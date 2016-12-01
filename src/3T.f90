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
!                                                                         -
!        3T element                                                       -
!        Qi He,(2016)                                                     -
!        Tsinghua University                                              -
!                                                                         -
!                                                                         -
!--------------------------------------------------------------------------

SUBROUTINE ELEMENT_3T
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the 3T element subroutine            .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

USE GLOBALS
USE memAllocate

IMPLICIT NONE
INTEGER NUME, NUMMAT, MM, N101, N102, N103, N104, N105, N106

NUME = NPAR(2)
NUMMAT = NPAR(3)

! Allocate storage for element group data
  IF (IND == 1) THEN
      MM = 2*NUMMAT*ITWO + 7*NUME + 9*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: POISSON(NUMMAT)
! N103: LM(6,NUME)
! N104: XYZ(9,NUME)
! N105: MTAP(NUME)
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+6*NUME
  N105=N104+9*NUME*ITWO
  N106=N105+NUME
  NLAST=N106

  MIDEST=NLAST - NFIRST

  CALL ELEMENT_4Q_MAIN (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105))

  RETURN

END SUBROUTINE ELEMENT_3T


SUBROUTINE ELEMENT_3T_MAIN (ID,X,Y,Z,U,MHT,E,POISSON,LM,XYZ,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE

  
  INTERFACE
    FUNCTION NmatElast3T(eta,psi)
        IMPLICIT NONE
        REAL(8):: eta
        REAL(8):: psi
        REAL(8):: NmatElast3T(2,8)
    END FUNCTION
    FUNCTION BmatElast3T(eta,psi,C)
        IMPLICIT NONE
        REAL(8):: eta,psi
        REAL(8):: C(4,2)
        REAL(8):: BmatElast3T(3,8)
    END FUNCTION
  END INTERFACE
  
  INTEGER :: ID(3,NUMNP),LM(8,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),POISSON(NPAR(3)),  &
             XYZ(12,NPAR(2)),U(NEQ),DST(8,1)
  REAL(8) :: ETA,EPSILON

  INTEGER :: NPAR1, NUME, NUMMAT, ND, I1, I2, I3, I4, L, N, I, J
  INTEGER :: MTYPE, IPRINT
  REAL(8) :: GP(2),W(2),NMAT(2,8),BMAT(3,8),C(4,2),NA(1,4)
  REAL(8) :: KE(8,8),DETJ,D(3,3)
  REAL(8) :: X_GUASS(4,2),XY(1,2)
  REAL(8) :: STRESS_XX(NPAR(2),4),STRESS_YY(NPAR(2),4),STRESS_XY(NPAR(2),4),STRESS(3,1)
  COMMON DETJ
  
  !定义gauss积分常数
  GP(1)=-0.57735027
  GP(2)=0.57735027
  W(1)=1.0
  W(2)=1.0

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=8

! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.1, TRUSS ELEMENTS',/,      &
                   '     EQ.2, ELEMENTS CURRENTLY',/,  &
                   '     EQ.3, NOT AVAILABLE',//,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I5,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S     CROSS-SECTIONAL',/,  &
                   ' NUMBER     MODULUS',10X,'AREA',/,  &
                   15 X,'E',14X,'A')")

     DO I=1,NUMMAT
        READ (IIN,'(I5,2F10.0)') N,E(N),POISSON(N)  ! Read material information
        WRITE (IOUT,"(I5,4X,E12.5,2X,E14.6)") N,E(N),POISSON(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE       MATERIAL',/,   &
                      ' NUMBER-N      I        J       SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(7I5)') N,I1,I2,I3,I4,MTYPE  ! Read in element information

!       Save element information
        XYZ(1,N)=X(I1)  ! Coordinates of the element's 1st node
        XYZ(2,N)=Y(I1)
        XYZ(3,N)=Z(I1)

        XYZ(4,N)=X(I2)  ! Coordinates of the element's 2nd node
        XYZ(5,N)=Y(I2)
        XYZ(6,N)=Z(I2)
        
        XYZ(7,N)=X(I3)  ! Coordinates of the element's 3rd node
        XYZ(8,N)=Y(I3)
        XYZ(9,N)=Z(I3)
        
        XYZ(10,N)=X(I4) ! Coordinates of the element's 4th node
        XYZ(11,N)=Y(I4)
        XYZ(12,N)=Z(I4)

        MATP(N)=MTYPE  ! Material type

        DO L=1,8
           LM(L,N)=0
        END DO

        DO L=1,2
           LM(L,N)=ID(L,I1)     ! Connectivity matrix
           LM(L+2,N)=ID(L,I2)
           LM(L+4,N)=ID(L,I3)
           LM(L+6,N)=ID(L,I4)
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,4X,I5,4X,I5,4X,I5,7X,I5)") N,I1,I2,I3,I4,MTYPE

     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     DO N=1,NUME
        MTYPE=MATP(N)
              
        D(1,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/1D0,POISSON(MATP(N)),0D0/)
        D(2,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/POISSON(MATP(N)),1D0,0D0/)
        D(3,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/0D0,0D0,(1D0-POISSON(MATP(N)))/2D0/)
        
        C(1,:) = (/XYZ(1,N),XYZ(2,N)/)
        C(2,:) = (/XYZ(4,N),XYZ(5,N)/)
        C(3,:) = (/XYZ(7,N),XYZ(8,N)/)
        C(4,:) = (/XYZ(10,N),XYZ(11,N)/)
        
        KE = 0

        DO I=1,2
            DO J=1,2
                ETA = GP(I)
                EPSILON = GP(J)
                
                NMAT = NmatElast3T(ETA,EPSILON)
                BMAT = BmatElast3T(ETA,EPSILON,C)
                
                KE = KE + W(I)*W(J)*MATMUL(MATMUL(TRANSPOSE(BMAT),D),BMAT)*DETJ
                
            END DO
         END DO       

        CALL ADDBAN (DA(NP(3)),IA(NP(2)),KE,LM(1,N),ND)

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
                                           '  ELEMENT',9X,' X-CORRD',9X,'Y-CORRD',9X   &	
                                           'STRESS_XX',7X,'STRESS_YY',9X,'STRESS_XY')") NG
        MTYPE=MATP(N)
        
        D(1,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/1D0,POISSON(MATP(N)),0D0/)
        D(2,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/POISSON(MATP(N)),1D0,0D0/)
        D(3,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/0D0,0D0,(1D0-POISSON(MATP(N)))/2D0/)
        
        C(1,:) = (/XYZ(1,N),XYZ(2,N)/)
        C(2,:) = (/XYZ(4,N),XYZ(5,N)/)
        C(3,:) = (/XYZ(7,N),XYZ(8,N)/)
        C(4,:) = (/XYZ(10,N),XYZ(11,N)/)
        
        DO I=1,8
            IF (LM(I,N)==0) THEN
                DST(I,1)=0
            ELSE
                DST(I,1)=U(LM(I,N))
            END IF
        END DO
        
        DO I=1,2
            DO J=1,2
                ETA = GP(I)
                EPSILON = GP(J)
                
                NMAT = NmatElast3T(ETA,EPSILON)
                
                NA(1,:) = (/NMAT(1,1) , NMAT(1,3) , NMAT(1,5) , NMAT(1,7)/)
                
                XY = MATMUL(NA,C)
                X_GUASS(2*J+I-2,1) = XY(1,1)
                X_GUASS(2*J+I-2,2) = XY(1,2)
                

                BMAT = BmatElast3T(ETA,EPSILON,C)
                STRESS = MATMUL(D,MATMUL(BMAT,DST))
                
                STRESS_XX(N,2*J+I-2) = STRESS(1,1)
                STRESS_YY(N,2*J+I-2) = STRESS(2,1)
                STRESS_XY(N,2*J+I-2) = STRESS(3,1)
                
            WRITE (IOUT,"(1X,I5,4X,E13.6,4X,E13.6,11X,E13.6,4X,E13.6,4X,E13.6)") N,X_GUASS(2*J+I-2,1),X_GUASS(2*J+I-2,2),STRESS_XX(N,2*J+I-2),STRESS_YY(N,2*J+I-2),STRESS_XY(N,2*J+I-2)
                
            END DO
        END DO
                
     END DO

  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

END SUBROUTINE ELEMENT_3T_MAIN
    
FUNCTION NmatElast3T(eta,psi)
  IMPLICIT NONE
  REAL(8):: eta
  REAL(8):: psi
  REAL(8):: NmatElast3T(2,8),N(2,8)
  REAL(8):: N1,N2,N3,N4
  
  N1 = 0.25*(1.0-psi)*(1.0-eta)
  N2 = 0.25*(1.0+psi)*(1.0-eta) 
  N3 = 0.25*(1.0+psi)*(1.0+eta) 
  N4 = 0.25*(1.0-psi)*(1.0+eta)
  
  N(1,1)=N1
  N(1,2)=0.0
  N(1,3)=N2
  N(1,4)=0.0
  N(1,5)=N3
  N(1,6)=0.0
  N(1,7)=N4
  N(1,8)=0.0
  N(2,1)=0.0
  N(2,2)=N1
  N(2,3)=0.0
  N(2,4)=N2
  N(2,5)=0.0
  N(2,6)=N3
  N(2,7)=0.0
  N(2,8)=N4
  
  NmatElast3T=N
  
END FUNCTION NmatElast3T
  
FUNCTION BmatElast3T(X,Y)
    IMPLICIT NONE
    REAL(8):: X(3),Y(3)
    REAL(8):: X1,X2,X3,Y1,Y2,Y3
    REAL(8):: J(2,2)
    REAL(8):: DETJ
!   REAL(8):: INVJ(2,2),BB(2,3)
    REAL(8):: B1x,B2x,B3x,B1y,B2y,B3y
    REAL(8):: B(3,6),BmatElast3T(3,6)
    REAL(8):: ZERO
    COMMON DETJ
    
    X1 = X(1)
    X2 = X(2)
    X3 = X(3)
    Y1 = Y(1)
    Y2 = Y(2)
    Y3 = Y(3)
    
    J(1,:) = (/X1-X3,Y1-Y3/)
    J(2,:) = (/X2-X3,Y2-Y3/)
    
    DETJ = J(1,1)*J(2,2)-J(1,2)*J(2,1)
    
!   INVJ(1,:)=1/DETJ*(/J(2,2),-J(1,2)/)
!   INVJ(2,:)=1/DETJ*(/-J(2,1),J(1,1)/)
    
!   BB= MATMUL(INVJ,GN)
    
    B1x     = y2-y3
    B2x     = y3-y1
    B3x     = y1-y2
    B1y     = x3-x2
    B2y     = x1-x3
    B3y     = x2-x1
    ZERO    = 0.0
    
    B(1,:) = (/B1x   ,   ZERO   ,  B2x  ,   ZERO   ,   B3x  ,  ZERO   /)
    B(2,:) = (/  ZERO   ,  B1y  ,  ZERO  ,   B2y   ,   ZERO   ,  B3y  /)
    B(3,:) = (/B1y  ,   B1x  ,  B2y  ,  B2x  ,   B3y   , B3x /)
    
    BmatElast3T = B
END FUNCTION BmatElast3T


  