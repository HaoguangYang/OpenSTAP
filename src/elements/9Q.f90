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
!        9Q element                                                       -
!        Qi He,(2016)                                                     -
!        Tsinghua University                                              -
!                                                                         -
!                                                                         -
!--------------------------------------------------------------------------

SUBROUTINE ELEMENT_9Q
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the 4Q element subroutine            .
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
      MM = 2*NUMMAT*ITWO + 19*NUME + 18*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: POISSON(NUMMAT)
! N103: LM(18,NUME)
! N104: XY(18,NUME)
! N105: MTAP(NUME)
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+18*NUME
  N105=N104+18*NUME*ITWO
  N106=N105+NUME
  NLAST=N106

  MIDEST=NLAST - NFIRST

  CALL ELEMENT_9Q_MAIN (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105))

  RETURN

END SUBROUTINE ELEMENT_9Q


SUBROUTINE ELEMENT_9Q_MAIN (ID,X,Y,Z,U,MHT,E,POISSON,LM,XY,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE

  
  INTERFACE
    FUNCTION NmatElast9Q(eta,psi)
        IMPLICIT NONE
        REAL(8):: eta
        REAL(8):: psi
        REAL(8):: NmatElast9Q(1,9)
    END FUNCTION
    FUNCTION BmatElast9Q(eta,psi,C)
        IMPLICIT NONE
        REAL(8):: eta,psi
        REAL(8):: C(9,2)
        REAL(8):: BmatElast9Q(3,18)
    END FUNCTION
  END INTERFACE
  
  INTEGER :: ID(6,NUMNP),LM(18,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP), &
             E(NPAR(3)),POISSON(NPAR(3)), XY(18,NPAR(2)),U(NEQ),DST(18,1)
  REAL(8) :: ETA,EPSILON

  INTEGER :: NPAR1, NUME, NUMMAT, ND, I1, I2, I3, I4, I5, I6, I7, I8, I9, L, N, I, J
  INTEGER :: MTYPE, IPRINT
  INTEGER,PARAMETER :: GUASS_N=3
  REAL(8) :: NMAT(1,9),BMAT(3,18),C(9,2),NA(1,9)
  REAL(8) :: KE(18,18),DETJ,D(3,3),M(18,18), Rho, Density(NPAR(3))
  REAL(8) :: X_GUASS(4,2),XY0(1,2)
  REAL(8) :: STRESS_XX(NPAR(2),4),STRESS_YY(NPAR(2),4),STRESS_XY(NPAR(2),4),STRESS(3,1)
  REAL(8),ALLOCATABLE:: GP(:),W(:)
  COMMON DETJ
  
  !定义gauss积分常数
  
  ALLOCATE(GP(GUASS_N),W(GUASS_N))
  
  IF (GUASS_N == 2) THEN
      GP(1)=-0.5773502692
      GP(2)=0.5773502692
      W(1)=1.0
      W(2)=1.0
  ELSE IF (GUASS_N == 3) THEN
      GP(1)=-0.7745966692
      GP(2)=0.0
      GP(3)=0.7745966692
      W(1)=5.0/9
      W(2)=8.0/9
      W(3)=5.0/9
  END IF 


  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=18

! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I10,/,   &
                   '     EQ.3, 9Q ELEMENTS',//,     &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I10,/)") NPAR1,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I10,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S        POISSON',/,  &
                   ' NUMBER     MODULUS',9X,'RATIO',/,  &
                   15 X,'E',14X,'A')")

     DO I=1,NUMMAT
        READ (IIN,'(I10,2F10.0)') N,E(N),POISSON(N)  ! Read material information
        WRITE (IOUT,"(I10,4X,E12.5,2X,E14.6)") N,E(N),POISSON(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE       MATERIAL',/,   &
                      ' NUMBER-N      I        J       SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(11I10)') N,I1,I2,I3,I4,I5,I6,I7,I8,I9,MTYPE  ! Read in element information

!       Save element information
        XY(1,N)=X(I1)  ! Coordinates of the element's 1st node
        XY(2,N)=Y(I1)

        XY(3,N)=X(I2)  ! Coordinates of the element's 2nd node
        XY(4,N)=Y(I2)
        
        XY(5,N)=X(I3)  ! Coordinates of the element's 3th node
        XY(6,N)=Y(I3)
        
        XY(7,N)=X(I4)  ! Coordinates of the element's 4th node
        XY(8,N)=Y(I4)
        
        XY(9,N)=X(I5)  ! Coordinates of the element's 5th node
        XY(10,N)=Y(I5)
        
        XY(11,N)=X(I6)  ! Coordinates of the element's 6th node
        XY(12,N)=Y(I6)
        
        XY(13,N)=X(I7)  ! Coordinates of the element's 7th node
        XY(14,N)=Y(I7)
        
        XY(15,N)=X(I8)  ! Coordinates of the element's 8th node
        XY(16,N)=Y(I8)
        
        XY(17,N)=X(I9)  ! Coordinates of the element's 9th node
        XY(18,N)=Y(I9)

        MATP(N)=MTYPE  ! Material type

        DO L=1,8
           LM(L,N)=0
        END DO

        DO L=1,2
           LM(L,N)=ID(L,I1)     ! Connectivity matrix
           LM(L+2,N)=ID(L,I2)
           LM(L+4,N)=ID(L,I3)
           LM(L+6,N)=ID(L,I4)
           LM(L+8,N)=ID(L,I5)
           LM(L+10,N)=ID(L,I6)
           LM(L+12,N)=ID(L,I7)
           LM(L+14,N)=ID(L,I8)
           LM(L+16,N)=ID(L,I9)
        END DO

!       Update column heights and bandwidth
        if (.NOT. PARDISODOOR) CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I10,6X,I10,4X,I10,4X,I10,4X,I10,4X,I10,4X,I10,4X,I10,4X,I10,4X,I10,7X,I10)") N,I1,I2,I3,I4,I5,I6,I7,I8,I9,MTYPE

     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     DO N=1,NUME
        MTYPE=MATP(N)
              
        D(1,:)=E(MTYPE)/(1D0-POISSON(MTYPE)*POISSON(MTYPE))*(/1D0,POISSON(MTYPE),0D0/)
        D(2,:)=E(MTYPE)/(1D0-POISSON(MTYPE)*POISSON(MTYPE))*(/POISSON(MTYPE),1D0,0D0/)
        D(3,:)=E(MTYPE)/(1D0-POISSON(MTYPE)*POISSON(MTYPE))*(/0D0,0D0,(1D0-POISSON(MTYPE))/2D0/)
        
        DO I=1,9
            C(I,:)=(/XY(2*I-1,N),XY(2*I,N)/)
        END DO
        
        KE = 0
        DO I=1,GUASS_N
            DO J=1,GUASS_N
                ETA = GP(I)
                EPSILON = GP(J)
                
                BMAT = BmatElast9Q(ETA,EPSILON,C)
                
                KE = KE + W(I)*W(J)*MATMUL(MATMUL(TRANSPOSE(BMAT),D),BMAT)*abs(DETJ)
                
            END DO
         END DO       

        if(pardisodoor) then
            call pardiso_addban(DA(NP(3)),IA(NP(2)),IA(NP(5)),KE,LM(1,N),ND)
        else
            CALL ADDBAN (DA(NP(3)),IA(NP(2)),KE,LM(1,N),ND)
        end if
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
                                           '  ELEMENT',4X,' X-CORRD',9X,'Y-CORRD',17X   &	
                                           'STRESS_XX',8X,'STRESS_YY',8X,'STRESS_XY')") NG
        MTYPE=MATP(N)
        
        D(1,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/1D0,POISSON(MATP(N)),0D0/)
        D(2,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/POISSON(MATP(N)),1D0,0D0/)
        D(3,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/0D0,0D0,(1D0-POISSON(MATP(N)))/2D0/)
        
        DO I=1,9
            C(I,:)=(/XY(2*I-1,N),XY(2*I,N)/)
        END DO
        
        DO I=1,18
            IF (LM(I,N)==0) THEN
                DST(I,1)=0
            ELSE
                DST(I,1)=U(LM(I,N))
            END IF
        END DO
        
        DO I=1,GUASS_N
            DO J=1,GUASS_N
                ETA = GP(I)
                EPSILON = GP(J)
                
                NMAT = NmatElast9Q(ETA,EPSILON)
                
                XY0 = MATMUL(NMAT,C)
                BMAT = BmatElast9Q(ETA,EPSILON,C)
                STRESS = MATMUL(D,MATMUL(BMAT,DST))
                
                
            WRITE (IOUT,"(1X,I10,4X,E13.6,4X,E13.6,11X,E13.6,4X,E13.6,4X, E13.6)") N, XY0(1,1), &
                   XY0(1,2),STRESS(1,1),STRESS(2,1),STRESS(3,1)
                
            END DO
        END DO
                
     END DO

  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF
  deallocate(GP,W)
END SUBROUTINE ELEMENT_9Q_MAIN
    
FUNCTION NmatElast9Q(eta,psi)
  IMPLICIT NONE
  REAL(8):: eta
  REAL(8):: psi
  REAL(8):: NmatElast9Q(1,9),NMAT(1,9)
  REAL(8):: N00(3),N11(3)
  REAL(8):: N1,N2,N3,N4,N5,N6,N7,N8,N9
  REAL(8):: ZERO=0.0
  INTEGER:: I
  
  DO I=1,3
      N00(1)=0.5*eta*(eta-1.0)
      N00(2)=1-eta*eta
      N00(3)=0.5*eta*(eta+1.0)
  END DO
  
  DO I=1,3
      N11(1)=0.5*psi*(psi-1.0)
      N11(2)=1-psi*psi
      N11(3)=0.5*psi*(psi+1.0)
  END DO
  
  N1 = N00(1)*N11(1)
  N2 = N00(2)*N11(1)
  N3 = N00(3)*N11(1)
  N4 = N00(3)*N11(2)
  N5 = N00(3)*N11(3)
  N6 = N00(2)*N11(3)
  N7 = N00(1)*N11(3)
  N8 = N00(1)*N11(2)
  N9 = N00(2)*N11(2)
  
  NMAT(1,:) = (/  N1,N2,N3,N4,N5,N6,N7,N8,N9  /)

!  N(1,1)=N1
!  N(1,2)=ZERO
!  N(1,3)=N2
!  N(1,4)=ZERO
!  N(1,5)=N3
!  N(1,6)=ZERO
!  N(1,7)=N4
!  N(1,8)=ZERO
!  N(1,9)=N5
!  N(1,10)=ZERO
!  N(1,11)=N6
!  N(1,12)=ZERO
!  N(1,13)=N7
!  N(1,14)=ZERO
!  N(1,15)=N8
!  N(1,16)=ZERO
!  N(1,17)=N9
!  N(1,18)=ZERO
  
!  N(2,1)=ZERO
!  N(2,2)=N1
!  N(2,3)=ZERO
!  N(2,4)=N2
!  N(2,5)=ZERO
!  N(2,6)=N3
!  N(2,7)=ZERO
!  N(2,8)=N4
!  N(2,9)=ZERO
!  N(2,10)=N5
!  N(2,11)=ZERO
!  N(2,12)=N6
!  N(2,13)=ZERO
!  N(2,14)=N7
!  N(2,15)=ZERO
!  N(2,16)=N8
!  N(2,17)=ZERO
!  N(2,18)=N9
  
  NmatElast9Q=NMAT
  
END FUNCTION NmatElast9Q
  
FUNCTION BmatElast9Q(eta,psi,C)
    IMPLICIT NONE
    REAL(8):: eta,psi
    REAL(8):: C(9,2) ! 坐标
    REAL(8):: GN(2,9),J(2,2)
    REAL(8):: DETJ,INVJ(2,2)
    REAL(8):: BB(2,9)
    REAL(8):: B1x,B2x,B3x,B4x,B5x,B6x,B7x,B8x,B9x,B1y,B2y,B3y,B4y,B5y,B6y,B7y,B8y,B9y
    REAL(8):: B(3,18),BmatElast9Q(3,18)
    REAL(8),PARAMETER:: ZERO=0.0
    COMMON DETJ
    
    GN(1,:)=(/(psi-0.5)*0.5*eta*(eta-1.0),(1-eta*eta)*(psi-0.5),0.5*eta*(eta+1.0)*(psi-0.5), &
             0.5*eta*(eta+1.0)*(-2.0*psi),0.5*eta*(eta+1.0)*(psi+0.5),(1.0-eta*eta)*(psi+0.5),0.5*eta*(eta-1.0)*(psi+0.5), &
             0.5*eta*(eta-1.0)*(-2.0*psi),(1.0-eta*eta)*(-2.0*psi)/)
    GN(2,:)=(/(eta-0.5)*0.5*psi*(psi-1.0),(-2.0*eta)*0.5*psi*(psi-1.0),(eta+0.5)*0.5*psi*(psi-1.0), &
              (eta+0.5)*(1.0-psi*psi),(eta+0.5)*0.5*psi*(psi+1.0),(-2.0*eta)*0.5*psi*(psi+1.0),(eta-0.5)*0.5*psi*(psi+1.0), &
              (eta-0.5)*(1.0-psi*psi),(-2.0*eta)*(1.0-psi*psi)/)
    
    J = MATMUL(GN,C)
    DETJ = J(1,1)*J(2,2)-J(1,2)*J(2,1)
    
    INVJ(1,:)=1/DETJ*(/J(2,2),-J(1,2)/)
    INVJ(2,:)=1/DETJ*(/-J(2,1),J(1,1)/)
    
    BB= MATMUL(INVJ,GN)
    
    B1x     = BB(1,1) 
    B2x     = BB(1,2) 
    B3x     = BB(1,3) 
    B4x     = BB(1,4)
    B5x     = BB(1,5) 
    B6x     = BB(1,6)
    B7x     = BB(1,7) 
    B8x     = BB(1,8)
    B9x     = BB(1,9)
    
    B1y     = BB(2,1) 
    B2y     = BB(2,2) 
    B3y     = BB(2,3) 
    B4y     = BB(2,4)
    B5y     = BB(2,5) 
    B6y     = BB(2,6)
    B7y     = BB(2,7) 
    B8y     = BB(2,8)
    B9y     = BB(2,9)
        
    B(1,:) = (/B1x   ,   ZERO   ,  B2x  ,   ZERO   ,   B3x  ,  ZERO   ,   B4x   ,  ZERO   ,   B5x   ,&
               ZERO   ,   B6x   ,  ZERO   ,   B7x   ,  ZERO   ,   B8x   ,  ZERO   ,   B9x   ,  ZERO  /)
    B(2,:) = (/  ZERO   ,  B1y  ,  ZERO  ,   B2y   ,   ZERO   ,  B3y   ,  ZERO    ,  B4y   ,  ZERO  ,&
               B5y   ,  ZERO    ,  B6y   ,  ZERO    ,  B7y   ,  ZERO    ,  B8y   ,  ZERO    ,  B9y/)
    B(3,:) = (/B1y  ,   B1x  ,  B2y  ,  B2x  ,   B3y   , B3x  ,  B4y  ,   B4x  ,  B5y  ,   B5x  ,&
               B6y  ,   B6x  ,  B7y  ,   B7x  ,  B8y  ,   B8x  ,  B9y  ,   B9x/)
    
    BmatElast9Q = B
END FUNCTION BmatElast9Q


  
