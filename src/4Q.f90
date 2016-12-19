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
!        4Q element                                                       -
!        Qi He,(2016)                                                     -
!        Tsinghua University                                              -
!                                                                         -
!                                                                         -
!--------------------------------------------------------------------------

SUBROUTINE ELEMENT_4Q
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the 4Q element subroutine            .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

USE GLOBALS
USE memAllocate

IMPLICIT NONE
INTEGER NUME, NUMMAT, MM, N(7)

NUME = NPAR(2)
NUMMAT = NPAR(3)
NPAR(5) = 4

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: POISSON(NUMMAT)
! N103: LM(8,NUME)
! N104: XYZ(12,NUME)
! N105: MTAP(NUME)
  N(1)=0
  N(2)=N(1)+NUMMAT*ITWO
  N(3)=N(2)+NUMMAT*ITWO
  N(4)=N(3)+8*NUME
  N(5)=N(4)+12*NUME*ITWO
  N(6)=N(5)+NUME
  N(7)=N(6)+NPAR(5)*NPAR(2)
  
  MIDEST=N(7)
  if (IND .EQ. 1) then
        ! Allocate storage for element group data
        call MemAlloc(11,"ELEGP",MIDEST,1)
  end if
  NFIRST = NP(11)   ! Pointer to the first entry in the element group data array in the unit of single precision (corresponding to A)
  N(:) = N(:) + NFIRST
  NLAST=N(7)

  CALL ELEMENT_4Q_MAIN (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(7)),   &
       A(N(1)),A(N(2)),A(N(3)),A(N(4)),A(N(5)),A(N(6)))

  RETURN

END SUBROUTINE ELEMENT_4Q


SUBROUTINE ELEMENT_4Q_MAIN (ID,X,Y,Z,U,MHT,E,POISSON,LM,XYZ,MATP,Node)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE
  USE MathKernel

  IMPLICIT NONE

  
  INTERFACE
    FUNCTION NmatElast2D(eta,psi)
        IMPLICIT NONE
        REAL(8):: eta
        REAL(8):: psi
        REAL(8):: NmatElast2D(2,8)
    END FUNCTION
    FUNCTION BmatElast2D(eta,psi,C)
        IMPLICIT NONE
        REAL(8):: eta,psi
        REAL(8):: C(4,2)
        REAL(8):: BmatElast2D(3,8)
    END FUNCTION
  END INTERFACE
  
  INTEGER :: ID(3,NUMNP),LM(8,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),&
             E(NPAR(3)),POISSON(NPAR(3)),  &
             XYZ(12,NPAR(2)),U(NEQ),DST(8,1)
  REAL(8) :: ETA,EPSILON

  INTEGER :: NPAR1, NUME, NUMMAT, ND, Node(NPAR(2),NPAR(5)), L, N, I, J
  INTEGER :: MTYPE, IPRINT
  REAL(8) :: GP(2),W(2),NMAT(2,8),BMAT(3,8),C(4,2),NA(1,4)
  REAL(8) :: KE(8,8),DETJ,D(3,3)
  REAL(8) :: X_GUASS(4,2),XY(1,2), GaussianCollection(2,NPAR(2)*4)
  REAL(8) :: STRESS_XX(NPAR(2),4),STRESS_YY(NPAR(2),4),STRESS_XY(NPAR(2),4),STRESS(3,1)
  COMMON DETJ
  
  !定义gauss积分常数
  CALL GaussianMask(GP, W, 2)
  !GP(1)=-sqrt(3D0)/3.           !-0.57735027
  !GP(2)=sqrt(3D0)/3.            !0.57735027
  !W(1)=1.0
  !W(2)=1.0

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=8

! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.2, 4Q ELEMENTS',/,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I5,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S        POISSON',/,  &
                   ' NUMBER     MODULUS',9X,'RATIO',/,  &
                   15 X,'E',14X,'A')")

     DO I=1,NUMMAT
        READ (IIN,'(I5,2F10.0)') N,E(N),POISSON(N)  ! Read material information
        WRITE (IOUT,"(I5,4X,E12.5,2X,E14.6)") N,E(N),POISSON(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT    |------------- NODE -------------|    MATERIAL',/,   &
                      ' NUMBER-N      1        2        3        4      SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(7I5)') N, Node(N,1:NPAR(5)), MTYPE  ! Read in element information

!       Save element information
        XYZ(1:NPAR(5)*3-1:3,N)=X(Node(N,:))  ! Coordinates of the element's nodes
        XYZ(2:NPAR(5)*3  :3,N)=Y(Node(N,:))
        XYZ(3:NPAR(5)*3+1:3,N)=Z(Node(N,:))

        MATP(N)=MTYPE  ! Material type

        DO L=1,8
           LM(L,N)=0
        END DO

        DO L=1,2
           LM(L  ,N) = ID(L,Node(N,1))     ! Connectivity matrix
           LM(L+2,N) = ID(L,Node(N,2))
           LM(L+4,N) = ID(L,Node(N,3))
           LM(L+6,N) = ID(L,Node(N,4))
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,4X,I5,4X,I5,4X,I5,7X,I5)") N,Node(N,1:NPAR(5)),MTYPE
        write (VTKNodeTmp) NPAR(5), Node(N,:)-1

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
                
                NMAT = NmatElast2D(ETA,EPSILON)
                BMAT = BmatElast2D(ETA,EPSILON,C)
                
                KE = KE + W(I)*W(J)*MATMUL(MATMUL(TRANSPOSE(BMAT),D),BMAT)*abs(DETJ)
                
            END DO
         END DO       

        if(pardisodoor) then
            call pardiso_addban(DA(NP(8)),IA(NP(5)),IA(NP(6)),KE,LM(1,N),ND)
        else
            CALL ADDBAN (DA(NP(8)),IA(NP(2)),KE,LM(1,N),ND)
        end if
     END DO
        !write (*,*) "--------------------K--------------------"
        !write (*,*) DA(NP(3): NP(3)+NumberOfMatrixElements-1)
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
                
                NMAT = NmatElast2D(ETA,EPSILON)
                
                NA(1,:) = (/NMAT(1,1) , NMAT(1,3) , NMAT(1,5) , NMAT(1,7)/)
                
                XY = MATMUL(NA,C)
                X_GUASS(2*I+J-2,1) = XY(1,1)
                X_GUASS(2*I+J-2,2) = XY(1,2)
                

                BMAT = BmatElast2D(ETA,EPSILON,C)
                STRESS = MATMUL(D,MATMUL(BMAT,DST))
                
                STRESS_XX(N,2*I+J-2) = STRESS(1,1)
                STRESS_YY(N,2*I+J-2) = STRESS(2,1)
                STRESS_XY(N,2*I+J-2) = STRESS(3,1)
                
            WRITE (IOUT,"(1X,I5,4X,E13.6,4X,E13.6,11X,E13.6,4X,E13.6,4X, E13.6)") N,X_GUASS(2*I+J-2,1), &
            X_GUASS(2*I+J-2,2),STRESS_XX(N,2*I+J-2),STRESS_YY(N,2*I+J-2),STRESS_XY(N,2*I+J-2)
                
            END DO
        END DO
        
        I = (N-1)*4+1
        J = N*4
        GaussianCollection (:,I:J) = transpose(X_GUASS)
     END DO
     !call PostProcessor(NPAR(1), 2, XYZ((/1,2,4,5,7,8,10,11/),:), Node, 4, GaussianCollection, &
     !    TRANSPOSE(reshape((/transpose(STRESS_XX),transpose(STRESS_YY),transpose(STRESS_XY)/),(/4*NPAR(2),3/))), U)
  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

END SUBROUTINE ELEMENT_4Q_MAIN
    
FUNCTION NmatElast2D(eta,psi)
  IMPLICIT NONE
  REAL(8):: eta
  REAL(8):: psi
  REAL(8):: NmatElast2D(2,8),N(2,8)
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
  
  NmatElast2D=N
  
END FUNCTION NmatElast2D
  
FUNCTION BmatElast2D(eta,psi,C)
    IMPLICIT NONE
    REAL(8):: eta,psi
    REAL(8):: C(4,2)
    REAL(8):: GN(2,4),J(2,2)
    REAL(8):: DETJ,INVJ(2,2)
    REAL(8):: BB(2,4),B1x,B2x,B3x,B4x,B1y,B2y,B3y,B4y
    REAL(8):: B(3,8),BmatElast2D(3,8)
    REAL(8):: ZERO
    COMMON DETJ
    
    GN(1,:)=0.25*(/eta-1,1-eta,1+eta,-eta-1/)
    GN(2,:)=0.25*(/psi-1,-psi-1,1+psi,1-psi/)
    
    J = MATMUL(GN,C)
    DETJ = J(1,1)*J(2,2)-J(1,2)*J(2,1)
    
    INVJ(1,:)=1/DETJ*(/J(2,2),-J(1,2)/)
    INVJ(2,:)=1/DETJ*(/-J(2,1),J(1,1)/)
    
    BB= MATMUL(INVJ,GN)
    
    B1x     = BB(1,1) 
    B2x     = BB(1,2) 
    B3x     = BB(1,3) 
    B4x     = BB(1,4)
    B1y     = BB(2,1) 
    B2y     = BB(2,2) 
    B3y     = BB(2,3) 
    B4y     = BB(2,4) 
    ZERO    = 0.0
    
    B(1,:) = (/B1x   ,   ZERO   ,  B2x  ,   ZERO   ,   B3x  ,  ZERO   ,   B4x   ,  ZERO  /)
    B(2,:) = (/  ZERO   ,  B1y  ,  ZERO  ,   B2y   ,   ZERO   ,  B3y   ,  ZERO    ,  B4y/)
    B(3,:) = (/B1y  ,   B1x  ,  B2y  ,  B2x  ,   B3y   , B3x  ,  B4y  ,   B4x/)
    
    BmatElast2D = B
END FUNCTION BmatElast2D

