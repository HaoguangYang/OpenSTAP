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
!        Revision: Haoguang Yang  (12.2016)                               -
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
INTEGER NUME, NUMMAT, MM, N(7)

NUME = NPAR(2)
NUMMAT = NPAR(3)
NPAR(5) = 3

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: POISSON(NUMMAT)
! N103: LM(6,NUME)
! N104: XYZ(9,NUME)
! N105: MTAP(NUME)
  N(1)=0
  N(2)=N(1)+NUMMAT*ITWO
  N(3)=N(2)+NUMMAT*ITWO
  N(4)=N(3)+6*NUME
  N(5)=N(4)+9*NUME*ITWO
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

  CALL ELEMENT_3T_MAIN (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N(1)),A(N(2)),A(N(3)),A(N(4)),A(N(5)),A(N(6)))

  RETURN

END SUBROUTINE ELEMENT_3T


SUBROUTINE ELEMENT_3T_MAIN (ID,X,Y,Z,U,MHT,E,POISSON,LM,XYZ,MATP,Node)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE

  
  INTERFACE
    FUNCTION NmatElast3T(EPSILON,ETA)
        IMPLICIT NONE
        REAL(8):: EPSILON,ETA
        REAL(8):: NmatElast3T(1,3)
    END FUNCTION
    FUNCTION BmatElast3T(X,Y)
        IMPLICIT NONE
        REAL(8):: X(3,1),Y(3,1)
        REAL(8):: BmatElast3T(3,6)
    END FUNCTION
  END INTERFACE
  
  INTEGER :: ID(3,NUMNP),LM(12,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),POISSON(NPAR(3)),  &
             XYZ(9,NPAR(2)),U(NEQ),DST(6,1)
  REAL(8) :: ETA,EPSILON

  INTEGER :: NPAR1, NUME, NUMMAT, ND, Node(NPAR(2),NPAR(5)), L, N, I, J
  INTEGER :: MTYPE, IPRINT
  INTEGER,PARAMETER:: GUASS_N=3
  REAL(8),ALLOCATABLE:: GP1(:),GP2(:),W(:)
  REAL(8) :: NMAT(1,3),BMAT(3,6),C(3,2)
  REAL(8) :: KE(6,6),DETJ,D(3,3),XY(3,2), StressCollection(3,NPAR(2)*3), GaussianCollection(2,NPAR(2)*3)
  REAL(8),ALLOCATABLE:: STRESS_XX(:,:),STRESS_YY(:,:),STRESS_XY(:,:),STRESS(:,:)
  COMMON DETJ
  
  !定义gauss积分常数
  
  ALLOCATE(GP1(GUASS_N),GP2(GUASS_N),W(GUASS_N))
  ALLOCATE(STRESS(3,GUASS_N))
  
  IF (GUASS_N == 3) THEN
      GP1(1)=1./6.              !0.16666666666
      GP1(2)=2./3.              !0.66666666666
      GP1(3)=1./6.              !0.16666666666
      GP2(1)=1./6.              !0.16666666666
      GP2(2)=1./6.              !0.16666666666
      GP2(3)=2./3.              !0.66666666666
      W(1)=1./6.                !0.16666666666
      W(2)=1./6.                !0.16666666666
      W(3)=1./6.                !0.16666666666
  ELSE 
      WRITE(*,*) 'YOU NEED TO CHANGE THE PROGRAM FILE IN 3T.F90'
  END IF 

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=6

! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.8, 3T ELEMENTS',//,     & 
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
                      ' ELEMENT    |-------- NODE ---------|      MATERIAL',/,   &
                      ' NUMBER-N      1        2        3       SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(5I5)') N, Node(N,1:NPAR(5)) ,MTYPE  ! Read in element information

!       Save element information
        XYZ(1:NPAR(5)*3-1:3,N)=X(Node(N,:))  ! Coordinates of the element's nodes
        XYZ(2:NPAR(5)*3  :3,N)=Y(Node(N,:))
        XYZ(3:NPAR(5)*3+1:3,N)=Z(Node(N,:))
        
        MATP(N)=MTYPE  ! Material type

        DO L=1,6
           LM(L,N)=0
        END DO

        DO L=1,2
           LM(L,N)=ID(L,Node(N,1))     ! Connectivity matrix
           LM(L+2,N)=ID(L,Node(N,2))
           LM(L+4,N)=ID(L,Node(N,3))
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,4X,I5,4X,I5,7X,I5)") N,Node(N,1:NPAR(5)),MTYPE
        write (VTKNodeTmp) NPAR(5), Node(N,:)-1

     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     DO N=1,NUME
        MTYPE=MATP(N)
              
        D(1,:)=E(MTYPE)/(1D0-POISSON(MTYPE)*POISSON(MTYPE))*(/1D0,POISSON(MTYPE),0D0/)
        D(2,:)=E(MTYPE)/(1D0-POISSON(MTYPE)*POISSON(MTYPE))*(/POISSON(MTYPE),1D0,0D0/)
        D(3,:)=E(MTYPE)/(1D0-POISSON(MTYPE)*POISSON(MTYPE))*(/0D0,0D0,(1D0-POISSON(MTYPE))/2D0/)
        
        C(1,:) = (/XYZ(1,N),XYZ(2,N)/)
        C(2,:) = (/XYZ(4,N),XYZ(5,N)/)
        C(3,:) = (/XYZ(7,N),XYZ(8,N)/)
        
        KE = 0
        
        BMAT = BmatElast3T(C(:,1),C(:,2))
        KE = 1.0/2*MATMUL(MATMUL(TRANSPOSE(BMAT),D),BMAT)*DETJ
        
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
                                           '  ELEMENT',4X,' X-CORRD',9X,'Y-CORRD',17X   &	
                                           'STRESS_XX',8X,'STRESS_YY',8X,'STRESS_XY')") NG
        MTYPE=MATP(N)
        
        D(1,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/1D0,POISSON(MATP(N)),0D0/)
        D(2,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/POISSON(MATP(N)),1D0,0D0/)
        D(3,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/0D0,0D0,(1D0-POISSON(MATP(N)))/2D0/)
        
        C(1,:) = (/XYZ(1,N),XYZ(2,N)/)
        C(2,:) = (/XYZ(4,N),XYZ(5,N)/)
        C(3,:) = (/XYZ(7,N),XYZ(8,N)/)
        
        DO I=1,6
            IF (LM(I,N)==0) THEN
                DST(I,1)=0
            ELSE
                DST(I,1)=U(LM(I,N))
            END IF
        END DO
        
        DO I=1,GUASS_N
            ETA = GP1(I)
            EPSILON = GP2(I)
            
            NMAT = NmatElast3T(ETA,EPSILON)
                
            XY(I,:) = reshape(matmul(NMAT,C),(/2/))
            
            BMAT = BmatElast3T(C(:,1),C(:,2))
            STRESS(:,I) = reshape(MATMUL(D,MATMUL(BMAT,DST)),(/3/))
                
            WRITE (IOUT,"(1X,I5,4X,E13.6,4X,E13.6,11X,E13.6,4X,E13.6,4X,E13.6)") &
                    &  N , XY(I,1) , XY(I,2) &
                    &  , STRESS(1,I) , STRESS(2,I) , STRESS(3,I)
                
            
        END DO
        I = (N-1)*3+1
        J = N*3
        GaussianCollection (:,I:J) = transpose(XY)
        StressCollection (:,I:J) = Stress
     END DO
     !call PostProcessor(NPAR(1), 2, XYZ((/1,2,4,5,7,8/),:), &
     !                      Node, 3, GaussianCollection, StressCollection, U)

  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF
  deallocate (GP1, GP2, W, STRESS)
END SUBROUTINE ELEMENT_3T_MAIN
    
FUNCTION NmatElast3T(ETA,EPSILON)
  IMPLICIT NONE
  REAL(8):: EPSILON,ETA
  REAL(8):: X1,X2,X3,Y1,Y2,Y3
  REAL(8):: NmatElast3T(1,3),N(1,3)
  REAL(8):: N1,N2,N3
  
  X1=1.0
  X2=0.0
  X3=0.0
  Y1=0.0
  Y2=1.0
  Y3=0.0
  
  N1 = (x2*y3-x3*y2)+(y2-y3)*EPSILON+(x3-x2)*ETA
  N2 = (x3*y1-x1*y3)+(y3-y1)*EPSILON+(x1-x3)*ETA
  N3 = (x1*y2-x2*y1)+(y1-y2)*EPSILON+(x2-x1)*ETA
  
  N(1,:) = (/N1,N2,N3/)
  NmatElast3T=N
  
END FUNCTION NmatElast3T
  
FUNCTION BmatElast3T(X,Y)
    IMPLICIT NONE
    REAL(8):: X(3,1),Y(3,1)
    REAL(8):: X1,X2,X3,Y1,Y2,Y3
    REAL(8):: J(2,2)
    REAL(8):: DETJ
    REAL(8):: B1x,B2x,B3x,B1y,B2y,B3y
    REAL(8):: B(3,6),BmatElast3T(3,6)
    REAL(8):: ZERO
    COMMON DETJ
    
    X1 = X(1,1)
    X2 = X(2,1)
    X3 = X(3,1)
    Y1 = Y(1,1)
    Y2 = Y(2,1)
    Y3 = Y(3,1)
    
    J(1,:) = (/X1-X3,Y1-Y3/)
    J(2,:) = (/X2-X3,Y2-Y3/)
    
    DETJ = J(1,1)*J(2,2)-J(1,2)*J(2,1)
    
    B1x     = y2-y3
    B2x     = y3-y1
    B3x     = y1-y2
    B1y     = x3-x2
    B2y     = x1-x3
    B3y     = x2-x1
    ZERO    = 0.0
    
    B(1,:) =1/DETJ* (/B1x   ,  ZERO ,  B2x   ,   ZERO  ,   B3x   ,  ZERO /)
    B(2,:) =1/DETJ* (/ZERO  ,  B1y  ,  ZERO  ,   B2y   ,   ZERO  ,  B3y  /)
    B(3,:) =1/DETJ* (/B1y   ,  B1x  ,  B2y   ,   B2x   ,   B3y   ,  B3x  /)
    
    BmatElast3T = B
END FUNCTION BmatElast3T
