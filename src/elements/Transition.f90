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
!        Transition element                                                       -
!        Qi He,(2016)                                                     -
!        Tsinghua University                                              -
!                                                                         -
!                                                                         -
!--------------------------------------------------------------------------

SUBROUTINE ELEMENT_Transition
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the Transition element subroutine            .
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
      MM = 2*NUMMAT*ITWO + 11*NUME + 15*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: POISSON(NUMMAT)
! N103: LM(10,NUME)
! N104: XY(15,NUME)
! N105: MTAP(NUME)
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+10*NUME
  N105=N104+15*NUME*ITWO
  N106=N105+NUME
  NLAST=N106

  MIDEST=NLAST - NFIRST

  CALL ELEMENT_Transition_MAIN (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105))

  RETURN

END SUBROUTINE ELEMENT_Transition


SUBROUTINE ELEMENT_Transition_MAIN (ID,X,Y,Z,U,MHT,E,POISSON,LM,XYZ,MATP)
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
    FUNCTION NmatElastTr(eta,psi)
        IMPLICIT NONE
        REAL(8):: eta
        REAL(8):: psi
        REAL(8):: NmatElastTr(2,10)
    END FUNCTION
    FUNCTION BmatElastTr(eta,psi,C)
        IMPLICIT NONE
        REAL(8):: eta,psi
        REAL(8):: C(5,2)
        REAL(8):: BmatElastTr(3,10)
    END FUNCTION
  END INTERFACE
  
  INTEGER :: ID(6,NUMNP),LM(10,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),&
             E(NPAR(3)),POISSON(NPAR(3)),  &
             XYZ(15,NPAR(2)),U(NEQ),DST(10,1)
  REAL(8) :: ETA,EPSILON

  INTEGER :: NPAR1, NUME, NUMMAT, ND, Node(5), L, N, I, J
  INTEGER :: MTYPE, IPRINT
  INTEGER,PARAMETER :: GUASS_N=4
  REAL(8) :: NMAT(2,10),BMAT(3,10),C(5,2),NA(1,5)
  REAL(8) :: KE(10,10),DETJ,D(3,3)
  REAL(8) :: XY_GUASS(1,2)
  REAL(8) :: STRESS(3,1)
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
      W(1)=5.0/9.0
      W(2)=8.0/9.0
      W(3)=5.0/9.0
  ELSE IF (GUASS_N == 4) THEN
      GP(1) = -0.861136311594053
      GP(2) = -0.339981043584586
      GP(3) = 0.339981043584586
      GP(4) = 0.861136311594053
      W(1) = 0.347854845137454
      W(2) = 0.652145154837454
      W(3) = W(2)
      W(4) = W(1)
  END IF 

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=10

! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I10,/,   &
                   '     EQ.16, Transition ELEMENTS',/,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I10,/)") NPAR1,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I10,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S        POISSON      DENSITY',/,  &
                   ' NUMBER     MODULUS',9X,    'RATIO',/,  &
                   15 X,'E',14X,                  'v', 14X,    'ρ')")

     DO I=1,NUMMAT
        READ (IIN,'(I10,2F10.0)') N,E(N),POISSON(N)  ! Read material information
        WRITE (IOUT,"(I10,4X,E12.5,2X,E14.6)") N,E(N),POISSON(N)
     END DO
     
     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT    |------------- NODE -------------|    MATERIAL',/,   &
                      ' NUMBER-N      1        2        3        4      SET NUMBER')")     
     
     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(7I10)') N, Node(1:5), MTYPE  ! Read in element information

!       Save element information
        XYZ(1:13:3,N)=X(Node(:))  ! Coordinates of the element's nodes
        XYZ(2:14:3,N)=Y(Node(:))
        XYZ(3:15:3,N)=Z(Node(:))

        MATP(N)=MTYPE  ! Material type

        DO L=1,10
           LM(L,N)=0
        END DO

        DO L=1,2
           LM(L  ,N) = ID(L,Node(1))     ! Connectivity matrix
           LM(L+2,N) = ID(L,Node(2))
           LM(L+4,N) = ID(L,Node(3))
           LM(L+6,N) = ID(L,Node(4))
           LM(L+8,N) = ID(L,Node(5))
        END DO

!       Update column heights and bandwidth
        if (.NOT. PARDISODOOR) CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I10,6X,I10,4X,I10,4X,I10,4X,I10,4X,I10,7X,I10)") N,Node(1:5),MTYPE
!        write (VTKNodeTmp) NPAR(5), Node(N,:)-1

     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     DO N=1,NUME
        MTYPE=MATP(N)
              
        D(1,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/1D0,POISSON(MATP(N)),0D0/)
        D(2,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/POISSON(MATP(N)),1D0,0D0/)
        D(3,:)=E(MATP(N))/(1D0-POISSON(MATP(N))*POISSON(MATP(N)))*(/0D0,0D0,(1D0-POISSON(MATP(N)))/2D0/)
        
        DO I=1,5
            C(I,:)=(/XYZ(3*I-2,N),XYZ(3*I-1,N)/)
        END DO
        
        KE = 0
        DO I=1,GUASS_N
            DO J=1,GUASS_N
                ETA = GP(I)
                EPSILON = GP(J)
                BMAT = BmatElastTr(ETA,EPSILON,C)
                
                KE = KE + W(I)*W(J)*MATMUL(MATMUL(TRANSPOSE(BMAT),D),BMAT)*abs(DETJ)
                
            END DO
        END DO   
     

        !CALL ADDBAN (DA(NP(3)),IA(NP(2)),KE,LM(1,N),ND)
        if(pardisodoor) then
            call pardiso_addban(DA(NP(3)),IA(NP(2)),IA(NP(5)),KE,LM(1,N),ND)
        else
            CALL ADDBAN (DA(NP(3)),IA(NP(2)),KE,LM(1,N),ND)
        end if

     END DO
        !write (*,*) "--------------------K--------------------"
        !write (*,*) DA(NP(3): NP(3)+NumberOfMatrixElements-1)
     !RETURN

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
        
        DO I=1,5
            C(I,:)=(/XYZ(3*I-2,N),XYZ(3*I-1,N)/)
        END DO
        
        C(:,1) = (/0.0,1.0,1.0,0.0,1.0/)
        C(:,2) = (/0.0,0.0,1.0,1.0,0.5/)        
        
        DO I=1,10
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
                
                NMAT = NmatElastTr(ETA,EPSILON)
                
                NA(1,:) = (/NMAT(1,1) , NMAT(1,3) , NMAT(1,5) , NMAT(1,7) , NMAT(1,9)/)
                
                XY_GUASS = MATMUL(NA,C)

                BMAT = BmatElastTr(ETA,EPSILON,C)
                STRESS = MATMUL(D,MATMUL(BMAT,DST))
                
                
            WRITE (IOUT,"(1X,I10,4X,E13.6,4X,E13.6,11X,E13.6,4X,E13.6,4X, E13.6)") N, XY_GUASS(1,1), &
                   XY_GUASS(1,2),STRESS(1,1),STRESS(2,1),STRESS(3,1)
                
            END DO
        END DO
        
!        I = (N-1)*4+1
!        J = N*4
!        GaussianCollection (:,I:J) = transpose(X_GUASS)
     END DO
!     call PostProcessor(NPAR(1), 2, XYZ((/1,2,4,5,7,8,10,11/),:), Node, 4, GaussianCollection, &
!         TRANSPOSE(reshape((/transpose(STRESS_XX),transpose(STRESS_YY),transpose(STRESS_XY)/),(/4*NPAR(2),3/))), U)
  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

END SUBROUTINE ELEMENT_Transition_MAIN
    
FUNCTION NmatElastTr(ETA,PSI)
  IMPLICIT NONE
  REAL(8):: ETA
  REAL(8):: PSI
  REAL(8):: NmatElastTr(2,10),N(2,10)
  REAL(8):: N1,N2,N3,N4,N5
  
  N1 = 0.25*(1.0-ETA)*(1.0-PSI)
  N2 = -0.25*(1.0+ETA)*(1.0-PSI)*PSI
  N3 = 0.25*(1.0+ETA)*(1.0+PSI)*PSI
  N4 = 0.25*(1.0-ETA)*(1.0+PSI)
  N5 = (1.0-PSI**2)*0.5*(ETA+1.0)
  
  N(1,1)=N1
  N(1,2)=0.0
  N(1,3)=N2
  N(1,4)=0.0
  N(1,5)=N3
  N(1,6)=0.0
  N(1,7)=N4
  N(1,8)=0.0
  N(1,9)=N5
  N(1,10)=0.0  
  
  N(2,1)=0.0
  N(2,2)=N1
  N(2,3)=0.0
  N(2,4)=N2
  N(2,5)=0.0
  N(2,6)=N3
  N(2,7)=0.0
  N(2,8)=N4
  N(2,9)=0.0
  N(2,10)=N5
  
  NmatElastTr=N
  
END FUNCTION NmatElastTr
  
FUNCTION BmatElastTr(ETA,PSI,C)
    IMPLICIT NONE
    REAL(8):: ETA,PSI
    REAL(8):: C(5,2)
    REAL(8):: GN(2,5),J(2,2)
    REAL(8):: DETJ,INVJ(2,2)
    REAL(8):: BB(2,5),B1x,B2x,B3x,B4x,B1y,B2y,B3y,B4y,B5x,B5y
    REAL(8):: B(3,10),BmatElastTr(3,10)
    REAL(8):: ZERO
    COMMON DETJ
    
    GN(1,:)=(/-0.25*(1.0-PSI) , -0.25*(1.0-PSI)*PSI , 0.25*(1.0+PSI)*PSI , -0.25*(1.0+PSI) , (1.0-PSI**2)*0.5 /)
    GN(2,:)=(/-0.25*(1.0-ETA) , -0.25*(1.0+ETA)*(1.0-2*PSI) , 0.25*(1.0+ETA)*(1.0+2*PSI) , 0.25*(1.0-ETA) , (-2*PSI)*0.5*(ETA+1.0) /)
    
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
    B1y     = BB(2,1) 
    B2y     = BB(2,2) 
    B3y     = BB(2,3) 
    B4y     = BB(2,4)
    B5y     = BB(2,5)
    
    ZERO    = 0.0
    
    B(1,:) = (/B1x   ,   ZERO   ,  B2x  ,  ZERO   ,   B3x  ,  ZERO   ,   B4x   ,  ZERO  ,  B5x  ,  ZERO/)
    B(2,:) = (/ZERO   ,  B1y  ,  ZERO  ,   B2y   ,   ZERO   ,  B3y   ,  ZERO    ,  B4y  ,  ZERO  ,  B5y/)
    B(3,:) = (/B1y  ,   B1x  ,  B2y  ,  B2x  ,   B3y   , B3x  ,  B4y  ,   B4x  ,  B5y  ,  B5x/)
    
    BmatElastTr = B
END FUNCTION BmatElastTr

