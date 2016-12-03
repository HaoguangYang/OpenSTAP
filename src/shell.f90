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
    
SUBROUTINE SHELL
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the PLATE element subroutine         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N101, N102, N103, N104, N105, N106, N107

  NUME = NPAR(2)
  NUMMAT = NPAR(3)

! Allocate storage for element group data
! 此处材料要求每一种提供E, Possion
! 每个element需要
  IF (IND == 1) THEN
      MM = 2*NUMMAT*ITWO + 22*NUME + 12*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: POSSION(NUMMAT)
! N103: LM(20,NUME)
! N104: XYZ(12,NUME)
! N105: MTAP(NUME)
! N106: NLAST
  
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+20*NUME
  N105=N104+12*NUME*ITWO
  N106=N105+NUME
  N107=N106+NUME*ITWO
  NLAST=N107

  MIDEST=NLAST - NFIRST

  CALL SHELL4Q (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105), A(106))

  RETURN

END SUBROUTINE SHELL


SUBROUTINE SHELL4Q (ID,X,Y,Z,U,MHT,E,POSSION,LM,XYZ,MATP, THICK)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(5,NUMNP),LM(20,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),POSSION(NPAR(3)),  &
             XYZ(12,NPAR(2)),THICK(NPAR(2)),U(NEQ), DISP(20,1)=0

  INTEGER :: NPAR1, NUME, NUMMAT, ND, I, J, K, L, M, N
  INTEGER :: MTYPE, IPRINT

  REAL(8) :: Cb(3, 3),Cc(3, 3), Cs, Etemp, Ptemp, det
  REAL(8) :: GAUSS(2) = (/-0.5773502692,0.5773502692/)
  REAL(8) :: G1, G2, GN(2,4), Ja(2,2), Ja_inv(2,2), Bk(3,20),By(2,20),Bm(3,20), S(20,20), BB(2,4), NShape(1,4)
  REAL(8) :: X_Y(4, 2), XY_G(1,2), STR(3,1)
  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=20

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
                   ' NUMBER     MODULUS',10X,'AREA')")

     DO I=1,NUMMAT
        READ (IIN,'(I5,2F10.0)') N,E(N),POSSION(N)  ! Read material information
        WRITE (IOUT,"(I5,4X,E12.5,2X,E14.6)") N,E(N),POSSION(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE     NODE     NODE       MATERIAL',/,   &
                      ' NUMBER-N      I        J        K        L       SET NUMBER')")
     
     N=0
     LM = 0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(6I5, F10.0, I5)') N,I,J,K,L,MTYPE,THICK(N)  ! Read in element information

!       Save element information
        XYZ(1,N)=X(I)  
        XYZ(2,N)=Y(I)
        
        XYZ(4,N)=X(J) 
        XYZ(5,N)=Y(J)

        XYZ(7,N)=X(K) 
        XYZ(8,N)=Y(K)
        
        XYZ(10,N)=X(L)
        XYZ(11,N)=Y(L)

        MATP(N)=MTYPE  ! Material type

        DO M=1,5
           LM(M,N)=ID(M,I)     ! Connectivity matrix
           LM(M+5,N)=ID(M,J)
           LM(M+10,N)=ID(M,K)
           LM(M+15,N)=ID(M,L)
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,4X,I5,4X,I5,4X,I5,7X,I5)") N,I,J,K,L,MTYPE

     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN
    
     DO N=1,NUME
        MTYPE=MATP(N)
        Etemp = E(MTYPE)
        Ptemp = POSSION(MTYPE)
        DO L = 1,4
            X_Y(L,1) = XYZ(3*L-2, N)
            X_Y(L,2) = XYZ(3*L-1, N)
        END DO
! 计算D
        Cb(1,1) = 1
        Cb(1,2) = Ptemp
        Cb(1,3) = 0
        Cb(2,1) = Ptemp
        Cb(2,2) = 1
        Cb(2,3) = 0
        Cb(3,1) = 0
        Cb(3,2) = 0
        Cb(3,3) = (1-Ptemp)/2
        Cc = Cb*Etemp/12.0/(1-Ptemp*Ptemp)*5.0/6.0
        Cb = Cb*Etemp/12.0/(1-Ptemp*Ptemp)*5.0/6.0

        Cs = Etemp/(2*(1+Ptemp))
! Gauss 积分常数
        S = 0
        DO L=1,2
            DO M=1,2
                G1 = GAUSS(L)
                G2 = GAUSS(M)
! 计算Jacobian
                GN = reshape((/G2-1,G1-1, 1-G2,-G1-1, 1+G2,1+G1, -G2-1,1-G1/), shape(GN))/4
                Ja = matmul(GN,X_Y)
                det = Ja(1,1)*Ja(2,2) - Ja(1,2)*Ja(2,1)
                Ja_inv(1,1) = Ja(2,2)
                Ja_inv(2,1) = -Ja(2,1)
                Ja_inv(1,2) = -Ja(1,2)
                Ja_inv(2,2) = Ja(1,1)
                Ja_inv = Ja_inv/det
                BB = matmul(Ja_inv, GN)
! 为弯曲部分的Bk赋值，改成循环
                Bk = 0
                DO K = 1,4
                    Bk(1,5*K-3) = BB(1,K)
                    Bk(2,5*K-2)   = BB(2,K)
                    Bk(3,5*K-3) = BB(2,K)
                    Bk(3,5*K-2)   = BB(1,K)
                END DO
! 为剪切部分的By赋值。改成循环
                By = 0
                DO K = 1,4
                    By(1,5*K-4) = BB(1,K)
                    By(1,5*K-3) = -1
                    By(2,5*K-4) = BB(2,K)
                    By(2,5*K-2)   = -1
                END DO
! 为平面部分的Bm赋值。改成循环
                Bm = 0
                DO K = 1,4
                    Bm(1,5*K-1) = BB(1,K)
                    Bm(2,5*K) = BB(2,K)
                    Bm(3,5*K-1) = BB(2,K)
                    Bm(3,5*K)   = BB(1,K)
                END DO
! 这里不要忘了还要乘上z方向积分
                S = S + (matmul(matmul(transpose(Bk), Cb), Bk)/12.0 + 5.0/6.0*Cs*matmul(transpose(By), By) + matmul(matmul(transpose(Bm), Cc), Bm))*abs(det)
            END DO
        END DO

        CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)  ! 这里要输出的S就是制作好了的local stiffness matrix
        
     END DO

     RETURN

! Stress calculations
  ELSE IF (IND .EQ. 3) THEN
    
  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

END SUBROUTINE SHELL4Q