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
    
SUBROUTINE PLATE8Q
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the PLATE element subroutine         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N(8)

  NUME = NPAR(2)
  NUMMAT = NPAR(3)
  NPAR(5) = 8

! Allocate storage for element group data
! 此处材料要求每一种提供E, Possion
! 每个element需要

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: POSSION(NUMMAT)
! N103: LM(24,NUME)
! N104: XYZ(24,NUME)
! N105: MTAP(NUME)
! N106: THICK(NUME)
! N107: NLAST

  N(1)=0
  N(2)=N(1)+NUMMAT*ITWO
  N(3)=N(2)+NUMMAT*ITWO
  N(4)=N(3)+24*NUME
  N(5)=N(4)+24*NUME*ITWO
  N(6)=N(5)+NUME
  N(7)=N(6)+NUME*ITWO
  N(8)=N(7)+NPAR(5)*NPAR(2)
  
  MIDEST=N(8)
  if (IND .EQ. 1) then
        ! Allocate storage for element group data
        call MemAlloc(11,"ELEGP",MIDEST,1)
  end if
  NFIRST = NP(11)   ! Pointer to the first entry in the element group data array in the unit of single precision (corresponding to A)
  N(:) = N(:) + NFIRST
  NLAST=N(8)

  CALL PLATE8 (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N(1)),A(N(2)),A(N(3)),A(N(4)),A(N(5)),A(N(6)),A(N(7)))

  RETURN

END SUBROUTINE PLATE8Q


SUBROUTINE PLATE8 (ID,X,Y,Z,U,MHT,E,POSSION,LM,XYZ,MATP,THICK, Node)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE
  USE MathKernel

  IMPLICIT NONE
  INTEGER :: ID(6,NUMNP),LM(24,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),POSSION(NPAR(3)),  &
             XYZ(24,NPAR(2)),THICK(NPAR(2)),U(NEQ)

  REAL(8) :: DE(24,1)
  INTEGER :: NPAR1, NUME, NUMMAT, ND, I, J, K, L, I1,J1,K1,L1, M, N
  INTEGER :: MTYPE, IPRINT, Node(NPAR(2),NPAR(5))

  REAL(8) :: Cb(3, 3), Cs, Etemp, Ptemp, detJ, StressCollection(6, NPAR(2)*9), GaussianCollection(3, NPAR(2)*9)
  REAL(8) :: GAUSS(3), GAUSS_COF(3)
  REAL(8) :: G1, G2, GN(2,4), GN8(2,8), Ja(2,2), Ja_inv(2,2), Bk(3,24),By(2,24), S(24,24), BB(2,8)
  REAL(8) :: X_Y(4, 2), STR1(3,1), STR2(2,1), NN(1,8)
  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=24
  call GaussianMask(GAUSS, GAUSS_COF, 3)
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
                      ' ELEMENT   |------------------- NODE --------------------|      MATERIAL',/,   &
                      ' NUMBER-N     1    2     3     4      5     6     7     8     SET NUMBER')")
     
     N=0
     LM = 0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(10I5, F10.0, I5)') N,Node(N,1:NPAR(5)),MTYPE,THICK(N)  ! Read in element information

!       Save element information
        XYZ(1:NPAR(5)*3-1:3,N)=X(Node(N,:))     ! Coordinates of the element's nodes
        XYZ(2:NPAR(5)*3  :3,N)=Y(Node(N,:))
        XYZ(3:NPAR(5)*3+1:3,N)=Z(Node(N,:))
        MATP(N)=MTYPE  ! Material type

        DO M=1,3
           LM(M:M+22:3,N)=ID(M+2,Node(N,:))     ! Connectivity matrix
           !LM(M,N)=ID(M+2,I)
           !LM(M+3,N)=ID(M+2,J)
           !LM(M+6,N)=ID(M+2,K)
           !LM(M+9,N)=ID(M+2,L)
           !LM(M+12,N)=ID(M+2,I1)     
           !LM(M+15,N)=ID(M+2,J1)
           !LM(M+18,N)=ID(M+2,K1)
           !LM(M+21,N)=ID(M+2,L1)
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,8(I5,1X),4X,I5)") N,Node(N,1:NPAR(5)),MTYPE
        write (VTKNodeTmp) NPAR(5), Node(N,:)-1

     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN
    
     DO N=1,NUME
        MTYPE=MATP(N)
        Etemp = E(MTYPE)
        Ptemp = POSSION(MTYPE)
        DO L = 1,4
            X_Y(L,1)  = XYZ(3*L-2, N)
            X_Y(L,2)  = XYZ(3*L-1, N)
            !DE(3*L-2) = U
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
        Cb = Cb*Etemp/12.0/(1-Ptemp*Ptemp)*THICK(N)*THICK(N)*THICK(N)

        Cs = Etemp/(2*(1+Ptemp))*THICK(N)*5/6
! Gauss 积分常数
        S = 0
        DO L=1,3
            DO M=1,3
                G1 = GAUSS(L)
                G2 = GAUSS(M)
! 计算Jacobian
                GN = reshape((/G2-1,G1-1, 1-G2,-G1-1, 1+G2,1+G1, -G2-1,1-G1/), shape(GN))/4
                Ja = matmul(GN,X_Y)
                detJ = Det(Ja,2)
                Ja_inv(1,1) = Ja(2,2)
                Ja_inv(2,1) = -Ja(2,1)
                Ja_inv(1,2) = -Ja(1,2)
                Ja_inv(2,2) = Ja(1,1)
                Ja_inv = Ja_inv/detJ
                
            NN(1,1)=(1-G1)*(1-G2)*(-G1-G2-1)/4
            NN(1,2)=(1+G1)*(1-G2)*(G1-G2-1)/4
            NN(1,3)=(1+G1)*(1+G2)*(G1+G2-1)/4
            NN(1,4)=(1-G1)*(1+G2)*(-G1+G2-1)/4
            NN(1,5)=(1-G1*G1)*(1-G2)/2
            NN(1,6)=(1-G2*G2)*(1+G1)/2
            NN(1,7)=(1-G1*G1)*(1+G2)/2
            NN(1,8)=(1-G2*G2)*(1-G1)/2
! 因为可能写不成一行了，所以直接依次赋值了~
                GN8(1,1) = (1-G2)*(2*G1+G2)
                GN8(2,1) = (1-G1)*(G1+2*G2)
                GN8(1,2) = (1-G2)*(2*G1-G2)
                GN8(2,2) = (1+G1)*(2*G2-G1)
                GN8(1,3) = (1+G2)*(2*G1+G2)
                GN8(2,3) = (1+G1)*(2*G2+G1)
                GN8(1,4) = (1+G2)*(2*G1-G2)
                GN8(2,4) = (1-G1)*(2*G2-G1)
                GN8 = GN8/4
                GN8(1,5) = G1*(G2-1)
                GN8(2,5) = -(1-G1*G1)/2
                GN8(1,6) = (1-G2*G2)/2
                GN8(2,6) = -G2*(1+G1)
                GN8(1,7) = -G1*(1+G2)
                GN8(2,7) = (1-G1*G1)/2
                GN8(1,8) = (G2*G2-1)/2
                GN8(2,8) = G2*(G1-1)
                BB = matmul(Ja_inv, GN8)
! 为弯曲部分的Bk赋值，改成循环
                Bk = 0
                DO K = 1,8
                    Bk(1,3*K-1) = BB(1,K)
                    Bk(2,3*K)   = BB(2,K)
                    Bk(3,3*K-1) = BB(2,K)
                    Bk(3,3*K)   = BB(1,K)
                END DO
! 为剪切部分的By赋值。改成循环
                By = 0
                DO K = 1,8
                    By(1,3*K-2) = BB(1,K)
                    By(1,3*K-1) = -NN(1,K)
                    By(2,3*K-2) = BB(2,K)
                    By(2,3*K)   = -NN(1,K)
                END DO
 ! 这里不要忘了还要乘上z方向积分
                S = S + (matmul(matmul(transpose(Bk), Cb), Bk) + Cs*matmul(transpose(By), By)) &
                *abs(detJ)*GAUSS_COF(L)*GAUSS_COF(M)
            END DO
        END DO
        
        if(pardisodoor) then
            call pardiso_addban(DA(NP(3)),IA(NP(2)),IA(NP(5)),S,LM(1,N),ND)
        else
            CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)
        end if
        
     END DO

     RETURN

! Stress calculations
  ELSE IF (IND .EQ. 3) THEN
     WRITE (IOUT,"(//,' S T R E S S   I N F O R M A T I O N',//,  &
                  '           TAU_xx        TAU_yy        TAU_xy         TAU_yz       TAU_zx')")
     DO N=1,NUME
        WRITE (IOUT,"('ELEMENT', I3)") N
        MTYPE=MATP(N)
        Etemp = E(MTYPE)
        Ptemp = POSSION(MTYPE)
        DO L = 1,4
            X_Y(L,1) = XYZ(3*L-2, N)
            X_Y(L,2) = XYZ(3*L-1, N)
        END DO
        DO L = 1,ND
            IF(LM(L,N) == 0) THEN
                DE(L,1) = 0
            ELSE
                DE(L,1) = U(LM(L,N))
            ENDIF
        ENDDO
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
        Cb = Cb*Etemp/12.0/(1-Ptemp*Ptemp)*THICK(N)*THICK(N)*THICK(N)

        Cs = Etemp/(2*(1+Ptemp))*THICK(N)*5/6
! Gauss 积分常数
        S = 0
        DO L=1,3
            DO M=1,3
                G1 = GAUSS(L)
                G2 = GAUSS(M)
! 计算Jacobian
                GN = reshape((/G2-1,G1-1, 1-G2,-G1-1, 1+G2,1+G1, -G2-1,1-G1/), shape(GN))/4
                Ja = matmul(GN,X_Y)
                detJ = Det(Ja,2)
                Ja_inv(1,1) = Ja(2,2)
                Ja_inv(2,1) = -Ja(2,1)
                Ja_inv(1,2) = -Ja(1,2)
                Ja_inv(2,2) = Ja(1,1)
                Ja_inv = Ja_inv/detJ
                BB = matmul(Ja_inv, GN8)
                
                NN(1,1)=(1-G1)*(1-G2)*(-G1-G2-1)/4
                NN(1,2)=(1+G1)*(1-G2)*(G1-G2-1)/4
                NN(1,3)=(1+G1)*(1+G2)*(G1+G2-1)/4
                NN(1,4)=(1-G1)*(1+G2)*(-G1+G2-1)/4
                NN(1,5)=(1-G1*G1)*(1-G2)/2
                NN(1,6)=(1-G2*G2)*(1+G1)/2
                NN(1,7)=(1-G1*G1)*(1+G2)/2
                NN(1,8)=(1-G2*G2)*(1-G1)/2
! 因为可能写不成一行了，所以直接依次赋值了~
                GN8(1,1) = (1-G2)*(2*G1+G2)
                GN8(2,1) = (1-G1)*(G1+2*G2)
                GN8(1,2) = (1-G2)*(2*G1-G2)
                GN8(2,2) = (1+G1)*(2*G2-G1)
                GN8(1,3) = (1+G2)*(2*G1+G2)
                GN8(2,3) = (1+G1)*(2*G2+G1)
                GN8(1,4) = (1+G2)*(2*G1-G2)
                GN8(2,4) = (1-G1)*(2*G2-G1)
                GN8 = GN8/4
                GN8(1,5) = G1*(G2-1)
                GN8(2,5) = -(1-G1*G1)/2
                GN8(1,6) = (1-G2*G2)/2
                GN8(2,6) = -G2*(1+G1)
                GN8(1,7) = -G1*(1+G2)
                GN8(2,7) = (1-G1*G1)/2
                GN8(1,8) = (G2*G2-1)/2
                GN8(2,8) = G2*(G1-1)
                BB = matmul(Ja_inv, GN8)
! 为弯曲部分的Bk赋值，改成循环
                Bk = 0
                DO K = 1,8
                    Bk(1,3*K-1) = BB(1,K)
                    Bk(2,3*K)   = BB(2,K)
                    Bk(3,3*K-1) = BB(2,K)
                    Bk(3,3*K)   = BB(1,K)
                END DO
! 为剪切部分的By赋值。改成循环
                By = 0
                DO K = 1,8
                    By(1,3*K-2) = BB(1,K)
                    By(1,3*K-1) = -NN(1,K)
                    By(2,3*K-2) = BB(2,K)
                    By(2,3*K)   = -NN(1,K)
                END DO

            STR1 = -THICK(N)/2*matmul(Cb,matmul(Bk,DE))
            STR2(2:1:-1,:) = Cs*matmul(By, DE)
            WRITE (IOUT,"(5X,5E14.2)") STR1, STR2
            StressCollection(1:5, 9*N+3*L+M-12) = (/STR1, STR2/)
            GaussianCollection(1:3, 9*N+3*L+M-12) = reshape(matmul(reshape(XYZ(:,N), (/3,8/)),transpose(NN)), (/3/))
            END DO
        END DO
    END DO
    !call PostProcessor(NPAR(1), 2, XYZ((/((/3*k-2,3*k-1/),k=1,8)/),:), Node, 9, GaussianCollection(1:2,:), &
    !                   StressCollection, U)
  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

END SUBROUTINE PLATE8
