SUBROUTINE INFINITE4Q
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
  NPAR(5) = 4

! 此处材料要求每一种提供E, Possion

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: POSSION(NUMMAT)
! N103: LM(20,NUME)
! N104: XYZ(12,NUME)
! N105: MTAP(NUME)
! N106: THICK(NUME
! N107: NLAST
  N(1)=0
  N(2)=N(1)+NUMMAT*ITWO
  N(3)=N(2)+NUMMAT*ITWO
  N(4)=N(3)+8*NUME
  N(5)=N(4)+12*NUME*ITWO
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

  CALL INFINITE4 (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N(1)),A(N(2)),A(N(3)),A(N(4)),A(N(5)),A(N(6)),A(N(7)))

  RETURN

END SUBROUTINE INFINITE4Q


SUBROUTINE INFINITE4 (ID,X,Y,Z,U,MHT,E,POSSION,LM,XYZ,MATP, THICK, Node)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE
  USE MathKernel

  IMPLICIT NONE
  INTEGER :: ID(6,NUMNP),LM(8,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),POSSION(NPAR(3)),  &
             XYZ(12,NPAR(2)),THICK(NPAR(2)),U(NEQ)
  
  REAL(8) :: DE(8,1)
  INTEGER :: NPAR1, NUME, NUMMAT, ND, K, L, M, N
  INTEGER :: MTYPE, IPRINT, Node(NPAR(2),NPAR(5))

  REAL(8) :: Cb(3, 3),Cc(3, 3), Etemp, Ptemp, detJ, N1, N2, N3, N4
  REAL(8) :: GAUSS(2), W(2)
  REAL(8) :: G1, G2, GN(2,4), Ja(2,2), Ja_inv(2,2), Bm(3,8), S(8,8), BB(2,4), NShape(1,4)
  REAL(8) :: X_Y(4, 2), XY_G(1,2), STR1(3), STR2(2), GaussianCollection(3,NPAR(2)*4), StressCollection(6,NPAR(2)*4),NN0(1,4)
  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=8
  call GaussianMask(GAUSS, W, 2)
! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I10,/,   &
                   '     EQ.1, TRUSS ELEMENTS',/,      &
                   '     EQ.2, ELEMENTS CURRENTLY',/,  &
                   '     EQ.3, NOT AVAILABLE',//,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I10,/)") NPAR1,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I10,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S     CROSS-SECTIONAL',/,  &
                   ' NUMBER     MODULUS',10X,'AREA')")

     DO K =1,NUMMAT
        READ (IIN,'(I10,2F10.0)') N,E(N),POSSION(N)  ! Read material information
        WRITE (IOUT,"(I10,4X,E12.5,2X,E14.6)") N,E(N),POSSION(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE     NODE     NODE       MATERIAL',/,   &
                      ' NUMBER-N      I        J        K        L       SET NUMBER')")
     
     N=0
     LM = 0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(6I10, F10.0, I10)') N,Node(N,1:NPAR(5)),MTYPE,THICK(N)  ! Read in element information

!       Save element information
        XYZ(1:NPAR(5)*3-1:3,N)=X(Node(N,:))  ! Coordinates of the element's nodes
        XYZ(2:NPAR(5)*3  :3,N)=Y(Node(N,:))
        XYZ(3:NPAR(5)*3+1:3,N)=Z(Node(N,:))

        MATP(N)=MTYPE  ! Material type
        !DO L=1,5
        !   LM(L,N)=ID(L,I)     ! Connectivity matrix
        !   LM(L+5,N)=ID(L,J)
        !   LM(L+10,N)=ID(L,K)
        !   LM(L+15,N)=ID(L,L)
        !END DO
        ! 这里是为了让顺序

        DO L=1,2
           LM(L  ,N) = ID(L,Node(N,1))     ! Connectivity matrix
           LM(L+2,N) = ID(L,Node(N,2))
           LM(L+4,N) = ID(L,Node(N,3))
           LM(L+6,N) = ID(L,Node(N,4))
        END DO
        
!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I10,6X,I10,4X,I10,4X,I10,4X,I10,7X,I10)") N,Node(N,1:NPAR(5)),MTYPE
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
        Cc = Cb*THICK(N)*Etemp/(1-Ptemp*Ptemp)
! Gauss 积分常数
        S = 0
        DO L=1,2
            DO M=1,2
                G1 = GAUSS(L)
                G2 = GAUSS(M)
! 计算Jacobian
                !GN = reshape((/G2-1,G1-1, 1-G2,-G1-1, 1+G2,1+G1, -G2-1,1-G1/), shape(GN))/4
                GN(1,1) = -(1+G2)/(1-G1)**2
                GN(1,2) = -(1-G2)/(1-G1)**2
                GN(1,3) = (1-G2)/(1-G1)**2
                GN(1,4) = (1+G2)/(1-G1)**2
                GN(2,1) = -G1/(1-G1)
                GN(2,2) = G1/(1-G1)
                GN(2,3) = -(1+G1)/(2*(1-G1))
                GN(2,4) = (1+G1)/(2*(1-G1))
                Ja = matmul(GN,X_Y)
                detJ = Det(Ja,2)
                Ja_inv(1,1) = Ja(2,2)
                Ja_inv(2,1) = -Ja(2,1)
                Ja_inv(1,2) = -Ja(1,2)
                Ja_inv(2,2) = Ja(1,1)
                Ja_inv = Ja_inv/abs(detJ)
                BB = matmul(Ja_inv, GN)
                
!                NN0(1,1)=(1-G1)*(1-G2)/4
!                NN0(1,2)=(1+G1)*(1-G2)/4
!                NN0(1,3)=(1+G1)*(1+G2)/4
!                NN0(1,4)=(1-G1)*(1+G2)/4
!! 为弯曲部分的Bk赋值，改成循环
!                Bk = 0
!                DO K = 1,4
!                    Bk(1,5*K-3) = BB(1,K)
!                    Bk(2,5*K-2)   = BB(2,K)
!                    Bk(3,5*K-3) = BB(2,K)
!                    Bk(3,5*K-2)   = BB(1,K)
!                END DO
!! 为剪切部分的By赋值。改成循环
!                By = 0
!                DO K = 1,4
!                    By(1,5*K-4) = BB(1,K)
!                    By(1,5*K-3) = -NN0(1,K)
!                    By(2,5*K-4) = BB(2,K)
!                    By(2,5*K-2)   = -NN0(1,K)
!                END DO
                    !By(1,2) = -GN(1,1)
                    !By(1,7) = -GN(1,1)
                    !By(1,12) = -GN(1,4)
                    !By(1,17) = -GN(1,4)
                    !By(2,3) = -GN(2,1)
                    !By(2,8) = -GN(2,2)
                    !By(2,13) = -GN(2,2)
                    !By(2,18) = -GN(2,1)
! 为平面部分的Bm赋值。改成=循环
                Bm = 0
                DO K = 1,4
                    Bm(1,2*K-1) = BB(1,K)
                    Bm(2,2*K) = BB(2,K)
                    Bm(3,2*K-1) = BB(2,K)
                    Bm(3,2*K)   = BB(1,K)
                END DO
! 这里不要忘了还要乘上z方向积分
                S = S + W(L)*W(M)*(matmul(matmul(transpose(Bm), Cc), Bm))*abs(detJ)
            END DO
        END DO

        CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)  ! 这里要输出的S就是制作好了的local stiffness matrix
        
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
        Cb = Cb*Etemp/12.0/(1-Ptemp*Ptemp)*5.0/6.0

        Cc = Cb
        !
        !Cs = Etemp/(2*(1+Ptemp))
! Gauss 积分常数
        S = 0
        DO L=1,2
            DO M=1,2
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
                BB = matmul(Ja_inv, GN)
! 为弯曲部分的Bk赋值，改成循环
!                Bk = 0
!                DO K = 1,4
!                    Bk(1,3*K-1) = BB(1,K)
!                    Bk(2,3*K)   = BB(2,K)
!                    Bk(3,3*K-1) = BB(2,K)
!                    Bk(3,3*K)   = BB(1,K)
!                END DO
!! 为剪切部分的By赋值。改成循环
!                By = 0
!                DO K = 1,4
!                    By(1,5*K-4) = BB(1,K)
!                    By(1,5*K-3) = -NN0(1,K)
!                    By(2,5*K-4) = BB(2,K)
!                    By(2,5*K-2)   = -NN0(1,K)
!                END DO
! 为平面部分的Bm赋值。改成循环
                Bm = 0
                DO K = 1,4
                    Bm(1,2*K-1) = BB(1,K)
                    Bm(2,2*K) = BB(2,K)
                    Bm(3,2*K-1) = BB(2,K)
                    Bm(3,2*K)   = BB(1,K)
                END DO

                STR1 = RESHAPE(matmul(Cc,matmul(Bm, DE)), (/3/))
                !str2(2:1:-1) = reshape(cs*matmul(by, de), (/2/))
                WRITE (IOUT,"(5X,5E14.2)") STR1, STR2
                
                !N Matrix Elements
                N1 = 0.25*(1.0-G2)*(1.0-G1)
                N2 = 0.25*(1.0+G2)*(1.0-G1)
                N3 = 0.25*(1.0+G2)*(1.0+G1)
                N4 = 0.25*(1.0-G2)*(1.0+G1)
                
                GaussianCollection(1:3, N*4+2*L+M-6) = matmul(reshape(XYZ(:,N),(/3,4/)), (/N1, N2, N3, N4/))
                StressCollection(1:5, N*4+2*L+M-6) = (/STR1, STR2/)
            END DO
        END DO
    END DO
    
    call PostProcessor(NPAR(1), 2, XYZ, Node, 4, GaussianCollection, StressCollection, U)
  
  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

END SUBROUTINE INFINITE4
