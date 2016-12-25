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

SUBROUTINE TRUSS
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the truss element subroutine         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N(8)

  NUME = NPAR(2)
  NUMMAT = NPAR(3)
  NPAR(5) = 2

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: AREA(NUMMAT)
! N103: LM(6,NUME)
! N104: XYZ(6,NUME)
! N105: MTAP(NUME)
  N(1)=0
  N(2)=N(1)+NUMMAT*ITWO
  N(3)=N(2)+NUMMAT*ITWO
  if (DYNANALYSIS) then
        N(4) = N(3)+NPAR(3)*ITWO
  else
        N(4) = N(3)
  end if
  N(5)=N(4)+6*NUME
  N(6)=N(5)+6*NUME*ITWO
  N(7)=N(6)+NUME
  N(8)=N(7)+NPAR(5)*NPAR(2)
  
  MIDEST=N(8)
  if (IND .EQ. 1) then
        ! Allocate storage for element group data
        call MemAlloc(11,"ELEGP",MIDEST,1)
  end if
  NFIRST = NP(11)   ! Pointer to the first entry in the element group data array in the unit of single precision (corresponding to A)
  N(:) = N(:) + NFIRST
  NLAST=N(8)

  CALL RUSS (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
        A(N(1)),A(N(2)),A(N(3)),A(N(4)),A(N(5)),A(N(6)),A(N(7)))

  RETURN

END SUBROUTINE TRUSS


SUBROUTINE RUSS (ID,X,Y,Z,U,MHT,E,AREA,Density,LM,XYZ,MATP,Node)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(6,NUMNP),LM(6,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),AREA(NPAR(3)),  &
             XYZ(6,NPAR(2)),U(NEQ), M(6,6), Rho, Density(NPAR(3))
  REAL(8) :: S(6,6),ST(6),D(3), StressCollection(6,NPAR(2)), GaussianCollection(3,NPAR(2))

  INTEGER :: NPAR1, NUME, NUMMAT, ND, I, J, L, N, Node(NPAR(2),NPAR(5))
  INTEGER :: MTYPE, IPRINT
  REAL(8) :: XL2, XL, SQRT, XX, YY, STR, P

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=6

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

     WRITE (IOUT,"('  SET       YOUNG''S     CROSS-SECTIONAL   DENSITY',/,  &
                   ' NUMBER     MODULUS',10X,    'AREA',/,  &
                   15 X,'E',14X,                  'A',14X,       '¦Ñ')")

     if (DYNANALYSIS) then
        DO I=1,NUMMAT
            READ (IIN,'(I10,3F10.0)') N,E(N), AREA(N), Density(N)      ! Read Density
            WRITE (IOUT,"(I10,4X,E12.5,2(2X,E14.6))") N,E(N), AREA(N), Density(N)
        END DO
     else
        DO I=1,NUMMAT
            READ (IIN,'(I10,2F10.0)') N,E(N),AREA(N)  ! Read material information
            WRITE (IOUT,"(I10,4X,E12.5,2X,E14.6)") N,E(N),AREA(N)
        END DO
     end if

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE       MATERIAL',/,   &
                      ' NUMBER-N      I        J       SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(5I10)') N,Node(N,1:2),MTYPE  ! Read in element information

!       Save element information
        XYZ(1:NPAR(5)*3-1:3,N)=X(Node(N,:))  ! Coordinates of the element's nodes
        XYZ(2:NPAR(5)*3  :3,N)=Y(Node(N,:))
        XYZ(3:NPAR(5)*3+1:3,N)=Z(Node(N,:))
        MATP(N)=MTYPE  ! Material type

        DO L=1,6
           LM(L,N)=0
        END DO

        DO L=1,3
           LM(L,N)=ID(L,Node(N,1))     ! Connectivity matrix
           LM(L+3,N)=ID(L,Node(N,2))
        END DO

!       Update column heights and bandwidth
        if (.NOT. PARDISODOOR) CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I10,6X,I10,4X,I10,7X,I10)") N,Node(N,1:2),MTYPE
        write (VTKNodeTmp) NPAR(5), Node(N,:)-1

     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     DO N=1,NUME
        MTYPE=MATP(N)

        XL2=0.
        DO L=1,3
           D(L)=XYZ(L,N) - XYZ(L+3,N)
           XL2=XL2 + D(L)*D(L)
        END DO
        XL=SQRT(XL2)   ! Length of element N

        XX=E(MTYPE)*AREA(MTYPE)*XL   !  E*A*l

        DO L=1,3
           ST(L)=D(L)/XL2
           ST(L+3)=-ST(L)
        END DO

        DO J=1,ND
           YY=ST(J)*XX
           DO I=1,ND
              S(I,J)=ST(I)*YY
           END DO
        END DO
        
        if (DYNANALYSIS) then                           !Mass Matrix
            XX = Density(MTYPE)*AREA(MTYPE)*XL/3        !Density*A*l/3
            DO L=1,3
                ST(L)=D(L)/XL2                          ![ 2  1 ]
                ST(L+3)=ST(L)/2                         ![ 1  2 ]
            END DO

            DO J=1,ND
                YY=ST(J)*XX
                DO I=1,ND
                    M(I,J)=ST(I)*YY
                END DO
            END DO
        end if

        if(pardisodoor) then
<<<<<<< HEAD
            if(huge) then
                call pardiso_addban(stff,IA(NP(2)),columns,S,LM(1,N),ND)
            else
                call pardiso_addban(DA(NP(3)),IA(NP(2)),IA(NP(5)),S,LM(1,N),ND)
            end if
=======
            call pardiso_addban(DA(NP(3)),IA(NP(2)),IA(NP(5)),S,LM(1,N),ND)
            if (DYNANALYSIS) CALL pardiso_addban(DA(NP(10)),IA(NP(9)), IA(NP(8)),M,LM(:,N),ND)
>>>>>>> c85d1edd52a173423583b01cc6524492827c80cd
        else
            CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)
            IF (DYNANALYSIS) CALL ADDBAN (DA(NP(10)),IA(NP(2)),M,LM(:,N),ND)
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
                                           '  ELEMENT',13X,'FORCE',12X,'STRESS',/,'  NUMBER')") NG
        MTYPE=MATP(N)

        XL2=0.
        DO L=1,3
           D(L) = XYZ(L,N) - XYZ(L+3,N)
           XL2=XL2 + D(L)*D(L)
        END DO

        DO L=1,3
           ST(L)=(D(L)/XL2)*E(MTYPE)
           ST(L+3)=-ST(L)
        END DO

        STR=0.0
        DO L=1,3
           I=LM(L,N)
           IF (I.GT.0) STR=STR + ST(L)*U(I)

           J=LM(L+3,N)
           IF (J.GT.0) STR=STR + ST(L+3)*U(J)
        END DO

        P=STR*AREA(MTYPE)

        WRITE (IOUT,"(1X,I10,11X,E13.6,4X,E13.6)") N,P,STR
        GaussianCollection(:,N) = 0.5*(XYZ(4:6,N)+XYZ(1:3,N))
        StressCollection(1,N) = STR
     END DO
     StressCollection(2,:) = 0D0
     call PostProcessor(NPAR(1), 1, XYZ, &
                        Node, 1, GaussianCollection, StressCollection, U)
                      ! AT POINT 513 & 1370 WHERE TRUSSES MERGE, DO NOT USE LEAST SQUARE!
  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

END SUBROUTINE RUSS
