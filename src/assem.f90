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

SUBROUTINE COLHT (MHT,ND,LM)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To calculate column heights                                     .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : NEQ
  IMPLICIT NONE
  INTEGER :: ND, LM(ND),MHT(NEQ)
  INTEGER :: I, LS, II, ME

  LS=HUGE(1)   ! The largest integer number

  DO I=1,ND
     IF (LM(I) .NE. 0) THEN
        IF (LM(I)-LS .LT. 0) LS=LM(I)
     END IF
  END DO

  DO I=1,ND
     II=LM(I)
     IF (II.NE.0) THEN
        ME=II - LS
        IF (ME.GT.MHT(II)) MHT(II)=ME
     END IF
  END DO

  RETURN
END SUBROUTINE COLHT


SUBROUTINE ADDRES (MAXA,MHT)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To calculate addresses of diagonal elements in banded           .
! .   matrix whose column heights are known                           .
! .                                                                   .
! .   MHT  = Active column heights                                    .
! .   MAXA = Addresses of diagonal elements                           .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : NEQ, MK, NWK

  IMPLICIT NONE
  INTEGER :: MAXA(NEQ+1),MHT(NEQ)
  INTEGER :: NN, I

! Clear array maxa

  NN=NEQ + 1
  DO I=1,NN
     MAXA(I)=0.0
  END DO

  MAXA(1)=1
  MAXA(2)=2
  MK=0
  IF (NEQ.GT.1) THEN
     DO I=2,NEQ
        IF (MHT(I).GT.MK) MK=MHT(I)
        MAXA(I+1)=MAXA(I) + MHT(I) + 1
     END DO
  END IF
  MK=MK + 1
  NWK=MAXA(NEQ+1) - MAXA(1)

  RETURN
END SUBROUTINE ADDRES


SUBROUTINE ASSEM (AA)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To call element subroutines for assemblage of the               .
! .   structure stiffness matrix                                      .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : IELMNT, NUMEG, MIDEST, NPAR

  IMPLICIT NONE
  REAL :: AA(*)
  INTEGER :: N, I

  REWIND IELMNT
  DO N=1,NUMEG
     READ (IELMNT) MIDEST,NPAR,(AA(I),I=1,MIDEST)
     CALL ELEMNT
  END DO

  RETURN
END SUBROUTINE ASSEM


SUBROUTINE ADDBAN (A,MAXA,S,LM,ND)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To assemble element stiffness into compacted global stiffness   .
! .                                                                   .
! .      A = GLOBAL STIFFNESS (1D skyline storage)                    .
! .      S = ELEMENT STIFFNESS                                        .
! .      ND = DEGREES OF FREEDOM IN ELEMENT STIFFNESS                 .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  USE GLOBALS, ONLY : NWK, NEQ
  IMPLICIT NONE
  REAL(8) :: A(NWK),S(ND,ND)
  INTEGER :: MAXA(NEQ+1),LM(ND)
  INTEGER :: ND, I, J, II, JJ, KK
  
  DO J=1,ND
     JJ=LM(J)
     IF (JJ .GT. 0) THEN
        DO I=1,J
           II=LM(I)
           IF (II .GT. 0) THEN
              IF (JJ .GE. II) THEN
                 KK= MAXA(JJ) + JJ - II
              ELSE
                 KK= MAXA(II) + II - JJ
              END IF              
              A(KK)=A(KK) + S(I,J)
           END IF
        END DO
     END IF
  END DO

  RETURN
    END SUBROUTINE ADDBAN

subroutine pardiso_addban(A, rowIndex, columns, S, LM, ND)
  USE GLOBALS, ONLY : NWK, NEQ
  IMPLICIT NONE
  REAL(8) :: A(NWK),S(ND,ND)
  INTEGER :: rowIndex(NEQ+1),columns(nwk), LM(ND)
  INTEGER :: ND, I, J, II, JJ, k, KK, tempJ
  DO J=1,ND
     JJ=LM(J)
     IF (JJ .GT. 0) THEN
        DO I=J,ND
           II=LM(I)
           IF (II .GT. 0) THEN
              IF (JJ .GT. II) THEN ! 如果JJ > II，交换，让JJ永远小于II
                 tempJ = II
              else
                 tempJ = JJ                  
              END IF              
              loop1: do k = rowIndex(tempJ), rowIndex(tempJ+1)-1
                 if(columns(k) .EQ. II) then
                     A(k) = A(k) + S(I,J)
                     exit loop1
                 end if
             end do loop1
           END IF
        END DO
     END IF
  END DO

  RETURN
end subroutine
  
SUBROUTINE COLSOL (A,V,MAXA,NN,NWK,NNM,KKK)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To solve finite element static equilibrium equations in         .
! .   core, using compacted storage and column reduction scheme       .
! .                                                                   .
! .  - - Input variables - -                                          .
! .        A(NWK)    = Stiffness matrix stored in compacted form      .
! .        V(NN)     = Right-hand-side load vector                    .
! .        MAXA(NNM) = Vector containing addresses of diagonal        .
! .                    elements of stiffness matrix in a              .
! .        NN        = Number of equations                            .
! .        NWK       = Number of elements below skyline of matrix     .
! .        NNM       = NN + 1                                         .
! .        KKK       = Input flag                                     .
! .            EQ. 1   Triangularization of stiffness matrix          .
! .            EQ. 2   Reduction and back-substitution of load vector .
! .        IOUT      = UNIT used for output                           .
! .                                                                   .
! .  - - OUTPUT - -                                                   .
! .        A(NWK)    = D and L - Factors of stiffness matrix          .
! .        V(NN)     = Displacement vector                            .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : IOUT

  IMPLICIT NONE
  INTEGER :: MAXA(NNM),NN,NWK,NNM,KKK
  REAL(8) :: A(NWK),V(NN),C,B, D(10, 10)
  INTEGER :: N,K,KN,KL,KU,KH,IC,KLT,KI,J,ND,KK,L
  INTEGER :: MIN0

! Perform L*D*L(T) factorization of stiffness matrix

  IF (KKK == 1) THEN

      DO N=1,NN
         KN=MAXA(N)
         KL=KN + 1
         KU=MAXA(N+1) - 1
         KH=KU - KL

         IF (KH > 0) THEN
             K=N - KH
             IC=0
             KLT=KU
             DO J=1,KH
                IC=IC + 1
                KLT=KLT - 1
                KI=MAXA(K)
                ND=MAXA(K+1) - KI - 1
                IF (ND .GT. 0) THEN
                   KK=MIN0(IC,ND)
                   C=0.
                   DO L=1,KK
                      C=C + A(KI+L)*A(KLT+L)
                   END DO
                   A(KLT)=A(KLT) - C
                END IF
                K=K + 1
             END DO
         ENDIF

         IF (KH >= 0) THEN
             K=N
             B=0.
             DO KK=KL,KU
                K=K - 1
                KI=MAXA(K)
                C=A(KK)/A(KI)
                B=B + C*A(KK)
                A(KK)=C
             END DO
             A(KN)=A(KN) - B
         ENDIF

         IF (A(KN) .LE. 0) THEN
            WRITE (IOUT,"(//' STOP - STIFFNESS MATRIX NOT POSITIVE DEFINITE',//,  &
                            ' NONPOSITIVE PIVOT FOR EQUATION ',I8,//,' PIVOT = ',E20.12 )") N,A(KN)
            STOP
         END IF
      END DO

  ELSE IF (KKK == 2) THEN

! REDUCE RIGHT-HAND-SIDE LOAD VECTOR

       DO N=1,NN
         KL=MAXA(N) + 1
         KU=MAXA(N+1) - 1
         IF (KU-KL .GE. 0) THEN
            !K=N
            !C=0.
            !DO KK=KL,KU
            !   K=K - 1
            !   C=C + A(KK)*V(K)
            !END DO
            C = dot_product(A(KL:KU),V(N-1:N-(KU-KL)-1:-1))
            V(N)=V(N) - C
         END IF
      END DO

! BACK-SUBSTITUTE

      DO N=1,NN
         K=MAXA(N)
         V(N)=V(N)/A(K)
      END DO

      IF (NN.EQ.1) RETURN

      N=NN
      DO L=2,NN
         KL=MAXA(N) + 1
         KU=MAXA(N+1) - 1
         IF (KU-KL .GE. 0) THEN
            K=N
            DO KK=KL,KU
               K=K - 1
               V(K)=V(K) - A(KK)*V(N)
            END DO
         END IF
         N=N - 1
      END DO
       D = reshape((/ 3.43d8,1.24d8,0.d0,0.d0,0.d0,-2.09d8,-9.53d6,0.d0,0.d0,0.d0,   &
                1.24d8, 3.43d8, 0.d0, 0.d0, 0.d0, 9.53d6, 3.81d7, 0.d0, 0.d0, 0.d0,   &
                0.d0, 0.d0, 2.13d9, -1.60d9, -1.60d9, 0.d0, 0.d0, -5.34d8, -1.60d9, -1.60d9,  &
                0.d0, 0.d0, -1.60d9, 3.23d9, 1.03d7, 0.d0, 0.d0, 1.60d9, 3.18d9, -7.94d5 ,   &
                0.d0, 0.d0, -1.60d9, 1.03d7, 3.23d9, 0.d0, 0.d0, -1.60d9, 7.94d5, 3.20d9,   &
                -2.09d8, 9.53d6, 0.d0, 0.d0, 0.d0,3.43d8, -1.24d8, 0.d0, 0.d0, 0.d0, &
                -9.53d6, 3.81d7, 0.d0, 0.d0, 0.d0,-1.24d8, 3.43d8, 0.d0, 0.d0, 0.d0, &
                0.d0, 0.d0, -5.34d8, 1.60d9, -1.60d9, 0.d0, 0.d0, 2.13d9, 1.60d9, -1.60d9,&
                0.d0, 0.d0, -1.60d9, 3.18d9, 7.94d5, 0.d0, 0.d0, 1.60d9, 3.23d9, -1.03d7, &
                0.d0, 0.d0, -1.60d9, -7.94d5, 3.20d9, 0.d0, 0.d0, 1.60d9, -1.03d7, 3.23d9               /), shape(D))
  V = matmul(D, V)
  END IF

    END SUBROUTINE COLSOL

subroutine pardiso_solver(D,V, rowIndex, columns)
USE mkl_pardiso
USE GLOBALS, only : neq, nwk
IMPLICIT NONE
!.. Internal solver memory pointer 
TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
!.. All other variables
INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl, nnz
INTEGER error1
INTEGER, ALLOCATABLE :: iparm( : )
INTEGER, ALLOCATABLE :: ia( : )
INTEGER, ALLOCATABLE :: ja( : )
REAL(8), ALLOCATABLE :: a( : )
REAL(8), ALLOCATABLE :: b( : )
REAL(8), ALLOCATABLE :: x( : )
INTEGER i, idum(1)
REAL(8) ddum(1)
REAL(8) C(10,10)
REAL(8) D(nwk), V(neq)
INTEGER columns(nwk), rowIndex(neq+1)
!.. Fill all arrays containing matrix data.
n = 8 
nnz = 36

ALLOCATE(ia(n + 1))
ia = (/ 1, 9, 16, 22, 27, 31, 34, 36, 37 /)
ALLOCATE(ja(nnz))
ja = (/ 1, 2, 3, 4, 5, 6, 7, 8, &
           2, 3, 4, 5, 6, 7, 8, &
              3, 4, 5, 6, 7, 8, &
                 4, 5, 6, 7, 8, &
                    5, 6, 7, 8, &
                       6, 7, 8, &
                          7, 8, &
                             8 /)
ALLOCATE(a(nnz))
 C = reshape((/ 3.43d8,1.24d8,0.d0,0.d0,0.d0,-2.09d8,-9.53d6,0.d0,0.d0,0.d0,   &
                1.24d8, 3.43d8, 0.d0, 0.d0, 0.d0, 9.53d6, 3.81d7, 0.d0, 0.d0, 0.d0,   &
                0.d0, 0.d0, 2.13d9, -1.60d9, -1.60d9, 0.d0, 0.d0, -5.34d8, -1.60d9, -1.60d9,  &
                0.d0, 0.d0, -1.60d9, 3.23d9, 1.03d7, 0.d0, 0.d0, 1.60d9, 3.18d9, -7.94d5 ,   &
                0.d0, 0.d0, -1.60d9, 1.03d7, 3.23d9, 0.d0, 0.d0, -1.60d9, 7.94d5, 3.20d9,   &
                -2.09d8, 9.53d6, 0.d0, 0.d0, 0.d0,3.43d8, -1.24d8, 0.d0, 0.d0, 0.d0, &
                -9.53d6, 3.81d7, 0.d0, 0.d0, 0.d0,-1.24d8, 3.43d8, 0.d0, 0.d0, 0.d0, &
                0.d0, 0.d0, -5.34d8, 1.60d9, -1.60d9, 0.d0, 0.d0, 2.13d9, 1.60d9, -1.60d9,&
                0.d0, 0.d0, -1.60d9, 3.18d9, 7.94d5, 0.d0, 0.d0, 1.60d9, 3.23d9, -1.03d7, &
                0.d0, 0.d0, -1.60d9, -7.94d5, 3.20d9, 0.d0, 0.d0, 1.60d9, -1.03d7, 3.23d9               /), shape(C))


a = (/ 7.d0, 0.d0,  1.d0, 0.d0, 0.d0, 2.d0, 7.d0, 0.d0,  &
             -4.d0, 8.d0, 0.d0, 2.d0, 0.d0, 0.d0, 0.d0,  &
                    1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 5.d0,  &
                          7.d0, 0.d0, 0.d0, 9.d0, 0.d0,  &
                                5.d0, 1.d0, 5.d0, 0.d0,  &
                                     -1.d0, 0.d0, 5.d0,  &
                                           11.d0, 0.d0,  &
                                                  5.d0 /)
ALLOCATE(b(n))
ALLOCATE(x(neq))
!..
!.. Set up PARDISO control parameter
!..
ALLOCATE(iparm(64))

DO i = 1, 64
   iparm(i) = 0
END DO

iparm(1) = 1 ! no solver default
iparm(2) = 2 ! fill-in reordering from METIS
iparm(4) = 0 ! no iterative-direct algorithm
iparm(5) = 0 ! no user fill-in reducing permutation
iparm(6) = 0 ! =0 solution on the first n components of x
iparm(8) = 2 ! numbers of iterative refinement steps
iparm(10) = 13 ! perturb the pivot elements with 1E-13
iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
iparm(14) = 0 ! Output: number of perturbed pivots
iparm(18) = -1 ! Output: number of nonzeros in the factor LU
iparm(19) = -1 ! Output: Mflops for LU factorization
iparm(20) = 0 ! Output: Numbers of CG Iterations

error  = 0 ! initialize error flag
msglvl = 1 ! print statistical information
mtype  = -2 ! symmetric, indefinite
nrhs = 1 
maxfct = 1 
mnum = 1
!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

ALLOCATE (pt(64))
DO i = 1, 64
   pt(i)%DUMMY =  0 
END DO

!.. Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization

phase = 11 ! only reordering and symbolic factorization

CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, D, rowIndex, columns, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error)
    
WRITE(*,*) 'Reordering completed ... '
IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
END IF
WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)

!.. Factorization.
phase = 22 ! only factorization
CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, D, rowIndex, columns, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error)
WRITE(*,*) 'Factorization completed ... '
IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
ENDIF

!.. Back substitution and iterative refinement
iparm(8) = 2 ! max numbers of iterative refinement steps
phase = 33 ! only solving
DO i = 1, n
   b(i) = 1.d0
END DO
CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, D, rowIndex, columns, &
              idum, nrhs, iparm, msglvl, V, x, error)
WRITE(*,*) 'Solve completed ... '
IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
ENDIF
WRITE(*,*) 'The solution of the system is '
DO i = 1, neq
   WRITE(*,*) ' x(',i,') = ', x(i)
END DO
x = matmul(C, x)
DO i = 1, neq
   WRITE(*,*) ' x(',i,') = ', x(i)
END DO      
1000 CONTINUE
!.. Termination and release of memory
phase = -1 ! release internal memory
CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, ddum, idum, idum, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error1)

IF (ALLOCATED(ia))      DEALLOCATE(ia)
IF (ALLOCATED(ja))      DEALLOCATE(ja)
IF (ALLOCATED(a))       DEALLOCATE(a)
IF (ALLOCATED(b))       DEALLOCATE(b)
IF (ALLOCATED(x))       DEALLOCATE(x)
IF (ALLOCATED(iparm))   DEALLOCATE(iparm)

IF (error1 /= 0) THEN
   WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
   STOP 1
ENDIF

IF (error /= 0) STOP 1
end subroutine pardiso_solver