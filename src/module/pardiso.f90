subroutine pardiso_crop(A, rowIndex, columns)
    use GLOBALS, only : NWK, NEQ
    real(8) :: A(nwk)
    integer :: rowIndex(neq+1)
    integer :: columns(nwk)
    integer :: i, j, k, temp !循环变量
    ! 因为第一个元素肯定有数值，所以从第二个开始
    j = 0
    k = 2
    temp = 0
    do i = 1,nwk
        if( i == rowIndex(k)) then
            rowIndex(k) = rowIndex(k-1) + temp
            k = k+1 !换下一行
            temp = 0 !清零计数
        end if
        if( abs(A(i)) .gt. 1.0d-5) then ! 如果非零就赋值
            j = j+1
            A(j) = A(i)
            columns(j) = columns(i)
            temp = temp+1
        end if
    end do
    rowIndex(neq + 1) = rowIndex(neq)+temp
    nwk = j
    end subroutine
    
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
              IF (JJ .GT. II) THEN ! 如果JJ > II
                 loop1: do k = rowIndex(II), rowIndex(II+1)-1
                    if(columns(k) .EQ. JJ) then
                        A(k) = A(k) + S(I,J)
                        exit loop1
                    end if
                end do loop1
              else ! 如果II < JJ
                  loop2: do k = rowIndex(JJ), rowIndex(JJ+1)-1
                    if(columns(k) .EQ. II) then
                        A(k) = A(k) + S(I,J)
                        exit loop2
                    end if
                end do loop2                  
              END IF              
           END IF
        END DO
     END IF
  END DO

  RETURN
end subroutine   

subroutine pardiso_solver(K,V, rowIndex, columns)
USE mkl_pardiso
USE GLOBALS, only : neq, nwk
IMPLICIT NONE
!.. Internal solver memory pointer 
TYPE(MKL_PARDISO_HANDLE)  :: pt(64)
!.. All other variables
INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
INTEGER error1
INTEGER :: iparm(64)
REAL(8) x(neq)
INTEGER i, idum(1)
REAL(8) ddum(1)
REAL(8) K(nwk), V(neq)
INTEGER columns(nwk), rowIndex(neq+1)
!..
!.. Set up PARDISO control parameter
!..
iparm = 0

iparm(1) = 1 ! no solver default
iparm(2) = 3 ! fill-in reordering from METIS
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
mtype  = 2 ! symmetric, indefinite
nrhs = 1 
maxfct = 1 
mnum = 1
!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

DO i = 1, 64
   pt(i)%DUMMY =  0 
END DO

!.. Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization

phase = 11 ! only reordering and symbolic factorization

CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, K, rowIndex, columns, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error)
    
WRITE(*,*) 'Reordering completed ... '
IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
END IF
!WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
!WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)

!.. Factorization.
!phase = 22 ! only factorization
!CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, K, rowIndex, columns, &
!              idum, nrhs, iparm, msglvl, ddum, ddum, error)
!WRITE(*,*) 'Factorization completed ... '
!IF (error /= 0) THEN
!   WRITE(*,*) 'The following ERROR was detected: ', error
!   GOTO 1000
!ENDIF

!.. Back substitution and iterative refinement
!iparm(8) = 2 ! max numbers of iterative refinement steps
!phase = 33 ! only solving
!CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, K, rowIndex, columns, &
!              idum, nrhs, iparm, msglvl, V, x, error)
!WRITE(*,*) 'Solve completed ... '
!IF (error /= 0) THEN
!   WRITE(*,*) 'The following ERROR was detected: ', error
!   GOTO 1000
!ENDIF

iparm(8) = 2 ! max numbers of iterative refinement steps
phase = 23 ! only solving
CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, K, rowIndex, columns, &
              idum, nrhs, iparm, msglvl, V, x, error)

V = x
1000 CONTINUE
!.. Termination and release of memory
phase = -1 ! release internal memory
CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, ddum, idum, idum, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error1)

IF (error1 /= 0) THEN
   WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
   STOP 1
ENDIF

IF (error /= 0) STOP 1
end subroutine pardiso_solver