MODULE MathKernel

!=======================================================================================
! Common Math Operation Collections
! Includes:
!       Det : Detting determinant of a matrix
!       InvMat : General purpose matrix inversion
!       InvMat3 : Matrix inversion for 3x3 Matrix
!       InvMat2 : Matrix inversion for 2x2 Matrix
!       Cross : Vector cross product
!       IsPlanar : Judging if 4 points are on the same plane
!       GaussianMask : Generating gaussian points and weights for caritsian coordinates
!       AxialRotate : Rotate a group of points along a given axis for given angle
!       LeastSquare : Return the least-square coefficients for a given series of data
! 
! Author: Haoguang Yang
! 
!=======================================================================================

    
    implicit none

CONTAINS
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Function to find the determinant of a square matrix   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) FUNCTION Det(matrix, n)
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
    IMPLICIT NONE
    REAL(8), DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL(8) :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    
    select case (n)
    case (2)
        Det = matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
    case (3)
        Det = matrix(1,1)*matrix(2,2)*matrix(3,3)+ &
              matrix(1,2)*matrix(2,3)*matrix(3,1)+ &
              matrix(1,3)*matrix(2,1)*matrix(3,2)- &
              matrix(1,1)*matrix(2,3)*matrix(3,2)- &
              matrix(1,2)*matrix(2,1)*matrix(3,3)- &
              matrix(1,3)*matrix(2,2)*matrix(3,1)
    case (4:)
        l = 1
        !Convert to upper triangular form
        DO k = 1, n-1
            IF (matrix(k,k) == 0) THEN
                DetExists = .FALSE.
                DO i = k+1, n
                    IF (matrix(i,k) /= 0) THEN
                        DO j = 1, n
                            temp = matrix(i,j)
                            matrix(i,j)= matrix(k,j)
                            matrix(k,j) = temp
                        END DO
                        DetExists = .TRUE.
                        l=-l
                        EXIT
                    ENDIF
                END DO
                IF (DetExists .EQV. .FALSE.) THEN
                    Det = 0
                    return
                END IF
            ENDIF
            DO j = k+1, n
                m = matrix(j,k)/matrix(k,k)
                DO i = k+1, n
                    matrix(j,i) = matrix(j,i) - m*matrix(k,i)
                END DO
            END DO
        END DO
        !Calculate determinant by finding product of diagonal elements
        Det = l
        DO i = 1, n
            Det = Det * matrix(i,i)
        END DO
    END select
END FUNCTION Det


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     Functions To Calculate The Inverse Of A Small Matrix     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!2x2
SUBROUTINE InvMat2 (A, AINV, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(2,2), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(2,2), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(2,2) :: COFACTOR


      DET =   A(1,1)*A(2,2) - A(1,2)*A(2,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +A(2,2)
      COFACTOR(1,2) = -A(2,1)
      COFACTOR(2,1) = -A(1,2)
      COFACTOR(2,2) = +A(1,1)

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

END SUBROUTINE InvMat2
      
!3x3
SUBROUTINE InvMat3 (A, AINV, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

END SUBROUTINE InvMat3

!General inverse
subroutine InvMat(MatIn,MatOut,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision MatIn(n,n), MatOut(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=MatIn(i,k)/MatIn(k,k)
      L(i,k) = coeff
      do j=k+1,n
         MatIn(i,j) = MatIn(i,j)-coeff*MatIn(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = MatIn(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    MatOut(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine InvMat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      Function calculating cross product of 3D vectors        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION cross(a, b)
  implicit none
  real(8), DIMENSION(3) :: cross
  real(8), DIMENSION(3), INTENT(IN) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Subroutine Judging if four points are in the same place   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine IsPlanar(Judge, NormalVec, NPoints, PositionData)
    implicit none
    integer ::  NPoints
    real(8) ::  PositionData(NPoints*3), NormalVec(3), V(3,3), Criteria
    logical ::  Judge
    
    Judge = .False.
    if (NPoints == 4) then
        V=reshape((/PositionData(4:11:3)-PositionData(1)*(/1,1,1/), &
            PositionData(5:12:3)-PositionData(2)*(/1,1,1/), &
            PositionData(6:13:3)-PositionData(3)*(/1,1,1/)/),shape(V))              !problematic expression
        NormalVec = cross(V(2,:),V(3,:)-V(1,:))
        Criteria  = dot_product(NormalVec,V(1,:))
        
        if (abs(Criteria) .LT. 1.0D-10) then
            Judge = .True.
            return
        endif
    endif
end subroutine IsPlanar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Function Returning Weight Masks for Gaussian Quadrature    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GaussianMask (GaussianPoints, Weights, QuadratureOrder)
    implicit none
    !integer ::  EdgesPerNode, NumberOfEdges
    integer ::  QuadratureOrder
    real(8) ::  GaussianPoints (QuadratureOrder), Weights(QuadratureOrder)
    
    select case (QuadratureOrder)
    case (1)
        GaussianPoints=0.
        Weights=2.
    case (2)
        GaussianPoints=(/1./sqrt(3.), -1./sqrt(3.)/)
        Weights=(/1.,1./)
    case (3)
        GaussianPoints=(/sqrt(3./5.), 0., -sqrt(3./5.)/)
        Weights=(/5./9., 8./9., 5./9./)
    !case (4)
    !Other Higher Orders
    
    end select
end subroutine GaussianMask

!subroutine GaussianMaskSimplex (GaussianPts, Weights, QuadratureOrder)
!    implicit none
!    integer ::  QuadratureOrder
!    real(8) ::  GaussianPts(npg)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Subroutine performing spacial rotations along a given axis   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine AxialRotate(PtIn, PtOut, Axis, Theta)
	implicit none
	real(8) :: PtIn(3), PtOut(3), Axis(3), Theta
	real(8) :: C, S, R(3,3)
	
	C = cos(Theta)
	S = sin(Theta)
	Axis = Axis/sqrt(dot_product(Axis, Axis))
	R = reshape((/Axis(1)**2+(1-Axis(1)**2)*C, Axis(1)*Axis(2)*(1-C)-Axis(3)*S, Axis(1)*Axis(3)*(1-C)+Axis(2)*S, &
				  Axis(1)*Axis(2)*(1-C)+Axis(3)*S, Axis(2)**2+(1-Axis(2)**2)*C, Axis(2)*Axis(3)*(1-C)-Axis(1)*S, &
				  Axis(1)*Axis(3)*(1-C)-Axis(2)*S, Axis(2)*Axis(3)*(1-C)+Axis(1)*S, Axis(3)**2+(1-Axis(3)**2)*C/),(/3,3/))
	PtOut = matmul(R,PtIn)
end subroutine AxialRotate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!              Perform least-square approoximation              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LeastSquare (coeff, value, Ncoeff, Nval, sets)
    implicit none
    integer ::  Ncoeff, Nval, sets, i
    real(8) ::  coeff (Ncoeff, sets), value (Nval, Ncoeff+sets)
    real(8) ::  ATA (Ncoeff, Ncoeff), ATY (Ncoeff)
    
    ATA = matmul(transpose(value(:,1:Ncoeff)),value(:,1:Ncoeff))        !May waste a lot of calcullatons here...
    call InvMat(ATA, ATA, Ncoeff)
    do i = 1,sets
        ATY = matmul(transpose(value(:,1:Ncoeff)),value(:,Ncoeff+i))
        coeff(:,i) = matmul(ATA,ATY)
    end do
    
end subroutine LeastSquare
END MODULE MathKernel

