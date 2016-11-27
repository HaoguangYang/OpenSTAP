MODULE MathKernel
    
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
				  Axis(1)*Axis(3)*(1-C)-Axis(2)*S, Axis(2)*Axis(3)*(1-C)+Axis(1)*S, Axis(3)**2+(1-Axis(3)**2)*C/), shape(R))
	PtOut = matmul(R,PtIn)
end subroutine AxialRotate


END MODULE MathKernel

