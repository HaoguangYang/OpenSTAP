subroutine PostProcessor (ElementType, Dimen, PositionData, &
                          Node, QuadratureOrder, GaussianCollection, StressCollection, Displacement)
!Stress Recovery using SPR
    use globals
    use mathkernel
    
    implicit none
    
    integer :: Dimen, QuadratureOrder, ElementType, NumberOfStress
    integer :: NodeRelationFlag (NUMNP,QuadratureOrder**Dimen+2), Ncoeff, Nval, N, i, j, k, L, &
               ind0, ind1, ind2, Node(NPAR(2),NPAR(5))
    real(8) :: coeff(10,6), value(NPAR(5)*QuadratureOrder**Dimen,16), Stress(6,NUMNP), PositionData(3*NPAR(5), NPAR(2)), U(NEQ), &
               GaussianCollection(3, NPAR(2)*QuadratureOrder**Dimen), StressCollection(6,NPAR(2)*QuadratureOrder**Dimen)
    real(8) :: x, y, z, Displacement(NEQ)
    character(len=19) :: String
    
    write (IOUT,"(/,/)") 
    write (IOUT,*) "               S T R E S S   R E C O V E R Y   A T   N O D A L   P O I N T S"
    write (IOUT,*) " Node   |-----------------------------------------Stress----------------------------------------|"
    write (IOUT,*) "Number     Sigma_XX       Sigma_YY       Sigma_ZZ       Sigma_XY       Sigma_YZ       Sigma_ZX"
    NodeRelationFlag(:,:) = 0
    DO N = 1, NPAR(2)                                                                       !NumberOfElements
        do i = 1,QuadratureOrder**Dimen
            j = NodeRelationFlag(Node(N,i),9) + 1                                           !How many elements connected to the node
            NodeRelationFlag(Node(N,i),9) = j
            if (j .EQ. 1) then                                                              !Hint 1
                NodeRelationFlag(Node(N,i), 10) = i                                         !The first element connected to the node with its node i
            end if
            NodeRelationFlag(Node(N,i),j) = N  
            !write (*,*) NodeRelationFlag(Node(N,i),:)
        end do
    end do
    
    if (Dimen == 3) then
        do L =1, NUMNP
            coeff(:,:) = 0
            Nval = NodeRelationFlag(L,9) * QuadratureOrder**3
            ind0 = 1
            if (Nval .GE. 18) then
                Ncoeff = 10
            else
                Ncoeff = 4
            end if
            do ind2 = 1, NodeRelationFlag(L,9)
                N = NodeRelationFlag (L, ind2)
                
                ind1 = (N-1)*QuadratureOrder**3+1
                
                do i = 1, QuadratureOrder
                    do j = 1, QuadratureOrder
                        do k = 1, QuadratureOrder
                            x = GaussianCollection (1, ind1+mod(ind0-1,8))
                            y = GaussianCollection (2, ind1+mod(ind0-1,8))
                            z = GaussianCollection (3, ind1+mod(ind0-1,8))
                            Stress(:,L) = StressCollection (:,ind1+mod(ind0-1,8))
                            
                            if (Ncoeff .EQ. 10) &
                                value(ind0,1:Ncoeff+6) = reshape((/1D0, x, y, z, x*y, y*z, z*x, x**2, y**2, z**2, &
                                                                  Stress(1:6,L)/), (/Ncoeff+6/))
                            if (Ncoeff .EQ. 4) &
                                value(ind0,1:Ncoeff+6) = reshape((/1D0, x, y, z, Stress(1:6,L)/), (/Ncoeff+6/))
                            ind0 = ind0 + 1
                            !write (*,*) x, y, z
                            !write (*,*) Stress(:,L)
                            !write(*,*) Nval
                        end do
                    end do
                end do
            end do
            
            !sets = 6
            call LeastSquare (coeff(1:Ncoeff,:), value(1:Nval,1:Ncoeff+6), Ncoeff, Nval, 6)
            !write (*,*) coeff
            ind2 = NodeRelationFlag(L,10)
            x = PositionData(3*(ind2-1)+1,NodeRelationFlag(L,1))                                             !(L,1) relative to Hint 1
            y = PositionData(3*(ind2-1)+2,NodeRelationFlag(L,1))
            z = PositionData(3*(ind2-1)+3,NodeRelationFlag(L,1))
            if (Ncoeff .EQ. 10) Stress(:,L) = matmul(transpose(coeff),(/1D0, x, y, z, x*y, y*z, z*x, x**2, y**2, z**2/))
            if (Ncoeff .EQ. 4) Stress(:,L) = matmul(transpose(coeff(1:4,:)),(/1D0, x, y, z/))
            write (IOUT,"(I6, 3X, E13.6, 5(2X, E13.6))") L, Stress(1:6,L)
        end do
        
    else if (Dimen == 2) then
        do L =1, NUMNP
            coeff(:,:) = 0
            Nval = NodeRelationFlag(L,9) * QuadratureOrder**2           !Change if using triangular elements
            ind0 = 1
            if (Nval .GE. 8) then
                Ncoeff = 6
            else
                Ncoeff = 3
            end if
            do ind2 = 1, NodeRelationFlag(L,9)
                N = NodeRelationFlag (L, ind2)
                ind1 = (N-1)*QuadratureOrder**2+1
                    do j = 1, QuadratureOrder
                        do k = 1, QuadratureOrder
                            x = GaussianCollection (1, ind1+mod(ind0-1,4))
                            y = GaussianCollection (2, ind1+mod(ind0-1,4))
                            
                            Stress(1:3,L) = StressCollection (1:3,ind1+mod(ind0-1,4))
                            
                            if (Ncoeff .EQ. 6) &
                                value(ind0,1:Ncoeff+3) = reshape((/1D0, x, y, x*y, x**2, y**2, &
                                                                  Stress(1:3,L)/), (/Ncoeff+3/))
                            if (Ncoeff .EQ. 3) &
                                value(ind0,1:Ncoeff+3) = reshape((/1D0, x, y, Stress(1:3,L)/), (/Ncoeff+3/))
                            ind0 = ind0 + 1
                            !write (*,*) x, y
                            !write (*,*) Stress(:,L)
                            !write(*,*) Nval
                        end do
                    end do
            end do
            !sets = 3
            call LeastSquare (coeff(1:Ncoeff,:), value(1:Nval,1:Ncoeff+3), Ncoeff, Nval, 3)
            !write (*,*) coeff
            ind2 = NodeRelationFlag(L,10)
            x = PositionData(2*(ind2-1)+1,NodeRelationFlag(L,1))                                             !(L,1) relative to Hint 1
            y = PositionData(2*(ind2-1)+2,NodeRelationFlag(L,1))
            if (Ncoeff .EQ. 6) Stress(1:3,L) = matmul(transpose(coeff(1:6,1:3)),(/1D0, x, y, x*y, x**2, y**2/))
            if (Ncoeff .EQ. 3) Stress(1:3,L) = matmul(transpose(coeff(1:3,1:3)),(/1D0, x, y/))
            write (IOUT,"(I6, 3X, E13.6, 5(2X, E13.6))") L, Stress(1:2,L), 0, Stress(3,L), 0, 0
        end do
    end if
    NEL = NEL+NPAR(2)                               !renew total number of elements.
    NCONECT = NCONECT + NPAR(5)*(NPAR(2)+1)             !Renew total connectivity matrix element number
    
    if (Dimen == 3) NumberOfStress = 6
    if (Dimen == 2) NumberOfStress = 3
    write (String, "('Stress_Load_Case',I2.2)") CURLCASE
    write (VTKTmpFile) String, NumberOfStress, NUMNP
    do j = 1, NUMNP
        write (VTKTmpFile) Stress(1:NumberOfStress,j) !Stresses
    end do
    
end subroutine PostProcessor



subroutine VTKgenerate (Flag)
    use globals
    use memallocate
    implicit none
    integer :: i,j,Flag, Dat(NEL+100)
    character(len=25) :: string
    real(8) :: Dat1(NUMNP)
    
select case (Flag)
case (1)
    write (VTKFile,"(A26)") '# vtk DataFile Version 3.0'
    write (VTKFile,*) HED
    write (VTKFile,"(A5)") 'ASCII'
    write (VTKFile,"(A25)") 'DATASET UNSTRUCTURED_GRID'
    write (VTKFile,*) "POINTS",NUMNP,"double"
    do i = 1, NUMNP
        write (VTKFile,*) DA(NP(2)+i-1), DA(NP(3)+i-1), DA(NP(4)+i-1)   !X(i), Y(i), Z(i)
    end do
    
case (2)
    select case (NPAR(1))
     case (1)
        write (VTKElTypTmp) (3,I=1,NPAR(2))
     case (2)
        write (VTKElTypTmp) (9,I=1,NPAR(2))
     case (3)
        write (VTKElTypTmp) (28,I=1,NPAR(2))
     case (4)
        write (VTKElTypTmp) (12,I=1,NPAR(2))
     case (5)
        write (VTKElTypTmp) (3,I=1,NPAR(2))
     case (6)
        write (VTKElTypTmp) (9,I=1,NPAR(2))
     case (7)
        write (VTKElTypTmp) (9,I=1,NPAR(2))
     case (8)
        write (VTKElTypTmp) (23,I=1,NPAR(2))
     case (9)
        write (VTKElTypTmp) (23,I=1,NPAR(2))
     case (10)
        write (VTKElTypTmp) (5,I=1,NPAR(2))
     end select
     
case (3)
    write (VTKFile,*) "CELLS ", NEL, NCONECT          !Sum up all elements to generate a global picture
    rewind (VTKNodeTmp)
    do i = 1 , NEL
        read (VTKNodeTmp) Dat(1), Dat(2:1+Dat(1))
        write (VTKFile,*) Dat(1:1+Dat(1))
    end do
    write (VTKFile,*) "CELL_TYPES ", NEL
    rewind (VTKElTypTmp)
    read (VTKElTypTmp) Dat(1:NEL)
    write (VTKFile,*) Dat(1:NEL)
    write (VTKFile,*) "CELL_DATA", NEL
    write (VTKFile,*) "POINT_DATA", NUMNP
    write (VTKFile,*) "FIELD Result ", NLCASE*2
    rewind (VTKTmpFile)
    do i = 1 , NLCASE
        read (VTKTmpFile) string(1:25), Dat(1:2)                !Fetch Displacements of Load Cases
        write (VTKFile,*) string(1:25), Dat(1:2), "double"
        do j = 1, NUMNP
            read (VTKTmpFile) Dat1(1:Dat(1))
            write (VTKFile,*) Dat1(1:Dat(1))
        end do
        read (VTKTmpFile) string(1:19), Dat(1:2)                !Fetch Stress of Load Cases
        write (VTKFile,*) string(1:19), Dat(1:2), "double"
        do j = 1, NUMNP
            read (VTKTmpFile) Dat1(1:Dat(1))
            write (VTKFile,*) Dat1(1:Dat(1))
        end do
    end do
end select
    
    !do i = 1, NUMEG
    !    do i = 1,NPAR(2) !of ELEGP#i
    !        write (VTKNodeTmp,*) NPAR(5), Nodes(i,:)   !In each element case 1 before enddo.
    !    end do
    !end do
    
    
    !do i = 1, NUMEG
    !    select case (NPAR(1))
    !    case (1) write(VTKElTypTmp,*) ?
    !    case (2) write(VTKElTypTmp,*) ?
    !    case (3) write(VTKElTypTmp,*) ?    !In ELCAL.f90, subroutine Elcal
    !    end select
    !end do
    
    
end subroutine VTKgenerate
