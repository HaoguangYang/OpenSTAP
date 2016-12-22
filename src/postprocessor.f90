!===========================================================================================
! 
! POST PROCESSOR SUBROUTINES INCLUDING:
!    SPR Stress Recovery
!    VTK File Generation (use with specific sentences in stap90, noted in comments.
! 
! Author:    Haoguang Yang
! 
!===========================================================================================


subroutine PostProcessor (ElementType, Dimen, PositionData, &
                          Node, NGauss, GaussianCollection, StressCollection, Displacement)
!===============================================================================
!Function:  Stress Recovery using SPR
!Author:    Haoguang Yang
!Note:      Connectivity matrix required as Node input here.
!===============================================================================
    use globals
    use mathkernel
    
    implicit none
    
    integer :: Dimen, NGauss, ElementType, ref1, ref2, NStress
    integer :: NodeRelationFlag(NUMNP,NPAR(5)*2+12), Ncoeff, Nval, N, i, j, k, L, &
               ind0, ind1, ind2, Node(NPAR(2),NPAR(5))
    real(8) :: coeff(10,6), Stress(6,NUMNP), PositionData(Dimen*NPAR(5), NPAR(2)), U(NEQ), &
               GaussianCollection(Dimen, NPAR(2)*NGauss), StressCollection(6,NPAR(2)*NGauss)
    real(8) :: x, y, z, Displacement(NEQ)
    REAL(8), ALLOCATABLE :: value(:,:)
    character(len=19) :: String
    
    write (IOUT,"(/,/)") 
    write (IOUT,"(A81,/,A99,/,A96)") "               S T R E S S   R E C O V E R Y   A T   N O D A L   P O I N T S", &
                " Node   |-----------------------------------------Stress----------------------------------------|", &
                "Number     Sigma_XX       Sigma_YY       Sigma_ZZ       Sigma_XY       Sigma_YZ       Sigma_ZX"
    NodeRelationFlag(:,:) = 0
    ref1 = NPAR(5)*2+11                                             !Set Node-Connection Counter Position Reference
    ref2 = NPAR(5)*2+12                                             !Set First-Node Counter Position Reference
    Stress(:,:) = 0D0
    DO N = 1, NPAR(2)                                               !NumberOfElements
        do i = 1,NPAR(5)
            j = NodeRelationFlag(Node(N,i),ref1) + 1                !How many elements connected to the node
            NodeRelationFlag(Node(N,i),ref1) = j
            if (j .EQ. 1) then                                      !Hint 1
                NodeRelationFlag(Node(N,i), ref2) = i               !The first element connected to the node with its node i
            end if
            NodeRelationFlag(Node(N,i),j) = N  
            !write (*,*) 'Node Relation Flag',N,i
            !write (*,*) NodeRelationFlag(Node(N,i),:)
        end do
    end do
    
    i = maxval(NodeRelationFlag(:,ref1))*NGauss                         !How many gaussian points will contribute to the SPR
    if (Dimen == 3) then
        NStress = 6
        allocate (value(i,16))
        do L =1, NUMNP
            coeff(:,:) = 0
            Nval = NodeRelationFlag(L,ref1) * NGauss
            if ( Nval .NE. 0) then
                ind0 = 1
                if (Nval .GE. 25) then                                  !Chooose whether to use quadratic or linear interplotation
                    Ncoeff = 7
                else if (Nval .GE. 6) then
                    Ncoeff = 4
                else
                    Ncoeff = 0
                end if
                !Ncoeff = 0
                do ind2 = 1, NodeRelationFlag(L,ref1)                   !Must Run Serial!
                    N = NodeRelationFlag (L, ind2)
                    ind1 = (N-1)*NGauss+1
                    do j = 1, NGauss
                        ind0 = (ind2-1)*NGauss+j
                        if (Ncoeff .GT. 0) then
                            x = GaussianCollection (1, ind1+mod(ind0-1,NGauss))
                            y = GaussianCollection (2, ind1+mod(ind0-1,NGauss))
                            z = GaussianCollection (3, ind1+mod(ind0-1,NGauss))
                        end if
                        
                        Stress(:,L) = StressCollection (:,ind1+mod(ind0-1,NGauss))
                    
                        select case (Ncoeff)
                        case (7)
                            value(ind0,1:Ncoeff+6) = reshape((/1D0, x, y, z, x*y, y*z, z*x, &                       !x**2, y**2, z**2, &
                                                                      Stress(1:6,L)/), (/Ncoeff+6/))
                        case (4)
                            value(ind0,1:Ncoeff+6) = reshape((/1D0, x, y, z, Stress(1:6,L)/), (/Ncoeff+6/))
                        case (0)
                            value(ind0,1:6) = Stress(1:6,L)
                        end select
                        !ind0 = ind0 + 1
                        !write(*,*) 'Nval= ',Nval,'-----------------------------------------------------'
                        !write (*,*)'x,y,z= ', x, y, z
                        !write (*,*)'Stress = ', Stress(:,L)
                    end do
                end do
            
                !sets = 6
                if (Ncoeff .GT. 0) THEN
                    call LeastSquare (coeff(1:Ncoeff,:), value(1:Nval,1:Ncoeff+6), Ncoeff, Nval, 6)
                    !write (*,*) coeff
                    ind2 = NodeRelationFlag(L,ref2)
                    x = PositionData(3*(ind2-1)+1,NodeRelationFlag(L,1))    !(L,1) relative to Hint 1
                    y = PositionData(3*(ind2-1)+2,NodeRelationFlag(L,1))
                    z = PositionData(3*(ind2-1)+3,NodeRelationFlag(L,1))
                    if (Ncoeff .EQ. 7) Stress(:,L) = matmul(transpose(coeff(1:7,:)),(/1D0, x, y, z, x*y, y*z, z*x/))              !, x**2, y**2, z**2/))
                    if (Ncoeff .EQ. 4) Stress(:,L) = matmul(transpose(coeff(1:4,:)),(/1D0, x, y, z/))
                else
                    Stress(:,L) = (sum(value(1:Nval, 1:6),DIM = 1))/NVal
                end if
                write (IOUT,"(I6, 3X, E13.6, 5(2X, E13.6))") L, Stress(1:6,L)
            end if
        end do
        
    else if (Dimen == 2) then
        if (ElementType>=6 .AND. ElementType<=9) then
            NStress = 5                                                 !Shell
        else
            NStress = 3
        end if
        allocate (value(i,6+NStress))
        do L =1, NUMNP
            coeff(:,:) = 0
            Nval = NodeRelationFlag(L,ref1) * NGauss
            if (Nval .NE. 0) then
                ind0 = 1
                if (Nval .GE. 12) then
                    Ncoeff = 4
                else if (Nval .GE. 4) then
                    Ncoeff = 3
                else
                    Ncoeff = 0
                end if
                !Ncoeff = 0
                do ind2 = 1, NodeRelationFlag(L,ref1)                   !Must Run Serial!
                    N = NodeRelationFlag (L, ind2)
                    ind1 = (N-1)*NGauss+1
                    do j = 1, NGauss
                        ind0 = (ind2-1)*NGauss+j
                        if (Ncoeff .GT. 0) then
                            x = GaussianCollection (1, ind1+mod(ind0-1,NGauss))
                            y = GaussianCollection (2, ind1+mod(ind0-1,NGauss))
                        end if
                                
                        Stress(1:NStress,L) = StressCollection (1:NStress,ind1+mod(ind0-1,NGauss))
                                
                        select case (Ncoeff)
                        case (4)
                            value(ind0,1:Ncoeff+NStress) = reshape((/1D0, x, y, x*y, &                          !x**2, y**2, &
                                                                      Stress(1:NStress,L)/), (/Ncoeff+NStress/))
                        case (3)
                            value(ind0,1:Ncoeff+NStress) = reshape((/1D0, x, y, Stress(1:NStress,L)/), (/Ncoeff+NStress/))
                        case (0)
                            value(ind0,1:NStress) = Stress(1:NStress,L)
                        end select
                        !ind0 = ind0 + 1
                        !write(*,*) 'Nval= ',Nval,'-----------------------------------------------------'
                        !write (*,*)'x,y,z= ', x, y, z
                        !write (*,*)'Stress = ', Stress(:,L)
                    end do
                end do
                !sets = 3
                if (Ncoeff .GT. 0) then
                    call LeastSquare (coeff(1:Ncoeff,:), value(1:Nval,1:Ncoeff+NStress), Ncoeff, Nval, NStress)
                    !write (*,*) coeff
                    !write (*,*) value
                    ind2 = NodeRelationFlag(L,ref2)
                    x = PositionData(2*(ind2-1)+1,NodeRelationFlag(L,1))    !(L,1) relative to Hint 1
                    y = PositionData(2*(ind2-1)+2,NodeRelationFlag(L,1))
                    if (Ncoeff .EQ. 4) Stress(1:NStress,L) = matmul(transpose(coeff(1:4,1:NStress)),(/1D0, x, y, x*y/))         !, x**2, y**2/))
                    if (Ncoeff .EQ. 3) Stress(1:NStress,L) = matmul(transpose(coeff(1:3,1:NStress)),(/1D0, x, y/))
                else
                    Stress(1:NStress,L) = (sum(value(1:Nval, 1:NStress),DIM = 1))/NVal
                end if
                if (NStress == 3) then
                    write (IOUT,"(I6, 3X, E13.6, 2X, E13.6, A13, 2X, E13.6, 2(2X, A13))") &
                                                                    L, Stress(1:2,L), "---", Stress(3,L), "---", "---"
                else if (NStress == 5) then
                    write (IOUT,"(I6, 3X, E13.6, 2X, E13.6, A13, 3(2X, E13.6))") &
                                                                    L, Stress(1:2,L), "---", Stress(3:5,L)
                end if
            end if
        end do
        
    else if (Dimen == 1) then                                           !Truss or beam
        if (ElementType ==1) then                                       !Truss
            Nstress = 1
        else if (ElementType == 5) then                                 !Beam
            NStress = 6
        end if
        allocate (value(i,6+NStress))
        do L =1, NUMNP
            coeff(:,:) = 0
            Nval = NodeRelationFlag(L,ref1) * NGauss
            if (Nval .NE. 0) then
                ind0 = 1
                do ind2 = 1, NodeRelationFlag(L,ref1)                   !Must Run Serial!
                    N = NodeRelationFlag (L, ind2)
                    ind1 = (N-1)*NGauss+1
                    do j = 1, NGauss
                        ind0 = (ind2-1)*NGauss+j
                        Stress(1:NStress,L) = StressCollection (1:NStress,ind1+mod(ind0-1,NGauss))
                        value(ind0,1:NStress) = Stress(1:NStress,L)
                        !ind0 = ind0 + 1
                        !write (*,*) x, y
                        !write (*,*) Stress(:,L)
                        !write(*,*) Nval
                    end do
                end do
                !sets = 1
                Stress(1:NStress,L) = (sum(value(1:Nval, 1:NStress),DIM = 1))/NVal
                if (NStress == 1) then
                    write (IOUT,"(I6, 3X, E13.6, 5(2X, A13))") L, Stress(1,L), "---", "---", "---", "---", "---"
                else 
                    write (IOUT,"(I6, 3X, E13.6, 5(2X, E13.6))") L, Stress(1:6,L)
                end if
            end if
        end do
    
    end if
    deallocate (value)
    NEL = NEL+NPAR(2)                                               !renew total number of elements.
    NCONECT = NCONECT + NPAR(2)*(NPAR(5)+1)                         !Renew total connectivity matrix element number
    
    write (String, "('Stress_Load_Case',I2.2)") CURLCASE
    write (VTKTmpFile) String, NStress, NUMNP
    do j = 1, NUMNP
        write (VTKTmpFile) Stress(1:NStress,j)                      !Stresses at each nodal point
    end do
    write (String, "('Mises_Load_Case',I2.2)") CURLCASE
    write (VTKTmpFile) String, 1, NUMNP
    do j = 1, NUMNP
        write (VTKTmpFile) sqrt(((Stress(1,j)-Stress(2,j))**2   &
                                +(Stress(2,j)-Stress(3,j))**2   &
                                +(Stress(3,j)-Stress(1,j))**2   &
                                +6*(dot_product(Stress(3:6,j),Stress(3:6,j))))/2)
        !von-Mises Stress at each nodal point
    end do
    
    
end subroutine PostProcessor



subroutine VTKgenerate (Flag)
!===============================================================================
!Function:  Generating VTK Files for ParaView Visualization
!Author:    Haoguang Yang
!Note:      Called at places written in comments.
!===============================================================================
    use globals
    use memallocate
    implicit none
    integer :: i,j,k,Flag, Dat(NEL+100)
    character(len=25) :: string
    real(8) :: Dat1(6), Stress(7,NUMNP)                             !Stress component at xx, yy, zz, xy, yz, zx, von-Mises
    
select case (Flag)
case (1)                                                            !Called in solution phase IND=1
    NCONECT = 0
    write (VTKFile,"(A26)") '# vtk DataFile Version 3.0'
    write (VTKFile,*) HED
    write (VTKFile,"(A5)") 'ASCII'
    write (VTKFile,"(A25)") 'DATASET UNSTRUCTURED_GRID'
    write (VTKFile,*) "POINTS",NUMNP,"double"
    do i = 1, NUMNP
        write (VTKFile,*) DA(NP(2)+i-1), DA(NP(3)+i-1), DA(NP(4)+i-1)           !X(i), Y(i), Z(i)
    end do
    
case (2)                                                            !Called in elcal.f90, subroutine ELCAL
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
     
case (3)                                                            !Called in stap.f90, STAP at solution phase IND=3
    write (VTKFile,*) "CELLS ", NEL, NCONECT                        !Sum up all elements to generate a global picture
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
    write (VTKFile,*) "FIELD Result ", NLCASE*3
    rewind (VTKTmpFile)
    do i = 1 , NLCASE
        read (VTKTmpFile) string(1:25), Dat(1:2)                    !Fetch Displacements of Load Cases
        write (VTKFile,*) string(1:25), Dat(1:2), "double"          !Displacements are written into temp files in subroutine writd, stap.f90
        do j = 1, NUMNP
            read (VTKTmpFile) Dat1(1:Dat(1))
            write (VTKFile,*) Dat1(1:Dat(1))
        end do
        do k = 1, NUMEG
            read (VTKTmpFile) string(1:19), Dat(1:2)                            !Fetch Stress of Load Cases
            if (k == 1) write (VTKFile,*) string(1:19), 6, Dat(2), "double"     !Dat(1) should be <=6 for no more than 6 stress components
            do j = 1, NUMNP
                read (VTKTmpFile) Dat1(1:Dat(1))
                if ((Stress(1,j) .NE. 0.) .OR.(Stress(2,j) .NE. 0.) .OR.(Stress(3,j) .NE. 0.) .OR. &
                    (Stress(4,j) .NE. 0.) .OR.(Stress(5,j) .NE. 0.) .OR.(Stress(6,j) .NE. 0.)) then
                    Stress(1:Dat(1),j) = 0.5*(Stress(1:Dat(1),j)+Dat1(1:Dat(1)))
                else
                    Stress(1:Dat(1),j) = Dat1(1:Dat(1))
                end if
            end do
            read (VTKTmpFile) string(1:18), Dat(1:2)                            !Fetch Mises Stress of Load Cases
                                                                                !Dat(1) should be 1
            do j = 1, NUMNP
                read (VTKTmpFile) Dat1(1)
                if (Stress(7,j) .NE. 0) then
                    Stress(7,j) = 0.5*(Stress(7,j)+Dat1(1))
                else
                    Stress(7,j) = Dat1(1)
                end if
            end do
        end do
        do j = 1,NUMNP
            do k = 1,6
                if (isnan(Stress(k,j))) Stress(k,j) = 0.
            end do
            write (VTKFile,*) Stress(1:6, j)
        end do
        write (VTKFile,*) string(1:18), 1, NUMNP, "double"
        write (VTKFile,*) Stress(7,1:NUMNP)
    end do
end select
    
    !do i = 1, NUMEG                                            ----|
    !    do i = 1,NPAR(2) !of ELEGP#i                               |
    !        write (VTKNodeTmp,*) NPAR(5), Nodes(i,:)               !This part called in each element case 1 before enddo.
    !    end do                                                     |
    !end do                                                     ----|
    
end subroutine VTKgenerate
