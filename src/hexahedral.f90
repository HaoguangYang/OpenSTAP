subroutine hexahedral
    use globals
    use memallocate
    
    implicit none
    integer :: NumberOfElements, NumberOfMaterials, ElementGroupSize, QuadratureOrder=2
    integer :: N(11) !Pointers
    
    NPAR(5) = 8
    NPAR(6) = 8
    
    NumberOfElements = NPAR(2)
    NumberOfMaterials = NPAR(3)
    
    !Preallocate Memory
    
!pointer lists
! Calculate the pointer to the arrays in the element group data
! N101: E(NumberOfMaterials)
! N102: v(NumberOfMaterials)
! N103: Density(NumberOfMaterials)
! N104: Gravity(NumberOfMaterials)
! N105: LM(3*NPAR(5),NumberOfElements)
! N106: PositionData(3*NPAR(5),NumberOfElements)
! N107: MaterialData(NumberOfElements)
! N108: BMat (NumberOfElements*18*NPAR(5)*8)
! N109: Jacobian(NumberOfElements*QuadratureOrder^3)

    N(1) = 0
    N(2) = N(1)+NumberOfMaterials*ITWO
    N(3) = N(2)+NumberOfMaterials*ITWO
    
    if (NPAR(4) .GT. 0) then
        N(4) = N(3)+NumberOfMaterials*ITWO
        N(5) = N(4)+NumberOfMaterials*ITWO
    else
        N(4) = N(3)
        N(5) = N(4)
    end if
    
    N(6) = N(5) + 3*NPAR(5)*NumberOfElements
    N(7) = N(6) + 3*NPAR(5)*NumberOfElements*ITWO
    N(8) = N(7) + NumberOfElements
    N(9) = N(8) + 18*NPAR(5)*NumberOfElements*QuadratureOrder**3*2
    N(10)= N(9) + NumberOfElements*2*QuadratureOrder**3
    N(11)= N(10)+ NPAR(5)*NPAR(2)
    
    MIDEST = N(11)
    
    if (IND .EQ. 1) then
        call MemAlloc(11,"ELEGP",MIDEST,1)
    end if
    NFIRST = NP(11)
    N(:) = N(:) + NFIRST
    NLAST  = N(11)
    
    if ((NPAR(5) .EQ. 8) .AND. (NPAR(6) .EQ. 8)) call &
        HexEight (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
                  A(N(1)),A(N(2)),A(N(3)),A(N(4)),A(N(5)),A(N(6)),A(N(7)),A(N(8)),A(N(9)),A(N(10)))
    
    !Reuse DA(NP(4)) at Solution Phase 3 as displacement U    
    return
end subroutine hexahedral

subroutine HexEight (ID,X,Y,Z,U,MHT,E, PoissonRatio, Density, Gravity, LM, PositionData, MaterialData, BMat, Jacobian, Node)
    use globals
    use MemAllocate
    use MathKernel
    
    implicit none
    integer ::  ElementShapeNodes, NumberOfMaterials, NumberOfElements
    integer ::  QuadratureOrder     = 2                     ! NPAR(6) -- Element Load Nodes
    INTEGER ::  ID(3,NUMNP), MHT(NEQ), MaterialData(NPAR(2))
    integer ::  MaterialType, MaterialComp, ND, L, N, i, j, LM(24,NPAR(2)), ElementType, ind0, iprint, k, &
                NodeRelationFlag(NUMNP,2**3+2), ind1, ind2, Ncoeff, Nval, Node(NPAR(2),NPAR(5))
    real(8) ::  X(NUMNP), Y(NUMNP), Z(NUMNP), U(NEQ), &
                DetJ(2,2,2), E(NPAR(3)), PoissonRatio(NPAR(3)), BMat(NPAR(2)*18*NPAR(5)*2**3), ElementDisp(24)
    real(8) ::  BMatrix(6, 3*NPAR(5)), PositionData(3*NPAR(5), NPAR(2)), DMatrix(6,6), x1, y1, z1, &
                Transformed(3), W(2), Weight(2,2,2), GaussianPts(2), Jacobian(NPAR(2)*2**3), coeff(10,6), &
                value(8*2**3,16)
    real(8) ::  Young, v, S(3*NPAR(5),3*NPAR(5)), GaussianPtsPosit(3,2**3), Strain(6,2**3), Stress(6,2**3), &
                Density, Gravity, NMatrix(3,3*NPAR(5)), NormalVec(3), Point(3*NPAR(5),3*NPAR(5))
                
    ElementType         = NPAR(1)
    NumberOfElements    = NPAR(2)
    NumberOfMaterials   = NPAR(3)
    ElementShapeNodes   = NPAR(5)                           !NPAR(5)=8
    ND = 24
    
    SELECT CASE (IND)
    CASE(1)
        WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.1, TRUSS ELEMENTS',/,   &
                   '     EQ.2, 4Q ELEMENTS',/,      &
                   '     EQ.3, 9Q ELEMENTS',//,     &
                   '     EQ.4, 8H ELEMENTS',//,     &
                   '     EQ.5, 3T ELEMENTS',//,     & 
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") ElementType, NumberOfElements
        IF (NumberOfMaterials.EQ.0) NumberOfMaterials=1
        WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I5,/)") NumberOfMaterials
        
        WRITE (IOUT,"('  SET       YOUNG''S    POISSON',/,  & 
                            ' NUMBER     MODULUS      RATIO',/,  &
                   15 X,'E', 12 X, 'v')")

        DO I=1,NumberOfMaterials
            READ (IIN,'(I5,2F10.0)') N,E(N), PoissonRatio(N)      ! Read material information for 3D Homogeneous
            WRITE (IOUT,"(I5,4X,E12.5,2X,E14.6)") N,E(N), PoissonRatio(N)
        END DO
        WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT        |------------------------- NODES -------------------------|       MATERIAL',/,   &
                      ' NUMBER-N        1       2       3       4       5       6       7       8       SET NUMBER')")
        N=0
        
        CALL GaussianMask(GaussianPts, W, QuadratureOrder)
        
        DO WHILE (N .NE. NumberOfElements)
            READ (IIN,'(11I5)') N,Node(N,1:8),MaterialType          ! Read in element information
            
    !       Save element information
            PositionData(1:ElementShapeNodes*3-1:3,N)=X(Node(N,:))        ! Coordinates of the element's nodes
            PositionData(2:ElementShapeNodes*3  :3,N)=Y(Node(N,:))
            PositionData(3:ElementShapeNodes*3+1:3,N)=Z(Node(N,:))
            
            MaterialData(N) = MaterialType                              ! Material type

            DO L=1,ND
                LM(L,N)=0
            END DO
            DO L=1,3
                LM(L:ElementShapeNodes*3+L-2:3,N) = ID(L,Node(N,:))       ! Connectivity matrix
            END DO
            DetJ(:, :, :)=0
            ind0 = 1
            do i = 1, QuadratureOrder
                do j = 1, QuadratureOrder
                    do k = 1, QuadratureOrder
                        Transformed = (/GaussianPts(i),GaussianPts(j),GaussianPts(k)/) 
                        CALL HexB(BMatrix, DetJ(i,j,k), ElementShapeNodes, Transformed, &
                                  (/X(Node(N,:)), Y(Node(N,:)), Z(Node(N,:))/)) !Needs reshaping
                        BMat(((n-1)*QuadratureOrder**3+ind0-1)*18*ElementShapeNodes+1 : &
                             ((n-1)*QuadratureOrder**3+ind0)*18*ElementShapeNodes) = reshape(BMatrix, (/18*ElementShapeNodes/))
                        ind0 = ind0 +1
                        
                        !write (*,*) ((n-1)*QuadratureOrder**3+ind0-2)*18*ElementShapeNodes+1,'----------------------',& 
                        !            ((n-1)*QuadratureOrder**3+ind0-1)*18*ElementShapeNodes
                        !write (*,*) BMatrix
                        
                    end do
                end do
            end do
            Jacobian((n-1)*QuadratureOrder**3+1 : n*QuadratureOrder**3)=reshape(DetJ, (/QuadratureOrder**3/))
            
            !write (*,*) "DetJ",(n-1)*QuadratureOrder**3+1,'---',n*QuadratureOrder**3, DetJ
            !BMat((n-1)*6*ElementShapeNodes+1 : n*6*ElementShapeNodes)=reshape(BMatrix, (/6*ElementShapeNodes/))
            
            CALL COLHT (MHT,ND,LM(:,N))
            WRITE (IOUT,"(I7,5X,7(I7,1X),I7,4X,I5)") N,Node(N,1:ElementShapeNodes),MaterialType
            
            !write (IOUT,*) 'MHT',MHT
            !write (IOUT,*) 'NormalVec',NormalVec
            !write (IOUT,*) 'Transformed',Transformed
            !write (IOUT,*) 'Flattened',Flattened
            !write (IOUT,*) 'BMat',BMat
            !write (IOUT,*) 'DetJ',Jacobian
            
        enddo
        return
    
    CASE (2)
        CALL GaussianMask(GaussianPts, W, QuadratureOrder)
        do i = 1, QuadratureOrder
            do j = 1, QuadratureOrder
                do k = 1, QuadratureOrder
                    Weight(i,j,k)  = W(i)*W(j)*W(k)
                end do
            end do
        end do
        
        MaterialComp = -1
        do n = 1, NumberOfElements
            S(:,:) = 0
            MaterialType = MaterialData(N)
            if (MaterialType .NE. MaterialComp) then
                Young        = E(MaterialType)
                v            = PoissonRatio(MaterialType)
                call GetDMat(Young, v, DMatrix)
                
                MaterialComp = MaterialType
            end if
            
            DetJ = reshape(Jacobian((n-1)*QuadratureOrder**3+1 : n*QuadratureOrder**3), &
                           (/QuadratureOrder, QuadratureOrder, QuadratureOrder/))
            ind0  = 1
            do i = 1, QuadratureOrder
                do j = 1, QuadratureOrder
                    do k = 1, QuadratureOrder
                        BMatrix = reshape(BMat(((n-1)*QuadratureOrder**3+ind0-1)*18*ElementShapeNodes+1 : &
                                          ((n-1)*QuadratureOrder**3+ind0)*18*ElementShapeNodes),(/6, 3*ElementShapeNodes/))
                        
                        !write (*,*) ((n-1)*QuadratureOrder**3+ind0-2)*18*ElementShapeNodes+1,'----------------------',& 
                        !            ((n-1)*QuadratureOrder**3+ind0-1)*18*ElementShapeNodes
                        !write (*,*) BMatrix
                        
                        Point   = matmul(matmul(transpose(BMatrix),DMatrix),BMatrix)
                        S = S + (Weight(i,j,k)*DetJ(i,j,k))*Point
                        ind0 = ind0 + 1
                    end do
                end do
            end do
            
            !write(*,*) "S",S
            
            CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(:,N),ND)
        end do
        
    CASE (3)
        
        IPRINT=0
        call GaussianMask(GaussianPts, W, QuadratureOrder)
        DO N=1,NumberOfElements
            IPRINT=IPRINT + 1
            IF (IPRINT.GT.50) IPRINT=1
            IF (IPRINT.EQ.1) WRITE (IOUT,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  ',  &
                                           'E L E M E N T  G R O U P',I4,//,   &
            '  ELEMENT',13X, 'COORDINSTES',19X, 'Sigma_xx',9X, 'Sigma_yy',9X,'Sigma_zz',9X, &
            'Sigma_xy',9X,'Sigma_yz',9X,'Sigma_xz', /, &
            '  NUMBER', 8X,'X',10X,'Y',10X,'Z')") NG
            MaterialType = MaterialData(N)
            Young        = E(MaterialType)
            v            = PoissonRatio(MaterialType)
            call GetDMat(Young, v, DMatrix)
            
            ElementDisp(:) = 0
            
            do i = 1,ND
                if (LM(i,N) .NE. 0) ElementDisp(i)  =  U(LM(i,N))
            end do
            
            ind0 = 1
            do i = 1,QuadratureOrder
                do j = 1,QuadratureOrder
                    do k = 1,QuadratureOrder
                        Transformed   = (/GaussianPts(i), GaussianPts(j), GaussianPts(k)/)
                        call HexN (NMatrix, ElementShapeNodes, Transformed)
                        GaussianPtsPosit(:,ind0) = matmul(reshape(PositionData(:,N), (/3,ElementShapeNodes/)), &
                                                         NMatrix(1, 1:3*ElementShapeNodes:3))
                        ind1 = ((n-1)*QuadratureOrder**3+ind0-1)*18*ElementShapeNodes+1
                        ind2 = ((n-1)*QuadratureOrder**3+ind0)*18*ElementShapeNodes
                        BMatrix = reshape(BMat(ind1 : ind2),(/6, 3*ElementShapeNodes/))
                        Strain(:,ind0) = matmul(BMatrix, ElementDisp)             
                        Stress(:,ind0) = matmul(DMatrix, Strain(:,ind0))
                        
                        !Change the temporatory BMatrix Array to Gaussian Pts Array and Global Stress Array.
                        BMat (ind1+9*(ind0-1) : ind1+9*(ind0-1)+2) = GaussianPtsPosit(:,ind0)
                        BMat (ind1+9*(ind0-1)+3 : ind1+9*(ind0)-1) = Stress(:,ind0)
                        
                        ind0 = ind0 + 1
                    end do
                end do
            end do
            
            write (IOUT,"(I6,3(3X, F10.4),6(4X, E13.6),/,7(6X, 3(3X, F10.4),6(4X, E13.6),/))") &
                               N, (GaussianPtsPosit(:,I), Stress(:,I), I=1,QuadratureOrder**3)
        END DO

        !Stress Recovery using SPR
        write (IOUT,"(/,/)") 
        write (IOUT,*) "               S T R E S S   R E C O V E R Y   A T   N O D A L   P O I N T S"
        write (IOUT,*) " Node   |-----------------------------------------Stress----------------------------------------|"
        write (IOUT,*) "Number     Sigma_XX       Sigma_YY       Sigma_ZZ       Sigma_XY       Sigma_YZ       Sigma_ZX"
        NodeRelationFlag(:,:) = 0
        DO N = 1, NumberOfElements
            do i = 1,8
                j = NodeRelationFlag(Node(N,i),9) + 1
                NodeRelationFlag(Node(N,i),9) = j
                if (j .EQ. 1) then                                                              !Hint 1
                    NodeRelationFlag(Node(N,i), 10) = i
                end if
                NodeRelationFlag(Node(N,i),j) = N         
            end do
        end do
        do L =1, NUMNP
            coeff(:,:) = 0
            Nval = NodeRelationFlag(L,9) * 2**3
            ind0 = 1
            if (NodeRelationFlag(L,9) .GE. 4) then
                Ncoeff = 10
            else
                Ncoeff = 4
            end if
            do ind1 = 1, NodeRelationFlag(L,9)
                N = NodeRelationFlag (L, ind1)
                do i = 1, QuadratureOrder
                    do j = 1, QuadratureOrder
                        do k = 1, QuadratureOrder
                            ind2 = ((n-1)*QuadratureOrder**3+mod(ind0-1,8))*18*ElementShapeNodes+1
                            x1 = BMat (ind2+9*mod(ind0-1,8))
                            y1 = BMat (ind2+9*mod(ind0-1,8)+1)
                            z1 = BMat (ind2+9*mod(ind0-1,8)+2)
                            Stress(:,1) = BMat (ind2+9*mod(ind0-1,8)+3 : ind2+9*(mod(ind0-1,8)+1)-1)
                            if (Ncoeff .EQ. 10) &
                                value(ind0,1:Ncoeff+6) = reshape((/1D0, x1, y1, z1, x1*y1, y1*z1, z1*x1, x1**2, y1**2, z1**2, &
                                                                  Stress(1:6,1)/), (/Ncoeff+6/))
                            if (Ncoeff .EQ. 4) &
                                value(ind0,1:Ncoeff+6) = reshape((/1D0, x1, y1, z1, Stress(1:6,1)/), (/Ncoeff+6/))
                            ind0 = ind0 + 1
                            !write (*,*) x1, y1, z1
                            !write (*,*) Stress(:,1)
                            !write(*,*) Nval
                        end do
                    end do
                end do
            end do
            
            !sets = 6
            call LeastSquare (coeff(1:Ncoeff,:), value(1:Nval,1:Ncoeff+6), Ncoeff, Nval, 6)
            !write (*,*) coeff
            ind2 = NodeRelationFlag(L,10)
            x1 = PositionData(3*(ind2-1)+1,NodeRelationFlag(L,1))                                             !(L,1) relative to Hint 1
            y1 = PositionData(3*(ind2-1)+2,NodeRelationFlag(L,1))
            z1 = PositionData(3*(ind2-1)+3,NodeRelationFlag(L,1))
            if (Ncoeff .EQ. 10) Stress(:,1) = matmul(transpose(coeff),(/1D0, x1, y1, z1, x1*y1, y1*z1, z1*x1, x1**2, y1**2, z1**2/))
            if (Ncoeff .EQ. 4) Stress(:,1) = matmul(transpose(coeff(1:4,:)),(/1D0, x1, y1, z1/))
            write (IOUT,"(I6, 3X, E13.6, 5(2X, E13.6))") L, Stress(1:6,1)
        end do
        
    END SELECT

end subroutine HexEight


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!              Calculate Shape Function Matrix             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HexN (NMatrix, ElementShapeNodes, Transformed)
    implicit none
    integer ::  ElementShapeNodes, i
    real(8) ::  NMatrix(3, 3*ElementShapeNodes), Transformed(3), N(ElementShapeNodes)
    
    select case (ElementShapeNodes)
    case (8)
        N(1) = (1-Transformed(1))*(1-Transformed(2))*(1-Transformed(3))/8
        N(2) = (1+Transformed(1))*(1-Transformed(2))*(1-Transformed(3))/8
        N(3) = (1+Transformed(1))*(1+Transformed(2))*(1-Transformed(3))/8
        N(4) = (1-Transformed(1))*(1+Transformed(2))*(1-Transformed(3))/8
        N(5) = (1-Transformed(1))*(1-Transformed(2))*(1+Transformed(3))/8
        N(6) = (1+Transformed(1))*(1-Transformed(2))*(1+Transformed(3))/8
        N(7) = (1+Transformed(1))*(1+Transformed(2))*(1+Transformed(3))/8
        N(8) = (1-Transformed(1))*(1+Transformed(2))*(1+Transformed(3))/8
        NMatrix = reshape((/((/N(i), 0D0, 0D0, 0D0, N(i), 0D0, 0D0, 0D0, N(i)/),i=1,ElementShapeNodes)/), &
                          shape(NMatrix))
    end select
end subroutine HexN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       Calculate Gradiient of Shape Function Matrix       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Still Problematic.
subroutine HexB (BMatrix, DetJ, ElementShapeNodes, Transformed, Original)
    USE MathKernel
    implicit none
    integer, intent(in) ::  ElementShapeNodes
    real(8) ::  BMatrix(6, 3*ElementShapeNodes), DetJ, Transformed(3), Original(ElementShapeNodes, 3)
    real(8) ::  GradN(3,ElementShapeNodes), J(3,3), InvMatJ(3,3), DerivN(3,ElementShapeNodes)
    real(8) ::  Bx(ElementShapeNodes), By(ElementShapeNodes), Bz(ElementShapeNodes), xi, eta, zta
    integer ::  i
    logical ::  OK_Flag
    
    select case (ElementShapeNodes)
    case (8)
        xi      = Transformed(1)
        eta     = Transformed(2)
        zta     = Transformed(3)
        GradN   = 0.125*reshape((/-(1-eta)*(1-zta), -(1-xi)*(1-zta), -(1-xi)*(1-eta), &
                                   (1-eta)*(1-zta), -(1+xi)*(1-zta), -(1+xi)*(1-eta), &
                                   (1+eta)*(1-zta),  (1+xi)*(1-zta), -(1+xi)*(1+eta), &
                                  -(1+eta)*(1-zta),  (1-xi)*(1-zta), -(1-xi)*(1+eta), &
                                  -(1-eta)*(1+zta), -(1-xi)*(1+zta),  (1-xi)*(1-eta), &
                                   (1-eta)*(1+zta), -(1+xi)*(1+zta),  (1+xi)*(1-eta), &
                                   (1+eta)*(1+zta),  (1+xi)*(1+zta),  (1+xi)*(1+eta), &
                                  -(1+eta)*(1+zta),  (1-xi)*(1+zta),  (1-xi)*(1+eta)/), &
                                shape(GradN))
        J       = matmul(GradN, Original)
        DetJ    = Det(J, 3)
        call InvMat3(J, InvMatJ, OK_Flag)
        if (OK_Flag .EQV. .FALSE.) STOP "***ERROR*** Derivative Of Shape Function Is SINGULAR"
        DerivN  = matmul(InvMatJ,GradN)
        Bx(:)   = DerivN(1,:)
        By(:)   = DerivN(2,:)
        Bz(:)   = DerivN(3,:)
        BMatrix = reshape((/((/Bx(i), 0D0, 0D0, By(i), 0D0, Bz(i), &
                             0D0, By(i), 0D0, Bx(i), Bz(i), 0D0, &
                             0D0, 0D0, Bz(i), 0D0, By(i), Bx(i)/), &
                             i = 1, ElementShapeNodes)/), Shape(BMatrix))
    end select
end subroutine HexB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                Calculate Stiffness Matrix                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetDMat (Young, Poisson, DMatrix)
    implicit none
    real(8),dimension(6,6)  ::  DMatrix
    real(8), intent(in)     ::  Young, Poisson
    real(8)                 ::  G, lambda
    integer                 ::  i
    
    DMatrix(:,:) = 0
    G = Young/(2*(1+Poisson))
    lambda = Poisson*Young/(1+Poisson)/(1-2*Poisson)
    DMatrix(1:3,1:3) = lambda
    do i = 1,6
        DMatrix(i,i) = DMatrix(i,i) + 2*G
    end do
    return
end subroutine GetDMat
