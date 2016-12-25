subroutine EIGENVAL(Stiff, Mass, MAXA, NN, NWK2, NNM, NRoot)
    USE GLOBALS, ONLY : IOUT, PARDISODOOR, NWK
    use memallocate
    implicit none
    real(8) :: Stiff(NWK2), Mass(NWK2)
    integer :: MAXA(NNM)            !NNM = NN+1
    integer :: NN, NWK2, NNM, NRoot
    real(8), parameter :: RTol = 1.0D-6
    real(8), allocatable :: EignVec(:,:),EignVal(:)
    !LANCZOS variables
    integer, parameter :: StiffTmp = 16
    integer :: NC, NNC, NRestart
    real(8) :: TT(NN),W(NN)
    REAL(8),allocatable :: VEC(:,:),D(:),RTOLV(:),BUP(:),BLO(:)
    REAL(8),allocatable :: BUPC(:),AR(:),BR(:),Q(:,:)
    logical :: IFSS, IFPR
    !dfeast_scsrgv variables
    character :: uplo
    integer :: fpm(128), info, loop
    real(8) :: epsout, emin, emax
    real(8),allocatable :: res(:)
    
    write(IOUT,*)'-------------------------------------------------------------------------------------'
    write(IOUT,*)'           E I G E N   V A L U E   C A L C U L A T I O N   R E S U L T S'
    
    NC = minval((/2*NRoot, NRoot+8, NWK2/))
    allocate (EignVec(NN,NC),EignVal(NC))
    if (.NOT. PARDISODOOR) then
        open(StiffTmp, FILE = "Stff.tmp", FORM = "UNFORMATTED", STATUS = "Replace")
        REWIND StiffTmp
        WRITE (StiffTmp) Stiff(1:NWK2)
        NNC=NC*(NC+1)/2
        allocate (AR(NNC),BR(NNC),VEC(NC,NC),D(NC),RTOLV(NC),BUP(NC),BLO(NC),BUPC(NC),Q(NN,NC))
        TT = 0
        W = 0
        AR = 0
        BR = 0
        VEC = 0
        D = 0
        RTOLV = 0
        BUP = 0
        BLO = 0
        BUPC = 0
        Q = 0
        EignVec = 0
        EignVal = 0
        !THE PARAMETERS NC AND/OR NRestart MUST BE INCREASED IF A SOLUTION HAS NOT CONVERGED
        NRestart = NNC
        IFSS = .TRUE.
        IFPR = .FALSE.
        write (*,*) Stiff
        write (*,*) Mass
        write (*,*) MAXA
        !write(*,*) EignVec, EignVal,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,BUPC, NN, NNM, NWK2, NWK2, NRoot, RTol, NC, NNC, NRestart, IFSS, IFPR, StiffTmp, IOUT
        call LANCZOS(Stiff, Mass, MAXA, EignVec, EignVal,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,BUPC, NN, NNM, NWK2, NWK2, NRoot, RTol, NC, NNC, NRestart, IFSS, IFPR, StiffTmp, IOUT,Q)
        REWIND StiffTmp
        READ (StiffTmp) Stiff
        close(StiffTmp)
        deallocate (AR,BR,VEC,D,RTOLV,BUP,BLO,BUPC,Q)
    else
        !uplo = 'U'
        allocate (res(NC))
        !IA(NP(9))      --  Mass Row Index
        !IA(NP(8))      --  Mass Column Indicator
        !IA(NP(5))      --  Stiffness Column Indicator
        !IA(NP(2))      --  Stiffness Row Index
        !DA(NP(10))     --  Mass Matrix
        !DA(NP(3))      --  Stiffness Matrix
        emin = 1.
        emax = 400.
        call feastinit (fpm)
        call pardiso_crop(DA(NP(10)), IA(NP(9)), IA(NP(8)))
        NWK2 = NWK                                          !Renew NWK
        call pardiso_crop(DA(NP(3)), IA(NP(9)), IA(NP(8)))
        call dfeast_scsrgv('u', NN, Stiff(1:NWK), IA(NP(2)), IA(NP(5)), Mass(1:NWK2), IA(NP(9)), IA(NP(8)), fpm, epsout, loop, &
                           emin, emax, NC, EignVal, EignVec, NRoot, res, info)
        write(*,*) "Eigen Values:",EignVal
        write(*,*) "Eigen Vectors:",EignVec
        write(IOUT,*) "Eigen Values:",EignVal
        write(IOUT,*) "Eigen Vectors:",EignVec
        write(IOUT,*) "Residual:",res
        deallocate (res)
    end if
    write(IOUT,*)'-------------------------------------------------------------------------------------'
    deallocate (EignVec, EignVal)
end subroutine EIGENVAL
    
subroutine prepare_MassMatrix
    use globals
    use memallocate
    implicit none

    CALL MEMALLOC(10,"M    ",NWK,ITWO)
    if (PARDISODOOR) then
        call memalloc(9,"MrInd",NEQ+1,1)
        call memalloc(8,"Mcolm",NWK,1)
        IA(NP(8):NP(8)+NWK) = IA(NP(5):NP(5)+NWK)
        IA(NP(9):NP(9)+NEQ+1) = IA(NP(2):NP(2)+NEQ+1)
    end if
end subroutine prepare_MassMatrix

