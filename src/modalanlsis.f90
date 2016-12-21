subroutine EIGENVAL(Stiff, Mass, MAXA, NN, NWK, NNM, NRoot)
    USE GLOBALS, ONLY : IOUT, PARDISODOOR
    implicit none
    real(8) :: Stiff(NWK), Mass(NWK)
    integer :: MAXA(NNM)            !NNM = NN+1
    integer :: NN, NWK, NNM, NRoot
    real(8), parameter :: RTol = 1.0D-6
    real(8), allocatable :: EignVec(:,:),EignVal(:)
    !LANCZOS variables
    integer, parameter :: StiffTmp = 16
    integer :: NC, NNC, NRestart
    logical :: IFSS, IFPR
    !dfeast_scsrgv variables
    character(1) :: uplo
    integer :: fpm(128), info, loop
    real(8) :: epsout
    real(8),allocatable :: res(:)
    
    write(IOUT,*)'-------------------------------------------------------------------------------------'
    write(IOUT,*)'           E I G E N   V A L U E   C A L C U L A T I O N   R E S U L T S'
    
    NC = minval((/2*NRoot, NRoot+8, NWK/))
    allocate (EignVec(NN,NC),EignVal(NC))
    if (.NOT. PARDISODOOR) then
        open(StiffTmp, FILE = "Stff.tmp", FORM = "UNFORMATTED", STATUS = "Replace")
        REWIND StiffTmp
        WRITE (StiffTmp) Stiff(1:NWK)
        NNC=NC*(NC+1)/2
        !THE PARAMETERS NC AND/OR NRestart MUST BE INCREASED IF A SOLUTION HAS NOT CONVERGED
        call LANCZOS(Stiff, Mass, MAXA, EignVec, EignVal, NN, NNM, NWK, NWK, NRoot, RTol, NC, NNC, NRestart, IFSS, IFPR, StiffTmp, IOUT)
        REWIND StiffTmp
        READ (StiffTmp) Stiff
        close(StiffTmp)
    else
        uplo = 'U'
        allocate (res(NC))
        !call dfeast_scsrgv(uplo, NN, Stiff, iStiff, jStiff, Mass, iMass, jMass, fpm, epsout, loop, emin, emax, NC, EignVal, EignVec, NRoot, res, info)
        deallocate (res)
    end if
    write(IOUT,*)'-------------------------------------------------------------------------------------'
    deallocate (EignVec, EignVal)
end subroutine EIGENVAL

