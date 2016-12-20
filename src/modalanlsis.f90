subroutine EIGENVAL(Stiff, Mass, MAXA, NN, NWK, NNM, NRoot)
    USE GLOBALS, ONLY : IOUT
    implicit none
    real(8) :: Stiff(NWK), Mass(NWK)
    integer :: MAXA(NNM)
    integer :: NN, NWK, NNM
    real(8), parameter :: RTol = 1.0D-6
    integer, parameter :: StiffTmp = 16
    integer :: NC, NRoot, NNC, NRestart
    real(8), allocatable :: EignVec(:,:),EignVal(:)
    logical :: IFSS, IFPR
    
    write(IOUT,*)'------------------------------------------------------------------------------------'
    write(IOUT,*)'           E I G E N   V A L U E   C A L C U L A T I O N   R E S U L T S'
    

    NC = minval((/2*NRoot, NRoot+8, NWK/))
    allocate (EignVec(NN,NC),EignVal(NC))
    NNC=NC*(NC+1)/2
    open(StiffTmp, FILE = "Stff.tmp", FORM = "UNFORMATTED", STATUS = "Replace")
    !THE PARAMETERS NC AND/OR NRestart MUST BE INCREASED IF A SOLUTION HAS NOT CONVERGED
    call LANCZOS(Stiff, Mass, MAXA, EignVec, EignVal, NN, NNM, NWK, NWK, NRoot, RTol, NC, NNC, NRestart, IFSS, IFPR, StiffTmp, IOUT)
    write(IOUT,*)'------------------------------------------------------------------------------------'
    deallocate (EignVec, EignVal)
    REWIND StiffTmp
    READ (StiffTmp) Stiff
end subroutine EIGENVAL
