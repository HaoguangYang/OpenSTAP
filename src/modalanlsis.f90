module Modal
    integer, parameter :: StiffTmp = 16
contains
    
subroutine EIGENVAL(Stiff, Mass, MAXA, NN, NWK, NNM, NRoot)
    USE GLOBALS, ONLY : IOUT, PARDISODOOR
    implicit none
    real(8) :: Stiff(NWK), Mass(NWK)
    integer :: MAXA(NNM)
    integer :: NN, NWK, NNM
    real(8), parameter :: RTol = 1.0D-6
    integer :: NC, NRoot, NNC, NRestart
    real(8), allocatable :: EignVec(:,:),EignVal(:)
    logical :: IFSS, IFPR
    
    open(StiffTmp, FILE = "Stff.tmp", FORM = "UNFORMATTED", STATUS = "Replace")
    REWIND StiffTmp
    WRITE (StiffTmp) Stiff(1:NWK)
    write(IOUT,*)'-------------------------------------------------------------------------------------'
    write(IOUT,*)'           E I G E N   V A L U E   C A L C U L A T I O N   R E S U L T S'
    
    if (.NOT. PARDISODOOR) then
        NC = minval((/2*NRoot, NRoot+8, NWK/))
        allocate (EignVec(NN,NC),EignVal(NC))
        NNC=NC*(NC+1)/2
    
        !THE PARAMETERS NC AND/OR NRestart MUST BE INCREASED IF A SOLUTION HAS NOT CONVERGED
        call LANCZOS(Stiff, Mass, MAXA, EignVec, EignVal, NN, NNM, NWK, NWK, NRoot, RTol, NC, NNC, NRestart, IFSS, IFPR, StiffTmp, IOUT)
    !else
        !call dfeast_scsrgv(uplo, n, a, ia, ja, b, ib, jb, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)
    end if
    write(IOUT,*)'-------------------------------------------------------------------------------------'
    deallocate (EignVec, EignVal)
    if (.NOT. PARDISODOOR) then
        REWIND StiffTmp
        READ (StiffTmp) Stiff
        close(StiffTmp)
    end if
end subroutine EIGENVAL

end module Modal
