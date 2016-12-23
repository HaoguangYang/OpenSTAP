! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       .
! .                            S T A P 9 0                                .
! .                                                                       .
! .     AN IN-CORE SOLUTION STATIC ANALYSIS PROGRAM IN FORTRAN 90         .
! .     Adapted from STAP (KJ Bath, FORTRAN IV) for teaching purpose      .
! .                                                                       .
! .     Xiong Zhang, (2013)                                               .
! .     Computational Dynamics Group, School of Aerospace                 .
! .     Tsinghua Univerity                                                .
! .                                                                       .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .

! .  Define global variables

module GLOBALS

   integer, parameter :: IELMNT=1	! Unit storing element data 
   integer, parameter :: ILOAD=2	! Unit storing load vectors
   integer, parameter :: IIN=5		! Unit used for input
   integer, parameter :: IOUT=6		! Unit used for output
   integer, parameter :: VTKFile=7
   integer, parameter :: VTKTmpFile=8
   integer, parameter :: VTKNodeTmp=9
   integer, parameter :: VTKElTypTmp=10

   integer :: DIM       ! Dimension of the problem
   integer :: NUMNP		! Total number of nodal points
						! = 0 : Program stop
   integer :: NEQ		! Number of equations
   integer :: NWK		! Number of matrix elements 注意这个东西在skyline和pardiso的情况下不同
   integer :: MK		! Maximum half bandwidth

   integer :: IND		! Solution phase indicator
						!   1 - Read and generate element information
						!   2 - Assemble structure stiffness matrix
						!   3 - Stress calculations
   integer :: NPAR(10)	! Element group control data
						!   NPAR(1) - Element type
						!             1 : Truss element
                        !             2 : 4Q
                        !  
						!   NPAR(2) - Number of elements
						!   NPAR(3) - Number of different sets of material and 
						!             cross-sectional  constants
						!   NPAR(4) - Consider Gravity Or Not
						!   NPAR(5) - Numder of Element Shape Nodes
						!   NPAR(6) - Number of Element Load Interplotation Nodes
   integer :: NUMEG		! Total number of element groups, > 0
   integer :: NEL       ! Total number of elements (for post-processing .vtk)
   integer :: NCONECT   ! Total element numbers in connection matrix
   integer :: NLCASE
   integer :: CURLCASE  ! Current Load Case
   
   integer :: MODEX		! Solution mode: 0 - data check only;  1 -  execution                                   

   integer(8) :: TIM(10)		! Timing information
   character*80 :: HED	! Master heading information for use in labeling the output

   integer :: NFIRST
   integer :: NLAST
   integer :: MIDEST
   integer :: MAXEST

   integer :: NG
   
   logical :: BANDWIDTHOPT
   logical :: PARDISODOOR
   logical :: LOADANALYSIS
   logical :: DYNANALYSIS
   
   !FOR PLASTIC TRUSS ONLY
   !CONTROL MESSAGE
   LOGICAL :: PLASTICTRIAL =  .TRUE.
   LOGICAL :: PLASTICITERATION =  .FALSE.
   REAL(4) :: ITERATENUM = 0.0

end module GLOBALS
