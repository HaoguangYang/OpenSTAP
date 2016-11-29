关键变量名对应表
旧变量名-->新变量名                单元用到的变量

globals.f90
IELMNT  --> ElementTmpFile
ILOAD   --> LoadTmpFile
IIN     --> InputFile                   *
IOUT    --> OutputFile                  *
NUMNP   --> NumberOfNodalPoints         *
NEQ     --> NumberOfEquations           *
NWK     --> NumberOfMatrixElements
MK      --> MaxHalfBandwidth
IND     --> SolutionPhase               *
NUMEG   --> NumberOfElementGroups
MODEX   --> SolutionMode
TIM     --> Timing
HED     --> JobIdentifier
MIDEST  --> ElementGroupArraySize       *
MAXEST  --> MaxElementGroupArraySize
NG      --> CurrentElementGroup         *

memalloc.f90
NO CHANGE

assem.f90
MHT     --> ColumnHeight
ND      --> ElementDOF
LM      --> ElementLocationMatrix
LS      --> SmallestIndex
II      --> CurrentIndex
ME      --> IndexDisparity
AA      --> ElementStiffnessMatrix
A       --> SkylineK
II (addban) --> IndexI
JJ (addban) --> IndexJ
KK (addban) --> SkylineIndex
NN (colsol) --> NumberOfEquations
KKK         --> SolutionMode
V           --> LoadToDisplacement

datain.f90
ID      --> NodeIdentifier
NOD     --> NOdeOfLoad
IDIRN   --> LoadDirection
Added L as index for direct access of load temp file.
NLOAD   --> NumberOfConcentratedLoads
L is referenced later in LOADV as CurrentLoadCase.

elcal.f90
NPAR1   --> ElementType
AA      --> ElementGroupData

为方便组会讨论之用。

