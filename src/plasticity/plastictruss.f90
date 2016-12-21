SUBROUTINE PLASTICTRUSS
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the truss element subroutine         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N(9)

  NUME = NPAR(2)
  NUMMAT = NPAR(3)
  NPAR(5) = 2

! Calculate the pointer to the arrays in the element group data
! N1: E(NUMMAT)
! N2: AREA(NUMMAT)
! N3: YIELD STRESS (NUMMAT)
! N4: PLASTIC MODULES K (NUMMAT)
! N5: HISTORY OF THE ELEMENT //ELASTIC OR PLASTIC// (NUME)
! N6: LM(6,NUME)
! N7: XYZ(6,NUME)
! N8: MTAP(NUME) 
  N(1)=0
  N(2)=N(1)+NUMMAT*ITWO
  N(3)=N(2)+NUMMAT*ITWO
  N(4)=N(3)+NUMMAT*ITWO
  N(5)=N(4)+NUMMAT*ITWO
  N(6)=N(5)+NUME
  N(7)=N(6)+6*NUME
  N(8)=N(7)+6*NUME*ITWO
  N(9)=N(8)+NUME
  
  
  MIDEST=N(9)
  if (IND .EQ. 1) then
        ! Allocate storage for element group data
        call MemAlloc(11,"ELEGP",MIDEST,1)
  end if
  NFIRST = NP(11)   ! Pointer to the first entry in the element group data array in the unit of single precision (corresponding to A)
  N(:) = N(:) + NFIRST
  NLAST=N(9)

!  CALL PLASTICRUSS (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
!        A(N(1)),A(N(2)),A(N(3)),A(N(4)),A(N(5)),A(N(6)),A(N(7)))

  RETURN

END SUBROUTINE PLASTICTRUSS