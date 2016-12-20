!*******************************************************************************
!   Copyright(C) 2004-2013 Intel Corporation. All Rights Reserved.
!
!   The source code, information  and  material ("Material") contained herein is
!   owned  by Intel Corporation or its suppliers or licensors, and title to such
!   Material remains  with Intel Corporation  or its suppliers or licensors. The
!   Material  contains proprietary information  of  Intel or  its  suppliers and
!   licensors. The  Material is protected by worldwide copyright laws and treaty
!   provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way  without Intel's  prior  express written  permission. No  license
!   under  any patent, copyright  or  other intellectual property rights  in the
!   Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
!   implication, inducement,  estoppel or  otherwise.  Any  license  under  such
!   intellectual  property  rights must  be express  and  approved  by  Intel in
!   writing.
!
!   *Third Party trademarks are the property of their respective owners.
!
!   Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
!   this  notice or  any other notice embedded  in Materials by Intel or Intel's
!   suppliers or licensors in any way.
!
!*******************************************************************************
!   Content : MKL PARDISO Fortran-90 header file
!
!           Contains PARDISO routine definition.
!           For CDECL use only.
!
!*******************************************************************************
!DEC$ IF .NOT. DEFINED( __MKL_PARDISO_F90 )

!DEC$ DEFINE __MKL_PARDISO_F90

      MODULE MKL_PARDISO_PRIVATE

        TYPE MKL_PARDISO_HANDLE; INTEGER(KIND=8) DUMMY; END TYPE

        INTEGER, PARAMETER :: PARDISO_OOC_FILE_NAME = 1

      END MODULE MKL_PARDISO_PRIVATE

      MODULE MKL_PARDISO
        USE MKL_PARDISO_PRIVATE

!
! Subroutine prototype for PARDISO
!

      INTERFACE PARDISO
        SUBROUTINE PARDISO_S( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(INOUT) :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR

          REAL(KIND=4),     INTENT(IN)    :: A(*)
          REAL(KIND=4),     INTENT(INOUT) :: B(*)
          REAL(KIND=4),     INTENT(OUT)   :: X(*)
        END SUBROUTINE PARDISO_S

        SUBROUTINE PARDISO_D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(INOUT) :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR

          REAL(KIND=8),     INTENT(IN)    :: A(*)
          REAL(KIND=8),     INTENT(INOUT) :: B(*)
          REAL(KIND=8),     INTENT(OUT)   :: X(*)
        END SUBROUTINE PARDISO_D

        SUBROUTINE PARDISO_SC( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(INOUT) :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR

          COMPLEX(KIND=4),  INTENT(IN)    :: A(*)
          COMPLEX(KIND=4),  INTENT(INOUT) :: B(*)
          COMPLEX(KIND=4),  INTENT(OUT)   :: X(*)
        END SUBROUTINE PARDISO_SC

        SUBROUTINE PARDISO_DC( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(INOUT) :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=8),  INTENT(IN)    :: A(*)
          COMPLEX(KIND=8),  INTENT(INOUT) :: B(*)
          COMPLEX(KIND=8),  INTENT(OUT)   :: X(*)
        END SUBROUTINE PARDISO_DC


        SUBROUTINE PARDISO_S_2D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(INOUT) :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          REAL(KIND=4),     INTENT(IN)    :: A(*)
          REAL(KIND=4),     INTENT(INOUT) :: B(N,*)
          REAL(KIND=4),     INTENT(OUT)   :: X(N,*)
        END SUBROUTINE PARDISO_S_2D

        SUBROUTINE PARDISO_D_2D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(INOUT) :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          REAL(KIND=8),     INTENT(IN)    :: A(*)
          REAL(KIND=8),     INTENT(INOUT) :: B(N,*)
          REAL(KIND=8),     INTENT(OUT)   :: X(N,*)
        END SUBROUTINE PARDISO_D_2D

        SUBROUTINE PARDISO_SC_2D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(INOUT) :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=4),  INTENT(IN)    :: A(*)
          COMPLEX(KIND=4),  INTENT(INOUT) :: B(N,*)
          COMPLEX(KIND=4),  INTENT(OUT)   :: X(N,*)
        END SUBROUTINE PARDISO_SC_2D

        SUBROUTINE PARDISO_DC_2D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(INOUT) :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=8),  INTENT(IN)    :: A(*)
          COMPLEX(KIND=8),  INTENT(INOUT) :: B(N,*)
          COMPLEX(KIND=8),  INTENT(OUT)   :: X(N,*)
        END SUBROUTINE PARDISO_DC_2D
      END INTERFACE
!
! Subroutine prototype for PARDISO_64
!
! Note: The pardiso_64 interface is not supported on IA-32 architecture.
!       If called on IA-32, error = -12 is returned.
!
      INTERFACE PARDISO_64
        SUBROUTINE PARDISO_S_64( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(INOUT) :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          REAL(KIND=4),             INTENT(IN)    :: A(*)
          REAL(KIND=4),             INTENT(INOUT) :: B(*)
          REAL(KIND=4),             INTENT(OUT)   :: X(*)
        END SUBROUTINE PARDISO_S_64

        SUBROUTINE PARDISO_D_64( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(INOUT) :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          REAL(KIND=8),             INTENT(IN)    :: A(*)
          REAL(KIND=8),             INTENT(INOUT) :: B(*)
          REAL(KIND=8),             INTENT(OUT)   :: X(*)
        END SUBROUTINE PARDISO_D_64

        SUBROUTINE PARDISO_SC_64( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(INOUT) :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=4),          INTENT(IN)    :: A(*)
          COMPLEX(KIND=4),          INTENT(INOUT) :: B(*)
          COMPLEX(KIND=4),          INTENT(OUT)   :: X(*)
        END SUBROUTINE PARDISO_SC_64

        SUBROUTINE PARDISO_DC_64( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(INOUT) :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=8),          INTENT(IN)    :: A(*)
          COMPLEX(KIND=8),          INTENT(INOUT) :: B(*)
          COMPLEX(KIND=8),          INTENT(OUT)   :: X(*)
        END SUBROUTINE PARDISO_DC_64

        SUBROUTINE PARDISO_S_64_2D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(INOUT) :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          REAL(KIND=4),             INTENT(IN)    :: A(*)
          REAL(KIND=4),             INTENT(INOUT) :: B(N,*)
          REAL(KIND=4),             INTENT(OUT)   :: X(N,*)
        END SUBROUTINE PARDISO_S_64_2D

        SUBROUTINE PARDISO_D_64_2D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(INOUT) :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          REAL(KIND=8),             INTENT(IN)    :: A(*)
          REAL(KIND=8),             INTENT(INOUT) :: B(N,*)
          REAL(KIND=8),             INTENT(OUT)   :: X(N,*)
        END SUBROUTINE PARDISO_D_64_2D

        SUBROUTINE PARDISO_SC_64_2D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(INOUT) :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=4),          INTENT(IN)    :: A(*)
          COMPLEX(KIND=4),          INTENT(INOUT) :: B(N,*)
          COMPLEX(KIND=4),          INTENT(OUT)   :: X(N,*)
        END SUBROUTINE PARDISO_SC_64_2D

        SUBROUTINE PARDISO_DC_64_2D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(INOUT) :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=8),          INTENT(IN)    :: A(*)
          COMPLEX(KIND=8),          INTENT(INOUT) :: B(N,*)
          COMPLEX(KIND=8),          INTENT(OUT)   :: X(N,*)
        END SUBROUTINE PARDISO_DC_64_2D

      END INTERFACE

      INTERFACE

       SUBROUTINE PARDISOINIT(PT, MTYPE, IPARM)
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(OUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(OUT) :: IPARM(*)
       END SUBROUTINE PARDISOINIT

      END INTERFACE

      INTERFACE PARDISO_GET

       FUNCTION PARDISO_GETENV(PT, OptName, StrVal)
          USE MKL_PARDISO_PRIVATE
          INTEGER PARDISO_GETENV
          TYPE(MKL_PARDISO_HANDLE), INTENT(IN)  :: PT(*)
          INTEGER,                  INTENT(IN)  :: OptName
          CHARACTER(*),             INTENT(OUT)  :: StrVal
       END FUNCTION PARDISO_GETENV

      END INTERFACE

      INTERFACE PARDISO_SET

       FUNCTION PARDISO_SETENV(PT, OptName, StrVal)
          USE MKL_PARDISO_PRIVATE
          INTEGER PARDISO_SETENV
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,                  INTENT(IN)    :: OptName
          INTEGER,                  INTENT(IN)    :: StrVal(*)
       END FUNCTION PARDISO_SETENV

      END INTERFACE

      INTERFACE PARDISO_PIV
      
      FUNCTION MKL_PARDISO_PIVOT( AII, BII, EPS)
         REAL(KIND=8)   :: AII, BII, EPS
         INTEGER   MKL_PARDISO_PIVOT
         END
      END INTERFACE PARDISO_PIV

      INTERFACE PARDISO_GETDIAG

         SUBROUTINE PARDISO_GETDIAG_D(PT, DIAG_FACT, DIAG_A, MNUM, EPS)
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          REAL(KIND=8), INTENT(INOUT)             :: DIAG_FACT, DIAG_A
          INTEGER, INTENT(IN)            :: MNUM
          INTEGER, INTENT(INOUT)         :: EPS
         END

         SUBROUTINE PARDISO_GETDIAG_Z(PT, DIAG_FACT, DIAG_A, MNUM, EPS)
          USE MKL_PARDISO_PRIVATE
          TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
          COMPLEX(KIND=8), INTENT(INOUT)          :: DIAG_FACT, DIAG_A
          INTEGER, INTENT(IN)            :: MNUM
          INTEGER, INTENT(INOUT)         :: EPS
         END

      END INTERFACE PARDISO_GETDIAG




      END MODULE MKL_PARDISO

!DEC$ ENDIF
