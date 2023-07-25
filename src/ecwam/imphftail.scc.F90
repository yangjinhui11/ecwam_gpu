! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE IMPHFTAIL_SCC_MOD
  CONTAINS
  SUBROUTINE IMPHFTAIL_SCC (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1)
    ! ----------------------------------------------------------------------
    
    !**** *IMPHFTAIL* - IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM
    
    
    !*    PURPOSE.
    !     --------
    
    !     IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM ABOVE FREQUENCY INDEX MIJ
    
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *IMPHFTAIL (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1)
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - INDEX OF LAST GRIDPOINT
    !          *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
    !          *FLM*     - SPECTAL DENSITY MINIMUM VALUE
    !          *WAVNUM*  - WAVENUMBER
    !          *XK2CG*   - (WAVNUM)**2 * GROUP SPEED
    !          *FL1*     - SPECTRUM (INPUT AND OUTPUT).
    
    !     METHOD.
    !     -------
    
    !     EXTERNALS.
    !     ---------
    
    !     REFERENCE.
    !     ----------
    
    ! ----------------------------------------------------------------------
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPARAM, ONLY: NANG, NFRE
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(KIJL) :: MIJ
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: FLM
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, XK2CG
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FL1
    
    
    INTEGER(KIND=JWIM) :: IJ, K, M
    
    REAL(KIND=JWRB) :: AKM1, TFAC
    REAL(KIND=JWRB) :: TEMP1, TEMP2
!$acc routine vector
!$acc data present( MIJ, FLM, WAVNUM, XK2CG, FL1 )
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      !*    DIAGNOSTIC TAIL.
      !     ----------------
      
      TEMP1 = 1.0_JWRB / XK2CG(IJ, MIJ(IJ)) / WAVNUM(IJ, MIJ(IJ))
      
!$acc loop seq
      DO M=MIJ(IJ) + 1,NFRE
        TEMP2 = 1.0_JWRB / XK2CG(IJ, M) / WAVNUM(IJ, M)
        TEMP2 = TEMP2 / TEMP1
        
        !*    MERGE TAIL INTO SPECTRA.
        !     ------------------------
!$acc loop seq
        DO K=1,NANG
          TFAC = FL1(IJ, K, MIJ(IJ))
          FL1(IJ, K, M) = MAX(TEMP2*TFAC, FLM(IJ, K))
        END DO
      END DO
      
      ! ----------------------------------------------------------------------
      
      
    END DO
!$acc end data
  END SUBROUTINE IMPHFTAIL_SCC
END MODULE IMPHFTAIL_SCC_MOD
