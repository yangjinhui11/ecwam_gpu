! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDIWBK_SCC_MOD
  CONTAINS
  SUBROUTINE SDIWBK_SCC (KIJS, KIJL, FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN)
    
    ! ----------------------------------------------------------------------
    
    !**** *SDIWBK* - COMPUTATION OF BOTTOM-INDUCED WAVE BREAKING DISSIPATION
    
    
    !*    PURPOSE.
    !     --------
    !       COMPUTE BOTTOM-INDUCED DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
    !       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
    !       OF DISSIPATION SOURCE FUNCTION.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SDIWBK (KIJS, KIJL, FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN)*
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - INDEX OF LAST GRIDPOINT
    !          *FL1*     - SPECTRUM.
    !          *FLD*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*      - TOTAL SOURCE FUNCTION ARRAY
    !          *DEPTH*   - WATER DEPTH
    !          *EMAXDPT* - MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
    !          *EMEAN*   - MEAN ENERGY DENSITY
    !          *F1MEAN*  - MEAN FREQUENCY BASED ON 1st MOMENT.
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCES.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    ! ----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPARAM, ONLY: NANG, NFRE, NFRE_RED
    USE YOWSTAT, ONLY: LBIWBK
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: DEPTH, EMAXDPT, EMEAN, F1MEAN
    
    
    INTEGER(KIND=JWIM) :: IJ, K, M, IC
    REAL(KIND=JWRB) :: ALPH, ARG, Q, Q_OLD, REL_ERR, EXPQ
    REAL(KIND=JWRB) :: SDS
    
    REAL, PARAMETER :: ALPH_B_J = 1.0_JWRB
    REAL, PARAMETER :: COEF_B_J = 2*ALPH_B_J
    REAL, PARAMETER :: DEPTHTRS = 50.0_JWRB
!$acc routine vector
!$acc data present( FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN )
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      !*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
      !*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
      !        --------------------------------------------------------------
      
      IF (LBIWBK) THEN
        !       (FOLLOWING BATTJES-JANSSEN AND BEJI)
        IF (DEPTH(IJ) < DEPTHTRS) THEN
          ALPH = 2.0_JWRB*EMAXDPT(IJ) / EMEAN(IJ)
          ARG = MIN(ALPH, 50.0_JWRB)
          Q_OLD = EXP(-ARG)
          !            USE NEWTON-RAPHSON METHOD
!$acc loop seq
          DO IC=1,15
            EXPQ = EXP(-ARG*(1.0_JWRB - Q_OLD))
            Q = Q_OLD - (EXPQ - Q_OLD) / (ARG*EXPQ - 1.0_JWRB)
            REL_ERR = ABS(Q - Q_OLD) / Q_OLD
            IF (REL_ERR < 0.00001_JWRB)             EXIT
            Q_OLD = Q
          END DO
          Q = MIN(Q, 1.0_JWRB)
          SDS = COEF_B_J*ALPH*Q*F1MEAN(IJ)
        END IF
        
!$acc loop seq
        DO M=1,NFRE_RED
!$acc loop seq
          DO K=1,NANG
            IF (DEPTH(IJ) < DEPTHTRS) THEN
              SL(IJ, K, M) = SL(IJ, K, M) - SDS*FL1(IJ, K, M)
              FLD(IJ, K, M) = FLD(IJ, K, M) - SDS
            END IF
          END DO
        END DO
        
      END IF
      
      
    END DO
!$acc end data
  END SUBROUTINE SDIWBK_SCC
END MODULE SDIWBK_SCC_MOD
