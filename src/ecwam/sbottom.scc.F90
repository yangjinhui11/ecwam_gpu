! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SBOTTOM_SCC_MOD
  CONTAINS
  SUBROUTINE SBOTTOM_SCC (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH)
    
    !SHALLOW
    ! ----------------------------------------------------------------------
    
    !**** *SBOTTOM* - COMPUTATION OF BOTTOM FRICTION.
    
    !     G.J.KOMEN AND Q.D.GAO
    !     OPTIMIZED BY L.F. ZAMBRESKY
    !     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTATION OF BOTTOM FRICTION DISSIPATION
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SBOTTOM (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH)
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - INDEX OF LAST GRIDPOINT
    !          *FL1*     - SPECTRUM.
    !          *FLD*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*      - TOTAL SOURCE FUNCTION ARRAY
    !          *WAVNUM*  - WAVE NUMBER
    !          *DEPTH*   - WATER DEPTH
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCES.
    
    !     REFERENCES.
    !     -----------
    
    !       HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)
    !       BOUWS AND KOMEN, JPO 13(1983)1653-1658
    
    ! ----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPARAM, ONLY: NANG, NFRE, NFRE_RED
    USE YOWPCONS, ONLY: GM1
    USE YOWSHAL, ONLY: BATHYMAX
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: DEPTH
    
    
    INTEGER(KIND=JWIM) :: IJ, K, M
    REAL(KIND=JWRB) :: CONST, ARG
    REAL(KIND=JWRB) :: SBO
!$acc routine vector
!$acc data present( FL1, FLD, SL, WAVNUM, DEPTH )
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      CONST = -2.0_JWRB*0.038_JWRB*GM1
!$acc loop seq
      DO M=1,NFRE_RED
        IF (DEPTH(IJ) < BATHYMAX) THEN
          ARG = 2.0_JWRB*DEPTH(IJ)*WAVNUM(IJ, M)
          ARG = MIN(ARG, 50.0_JWRB)
          SBO = CONST*WAVNUM(IJ, M) / SINH(ARG)
        ELSE
          SBO = 0.0_JWRB
        END IF
        
!$acc loop seq
        DO K=1,NANG
          SL(IJ, K, M) = SL(IJ, K, M) + SBO*FL1(IJ, K, M)
          FLD(IJ, K, M) = FLD(IJ, K, M) + SBO
        END DO
      END DO
      
      
    END DO
!$acc end data
  END SUBROUTINE SBOTTOM_SCC
END MODULE SBOTTOM_SCC_MOD
