! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDISSIP_SCC_MOD
  CONTAINS
  SUBROUTINE SDISSIP_SCC (KIJS, KIJL, FL1, FLD, SL, INDEP, WAVNUM, XK2CG, EMEAN, F1MEAN, XKMEAN, UFRIC, COSWDIF, RAORW)
    ! ----------------------------------------------------------------------
    
    !**** *SDISSIP* - COMPUTATION OF DEEP WATER DISSIPATION SOURCE FUNCTION.
    
    
    !*    PURPOSE.
    !     --------
    !       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
    !       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
    !       OF DISSIPATION SOURCE FUNCTION.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SDISSIP (KIJS, KIJL, FL1, FLD, SL, *
    !                        INDEP, WAVNUM, XK2CG,
    !                        EMEAN, F1MEAN, XKMEAN,*
    !                        UFRIC, COSWDIF, RAORW)*
    !         *KIJS* - INDEX OF FIRST GRIDPOINT
    !         *KIJL* - INDEX OF LAST GRIDPOINT
    !         *FL1*  - SPECTRUM.
    !         *FLD*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*  - TOTAL SOURCE FUNCTION ARRAY
    !       *INDEP*  - DEPTH INDEX
    !       *WAVNUM* - WAVE NUMBER
    !       *XK2CG*  - (WAVNUM)**2 * GROUP SPEED
    !        *EMEAN* - MEAN ENERGY DENSITY
    !       *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
    !       *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
    !       *UFRIC*  - FRICTION VELOCITY IN M/S.
    !       *RAORW*  - RATIO AIR DENSITY TO WATER DENSITY
    !       *COSWDIF*-  COS(TH(K)-WDWAVE(IJ))
    
    ! ----------------------------------------------------------------------
    USE SDISSIP_JAN_SCC_MOD, ONLY: SDISSIP_JAN_SCC
    USE SDISSIP_ARD_SCC_MOD, ONLY: SDISSIP_ARD_SCC
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWSTAT, ONLY: IPHYS
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL
    INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(KIJL) :: INDEP
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, XK2CG
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: EMEAN, F1MEAN, XKMEAN
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: UFRIC, RAORW
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF
    
!$acc routine vector
!$acc data present( FL1, FLD, SL, INDEP, WAVNUM, XK2CG, EMEAN, F1MEAN, XKMEAN, UFRIC, COSWDIF, RAORW )
    
    ! ----------------------------------------------------------------------
    
    
    SELECT CASE (IPHYS)
    CASE (0)
      CALL SDISSIP_JAN_SCC(KIJS, KIJL, FL1, FLD, SL, WAVNUM, EMEAN, F1MEAN, XKMEAN)
      
    CASE (1)
      CALL SDISSIP_ARD_SCC(KIJS, KIJL, FL1, FLD, SL, INDEP, WAVNUM, XK2CG, UFRIC, COSWDIF, RAORW)
    END SELECT
    
    
!$acc end data
  END SUBROUTINE SDISSIP_SCC
END MODULE SDISSIP_SCC_MOD
