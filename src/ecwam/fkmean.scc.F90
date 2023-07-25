! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE FKMEAN_SCC_MOD
  CONTAINS
  SUBROUTINE FKMEAN_SCC (KIJS, KIJL, FL1, WAVNUM, EM, FM1, F1, AK, XK)
    
    ! ----------------------------------------------------------------------
    
    !**** *FKMEAN* - COMPUTATION OF MEAN FREQUENCIES AT EACH GRID POINT
    !                AND MEAN WAVE NUMBER (based in  sqrt(k)*F moment) .
    !                COMPUTATION OF THE MEAN WAVE ENERGY WAS ALSO
    !                ADDED SUCH THAT A CALL TO FKMEAN DOES NOT NEED
    
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTE MEAN FREQUENCIES AND WAVE NUMBER AT EACH GRID POINT.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *FKMEAN (KIJS, KIJL, FL1, WAVNUM, EM, FM1, F1, AK, XK)*
    !              *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
    !              *KIJL*    - LOCAL INDEX OF LAST GRIDPOINT
    !              *FL1*     - SPECTRUM.
    !              *WAVNUM*  - WAVE NUMBER.
    !              *EM*      - MEAN WAVE ENERGY
    !              *FM1*     - MEAN WAVE FREQUENCY BASED ON (1/f)*FL1 INTEGRATION
    !              *F1*      - MEAN WAVE FREQUENCY BASED ON f*FL1 INTEGRATION
    !              *AK*      - MEAN WAVE NUMBER  BASED ON sqrt(1/k)*FL1 INTGRATION
    !                          ONLY FOR SHALLOW WATER RUNS.
    !!!                        AK IS STILL NEEDED IN SNONLIN !!!!
    !!!                        IF THE OLD FORMULATION IS USED.
    !              *XK*      - MEAN WAVE NUMBER  BASED ON sqrt(k)*FL1 INTEGRATION
    !                          ONLY FOR SHALLOW WATER RUNS.
    
    !     METHOD.
    !     -------
    
    !       NONE.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       NONE.
    
    ! ----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: FR, DFIM, DFIMOFR, DFIMFR, DELTH, WETAIL, FRTAIL, WP1TAIL
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: G, ZPI, EPSMIN
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM
    
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: EM, FM1, F1, AK, XK
    
    
    INTEGER(KIND=JWIM) :: IJ, M, K
    REAL(KIND=JWRB) :: DELT25, COEFM1, COEF1, COEFA, COEFX, SQRTK
    REAL(KIND=JWRB) :: TEMPA, TEMPX, TEMP2
!$acc routine vector
!$acc data present( FL1, WAVNUM, EM, FM1, F1, AK, XK )
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      
      !*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
      !        ------------------------------------------------
      
      EM(IJ) = EPSMIN
      FM1(IJ) = EPSMIN
      F1(IJ) = EPSMIN
      AK(IJ) = EPSMIN
      XK(IJ) = EPSMIN
      
      DELT25 = WETAIL*FR(NFRE)*DELTH
      COEFM1 = FRTAIL*DELTH
      COEF1 = WP1TAIL*DELTH*FR(NFRE)**2
      COEFA = COEFM1*SQRT(G) / ZPI
      COEFX = COEF1*(ZPI / SQRT(G))
      
      !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
      !        ------------------------------------------
      
      !*    2.2 SHALLOW WATER INTEGRATION.
      !         --------------------------
      
!$acc loop seq
      DO M=1,NFRE
        SQRTK = SQRT(WAVNUM(IJ, M))
        TEMPA = DFIM(M) / SQRTK
        TEMPX = SQRTK*DFIM(M)
        K = 1
        TEMP2 = FL1(IJ, K, M)
!$acc loop seq
        DO K=2,NANG
          TEMP2 = TEMP2 + FL1(IJ, K, M)
        END DO
        EM(IJ) = EM(IJ) + DFIM(M)*TEMP2
        FM1(IJ) = FM1(IJ) + DFIMOFR(M)*TEMP2
        F1(IJ) = F1(IJ) + DFIMFR(M)*TEMP2
        AK(IJ) = AK(IJ) + TEMPA*TEMP2
        XK(IJ) = XK(IJ) + TEMPX*TEMP2
      END DO
      
      !*      ADD TAIL CORRECTION TO MEAN FREQUENCY AND
      !*      NORMALIZE WITH TOTAL ENERGY.
      EM(IJ) = EM(IJ) + DELT25*TEMP2
      FM1(IJ) = FM1(IJ) + COEFM1*TEMP2
      FM1(IJ) = EM(IJ) / FM1(IJ)
      F1(IJ) = F1(IJ) + COEF1*TEMP2
      F1(IJ) = F1(IJ) / EM(IJ)
      AK(IJ) = AK(IJ) + COEFA*TEMP2
      AK(IJ) = (EM(IJ) / AK(IJ))**2
      XK(IJ) = XK(IJ) + COEFX*TEMP2
      XK(IJ) = (XK(IJ) / EM(IJ))**2
      
      
    END DO
!$acc end data
  END SUBROUTINE FKMEAN_SCC
END MODULE FKMEAN_SCC_MOD
