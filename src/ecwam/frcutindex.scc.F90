! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE FRCUTINDEX_SCC_MOD
  CONTAINS
  SUBROUTINE FRCUTINDEX_SCC (KIJS, KIJL, FM, FMWS, UFRIC, CICOVER, MIJ, RHOWGDFTH)
    
    ! ----------------------------------------------------------------------
    
    !**** *FRCUTINDEX* - RETURNS THE LAST FREQUENCY INDEX OF
    !                    PROGNOSTIC PART OF SPECTRUM.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *FRCUTINDEX (KIJS, KIJL, FM, FMWS, CICOVER, MIJ, RHOWGDFTH)
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FM*     - MEAN FREQUENCY
    !          *FMWS*   - MEAN FREQUENCY OF WINDSEA
    !          *UFRIC*  - FRICTION VELOCITY IN M/S
    !          *CICOVER*- CICOVER
    !          *MIJ*    - LAST FREQUENCY INDEX for imposing high frequency tail
    !          *RHOWGDFTH - WATER DENSITY * G * DF * DTHETA
    !                       FOR TRAPEZOIDAL INTEGRATION BETWEEN FR(1) and FR(MIJ)
    !                       !!!!!!!!  RHOWGDFTH=0 FOR FR > FR(MIJ)
    
    
    !     METHOD.
    !     -------
    
    !*    COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
    !*    FREQUENCIES LE 2.5*MAX(FMWS,FM).
    
    
    !!! be aware that if this is NOT used, for iphys=1, the cumulative dissipation has to be
    !!! re-activated (see module yowphys) !!!
    
    
    ! ----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: FR, DFIM, FRATIO, FLOGSPRDM1, ZPIFR, DELTH, RHOWG_DFIM, FRIC
    USE YOWICE, ONLY: CITHRSH_TAIL
    USE YOWPARAM, ONLY: NFRE
    USE YOWPCONS, ONLY: G, EPSMIN
    USE YOWPHYS, ONLY: TAILFACTOR, TAILFACTOR_PM
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(KIJL)
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: FM, FMWS, UFRIC, CICOVER
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NFRE) :: RHOWGDFTH
    
    
    INTEGER(KIND=JWIM) :: IJ, M
    
    REAL(KIND=JWRB) :: FPMH, FPPM, FM2, FPM, FPM4
!$acc routine vector
!$acc data present( FM, FMWS, UFRIC, CICOVER, MIJ, RHOWGDFTH )
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      !*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
      !*    FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
      !*    WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
      !*    VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
      !     ------------------------------------------------------------
      
      FPMH = TAILFACTOR / FR(1)
      FPPM = TAILFACTOR_PM*G / (FRIC*ZPIFR(1))
      
      IF (CICOVER(IJ) <= CITHRSH_TAIL) THEN
        FM2 = MAX(FMWS(IJ), FM(IJ))*FPMH
        FPM = FPPM / MAX(UFRIC(IJ), EPSMIN)
        FPM4 = MAX(FM2, FPM)
        MIJ(IJ) = NINT(LOG10(FPM4)*FLOGSPRDM1) + 1
        MIJ(IJ) = MIN(MAX(1, MIJ(IJ)), NFRE)
      ELSE
        MIJ(IJ) = NFRE
      END IF
      
      !     SET RHOWGDFTH
!$acc loop seq
      DO M=1,MIJ(IJ)
        RHOWGDFTH(IJ, M) = RHOWG_DFIM(M)
      END DO
      IF (MIJ(IJ) /= NFRE)       RHOWGDFTH(IJ, MIJ(IJ)) = 0.5_JWRB*RHOWGDFTH(IJ, MIJ(IJ))
!$acc loop seq
      DO M=MIJ(IJ) + 1,NFRE
        RHOWGDFTH(IJ, M) = 0.0_JWRB
      END DO
      
      
    END DO
!$acc end data
  END SUBROUTINE FRCUTINDEX_SCC
END MODULE FRCUTINDEX_SCC_MOD
