! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE Z0WAVE_SCC_MOD
  CONTAINS
  SUBROUTINE Z0WAVE_SCC (KIJS, KIJL, US, TAUW, UTOP, Z0, Z0B, CHRNCK)
    
    ! ----------------------------------------------------------------------
    
    !**** *Z0WAVE* - DETERMINE THE SEA STATE DEPENDENT ROUGHNESS LENGTH.
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTE ROUGHNESS LENGTH.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *Z0WAVE (KIJS, KIJL, US, TAUW, UTOP, Z0, Z0B, CHRNCK)
    !          *KIJS* - INDEX OF FIRST GRIDPOINT.
    !          *KIJL* - INDEX OF LAST GRIDPOINT.
    !          *US*   - OUTPUT BLOCK OF SURFACE STRESSES.
    !          *TAUW* - INPUT BLOCK OF WAVE STRESSES.
    !          *UTOP* - WIND SPEED.
    !          *Z0*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.
    !          *Z0B*  - BACKGROUND ROUGHNESS LENGTH.
    !          *CHRNCK- CHARNOCK COEFFICIENT
    
    !     METHOD.
    !     -------
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ---------
    
    !       NONE.
    
    ! ----------------------------------------------------------------------
    
    USE CHNKMIN_SCC_MOD, ONLY: CHNKMIN_SCC
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: LLCAPCHNK
    USE YOWPCONS, ONLY: G, GM1
    USE YOWPHYS, ONLY: ALPHA
    USE YOWTABL, ONLY: EPS1
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: US, TAUW, UTOP
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: Z0, Z0B, CHRNCK
    
    
    INTEGER(KIND=JWIM) :: IJ
    REAL(KIND=JWRB) :: UST2, UST3, ARG
    REAL(KIND=JWRB) :: ALPHAOG
!$acc routine vector
!$acc data present( US, TAUW, UTOP, Z0, Z0B, CHRNCK )
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      IF (LLCAPCHNK) THEN
        ALPHAOG = CHNKMIN_SCC(UTOP(IJ))*GM1
      ELSE
        ALPHAOG = ALPHA*GM1
      END IF
      
      UST2 = US(IJ)**2
      UST3 = US(IJ)**3
      ARG = MAX(UST2 - TAUW(IJ), EPS1)
      Z0(IJ) = ALPHAOG*UST3 / SQRT(ARG)
      Z0B(IJ) = ALPHAOG*UST2
      CHRNCK(IJ) = G*Z0(IJ) / UST2
      
      
    END DO
!$acc end data
  END SUBROUTINE Z0WAVE_SCC
END MODULE Z0WAVE_SCC_MOD
