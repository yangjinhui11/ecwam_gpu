! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE AIRSEA_SCC_MOD
  CONTAINS
  SUBROUTINE AIRSEA_SCC (KIJS, KIJL, HALP, U10, U10DIR, TAUW, TAUWDIR, RNFAC, US, Z0, Z0B, CHRNCK, ICODE_WND, IUSFG)
    
    ! ----------------------------------------------------------------------
    
    !**** *AIRSEA* - DETERMINE TOTAL STRESS IN SURFACE LAYER.
    
    !     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
    !     JEAN BIDLOT         ECMWF     FEBRUARY 1999 : TAUT is already
    !                                                   SQRT(TAUT)
    !     JEAN BIDLOT         ECMWF     OCTOBER 2004: QUADRATIC STEP FOR
    !                                                 TAUW
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTE TOTAL STRESS.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *AIRSEA (KIJS, KIJL, FL1, WAVNUM,
    !                       HALP, U10, U10DIR, TAUW, TAUWDIR, RNFAC,
    !                       US, Z0, Z0B, CHRNCK, ICODE_WND, IUSFG)*
    
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT.
    !          *KIJL*    - INDEX OF LAST GRIDPOINT.
    !          *FL1*     - SPECTRA
    !          *WAVNUM*  - WAVE NUMBER
    !          *HALP*    - 1/2 PHILLIPS PARAMETER
    !          *U10*     - WINDSPEED U10.
    !          *U10DIR*  - WINDSPEED DIRECTION.
    !          *TAUW*    - WAVE STRESS.
    !          *TAUWDIR* - WAVE STRESS DIRECTION.
    !          *RNFAC*   - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !          *US*      - OUTPUT OR OUTPUT BLOCK OF FRICTION VELOCITY.
    !          *Z0*      - OUTPUT BLOCK OF ROUGHNESS LENGTH.
    !          *Z0B*     - BACKGROUND ROUGHNESS LENGTH.
    !          *CHRNCK*  - CHARNOCK COEFFICIENT
    !          *ICODE_WND* SPECIFIES WHICH OF U10 OR US HAS BEEN FILED UPDATED:
    !                     U10: ICODE_WND=3 --> US will be updated
    !                     US:  ICODE_WND=1 OR 2 --> U10 will be updated
    !          *IUSFG*   - IF = 1 THEN USE THE FRICTION VELOCITY (US) AS FIRST GUESS in TAUT_Z0
    !                           0 DO NOT USE THE FIELD US
    
    
    ! ----------------------------------------------------------------------
    
    USE Z0WAVE_SCC_MOD, ONLY: Z0WAVE_SCC
    USE TAUT_Z0_SCC_MOD, ONLY: TAUT_Z0_SCC
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPHYS, ONLY: XKAPPA, XNLEV
    USE YOWTEST, ONLY: IU06
    USE YOWWIND, ONLY: WSPMIN
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL, ICODE_WND, IUSFG
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: HALP, U10DIR, TAUW, TAUWDIR, RNFAC
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: U10, US
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: Z0, Z0B, CHRNCK
    
    INTEGER(KIND=JWIM) :: IJ, I, J
    
    REAL(KIND=JWRB) :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, UST2, ARG, SQRTCDM1
    REAL(KIND=JWRB) :: XKAPPAD, XLOGLEV
    REAL(KIND=JWRB) :: XLEV
!$acc routine vector
!$acc data present( HALP, U10, U10DIR, TAUW, TAUWDIR, RNFAC, US, Z0, Z0B, CHRNCK )
    
    ! ----------------------------------------------------------------------
    
    !*    2. DETERMINE TOTAL STRESS (if needed)
    !        ----------------------------------
    
    IF (ICODE_WND == 3) THEN
      
      CALL TAUT_Z0_SCC(KIJS, KIJL, IUSFG, HALP, U10, U10DIR, TAUW, TAUWDIR, RNFAC, US, Z0, Z0B, CHRNCK)
      
    ELSE IF (ICODE_WND == 1 .or. ICODE_WND == 2) THEN
      
      !*    3. DETERMINE ROUGHNESS LENGTH (if needed).
      !        ---------------------------
      
      CALL Z0WAVE_SCC(KIJS, KIJL, US, TAUW, U10, Z0, Z0B, CHRNCK)
      
!$acc loop vector
      DO IJ=KIJS,KIJL
        
        !*    3. DETERMINE U10 (if needed).
        !        ---------------------------
        
        XKAPPAD = 1.0_JWRB / XKAPPA
        XLOGLEV = LOG(XNLEV)
        
        U10(IJ) = XKAPPAD*US(IJ)*(XLOGLEV - LOG(Z0(IJ)))
        U10(IJ) = MAX(U10(IJ), WSPMIN)
        
      END DO
    END IF
    
    
!$acc end data
  END SUBROUTINE AIRSEA_SCC
END MODULE AIRSEA_SCC_MOD
