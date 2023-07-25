! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE HALPHAP_SCC_MOD
  CONTAINS
  SUBROUTINE HALPHAP_SCC (KIJS, KIJL, WAVNUM, COSWDIF, FL1, HALP)
    
    ! ----------------------------------------------------------------------
    
    !**** *HALPHAP* - COMPUTATION OF 1/2 PHILLIPS PARAMETER
    
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *HALPHAP(KIJS, KIJL, WAVNUM, UDIR, FL1, HALP)
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *WAVNUM* - WAVE NUMBER
    !          *COSWDIF*- COSINE ( WIND SPEED DIRECTION - WAVE DIRECTIONS)
    !          *FL1*    - SPECTRA
    !          *HALP*   - 1/2 PHILLIPS PARAMETER
    
    !     METHOD.
    !     -------
    
    ! ----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: FR, TH, FR5, DELTH, WETAIL, FRTAIL, DFIM, DFIMOFR
    USE YOWPARAM, ONLY: NANG, NFRE, NANG_PARAM
    USE YOWPCONS, ONLY: G, ZPI, ZPI4GM2, EPSMIN
    USE YOWPHYS, ONLY: ALPHAPMAX
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: HALP
    
    INTEGER(KIND=JWIM) :: IJ, K, M
    
    REAL(KIND=JWRB) :: ZLNFRNFRE
    REAL(KIND=JWRB) :: DELT25, DELT2, DEL2
    REAL(KIND=JWRB) :: TEMP1, TEMP2
    REAL(KIND=JWRB) :: ALPHAP
    REAL(KIND=JWRB) :: XMSS, EM, FM, F1D
    REAL(KIND=JWRB), DIMENSION(NANG_PARAM) :: FLWD
!$acc routine vector
!$acc data present( WAVNUM, COSWDIF, FL1, HALP )
    
!$acc loop vector private( FLWD )
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      ZLNFRNFRE = LOG(FR(NFRE))
      
      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH
      
      ! Find spectrum in wind direction
!$acc loop seq
      DO M=1,NFRE
!$acc loop seq
        DO K=1,NANG
          FLWD(K) = FL1(IJ, K, M)*0.5_JWRB + 0.5_JWRB*SIGN(1.0_JWRB, COSWDIF(IJ, K))
        END DO
        
        XMSS = 0._JWRB
        TEMP1 = DFIM(M)*WAVNUM(IJ, M)**2
        TEMP2 = 0.0_JWRB
!$acc loop seq
        DO K=1,NANG
          TEMP2 = TEMP2 + FLWD(K)
        END DO
        XMSS = XMSS + TEMP1*TEMP2
        
        K = 1
        EM = 0._JWRB
        FM = 0._JWRB
        TEMP2 = MAX(FLWD(K), EPSMIN)
!$acc loop seq
        DO K=2,NANG
          TEMP2 = TEMP2 + MAX(FLWD(K), EPSMIN)
        END DO
        EM = EM + TEMP2*DFIM(M)
        FM = FM + DFIMOFR(M)*TEMP2
      END DO
      
!$acc loop seq
      DO K=1,NANG
        FLWD(K) = FL1(IJ, K, NFRE)*0.5_JWRB + 0.5_JWRB*SIGN(1.0_JWRB, COSWDIF(IJ, K))
      END DO
      
      EM = EM + DELT25*TEMP2
      FM = FM + DELT2*TEMP2
      FM = EM / FM
      FM = MAX(FM, FR(1))
      
      IF (EM > 0.0_JWRB .and. FM < FR(NFRE - 2)) THEN
        ALPHAP = XMSS / (ZLNFRNFRE - LOG(FM))
        IF (ALPHAP > ALPHAPMAX) THEN
          ! some odd cases, revert to tail value
          F1D = 0.0_JWRB
!$acc loop seq
          DO K=1,NANG
            F1D = F1D + FLWD(K)*DELTH
          END DO
          ALPHAP = ZPI4GM2*FR5(NFRE)*F1D
        END IF
      ELSE
        F1D = 0.0_JWRB
!$acc loop seq
        DO K=1,NANG
          F1D = F1D + FLWD(K)*DELTH
        END DO
        ALPHAP = ZPI4GM2*FR5(NFRE)*F1D
      END IF
      
      !     1/2 ALPHAP:
      HALP(IJ) = 0.5_JWRB*MIN(ALPHAP, ALPHAPMAX)
      
      
    END DO
!$acc end data
  END SUBROUTINE HALPHAP_SCC
END MODULE HALPHAP_SCC_MOD
