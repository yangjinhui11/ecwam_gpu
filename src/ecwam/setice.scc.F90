! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!-----------------------------------------------------------------------
MODULE SETICE_SCC_MOD
  CONTAINS
  SUBROUTINE SETICE_SCC (KIJS, KIJL, FL1, CICOVER, COSWDIF)
    
    !-----------------------------------------------------------------------
    
    !**** *SETICE* ROUTINE TO SET SPECTRA ON ICE TO NOISE LEVEL.
    
    !     R.PORTZ      MPI         OKT.1992
    !     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING
    
    !     PURPOSE.
    !     -------
    
    !          *SETICE* SET ICE SPECTRA (FL1) TO NOISE LEVEL
    
    !**   INTERFACE.
    !     ----------
    
    !         *CALL* *SETICE(KIJS, KIJL, FL1, CICOVER, WSWAVE, COSWDIF)*
    !          *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - LOCAL  INDEX OF LAST GRIDPOINT
    !          *FL1*     - SPECTRA
    !          *CICOVER* - SEA ICE COVER
    !          *WSWAVE*  - WIND SPEED.
    !          *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
    
    ! ----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWICE, ONLY: FLMIN, CITHRSH
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: EPSMIN
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CICOVER
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF
    
    
    INTEGER(KIND=JWIM) :: IJ, M, K
    
    REAL(KIND=JWRB) :: CIREDUC, TEMP, ICEFREE
!$acc routine vector
!$acc data present( FL1, CICOVER, COSWDIF )
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      ! ----------------------------------------------------------------------
      
      
      !*    1. SET SPECTRA TO NOISE LEVEL OVER ICE POINTS.
      !     ----------------------------------------------
      
      IF (CICOVER(IJ) > CITHRSH) THEN
        CIREDUC = MAX(EPSMIN, (1.0_JWRB - CICOVER(IJ)))
        ICEFREE = 0.0_JWRB
      ELSE
        CIREDUC = 0.0_JWRB
        ICEFREE = 1.0_JWRB
      END IF
      
      TEMP = CIREDUC*FLMIN
!$acc loop seq
      DO M=1,NFRE
!$acc loop seq
        DO K=1,NANG
          FL1(IJ, K, M) = FL1(IJ, K, M)*ICEFREE + TEMP*MAX(0.0_JWRB, COSWDIF(IJ, K))**2
        END DO
      END DO
      
      
    END DO
!$acc end data
  END SUBROUTINE SETICE_SCC
END MODULE SETICE_SCC_MOD
