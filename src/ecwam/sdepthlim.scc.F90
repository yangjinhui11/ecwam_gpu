! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDEPTHLIM_SCC_MOD
  CONTAINS
  SUBROUTINE SDEPTHLIM_SCC (KIJS, KIJL, EMAXDPT, FL1)
    ! ----------------------------------------------------------------------
    !     J. BIDLOT    ECMWF   NOVEMBER 2017
    
    !*    PURPOSE.
    !     --------
    !     LIMITS THE SPECTRAL VARIANCE SUCH THAT THE TOTAL VARIANCE
    !     DOES NOT EXCEED THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
    
    !**   INTERFACE.
    !     ----------
    !     *CALL* *SDEPTHLIM((KIJS, KIJL, EMAXDPT, FL1)
    !          *KIJS*   - LOCAL INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - LOCAL  INDEX OF LAST GRIDPOIN
    !          *EMAXDPT - MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
    !          *FL1*    - SPECTRUM.
    
    
    !     METHOD.
    !     -------
    
    !     EXTERNALS.
    !     ----------
    
    !     REFERENCE.
    !     ----------
    !     NONE
    
    ! ----------------------------------------------------------------------
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: EPSMIN
    USE YOWFRED, ONLY: FR, DFIM, DELTH, WETAIL
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: EMAXDPT
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FL1
    
    INTEGER(KIND=JWIM) :: IJ, K, M
    REAL(KIND=JWRB) :: DELT25
    REAL(KIND=JWRB) :: EM
    REAL(KIND=JWRB) :: TEMP
    LOGICAL :: LLEPSMIN
!$acc routine vector
!$acc data present( EMAXDPT, FL1 )
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      EM = EPSMIN
!$acc loop seq
      DO M=1,NFRE
        K = 1
        TEMP = FL1(IJ, K, M)
!$acc loop seq
        DO K=2,NANG
          TEMP = TEMP + FL1(IJ, K, M)
        END DO
        EM = EM + DFIM(M)*TEMP
      END DO
      ! ----------------------------------------------------------------------
      
      !*    3. ADD TAIL ENERGY.
      !        ----------------
      
      DELT25 = WETAIL*FR(NFRE)*DELTH
      EM = EM + DELT25*TEMP
      
      EM = MIN(EMAXDPT(IJ) / EM, 1.0_JWRB)
      
!$acc loop seq
      DO M=1,NFRE
!$acc loop seq
        DO K=1,NANG
          FL1(IJ, K, M) = MAX(FL1(IJ, K, M)*EM, EPSMIN)
        END DO
      END DO
      
      
    END DO
!$acc end data
  END SUBROUTINE SDEPTHLIM_SCC
END MODULE SDEPTHLIM_SCC_MOD
