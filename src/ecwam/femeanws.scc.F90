! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE FEMEANWS_SCC_MOD
  CONTAINS
  SUBROUTINE FEMEANWS_SCC (KIJS, KIJL, FL1, XLLWS, FM, EM)
    
    ! ----------------------------------------------------------------------
    
    !**** *FEMEANWS* - COMPUTATION OF MEAN ENERGY, MEAN FREQUENCY
    !                  FOR WINDSEA PART OF THE SPECTRUM AS DETERMINED
    !                  BY XLLWS
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTE MEAN FREQUENCY AT EACH GRID POINT FOR PART OF THE
    !       SPECTRUM WHERE XLLWS IS NON ZERO.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *FEMEANWS (KIJS, KIJL, FL1, XLLWS, EM, FM)*
    !              *KIJS*   - INDEX OF FIRST GRIDPOINT
    !              *KIJL*   - INDEX OF LAST GRIDPOINT
    !              *FL1*    - SPECTRUM.
    !              *XLLWS* - TOTAL WINDSEA MASK FROM INPUT SOURCE TERM
    !              *EM*     - MEAN WAVE ENERGY (OUTPUT)
    !              *FM*     - MEAN WAVE FREQUENCY (OUTPUT)
    
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
    
    USE YOWFRED, ONLY: FR, DFIM, DFIMOFR, DELTH, WETAIL, FRTAIL
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: EPSMIN
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1, XLLWS
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: FM
    REAL(KIND=JWRB), OPTIONAL, INTENT(OUT), DIMENSION(KIJL) :: EM
    
    
    INTEGER(KIND=JWIM) :: IJ, M, K
    
    REAL(KIND=JWRB) :: DELT25, DELT2
    REAL(KIND=JWRB) :: TEMP2, EM_LOC
!$acc routine vector
!$acc data present( FL1, XLLWS, FM, EM )
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      !*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
      !        ------------------------------------------------
      
      EM_LOC = EPSMIN
      FM(IJ) = EPSMIN
      
      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH
      
      
      !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
      !        ------------------------------------------
      
!$acc loop seq
      DO M=1,NFRE
        TEMP2 = 0.0_JWRB
!$acc loop seq
        DO K=1,NANG
          TEMP2 = TEMP2 + XLLWS(IJ, K, M)*FL1(IJ, K, M)
        END DO
        EM_LOC = EM_LOC + DFIM(M)*TEMP2
        FM(IJ) = FM(IJ) + DFIMOFR(M)*TEMP2
      END DO
      
      !*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
      !*       NORMALIZE WITH TOTAL ENERGY.
      !        ------------------------------------------
      
      EM_LOC = EM_LOC + DELT25*TEMP2
      FM(IJ) = FM(IJ) + DELT2*TEMP2
      FM(IJ) = EM_LOC / FM(IJ)
      
      IF (PRESENT(EM)) THEN
        EM(IJ) = EM_LOC
      END IF
      
      
    END DO
!$acc end data
  END SUBROUTINE FEMEANWS_SCC
END MODULE FEMEANWS_SCC_MOD
