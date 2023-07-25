interface
SUBROUTINE CTUWINI (KIJS, KIJL, NINF, NSUP, BLK2GLO, COSPHM1_EXT,   &
 &                  WLATM1, WCORM1, DP)
use parkind_wave, only:&
 & jwim,&
 & jwrb
USE YOWDRVTYPE , ONLY : WVGRIDGLO
use yowparam , only:&
 & niblo
INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL  ! GRID POINT INDEXES
INTEGER(KIND=JWIM), INTENT(IN) :: NINF, NSUP  ! HALO EXTEND NINF:NSUP+1
TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO  ! BLOCK TO GRID TRANSFORMATION
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: COSPHM1_EXT ! 1/COSPH
REAL(KIND=JWRB), DIMENSION(NINF:NSUP,2), INTENT(OUT) :: WLATM1 ! 1 - WLAT
REAL(KIND=JWRB), DIMENSION(NINF:NSUP,4), INTENT(OUT) :: WCORM1 ! 1 - WCOR
REAL(KIND=JWRB), DIMENSION(NINF:NSUP,2), INTENT(OUT) :: DP     ! COS PHI FACTOR
END SUBROUTINE CTUWINI
end interface