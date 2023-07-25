interface
SUBROUTINE PROPAGS1 (F1, F3, NINF, NSUP, KIJS, KIJL,           &
 &                   BLK2GLO,                                  &
 &                   DEPTH_EXT,                                &
 &                   CGROUP_EXT, OMOSNH2KD_EXT,                &
 &                   DELLAM1_EXT, COSPHM1_EXT,                 &
 &                   U_EXT, V_EXT,                             &
 &                   L1STCALL)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWDRVTYPE , ONLY : WVGRIDGLO
 use yowparam , only:&
 & niblo,&
 & nang,&
 & nfre_red
REAL(KIND=JWRB),DIMENSION(NINF:NSUP+1, NANG, NFRE_RED), INTENT(IN) :: F1
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP+1, NANG, NFRE_RED), INTENT(OUT) :: F3
      INTEGER(KIND=JWIM), INTENT(IN) :: NINF, NSUP
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN):: DEPTH_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: DELLAM1_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: COSPHM1_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(IN) :: CGROUP_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(IN) :: OMOSNH2KD_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: U_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: V_EXT
      LOGICAL, INTENT(IN) :: L1STCALL
END SUBROUTINE PROPAGS1
end interface
