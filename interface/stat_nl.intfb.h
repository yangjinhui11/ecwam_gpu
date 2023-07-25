interface
 SUBROUTINE STAT_NL(KIJS, KIJL, &
 & XM0, XK0, BF2, XNU, SIG_TH, DPTH, &
 & C3, C4, ETA_M, R, C4_B, C4_DYN)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: XM0, XK0, BF2, XNU, SIG_TH, DPTH
 REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(OUT) :: C3, C4, ETA_M, R, C4_B, C4_DYN
 END SUBROUTINE STAT_NL
end interface
