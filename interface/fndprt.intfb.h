interface
 SUBROUTINE FNDPRT (KIJS, KIJL, NPMAX, &
 & NPEAK, MIJ, NTHP, NFRP, &
 & FLLOW, LLCOSDIFF, FLNOISE, &
 & FL1, SWM, &
 & ENE, DIR, PER)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWPARAM , ONLY : NANG ,NFRE
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL, NPMAX
 INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(KIJS:KIJL) :: MIJ
 INTEGER(KIND=JWIM), INTENT(INOUT), DIMENSION(KIJS:KIJL) :: NPEAK
 INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(KIJS:KIJL,NPMAX) :: NTHP, NFRP
 REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJS:KIJL) :: FLNOISE
 REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJS:KIJL,NANG,NFRE) :: FLLOW, FL1
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: SWM
 REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJS:KIJL,0:NPMAX) :: DIR, PER, ENE
 LOGICAL, INTENT(IN), DIMENSION(KIJS:KIJL,NANG) :: LLCOSDIFF
 END SUBROUTINE FNDPRT
end interface
