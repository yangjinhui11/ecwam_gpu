interface
 SUBROUTINE CHECKCFL (KIJS, KIJL, DEPTH, DTC, &
 & CFLEA,CFLWE,CFLNO,CFLSO,CFLNO2, &
 & CFLSO2, CFLTP,CFLTM,CFLOP,CFLOM)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH
 REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: DTC
 REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: CFLEA,CFLWE,CFLNO,CFLSO,CFLNO2
 REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: CFLSO2,CFLTP,CFLTM,CFLOP,CFLOM
 END SUBROUTINE CHECKCFL
end interface
