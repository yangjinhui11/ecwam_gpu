interface
 SUBROUTINE SPR (NANG, THETAQ, THETA, ST)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 INTEGER(KIND=JWIM), INTENT(IN) :: NANG
 REAL(KIND=JWRB), INTENT(IN) :: THETAQ
 REAL(KIND=JWRB), DIMENSION(NANG), INTENT(IN) :: THETA
 REAL(KIND=JWRB), DIMENSION(NANG), INTENT(OUT) :: ST
 END SUBROUTINE SPR
end interface
