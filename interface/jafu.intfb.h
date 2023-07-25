interface
INTEGER(KIND=JWIM) FUNCTION JAFU (CL, J, IAN)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 REAL(KIND=JWRB), INTENT(IN):: CL
 INTEGER(KIND=JWIM), INTENT(IN) :: J, IAN
END FUNCTION JAFU
end interface
