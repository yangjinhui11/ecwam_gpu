interface
REAL(KIND=JWRB) FUNCTION TRANSF(XK,D)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB), INTENT(IN) :: XK,D
END FUNCTION TRANSF
end interface
