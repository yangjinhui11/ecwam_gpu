interface
REAL(KIND=JWRB) FUNCTION AKI_ICE (G,XK,DEPTH,RHOW,CITH)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB), INTENT(IN) :: G, XK, DEPTH, RHOW, CITH
END FUNCTION AKI_ICE
end interface
