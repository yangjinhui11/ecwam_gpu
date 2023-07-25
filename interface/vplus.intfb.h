interface
REAL(KIND=JWRB) FUNCTION VPLUS(XI,XJ,XK,THI,THJ,THK)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB), INTENT(IN) :: XI, XJ, XK, THI, THJ, THK
END FUNCTION VPLUS
end interface
