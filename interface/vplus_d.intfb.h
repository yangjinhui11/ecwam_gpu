interface
REAL(KIND=JWRB) FUNCTION VPLUS_D(XI,XJ,XK,XIJ,XIK,XJK,XOI,XOJ,XOK)
 use parkind_wave, only:&
 & jwrb
 REAL(KIND=JWRB), INTENT(IN) :: XI,XJ,XK,XIJ,XIK,XJK,XOI,XOJ,XOK
END FUNCTION VPLUS_D
end interface
