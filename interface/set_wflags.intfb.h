interface
 LOGICAL FUNCTION SET_WFLAGS(FLAG,N)
 use parkind_wave, only:&
 & jwim
 LOGICAL, INTENT(IN) :: FLAG(N)
 INTEGER(KIND=JWIM), INTENT(IN) :: N
 END FUNCTION SET_WFLAGS
end interface
