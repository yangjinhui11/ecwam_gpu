interface
SUBROUTINE MPDECOMP(NPR, MAXLEN, LLIRANK, LLWVENVI)
 use parkind_wave, only:&
 & jwim
 INTEGER(KIND=JWIM), INTENT(IN) :: NPR
 INTEGER(KIND=JWIM), INTENT(OUT) :: MAXLEN
 LOGICAL, INTENT(IN) :: LLIRANK
 LOGICAL, INTENT(IN) :: LLWVENVI
END SUBROUTINE MPDECOMP
end interface
