interface
 SUBROUTINE INITIALINT (IU06, NCA, NRA, &
 & NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO, &
 & AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC, &
 & ILONRGG, &
 & LLINTERPOL,DK1, DII1, DIIP1, KK, II, IIP)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 INTEGER(KIND=JWIM), INTENT(IN) :: IU06, NCA, NRA
 INTEGER(KIND=JWIM), INTENT(IN) :: NGX, NGY, KRGG
 INTEGER(KIND=JWIM), INTENT(INOUT) :: IPERIODIC
 INTEGER(KIND=JWIM), DIMENSION(NGY), INTENT(IN) :: KLONRGG
 INTEGER(KIND=JWIM), DIMENSION(NGY), INTENT(INOUT) :: KK
 INTEGER(KIND=JWIM), DIMENSION(NRA), INTENT(INOUT) :: ILONRGG
 INTEGER(KIND=JWIM), DIMENSION(NGX, NGY), INTENT(INOUT) :: II, IIP
 REAL(KIND=JWRB), INTENT(IN) :: XDELLA
 REAL(KIND=JWRB), DIMENSION(NGY), INTENT(IN) :: ZDELLO
 REAL(KIND=JWRB), INTENT(IN) :: AMOWEP, AMOSOP, AMOEAP, AMONOP
 REAL(KIND=JWRB), DIMENSION(NGY), INTENT(INOUT) :: DK1
 REAL(KIND=JWRB), DIMENSION(NGX, NGY), INTENT(INOUT) :: DII1, DIIP1
 LOGICAL, INTENT(INOUT) :: LLINTERPOL
 END SUBROUTINE INITIALINT
end interface