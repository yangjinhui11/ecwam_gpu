interface
 SUBROUTINE TOPOAR (IU01, BATHY)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWPARAM , ONLY : NGX ,NGY
 INTEGER(KIND=JWIM), INTENT(IN) :: IU01
 REAL(KIND=JWRB), DIMENSION(NGX,NGY), INTENT(OUT) :: BATHY
 END SUBROUTINE TOPOAR
end interface
