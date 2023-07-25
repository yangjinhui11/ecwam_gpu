interface
 SUBROUTINE MAKEGRID (BLOCK, GRID, PMISS)
 use parkind_wave, only:&
 & jwrb
 use yowparam , only:&
 & ngx,&
 & ngy
 REAL(KIND=JWRB), INTENT(INOUT) :: BLOCK(:)
 REAL(KIND=JWRB), INTENT(IN) :: PMISS
 REAL(KIND=JWRB), INTENT(OUT) :: GRID(NGX,NGY)
 END SUBROUTINE MAKEGRID
end interface
