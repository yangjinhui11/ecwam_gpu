interface
 SUBROUTINE MBLOCK (BATHY, KA, KE, IPP)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 use yowparam , only:&
 & ngx,&
 & ngy
 INTEGER(KIND=JWIM) :: KA, KE
 INTEGER(KIND=JWIM) :: IPP(NGY)
 REAL(KIND=JWRB) :: BATHY(NGX, NGY)
 END SUBROUTINE MBLOCK
end interface
